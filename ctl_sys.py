from sympy import (collect, diff, sympify, ceiling, zeros, reduce_inequalities, solve, nan, ln,
                   laplace_transform, inverse_laplace_transform, oo, summation, sqrt, exp, pi, limit)
from sympy import Basic, Symbol, Function, Add, Eq, Matrix, Mul, Expr, Poly, Abs, Heaviside, DiracDelta

t: Symbol = Symbol('t', positive=True, real=True)
s: Symbol = Symbol('s', positive=True)

n: Symbol = Symbol('n', positive=True, real=True)
z: Symbol = Symbol('z')

x: Function = Function('x')
y: Function = Function('y')

xt = x(t)
yt = y(t)

kp: Symbol = Symbol('K_p', real=True)
ki: Symbol = Symbol('K_i', real=True)
kd: Symbol = Symbol('K_d', real=True)
k: Symbol = Symbol('k', real=True)

avail_funcs = [x, y]


def diff_n(f: Function, n: int) -> Function:
    result = f
    for i in range(n):
        result = diff(result, t)
    return result


def tdom_to_sdom(sys: Function) -> Add:
    terms = []
    if sys.is_Equality:
        terms.extend(sys.args[0].args)
        terms.extend([-arg for arg in sys.args[1].args])
    else:
        terms = sys.args

    s_terms = []
    for term in terms:
        coeff = sympify(1)
        degree = sympify(0)

        if term.is_Mul:
            coeff, term = term.args

        if term.is_Derivative:
            degree = term.args[1][1]
            term = term.args[0]

        control = term.subs(t, s)

        s_terms.append(coeff * control * s ** degree)

    return Add(*s_terms)


def sdom_to_tf(s_func: Function) -> Function:
    collected_sdom: Function = collect(s_func, [v(s) for v in avail_funcs])
    if not collected_sdom.is_Add:
        raise ValueError()

    # Extract the Y and X coefficients, then return Y/X
    left, right = collected_sdom.args
    if not left.is_Mul or not right.is_Mul:
        raise ValueError()

    left_s, left_control = left.args
    right_s, right_control = right.args

    # Undo cross-multiplication
    if right_control == y(s):
        Ys = left_s
        Xs = right_s
    else:
        Ys = right_s
        Xs = left_s

    Xs *= -1    # Bring X to other side

    return Ys / Xs


def tdom_to_tf(tdom_func: Function) -> Function:
    s_func = tdom_to_sdom(tdom_func)
    tf = sdom_to_tf(s_func)
    return tf


def get_numerator_and_denominator(func: Expr) -> tuple[Symbol, Symbol]:
    func = func.expand(numer=True, denom=True)
    one = sympify(1)

    if func.is_Pow:
        b, p = func.args
        if p > 0:
            return b, one
        else:
            return one,  b
    if not func.is_Mul:
        return func, one

    numerator = one
    denominator = one
    for arg in func.args:
        if arg.is_Pow and arg.args[1] < 0:
            denominator *= arg.args[0]
        else:
            numerator *= arg

    return numerator.expand(), denominator.expand()


def tf_to_tdom(tf: Function) -> Eq:
    numerator, denominator = get_numerator_and_denominator(tf)

    def _sdom_to_tdom(term: Basic, control: Function, tdom_terms: list[Basic]):
        if term.is_Add:
            for subterm in term.args:
                _sdom_to_tdom(subterm, control, tdom_terms)
        else:
            coeff = sympify(1)
            degree = sympify(0)
            if term.is_Number:
                coeff = term
            else:
                if term.is_Mul:
                    coeff, term = term.args
                    degree = 1
                if term.is_Pow:
                    degree = term.args[1]

            tdom_term = coeff * diff_n(control, degree)
            tdom_terms.append(tdom_term)

        return tdom_terms

    y_terms = []
    _sdom_to_tdom(denominator, yt, y_terms)
    y_t = Add(*y_terms)

    x_terms = []
    _sdom_to_tdom(numerator, xt, x_terms)
    x_t = Add(*x_terms)

    sys_t = Eq(x_t, y_t)
    return sys_t


def sdom_to_zdom(sdom_func: Expr) -> Expr:
    z2s = (s+1) / (s-1)
    return sdom_func.subs(z, z2s).simplify()


def zdom_to_sdom(zdom_func: Expr) -> Expr:
    s2z = (z-1) / (z+1)
    return zdom_func.subs(s, s2z).simplify()


def routhhurwitz_table(char_eq_coeffs: list[Basic]) -> Matrix:
    degree = len(char_eq_coeffs)

    num_columns = int(ceiling(degree / 2))
    if degree % 2 == 0:
        num_columns += 1

    table: Matrix = zeros(degree, num_columns)

    # Populate first two rows
    for c in range(degree):
        col_idx = c // 2
        coeff = char_eq_coeffs[c]
        if c % 2 == 0:
            table[0, col_idx] = coeff
        else:
            table[1, col_idx] = coeff

    # Perform criterion
    for r in range(2, degree):
        divisor = table[r - 1, 0]

        tl = table[r - 2, 0]
        bl = divisor

        for c in range(num_columns):
            if c >= num_columns - 1:
                tr = 0
                br = 0
            else:
                tr = table[r - 2, c + 1]
                br = table[r - 1, c + 1]
            m = Matrix([[tl, tr], [bl, br]])
            table[r, c] = -m.det() / divisor

    table.simplify()
    return table


def routhhurwitz_criterion(table: Matrix) -> any:
    conditions = []
    for expr in table.col(0):
        if expr == nan:
            raise ValueError("Marginally stable")
        condition = expr > 0
        try:
            bool(condition)
        except:
            conditions.append(condition)

    return reduce_inequalities(conditions, [kp, ki, kd, k]).simplify()


def routhhurwitz_complete(char_eq_coeffs: list[Basic]) -> any:
    table = routhhurwitz_table(char_eq_coeffs)
    return routhhurwitz_criterion(table)


def tf_poles_zeros(tf: Function) -> tuple[list[Basic], list[Basic]]:
    numerator, denominator = get_numerator_and_denominator(tf)
    return solve(numerator), solve(denominator)


def zdom_bibo_stable(zdom_func: Expr) -> Expr:
    sdom_func = zdom_to_sdom(zdom_func)
    _, denominator = get_numerator_and_denominator(sdom_func)

    # Extract coefficients of characteristic equation
    char_eq = Poly(denominator)
    char_eq_coeffs = char_eq.all_coeffs()

    return routhhurwitz_complete(char_eq_coeffs)


def zdom_stable(zdom_func: Expr) -> bool | Expr:
    _, denominator = get_numerator_and_denominator(zdom_func)
    print(denominator)
    char_eq = Poly(denominator)
    poles = solve(char_eq, z)

    conditions = []
    for expr in poles:
        condition = Abs(expr) < 1
        try:
            if not bool(condition):
                return False
        except:
            conditions.append(condition)

    return reduce_inequalities(conditions, [kp, ki, kd, k]).simplify()


def laplace(f: Expr) -> Expr:
    return laplace_transform(f, t, s, noconds=True)


def ilaplace(f: Expr) -> Expr:
    return inverse_laplace_transform(f, s, t, noconds=True)


def open_to_closed_loop(GOL: Expr) -> Expr:
    return GOL / (1 + GOL)


def build_closed_loop(Cs: Expr, Ps: Expr) -> Expr:
    return open_to_closed_loop(Cs * Ps).simplify()


def build_feedback(Fs: Expr, Hs: Expr) -> Expr:
    return Fs / (1 + Fs * Hs)


def z_transform(ndom_func: Expr) -> Expr:
    return summation(ndom_func * z**-n, (n, -oo, oo)).doit()


def heaviside(delay: Expr = 0) -> Heaviside:
    """
    Creates a time-domain unit step function.
    """
    return Heaviside(t - delay)


def get_damping_form(tf: Expr) -> tuple[Expr, Expr]:
    b, denominator = get_numerator_and_denominator(tf)
    poly: Poly = Poly(denominator, s)
    if poly.degree() != 2:
        raise ValueError

    a2, a1, a0 = poly.coeffs()

    omega_n = sqrt(sympify(a0) / sympify(a2))
    zeta_n = (a1 / (2 * a2 * omega_n))
    return omega_n, zeta_n


def poles_from_damping(omega_n: Expr, zeta_n: Expr) -> tuple[Expr, Expr]:
    a = -zeta_n * omega_n
    b = omega_n * sqrt(zeta_n**2 - 1)
    return a+b, a-b


def max_overshoot(zeta_n: Expr) -> Expr:
    return exp(-pi * zeta_n / sqrt(1 - zeta_n**2)).simplify()


def delay_time(omega_n: Expr, zeta_n: Expr) -> Expr:
    return ((1 + 7*zeta_n/10) / omega_n).simplify()


def rise_time(omega_n: Expr, zeta_n: Expr) -> Expr:
    return ((8/10 + 5*zeta_n/2) / omega_n).simplify()


def settling_time(omega_n: Expr, zeta_n: Expr, delta_E: Expr = sympify(0.05)) -> Expr:
    return (-ln(delta_E * sqrt(1 - zeta_n**2)) / (omega_n * zeta_n)).simplify()


def steady_state_error(Gs: Expr) -> Expr:
    return final_value_theorem(Gs, 1/s)


def final_value_theorem(Gs: Expr, Xs: Expr) -> Expr:
    Ys = Gs * Xs
    return limit(s * Ys, s, 0)


def initial_value_theorem(Gs: Expr, Xs: Expr) -> Expr:
    Ys = Gs * Xs
    return limit(s * Ys, s, oo)