from sympy import collect, diff, sympify
from sympy import Basic, Symbol, Function, Add, Eq

t: Symbol = Symbol('t')
s: Symbol = Symbol('s')
z: Symbol = Symbol('z')

x: Function = Function('x')
y: Function = Function('y')

xt = x(t)
yt = y(t)

kp: Symbol = Symbol('K_p')
ki: Symbol = Symbol('K_i')
kd: Symbol = Symbol('K_d')
k: Symbol = Symbol('k')

avail_funcs = [x, y]


def diff_n(f: Function, n: int) -> Function:
    result = f
    for i in range(n):
        result = diff(result, t)
    return result


def tdom_to_sdom(sys: Function) -> Add:
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

    y_term = left_s
    x_term = right_s
    if right_control == y(s):
        y_term, x_term = x_term, y_term

    return y_term / x_term


def tdom_to_tf(tdom_func: Function) -> Function:
    s_func = tdom_to_sdom(tdom_func)
    tf = sdom_to_tf(s_func)
    return tf


def tf_to_tdom(tf: Function) -> Eq:
    if not tf.is_Mul:
        raise ValueError()

    a, b = tf.args

    if a.is_Pow and a.args[1] < 0:
        denominator = 1 / a
        numerator = b
    else:
        denominator = 1 / b
        numerator = a

    numerator = numerator.expand()
    denominator = denominator.expand()

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