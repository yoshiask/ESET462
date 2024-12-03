from sympy import init_printing, symbols, pprint, zeros, reduce_inequalities, collect
from sympy import Matrix, Symbol
from numpy import ceil

init_printing()

Kp: Symbol = symbols('K_p')
Ki: Symbol = symbols('K_i')
Kd: Symbol = symbols('K_d')
K: Symbol = symbols('k')

d_s = [2*Kp+5, 1-Kp, -Kp]

degree = len(d_s)

num_columns = int(ceil(degree / 2))
if degree % 2 == 0:
    num_columns += 1

table: Matrix = zeros(degree, num_columns)

# Populate first two rows
for c in range(degree):
    col_idx = c // 2
    coeff = d_s[c]
    if c % 2 == 0:
        table[0, col_idx] = coeff
    else:
        table[1, col_idx] = coeff

pprint(table)

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
print()
pprint(table)

conditions = []
for r in range(degree):
    expr = table[r, 0]
    condition = expr > 0
    try:
        bool(condition)
    except:
        conditions.append(condition)

print()
pprint(conditions)

print()
pprint(reduce_inequalities(conditions, Ki).simplify())
