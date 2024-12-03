from sympy import init_printing, pprint
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import det
import sys

init_printing()
np.set_printoptions(suppress=True)

d_s = [250, -500, 250]

degree = len(d_s)

num_columns = int(np.ceil(degree / 2))
if degree % 2 == 0:
    num_columns += 1

rows = np.zeros([degree, num_columns])

# Populate first two rows
for c in range(degree):
    col_idx = c // 2
    coeff = d_s[c]
    if c % 2 == 0:
        rows[0][col_idx] = coeff
    else:
        rows[1][col_idx] = coeff

pprint(rows)

for r in range(2, degree):
    divisor = rows[r - 1][0]
    
    tl = rows[r - 2][0]
    bl = divisor

    for c in range(num_columns):
        if c >= num_columns - 1:
            tr = 0
            br = 0
        else:
            tr = rows[r - 2][c + 1]
            br = rows[r - 1][c + 1]
        rows[r][c] = -det([[tl, tr], [bl, br]]) / divisor

print()
pprint(rows)

# Count sign changes
num_sign_changes = 0
for r in range(1, degree):
    is_previous_pos = rows[r - 1][0] >= 0
    is_current_pos = rows[r][0] >= 0
    if is_current_pos != is_previous_pos:
        num_sign_changes += 1

print(num_sign_changes)