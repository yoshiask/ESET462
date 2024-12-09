from sympy import init_printing, pprint, srepr, Heaviside, symbols
from ctl_sys import *

init_printing()

Gs = 1 / (s**6 + 2*s**5 + 8*s**4 + 20*s**2 + 16*s + 16)
coeffs = char_eq_coeffs_from_tf(Gs)
pprint(coeffs)

table = routhhurwitz_table(coeffs)
pprint(table)