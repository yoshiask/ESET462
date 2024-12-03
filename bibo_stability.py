from sympy import init_printing, pprint
from ctl_sys import *

init_printing()

d_s = [2*kp+5, 1-kp, -kp]

pprint(routhhurwitz_complete(d_s))
