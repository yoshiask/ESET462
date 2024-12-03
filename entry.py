from sympy import init_printing, pprint
from ctl_sys import *

init_printing()

tf = (2*s + 1) / (s**3 + 2*s**2 + 5)
pprint(tf)

tdom_func = tf_to_tdom(tf)
print()
pprint(tdom_func)
