from sympy import init_printing, pprint, srepr
from ctl_sys import *

init_printing()

tf = (s + 4) / ((s**2+1) * (s+2) * (s+3))
pprint(tf)

zeros, poles = tf_poles_zeros(tf)

print("\r\nzeros:")
pprint(zeros)
print("\r\npoles:")
pprint(poles)
