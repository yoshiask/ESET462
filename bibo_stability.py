from sympy import init_printing, pprint
from ctl_sys import *

init_printing()

d_s = [2*kp+5, 1-kp, -kp]

#pprint(routhhurwitz_complete(d_s))

#Gz = (2 - z**-1) / (3 + 5*z**-1 + 6*z**-2 + z**-3)
Gz = (2*z - 1) / ((3*z+1) * (2*z+k))
pprint(Gz)

cond = zdom_stable(Gz)
print()
pprint(cond)

cond = zdom_bibo_stable(Gz)
print()
pprint(cond)
