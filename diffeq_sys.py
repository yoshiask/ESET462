from sympy import init_printing, pprint
from ctl_sys import *

init_printing()

sys_l = diff_n(yt, 3) + 2*diff_n(yt, 2) + 3*diff_n(yt, 1) + 2*yt
sys_r = diff_n(xt, 1) + xt

sys_c = sys_l - sys_r
print("System (differential eq.):")
pprint(sys_c)

print("\r\nSystem transfer function (s-domain):")
tf = diffeq_to_tf(sys_c)
pprint(tf)