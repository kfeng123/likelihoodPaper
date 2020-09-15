from sympy import *
from sympy.functions.combinatorial.numbers import nC

init_printing(use_unicode = True)

x, w, p1, p2, p = symbols("x w p1 p2 p")

full_expr = ( (1-w) * (1- p1)**(3-x) * p1**x + w * (1-p2)**(3-x) * p2**x )

null_expr = full_expr.subs([(p1, p), (p2, p) , (w, 0) ])

sqrt_full_times_null = sqrt(full_expr * null_expr)

for i in range(4):
    if i == 0:
        theSum = sqrt_full_times_null.subs([(x, i)])
    else:
        theSum = theSum + sqrt_full_times_null.subs([(x, i)])

print("should be zero:")
pprint(simplify(theSum))
