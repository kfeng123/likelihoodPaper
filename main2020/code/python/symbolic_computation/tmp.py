from sympy import *

init_printing(use_unicode = True)

u, u0 = symbols("u u0")
expr = u* (1-u) * ( sqrt((1-u) * u0) -sqrt((1-u0) * u ) ) ** 2
a1 = diff(expr, u)
a2 = diff(a1, u)

shouldBeZero = a2.subs([ (u,u0)])

print("should be zero:")
pprint(simplify(shouldBeZero))
