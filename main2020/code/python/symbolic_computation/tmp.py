import sympy as sp
from sympy import *

init_printing(use_unicode = True)

u, u0 = symbols("u u0")
expr = u* (1-u) * ( sqrt((1-u) * u0) -sqrt((1-u0) * u ) ) ** 2
a1 = diff(expr, u)
a2 = diff(a1, u)

shouldBeZero = a2.subs([ (u,u0)])

print("should be zero:")
pprint(simplify(shouldBeZero))



z, t, M = symbols("z t M")
expr = 2**t * (1 + sp.exp(z * t * M)) / (1 + sp.exp(z * M ))**t
expr = expr.subs([(M, 1)])
a1 = diff(expr, z)
a2 = diff(a1, z)
a3 = diff(a2, z)
a4 = diff(a3, z)


print("first derivative:")
pprint(simplify(a1))
pprint(simplify(a1.subs([(z, 0), (M, 1)])))
print("second derivative:")
pprint(simplify(a2))
pprint(simplify(a2.subs([(z, 0), (M, 1)])))
print("third derivative:")
pprint(simplify(a3))
pprint(factor(a3))
pprint(simplify(a3.subs([(z, 0), (M, 1)])))

print("fourth derivative:")
pprint(simplify(a4))
pprint(factor(a4))
pprint(simplify(a4.subs([(z, 0), (M, 1)])))
