from sympy import *

init_printing(use_unicode = True)

n, h1, h2, x = symbols("n h1 h2 x")
expr = -h2/(2 * sqrt(n)) * x**2 + (1+h2/sqrt(n)) * h1/sqrt(n) * x - 1/2 * (1+ h2/sqrt(n)) * h1**2 / n

print("original:")
pprint(expr)

expr2 = 1+expr + expr**2 / 2

expr2 = (1 +  expr2 * ( 1+ h2/(2*sqrt(n)) - h2 ** 2 / (8*n) ) )/2

expr2 = expr2 - 1 - (expr2 - 1)**2 / 2
expr2 = expand(expr2)
expr2 = collect(expr2, x)
expr2 = collect(expr2, n)

print("result:")
pprint(expr2)


k, p, p0, x = symbols("k p p0 x")
expr = 2**k * (1-p)**(k-x) * p ** x

print("result:")
pprint(expr)
expr = expr.series(p, S(1)/2, 3).removeO()
expr = simplify(expr)


print("result:")
pprint(expr)
expr = 1 + (4* x - 2* k) * (p - S(1)/2) + 2 * ((2*x -k)**2 -k ) * (p-S(1)/2)**2
print("result:")
pprint(expr)
