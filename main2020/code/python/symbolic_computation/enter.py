from sympy import *

init_printing(use_unicode = True)

x, k, u, w, t1, t2, t3 = symbols("x k u w t1 t2 t3")
expr1 = (k-x) * ln(1-t2/(1-u))
expr2 = x * ln(1 + t2/u)
expr3 = ln(1+(w+t1) * ( (1-t3/(1-u-t2))**(k-x) * (1+t3/(u + t2)) ** x - 1 ))
expr = expr1 + expr2 + expr3
a1 = diff(expr, t1)
a2 = diff(expr, t2)
a3 = diff(expr, t3)
print("a1")
pprint(a1.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a2")
pprint(a2.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a3")
pprint(a3.subs([ (t1,0), (t2,0), (t3,0) ]))


a11 = diff(a1, t1)
a12 = diff(a1, t2)
a13 = diff(a1, t3)


a21 = diff(a2, t1)
a22 = diff(a2, t2)
a23 = diff(a2, t3)

a31 = diff(a3, t1)
a32 = diff(a3, t2)
a33 = diff(a3, t3)


print("a11:")
pprint(a11.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a12:")
pprint(a12.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a13:")
pprint(a13.subs([ (t1,0), (t2,0), (t3,0) ]))

print("a21:")
pprint(a21.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a22:")
pprint(a22.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a23:")
pprint(a23.subs([ (t1,0), (t2,0), (t3,0) ]))

print("a31:")
pprint(a31.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a32:")
pprint(a32.subs([ (t1,0), (t2,0), (t3,0) ]))
print("a33:")
pprint(a33.subs([ (t1,0), (t2,0), (t3,0) ]))

myA33 = w * ( ((1-w)*(x-u*k)**2 + (2*u-1) * x - k* u**2 ) / (u**2 * (1-u) ** 2) )

shouldBeZero = a33.subs([ (t1,0), (t2,0), (t3,0) ]) - myA33

print("should be zero:")
pprint(simplify(shouldBeZero))
