from sympy import *
import math
import sys

N = int(sys.argv[1])

A = [0, 1, 1]
for i in range(2, N + 1):
    A.append(simplify(A[i] + sympify(f"x ** 2 / (4 * (4 * {i} ** 2 - 1))") * A[i - 1]))

D = reversed(Poly(A[-1]).coeffs())

x = Symbol("x", real=True)
k = Symbol("k", integer=True)

base = x - k + Rational(N + 1, 2)
poly = Piecewise((base**N, base >= 0), (0, True))

term = (-1) ** k * binomial(N + 1, k) * poly
bspline = simplify(sum(term.subs(k, i) for i in range((N + 2) // 2)) / factorial(N))

omoms = 0
for i in D:
    omoms += i * bspline
    bspline = diff(diff(bspline, x), x)

omoms = expand(simplify(omoms.subs(x, -x)))

C = 1
for expr, _ in omoms.args:
    try:
        C = math.lcm(C, *[denom(x) for x in Poly(expr).coeffs()])
    except:
        pass


def format(expr, cond, last):
    expr = horner(C * expr)
    expr = str(expr.xreplace({n: float(n) for n in expr.atoms(Number)}))
    if not last:
        print(f"    if ({cond.xreplace({n: float(n) for n in cond.atoms(Number)})})")
        print(f"        return {expr.replace('*', ' * ')};")
    else:
        print(f"    return {expr.replace('*', ' * ')};")


print(f"// Normalized to {C}.0")
print(f"double OMOMS{N}(double x) {{")
for expr, cond in omoms.args[:-1]:
    format(expr, cond, expr == omoms.args[-2][0])
print("}")

mpstr = []
for i in range(N // 2 + 1):
    v = simplify(omoms.subs(x, i))
    mpstr.append(f"mp.mpf({numer(v)}) / {denom(v)}")

mpstr = ", ".join(list(reversed(mpstr[1:])) + mpstr)
print(f"\nD = [{mpstr}] # omoms{N}")
