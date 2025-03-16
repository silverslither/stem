from fractions import Fraction
from itertools import combinations
from math import isqrt
import sys

sys.set_int_max_str_digits(65536)

# returns isqrt, with optional check for exactness; if not exact, returns -1
def fraction_isqrt(x: Fraction, exact = False):
    y = Fraction(isqrt(x.numerator), isqrt(x.denominator))
    if exact and y * y != x:
        return Fraction(-1)
    return y

# returns [r, C, (x, y, z)] or None if input is not a valid Heronian triangle
def tangent_parameterization(sides: list[Fraction]):
    [a, b, c] = sides
    if a == 0 or b == 0 or c == 0 or a + b <= c or b + c <= a or c + a <= b:
        return None

    s = (a + b + c) / 2
    A = fraction_isqrt(s * (s - a) * (s - b) * (s - c), True)
    if A == -1:
        return None

    r = A / s
    C = s * s / A

    x = fraction_isqrt(-((a + b - c) * (a + b + c)) / ((a - b - c) * (a - b + c)))
    y = fraction_isqrt(-((b + c - a) * (b + c + a)) / ((b - c - a) * (b - c + a)))
    z = fraction_isqrt(-((c + a - b) * (c + a + b)) / ((c - a - b) * (c - a + b)))

    return [r, C, tuple(sorted([x, y, z]))]

# implicit differentiation of (x + y)(1 + 1/(xy - 1)) = C gives dy/dx = (y(x(2 - xy) + y))/(x(x(y^2 - 1) - 2y))
# returns (a, b) where the equation of tangent line is y = ax + b or None if slope is not valid
def find_tangent_expr(p: tuple[Fraction, Fraction]):
    (x, y) = p
    m = y * (x * (2 - x * y) + y)
    n = x * (x * (y * y - 1) - 2 * y)
    if m == 0 or n == 0 or m == -n:
        return None
    a = m / n
    return (a, -a * x + y)

# returns (a, b) where the equation of secant line is y = ax + b or None if slope is not valid
def find_secant_expr(p1: tuple[Fraction, Fraction], p2: tuple[Fraction, Fraction]):
    (x1, y1) = p1
    (x2, y2) = p2
    m = y2 - y1
    n = x2 - x1
    if m == 0 or n == 0 or m == -n:
        return None
    a = m / n
    return (a, y1 - a * x1)

# if (x + y)(1 + 1/(xy - 1)) = C and y = ax + b, (a^2 + a)x^3 + (2ab - aC + b)x^2 + (b^2 - bC)x + C = 0
# hence if solutions p, q (not necessarily unique) are known, the new solution is x = -C/(pqa(a + 1)), y = -C/(pq(a + 1)) + b
def get_new_solution(C: Fraction, line: tuple[Fraction, Fraction], solutions: tuple[Fraction, Fraction]):
    temp = -C / (solutions[0] * solutions[1] * (line[0] + 1))
    return (temp / line[0], temp + line[1])

def main() -> int:
    r = Fraction(-1)
    C = Fraction(-1)
    basis = []
    MAX_DEPTH = 18

    for i in range(1, len(sys.argv)):
        if "," in sys.argv[i]:
            params = tangent_parameterization([Fraction(x) for x in sys.argv[i].split(",")])
            if params == None:
                print(f"invalid Heronian triangle '{sys.argv[i]}', exiting")
                return 1
            [_r, _C, ratios] = params
            basis.extend([(ratios[0], ratios[1]), (ratios[0], ratios[2])])
            if (r != -1 and r != _r) or (C != -1 and C != _C):
                print(f"triangle '{sys.argv[i]}' has mismatched properties, exiting")
                return 1
            r = _r
            C = _C
        else:
            MAX_DEPTH = int(sys.argv[i])

    if r == -1 or C == -1:
        print("usage: solver.py <a,b,c> [d,e,f] ... [depth]")
        return 1

    solutions = [basis]
    visited = set(basis)

    # tangent step is used to reduce runtime, as it reduces required number of basis points
    if MAX_DEPTH > 0:
        iteration = []
        for node in basis:
            tangent = find_tangent_expr(node)
            if tangent == None:
                continue
            solution = get_new_solution(C, tangent, (node[0], node[0]))

            if solution[0] > solution[1]:
                solution = (solution[1], solution[0])
            if solution in visited:
                continue
            visited.add(solution)

            iteration.append(solution)
        for (node, mul) in combinations(basis, 2):
            secant = find_secant_expr(node, mul)
            if secant == None:
                continue
            solution = get_new_solution(C, secant, (node[0], mul[0]))

            if solution[0] > solution[1]:
                solution = (solution[1], solution[0])
            if solution in visited:
                continue
            visited.add(solution)

            iteration.append(solution)
        solutions.append(iteration)

    for i in range(1, MAX_DEPTH):
        iteration = []
        for node in solutions[i]:
            for mul in basis:
                secant = find_secant_expr(node, mul)
                if secant == None:
                    continue
                solution = get_new_solution(C, secant, (node[0], mul[0]))

                if solution[0] > solution[1]:
                    solution = (solution[1], solution[0])
                if solution in visited:
                    continue
                visited.add(solution)

                iteration.append(solution)
        solutions.append(iteration)

    final = set()
    for temp in solutions:
        for node in temp:
            (x, y) = node
            z = (x + y) / (x * y - 1)
            if x > 0 and y > 0 and z > 0:
                x *= r
                y *= r
                z *= r
                final.add(tuple(sorted([x + y, y + z, z + x])))

    triangles = list(final)
    triangles.sort(key = lambda x: x[0].denominator + x[1].denominator + x[2].denominator)
    for node in triangles:
        print(", ".join([f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator) for x in node]))
    print(f"count: {len(triangles)}")

    return 0

if __name__ == "__main__":
    exit(main())
