from collections.abc import Iterable, Generator
from fractions import Fraction
from math import isqrt, gcd, lcm
import time
import sys

sys.set_int_max_str_digits(65536)

def height(solution: tuple[Fraction, Fraction, Fraction]):
    return (Fraction(lcm(*[x.denominator for x in solution]), gcd(*[x.numerator for x in solution])), solution[0])

# returns (s, A, (x, y, z)) or None if input is not a valid triangle
def variable_change(sides: list[Fraction]):
    if len(sides) != 3:
        return None

    [a, b, c] = sides
    if a <= 0 or b <= 0 or c <= 0 or a + b <= c or b + c <= a or c + a <= b:
        return None

    s = (a + b + c) / 2
    A = s * (s - a) * (s - b) * (s - c)

    x = A / ((s - a) * (s - b))
    y = A / ((s - b) * (s - c))
    z = A / ((s - c) * (s - a))

    variables = sorted([x, y, z])

    return (s, A, (variables[0], variables[1], variables[2]))

# generate ordered coprime pairs using the tree (a, b) => (a + b, b), (a, a + b), starting at (1, 2)
def generate_pairs(max_depth: int) -> Generator[tuple[int, int], None, None]:
    yield (1, 1)
    if max_depth <= 0:
        return
    stack = [(1, 2, 0)]
    while len(stack):
        (a, b, depth) = stack.pop()
        if a > b:
            yield (b, a)
        else:
            yield (a, b)
        depth += 1
        if depth == max_depth:
            continue
        stack.append((a + b, b, depth))
        stack.append((a, a + b, depth))

# return a rational solution for y if (x + y)(1 + 1/(xy/A - 1)) = C if it exists, given x
# otherwise return -1 if no solution exists, or 0 if irrational solution(s) exist
# if (x + y)(1 + 1/(xy/A - 1)) = C, y = (+/-sqrt(x(x(C - x)^2 - 4AC)) + Cx - x^2)/(2x)
def get_rational_y(A: Fraction, C: Fraction, x: Fraction):
    d = x * (x * (C - x) ** 2 - 4 * A * C)
    if d < 0:
        return -1
    if d == 0:
        return (C * x - x * x)/(2 * x)

    m = isqrt(d.numerator)
    n = isqrt(d.denominator)
    if m * m != d.numerator or n * n != d.denominator:
        return 0
    return (-Fraction(m, n) + C * x - x * x)/(2 * x)

# searcher step, which tests up to 2^depth rational numbers between 0 and C
def searcher(C: Fraction, A: Fraction, depth: int):
    solutions: set[tuple[Fraction, Fraction, Fraction]] = set()
    upper_bound = C
    lower_bound = 0
    for (m, n) in generate_pairs(depth):
        x = Fraction(m, n)
        if x < lower_bound or x > upper_bound:
            continue

        y = get_rational_y(A, C, x * C)
        if y == -1:
            if x > Fraction(1, 3):
                upper_bound = x
            else:
                lower_bound = x
            continue
        if y == 0:
            continue

        x *= C
        z = (x + y) / (x * y / A - 1)
        variables = sorted([x, y, z])
        solutions.add((variables[0], variables[1], variables[2]))
    return solutions

# implicit differentiation of (x + y)(1 + 1/(xy/A - 1)) = C gives dy/dx = -(y(y(x^2 - A) - 2Ax))/(x(x(y^2 - A) - 2Ay))
# returns (a, b) where the equation of tangent line is y = ax + b or None if slope is not valid
def find_tangent_expr(point: tuple[Fraction, Fraction], A: Fraction):
    (x, y) = point
    m = y * (y * (x * x - A) - 2 * A * x)
    n = x * (x * (y * y - A) - 2 * A * y)
    if m == 0 or n == 0 or m == -n:
        return None
    a = -m / n
    return (a, -a * x + y)

# returns (a, b) where the equation of secant line is y = ax + b or None if slope is not valid
def find_secant_expr(point1: tuple[Fraction, Fraction], point2: tuple[Fraction, Fraction], A: Fraction):
    if point1 == point2: # handle tangent case automatically
        return find_tangent_expr(point1, A)
    (x1, y1) = point1
    (x2, y2) = point2
    m = y2 - y1
    n = x2 - x1
    if m == 0 or n == 0 or m == -n:
        return None
    a = m / n
    return (a, y1 - a * x1)

# return a new (reflected) solution to (x + y)(1 + 1/(xy/A - 1)) = C, given two existing ones
# if (x + y)(1 + 1/(xy/A - 1)) = C and y = ax + b, (a^2 + a)x^3 + (2ab - aC + b)x^2 + (b^2 - bC)x + AC
# hence if solutions p, q (not necessarily unique) are known, the new solution is x = -AC/(pqa(a + 1)), y = -AC/(pq(a + 1)) + b
def get_new_solution(C: Fraction, A: Fraction, line: tuple[Fraction, Fraction], solutions: tuple[Fraction, Fraction]):
    temp = -C * A / (solutions[0] * solutions[1] * (line[0] + 1))
    return (temp + line[1], temp / line[0])

def solver_recurse(C: Fraction, A: Fraction, basis: list[tuple[Fraction, Fraction]], depth: int, start: tuple[Fraction, Fraction] | None) -> Generator[tuple[Fraction, Fraction], None, None]:
    if len(basis) == 0 or depth == 0:
        if start != None:
            yield start
        return

    new_basis = basis[1:]
    yield from solver_recurse(C, A, new_basis, depth, start)

    for i in range(1, depth + 1):
        if start == None:
            start = basis[0]
        else:
            line = find_secant_expr(start, basis[0], A)
            if line == None:
                return
            start = get_new_solution(C, A, line, (start[0], basis[0][0]))
        yield from solver_recurse(C, A, new_basis, depth - i, start)

# solver step using chord and tangent methods up to a maximum Manhattan distance of depth in Z^n, with n basis points
def solver(C: Fraction, A: Fraction, basis: Iterable[tuple[Fraction, Fraction]], depth: int):
    depth += 1
    triples: set[tuple[Fraction, Fraction, Fraction]] = set()
    for (x, y) in solver_recurse(C, A, list(basis), depth, None):
        z = (x + y) / (x * y / A - 1)
        if x <= 0 or y <= 0 or z <= 0:
            continue
        variables = sorted([x, y, z])
        triples.add((variables[0], variables[1], variables[2]))
    return triples

# reducer step using solver function with a maximum Manhattan distance of depth
def reducer(C: Fraction, A: Fraction, inputs: Iterable[tuple[Fraction, Fraction, Fraction]], depth: int):
    sorted_inputs = sorted(inputs, key=height)

    superset = sorted_inputs
    for i in range(1, len(sorted_inputs)):
        result = solver(C, A, [(x[0], x[2]) for x in sorted_inputs[:i]], depth)
        if all([x in result for x in inputs]):
            superset = sorted_inputs[:i]
            break

    flag = True
    while flag:
        for i in range(len(superset) - 1, -1, -1):
            generators = superset.copy()
            del generators[i]
            result = solver(C, A, [(x[0], x[2]) for x in generators], depth)
            if all([x in result for x in inputs]):
                superset = generators
                flag = True
                break
        flag = False

    return superset

# convert elliptic curve solutions to triangles, and print them
def print_solutions(s: Fraction, solutions: Iterable[tuple[Fraction, Fraction, Fraction]], command_max: int = -1):
    triangles = []
    for (x, y, z) in solutions:
        x /= s
        y /= s
        z /= s
        triangles.append((x + y, x + z, y + z))
    triangles.sort(key=height)
    triangles = [[str(y) for y in x] for x in triangles]
    for point in triangles:
        print(", ".join(point))
    print(f"count: {len(triangles)}")
    if len(triangles) <= command_max:
        print("solver.py " + " ".join([",".join(x) for x in triangles]))

def main() -> int:
    start = time.perf_counter_ns()
    previous = start
    s = None
    A = None
    inputs: set[tuple[Fraction, Fraction, Fraction]] = set()
    depths = [12, 8, 8]
    verbose_printing = False

    i = 1
    while i < len(sys.argv) and "," in sys.argv[i]:
        sides = [Fraction(x) for x in sys.argv[i].split(",")]
        variables = variable_change(sides)
        if variables == None:
            print(f"invalid triangle '{sys.argv[i]}', exiting")
            return 1
        [_s, _A, variables] = variables
        inputs.add(variables)
        if (s != None and s != _s) or (A != None and A != _A):
            print(f"triangle '{sys.argv[i]}' has mismatched properties, exiting")
            return 1
        s = _s
        A = _A
        i += 1

    j = 0
    while i + j < len(sys.argv):
        if sys.argv[i + j] == "-v":
            verbose_printing = True
        elif j < 3:
            depths[j] = int(sys.argv[i + j])
        j += 1
    [searcher_depth, reducer_depth, solver_depth] = depths

    if s == None or A == None:
        print("usage: solver.py <a,b,c [d,e,f] ...> [search depth = 12] [reducer depth = 10] [solver depth = 8] [-v]")
        return 1
    C = s * s
    
    inputs = inputs.union(searcher(C, A, searcher_depth))
    if verbose_printing:
        print_solutions(s, inputs, 32)
        elapsed = time.perf_counter_ns() - previous
        previous += elapsed
        print(f"searcher finished in {elapsed // 1000000} ms\n")

    generators = reducer(C, A, inputs, reducer_depth)
    if verbose_printing:
        print_solutions(s, generators, 32)
        elapsed = time.perf_counter_ns() - previous
        previous += elapsed
        print(f"reducer finished in {elapsed // 1000000} ms\n")

    basis = [(x[0], x[2]) for x in generators]
    solutions = solver(C, A, basis, solver_depth)
    print_solutions(s, solutions)
    if verbose_printing:
        elapsed = time.perf_counter_ns() - previous
        previous += elapsed
        print(f"solver finished in {elapsed // 1000000} ms")
        elapsed = previous - start
        print(f"script finished in {elapsed // 1000000} ms")

    return 0

if __name__ == "__main__":
    exit(main())
