from collections.abc import Iterable, Generator
from fractions import Fraction
from math import isqrt, gcd, lcm
import time
import sys

sys.set_int_max_str_digits(65536)


def height(solution: tuple[Fraction, Fraction, Fraction]):
    return (Fraction(lcm(*[x.denominator for x in solution]), gcd(*[x.numerator for x in solution])), solution[0])


# returns (s, r, (x, y, z)) or None if input is not a valid triangle
def variable_change(sides: list[Fraction]):
    if len(sides) != 3:
        return None

    [a, b, c] = sides
    if a <= 0 or b <= 0 or c <= 0 or a + b <= c or b + c <= a or c + a <= b:
        return None

    s = (a + b + c) / 2
    x = s - a
    y = s - b
    z = s - c

    A = s * x * y * z
    r = A / (s * s)

    variables = sorted([x, y, z])

    return s, r, (variables[0], variables[1], variables[2])


# generate ordered coprime pairs using the tree (a, b) => (a + b, b), (a, a + b), starting at (1, 2)
def generate_pairs(max_depth: int) -> Generator[tuple[int, int], None, None]:
    yield (1, 1)
    if max_depth <= 0:
        return
    stack = [(1, 2, 0)]
    while len(stack):
        a, b, depth = stack.pop()
        if a > b:
            yield (b, a)
        else:
            yield (a, b)
        depth += 1
        if depth == max_depth:
            continue
        stack.append((a + b, b, depth))
        stack.append((a, a + b, depth))


# return a rational solution for y if (x + y)(1 + 1/(xy/r - 1)) = s if it exists, given x
# otherwise return -1 if no solution exists, or 0 if irrational solution(s) exist
# if (x + y)(1 + 1/(xy/r - 1)) = s, y = (+/-sqrt(x(x(s - x)^2 - 4rs)) + sx - x^2)/(2x)
def get_rational_y(s: Fraction, r: Fraction, x: Fraction):
    d = x * (x * (s - x) ** 2 - 4 * r * s)
    if d < 0:
        return -1
    if d == 0:
        return (s * x - x * x) / (2 * x)

    m = isqrt(d.numerator)
    n = isqrt(d.denominator)
    if m * m != d.numerator or n * n != d.denominator:
        return 0
    return (-Fraction(m, n) + s * x - x * x) / (2 * x)


# searcher step, which tests up to 2^depth rational numbers between 0 and s
def searcher(s: Fraction, r: Fraction, depth: int):
    solutions: set[tuple[Fraction, Fraction, Fraction]] = set()
    upper_bound = s
    lower_bound = 0
    for m, n in generate_pairs(depth):
        x = Fraction(m, n)
        if x < lower_bound or x > upper_bound:
            continue

        y = get_rational_y(s, r, x * s)
        if y == -1:
            if x > Fraction(1, 3):
                upper_bound = x
            else:
                lower_bound = x
            continue
        if y == 0:
            continue

        x *= s
        z = (x + y) / (x * y / r - 1)
        variables = sorted([x, y, z])
        solutions.add((variables[0], variables[1], variables[2]))
    return solutions


# implicit differentiation of (x + y)(1 + 1/(xy/r - 1)) = s gives dy/dx = -(y(y(x^2 - r) - 2rx))/(x(x(y^2 - r) - 2ry))
# returns (a, b) where the equation of tangent line is y = ax + b or None if slope is not valid
def find_tangent_expr(point: tuple[Fraction, Fraction], r: Fraction):
    x, y = point
    m = y * (y * (x * x - r) - 2 * r * x)
    n = x * (x * (y * y - r) - 2 * r * y)
    if m == 0 or n == 0 or m == -n:
        return None
    a = -m / n
    return (a, -a * x + y)


# returns (a, b) where the equation of secant line is y = ax + b or None if slope is not valid
def find_secant_expr(point1: tuple[Fraction, Fraction], point2: tuple[Fraction, Fraction], r: Fraction):
    if point1 == point2:  # handle tangent case automatically
        return find_tangent_expr(point1, r)
    x1, y1 = point1
    x2, y2 = point2
    m = y2 - y1
    n = x2 - x1
    if m == 0 or n == 0 or m == -n:
        return None
    a = m / n
    return (a, y1 - a * x1)


# return a new (reflected) solution to (x + y)(1 + 1/(xy/r - 1)) = s, given two existing ones
# if (x + y)(1 + 1/(xy/r - 1)) = s and y = ax + b, (a^2 + a)x^3 + (2ab - as + b)x^2 + (b^2 - bs)x + rs = 0
# hence if solutions p, q (not necessarily unique) are known, the new solution is x = -rs/(pqa(a + 1)), y = -rs/(pq(a + 1)) + b
def get_new_solution(s: Fraction, r: Fraction, line: tuple[Fraction, Fraction], solutions: tuple[Fraction, Fraction]):
    temp = -s * r / (solutions[0] * solutions[1] * (line[0] + 1))
    return (temp + line[1], temp / line[0])


def solver_recurse(s: Fraction, r: Fraction, basis: list[tuple[Fraction, Fraction]], depth: int, start: tuple[Fraction, Fraction] | None) -> Generator[tuple[Fraction, Fraction], None, None]:
    if len(basis) == 0 or depth == 0:
        if start != None:
            yield start
        return

    new_basis = basis[1:]
    yield from solver_recurse(s, r, new_basis, depth, start)

    for i in range(1, depth + 1):
        if start == None:
            start = basis[0]
        else:
            line = find_secant_expr(start, basis[0], r)
            if line == None:
                return
            start = get_new_solution(s, r, line, (start[0], basis[0][0]))
        yield from solver_recurse(s, r, new_basis, depth - i, start)


# solver step using chord and tangent methods up to a maximum Manhattan distance of depth in Z^n, with n basis points
def solver(s: Fraction, r: Fraction, basis: Iterable[tuple[Fraction, Fraction]], depth: int):
    depth += 1
    triples: set[tuple[Fraction, Fraction, Fraction]] = set()
    for x, y in solver_recurse(s, r, list(basis), depth, None):
        z = (x + y) / (x * y / r - 1)
        if x <= 0 or y <= 0 or z <= 0:
            continue
        variables = sorted([x, y, z])
        triples.add((variables[0], variables[1], variables[2]))
    return triples


# reducer step using solver function with a maximum Manhattan distance of depth
def reducer(s: Fraction, r: Fraction, inputs: Iterable[tuple[Fraction, Fraction, Fraction]], depth: int):
    sorted_inputs = sorted(inputs, key=height)

    superset = sorted_inputs
    for i in range(1, len(sorted_inputs)):
        result = solver(s, r, [(x[0], x[2]) for x in sorted_inputs[:i]], depth)
        if all([x in result for x in inputs]):
            superset = sorted_inputs[:i]
            break

    flag = True
    while flag:
        for i in range(len(superset) - 1, -1, -1):
            generators = superset.copy()
            del generators[i]
            result = solver(s, r, [(x[0], x[2]) for x in generators], depth)
            if all([x in result for x in inputs]):
                superset = generators
                flag = True
                break
        flag = False

    return superset


# convert elliptic curve solutions to triangles, and print them
def print_solutions(solutions: Iterable[tuple[Fraction, Fraction, Fraction]], command_max: int = -1):
    triangles: list[tuple[Fraction, Fraction, Fraction]] = []
    for x, y, z in solutions:
        triangles.append((x + y, x + z, y + z))
    triangles.sort(key=height)
    strings = [[str(y) for y in x] for x in triangles]
    for point in strings:
        print(", ".join(point))
    print(f"count: {len(strings)}")
    if len(strings) <= command_max:
        print("solver.py " + " ".join([",".join(x) for x in strings]))


def main() -> int:
    start = time.perf_counter_ns()
    previous = start
    s = None
    r = None
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
        [_s, _r, variables] = variables
        inputs.add(variables)
        if (s != None and s != _s) or (r != None and r != _r):
            print(f"triangle '{sys.argv[i]}' has mismatched properties, exiting")
            return 1
        s = _s
        r = _r
        i += 1

    j = 0
    while i + j < len(sys.argv):
        if sys.argv[i + j] == "-v":
            verbose_printing = True
        elif j < 3:
            depths[j] = int(sys.argv[i + j])
        j += 1
    [searcher_depth, reducer_depth, solver_depth] = depths

    if s == None or r == None:
        print("usage: solver.py <a,b,c [d,e,f] ...> [search depth = 12] [reducer depth = 10] [solver depth = 8] [-v]")
        return 1

    inputs = inputs.union(searcher(s, r, searcher_depth))
    if verbose_printing:
        print_solutions(inputs, 32)
        elapsed = time.perf_counter_ns() - previous
        previous += elapsed
        print(f"searcher finished in {elapsed // 1000000} ms\n")

    generators = reducer(s, r, inputs, reducer_depth)
    if verbose_printing:
        print_solutions(generators, 32)
        elapsed = time.perf_counter_ns() - previous
        previous += elapsed
        print(f"reducer finished in {elapsed // 1000000} ms\n")

    basis = [(x[0], x[2]) for x in generators]
    solutions = solver(s, r, basis, solver_depth)
    print_solutions(solutions)
    if verbose_printing:
        elapsed = time.perf_counter_ns() - previous
        previous += elapsed
        print(f"solver finished in {elapsed // 1000000} ms")
        elapsed = previous - start
        print(f"script finished in {elapsed // 1000000} ms")

    return 0


if __name__ == "__main__":
    exit(main())
