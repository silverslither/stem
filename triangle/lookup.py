from itertools import product, combinations
from fractions import Fraction
from math import gcd
import time
import sys

sys.set_int_max_str_digits(65536)

# generate ordered coprime pairs using the tree (a, b) => (a + b, b), (a, a + b), starting at (1, 2)
def generate_pairs(max_depth: int) -> list[tuple[int, int]]:
    pairs = [(1, 1)]
    if max_depth <= 0:
        return pairs
    stack = [(1, 2, 0)]
    while len(stack):
        (a, b, depth) = stack.pop()
        if a > b:
            pairs.append((b, a))
        else:
            pairs.append((a, b))
        depth += 1
        if depth == max_depth:
            continue
        stack.append((a + b, b, depth))
        stack.append((a, a + b, depth))
    return pairs

# return sorted tuple given a <= b
def sift(a, b, c):
    if c >= b:
        return (a, b, c)
    if c >= a:
        return (a, c, b)
    return (c, a, b)

# generate primitive triangles using the generator (a, b), (c, d) => (ac, bc, d), (ac, ad, b), (c, ad, bd), (a, bc, bd), from list of ordered coprime pairs
def generate_triangles(pairs: list[tuple[int, int]]):
    triangles: set[tuple[int, int, int]] = set()
    for ((a, b), (c, d)) in combinations(pairs, 2):
        ad = a * d
        bc = b * c
        ac = a * c
        bd = b * d
        if bd < c + ad:
            triangles.add((c, ad, bd))
        if bd < a + bc:
            triangles.add((a, bc, bd))
        _ = sift(ac, bc, d)
        if _[2] < _[0] + _[1]:
            triangles.add(_)
        _ = sift(ac, ad, b)
        if _[2] < _[0] + _[1]:
            triangles.add(_)

    map: dict[Fraction, list[tuple[int, int, int]]] = {}
    for triangle in triangles:
        (a, b, c) = triangle
        p = a + b + c
        C = Fraction(p ** 3, (p - 2 * a) * (p - 2 * b) * (p - 2 * c))
        if not C in map:
            map[C] = []
        map[C].append(triangle)

    return map

def from_params(w: int, x: int, y: int, z: int):
    if w * y <= x * z:
        return None
    a = w * x * (y * y + z * z)
    b = y * z * (w * w + x * x)
    c = (w * z + x * y) * (w * y - x * z)

    A = w * x * y * z * c
    p = a + b + c
    C = Fraction(p * p, A)

    d = gcd(a, b, c)
    a //= d
    b //= d
    c //= d
    [a, b, c] = sorted([a, b, c])
    
    return (C, (a, b, c))


# generates Heronian triangles with integer sides.
def generate_heronian_triangles(pairs: list[tuple[int, int]]):
    map: dict[Fraction, set[tuple[int, int, int]]] = dict()

    for ((x, w), (z, y)) in product(pairs, pairs):
        kvpair = from_params(w, x, y, z)
        if kvpair == None:
            continue
        (C, triangle) = kvpair
        if not C in map:
            map[C] = set()
        map[C].add(triangle)

    for ((w, x), (z, y)) in product(pairs, pairs):
        kvpair = from_params(w, x, y, z)
        if kvpair == None:
            continue
        (C, triangle) = kvpair
        if not C in map:
            map[C] = set()
        map[C].add(triangle)

    return map

def main() -> int:
    start = time.perf_counter_ns()
    if len(sys.argv) < 3:
        print("usage: lookup.py <depth> <threshold> [-h]")
        return 1

    depth = int(sys.argv[1])
    threshold = int(sys.argv[2])
    heronian = False
    if len(sys.argv) > 3 and sys.argv[3] == "-h":
        heronian = True

    map = None
    if heronian:
        map = generate_heronian_triangles(generate_pairs(depth))
    else:
        map = generate_triangles(generate_pairs(depth))

    for key in sorted(map.keys(), key=lambda x: x.numerator * x.denominator):
        if len(map[key]) >= threshold:
            print(f"{key}:", end="")
            for triangle in sorted(map[key], key=lambda x: (sum(x), x[0])):
                print(" " + ",".join([f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator) for x in triangle]), end="")
            print("")
    elapsed = time.perf_counter_ns() - start
    print(f"script finished in {elapsed // 1000000} ms")

    return 0

if __name__ == "__main__":
    exit(main())
