from fractions import Fraction
from math import isqrt, gcd
import sys

sys.set_int_max_str_digits(65536)

# returns isqrt, with optional check for exactness; if not exact, returns -1
def fraction_isqrt(x: Fraction, exact = False):
    y = Fraction(isqrt(x.numerator), isqrt(x.denominator))
    if exact and y * y != x:
        return Fraction(-1)
    return y

# generates Heronian triangles with integer sides.
# returns { C: set[(a, b, c)] }
def generate_triangles(MAX: int):
    map: dict[Fraction, set[tuple[int, int, int]]] = dict()
    for m in range(1, MAX + 1):
        for n in range(1, MAX + 1):
            for p in range(1, MAX + 1):
                for q in range(1, MAX + 1):
                    if m * p <= n * q:
                        continue
                    a = m * n * (p * p + q * q)
                    b = p * q * (m * m + n * n)
                    c = (m * q + n * p) * (m * p - n * q)
                    
                    s = Fraction(a + b + c, 2)
                    A = Fraction(m * n * p * q * c)
                    C = s * s / A

                    d = gcd(a, b, c)
                    a //= d
                    b //= d
                    c //= d

                    if a > b:
                        a, b = b, a
                    if b > c:
                        b, c = c, b
                    if a > b:
                        a, b = b, a

                    if not C in map:
                        map[C] = set()
                    map[C].add((a, b, c))
    return map

def main() -> int:
    if len(sys.argv) < 2:
        print("usage: lookup.py <a,b,c> [max] [threshold]")
        return 1
    
    i = 1
    s = None
    C = None
    if "," in sys.argv[1]:
        [a, b, c] = [Fraction(x) for x in sys.argv[1].split(",")]
        if a == 0 or b == 0 or c == 0 or a + b <= c or b + c <= a or c + a <= b:
            print(f"invalid Heronian triangle, exiting")
            return 1

        s = (a + b + c) / 2
        A = fraction_isqrt(s * (s - a) * (s - b) * (s - c), True)
        if A == -1:
            print(f"invalid Heronian triangle, exiting")
            return 1

        C = s * s / A
        i += 1

    MAX = 20
    if len(sys.argv) > i:
        MAX = int(sys.argv[i])
    THRESHOLD = None
    if len(sys.argv) > i + 1:
        THRESHOLD = int(sys.argv[i + 1])

    map = generate_triangles(MAX)
    if THRESHOLD != None:
        for key in map:
            if len(map[key]) >= THRESHOLD:
                print(key, map[key])

    if s == None or C == None:
        return 0

    if not C in map:
        print("count: 0")
        return 0

    final: list[list[Fraction]] = []
    for sides in map[C]:
        _s = Fraction(sides[0] + sides[1] + sides[2], 2)
        ratio = s / _s
        final.append([ratio * x for x in sides])

    final.sort(key = lambda x: x[0])
    for scaled in final:
        print(", ".join([f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator) for x in scaled]))
    print(f"count: {len(final)}")

    return 0

if __name__ == "__main__":
    exit(main())
