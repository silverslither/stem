# vi: syntax=python

P_MAX = 100
S_MAX = 50

solutions: set[tuple[int, ...]] = set()

for m in range(1, P_MAX + 1):
    for n in range(1, P_MAX + 1):
        a = (m + n) * (2 * m^3 + m^2 * n + 2 * m * n^2 - n^3)
        b = (m - n) * (2 * m^3 - m^2 * n + 2 * m * n^2 + n^3)
        c = 2 * m^2 * (m^2 + 3 * n^2)
        if a <= 0 or b <= 0 or c <= 0:
            continue
        d = gcd(gcd(a, b), c)
        a //= d
        b //= d
        c //= d

        solutions.add(tuple(sorted([a, b, c])))

R.<x,y,z> = QQ[]
for sides in sorted(solutions, key=lambda x: sum(x))[:S_MAX]:
    print(", ".join([str(x) for x in sides]))

    [a, b, c] = sides
    s = (a + b + c) / 2
    AA = s * (s - a) * (s - b) * (s - c)

    E = EllipticCurve(AA / s * z^3 + (x + y - s * z) * x * y, [s - a, s - b, 1])
    print("rank =", E.rank(only_use_mwrank=False))
