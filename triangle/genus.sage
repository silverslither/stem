# vi: syntax=python

R.<x, y> = QQ[]

T(p, q) = (
    (q * (p ** 2 * q - 2 * p - q) ** 2) / ((p - q) * (p * q - 1) * (p * q ** 2 - p - 2 * q)),
    (p * (p * q ** 2 - p - 2 * q) ** 2) / ((q - p) * (p * q - 1) * (p ** 2 * q - 2 * p - q))
)

def simplify(c):
    c = R(c.full_simplify().numerator())
    c = list(c.factor())[-1][0]
    return c

x2, y2 = T(x, y)
c = simplify(x2 - y2)

print(c)
print("deg =", c.degree())
print("g =", Curve(c).genus())
