import subprocess
from itertools import combinations
from math import gcd

p = 0
s: set[tuple[int, ...]] = set()
for t in combinations(list(range(100)), 3):
    if t[0] + t[1] <= t[2] or t[1] + t[2] <= t[0] or t[0] + t[2] <= t[1]:
        continue
    if t[0] == t[1] or t[1] == t[2] or t[2] == t[0]:
        continue
    d = gcd(*t)
    t = tuple(v // d for v in sorted(t))
    if t in s:
        continue
    s.add(t)

    if t[0] != p:
        print("f", t[0])
        p = t[0]

    if t[0] < 62:
        continue

    t = ",".join([str(x) for x in t])
    r = subprocess.run(["python", "solver.py", t, "0", "0", "2"], capture_output=True)
    l = r.stdout.splitlines()
    if len(l) >= 3:
        continue
    print(l)
