from mpmath import mp

mp.dps = 128

# D = [mp.mpf(1) / 8, mp.mpf(3) / 4, mp.mpf(1) / 8]  # bspline2
# D = [mp.mpf(1) / 6, mp.mpf(2) / 3, mp.mpf(1) / 6]  # bspline3
# D = [mp.mpf(4) / 21, mp.mpf(13) / 21, mp.mpf(4) / 21]  # omoms3
# D = [mp.mpf(346) / 675675, mp.mpf(6101) / 200200, mp.mpf(1202) / 5005, mp.mpf(247409) / 540540, mp.mpf(1202) / 5005, mp.mpf(6101) / 200200, mp.mpf(346) / 675675] # omoms7
D = [mp.mpf(302399) / 1556845012800, mp.mpf(673972423) / 6227380051200, mp.mpf(2711965633) / 518948337600, mp.mpf(31678817987) / 518948337600, mp.mpf(3006012673) / 12355912800, mp.mpf(56442621569) / 148270953600, mp.mpf(3006012673) / 12355912800, mp.mpf(31678817987) / 518948337600, mp.mpf(2711965633) / 518948337600, mp.mpf(673972423) / 6227380051200, mp.mpf(302399) / 1556845012800] # omoms11
d = 128

def lu_decomp(mat):
    n = len(mat)
    L = [[mp.mpf(0) for _ in range(n)] for _ in range(n)]
    U = [[mp.mpf(0) for _ in range(n)] for _ in range(n)]

    for i in range(n):
        for j in range(i, n):
            U[i][j] = mp.mpf(mat[i][j]) - sum(L[i][k] * U[k][j] for k in range(i))

        for j in range(i, n):
            if i == j:
                L[i][i] = mp.mpf(1)
            else:
                L[j][i] = (mp.mpf(mat[j][i]) - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

    return L, U


def pretty_zpadded(arr):
    n = len(arr)
    start = 0
    end = n
    for i in range(n):
        if arr[i] == 0:
            start = i + 1
        else:
            break
    for i in range(n - 1, -1, -1):
        if arr[i] == 0:
            end = i
        else:
            break
    if end < start:
        return f"({n})"
    return f'{f"({start})" : <5}{str(arr[start:end])} ({n - end})'

M = []
H = len(D) // 2
for i in range(d):
    o = i - H
    n = sum(D)
    m = []
    for j in range(-H, H + 1):
        if i + j < 0:
            o += 1
            n -= D[j + H]
        elif i + j >= d:
            n -= D[j + H]
        else:
            m.append(D[j + H])
    m = [x / n for x in m]
    m = [0] * o + m
    m = m + [0] * (d - len(m))
    M.append(m)

L, U = lu_decomp(M)
N = 1 / U[-H - 1][-1]

for v in M:
    print("M", pretty_zpadded([float(x) for x in v]))
print("\nN", float(N), "\n")
for v in L:
    print("L", pretty_zpadded([float(x) for x in v]))
print("")
for v in U:
    print("U", pretty_zpadded([float(N * x) for x in v]))
C = []
for _ in range(1, H + 1):
    C.append([])
    print("\nL CONSTANTS")
    for i in range(_, len(L)):
        print(float(L[i][i - _]), end=",\n")
        if len(L) - i <= H:
            C[-1].append(str(float(L[i][i - _] / L[i - H][i - _ - H])))

print("\nMULT CONSTANTS")
print(",\n".join([y for x in C for y in reversed(x)]))
