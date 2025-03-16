from mpmath import mp
mp.dps = 64

D = [mp.mpf(1)/8, mp.mpf(3)/4, mp.mpf(1)/8] # bspline2
#D = [mp.mpf(1)/6, mp.mpf(2)/3, mp.mpf(1)/6] # bspline3
#D = [mp.mpf(3)/21, mp.mpf(13)/21, mp.mpf(4)/21] # omoms3
d = 20

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
    return f"{f"({start})" : <5}{str(arr[start:end]) : <40}({n - end})"


D += [0] * (d - len(D))

C = 1 / (D[0] + D[1])

M = [[D[1] * C, D[0] * C] + [0] * (d - 2)]
for i in range(d - 2):
    M.append(D[-i:] + D[:-i])
M.append([0] * (d - 2) + [D[0] * C, D[1] * C])

L, U = lu_decomp(M)
N = C / U[0][1];

print("N", float(N))
print("C", float(C), "\n")
for v in L:
    print("L", pretty_zpadded([float(x) for x in v]))
print("")
for v in U:
    print("U", pretty_zpadded([float(N * x) for x in v]))
print("\nL CONSTANTS")
for i in range(1, len(L)):
    print(float(L[i][i - 1]))
