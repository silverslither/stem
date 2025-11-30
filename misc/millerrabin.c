#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define assert(cond, ...)             \
    if (!(cond)) {                    \
        fprintf(stderr, __VA_ARGS__); \
        exit(1);                      \
    }                                 \
    (void)0

typedef uint32_t u32;
typedef uint64_t u64;

// compute a^b (mod m)
u64 expmod(u64 a, u64 b, u64 m) {
    u64 r = 1;
    for (int i = 0;; i++) {
        if (b & 1)
            r = (r * a) % m;
        if (i >= 31)
            break;
        b >>= 1;
        a = (a * a) % m;
    }
    return r;
}

// checks if w^n (mod m) = 1 or w^(2^i * n) (mod m) = -1 for 0 <= i < a
// returns 1 for a probable prime, 0 for a composite
int miller_rabin_witness(u64 w, u64 n, int a, u64 m) {
    if (w == m) // k would always be zero, so return probable prime
        return 1;
    u64 k = expmod(w, n, m);
    if (k == 1 || k == m - 1)
        return 1;
    for (int i = 1; i < a; i++) {
        k = (k * k) % m;
        if (k == m - 1)
            return 1;
    }
    return 0;
}

// driver logic for miller-rabin 
// returns 1 for a prime, 0 for a composite
int miller_rabin(u32 n) {
    if (n == 2)
        return 1;
    if ((n & 1) == 0 || n == 1)
        return 0;

    u32 m = n;
    n--;

    int a = 0;
    while ((n & 1) == 0) {
        n >>= 1;
        a++;
    }

    return miller_rabin_witness(2, n, a, m) &&
           miller_rabin_witness(7, n, a, m) &&
           miller_rabin_witness(61, n, a, m);
}

int main() {
    u32 n;
    assert(scanf("%u", &n) == 1, "invalid input");
    if (miller_rabin(n))
        printf("prime\n");
    else
        printf("not prime\n");
    return 0;
}
