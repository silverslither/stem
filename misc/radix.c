#include <stdint.h>
#include <stdlib.h>

typedef int32_t s32;
typedef uint32_t u32;

void rsort(u32 *arr, size_t len) {
    size_t i;
    u32 *aux = malloc(len << 2);

    size_t buf[0x1400] = { 0 };
    size_t *c0 = buf;
    size_t *c1 = buf + 0x800;
    size_t *c2 = buf + 0x1000;

    for (i = 0; i < len; i++) {
        c0[arr[i] & 0x7ff]++;
        c1[(arr[i] >> 11) & 0x7ff]++;
        c2[arr[i] >> 22]++;
    }

    for (i = 1; i < 0x800; i++) {
        c0[i] += c0[i - 1];
        c1[i] += c1[i - 1];
    }

    for (i = 1; i < 0x400; i++)
        c2[i] += c2[i - 1];

    for (i = len; i-- > 0;)
        aux[--c0[arr[i] & 0x7ff]] = arr[i];

    for (i = len; i-- > 0;)
        arr[--c1[(aux[i] >> 11) & 0x7ff]] = aux[i];

    for (i = len; i-- > 0;)
        aux[--c2[arr[i] >> 22]] = arr[i];

    for (i = 0; i < len; i++)
        arr[i] = aux[i];

    free(aux);
}

static inline u32 s2u_u2s(s32 v) {
    return v + 0x80000000;
}

void rsorts(s32 *arr, size_t len) {
    size_t i;
    u32 *u_arr = (u32 *)arr;
    for (i = 0; i < len; i++)
        u_arr[i] = s2u_u2s(arr[i]);
    rsort(u_arr, len);
    for (i = 0; i < len; i++)
        arr[i] = s2u_u2s(u_arr[i]);
}

static inline u32 f2u(float v) {
    u32 r = *(u32 *)&v;
    u32 mask = ((s32)r >> 31) | 0x80000000;
    return r ^ mask;
}

static inline float u2f(u32 r) {
    u32 mask = ((s32)(~r) >> 31) | 0x80000000;
    u32 v = r ^ mask;
    return *(float *)&v;
}

void rsortf(float *arr, size_t len) {
    size_t i;
    u32 *u_arr = (u32 *)arr;
    for (i = 0; i < len; i++)
        u_arr[i] = f2u(arr[i]);
    rsort(u_arr, len);
    for (i = 0; i < len; i++)
        arr[i] = u2f(u_arr[i]);
}
