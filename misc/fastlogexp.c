#include <math.h>
#include <stdint.h>

double q_exp2(double x) {
    if (x != x)
        return x;
    if (x < -1022.0)
        return 0.0;
    if (x >= 1024.0)
        return INFINITY;
    const double n = __builtin_roundeven(x);
    const double z = x - n;
    double v = 7.074194542204488e-09 + 4.457533149527563e-10 * z;
    v = 1.0178045295522382e-07 + v * z;
    v = 0.000001321543253438167 + v * z;
    v = 0.000015252733871037077 + v * z;
    v = 0.00015403530463727703 + v * z;
    v = 0.0013333558146374227 + v * z;
    v = 0.009618129107587253 + v * z;
    v = 0.05550410866482177 + v * z;
    v = 0.24022650695910158 + v * z;
    v = 0.6931471805599453 + v * z;
    v = 1.0 + v * z;
    const int64_t v_ = *(int64_t *)&v + ((int64_t)n << 52);
    return *(const double *)&v_;
}

/*
#include <stdio.h>

const double INC = 1.0 / ((uint64_t)1 << 0);
int main() {
    double maxerr = 0;
    for (double i = -1022; i <= 1025; i += INC) {
        double exact = exp2(i);
        double approx = q_exp2(i);
        if (exact != 0.0) {
            double err = fabs((double)(*(int64_t*)&approx - *(int64_t*)&exact));
            maxerr = err > maxerr ? err : maxerr;
        }
        printf("%-24.16lg %-24.16lg %-24.16lg\n", i, exact, approx);
    }
    printf("MAXERR %g ulp\n", maxerr);
    return 0;
}
/*/
int main() {
    volatile double x;
    for (double i = 0; i < 16; i += 1.0 / (1 << 24)) {
        x = exp2(i);
    }
}
//*/
