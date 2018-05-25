#include <math.h>
extern "C" {
double fitFun(double kT2, double x, const double *p) {
    return (exp(p[0]*log(kT2) + p[1]*log(kT2)*log(kT2)) * pow(1-x,p[2]));
}
}
