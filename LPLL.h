#ifndef LPLL_H

#define LPLL_H

#define _USE_MATH_DEFINES
#include <math.h>

struct LPLL_SS
{
    double A[2][2];
    double B[2];
    double C[2][2];
    double D[2];
};


struct LPLL_SS LPLL_filterDesign(double f0, double BW, double dt);
void LPLL_step(const struct LPLL_SS *ss, double x[], double u, double y[]);

#endif


