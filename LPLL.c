#include "LPLL.h"

__attribute__((const)) struct LPLL_SS LPLL_filterDesign(double f0, double BW, double dt)
{
    struct LPLL_SS ss;
    const double w0 = 2*M_PI*f0;
    const double Q = f0/BW;
    const double A[2][2] = { { 0.0, -w0 }, 
                             { w0, -w0/Q } };
    const double B[2] = { 0, w0/Q };
    const double a = w0/tan(dt*w0/2.0);
    const double bcg = sqrt(2.0/a);

    const double id[2][2] = {   { 1.0,            -A[0][1]/a     }, 
                                { -A[1][0]/a,   -A[1][1]/a + 1.0 } };

    const double d = id[0][0]*id[1][1] - id[0][1]*id[1][0];

    const double I[2][2] = { { id[1][1]/d,     -id[0][1]/d   }, 
                             { -id[1][0]/d,   id[0][0]/d     } };

    const double P[2][2] =  { { 1.0,        A[0][1]/a         }, 
                            { A[1][0]/a,  A[1][1]/a + 1.0   } };


    ss.A[0][0] = I[1][0]*P[0][1] + I[0][0]*P[0][0];
    ss.A[0][1] = I[1][1]*P[0][1] + I[0][1]*P[0][0];
    ss.A[1][0] = I[1][0]*P[1][1] + I[0][0]*P[1][0];
    ss.A[1][1] = I[1][1]*P[1][1] + I[1][0]*P[0][1];

    ss.B[0] = ( B[1]*I[0][1] + B[0]*I[0][0] )*bcg;
    ss.B[1] = ( B[1]*I[1][1] + B[0]*I[1][0] )*bcg;

    ss.C[0][0] = I[0][0]*bcg;
    ss.C[0][1] = I[0][1]*bcg;
    ss.C[1][0] = I[1][0]*bcg;
    ss.C[1][1] = I[1][1]*bcg;

    ss.D[0] = ( B[1]*I[0][1] + B[0]*I[0][0] )/a;
    ss.D[1] = ( B[1]*I[1][1] + B[0]*I[1][0] )/a;

    return ss; 
}

void LPLL_step(const struct LPLL_SS *ss, double x[], double u, double y[])
{
    double x_next[2];

    y[0] = ss->D[0]*u + ss->C[0][1]*x[1] + ss->C[0][0]*x[0];
    y[1] = ss->D[1]*u + ss->C[1][1]*x[1] + ss->C[1][0]*x[0];

    x_next[0] = ss->B[0]*u + ss->A[0][1]*x[1] + ss->A[0][0]*x[0];
    x_next[1] = ss->B[1]*u + ss->A[1][1]*x[1] + ss->A[1][0]*x[0];

    x[0] = x_next[0];
    x[1] = x_next[1];
}