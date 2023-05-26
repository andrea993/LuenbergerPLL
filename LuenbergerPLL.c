#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


struct SS
{
    double A[2][2];
    double B[2];
    double C[2][2];
    double D[2];
};


struct SS FilterDesign(double f0, double BW, double dt)
{
    struct SS ss;
    double w0 = 2*M_PI*f0;
    double Q = f0/BW;
    double A[2][2] = {{ 0.0, -w0 }, { w0, -w0/Q }};
    double B[2] = { 0, w0/Q };
    double a = w0/tan(dt*w0/2.0);
    double bcg = sqrt(2.0/a);

    double id[2][2] = { { 1.0,            -A[0][1]/a     }, 
                        { -A[1][0]/a,   -A[1][1]/a + 1.0 } };

    double d = id[0][0]*id[1][1] - id[0][1]*id[1][0];

    double I[2][2] = {  { id[1][1]/d,     -id[0][1]/d   }, 
                        { -id[1][0]/d,   id[0][0]/d     } };

    double P[2][2] =  { { 1.0,        A[0][1]/a         }, 
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

void Step(const struct SS *ss, double x[], double u, double y[], double x_next[])
{
    y[0] = ss->D[0]*u + ss->C[0][1]*x[1] + ss->C[0][0]*x[0];
    y[1] = ss->D[1]*u + ss->C[1][1]*x[1] + ss->C[1][0]*x[0];

    x_next[0] = ss->B[0]*u + ss->A[0][1]*x[1] + ss->A[0][0]*x[0];
    x_next[1] = ss->B[1]*u + ss->A[1][1]*x[1] + ss->A[1][0]*x[0];
}


/* Test on a signal */
int main()
{
    double dt = 0.01; /* sec sampling time */
    double f0 = 120; /* Hz */
    double x[2] = { 0, 1 }; /* filter state, please avoid (0,0) point */
    double x_next[2]; /* filter predicted state */
    double y[2]; /* filter output */

    struct SS ss = FilterDesign(f0, 8, dt);

    /* this example will generate a sin wave with additive noise and follow its characteristics */
    #define L 10240

    double M = 15; /* wave amplitude */
    double phi_0 = 2.0/3.0*M_PI; 
    double M_est[L]; /* estimated amplitude */
    double phi_est[L];
    
    double n;
    double u;
    int i;



    for (i = 0; i < L ; i++)
    {
        n = rand()/(double)RAND_MAX;
        u = M*cos(2*M_PI*f0*i*dt + phi_0) + n;

        Step(&ss, x, u, y, x_next);
        
        /* State update */
        x[0] = x_next[0];
        x[1] = x_next[1];

        M_est[i] = sqrt(y[0]*y[0] + y[1]*y[1]);
        /* M_est[i] = sqrt((y[0]*y[0] + y[1]*y[1])/2.0) for RMS value */
        phi_est[i] = atan2(y[1], y[0]);

        printf("%f ", M_est[i]);
    }

    

    
    return  EXIT_SUCCESS;
}

