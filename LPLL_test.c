#include <stdio.h>
#include <stdlib.h>

#include "LPLL.h"


int main()
{
    double dt = 1e-4; /* sec sampling time */
    double f0 = 120; /* Hz */
    double x[2] = { 0, 1 }; /* filter state, please avoid (0,0) point */
    double x_next[2]; /* filter predicted state */
    double y[2]; /* filter output */

    struct LPLL_SS ss = LPLL_filterDesign(f0, 8, dt);

    /* this example will generate a sin wave with additive noise and follow its characteristics */
    #define L 10240

    double M = 15; /* wave amplitude */
    double phi_0 = -2.0/3.0*M_PI; 
    double M_est; /* estimated amplitude */
    double phi_est; /* estimated phase */
    
    double n;
    double u;
    int i;



    for (i = 0; i < L ; i++)
    {
        n = (rand()/(double)RAND_MAX - 0.5)*(M/2);
        u = M*cos(2*M_PI*f0*i*dt + phi_0) + n;

        LPLL_step(&ss, x, u, y);

        M_est = sqrt(y[0]*y[0] + y[1]*y[1]);
        /* M_est[i] = sqrt((y[0]*y[0] + y[1]*y[1])/2.0) for RMS value */
        phi_est = atan2(y[1], y[0]);

        printf("%f %f %f %f\n", u, y[1], M_est, phi_est);
    }

    

    
    return  EXIT_SUCCESS;
}

