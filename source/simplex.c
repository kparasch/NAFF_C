#include "simplex.h"

void simplex_swap(double* f1, double *f2, double _Complex *x1, double _Complex *x2)
{
    double temp = *f1;
    *f1 = *f2;
    *f2 = temp;

    double _Complex tempc = *x1;
    *x1 = *x2;
    *x2 = tempc;
    
    return;
}

void simplex_order(double *f1, double *f2, double *f3, double _Complex *x1, double _Complex *x2, double _Complex *x3)
{
    if( *f1 > *f3 )
        simplex_swap(f1, f3, x1, x3);
    if( *f1 > *f2 )
        simplex_swap(f1, f2, x1, x2);
    if( *f2 > *f3 )
        simplex_swap(f2, f3, x2, x3);
    return;
}

double _Complex simplex_minimize(double (*minfunc)(double _Complex, const merit_args*), double _Complex x1, double _Complex x2, double _Complex x3, const merit_args* S )
{
    //const int max_iter = 10000;
    const int max_iter = 100;
    const double tolerance = 1.490116e-8;
    //const double tolerance = 1.490116e-15;

    const double alpha = 1.;
    const double gamma = 2.;
    const double rho = 0.5;
    const double sigma = 0.5;
    
    double f1 = minfunc(x1,S);
    double f2 = minfunc(x2,S);
    double f3 = minfunc(x3,S);

    for(int i = max_iter; i--;)
    {
        //1 order
        simplex_order(&f1, &f2, &f3, &x1, &x2, &x3);
        printf("%lf, %lf\n",creal(x1), cimag(x1));
        if( f3 - f1 < tolerance )
            return x1;
    
        //2 centroid
        double _Complex x0 = 0.5*(x1+x2);
    
        //3 reflection
        double _Complex xr = x0 + alpha*(x0-x3);
        double fr = minfunc(xr,S);
        if( f1 <= fr && fr < f3 )
        {
            x3 = xr;
            f3 = fr;
            continue;
        }

        //4 expansion
        if( fr < f1 )
        {
            double _Complex xe = x0 + gamma*(xr-x0);
            double fe = minfunc(xe,S);
            if( fe < fr )
            {
                x3 = xe;
                f3 = fe;
                continue;
            }
            else
            {
                x3 = xr;
                f3 = fr;
                continue;
            }
        }

        //5 contraction
        double _Complex xc = x0 + rho*(x3 - x0);
        double fc = minfunc(xc,S);
        if( fc < f3 )
        {
            x3 = xc;
            f3 = fc;
            continue;
        }

        //6 shrink
        x2 = x1 + sigma*(x2 - x1);
        x3 = x1 + sigma*(x3 - x1);
    }

    printf("***WARNING***: Simplex exceeded maximum number of calls.\n");
    printf("***WARNING***: Simplex convergence is not guaranteed.\n");
    return x1;
}

double rosenbrock(double _Complex z, const merit_args* S)
{
    double x = creal(z);
    double y = cimag(z);
    const double a = 1.;
    const double b = 100.;

    return (a-x)*(a-x) + b*(y-x*x)*(y-x*x);
}

void test_simplex()
{
    merit_args* S = NULL;
    double (*minfunc)(double _Complex, const merit_args*) = rosenbrock;
    double _Complex x1 = -2.1 -I*2.1;
    double _Complex x2 = -2.0 -I*2.0;
    double _Complex x3 = -2.1 -I*2.0;
    double _Complex minres = simplex_minimize(minfunc, x1, x2, x3, S);
    printf("Rosenbrock Simplex Minimization result:  %1.15lf + i %1.15lf\n",creal(minres),cimag(minres));
    printf("          True minimum lies at        :  %1.15lf + i %1.15lf\n",1.,1.);
    return;

}
