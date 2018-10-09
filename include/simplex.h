#ifndef __NAFFLIB_SIMPLEX_H__
#define __NAFFLIB_SIMPLEX_H__

#include <complex.h>
#include <frequency.h>

void simplex_swap(double* f1, double *f2, double _Complex *x1, double _Complex *x2);
void simplex_order(double *f1, double *f2, double *f3, double _Complex *x1, double _Complex *x2, double _Complex *x3);
double _Complex simplex_minimize(double (*minfunc)(double _Complex, const merit_args*), double _Complex x1, double _Complex x2, double _Complex x3, const merit_args* S );
double rosenbrock(double _Complex z, const merit_args* S);
void test_simplex();

#endif
