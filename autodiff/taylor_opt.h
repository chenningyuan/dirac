#ifndef TAYLOR_OPT_H
#define TAYLOR_OPT_H
#include "taylor.h"


// Multiply the Hessian of the taylor expansion t by the vector x,
// put the result in y.
template<int Nvar>
static void hessvec(tnum_t y[Nvar], const taylor<Nvar, 2> &t, const tnum_t x[Nvar])
{
    int k = Nvar+1;
    for (int i=0;i<Nvar;i++)
    {
	y[i] = 0;
	for (int j=i;j<Nvar;j++)
	    y[i] += x[j]*t[k++];
    }
    k = Nvar + 1;
    for (int j=0;j<Nvar;j++)
	for (int i=j;i<Nvar;i++)
	    y[i] += x[j]*t[k++];
}

template<int N>
static tnum_t dot(const tnum_t x[], const tnum_t y[])
{
    tnum_t sum = 0;
    for (int i=0;i<N;i++)
	sum += x[i]*y[i];
    return sum;
}

// Find the minimum x of a second order Taylor expansion. The value of x on
// input will be used as a starting guess for the conjugate gradient method.
// The expansion must be positive definite to guarantee a solution.
template<int Nvar>
int taylor_opt(tnum_t x[Nvar], const taylor<Nvar, 2> &t, tnum_t thres)
{
    int n = 1;
    tnum_t r[Nvar], d[Nvar], q[Nvar];
    tnum_t delta_new, delta_old, stopping = -1, alpha;
 restart:
    hessvec(r,t,x);
    for (int i=0;i<Nvar;i++)
    {
	r[i] = -t[i+1] - r[i];
	d[i] = r[i];
    }
    delta_new = dot<Nvar>(r,r);
    if (stopping < 0)
	stopping = delta_new*thres*thres;
    while (n < 2*Nvar and delta_new > stopping)
    {
	hessvec(q,t,d);
	alpha = delta_new/dot<Nvar>(d,q);
	for (int i=0;i<Nvar;i++)
	    x[i] += alpha*d[i];
	if (n % 50 != 0)
	    for (int i=0;i<Nvar;i++)
		r[i] -= alpha*q[i];
	else
	    goto restart;
	delta_old = delta_new;
	delta_new = dot<Nvar>(r,r);
	cout << delta_new << endl;
	for (int i=0;i<Nvar;i++)
	    d[i] = r[i]+(delta_new/delta_old)*d[i];
	n++;
    }
    //cout << "delta new at end: " << delta_new << " n = " << n << endl;
    if (delta_new > stopping)
	return -1;
    else
	return 0;
}

#endif
