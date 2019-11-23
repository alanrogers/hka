#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eprintf.h"
#include "hka.h"

/** Constructor */
HKA *
HKA_construct(unsigned nLoci, const unsigned * nA, const unsigned * nB,
	      const unsigned * nsA, const unsigned * nsB,
	      const unsigned * nsAB, const double * SA, const double * SB,
	      const double * D) {
    unsigned bytes;
    char * base;

    HKA * h;

    h = (HKA *) malloc(sizeof(HKA));
    if(h == NULL)
	return NULL;

    h->nLoci = nLoci;

    /* 6 arrays of unsigned and 9 of double */
    bytes = nLoci*(6*sizeof(unsigned) + 9*sizeof(double));

    /* All arrays are allocated w/ a single call to malloc */
    base = (char *) malloc( bytes );
    if(base == NULL) {
	free(h);
	return NULL;
    }

    h->nsA = (unsigned *) base;          /* on destruction, free nsA */
    base += nLoci * sizeof(unsigned);
    
    h->nsA = (unsigned *) base;
    base += nLoci * sizeof(unsigned);

    h->nsB = (unsigned *) base;
    base += nLoci * sizeof(unsigned);

    h->nsAB = (unsigned *) base;
    base += nLoci * sizeof(unsigned);

    h->nA = (unsigned *) base;
    base += nLoci * sizeof(unsigned);

    h->nB = (unsigned *) base;
    base += nLoci * sizeof(unsigned);

    h->SA = (double *) base;
    base += nLoci * sizeof(double);

    h->SB = (double *) base;
    base += nLoci * sizeof(double);

    h->D = (double *) base;
    base += nLoci * sizeof(double);

    h->aA = (double *) base;
    base += nLoci * sizeof(double);

    h->aB = (double *) base;
    base += nLoci * sizeof(double);

    h->bA = (double *) base;
    base += nLoci * sizeof(double);

    h->bB = (double *) base;
    base += nLoci * sizeof(double);

    h->theta = (double *) base;
    base += nLoci * sizeof(double);

    h->tmp = (double *) base;

    HKA_reset(nA, nB, nsA, nsB, nsAB, SA, SB, D, h);

    return h;
}

/** Free memory allocated by h */
void
HKA_destroy(HKA * h) {
    free(h->nsA);
    free(h);
}

/** Put data into the arrays and estimate. */
void
HKA_reset(const unsigned * nA, const unsigned * nB,
      const unsigned * nsA, const unsigned * nsB, const unsigned * nsAB,
      const double * SA, const double * SB, const double * D,
      HKA * h) {
    unsigned i;
    
    memcpy(h->nA, nA, h->nLoci * sizeof(unsigned));
    memcpy(h->nB, nB, h->nLoci * sizeof(unsigned));
    memcpy(h->nsA, nsA, h->nLoci * sizeof(unsigned));
    memcpy(h->nsB, nsB, h->nLoci * sizeof(unsigned));
    memcpy(h->nsAB, nsAB, h->nLoci * sizeof(unsigned));
    memcpy(h->SA, SA, h->nLoci * sizeof(double));
    memcpy(h->SB, SB, h->nLoci * sizeof(double));
    memcpy(h->D, D, h->nLoci * sizeof(double));
    memset(h->theta, 0, h->nLoci * sizeof(double));

    for(i=0; i < h->nLoci; ++i) {
	sumRecip(&(h->aA[i]), &(h->bA[i]), h->nA[i]);
	sumRecip(&(h->aB[i]), &(h->bB[i]), h->nB[i]);
    }
    HKA_estimate(h);
}

/*
 * Calculate the X^2 goodness-of-fit statistic.  (See p 154 of
 * Hudson, Kreitman, and Aguade 1987.
 */
double
HKA_X2(HKA * h) {
    double x=0;
    unsigned i;

    for(i=0; i < h->nLoci; ++i) {
	double k = h->theta[i] * h->nsA[i];
	double m = k * h->aA[i];               /* mean */
	double v = m + k * k * h->bA[i];       /* var  */
	double diff = h->SA[i] - m;            /* obs - mean */

	x += diff*diff/v;                        /* X^2 for SA */

	if(h->nB[i] > 1) {
	    k = h->theta[i] * h->nsB[i];
	    m = k * h->aB[i];               /* mean */
	    v = m + k*k* h->bB[i];          /* var  */
	    diff = h->SB[i] - m;            /* obs - mean */
	    x += diff*diff/v;                 /* X^2 for SB */

	    k = h->nsAB[i] * h->theta[i];
	    m = k*(h->T + (1.0 + h->f)/2.0);
	    k *= (1.0 + h->f)/2.0;
	    v = m + k*k;
	}else {
	    k = h->nsAB[i] * h->theta[i];
	    m = k*(h->T + 1.0);
	    v = m + k*k;
	}
	diff = h->D[i] - m;
	x += diff*diff/v;                     /* X^2 for D */
    }
    return x;
}

/*
 * Estimate parameters theta_i, f, and T.  The function uses the
 * EM algorithm to solve equations 5, on p 154 of Hudson,
 * Kreitman, and Aguade 1987. 
 *
 * The system of equations must be re-written to account for variation in
 * numbers of sites.  There are several ways to do this.  This function
 * used the EM algorithm to solve the following system (which is modeled on
 * eqn 6 of Hudson, Kreitman and Aguade):
 *
 * nLoci       nLoci
 * Sum SA[i] = Sum theta[i]*nsA[i]*aA[i]
 * i=1         i=1
 *
 * nLoci           nLoci
 * Sum SB[i] = f * Sum theta[i] * nsB[i] * aB[i]
 * i=1              i=1
 *
 * nLoci                      nLoci
 * Sum D[i] = (T + (1+f)/2) * Sum theta[i] * nsAB[i]
 * i=1                        i=1
 *
 * SA[i] + SB[i] + D[i] 
 *         = theta[i]*(aA[i]*nsA[i] + f*aB[i]*nsB[i]
 *           + (T + (1+f)/2))*nsAB[i]
 *                                                       i=1,2,...,nLoci-1.
 *
 * Loci at which nB=1 do not contribute to the SB sum.
 */
void
HKA_estimate(HKA * h) {
    double hold_T, hold_f;
    double tol = 1e-6, delta, x, y, sumA, sumB, sumD;
    unsigned i;

    sumA = sumB = sumD = x = 0.0;
    for(i=0; i < h->nLoci; ++i) {
	sumA += h->SA[i];
	sumB += h->SB[i];
	sumD += h->D[i];
	x += h->nsA[i] * h->aA[i];
    }

    /* initialize */
    hold_T = hold_f = -9999.9;
    x = sumA/x;
    for(i=0; i < h->nLoci; ++i)
	h->theta[i] = x;

    /* loop until EM algorithm converges */
    do{
	hold_T = h->T;
	hold_f = h->f;
	memcpy(h->tmp, h->theta, (h->nLoci) * sizeof(double));

	h->use_f = 0;

	// Find f, T
	h->f = x = y = 0.0;
	for(i=0; i < h->nLoci; ++i) {
	    if(h->nB[i] > 1) {
		h->use_f = 1;
		x += h->theta[i] * h->nsB[i] * h->aB[i];
	    }
	    y += h->theta[i] * h->nsAB[i];
	}
	if(h->use_f) {
	    h->f = sumB / x;
	    h->T = sumD/y - (1.0 + h->f)/2.0;
	}else
	    h->T = sumD/y - 1.0;

	/* Find theta[i] */
	for(i=0; i < h->nLoci; ++i) {
	    x = h->SA[i] + h->SB[i] + h->D[i];
	    if(h->use_f) {
		y = h->aA[i] * h->nsA[i]
		    + h->f * h->aB[i] * h->nsB[i]
		    + (h->T + (1.0 + h->f)/2.0) * h->nsAB[i];
	    }else {
		y = h->aA[i] * h->nsA[i]
		    + (h->T + 1.0) * h->nsAB[i];
	    }
	    h->theta[i] = x/y;
	}

	/* delta measures change in this iteration */
	delta = fabs(h->T - hold_T) + fabs(h->f - hold_f);
	for(i=0; i < h->nLoci; ++i)
	    delta += fabs(h->theta[i] - h->tmp[i]);

	/* stop when delta is small */
    }while(delta > tol);
}

/* Report the value of parameter T */
double
HKA_T(HKA * h) { return h->T;}

/* Report the value of parameter f */
double
HKA_f(HKA * h) { return h->f; }

/* Report the value of parameter theta[i] */
double
HKA_theta(unsigned i, HKA * h) { return h->theta[i]; }

/* Report degrees of freedom */
unsigned
HKA_degreesOfFreedom(HKA * h) {
    unsigned i, df=0;

    /* count equations */
    for(i=0; i < h->nLoci; ++i) {
	if(h->nB[i] > 1)
	    df += 3;
	else
	    df += 2;
    }

    /* subtract parameters */
    df -= h->nLoci + 1;
    if(h->use_f)
	--df;

    return df;
}

/*
 * *a =  1 + 1/2 + 1/3 + ... + 1/(n-1)
 * *b =  1 + 1/(2*2) + 1/(3*3) + ... + 1/((n-1)*(n-1))
 */
void
sumRecip(double * a, double * b, unsigned n) {
    unsigned i;

    *a = *b = 0.0;
    for(i=1; i < n; ++i) {
	*a += 1.0/i;
	*b += 1.0/(i*i);
    }
}
