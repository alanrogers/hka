#ifndef ARR_HKA_H

/** Structure containing the data used in HKA calculations */
typedef struct {
    unsigned   nLoci; /* number of loci */
    unsigned * nsA;   /* # sites in within-A comparisons */
    unsigned * nsB;   /* # sites in within-B comparisons */
    unsigned * nsAB;  /* # sites in between-pop comparisons */
    unsigned * nA;    /* haploid sample sizes from group A */
    unsigned * nB;    /* haploid sample sizes from group B */
    double   * SA;    /* segregating sites at  group A */
    double   * SB;    /* segregating sites at group B */
    double   * D;     /* mean difference between A and B */
    double   * aA;    /* Sum 1/i for A */
    double   * aB;    /* Sum 1/i for B */
    double   * bA;    /* Sum 1/(i*i) for A */
    double   * bB;    /* Sum 1/(i*i) for B */
    double   * theta; /* parameter: theta_i = 4*N_i*u_i */
    double   * tmp;   /* temporary storage */
    double     T;     /* parameter: time since species separated */
    double     f;     /* parameter: _f = (pop size B)/(pop size A) */
    int        use_f; /* true if _f is estimated */
} HKA;

HKA *
HKA_construct(unsigned nLoci, const unsigned * nA, const unsigned * nB,
	      const unsigned * nsA, const unsigned * nsB,
	      const unsigned * nsAB, const double * SA, const double * SB,
	      const double * D);

/** Free memory allocated by h */
void
HKA_destroy(HKA * h);

/** Put data into the arrays and estimate. */
void
HKA_reset(const unsigned * nA, const unsigned * nB,
      const unsigned * nsA, const unsigned * nsB, const unsigned * nsAB,
      const double * SA, const double * SB, const double * D,
	  HKA * h);

/*
 * Calculate the X^2 goodness-of-fit statistic.  (See p 154 of
 * Hudson, Kreitman, and Aguade 1987.
 */
double
HKA_X2(HKA * h);

void
HKA_estimate(HKA * h);

/* Report the value of parameter T */
double
HKA_T(HKA * h);

/* Report the value of parameter f */
double
HKA_f(HKA * h);

/* Report the value of parameter theta[i] */
double
HKA_theta(unsigned i, HKA * h);

/* Report degrees of freedom */
unsigned
HKA_degreesOfFreedom(HKA * h);

/*
 * *a =  1 + 1/2 + 1/3 + ... + 1/(n-1)
 * *b =  1 + 1/(2*2) + 1/(3*3) + ... + 1/((n-1)*(n-1))
 */
void
sumRecip(double * a, double * b, unsigned n);

#endif /* ARR_HKA_H */
