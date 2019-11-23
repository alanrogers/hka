#ifndef ARR_HKA_HPP_INCL
#include <boost/smart_ptr.hpp>
#include <cassert>
#include <cstring>
#include <cmath>
#include <iostream>
#endif // ARR_HKA_HPP_INCL

#ifndef ARR_HKA_HPP

class HKA {
    typedef boost::scoped_array<unsigned> unsigned_array_t;
    typedef boost::scoped_array<double>   double_array_t;
public:
    /// Constructor
    HKA(unsigned nLoci, const unsigned * nA, const unsigned * nB,
	const unsigned * nsA, const unsigned * nsB, const unsigned * nsAB,
	const double * SA, const double * SB, const double * D)
	: _nLoci(nLoci),
	  _nA(new unsigned[_nLoci]),
	  _nB(new unsigned[_nLoci]),
	  _nsA(new unsigned[_nLoci]),
	  _nsB(new unsigned[_nLoci]),
	  _nsAB(new unsigned[_nLoci]),
	  _SA(new double[_nLoci]),
	  _SB(new double[_nLoci]),
	  _D(new double[_nLoci]),
	  _aA(new double[_nLoci]),
	  _aB(new double[_nLoci]),
	  _bA(new double[_nLoci]),
	  _bB(new double[_nLoci]),
	  _theta(new double[_nLoci]) {

	// Values that never change are set here
	memcpy(&_nA[0], nA, _nLoci * sizeof(unsigned));
	memcpy(&_nB[0], nB, _nLoci * sizeof(unsigned));
	memcpy(&_nsA[0], nsA, _nLoci * sizeof(unsigned));
	memcpy(&_nsB[0], nsB, _nLoci * sizeof(unsigned));
	memcpy(&_nsAB[0], nsAB, _nLoci * sizeof(unsigned));

	for(unsigned i=0; i < _nLoci; ++i) {
	    sumRecip(&_aA[i], &_bA[i], _nA[i]);
	    sumRecip(&_aB[i], &_bB[i], _nB[i]);
	}

	// Values that change with each simulated tree are
	// set by "reset".
	reset(SA, SB, D);
    }

    /// Put data into the arrays and estimate.
    void
    reset(const double * SA, const double * SB, const double * D) {
	memcpy(&_SA[0], SA, _nLoci * sizeof(double));
	memcpy(&_SB[0], SB, _nLoci * sizeof(double));
	memcpy(&_D[0], D, _nLoci * sizeof(double));
	memset(&_theta[0], 0, _nLoci * sizeof(double));

	estimate();
    }

    /// Calculate the X^2 goodness-of-fit statistic.  (See p 154 of
    /// Hudson, Kreitman, and Aguade 1987.
    double
    X2() {
	double x=0;

	for(unsigned i=0; i < _nLoci; ++i) {
	    // Loci w/ no wi-group variation contribute nothing
	    if(_SA[i] == 0 && _SB[i] == 0 && _D[i] == 0.0)
		continue;
	    
	    double k = _theta[i]*_nsA[i];
	    double m = k * _aA[i];                   // mean
	    double v = m + k*k*_bA[i];               // var
	    double diff = _SA[i] - m;                // mean - var

	    x += diff*diff/v;                        // X^2 for SA
	    if(_nB[i] > 1) {
		k = _theta[i]*_nsB[i];
		m = k * _aB[i];                      // mean
		v = m + k*k*_bB[i];                  // var
		diff = _SB[i] - m;                   // mean - var
		x += diff*diff/v;                    // X^2 for SB

		k = _nsAB[i] * _theta[i];
		m = k*(_T + (1.0+_f)/2.0);
		k *= (1.0+_f)/2.0;
		v = m + k*k;
	    }else {
		k = _nsAB[i]*_theta[i];
		m = k*(_T + 1.0);
		v = m + k*k;
	    }
	    diff = _D[i] - m;
	    x += diff*diff/v;                        // X^2 for D
	}
	return x;
    }

    /// Estimate parameters theta_i, f, and T.  The function uses the
    /// EM algorithm to solve equations 5, on p 154 of Hudson,
    /// Kreitman, and Aguade 1987. 
    ///
    /// The system of equations must be re-written to account for variation in
    /// numbers of sites.  There are several ways to do this.  This function
    /// used the EM algorithm to solve the following system (which is modeled on
    /// eqn 6 of Hudson, Kreitman and Aguade):
    ///
    /// _nLoci       _nLoci
    /// Sum _SA[i] = Sum _theta[i]*_nsA[i]*_aA[i]
    /// i=1          i=1
    ///
    /// _nLoci            _nLoci
    /// Sum _SB[i] = _f * Sum _theta[i] * _nsB[i] * _aB[i]
    /// i=1               i=1
    ///
    /// _nLoci                        _nLoci
    /// Sum _D[i] = (_T + (1+_f)/2) * Sum _theta[i] * _nsAB[i]
    /// i=1                           i=1
    ///
    /// _SA[i] + _SB[i] + _D[i] 
    ///         = _theta[i]*(_aA[i]*_nsA[i] + _f*_aB[i]*_nsB[i]
    ///           + (_T + (1+_f)/2))*_nsAB[i]
    ///                                                       i=1,2,...,_nLoci-1.
    ///
    /// Loci at which _nB=1 do not contribute to the _SB sum.
    void
    estimate() {
	double_array_t hold_theta(new double[_nLoci]);
	double hold_T, hold_f;
	double tol = 1e-6, delta, x, y, sumA, sumB, sumD;
	unsigned i;

	sumA = sumB = sumD = x = 0.0;
	for(i=0; i < _nLoci; ++i) {
	    sumA += _SA[i];
	    sumB += _SB[i];
	    sumD += _D[i];
	    x += _nsA[i] * _aA[i];
	}
	x = sumA/x;
	for(i=0; i < _nLoci; ++i)
	    _theta[i] = x;         // initialize EM algorithm

	// loop until EM algorithm converges
	do{
	    hold_T = _T;
	    hold_f = _f;
	    memcpy(&hold_theta[0], &_theta[0], _nLoci * sizeof(double));

	    _use_f = false;

	    // Find f, T
	    _f = x = y = 0.0;
	    for(i=0; i < _nLoci; ++i) {
		if(_nB[i] > 1) {
		    _use_f = true;
		    x += _theta[i] * _nsB[i] * _aB[i];
		}
		y += _theta[i] * _nsAB[i];
	    }
	    if(_use_f) {
		_f = sumB / x;
		_T = sumD/y - (1.0+_f)/2.0;
	    }else
		_T = sumD/y - 1.0;

	    // Find theta[i]
	    for(i=0; i < _nLoci; ++i) {
		x = _SA[i] + _SB[i] + _D[i];
		if(_use_f) {
		    y = _aA[i]*_nsA[i] + _f*_aB[i]*_nsB[i]
			+ (_T + (1.0+_f)/2.0)*_nsAB[i];
		}else {
		    y = _aA[i]*_nsA[i] + (_T + 1.0)*_nsAB[i];
		}
		_theta[i] = x/y;
	    }

	    // delta measures change in this iteration
	    delta = fabs(_T - hold_T) + fabs(_f - hold_f);
	    for(i=0; i<_nLoci; ++i)
		delta += fabs(_theta[i] - hold_theta[i]);

	    // Algorithm ends when delta is small
	}while(delta > tol);
    }

    /// Report the value of parameter T
    double
    T() const { return _T;}

    /// Report the value of parameter f
    double
    f() const { return _f; }

    /// Report the value of parameter theta[i]
    double
    theta(unsigned i) const { return _theta[i]; }

    /// Return true if f has been estimated
    bool
    fEstimated() const { return _use_f; }

    /// Report degrees of freedom
    int
    degreesOfFreedom() const {
	int df=0;

	// count equations
	for(unsigned i=0; i < _nLoci; ++i) {
	    if(_SA[i] == 0 && _SB[i] == 0 && _D[i] == 0.0)
		continue;
	    if(_nB[i] > 1)
		df += 2; // 3 equations - 1 parameter
	    else
		df += 1; // 2 equations - 1 parameter
	}

	--df;           // subtract one df for T
	if(_use_f)
	    --df;       // subtract one df for f

	return df;
    }

private:

    /// *a =  1 + 1/2 + 1/3 + ... + 1/(n-1)
    /// *a2 =  1 + 1/(2*2) + 1/(3*3) + ... + 1/((n-1)*(n-1))
    void
    sumRecip(double * a, double * a2, unsigned n) {
	*a = *a2 = 0.0;
	for(unsigned i=1; i < n; ++i) {
	    *a += 1.0/i;
	    *a2 += 1.0/(i*i);
	}
    }

    unsigned _nLoci;   // number of loci
    unsigned_array_t _nA;    // haploid sample sizes from group A 
    unsigned_array_t _nB;    // haploid sample sizes from group B
    unsigned_array_t _nsA;   // # sites in within-A comparisons
    unsigned_array_t _nsB;   // # sites in within-B comparisons
    unsigned_array_t _nsAB;  // # sites in between-pop comparisons
    double_array_t _SA;      // segregating sites at  group A 
    double_array_t _SB;      // segregating sites at group B 
    double_array_t   _D;     // mean difference between A and B
    double_array_t   _aA;    // Sum 1/i for A
    double_array_t   _aB;    // Sum 1/i for B
    double_array_t   _bA;    // Sum 1/(i*i) for A
    double_array_t   _bB;    // Sum 1/(i*i) for B
    double_array_t   _theta; // parameter: theta_i = 4*N_i*u_i
    double           _T;     // parameter: time since species separated
    double           _f;     // parameter: _f = (pop size B)/(pop size A)
    bool             _use_f; // true if _f is estimated
};

    

#endif // ARR_HKA_HPP
