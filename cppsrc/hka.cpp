/**
   @file hka.cpp
   @brief Perform HKA test of selective neutrality
   @internal
   Copyright (C) 2003 Alan R. Rogers

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

   Alan R. Rogers, Department of Anthropology, University of Utah,
   Salt Lake City, UT 84112. rogers@anthro.utah.edu
**/
#include "hka.hpp"
#include "wibtw.hpp"
#include <boost/scoped_array.hpp>
#include <gtree/epoch.hpp>
#include <gtree/pophist.hpp>
#include <gtree/gtexcept.hpp>
#include <gtree/isimtree.hpp>
#include <arrin/ioexcept.hpp>
#include <arrin/icistream.hpp>
#include <ctime>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include  "hkaversion.hpp"
#define NITERATIONS 0

using namespace boost;
using namespace std;
using namespace arr_gtree;

// convenience typedef
typedef Epoch<double,IslandModel> EpochType;
typedef PopHist<EpochType> PHist;
typedef PHist * PHist_ptr;
typedef vector<unsigned> u_vector;
typedef u_vector * u_vector_ptr;

void    usage(void)
{
    cout << "usage: hka [options] input_file\n"
	 << "  where options may include:\n"
	 << "  -v     verbose\n"
	 << "  -i <n> do <n> iterations (default: "
	 << NITERATIONS << ")\n"
	 << endl;
    exit(1);
}

enum optndx {opt_niterations, opt_verbose};
typedef map<string, enum optndx> optmap;

int     main(unsigned argc, char **argv)
{
    unsigned locus, nIterations = NITERATIONS;
    bool verbose = false;
    
    // Process command line arguments
    optmap option;
    option["-i"] =  opt_niterations;
    option["-v"] =  opt_verbose;

    // pointer to file name
    char * fname = 0;

    for (unsigned i = 1; i < argc; ++i) {
	if(argv[i][0] == '-') {
	    optmap::iterator pos = option.find(argv[i]);
	    if(pos == option.end()) {
		cout << "unknown command-line arg: " << argv[i] << endl;
		usage();
	    }
	    switch(pos->second) {
		case opt_verbose:
		    verbose = true;
		    break;
		case opt_niterations:
		    if (++i >= argc)
			usage();
		    nIterations = strtol(argv[i], 0, 10);
		    break;
		default:
		    usage();
	    }
	}else {
	    if(fname != 0)
		usage();   // only one file name allowed
	    fname = strdup(argv[i]);
	}
    }

    if(fname == 0)
	usage();
	
    try{
	ifstream phis(fname); // input file stream

	if(!phis.good()) {
	    cerr << "Couldn't open " << fname << " for input" << endl;
	    exit(1);
	}

	arr_in::ICstreambuf phissb(phis.rdbuf()); // strip comments
        arr_in::ICistream in(&phissb);        // comment-free input file stream


	unsigned nLoci;
	
	
	// Read Input
	in >> nLoci;

	// Allocate arrays
	scoped_array<string> locusName(new string[nLoci]);
	scoped_array<unsigned> nsA(new unsigned[nLoci]);  // number of sites A
	scoped_array<unsigned> nsB(new unsigned[nLoci]);  // number of sites B
	scoped_array<unsigned> nsAB(new unsigned[nLoci]); // number of sites AB
	scoped_array<unsigned> nA(new unsigned[nLoci]);   // haploid sample size A
	scoped_array<unsigned> nB(new unsigned[nLoci]);   // haploid sample size B
	scoped_array<double> SA(new double[nLoci]);       // polymorphic sites A
	scoped_array<double> SB(new double[nLoci]);       // polymorphic sites B
	scoped_array<double> D(new double[nLoci]);        // mean diff AB

	// Finish reading input
	for(locus=0; locus < nLoci; ++locus) {
	    in >> locusName[locus];
	    in >> nsA[locus];
	    in >> nsB[locus];
	    in >> nsAB[locus];
	    in >> nA[locus];
	    in >> nB[locus];
	    in >> SA[locus];
	    in >> SB[locus];
	    in >> D[locus];
	}

	cout << "\n########## HKA Version "
	     << HKA_VERSION
	     << " ##########" << endl;
	// Echo command line
	cout << "# Command line:";
	for(unsigned i=0; i < argc; ++i)
	    cout << " " << argv[i];
	cout << "\n" << endl;

	// Echo data
	unsigned w = 6;
	cout << setw(w) << nLoci << " # loci" << endl;
	cout << setw(10) << "#    Locus"
	     << " " << setw(w) << "nsA"
	     << " " << setw(w) << "nsB"
	     << " " << setw(w) << "nsAB"
	     << " " << setw(w) << "nA"
	     << " " << setw(w) << "nB"
	     << " " << setw(w) << "SA"
	     << " " << setw(w) << "SB"
	     << " " << setw(w) << "D" << endl;
	for(locus=0; locus < nLoci; ++locus) {
	    cout << setw(10) << locusName[locus]
		 << " " << setw(w) << nsA[locus]
		 << " " << setw(w) << nsB[locus]
		 << " " << setw(w) << nsAB[locus]
		 << " " << setw(w) << nA[locus]
		 << " " << setw(w) << nB[locus]
		 << " " << setw(w) << SA[locus]
		 << " " << setw(w) << SB[locus]
		 << " " << setw(w) << D[locus]
		 << endl;
	}
	cout << "#\n"
	     << "# nsA  Number of sites used w/i species A\n"
	     << "# nsB  w/i species B\n"
	     << "# nsAB in comparison btw A and B\n"
	     << "# nA   Haploid sample size for A\n"
	     << "# nB   for B\n"
	     << "# SA   Polymorphic sites w/i A\n"
	     << "# SB   w/i B\n"
	     << "# D    Mean number of site differences btw A and B"
	     << endl;

	// Create object of type HKA
	HKA hka(nLoci, &nA[0], &nB[0], &nsA[0], &nsB[0], &nsAB[0],
		&SA[0], &SB[0], &D[0]);

	// observed X2 value
	double obsX2 = hka.X2();

	// Print results
	w = 15;
	cout << "\n## RESULTS ##" << endl;
	cout << "#        T = " << hka.T() << endl;
	cout << "#        f = " << hka.f() << endl;
	for(unsigned i=0; i < nLoci; ++i)
	    cout << "# theta[" << i << "] = " << hka.theta(i) << endl;
	cout << "#       df = " << hka.degreesOfFreedom() << endl;
	cout << "#      X^2 = " << obsX2 << endl;

	// Set up array of sample size (ss) vectors, one per locus
	u_vector_ptr * ss = new u_vector_ptr[nLoci];
	for(locus=0; locus < nLoci; ++locus) {
	    ss[locus] = new u_vector(2);
	    ss[locus]->clear();
	    ss[locus]->push_back(nA[locus]);
	    ss[locus]->push_back(nB[locus]);
	}

	// Set up PHist array.  Each entry contains the population
	// history parameters for one locus.
	EpochType ep(2);                        // epoch with 2 Demes
	PHist_ptr * ph = new PHist_ptr[nLoci];  // array of ptrs to PHist
	for(locus=0; locus < nLoci; ++locus) {
	    ph[locus] = new PHist(2);               // Each PHist has 2 demes

	    // define epoch 0
	    double t=0.0;
	    double theta0 = hka.theta(locus) * nsA[locus];
	    double theta1 = theta0;
	    double length = hka.T() * theta0;
	    if(hka.fEstimated())
		theta1 *= hka.f();
	    ep.t( t );             // backwards time at which epoch begins
	    ep.length( length );   // length of epoch
	    ep.theta(0, theta0 );  // set theta for pop 0
	    ep.theta(1, theta1);   // set theta for pop 1

	    // push epoch 0 onto ph[locus]
	    ph[locus]->clear();
	    ph[locus]->push_back(ep);

	    // define epoch 1
	    t= length;
	    length = HUGE_VAL;
	    theta1 = 0.0;
	    if(hka.fEstimated())
		theta0 *= 0.5*(1.0 + hka.f());

	    ep.t( t );             // backwards time at which epoch begins
	    ep.length( length );   // length of epoch
	    ep.theta(0, theta0 );  // set theta for pop 0
	    ep.theta(1, theta1);   // set theta for pop 1

	    ph[locus]->push_back(ep);

	    if(verbose) {
		cout << "# Population History for locus "
		     << locus << ":\n"
		     << *ph[locus] << endl;
	    }
	}

	// set up random number generator
	BaseGeneratorType r;
	// read clock to set seed
	r.seed(time(0));
	Uniform01 u(r);
	Poisson poisson(r);

	// An array ow WiBtw objects, which are used to measure
	// wi- and btw-pop differences
	WiBtw ** wibtw = new (WiBtw *)[nLoci];
	for(locus=0; locus < nLoci; ++locus)
	    wibtw[locus] = new WiBtw(2, &(*ss[locus])[0]);

	// Count number of iterations in which X2 > obsX2
	unsigned long upperTail = 0;

	if(nIterations == 0)
	    return 0;

	cout.flush();
	cerr << "Simulating " << nIterations << " data sets:" << endl;

	// simulation
	ISIMtree * tPtr = 0;
	for(unsigned iteration=0; iteration < nIterations; ++iteration) {
	    if(!verbose && iteration%10 == 0)
		cerr << '.' << flush;
	    for(locus=0; locus < nLoci; ++locus) {
		tPtr = new ISIMtree(*ss[locus], *ph[locus], poisson, u);
		tPtr->sanityCheck();
		wibtw[locus]->calcDiff(tPtr->root());

		// Put data from tree into SA, SB, and D
		SA[locus] = wibtw[locus]->wi(0);
		SB[locus] = wibtw[locus]->wi(1);
		D[locus]  = wibtw[locus]->btw(0,1);

		delete tPtr;
	    }
	    hka.reset(&SA[0], &SB[0], &D[0]);

	    double simX2 = hka.X2();

	    if(verbose) {
		cout << iteration << ": X2 = " << simX2 << endl;
	    }
	    
	    if(simX2 >= obsX2)
		++upperTail;
	}
	cerr << endl;

	// Probability the X2 >= observed value
	double p = static_cast<double>(upperTail)/nIterations;
	cout << "#    sim p = " << p << endl;

    }catch(invalid_argument & e){
	cerr << "invalid_argument: " << e.what() << endl;
	exit(2);
    }catch(std::bad_alloc & e) {
	cerr << "Ran out of memory" << endl;
	exit(2);
    }catch(exception &e){
	cerr << "Generic exception: " << e.what() << endl;
	exit(2);
    }catch(...){
	cerr << "Unknown exception" << endl;
	exit(2);
    }
    return 0;
}
