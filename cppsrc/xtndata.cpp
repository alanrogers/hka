#include "gtree/epoch.hpp"
#include "gtree/pophist.hpp"
#include <arrin/icistream.hpp>
#include "tndata.hpp"
#include <fstream>
#undef NDEBUG
using namespace arr_gtree;
using namespace arr_in;
using namespace std;

// convenience typedef
typedef Epoch<double,IslandModel> EpochType;
typedef PopHist<EpochType> PHist;

static const char * testInput = 
"2 # demes\n\
#len theta0 theta1 M\n\
10      1     1     0\n\
Inf     1     0       0\n";

static const char * stationaryTestInput = 
"1 # deme\n\
#len theta0 M\n\
Inf  1	    0\n";

static const char * fname = "xisimtree.tmp";

void    usage(void)
{
    cout << "usage: xisimtree [options]\n";
    cout << "  where options may include:\n";
    cout << "   -i nIterations\n";
    cout << "   -v <for verbose>\n";
    cout << "   -vv <for very verbose>\n";
    exit(1);
}

int     main(unsigned argc, char **argv)
{
    try{
	unsigned i;
	int      nIterations = 3;
	bool     verbose = false;
	bool     veryVerbose = false;

	for (i = 1; i < argc; i++) {
	    if (argv[i][0] == '-') {
		if (strcmp(argv[i] + 1, "i") == 0) {
		    if (++i >= argc)
			usage();
		    nIterations = strtol(argv[i], 0, 10);
		} else if (strcmp(argv[i] + 1, "v") == 0) {
		    verbose = true;
		} else if (strcmp(argv[i] + 1, "vv") == 0) {
		    veryVerbose = true;
		} else
		    usage();
	    } else
		usage();
	}

	if(veryVerbose)
	    verbose = true;

	ofstream ofs(fname);    // output file stream
	ofs << testInput ;
	ofs.close();

	ifstream ifs(fname); // input file stream
	ICstreambuf icsb(ifs.rdbuf()); // strip comments
        ICistream cfifs(&icsb);        // comment-free input file stream

	PHist ph;            // empty PopHist
	cfifs >> ph;         // Read from input stream
	ifs.close();         // avoid segmentation fault
	if(verbose) {
	    cout << "Read: " << ph << endl;
	}

	vector<unsigned> sampleSize;

	// The number of demes
	unsigned nDemes = ph.begin()->nDemes();

	// Construct a vector of sample sizes.  There must be as many
	// samples as there are demes in Epoch ph.begin().
	assert(sampleSize.size() == 0);
	for(i=1; i <= nDemes; ++i)
	    sampleSize.push_back(i * 10);

	assert(nDemes == sampleSize.size());

	unsigned nLeaves = 0;
	for(vector<unsigned>::const_iterator itr=sampleSize.begin();
	    itr != sampleSize.end();
	    ++itr)
	    nLeaves += *itr;

	double  theta = ph[0].theta(0);

	BaseGeneratorType r;
	Uniform01 u(r);
	Poisson poisson(r);

	// Upper-triangular matrix of segregating-site counts
	arr_matrix::UTriMatrix<unsigned> seg(nDemes);
	seg.zero();

	// Basic tests
	GTree *tPtr = 
	    new GTree(sampleSize, ph, poisson, u);
	tPtr->sanityCheck();

	assert(tPtr->nLeaves() == nLeaves);

	setLeafward(tPtr->root());
	setRootward(tPtr->root());
	segregatingSites(tPtr->root(), seg);

	cout << "Segregating site matrix:\n" << seg << endl;

	cout << "Segregating sites overall: "
	     << tPtr->root()->treeMutations() << endl;

	GTree t2(*tPtr);
	t2.sanityCheck();

	GTree *t3Ptr = 
	    new GTree(sampleSize, ph, poisson, u);

	t3Ptr->sanityCheck();

	GeneTree gt1(*tPtr), gt2(*t3Ptr);
	gt1.sanityCheck();
	gt2.sanityCheck();

	gt1 = gt2;
	gt1.sanityCheck();

	gt1.mutate(poisson, u);
	gt1.sanityCheck();

	if(veryVerbose)
	    cout << *tPtr << endl;

	delete tPtr;
	tPtr = 0;

	// Check statistical properties with stationary random-mating
	// coalescent.

	// For this we need a different PopHist.
	ofs.open(fname);
	ofs << stationaryTestInput ;
	ofs.close();

	ifs.open(fname); // input file stream
	ICstreambuf icsb2(ifs.rdbuf()); // strip comments
        ICistream cfifs2(&icsb2);       // comment-free input file stream

	ph.clear();           // empty PopHist
	cfifs2 >> ph;         // Read from input stream
	ifs.close();          // avoid segmentation fault
	ph[0].theta(0,theta); // set value of theta
	if(verbose) {
	    cout << "New history: " << ph << endl;
	}

	// Sample size vector must contain only 1 entry
	sampleSize.resize(1);
	sampleSize[0] = nLeaves;
	
	double a=0.0, b=0.0;
	for (i = 1; i < nLeaves; i++) {
	    a += 1.0 / i;	    // sum of 1/i
	    b += 1.0 / (i * i);     // sum of 1/(i*i)
	}
	double meanDTheory = 2 * theta * (1.0 - 1.0 / nLeaves);
	double varDTheory = 0.0;
	double varLTheory = 0.0;
	double x;
	for (i = 2; i <= nLeaves; i++) {
	    x = 2.0 * theta / (i * (i - 1));
	    varDTheory += x * x;
	    varLTheory += x * x * i * i;
	}
	double meanLTheory = 2.0 * theta * a;

	if(verbose)
	    cout << " nIterations=" << nIterations << endl;

	double  meanD = 0.0, meanL = 0.0, varD = 0.0, varL = 0.0;

	// Build a series of GTrees, measure depth and length
	for (int iteration = 0; iteration < nIterations; iteration++) {
	    if (veryVerbose)
		cout << "Iteration " << iteration << endl;

	    tPtr = new GTree(sampleSize, ph, poisson, u);
	    tPtr->sanityCheck();

	    if (veryVerbose) 
		cout << *tPtr << endl;  // print GTree

	    x = tPtr->depth();          // depth of GTree
	    meanD += x;
	    varD += x * x;

	    x = tPtr->length();         // length of all branches
	    meanL += x;
	    varL += x * x;

	    delete tPtr;                // free memory
	    tPtr = 0;
	}				// end iteration
	meanD /= nIterations;
	meanL /= nIterations;
	varD = (varD - nIterations * meanD * meanD) / (nIterations - 1);
	varL = (varL - nIterations * meanL * meanL) / (nIterations - 1);

	if(verbose) {
	    cout << "Did " << nIterations << " iterations. a=" << a
		 << " b=" << b << "\n";
	    cout << "meanD = " << meanD << "; should be " << meanDTheory
		 << endl;
	    cout << "varD  = " << varD << "; should be " << varDTheory
		 << endl;
	    cout << "meanL = " << meanL << "; should be " << meanLTheory
		 << endl;
	    cout << "varL  = " << varL << "; should be " << varLTheory
		 << endl << endl;
	}
    }catch(invalid_argument & e){
	cerr << "invalid_argument: " << e.what() << endl;
	exit(2);
    }catch(arr_in::IOException & e){
	cerr << "IOException: " << e.what() << endl;
	exit(2);
    }catch(arr_gtree::NoConvergenceException & e){
	cerr << "NoConvergenceException: " << e.what() << endl;
	exit(2);
    }catch(exception &e){
	cerr << "Generic exception: " << e.what() << endl;
	exit(2);
    }catch(...){
	cerr << "Unknown exception" << endl;
	exit(2);
    }

    cout << "GTree OK" << endl;
    return 0;
}

