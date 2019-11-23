#include "wibtw.hpp"
#include <gtree/isimtree.hpp>
#include <arrin/icistream.hpp>
#include <fstream>
#undef NDEBUG
using namespace arr_gtree;
using namespace arr_in;
using namespace std;

// convenience typedef
typedef Epoch<double,IslandModel> EpochType;
typedef PopHist<EpochType> PHist;

// Input specifies an enormous migration rate so that the
// population is effectively mating at random.
static const char * testInput = 
"2 # demes\n\
#len theta0 theta1 M\n\
Inf       1      1 1000\n";

static const char * fname = "xwibtw.tmp";

void    usage(void)
{
    cout << "usage: xwibtw [options]\n";
    cout << "  where options may include:\n";
    cout << "   -i nIterations\n";
    cout << "   -v <for verbose>\n";
    cout << "   -vv <for very verbose>\n";
    exit(1);
}

int     main(unsigned argc, char **argv)
{
    try{
	unsigned i, j;
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


	// The number of demes
	unsigned nDemes = ph.begin()->nDemes();
	unsigned demeSize = 10;  // haploid sample sizes

	// A vector with nDemes entries, each equal to demeSize
	vector<unsigned> sampleSize(nDemes, demeSize);

	// make sure
	assert(nDemes == sampleSize.size());
	for(i=0; i < nDemes; ++i)
	    assert(sampleSize[i] == demeSize);

	unsigned nLeaves = nDemes * demeSize;
	double  theta = ph[0].theta(0);

	BaseGeneratorType r;
	Uniform01 u(r);
	Poisson poisson(r);

	// Construct a tree
	ISIMtree *tPtr = new ISIMtree(sampleSize, ph, poisson, u);
	tPtr->sanityCheck();

	assert(tPtr->nLeaves() == nLeaves);

	if(verbose) {
	    cout << "nLeaves = " << nLeaves
		 << " tPtr->nLeaves() = " << tPtr->nLeaves()
		 << endl;
	}

	// Construct a WiBtw object and use it to measure differences
	// in tree just constructed.
	WiBtw wibtw(nDemes, &sampleSize[0]);
	wibtw.calcDiff(tPtr->root());

	cout << "S:" << flush;
	for(i=0; i < nDemes; ++i)
	    cout << " " << wibtw.wi(i) << flush;
	cout << endl;

	cout << "Mean pairwise diffs btw groups:" << endl;
	for(i=0; i < nDemes; ++i) {
	    for(j = 0; j < nDemes; ++j) {
		if(j <= i) {
		    cout << setw(10) << 0;
		}else
		    cout << setw(10) << wibtw.btw(i,j);
	    }
	    cout << endl;
	}

	delete tPtr;
	tPtr = 0;

	double a=0.0;
	for (i = 1; i < demeSize; i++) {
	    a += 1.0 / i;	    // sum of 1/i
	}

	// The theoretical values are those for random mating in a
	// population in which 2*N*u = nDemes*theta. This is
	// appropriate because the migration rate assumed here is very
	// large.  See the definition of testInput at the top of this
	// file.
	double meanSTheory = nDemes * theta * a;
	double meanBTheory = nDemes * theta;

	if(verbose)
	    cout << " nIterations=" << nIterations << endl;

	double  meanS = 0.0, meanB = 0.0;

	// Build a series of ISIMtrees
	for (int iteration = 0; iteration < nIterations; iteration++) {
	    if (veryVerbose)
		cout << "Iteration " << iteration << endl;

	    tPtr = new ISIMtree(sampleSize, ph, poisson, u);
	    tPtr->sanityCheck();

	    if (veryVerbose) 
		cout << *tPtr << endl;    // print ISIMtree

	    wibtw.calcDiff(tPtr->root()); // calculate differences
	    
	    for(i=0; i < nDemes; ++i) {   // keep track of sums
		meanS += wibtw.wi(i);
		for(j=i+1; j < nDemes; ++j)
		    meanB += wibtw.btw(i,j);
	    }
	    
	    delete tPtr;                  // free memory
	    tPtr = 0;
	}				  // end loop
	meanS /= nIterations * nDemes;
	meanB /= nIterations * nDemes * (nDemes-1) * 0.5;

	if(verbose) {
	    cout << "Did " << nIterations << " iterations. a="
		 << a << "\n";
	    cout << "meanS = " << meanS << "; should be " << meanSTheory
		 << endl;
	    cout << "meanB = " << meanB << "; should be " << meanBTheory
		 << endl;
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

    cout << "WiBtw OK" << endl;
    return 0;
}

