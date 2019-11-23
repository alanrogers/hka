/**
   @file wibtw.hpp
   
   @brief Header for class WiBtw, which measures the differences within
   and between populations on a tree.

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
#ifndef WIBTW_HPP_INCL
#define WIBTW_HPP_INCL

#include <iosfwd>
#include <gtree/rand.hpp>
#include <gtree/epoch.hpp>
#include <gtree/pophist.hpp>
#include <gtree/treenode.hpp>
#include <gtree/genetree.hpp>
#include <gtree/infsitedf.hpp>
#include <gtree/imcoales.hpp>
#include <arrmatrix/trimat.hpp>

#include <boost/scoped_array.hpp>

#endif

#ifndef WIBTW_HPP
#define WIBTW_HPP

namespace arr_gtree {
    /**
       WiBtw: Calculate differences within and between populations
       on a gene tree.
    **/
    class WiBtw {
    private:
	// _nPops is the number of pops.  I call them pops rather than demes
	// within this class because no assumption is made here about
	// random mating.
	unsigned _nPops;
	
	// _hapN[i] is haploid sample size of sample from pop i.
	boost::scoped_array<unsigned> _hapN;

	// Size of _work vector
	unsigned _workSize;
	
	// temporary storage used in calculations
	boost::scoped_array<unsigned> _work;

	// _btw[i][j] (where j > i ) is the mean difference between pops
	// i and j.
	arr_matrix::UTriMatrix<double> _btw;

	// _S[i] is the number of segregating (i.e. polymorphic) sites
	// within pop i.
	boost::scoped_array<unsigned> _S;

	////////////////////////////////////////////////////////////////

	/// Copy constructor is never defined
	WiBtw(const WiBtw & rhs);

	/// Assignment operator is never defined
	WiBtw &
	operator=(const WiBtw & rhs);

    public:

	/// Return the mean pairwise difference between pops i and j.
	double
	btw(unsigned i, unsigned j) {

	    assert(i < _nPops);
	    assert(j < _nPops);
	    assert(i != j);

	    if(j > i)
		return _btw[i][j];

	    return _btw[j][i];
		
	}

	/// Return the number of polymorphic sites within pop i.
	unsigned
	wi(unsigned i) {  return _S[i];	}

	/// Calculate the differences wi and btw pops in a given tree
	void
	calcDiff(TreeNode * root) {
	    unsigned i, j;
	    
	    memset(&_work[0], 0, _nPops * sizeof(unsigned));
	    memset(&_S[0], 0, _nPops * sizeof(unsigned));
	    _btw.zero();

	    traverse(root, &_work[0], &_work[_nPops]);

#ifndef NDEBUG
	    // On return from traverse, _work should equal _hapN.
	    for(i=0; i < _nPops; ++i) {
		assert(_work[i] == _hapN[i]);
	    }
#endif	    

	    // After traverse, _btw[i][j] contains the number of
	    // differences between pairs of chromosomes from pops i
	    // and j.  We need to divide the the number of pairs
	    // in order to make this into a mean.
	    for(i=0; i < _nPops; ++i)
		for(j=i+1; j < _nPops; ++j)
		    _btw[i][j] /= _hapN[i]*_hapN[j];
	}

	/// Traverse tree
	void
	traverse(TreeNode * node, unsigned * n, unsigned * w) {

	    unsigned i, j;

	    assert(node != 0);

	    // At leaf
	    if(node->left() == 0) {
		assert(node->right() == 0);
		++n[node->pop()];
		return;
	    }

	    // At internal node
	    memset(w, 0, _nPops * sizeof(unsigned));
	    traverse(node->left(), w, w + _nPops);
	    traverse(node->right(), w, w + _nPops);

	    // Now w[i] is the number of this node's descendants that
	    // are in pop i.
	    for(i=0; i < _nPops; ++i) {

		// count mutations within each pop
		if(w[i] > 0 && w[i] < _hapN[i])
		    _S[i] += node->nMutations();

		// Count between-pop differences.
		for(j=i+1; j < _nPops; ++j) {
		    _btw[i][j] += w[i]*(_hapN[j]-w[j])*node->nMutations();
		    _btw[i][j] += (_hapN[i]-w[i])*w[j]*node->nMutations();
		}
	    }

	    // Increment n, which carries the counts below this node
	    // back up to the parent node.
	    for(i=0; i < _nPops; ++i)
		n[i] += w[i];
	}
	
	/// Construct from nPops, the number of pops, and
	/// hapN, a vector of haploid pop sample sizes.
	explicit
	WiBtw(unsigned nPops, unsigned * hapN)
	    : _nPops(nPops),
	      _hapN(new unsigned[nPops]),
	      _btw(nPops),
	      _S(new unsigned[nPops])  {
	    memcpy(&_hapN[0], hapN, _nPops * sizeof(unsigned));

	    // Calculate total sample size
	    unsigned n=0;
	    for(unsigned i=0; i < _nPops; ++i) {
		n += _hapN[i];
	    }

	    // Allocate _work vector.  The _work vector is a concatenated
	    // list of vectors, each of size _nPops.  At worst, we
	    // need one _nPop-vector for each internal node and one
	    // for the last leaf visited.  Since there are n-1 internal
	    // nodes, we need n vectors each of length _nPops.
	    _workSize = n * _nPops;
	    _work.reset( new unsigned[_workSize] );
	}

	/** Destructor */
	~WiBtw(){}
    };
}
#endif // WIBTW_HPP


