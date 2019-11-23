/**
   @file tndata.hpp
   
   @brief Header for class TNData, which describes the
   state of a TreeNode under the model of Infinite Sites.

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
#ifndef ISTNDATA_HPP_INCL
#define ISTNDATA_HPP_INCL

#include <iosfwd>
#include <set>
#include <gtree/rand.hpp>
#include <gtree/epoch.hpp>
#include <gtree/pophist.hpp>
#include <gtree/treenode.hpp>
#include <gtree/genetree.hpp>
#include <gtree/infsitedf.hpp>
#include <gtree/imcoales.hpp>
#include <arrmatrix/trimat.hpp>

#endif

#ifndef ISTNDATA_HPP
#define ISTNDATA_HPP

namespace arr_gtree {
    /**
       ISTNData: Describes the mutational change along a single
       branch within a infinite sites GeneTree.
    **/
    class ISTNData : public InfSiteDiff {
    private:
	std::set<int> _leafward;
	std::set<int> _rootward;

	friend void
	setRootward(TreeNode * node);

	friend void
	setLeafward(TreeNode * node);

	friend void
	segregatingSites(TreeNode * node,
			 arr_matrix::UTriMatrix<unsigned> & s);

	friend bool
	polymorphic(int i, ISTNData & data);

    public:
	typedef std::set<int>::iterator iterator;
	
	/**
	   Construct from the count of mutations.  Used by
	    operator+
	**/
	explicit
	ISTNData(int nMutations = 0)
	    : InfSiteDiff(nMutations) {
	}

	/** Copy constructor */
	ISTNData(const ISTNData & rhs)
	    : InfSiteDiff(rhs),
	      _leafward(rhs._leafward),
	      _rootward(rhs._rootward) {}

	/** Destructor */
	~ISTNData(){}

	ISTNData &
	operator=(const ISTNData & rhs){
	    if(this != &rhs) {
		InfSiteDiff::operator=(rhs);
		_leafward = rhs._leafward;
		_rootward = rhs._rootward;
	    }
	    return *this;
	}

	/// virtual assignment
	virtual TreeNodeData &
	assign(const TreeNodeData & rhs) {
	    return operator=(dynamic_cast<const ISTNData &>(rhs));
	}

	/** @return a pointer to a newly allocated copy of *this. */
	virtual ISTNData *
	clone() const { return ( new ISTNData(*this) ); }
    };

    void
    setLeafward(TreeNode * node) {

	ISTNData * ego;
	ISTNData * child;
	
	// Assume w/o check that node is a valid object.
	ego = dynamic_cast<ISTNData *>(node->data());

	// At a leaf
	if(node->left() == 0) {
	    assert(node->right() == 0);
	    ego->_leafward.insert( node->pop() );
	    return;
	}

	// At an internal node
	setLeafward(node->left());
	child = dynamic_cast<ISTNData *>(node->left()->data());
	ego->_leafward.insert( child->_leafward.begin(),
			       child->_leafward.end());

	setLeafward(node->right());
	child = dynamic_cast<ISTNData *>(node->right()->data());
	ego->_leafward.insert( child->_leafward.begin(),
			       child->_leafward.end());
	
    }

    // At each internal node, add the list of descendants of the left
    // node to the list of ancestors of the right node and vice versa.
    // setLeafward should be called before setRootward.
    void
    setRootward(TreeNode * node) {

	ISTNData * right;
	ISTNData * left;
	
	// Ignore leaves
	if(node->left() == 0)
	    return;

	// At internal node
	left  = dynamic_cast<ISTNData *>(node->left()->data());
	right = dynamic_cast<ISTNData *>(node->right()->data());

	left->_rootward.insert(right->_leafward.begin(),
			       right->_leafward.end());

	right->_rootward.insert(left->_leafward.begin(),
				left->_leafward.end());

	setRootward(node->left());
	setRootward(node->right());
    }

    /// Is population i present in set v?
    bool
    inSet(int i, const std::set<int> & v) {
	return (v.find(i) != v.end());
    }

    /// Are the mutations associated with this branch polymorphic
    /// within population i?
    bool
    polymorphic(int i, ISTNData & data) {
	return (inSet(i, data._leafward) && inSet(i, data._rootward));
    }

    /// On return, s[i][j] is the number of sites that differ between
    /// populations i and j, and s[i][i] is the number that are
    /// polymorphic within population i.
    void
    segregatingSites(TreeNode * node,
		     arr_matrix::UTriMatrix<unsigned> & s) {

	if(node == 0)
	    return;

	ISTNData & data = *dynamic_cast<ISTNData *>(node->data());
	unsigned dim = s.nRows();
	
	for(unsigned i=0; i < s.nRows(); ++i) {
	    if( polymorphic(i, data) )
		s[i][i] += node->nMutations();

	    for(unsigned j=i+1; j<dim; ++j) {
		// This code counts btw-group differences even if
		// they are also polymorphic within populations.
		if( (inSet(i,data._leafward) && inSet(j,data._rootward))
		    || (inSet(j,data._leafward) && inSet(i,data._rootward)))
		    s[i][j] += node->nMutations();
	    }
	}
	    
	segregatingSites(node->left(), s);
	segregatingSites(node->right(), s);
    }

    /** Factory producing ISTNData objects */
    class ISTNDataFactory : public TNDFactory {
    public:
	ISTNData *
	operator()() const{
	    return new ISTNData();
	}
	
	ISTNDataFactory *
	clone() const {
	    return new ISTNDataFactory(*this);
	}
    };

    /**
       A concrete class of gene tree.  Uses the ISTNData class defined
       just above and constructs itself using the Island Model
       coalescent algorithm  

       @param _sampleSize A vector whose i'th entry is the number of
       leaves in subdivision i.
	   
       @param ph A reference to a PopHist object
	   
       @param u reference to a random number generator
    **/
    class GTree : public GeneTree {
    public:
        typedef Epoch<double,IslandModel> epoch_t;
        typedef const Epoch<double,IslandModel> const_epoch_t;

	/// Construct
	GTree(const std::vector<unsigned> & sampleSize,
		 const PopHist<epoch_t> & ph,
		 Poisson & poisson,
		 Uniform01 & u)
	    : GeneTree(
		new IMCoalescent<epoch_t>(sampleSize, ph, u,
					  ISTNDataFactory()
		    ) ) { mutate(poisson, u); }
    };
}
#endif // ISTNDATA_HPP


