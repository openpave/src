/**************************************************************************

	TREE.H - Binary search tree classes

	The contents of this file are subject to the Academic Development
	and Distribution License Version 1.0 (the "License"); you may not
	use this file except in compliance with the License.  You should
	have received a copy of the License with this file.  If you did not
	then please contact whoever distributed this file too you, since
	they may be in violation of the License, and this may affect your
	rights under the License.

	Software distributed under the License is distributed on an "AS IS"
	basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See
	the License for the specific language governing rights and
	limitations under the License.

	The Initial Developer of the Original Software is Jeremy Lea.

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements a templated C++ class for left leaning
		red-black trees, which are a special variant of binary search
		trees.

	History:
		2008/03/27 - Created a basic implementation.

**************************************************************************/

#ifndef __TREE_H
#define __TREE_H

#ifdef TEST_TREES
#include <stdio.h>
#endif

/*
 * struct BST
 *
 * A basic implementation of a left-leaning red/black binary search tree,
 * using individual allocations for each node.  Not very fancy, mostly here
 * to understand how the other trees work...  Will be removed.
 *
 * http://www.cs.princeton.edu/~rs/talks/LLRB/RedBlack.pdf
 */
template<class K, class V>
struct BST
{
	BST()
	  : root(nullptr) {
	}
	~BST() {
		delete root;
	}
	V * get(const K & k) {
		node * x = root;
		while (x != nullptr) {
			if (k == x->key)
				return &(x->value);
			else if (k < x->key)
				x = x->left;
			else
				x = x->right;
		}
		return nullptr;
	}
	void insert(const K & k, const V & v) {
		bool add = false;
		root = insert(root,k,v,&add);
		root->red = false;
	}
	void remove(const K & k) {
		if (root == nullptr)
			return;
		root = remove(root,k);
		if (root != nullptr)
			root->red = false;
	}
#ifdef TEST_TREES
	void print() {
		if (root == nullptr)
			printf("Empty!\n");
		else
			print(0,0,root);
	}
#endif

private:
	struct node {
		K key;
		V value;
		unsigned order;
		node * left, * right;
		bool red;
		node(const K & k, const V & v)
		  : key(k), value(v), order(0), left(nullptr), right(nullptr),
		    red(true) {
		}
		~node() {
			delete left;
			delete right;
		}
	} * root;

	node * insert(node * h, const K & k, const V & v, bool * o) {
		// If we're zero that means we need to make a new node...
		if (h == nullptr) {
			*o = true;
			return new node(k,v);
		}
		// Split double reds on the way down.
		if (h->left != nullptr && h->left->red) {
			if (h->left->left != nullptr && h->left->left->red) {
				h = rotR(h);
				h->left->red = false;
			}
		}
		if (k == h->key)
			// Key already exists, so replace.
			h->value = v;
		else if (k < h->key) {
			// Re-root insert into left tree.
			h->left = insert(h->left,k,v,o);
			if (*o)
				h->order++;
		} else {
			// Or right tree.
			h->right = insert(h->right,k,v,o);
		}
		if (h->right != nullptr && h->right->red)
			h = leanLeft(h);
		return h;
	}
	node * rotL(node * h) {
		node * x = h->right;
		h->right = x->left;
		x->left = h;
		x->order += h->order + 1;
		return x;
	}
	node * rotR(node * h) {
		node * x = h->left;
		h->left = x->right;
		x->right = h;
		h->order -= x->order + 1;
		return x;
	}
	node * leanLeft(node * h) {
		h = rotL(h);
		h->red = h->left->red;
		h->left->red = true;
		return h;
	}
	node * leanRight(node * h) {
		h = rotR(h);
		h->red = h->right->red;
		h->right->red = true;
		return h;
	}
	node * remove(node * h, const K & k) {
		if (k < h->key) {
			// If K is missing do nothing.
			if (h->left == nullptr)
				return h;
			// move red left is needed.
			if (!h->left->red
			 && (h->left->left == nullptr || !h->left->left->red)) {
				h->red = false;
				h->left->red = true;
				if (h->right->left != nullptr && h->right->left->red) {
					h->right = rotR(h->right);
					h = rotL(h);
				} else
					h->right->red = true;
			}
			// Re-root delete in left branch.
			h->left = remove(h->left,k);
		} else {
			// Lean red links right going down.
			if (h->left != nullptr && h->left->red)
				h = leanRight(h);
			// We've found our node, and it's a leaf.
			if (k == h->key && h->right == nullptr) {
				delete h;
				return nullptr;
			}
			// Push red right going down.
			if (!h->right->red
			 && (h->right->left == nullptr || !h->right->left->red)) {
				h->red = false;
				h->right->red = true;
				if (h->left->left != nullptr && h->left->left->red) {
					h = rotR(h);
					h->red = true;
					h->left->red = false;
				} else
					h->left->red = true;
			}
			// We've found our node, but it's not a leaf...
			// Replace it with the minimum right node and
			// then remove that node...
			if (k == h->key) {
				node * x = h->right;
				while (x->left != nullptr)
					x = x->left;
				h->key = x->key;
				h->value = x->value;
				h->right = remove(h->right,x->key);
			} else
				// Look on the right.
				h->right = remove(h->right,k);
		}
		// fixed up right leaning reds on the way up.
		if (h->right != nullptr && h->right->red)
			h = leanLeft(h);
		return h;
	}
#ifdef TEST_TREES
	void print(int level, unsigned order, node * h) {
		if (h->right != nullptr)
			print(level+1,order+h->order+1,h->right);
		printf("(%d %d %d) ",order,h->order,order+h->order);
		for (int i = 0; i < level; i++)
			printf(" ");
		printf("%d: %f %s\n",h->key.i,h->value.d,(h->red?"RED  ":"BLACK"));
		if (h->left != nullptr)
			print(level+1,order,h->left);
	}
#endif
};

/*
 * class ktree_avl
 *
 * A key/value tree, like ksset<K,V>, using the AVL implementation since
 * it is relatively simple.  The common code between this any the LLRB tree
 * needs to be abstracted out into a common base class.
 *
 * This does not shrink yet.
 *
 * This uses index values rather than pointers, and a fixed storage array
 * underneath, to minimize the use of small chunks of memory.  The key class
 * K must implement compare().
 *
 * See https://en.wikipedia.org/wiki/AVL_tree 
 */
template <class K, class V>
class ktree_avl {
public:
	// Make one...
	inline explicit ktree_avl()
	  : size(0), buffer(0), block(DFLT_BLK), root(UINT_MAX), value(nullptr) {
	}
	// Clean up.
	inline ~ktree_avl() {
		if (value)
			free(value);
	}
	// The length. Nice for lots of things...
	inline unsigned length() const {
		return size;
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned haskey(const K & k) const {
		unsigned x = root;
		while (x != UINT_MAX) {
			int cmp = k.compare(static_cast<const K &>(value[x]._v));
			if (cmp == 0)
				break;
			else if (cmp < 0)
				x = value[x].left;
			else
				x = value[x].right;
		}
		return x;
	}
	// We can use node numbers to find our nodes...
	inline V & operator[] (const unsigned p) const {
		return value[p]._v;
	}
	// Add node.
	inline void add(const V & v) {
		bool grew = false;
		
		allocate(size+1);
		append(root,v,&grew);
		if (grew)
			rebalance(root);
	}
#ifdef TEST_TREES
	void print() {
		if (root == UINT_MAX)
			printf("Empty!\n");
		else
			print(0,root);
	}
#endif

protected:
	unsigned size;             // The size of the set...
	unsigned buffer;           // The allocated buffer size...
	unsigned block;            // The minimum block size.
	unsigned root;             // Root of the tree.
	struct _V {
	private:
		V _v;
		
		friend class ktree_avl;

		int weight;            // Number of nodes under this one (including it).
		unsigned left, right;  // Left and right node numbers

		explicit _V(const V & v)
		  : _v(v), 
		    weight(1), left(UINT_MAX), right(UINT_MAX) {
		}
		// Placement new to support inplace init in the list.
		void * operator new(size_t, void * p) {
			return p;
		}
		void operator delete(void *, void *) {
		}
	} * value;                 // Take a guess...

	// Make some space...
	void allocate(const unsigned s) {
		while (s > 8*block)
			block *= 8;
		unsigned b = block*(s/block+(s%block?1:0));
		if (b == buffer)
			return;
		_V * temp = static_cast<_V *>(realloc(value,b*sizeof(_V)));
		if (temp == nullptr)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void append(unsigned & r, const V & v, bool * grew) {
		// If we're UINT_MAX that means we need to make a new node...
		if (r == UINT_MAX) {
			new(&value[size]) _V(v);
			r = size++;
			*grew = true;
			return;
		}
#ifdef TEST_TREES
		printf("Append @%d before:\n",r);
		print();
#endif
		int cmp = static_cast<const K &>(v).compare(value[r]._v);
		if (cmp == 0) {
			value[r]._v = v;
			return;
		} else if (cmp < 0) {
			// Re-root insert into left tree.
			append(value[r].left,v,grew);
			if (*grew)
				value[r].weight++;
		} else {
			append(value[r].right,v,grew);
			if (*grew)
				value[r].weight++;
		}
#ifdef TEST_TREES
		printf("Append @%d after:\n",r);
		print();
#endif
	}
	int balance(unsigned r) const {
		return (value[r].right != UINT_MAX ? value[value[r].right].weight : 0)
			   -(value[r].left != UINT_MAX ? value[value[r].left].weight : 0);
	}
	void rebalance(unsigned & r) {
		unsigned x;
		int b;
		
		if (r == UINT_MAX || value[r].weight <= 2)
			return;
#ifdef TEST_TREES
		printf("Rebalance @%d before:\n",r);
		print();
#endif
		b = balance(r);
		if (b < -1) {
			x = value[r].left;
			while (value[x].right != UINT_MAX) {
				rotate_left(x);
				value[r].left = x;
			}
			rotate_right(r);
		} else if (b > 1) {
			x = value[r].right;
			while (value[x].left != UINT_MAX) {
				rotate_right(x);
				value[r].right = x;
			}
			rotate_left(r);
		}
		rebalance(value[r].right);
		rebalance(value[r].left);
#ifdef TEST_TREES
		printf("Rebalance @%d after:\n",r);
		print();
#endif
	}
	void rotate_left(unsigned & r) {
		unsigned x = value[r].right;
		value[r].weight -= value[x].weight;
		value[r].right = value[x].left;
		value[r].weight += (value[x].left != UINT_MAX ? value[value[x].left].weight : 0);
		value[x].weight -= (value[x].left != UINT_MAX ? value[value[x].left].weight : 0);
		value[x].left = r;
		value[x].weight += value[r].weight;
		r = x;
	}
	void rotate_right(unsigned & r) {
		unsigned x = value[r].left;
		value[r].weight -= value[x].weight;
		value[r].left = value[x].right;
		value[r].weight += (value[x].right != UINT_MAX ? value[value[x].right].weight : 0);
		value[x].weight -= (value[x].right != UINT_MAX ? value[value[x].right].weight : 0);
		value[x].right = r;
		value[x].weight += value[r].weight;
		r = x;
	}
#ifdef TEST_TREES
	void print(int level, unsigned r) {
		if (value[r].right != UINT_MAX)
			print(level+1,value[r].right);
		printf("%2d (%2i %+i): ",r,value[r].weight,balance(r));
		for (int i = 0; i < level; i++)
			printf(" ");
		printf("%2d: %2.1f\n",value[r]._v.i,value[r]._v.d);
		if (value[r].left != UINT_MAX)
			print(level+1,value[r].left);
	}
#endif
};

/*
 * class ktree_llrb
 *
 * A key/value tree, like ksset<K,V>, using the left-leaning red/black
 * implementation since it is relatively simple.  The common code between
 * this any the AVL tree needs to be abstracted out into a common base
 * class.
 *
 * This uses index values rather than pointers, and a fixed storage array
 * underneath, to minimize the use of small chunks of memory.  The key class
 * K must implement compare().
 *
 * http://www.cs.princeton.edu/~rs/talks/LLRB/RedBlack.pdf
 */
template <class K, class V>
class ktree_llrb {
public:
	// Make one...
	inline explicit ktree_llrb()
	  : size(0), buffer(0), block(DFLT_BLK), root(UINT_MAX), value(nullptr) {
	}
	// Clean up.
	inline ~ktree_llrb() {
		if (value)
			free(value);
	}
	// The length. Nice for lots of things...
	inline unsigned length() const {
		return size;
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned haskey(const K & k) const {
		unsigned x = root;
		while (x != UINT_MAX) {
			int cmp = k.compare(value[x]._v);
			if (cmp == 0)
				break;
			else if (cmp < 0)
				x = value[x].left;
			else
				x = value[x].right;
		}
		return x;
	}
	// We can use node numbers to find our nodes...
	inline V & operator[] (const unsigned p) const {
		return value[p]._v;
	}
	// Allow sorted access.
	inline V & getindex(const unsigned i) const {
		unsigned x = root, order = 0;
		while (x != UINT_MAX) {
			if (i == order + value[x].order)
				break;
			else if (i < order + value[x].order)
				x = value[x].left;
			else {
				order += value[x].order + 1;
				x = value[x].right;
			}
		}
		// XXX throw if x == UINT_MAX
		return value[x]._v;
	}
	// Get the position of an element in the sort.
	inline unsigned getorder(const unsigned i) const {
		const K & k = static_cast<K &>(value[i]._v);
		unsigned x = root, order = 0;
		while (x != UINT_MAX) {
			int cmp = k.compare(value[x]._v);
			if (cmp == 0) {
				order += value[x].order;
				break;
			} else if (cmp < 0) {
				x = value[x].left;
			} else {
				order += value[x].order + 1;
				x = value[x].right;
			}
		}
		return order;
	}
	// Add node.
	inline void add(const V & v) {
		bool grew = false;
		
		allocate(size+1);
		append(root,v,&grew);
		allocate(size);
		value[root].red = false;
	}
#ifdef TEST_TREES
	void print() {
		if (root == UINT_MAX)
			printf("Empty!\n");
		else
			print(0,0,root);
	}
#endif

protected:
	unsigned size;             // The size of the set...
	unsigned buffer;           // The allocated buffer size...
	unsigned block;            // The minimum block size.
	unsigned root;             // Root of red-black tree.
	struct _V {
	private:
		V _v;
		
		friend class ktree_llrb;

		unsigned order;        // Number of nodes on left
		unsigned left, right;  // Left and right node numbers
		bool red;              // Colour of link to parent

		explicit _V(const V & v)
		  : _v(v), 
		    order(0), left(UINT_MAX), right(UINT_MAX), red(true) {
		}
		// Placement new to support inplace init in the list.
		void * operator new(size_t, void * p) {
			return p;
		}
		void operator delete(void *, void *) {
		}
	} * value;            // Take a guess...

	// Make some space...
	void allocate(const unsigned s) {
		while (s > 8*block)
			block *= 8;
		unsigned b = block*(s/block+(s%block?1:0));
		if (b == buffer)
			return;
		_V * temp = static_cast<_V *>(realloc(value,
				b*sizeof(_V)));
		if (temp == nullptr)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void append(unsigned & r, const V & v, bool * grew) {
		unsigned x;
		
		// If we're UINT_MAX that means we need to make a new node...
		if (r == UINT_MAX) {
			new(&value[size]) _V(v);
			r = size++;
			*grew = true;
			return;
		}
		// Split double reds on the way down.
		x = value[r].left;
		if (x != UINT_MAX && value[x].red) {
			if (value[x].left != UINT_MAX && value[value[x].left].red) {
				value[r].left = value[x].right;
				value[x].right = r;
				value[r].order -= value[x].order + 1;
				r = x;
				value[value[r].left].red = false;
			}
		}
		int cmp = static_cast<const K &>(v).compare(value[r]._v);
		if (cmp == 0) {
			value[r]._v = v;
			return;
		} else if (cmp < 0) {
			// Re-root insert into left tree.
			append(value[r].left,v,grew);
			if (*grew)
				value[r].order++;
		} else {
			// Or right tree.  Same issue here as above.
			append(value[r].right,v,grew);
		}
		x = value[r].right;
		if (x != UINT_MAX && value[x].red) {
			value[r].right = value[x].left;
			value[x].left = r;
			value[x].order += value[r].order + 1;
			r = x;
			value[r].red = value[value[r].left].red;
			value[value[r].left].red = true;
		}
	}
	unsigned remove(unsigned r, const K & k) {
		int cmp = k.compare(value[r]._v);
		if (cmp < 0) {
			// If k is missing do nothing.
			if (value[r].left == UINT_MAX)
				return r;
			// move red left is needed.
			if (!value[value[r].left].red
			 && (value[value[r].left].left == UINT_MAX
			  || !value[value[value[r].left].left].red)) {
				value[r].red = false;
				value[value[r].left].red = true;
				if (value[value[r].right].left != UINT_MAX
				 && value[value[value[r].right].left].red) {
					unsigned x = value[value[r].right].left;
					value[r].left = value[x].right;
					value[x].right = r;
					value[r].order -= value[x].order + 1;
					value[r].right = value[x].left;
					value[x].left = r;
					value[x].order += value[r].order + 1;
					r = x;
				} else
					value[value[r].right].red = true;
			}
			// Re-root delete in left branch.
			value[r].left = remove(value[r].left,k);
		} else {
			// Lean red links right going down.
			if (value[value[r].left != UINT_MAX && value[r].left].red) {
				unsigned x = value[r].left;
				value[r].left = value[x].right;
				value[x].right = r;
				value[r].order -= value[x].order + 1;
				r = x;
				value[r].red = value[value[r].right].red;
				value[value[r].right].red = true;
			}
			// We've found our node, and it's a leaf.
			if (cmp == 0 && value[r].right == UINT_MAX) {
				value[r].~_V();
				return UINT_MAX;
			}
			// Push red right going down.
			if (!value[value[r].right].red
			 && (value[value[r].right->left == UINT_MAX
			  || !value[r].right->left].red)) {
				value[r].red = false;
				value[value[r].right].red = true;
				if (value[value[r].left->left != UINT_MAX
				 && value[r].left->left].red) {
					unsigned x = value[r].left;
					value[r].left = value[x].right;
					value[x].right = r;
					value[r].order -= value[x].order + 1;
					r = x;
					value[r].red = true;
					value[value[r].left].red = false;
				} else
					value[value[r].left].red = true;
			}
			// We've found our node, but it's not a leaf...
			// Replace it with the minimum right node and
			// then remove that node...
			if (cmp == 0) {
				unsigned x = value[r].right;
				while (value[x].left != UINT_MAX)
					x = value[x].left;
				value[r]._v = value[x]._v;
				value[r].right = remove(value[r].right,value[x]._v);
			} else
				// Look on the right.
				value[r].right = remove(value[r].right,k);
		}
		// fixed up right leaning reds on the way up.
		if (value[r].right != 0 && value[r].right->red) {
			unsigned x = value[r].right;
			value[r].right = value[x].left;
			value[x].left = r;
			value[x].order += value[r].order + 1;
			r = x;
			value[r].red = value[value[r].left].red;
			value[value[r].left].red = true;
		}
		return r;
	}
#ifdef TEST_TREES
	void print(int level, unsigned order, unsigned r) {
		if (value[r].right != UINT_MAX)
			print(level+1,order+value[r].order+1,value[r].right);
		printf("%d: (%d %d %d) ",r,order,value[r].order,order+value[r].order);
		for (int i = 0; i < level; i++)
			printf(" ");
		printf("%02d: %f %s\n",value[r]._v.i,value[r]._v.d,(value[r].red?"RED  ":"BLACK"));
		if (value[r].left != UINT_MAX)
			print(level+1,order,value[r].left);
	}
#endif
};

#endif // TREE_H
