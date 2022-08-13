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

	Portions Copyright (C) 2006-2022 OpenPave.org.

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

#if !defined(DFLT_BLK)
#define DFLT_BLK    64
#endif

#include <cstring>
#include <stdexcept>
#ifdef TEST_TREES
#include <cstdio>
#endif
#include "mathplus.h"
#include "hascompare.h"

namespace OP {

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
	BST() noexcept
	  : root(nullptr) {
	}
	BST(const BST &) = delete;
	BST & operator = (const BST &) = delete;
	~BST() {
		delete root;
	}
	V * get(const K & k) noexcept {
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
	void remove(const K & k) noexcept {
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
		node(const K & k, const V & v) noexcept
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
	node * rotL(node * h) noexcept {
		node * x = h->right;
		h->right = x->left;
		x->left = h;
		x->order += h->order + 1;
		return x;
	}
	node * rotR(node * h) noexcept {
		node * x = h->left;
		h->left = x->right;
		x->right = h;
		h->order -= x->order + 1;
		return x;
	}
	node * leanLeft(node * h) noexcept {
		h = rotL(h);
		h->red = h->left->red;
		h->left->red = true;
		return h;
	}
	node * leanRight(node * h) noexcept {
		h = rotR(h);
		h->red = h->right->red;
		h->right->red = true;
		return h;
	}
	node * remove(node * h, const K & k) noexcept {
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
 * class tree - Tree base class
 *
 * This class provides the basic working for the rest of the tree classes,
 * including the calculation of the buffer size.
 *
 * This uses index values rather than pointers, and a fixed storage array
 * underneath, to minimize the use of small chunks of memory.  The key class
 * K must implement compare().
 */
template<typename K, typename V, template<typename,typename> class N>
class tree
{
protected:
	template<typename KK, typename VV>
	struct merged_v
	{
		using type = typename std::conditional<
		std::is_base_of<KK,VV>::value || std::is_same<KK,VV>::value
			,std::true_type,std::false_type>::type;
		static constexpr bool value = type::value;
	};

public:
	// The length. Nice for lots of things...
	unsigned length() const noexcept {
		return size;
	}
	// Check if an index is within the set...
	bool inbounds(unsigned p) const noexcept {
		return (p < size ? true : false);
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	bool haskey(const K & k) const noexcept {
		return (getposition(k) != UINT_MAX);
	}
	// Allow acces to the keys.
	const K & getkey(unsigned p) const {
		if (!inbounds(p))
			throw std::out_of_range("key index out of bounds!");
		return value[p];
	}
	// We can use node numbers to find our nodes...
	V & operator [] (unsigned p) const {
		return getatposition(p);
	}
	V & getatposition(unsigned p) const {
		if (!inbounds(p))
			throw std::out_of_range("unordered index out of bounds!");
		return value[p];
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	unsigned getposition(const K & k) const noexcept {
		unsigned x = root;
		while (x != UINT_MAX) {
			int cmp = value[x].compare(k);
			if (cmp == 0)
				break;
			else if (cmp < 0)
				x = value[x].left;
			else
				x = value[x].right;
		}
		return x;
	}
	V & get(const K & k) const {
		return getatposition(getposition(k));
	}
	// Allow sorted access by returning the in-order position in the tree.
	V & getatorder(unsigned i) const {
		unsigned x = root, o = 0;
		while (x != UINT_MAX) {
			if (i == o + order(x))
				break;
			else if (i < o + order(x))
				x = value[x].left;
			else {
				o += order(x) + 1;
				x = value[x].right;
			}
		}
		if (x == UINT_MAX)
			throw std::out_of_range("ordered index out of bounds!");
		return value[x];
	}
	// Get the ordered position of an element in the sort.
	unsigned getorderof(const K & k) const noexcept {
		unsigned x = root, o = 0;
		while (x != UINT_MAX) {
			int cmp = value[x].compare(k);
			if (cmp == 0) {
				o += order(x);
				break;
			} else if (cmp < 0) {
				x = value[x].left;
			} else {
				o += order(x) + 1;
				x = value[x].right;
			}
		}
		if (x == UINT_MAX)
			o = UINT_MAX;
		return o;
	}
	void empty() {
		allocate(0);
	}

protected:
	unsigned size = 0;         // The size of the tree...
	unsigned buffer = 0;       // The allocated buffer size...
	unsigned root = UINT_MAX;  // Root of the tree.
	struct _node_base : N<K,V>::node_base {
		_node_base() noexcept
		  : N<K,V>::node_base(), left(UINT_MAX), right(UINT_MAX) {
		}
		unsigned weight = 1;   // Number of nodes below (including this)
		unsigned left, right;  // Left and right node numbers
	};
	template<typename KK, typename VV, typename = void>
	struct _node_val {
		V _v;
		explicit _node_val(const V & v) noexcept
		  : _v(v) {
		}
		explicit _node_val(V && v) noexcept
		  : _v(std::move(v)) {
		}
		operator const V & () const noexcept {
			return _v;
		}
		operator V & () noexcept {
			return _v;
		}
	};
	template<typename KK, typename VV>
	struct _node_val<KK,VV,typename std::enable_if<merged_v<KK,VV>::value
			&& !std::is_same<KK,VV>::value>::type> {
		V _v;
		explicit _node_val(const V & v) noexcept
		  : _v(v) {
		}
		explicit _node_val(V && v) noexcept
		  : _v(std::move(v)) {
		}
		operator const K & () const noexcept {
			return _v;
		}
		operator const V & () const noexcept {
			return _v;
		}
		operator V & () noexcept {
			return _v;
		}
	};
	template<typename KK, typename VV>
	struct _node_val<KK,VV,typename std::enable_if<
			!merged_v<KK,VV>::value>::type> {
		K _k;
		V _v;
		explicit _node_val(const K & k, const V & v)
		  : _k(k), _v(v) {
		}
		explicit _node_val(K && k, V && v) noexcept
		  : _k(std::move(k)), _v(std::move(v)) {
		}
		operator const K & () const noexcept {
			return _k;
		}
		operator const V & () const noexcept {
			return _v;
		}
		operator V & () noexcept {
			return _v;
		}
	};
	struct _V : _node_base, _node_val<K,V> {
		template<typename KK = K, typename VV = V>
		_V(const V & v, typename std::enable_if<
				merged_v<KK,VV>::value>::type * = nullptr) noexcept
		  : _node_base(), _node_val<K,V>(v) {
		}
		template<typename KK = K, typename VV = V>
		_V(const K & k, const V & v, typename std::enable_if<
				!merged_v<KK,VV>::value>::type * = nullptr)
		  : _node_base(), _node_val<K,V>(k,v) {
		}
		template<typename KK = K, typename VV = V>
		_V(V && v, typename std::enable_if<
				merged_v<KK,VV>::value>::type * = nullptr)
		  : _node_base(), _node_val<K,V>(std::move(v)) {
		}
		template<typename KK = K, typename VV = V>
		_V(K && k, V && v, typename std::enable_if<
				!merged_v<KK,VV>::value>::type * = nullptr)
		  : _node_base(), _node_val<K,V>(std::move(k),std::move(v)) {
		}
		_V(const _V & v) noexcept
		  : _node_base(), _node_val<K,V>(v) {
		}
		_V(_V && v) noexcept
		  : _node_base(), _node_val<K,V>(std::move(v)) {
		}
		// Placement new to support in-place initialization in the list.
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void *, void *) noexcept {
		}
		// note these are reversed for historical reasons.
		// use the compare function if it has one
		template<typename T = K>
		typename std::enable_if<has_compare<T>::value,int>::type
		compare(const K & k) const noexcept {
			return k.compare(*this);
		}
		// else resort to operators
		template<typename T = K>
		typename std::enable_if<!has_compare<T>::value,int>::type
		compare(const K & k) const noexcept {
			return (*this < k ? 1 :	*this == k ? 0 : -1);
		}
	};
	struct _V * value = nullptr;       // Take a guess...

	// Simple constructor...
	explicit tree(unsigned b = DFLT_BLK) noexcept
	  : block(b > 1 ? b : 1) {
	}
	// Copy constructor.
	explicit tree(const tree & t)
	  : block(t.block) {
		allocate(t.size);
		// The other real tree is responsible for the copy.
	}
	// Move constructor.
	explicit tree(tree && t) noexcept
	  : size(t.size), buffer(t.buffer), root(t.root), value(t.value),
		block(t.block) {
		t.size = 0;
		t.buffer = 0;
		t.root = UINT_MAX;
		t.value = nullptr;
	}
	tree & operator = (const tree &) = delete;
	tree & operator = (tree && t) noexcept {
		std::swap(size,t.size);
		std::swap(buffer,t.buffer);
		std::swap(root,t.root);
		std::swap(value,t.value);
		std::swap(block,t.block);
		return *this;
	};
	~tree() {
		allocate(0);
	}
	// Make some space...
	void allocate(unsigned s) {
		unsigned b = bufsize(s);
		if (b == buffer)
			return;
		if (b == 0) {
			for (unsigned i = 0; i < size; i++)
				value[i].~_V();
			free(value);
			value = nullptr;
			buffer = 0;
			size = 0;
			root = UINT_MAX;
			return;
		}
		_V * temp = static_cast<_V *>(realloc(value,b*sizeof(_V)));
		if (temp == nullptr)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	// Insert an element
	unsigned insert(const _V & v) noexcept {
		new(&value[size]) _V(v);
		return size++;
	}
	unsigned insert(_V && v) noexcept {
		new(&value[size]) _V(std::move(v));
		return size++;
	}
	void expunge(unsigned p) noexcept {
		// Do the ugly work of removing the element and adjusting the
		// offsets into value[].  This would be better if deferred to a
		// compact() function of some form and used the nodes to make a
		// free list.
		value[p].~_V();
		if (--size > p)
			std::memmove(&value[p],&value[p+1],(size-p)*sizeof(_V));
		if (root != UINT_MAX && root > p)
			root--;
		for (unsigned i = 0; i < size; i++) {
			if (value[i].left != UINT_MAX && value[i].left > p)
				value[i].left--;
			if (value[i].right != UINT_MAX && value[i].right > p)
				value[i].right--;
		}
	}
	// This is the number of nodes in this sub-tree including the root
	unsigned weight(unsigned r) const noexcept {
		return (r != UINT_MAX ? value[r].weight : 0);
	}
	// This is the number of nodes in the left sub-tree
	unsigned order(unsigned r) const noexcept {
		return (r != UINT_MAX ? weight(value[r].left) : 0);
	}
	// Get the new weight if we have disturbed the tree
	unsigned new_weight(unsigned r) const noexcept {
		return weight(value[r].left)+weight(value[r].right)+1;
	}

private:
	unsigned block;            // The minimum block size.
	// Calculate the buffer size.
	unsigned bufsize(unsigned s) noexcept {
		while (s > 8*block)
			block *= 8;
		//while (64*s < block)
		//	block /= 8;
		return block*(s/block+(s%block?1:0));
	}
};

/*
 * class ktree_avl
 *
 * A key/value tree, like ksset<K,V>, using the AVL implementation since
 * it is relatively simple.
 *
 * See https://en.wikipedia.org/wiki/AVL_tree
 */
template <typename K, typename V = K>
class ktree_avl : public tree<K,V,ktree_avl>
{
protected:
	template<typename KK, typename VV>
	using merged_v = typename tree<K,V,OP::ktree_avl>::template merged_v<KK,VV>;
	using _V = typename tree<K,V,OP::ktree_avl>::_V;
	using tree<K,V,OP::ktree_avl>::value;
	using tree<K,V,OP::ktree_avl>::root;
	using tree<K,V,OP::ktree_avl>::size;
	using tree<K,V,OP::ktree_avl>::allocate;
	using tree<K,V,OP::ktree_avl>::insert;
	using tree<K,V,OP::ktree_avl>::expunge;
	using tree<K,V,OP::ktree_avl>::weight;
	using tree<K,V,OP::ktree_avl>::new_weight;

public:
	// Make one...
	ktree_avl(unsigned b = DFLT_BLK) noexcept
	  : tree<K,V,OP::ktree_avl>(b) {
	}
	// Copy constructor.
	explicit ktree_avl(const ktree_avl & t)
	  : tree<K,V,OP::ktree_avl>(t) {
		unsigned p = UINT_MAX;
		for (unsigned i = 0; i < t.size; i++)
			append(root,_V(t.value[i]),&p);
	}
	// Move constructor.
	explicit ktree_avl(ktree_avl && t) noexcept
	  : tree<K,V,OP::ktree_avl>(std::move(t)) {
	}
	ktree_avl & operator = (const ktree_avl & t) {
		allocate(0);
		allocate(t.size);
		unsigned p = UINT_MAX;
		for (unsigned i = 0; i < t.size; i++)
			append(root,_V(t.value[i]),&p);
	}
	ktree_avl & operator = (ktree_avl && t) noexcept {
		tree<K,V,OP::ktree_avl>::operator=(std::move(t));
		return *this;
	}
	// Clean up
	~ktree_avl() {
	}
	// Add a node.  Returns the new position in value (usually == size)
	template<typename KK = K, typename VV = V,
		typename = typename std::enable_if<merged_v<KK,VV>::value>::type>
	void add(const V & v) {
		unsigned p = UINT_MAX;
		allocate(size+1);
		append(root,_V(v),&p);
		allocate(size);
#ifdef TEST_TREES
		assert_avl();
#endif
	}
	template<typename KK = K, typename VV = V,
		typename = typename std::enable_if<!merged_v<KK,VV>::value>::type>
	void add(const K & k, const V & v) {
		unsigned p = UINT_MAX;
		allocate(size+1);
		append(root,_V(k,v),&p);
		allocate(size);
#ifdef TEST_TREES
		assert_avl();
#endif
	}
	// Add a node.  Returns the new position in value (usually == size)
	template<typename KK = K, typename VV = V,
		typename = typename std::enable_if<merged_v<KK,VV>::value>::type>
	void add(V && v) {
		unsigned p = UINT_MAX;
		allocate(size+1);
		append(root,_V(std::move(v)),&p);
		allocate(size);
#ifdef TEST_TREES
		assert_avl();
#endif
	}
	template<typename KK = K, typename VV = V,
		typename = typename std::enable_if<!merged_v<KK,VV>::value>::type>
	void add(K && k, V && v) {
		unsigned p = UINT_MAX;
		allocate(size+1);
		append(root,_V(std::move(k),std::move(v)),&p);
		allocate(size);
#ifdef TEST_TREES
		assert_avl();
#endif
	}
	// Remove a node.  Returns the old position in value
	void remove(const K & k) {
		unsigned p = UINT_MAX;
		remove(root,k,&p);
		if (p == UINT_MAX)
			return;
		expunge(p);
		allocate(size);
#ifdef TEST_TREES
		assert_avl();
#endif
	}
#ifdef TEST_TREES
	void assert_avl() {
		int h;
		unsigned s = 0;
		assert_avl(root,&h,&s);
		if (s != size)
			throw std::runtime_error("oops");
	}
	void print() {
		if (root == UINT_MAX)
			printf("Empty!\n");
		else
			print(0,root);
	}
#endif

private:
	friend class tree<K,V,OP::ktree_avl>;

	// For AVL trees the only additional node information is the height
	struct node_base {
		int height = 1;
	};

	int height(unsigned r) const noexcept {
		return (r != UINT_MAX ? value[r].height : 0);
	}
	int balance(unsigned r) const noexcept {
		return height(value[r].right)-height(value[r].left);
	}
	int new_height(unsigned r) const noexcept {
		return MAX(height(value[r].left), height(value[r].right)) + 1;
	}
	void append(unsigned & r, _V && v, unsigned * p) {
		// If we're UINT_MAX that means we need to make a new node...
		if (r == UINT_MAX) {
			*p = r = insert(std::move(v));
			return;
		}
#ifdef TEST_TREES
		printf("Append @%d before:\n",r);
		print();
#endif
		int cmp = value[r].compare(v);
		if (cmp == 0) {
			value[*p = r]._v = std::move(v._v);
			return;
		} else if (cmp < 0)
			append(value[r].left,std::move(v),p);
		else
			append(value[r].right,std::move(v),p);
		if (*p != UINT_MAX) {
			value[r].weight = new_weight(r);
			value[r].height = new_height(r);
		}
		rebalance(r);
#ifdef TEST_TREES
		printf("Append @%d after:\n",r);
		print();
#endif
	}
	void remove(unsigned & r, const K & k, unsigned * p) {
		if (r == UINT_MAX)
			return; // XXX throw?
#ifdef TEST_TREES
		printf("Remove @%d before:\n",r);
		print();
#endif
		int cmp = value[r].compare(k);
		if (cmp == 0) {
			// We've found our node, and it's a leaf.
			if (value[r].right == UINT_MAX) {
				*p = r, r = value[r].left, value[*p].left = UINT_MAX;
				return;
			}
			// We've found our node, but it's not a leaf. Swap it with
			// the minimum right node and then remove that node...
			unsigned x = value[r].right;
			while (value[x].left != UINT_MAX)
				x = value[x].left;
			swap(value[r]._v,value[x]._v);
			remove(value[r].right,k,p);
		} else if (cmp < 0)
			remove(value[r].left,k,p);
		else
			remove(value[r].right,k,p);
		if (*p != UINT_MAX) {
			value[r].weight = new_weight(r);
			value[r].height = new_height(r);
		}
#ifdef TEST_TREES
		printf("Remove @%d after:\n",r);
		print();
#endif
		rebalance(r);
		return;
	}
	void rebalance(unsigned & r) noexcept {
		if (r == UINT_MAX || value[r].height < 2)
			return;
#ifdef TEST_TREES
		printf("Rebalance @%d before:\n",r);
		print();
#endif
		int b = balance(r);
		if (b < -1) {
			if (balance(value[r].left) > 0)
				rotate_left(value[r].left);
			rotate_right(r);
		} else if (b > 1) {
			if (balance(value[r].right) < 0)
				rotate_right(value[r].right);
			rotate_left(r);
		}
#ifdef TEST_TREES
		printf("Rebalance @%d after:\n",r);
		print();
		int h;
		unsigned s = 0;
		assert_avl(r,&h,&s);
#endif
		return;
	}
	void rotate_left(unsigned & r) noexcept {
		unsigned x = value[r].right;
		value[r].right = value[x].left;
		value[r].weight -= weight(x)-weight(value[x].left);
		value[r].height = new_height(r);
		value[x].weight += weight(r)-weight(value[x].left);
		value[x].left = r;
		value[x].height = new_height(x);
		r = x;
	}
	void rotate_right(unsigned & r) noexcept {
		unsigned x = value[r].left;
		value[r].left = value[x].right;
		value[r].weight -= weight(x)-weight(value[x].right);
		value[r].height = new_height(r);
		value[x].weight += weight(r)-weight(value[x].right);
		value[x].right = r;
		value[x].height = new_height(x);
		r = x;
	}
#ifdef TEST_TREES
	void assert_avl(unsigned r, int * h, unsigned * s) {
		int lh, rh;
		unsigned ls = 0, rs = 0;
		if (r == UINT_MAX) {
			*h = 0;
			return;
		}
		assert_avl(value[r].left,&lh,&ls);
		assert_avl(value[r].right,&rh,&rs);
		*h = MAX(lh,rh)+1;
		*s += ls+rs+1;
		if (value[r].height != new_height(r))
			throw std::runtime_error("oops");
		if (value[r].height != *h)
			throw std::runtime_error("oops");
		if (abs(balance(r)) > 1)
			throw std::runtime_error("oops");
		if (balance(r) != (rh-lh))
			throw std::runtime_error("oops");
		if (weight(r) != (ls+rs+1))
			throw std::runtime_error("oops");
	}
	void print(int level, unsigned r) {
		if (value[r].right != UINT_MAX)
			print(level+1,value[r].right);
		printf("%2d (%2i %2i %+i): ",r,weight(r),height(r),balance(r));
		for (int i = 0; i < level; i++)
			printf(" ");
		printf("%2d: %2.1f %s\n",value[r]._v.i,value[r]._v.d,
			(value[r].height != new_height(r) ? "***" : ""));
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
template<class K, class V>
class ktree_llrb {
public:
	// Make one...
	ktree_llrb() noexcept
	  : size(0), buffer(0), block(DFLT_BLK), root(UINT_MAX), value(nullptr) {
	}
	// Clean up.
	~ktree_llrb() {
		if (value)
			free(value);
	}
	// The length. Nice for lots of things...
	unsigned length() const noexcept {
		return size;
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	unsigned haskey(const K & k) const noexcept {
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
	V & operator [] (unsigned p) const {
		return value[p]._v;
	}
	// Allow sorted access.
	V & getindex(unsigned i) const noexcept {
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
	unsigned getorder(unsigned i) const noexcept {
		K & k = static_cast<K &>(value[i]._v);
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
	void add(const V & v) {
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
		V _v;
		unsigned order;        // Number of nodes on left
		unsigned left, right;  // Left and right node numbers
		bool red;              // Colour of link to parent
		explicit _V(const V & v) noexcept
		  : _v(v),
		    order(0), left(UINT_MAX), right(UINT_MAX), red(true) {
		}
		// Placement new to support inplace init in the list.
		void * operator new (size_t, void * p) noexcept {
			return p;
		}
		void operator delete (void *, void *) noexcept {
		}
	} * value;            // Take a guess...

	// Make some space...
	void allocate(unsigned s) {
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
	void append(unsigned & r, const V & v, bool * grew) noexcept {
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
	unsigned remove(unsigned r, const K & k) noexcept {
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

} // namespace OP

#endif // TREE_H
