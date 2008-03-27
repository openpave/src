/**************************************************************************

	TREE.H - Binary search tree classes

	$OpenPave$

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
		This header implements a templated C++ classes for left leaning
		red-black trees, which are a special variant of binary search
		trees.

	History:
		2008/03/27 - Created a basic implementation.

**************************************************************************/

#ifndef __TREE_H
#define __TREE_H

#include <stdio.h>

/*
 * class tree
 *
 * http://www.cs.princeton.edu/~rs/talks/LLRB/RedBlack.pdf
 */
template<class K, class V>
struct BST
{
	BST()
	  : root(0) {
	}
	~BST() {
		while (root != 0)
			remove(root->key);
	}
	V * get(const K & k) {
		node * x = root;
		while (x != 0) {
			if (k == x->key)
				return &(x->value);
			else if (k < x->key)
				x = x->left;
			else
				x = x->right;
		}
		return 0;
	}
	void insert(const K & k, const V & v) {
		root = insert(root,k,v);
		root->red = false;
	}
	void remove(const K & k) {
		if (root == 0)
			return;
		root = remove(root,k);
		if (root != 0)
			root->red = false;
	}
	void print() {
		printf("Dumping:\n");
		if (root == 0)
			printf("Empty!\n");
		else
			print(0,root);
		printf("\n");
	}

private:
	struct node {
		K key;
		V value;
		node * left, * right;
		bool red;
		node(const K & k, const V & v)
		  : key(k), value(v), red(true) {
		}
		~node() {
		}
	} * root;

	node * insert(node * h, const K & k, const V & v) {
		// If we're zero that means we need to make a new node...
		if (h == 0)
			return new node(k,v);
		// Split double reds on the way down.
		if (h->left != 0 && h->left->red) {
			if (h->left->left != 0 && h->left->left->red) {
				h = rotR(h);
				h->left->red = false;
			}
		}
		if (k == h->key)
			// Key already exists, so replace.
			h->value = v;
		else if (k < h->key)
			// Re-root insert into left tree.
			h->left = insert(h->left,k,v);
		else
			// Or right tree.
			h->right = insert(h->right,k,v);
		if (h->right != 0 && h->right->red)
			h = leanLeft(h);
		return h;
	}
	node * rotL(node * h) {
		node * x = h->right;
		h->right = x->left;
		x->left = h;
		return x;
	}
	node * rotR(node * h) {
		node * x = h->left;
		h->left = x->right;
		x->right = h;
		return x;
	}
	node * leanLeft(node * h) {
		h = rotL(h);
		h->red = h->left->red;
		h->left->red = true;
		return h;
	}
	node * remove(node * h, const K & k) {
		if (k < h->key) {
			// If K is missing do nothing.
			if (h->left == 0)
				return h;
			// move red left is needed.
			if (!h->left->red
			 && (h->left->left == 0 || !h->left->left->red)) {
				h->red = false;
				h->left->red = true;
				if (h->right->left != 0 && h->right->left->red) {
					h->right = rotR(h->right);
					h = rotL(h);
				} else
					h->right->red = true;
			}
			// Re-root delete in left branch.
			h->left = remove(h->left,k);
		} else {
			// Lean red links right going down.
			if (h->left != 0 && h->left->red) {
				h = rotR(h);
				h->red = h->right->red;
				h->right->red = true;
			}
			// We've found our node, and it's a leaf.
			if (k == h->key && h->right == 0) {
				delete h;
				return 0;
			}
			// Push red right going down.
			if (!h->right->red
			 && (h->right->left == 0 || !h->right->left->red)) {
				h->red = false;
				h->right->red = true;
				if (h->left->left != 0 && h->left->left->red) {
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
				while (x->left != 0)
					x = x->left;
				h->key = x->key;
				h->value = x->value;
				h->right = remove(h->right,x->key);
			} else
				// Look on the right.
				h->right = remove(h->right,k);
		}
		// fixed up right leaning reds on the way up.
		if (h->right != 0 && h->right->red)
			h = leanLeft(h);
		return h;
	}	
	void print(int level, node * h) {
		if (h->right != 0)
			print(level+1,h->right);
		for (int i = 0; i < level; i++)
			printf(" ");
		printf("%d: %f %s\n",h->key.i,h->value.d,(h->red?"RED  ":"BLACK"));
		if (h->left != 0)
			print(level+1,h->left);
	}
	
};

#endif // TREE_H
