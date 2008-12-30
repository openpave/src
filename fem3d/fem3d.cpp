/**************************************************************************

	FEM3D.CPP - A special 3D Finite element code for my PhD.

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

	Portions Copyright (C) 2007 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Design:
		This is a 3D Finite Element Method (FEM) code, which is designed
		for pavements problems.  Specifically, it is designed for my
		thesis work, which requires fast linear elastic 3D FEM analysis. 
		The code is fast.  It does most of the maths using expression
		templates, and to aid that uses small 3x3 matrices to store the
		DOF level data.  This works because the math for the FEM is really
		a tensor based math, with the matrix and vectors being a
		simplification.  In essence, this code solves it as two nested
		matrices, to make a forth order tensor.  This also has an
		advantage in the matrix storage, because we can manage less.
		
		The code follows the classic structure of finite elements for
		linear elasticity.  The solution is built of elements, which are
		connected to a global nodes.  These elements are integrated to
		obtain a element stiffness matrix, which is done by doing
		gaussian integration over the point stifness matrix, with the
		approriate transforms.  The element stiffness matrix is then
		assembled into a global stifness matrix.  A nodal force vector is
		assembled, and the bondary conditions are applied, and then the
		global stiffness matrix is inverted to find the nodal
		displacements.  These are then used to get the final results.
		
		The inversion is not done directly, but via a conjugate gradient
		itterative solver, with an incomplete Cholesky decompostion as a
		preconditioner.

	History:
		2007/09/07 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#define _EVENT_IMP
#define _PROGRESS_IMP

#include <stdlib.h>
#if defined(_MSC_VER)
#include <malloc.h>
#define alloca _alloca
#endif
#include "event.h"
#include "fixed.h"
#include "set.h"
#include "tmatrix.h"
#include "matrix.h"
#include "rng.h"
#include "pavement.h"

#if defined(__FreeBSD__)
#include <stdlib.h>
extern "C" {
extern const char * _malloc_options = "ajz";
};
#endif

#define NDOF 3
#define NDIM 3

/*
 * These are the block matrices used to store the DOF data.  These
 * typedefs just help to make the code more readable.
 */
typedef tmatrix<double,NDOF,NDOF> smatrix_dof;
typedef tmatrix<double,NDOF,1>    svector_dof;
typedef tmatrix<double,NDIM,NDOF> ematrix;

/*
 * Templated Kronecker delta
 */
template<unsigned I, unsigned J>
static inline double delta()
{
	return (I == J ? 1.0 : 0.0);
}

/*
 * struct meta_stiffness - template meta class to assign the stiffness
 * tensor based on lambda and mu.  DO NOT PLAY WITH THIS.
 */
template<unsigned I = NDIM, unsigned J = NDOF,
	unsigned K = NDIM, unsigned L = NDOF>
struct meta_stiffness
{
	static inline void subassignR(smatrix_dof & E,
			const double & l, const double & m) {
		meta_stiffness<I,J,K,L-1>::subassignR(E,l,m);
		E(J-1,L-1) = l*delta<I-1,J-1>()*delta<K-1,L-1>()
				+ m*(delta<I-1,K-1>()*delta<J-1,L-1>()
				   + delta<I-1,L-1>()*delta<J-1,K-1>());
	}
	static inline void subassign(smatrix_dof & E,
			const double & l, const double & m) {
		meta_stiffness<I,J-1,K,L>::subassign(E,l,m);
		meta_stiffness<I,J,K,L>::subassignR(E,l,m);
	}
	static inline void assignR(smatrix_dof (& E)[NDIM][NDIM],
			const double & l, const double & m) {
		meta_stiffness<I,J,K-1,L>::assignR(E,l,m);
		meta_stiffness<I,J,K,L>::subassign(E[I-1][K-1],l,m);
	}
	static inline void assign(smatrix_dof (& E)[NDIM][NDIM],
			const double & l, const double & m) {
		meta_stiffness<I-1,J,K,L>::assign(E,l,m);
		meta_stiffness<I,J,K,L>::assignR(E,l,m);
	}
	static inline void submulR(ematrix & s, const ematrix & e,
			const double & l, const double & m) {
		meta_stiffness<I,J,K,L-1>::submulR(s,e,l,m);
		s(I-1,J-1) += (l*delta<I-1,J-1>()*delta<K-1,L-1>()
				+ m*(delta<I-1,K-1>()*delta<J-1,L-1>()
				   + delta<I-1,L-1>()*delta<J-1,K-1>()))
				*e(K-1,L-1);
	}
	static inline void submul(ematrix & s, const ematrix & e,
			const double & l, const double & m) {
		meta_stiffness<I,J-1,K,L>::submul(s,e,l,m);
		meta_stiffness<I,J,K,L>::submulR(s,e,l,m);
	}
	static inline void mulR(ematrix & s, const ematrix & e,
			const double & l, const double & m) {
		meta_stiffness<I,J,K-1,L>::mulR(s,e,l,m);
		meta_stiffness<I,J,K,L>::submul(s,e,l,m);
	}
	static inline void mul(ematrix & s, const ematrix & e,
			const double & l, const double & m) {
		meta_stiffness<I-1,J,K,L>::mul(s,e,l,m);
		meta_stiffness<I,J,K,L>::mulR(s,e,l,m);
	}
};
template<unsigned I, unsigned J, unsigned K>
struct meta_stiffness<I,J,K,0>
{
	static inline void subassignR(smatrix_dof &,
			const double &, const double &) {}
	static inline void submulR(ematrix &, const ematrix &,
			const double &, const double &) {}
};
template<unsigned I, unsigned K, unsigned L>
struct meta_stiffness<I,0,K,L>
{
	static inline void subassign(smatrix_dof &,
			const double &, const double &) {}
	static inline void submul(ematrix &, const ematrix &,
			const double &, const double &) {}
};
template<unsigned I, unsigned J, unsigned L>
struct meta_stiffness<I,J,0,L>
{
	static inline void assignR(smatrix_dof (&)[NDIM][NDIM],
			const double &, const double &) {}
	static inline void mulR(ematrix &, const ematrix &,
			const double &, const double &) {}
};
template<unsigned J, unsigned K, unsigned L>
struct meta_stiffness<0,J,K,L>
{
	static inline void assign(smatrix_dof (&)[NDIM][NDIM],
			const double &, const double &) {}
	static inline void mul(ematrix &, const ematrix &,
			const double &, const double &) {}
};

/*
 * Calculate the principle stresses or strains.
 */
ematrix
principle(const ematrix & s)
{
	double t1, t2;

	assert(NDIM == 3 && NDOF == 3); 
	ematrix p(s);

	t1 = p(0,1)*p(0,1) + p(0,2)*p(0,2) + p(1,2)*p(1,2);
	while (t1 > 1e-30) {
		if (p(0,2) != 0.0) {
			t1 = p(2,2)-p(0,0);
			t2 = (t1 < 0.0 ? -2 : 2)*p(0,2);
			t1 = t2/(fabs(t1) + hypot(t1,t2));
			p(0,0) -= t1*p(0,2);
			p(2,2) += t1*p(0,2);
			t2 = 1.0/hypot(t1,1.0);
			p(1,2) *= t2;
			t2 = p(0,1) *= t2;
			p(0,1) -= t1*p(1,2);
			p(1,2) += t1*t2;
			p(0,2) = 0.0;
		}
		if (p(0,1) != 0.0) {
			t1 = p(1,1)-p(0,0);
			t2 = (t1 < 0.0 ? -2 : 2)*p(0,1);
			t1 = t2/(fabs(t1) + hypot(t1,t2));
			p(0,0) -= t1*p(0,1);
			p(1,1) += t1*p(0,1);
			p(1,2) /= hypot(t1,1.0);
			p(0,2) = -t1*p(1,2);
			p(0,1) = 0.0;
		}
		if (p(1,2) != 0.0) {
			t1 = p(2,2)-p(1,1);
			t2 = (t1 < 0.0 ? -2 : 2)*p(1,2);
			t1 = t2/(fabs(t1) + hypot(t1,t2));
			p(1,1) -= t1*p(1,2);
			p(2,2) += t1*p(1,2);
			p(0,2) /= hypot(t1,1.0);
			p(0,1) = -t1*p(0,2);
			p(1,2) = 0.0;
		}
		t1 = p(0,1)*p(0,1) + p(0,2)*p(0,2);
	}
	p(1,0) = p(0,1) = 0.0;
	p(2,0) = p(0,2) = 0.0;
	p(2,1) = p(1,2) = 0.0;
	for (unsigned i = 0; i < 2; i++) {
		for (unsigned j = i + 1; j < 3; j++) {
			if (p(i,i) > p(j,j))
				swap(p(i,i),p(j,j));
		}
	}

	return p;
}

/*
 * class material
 *
 * This class represents a material.  It will eventually grow to
 * support a wide variety of material types.  At the moment it
 * really is overkill...
 */
class material {
public:
	inline material(const double e, const double poissons)
	  : emod(e), v(poissons) {
	}
	// Get the point stiffness tensor.  This should be made to cache
	// the result.
	inline void pointstiffness(smatrix_dof (& E)[NDIM][NDIM]) const {
		double lambda = v*emod/(1+v)/(1-2*v);
		double mu = emod/2/(1+v);
		meta_stiffness<>::assign(E,lambda,mu);
	}
	// Get the point stress tensor, based on the strain.
	inline ematrix pointstress(const ematrix & e) const {
		double lambda = v*emod/(1+v)/(1-2*v);
		double mu = emod/2/(1+v);
		ematrix s(0.0);

		meta_stiffness<>::mul(s,e,lambda,mu);
		return s;
	}

//private:
	double emod;
	double v;
};

/*
 * struct coord3d
 *
 * A mesh point in 3D space.  This is not a point3d for two reasons:
 * 1. It sorts differently.
 * 2. We use a fixed point type to make sure the nodes line up.
 * 3. It has two z coordinates.  A real one and a relative one.
 */
struct coord3d {
	fixed<8> x, y, z;

	inline coord3d() {
	}
	inline coord3d(double px, double py, double pz)
	  : x(px), y(py), z(pz) {
	}
	inline coord3d(const coord3d & p)
	  : x(p.x), y(p.y), z(p.z) {
	}
	inline coord3d(const point3d & p)
	  : x(p.x), y(p.y), z(p.z) {
	}
	inline ~coord3d () {
	}
	int compare(const coord3d & p) const {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if (z > p.z) return -1;
		if (z < p.z) return  1;
		if (x < p.x) return -1;
		if (x > p.x) return  1;
		return (y < p.y ? -1 : 1);
	}
	bool operator == (const coord3d & p) const {
		return (x == p.x && y == p.y && z == p.z ? true : false);
	}
	bool operator != (const coord3d & p) const {
		return (x != p.x || y != p.y || z != p.z ? true : false);
	}
	inline bool operator > (const coord3d & p) const {
		return (compare(p) == 1 ? true : false);
	}
	inline bool operator >= (const coord3d & p) const {
		return (compare(p) != -1 ? true : false);
	}
	inline bool operator < (const coord3d & p) const {
		return (compare(p) == -1 ? true : false);
	}
	inline bool operator <= (const coord3d & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * Gauss Quadrature abscissa and weights.
 */
const double gp_2[2][2] = {{-1/sqrt(3.0), 1},
                           {+1/sqrt(3.0), 1}};
const double gp_3[3][2] = {{-sqrt(3/5.0), 5/9.0},
                           {           0, 8/9.0},
                           {+sqrt(3/5.0), 5/9.0}};
const double gp_4[4][2] = {{-sqrt(525+70*sqrt(30.0))/35, (18-sqrt(30.0))/36},
                           {-sqrt(525-70*sqrt(30.0))/35, (18+sqrt(30.0))/36},
                           {+sqrt(525-70*sqrt(30.0))/35, (18+sqrt(30.0))/36},
                           {+sqrt(525+70*sqrt(30.0))/35, (18-sqrt(30.0))/36}};
const double gp_A[4][2] = {{-(1+1/sqrt(3.0))/2, 1/2.0},
                           {-(1-1/sqrt(3.0))/2, 1/2.0},
                           {+(1-1/sqrt(3.0))/2, 1/2.0},
                           {+(1+1/sqrt(3.0))/2, 1/2.0}};

class mesh;                    // Forward declare.
class node_list;               // Forward declare.
double node_depth_callback(const coord3d &, const mesh *, const material *);
double node_emod_callback(const coord3d &, const mesh *, const material *);

/*
 * struct node3d
 *
 * A node in 3D, based on coord3d.  This stores the meshes neighbours,
 * so it can keep track of points for the variable node elements.
 */
struct node3d : public coord3d {
	// Neighbour in (x,y,z) plus/minus directions.
	// UINT_MAX if not there...
	struct {
		unsigned xm, xp, ym, yp, zm, zp;
	};
	// Final results.
	struct {
		double ux, uy, uz;
	};
	// Set new neighbours.  If this conflicts with what we know, we
	// complain (if we're debugging).
	void setneighbours(const unsigned x_m, const unsigned x_p,
			const unsigned y_m, const unsigned y_p,
			const unsigned z_m, const unsigned z_p) {
		if (x_m != UINT_MAX) {
			assert(xm == UINT_MAX || xm == x_m);
			xm = x_m;
		}
		if (x_p != UINT_MAX) {
			assert(xp == UINT_MAX || xp == x_p);
			xp = x_p;
		}
		if (y_m != UINT_MAX) {
			assert(ym == UINT_MAX || ym == y_m);
			ym = y_m;
		}
		if (y_p != UINT_MAX) {
			assert(yp == UINT_MAX || yp == y_p);
			yp = y_p;
		}
		if (z_m != UINT_MAX) {
			assert(zm == UINT_MAX || zm == z_m);
			zm = z_m;
		}
		if (z_p != UINT_MAX) {
			assert(zp == UINT_MAX || zp == z_p);
			zp = z_p;
		}
	}

private:
	friend class node_list;
	friend class mesh;
	
	unsigned order;            // Number of nodes on left
	unsigned left, right;      // Left and right node numbers
	bool red:1;                // Color of link to parent
	unsigned fixed:31;         // Bitmask of fixed DOFs

	explicit node3d(const coord3d & c)
	  : coord3d(c), xm(UINT_MAX), xp(UINT_MAX), ym(UINT_MAX), yp(UINT_MAX),
			zm(UINT_MAX), zp(UINT_MAX), ux(0.0), uy(0.0), uz(0.0),
			order(0), left(UINT_MAX), right(UINT_MAX), red(true), fixed(0) {
	}
	// Set our results.
	void setdisp(const tmatrix<double,NDIM,1> & t) {
		ux = t(0); uy = t(1); uz = t(2);
	}
	// Placement new to support inplace init in the list.
	void * operator new(size_t, void * p) {
		return p;
	} 
	void operator delete(void *, void *) {
	} 
};

/*
 * class node_list
 *
 * We use special set type to store nodes, since we have some interesting
 * requirements.  This is basically a left-leaning red-black tree, but
 * without using pointers...  It uses a fixed array to store the data
 * and array indices into that to find the nodes.
 *
 * Since we never delete nodes, this helps with some of the data management,
 * since we can skip a bunch of complexity...  We also have a very limited
 * set of functions, because we know our use case exactly.
 */
class node_list {
public:
	// Make one...
	inline explicit node_list()
	  : size(0), buffer(0), block(DFLT_BLK), root(UINT_MAX), value(0) {
	}
	// Clean up.
	inline ~node_list() {
		if (value)
			free(value);
	}
	// The length. Nice for lots of things...
	inline unsigned length() const {
		return size;
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned haskey(const coord3d & k) const {
		unsigned x = root;
		while (x != UINT_MAX) {
			int cmp = k.compare(value[x]);
			if (cmp == 0)
				break;
			else if (cmp < 0)
				x = value[x].left;
			else
				x = value[x].right;
		}
		return x;
	}
	// We use node numbers to find our nodes...
	inline node3d & operator[] (const unsigned p) const {
		return value[p];
	}
	// Allow sorted acess.
	inline node3d & getindex(const unsigned i) const {
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
		return value[x];
	}
	// Get the position of an element in the sort.
	inline unsigned getorder(const unsigned i) const {
		coord3d & k = static_cast<coord3d &>(value[i]);
		unsigned x = root, order = 0;
		while (x != UINT_MAX) {
			int cmp = k.compare(value[x]);
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
	inline bool add(const coord3d & v) {
		if (!allocate(size+1))
			return false;
		if (!append(root,v))
			return false;
		value[root].red = false;
		return true;
	}

protected:
	unsigned size;             // The size of the set...
	unsigned buffer;           // The allocated buffer size...
	unsigned block;            // The minimum block size.
	unsigned root;             // Root of red-black tree.
	node3d * value;            // Take a guess...

	// Make some space...
	bool allocate(const unsigned s) {
		while (s > 8*block)
			block *= 8;
		unsigned b = block*(s/block+(s%block?1:0));
		if (b == buffer)
			return true;
		node3d * temp = static_cast<node3d *>(realloc(value,
				b*sizeof(node3d)));
		if (temp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in node_list::allocate()!");
			return false;
		}
		value = temp;
		buffer = b;
		return true;
	}
	bool append(unsigned & r, const coord3d & c) {
		// If we're UINT_MAX that means we need to make a new node...
		if (r == UINT_MAX) {
			new(&value[size]) node3d(c);
			r = size++;
			return true;
		}
		// Split double reds on the way down.
		if (value[r].left != UINT_MAX && value[value[r].left].red) {
			if (value[value[r].left].left != UINT_MAX
					&& value[value[value[r].left].left].red) {
				unsigned x = value[r].left;
				value[r].left = value[x].right;
				value[x].right = r;
				value[r].order -= value[x].order + 1;
				r = x;
				value[value[r].left].red = false;
			}
		}
		int cmp = c.compare(value[r]);
		if (cmp == 0) {
			assert(false);
			return false;
		} else if (cmp < 0) {
			// Re-root insert into left tree.
			if (!append(value[r].left,c))
				return false;
			value[r].order++;
		} else {
			// Or right tree.  Same issue here as above.
			if (!append(value[r].right,c))
				return false;
		}
		if (value[r].right != UINT_MAX && value[value[r].right].red) {
			unsigned x = value[r].right;
			value[r].right = value[x].left;
			value[x].left = r;
			value[x].order += value[r].order + 1;
			r = x;
			value[r].red = value[value[r].left].red;
			value[value[r].left].red = true;
		}
		return true;
	}
};

/*
 * class smatrix_elem
 *
 * Element stiffness matrix, which is an upper triangular matrix,
 * so we only store one side.  This is only used to transfer the
 * information.
 */
class smatrix_elem {
public:
	explicit smatrix_elem(const unsigned n)
	  : nnd(n), K(0) {
		// calloc() to zero the memory...
		K = static_cast<smatrix_dof *>
				(calloc(size(),sizeof(smatrix_dof)));
		if (K == 0)
			event_msg(EVENT_ERROR,"Out of memory in smatrix_elem::smatrix_elem()!");
	}
	inline ~smatrix_elem() {
		if (K != 0)
			free(K);
	}
	inline smatrix_dof & operator() (unsigned i, unsigned j) const {
		assert(j >= i);
		return K[index(i,j)];
	}

private:
	friend class mesh;

	const unsigned nnd;
	smatrix_dof * K;

	// Size of the triangular matrix storage
	inline unsigned size() const {
		return nnd*(nnd+1)/2;
	}
	// Index into the triangular matrix storage
	inline unsigned index(const unsigned i, const unsigned j) const {
		return j*(j+1)/2 + i;
	}
};

/*
 * This is a specialised version of the general inv_mul_gauss function,
 * which does not pivot, and has fixed size, since we know that we never
 * pivot for the hessian, and it is always NDIMxNDIM.  It also takes the
 * transpose of the B matrix, since the column dimension is fixed.
 */
double
inv_mul_gauss(const unsigned nnd, double (& J)[NDIM][NDIM],
		double (* B)[NDIM])
{
	double det = 1.0;
	unsigned i, j, k;

	for (i = 0; i < NDIM; i++) {
		double pvt = J[i][i];
		if (fabs(pvt) < DBL_EPSILON) {
			event_msg(EVENT_ERROR,"Singular matrix in inv_mul_gauss()!");
			return 0.0;
		}
		det *= pvt;
		for (k = i+1; i < NDIM && k < NDIM; k++) {
			double tmp = J[k][i]/pvt;
			for (j = i+1; i < NDIM && j < NDIM; j++)
				J[k][j] -= tmp*J[i][j];
			for (j = 0; j < nnd; j++)
				B[j][k] -= tmp*B[j][i];
		}
	}
	for (i = NDIM; i > 0; i--) {
		for (j = i; j < NDIM; j++)
			for (k = 0; k < nnd; k++)
				B[k][i-1] -= J[i-1][j]*B[k][j];
		for (k = 0; k < nnd; k++)
			B[k][i-1] /= J[i-1][i-1];
	}
	return det;
}

/*
 * class element
 *
 * This abstract interface defines all of the things needed to access
 * a finite element.
 *
 */
class element : public listelement_o<mesh,element> {
public:
	// The types of elements.
	enum element_t {
		block8,
		block12,
		block16,
		block18,
		block27,
		block36,
		variable18,
		variable26,
		variable34,
		adaptor18,
		adaptor26,
		adaptor34,
		infinite8,
		infinite12,
		infinite16
	};
	// The types of shape functions.  The special numbering is so
	// we can use the enum to guess the node requirements.
	enum shape_t {
		linear     = 2,
		quadratic  = 3,
		cubic      = 4,
		absolute,
		inf_pos,
		inf_neg
	};

	element(mesh * o, const material & m)
	  : listelement_o<mesh,element>(o), mat(m), nnd(0) {
	}
	virtual ~element() {
	}
	// These are declared later, once we can access our mesh's node list,
	// without the problem of having the classes depend on one another.
	inline unsigned addnode(const coord3d & c) const;
	inline void updatenode(const node3d & n) const;
	inline const node3d & getnode(const unsigned i) const;
	
	virtual unsigned l2g(const unsigned i) const = 0;
	virtual smatrix_elem * stiffness() const = 0;
	virtual void results(const fset<point3d> &, fset<pavedata> &) const = 0;
	
//protected:
	friend class mesh;

	const material & mat;
	unsigned nnd;
};

/*
 * class element_base
 *
 * A finite element.  This base class does all the work for the
 * special element classes, which basically just exist to set things
 * up right...
 *
 * This class assumes a 3D rectangular geometry with shape functions in the
 * same shape functions in x and y, and a different function in the z
 * direction.
 */
template <element::shape_t SZ, unsigned NND>
class element_base : public element {
public:
	element_base(mesh * o, const shape_t x, const shape_t y,
			const material & m)
	  : element(o,m), sx(x), sy(y) {
	}
	virtual unsigned l2g(const unsigned i) const {
		return inel[i];
	}

protected:
	shape_t sx, sy;            // Shape functions in x and y.
	unsigned inel[NND];        // Local to Global node mapping.

	// This is the core routine for creating elements.  It is passed a set
	// of eight corners, and adds the required nodes.
	inline void setup_block(const fset<coord3d> & c) {
		assert(8 == c.length());
		assert(nnd == 0);
		coord3d cc[8];
		
		memcpy(cc,&c[0],8*sizeof(coord3d));
		coord_sort(cc);
		coord_add(sx,sy,cc);
		coord_neighbours();
	}
	// This is the core routine for creating variable elements.  It sets up
	// the mask of optional nodes when we are dropping nodes.
	inline void setup_mask(const fset<coord3d> & c, unsigned * mask) {
		assert(8 == c.length());
		assert(nnd == 0);
		assert(mask != 0);
		const unsigned nz = unsigned(SZ);
		coord3d cc[8];
		
		memcpy(cc,&c[0],8*sizeof(coord3d));
		coord_sort(cc);
		coord_add(linear,linear,cc);
		
		// If we have a mask, initialise it.
		for (unsigned i = 0; i < 5*nz; i++)
			mask[i] = UINT_MAX;
		// Loop over the nodes, setting the neighbours.
		for (unsigned i = 0; i < 4*nz; i++) {
			node3d n = getnode(inel[i]);
			unsigned x_m = UINT_MAX, x_p = UINT_MAX;
			unsigned y_m = UINT_MAX, y_p = UINT_MAX;
			unsigned z_m = UINT_MAX, z_p = UINT_MAX;
			(i%4 < 2 ? x_p : x_m) = inel[4*(i/4)+(i%4+2)%4];
			(i%2 < 1 ? y_p : y_m) = inel[2*(i/2)+(i%2+1)%2];
			if (i < 4*(nz-1))
				z_p = inel[i+4];
			if (i >= 4)
				z_m = inel[i-4];
			// Check for nodes in the middle, and set the mask.
			switch (i%4) {
			case 0:
				if (n.xp != UINT_MAX && n.xp != x_p) {
					mask[i+2] = x_p = n.xp;
				}
				if (n.yp != UINT_MAX && n.yp != y_p) {
					mask[i] = y_p = n.yp;
					// Check the center for the top and bottom.
					if ((i == 0 || i == 4*(nz-1))) {
						const node3d & mid = getnode(n.yp);
						if (mid.xp != UINT_MAX)
							mask[4*nz+(i/4)] = mid.xp;
					}
				}
				break;
			case 1:
				if (n.xp != UINT_MAX && n.xp != x_p)
					x_p = n.xp;
				if (n.ym != UINT_MAX && n.ym != y_m)
					y_m = n.ym;
				break;
			case 2:
				if (n.xm != UINT_MAX && n.xm != x_m)
					x_m = n.xm;
				if (n.yp != UINT_MAX && n.yp != y_p)
					y_p = n.yp;
				break;
			case 3:
				if (n.xm != UINT_MAX && n.xm != x_m) {
					mask[i-2] = x_m = n.xm;
				}
				if (n.ym != UINT_MAX && n.ym != y_m) {
					mask[i] = y_m = n.ym;
				}
				break;
			}
			n.setneighbours(x_m,x_p,y_m,y_p,z_m,z_p);
			updatenode(n);
		}
		// Add the extra nodes to our list and set the mask from the global
		// node number to the item number in the list.
		for (unsigned i = 0; i < 5*nz; i++) {
			if (mask[i] != UINT_MAX) {
				inel[nnd] = mask[i];
				mask[i] = nnd++;
			} else
				mask[i] = 0;
		}
	}
	// Build a matrix of the nodal coordinates.
	void buildxe(double (* xe)[NDIM]) const {
		assert(NDIM == 3);
		for (unsigned i = 0; i < nnd; i++) {
			const coord3d & c = getnode(inel[i]);
			xe[i][0] = c.x; xe[i][1] = c.y;
			xe[i][2] = node_depth_callback(c,owner,&mat);
		}
		if (sx == inf_pos) {
			xe[2][2] = xe[0][2]; xe[3][2] = xe[1][2];
			xe[6][2] = xe[4][2]; xe[7][2] = xe[5][2];
		}
		if (sx == inf_neg) {
			xe[0][2] = xe[2][2]; xe[1][2] = xe[3][2];
			xe[4][2] = xe[6][2]; xe[5][2] = xe[7][2];
		}
		if (sy == inf_pos) {
			xe[1][2] = xe[0][2]; xe[3][2] = xe[2][2];
			xe[5][2] = xe[4][2]; xe[7][2] = xe[6][2];
		}
		if (sy == inf_neg) {
			xe[0][2] = xe[1][2]; xe[2][2] = xe[3][2];
			xe[4][2] = xe[5][2]; xe[6][2] = xe[7][2];
		}
	}
	// Build a matrix of the emod delta values.
	void buildEe(double * Ee) const {
		for (unsigned i = 0; i < nnd; i++) {
			const coord3d & c = getnode(inel[i]);
			Ee[i] = node_emod_callback(c,owner,&mat) - mat.emod;
		}
		if (sx == inf_pos) {
			Ee[2] = Ee[0]; Ee[3] = Ee[1];
			Ee[6] = Ee[4]; Ee[7] = Ee[5];
		}
		if (sx == inf_neg) {
			Ee[0] = Ee[2]; Ee[1] = Ee[3];
			Ee[4] = Ee[6]; Ee[5] = Ee[7];
		}
		if (sy == inf_pos) {
			Ee[1] = Ee[0]; Ee[3] = Ee[2];
			Ee[5] = Ee[4]; Ee[7] = Ee[6];
		}
		if (sy == inf_neg) {
			Ee[0] = Ee[1]; Ee[2] = Ee[3];
			Ee[4] = Ee[5]; Ee[6] = Ee[7];
		}
	}
	// Build a matrix of the nodal deflections.
	void buildue(double (* ue)[NDOF]) const {
		assert(NDOF == 3);
		for (unsigned i = 0; i < nnd; i++) {
			const node3d & n = getnode(inel[i]);
			ue[i][0] = n.ux; ue[i][1] = n.uy; ue[i][2] = n.uz;
		}
	}
	// build the shape function vector and derivative matrix.
	virtual void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const = 0;
	// build the shape functions for linear and infinite elements
	// in plan view.
	void buildSF_block(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const {
		const unsigned ny = getn(sy),     my = 1;
		const unsigned nx = getn(sx),     mx = my*ny;
		const unsigned nz = unsigned(SZ), mz = mx*nx;
		assert(nnd == nz*mz);
		double Nx, Ny, Nz, dNdx, dNdy, dNdz;

		assert(NDIM == 3);
		// Loop over each level in the element.
		for (unsigned i = 0; i < nz; i++) {
			Nz = getN(SZ,i,gz,ismap); dNdz = getdN(SZ,i,gz,ismap);
			for (unsigned j = 0; j < nx; j++) {
				Nx = getN(sx,j,gx,ismap); dNdx = getdN(sx,j,gx,ismap);
				for (unsigned k = 0; k < ny; k++) {
					Ny = getN(sy,k,gy,ismap); dNdy = getdN(sy,k,gy,ismap);
					unsigned l = mz*i+mx*j+my*k;
					N[l] = Nx*Ny*Nz;
					dNdr[l][0] = dNdx*Ny*Nz;
					dNdr[l][1] = Nx*dNdy*Nz;
					dNdr[l][2] = Nx*Ny*dNdz;
				}
			}
		}
	}
	// Special shape function builder for variable node elements.
	void buildSF_variable(const double & gx, const double & gy,
			const double & gz, double * N, double (* dNdr)[NDIM],
			const unsigned * mask) const {
		assert(sx == absolute);
		assert(sy == absolute);
		const unsigned nz = unsigned(SZ);
		unsigned j;
		double Nx, Ny, Nz, dNdx, dNdy, dNdz;

		assert(NDIM == 3);
		// Loop over each level in the element.
		for (unsigned i = 0; i < nz; i++) {
			Nz = getN(SZ,i,gz); dNdz = getdN(SZ,i,gz);
			// Start with the corners.
			for (unsigned l = 0; l < 4; l++) {
				Nx = N_2(l/2,gx); dNdx = dN_2(l/2,gx);
				Ny = N_2(l%2,gy); dNdy = dN_2(l%2,gy);
				N[i*4+l] = Nx*Ny*Nz;
				dNdr[i*4+l][0] = dNdx*Ny*Nz;
				dNdr[i*4+l][1] = Nx*dNdy*Nz;
				dNdr[i*4+l][2] = Nx*Ny*dNdz;
				// Build the extra nodes at the half points.
				if ((j = mask[i*4+l]) == 0)
					continue;
				if (l == 0 || l == 3) {
					Ny = N_A(1,gy); dNdy = dN_A(1,gy);
				} else {
					Nx = N_A(1,gx); dNdx = dN_A(1,gx);
				}
				N[j] = Nx*Ny*Nz;
				dNdr[j][0] = dNdx*Ny*Nz;
				dNdr[j][1] = Nx*dNdy*Nz;
				dNdr[j][2] = Nx*Ny*dNdz;
			}
			// Build the top and bottom.
			if ((j = mask[4*nz+i]) == 0)
				continue;
			Nx = N_A(1,gx); dNdx = dN_A(1,gx);
			Ny = N_A(1,gy); dNdy = dN_A(1,gy);
			N[j] = Nx*Ny*Nz;
			dNdr[j][0] = dNdx*Ny*Nz;
			dNdr[j][1] = Nx*dNdy*Nz;
			dNdr[j][2] = Nx*Ny*dNdz;
		}
		// Correct the shape functions for the corners and sides,
		// depending on the existance of the extra nodes.
		for (unsigned l = 0; l < 4*nz; l++) {
			if ((j = mask[l]) != 0) {
				N[l] -= N[j]/2;
				dNdr[l][0] -= dNdr[j][0]/2;
				dNdr[l][1] -= dNdr[j][1]/2;
				dNdr[l][2] -= dNdr[j][2]/2;
			}
			const int Jxy[4] = {+2,-1,+1,-2};
			if ((j = mask[l+Jxy[l%4]]) != 0) {
				N[l] -= N[j]/2;
				dNdr[l][0] -= dNdr[j][0]/2;
				dNdr[l][1] -= dNdr[j][1]/2;
				dNdr[l][2] -= dNdr[j][2]/2;
			}
			if ((j = mask[4*nz+(l/4)]) != 0) {
				N[l] += N[j]/4;
				dNdr[l][0] += dNdr[j][0]/4;
				dNdr[l][1] += dNdr[j][1]/4;
				dNdr[l][2] += dNdr[j][2]/4;
			}
		}
		for (unsigned l = 0; l < 4*nz; l++) {
			if ((j = mask[4*nz+(l/4)]) != 0) {
				N[mask[l]] -= N[j]/2;
				dNdr[mask[l]][0] -= dNdr[j][0]/2;
				dNdr[mask[l]][1] -= dNdr[j][1]/2;
				dNdr[mask[l]][2] -= dNdr[j][2]/2;
			}
		}
	}
	// build the shape function vector and derivative matrix.
	void buildSF_adaptor(const double & gx, const double & gy,
			const double & gz, double * N, double (* dNdr)[NDIM],
			const unsigned * mask) const {
		assert(sx == quadratic);
		assert(sy == quadratic);
		const unsigned nz = unsigned(SZ);
		double Nx, Ny, Nz, dNdx, dNdy, dNdz;

		assert(NDIM == 3);
		// Loop over each level in the element.
		for (unsigned i = 0, j; i < nz; i++) {
			Nz = getN(SZ,i,gz); dNdz = getdN(SZ,i,gz);
			const int Jx[4] = {+2, 0, 0,-2};
			const int Jy[4] = { 0,-1,+1, 0};
			// Start with the corners.
			for (unsigned l = 0; l < 4; l++) {
				if (mask[i*4+l+Jx[l]] != 0) {
					Nx = N_3(2*(l/2),gx); dNdx = dN_3(2*(l/2),gx);
				} else {
					Nx = N_2(l/2,gx);     dNdx = dN_2(l/2,gx);
				}
				if (mask[i*4+l+Jy[l]] != 0) {
					Ny = N_3(2*(l%2),gy); dNdy = dN_3(2*(l%2),gy);
				} else {
					Ny = N_2(l%2,gy);     dNdy = dN_2(l%2,gy);
				}
				N[i*4+l] = Nx*Ny*Nz;
				dNdr[i*4+l][0] = dNdx*Ny*Nz;
				dNdr[i*4+l][1] = Nx*dNdy*Nz;
				dNdr[i*4+l][2] = Nx*Ny*dNdz;
				// Build the extra nodes at the half points.
				if ((j = mask[i*4+l]) == 0)
					continue;
				if (l == 0 || l == 3) {
					Ny = N_3(1,gy); dNdy = dN_3(1,gy);
				} else {
					Nx = N_3(1,gx); dNdx = dN_3(1,gx);
				}
				N[j] = Nx*Ny*Nz;
				dNdr[j][0] = dNdx*Ny*Nz;
				dNdr[j][1] = Nx*dNdy*Nz;
				dNdr[j][2] = Nx*Ny*dNdz;
			}
			if ((j = mask[4*nz+i]) == 0)
				continue;
			Nx = N_3(1,gx); dNdx = dN_3(1,gx);
			Ny = N_3(1,gy); dNdy = dN_3(1,gy);
			N[j] = Nx*Ny*Nz;
			dNdr[j][0] = dNdx*Ny*Nz;
			dNdr[j][1] = Nx*dNdy*Nz;
			dNdr[j][2] = Nx*Ny*dNdz;
		}
	}
	// Build the element stiffness matrix.
	smatrix_elem * buildKe() const {
		unsigned i, j, k, l;

		// Build the point stiffness tensor E_abcd
		smatrix_dof E[NDIM][NDIM];
		mat.pointstiffness(E);

		// Element nodal coords
		double (* xe)[NDIM] = static_cast<double (*)[NDIM]>
				(alloca(nnd*NDIM*sizeof(double)));
		buildxe(xe);
		double * Ee = static_cast<double *>(alloca(nnd*sizeof(double)));
		buildEe(Ee);

		// Build a list of gauss points.  This will need to be stored
		// in the element when we do plastic models with history.
		unsigned nx = 0, ny = 0, nz = 0;
		const double (* gx)[2] = 0, (* gy)[2] = 0, (* gz)[2] = 0;
		getgauss(sx,nx,gx);
		getgauss(sy,ny,gy);
		getgauss(SZ,nz,gz);

		smatrix_elem * K = new smatrix_elem(nnd);
		if (K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}

		double * N = static_cast<double *>(alloca(nnd*sizeof(double)));
		// This matrix is transposed.
		double (* dNdr)[NDIM] = static_cast<double (*)[NDIM]>
				(alloca(nnd*NDIM*sizeof(double)));
		double J[NDIM][NDIM];
		// Cheat on indenting a little to save screen width...
		for (unsigned gi = 0; gi < nx; gi++) {
		for (unsigned gj = 0; gj < ny; gj++) {
		for (unsigned gk = 0; gk < nz; gk++) {
			double gw = gx[gi][1]*gy[gj][1]*gz[gk][1];
			buildSF(true,gx[gi][0],gy[gj][0],gz[gk][0],N,dNdr);
			for (i = 0; i < NDIM; i++) {
				for (j = 0; j < NDIM; j++) {
					J[i][j] = 0.0;
					for (k = 0; k < nnd; k++)
						J[i][j] += dNdr[k][i]*xe[k][j];
				}
			}
			if (sx == inf_pos || sx == inf_neg
			 || sy == inf_pos || sy == inf_neg)
				buildSF(false,gx[gi][0],gy[gj][0],gz[gk][0],N,dNdr);
			double emod = 0.0;
			for (k = 0; k < nnd; k++)
				emod += N[k]*Ee[k];
			// This returns det(J);
			gw *= inv_mul_gauss(nnd,J,dNdr)*(1+emod/mat.emod);
			assert(gw > 0.0);
			for (i = 0; i < nnd; i++) {
				for (j = i; j < nnd; j++) {
					for (k = 0; k < NDIM; k++)
						for (l = 0; l < NDIM; l++)
							(*K)(i,j) += E[k][l]*(dNdr[i][k]*dNdr[j][l]*gw);
				}
			}
		}}}
		return K;
	}
	// Output results to a pavedata structure
	void builddata(const fset<point3d> & c, fset<pavedata> & data) const {
		unsigned n, i, j, k;
		double l[NDIM], g[NDIM];

		assert(NDIM == 3 && NDOF == 3);
		assert(c.length() == data.length());
		
		// Element nodal coords
		double (* xe)[NDIM] = static_cast<double (*)[NDIM]>
				(alloca(nnd*NDIM*sizeof(double)));
		buildxe(xe);
		// Element nodal deflections
		double (* ue)[NDOF] = static_cast<double (*)[NDOF]>
				(alloca(nnd*NDOF*sizeof(double)));
		buildue(ue);
		// Element nodal elastic modulus deltas
		double * Ee = static_cast<double *>(alloca(nnd*sizeof(double)));
		buildEe(Ee);

		double * N = static_cast<double *>(alloca(nnd*sizeof(double)));
		// This matrix is transposed.
		double (* dNdr)[NDIM] = static_cast<double (*)[NDIM]>
				(alloca(nnd*NDIM*sizeof(double)));
		double J[NDIM][NDIM];
		for (n = 0; n < c.length(); n++) {
			pavedata & d = data[n];
			// Stress and strain tensors.
			l[0] = c[n].x, l[1] = c[n].y, l[2] = c[n].z;
			buildSF(true,l[0],l[1],l[2],N,dNdr);
			for (j = 0; j < NDIM; j++) {
				g[j] = 0.0;
				for (k = 0; k < nnd; k++)
					g[j] += N[k]*xe[k][j];
			}
			d.x = g[0], d.y = g[1], d.z = g[2];
			for (i = 0; i < NDIM; i++) {
				for (j = 0; j < NDIM; j++) {
					J[i][j] = 0.0;
					for (k = 0; k < nnd; k++)
						J[i][j] += dNdr[k][i]*xe[k][j];
				}
			}
			if (sx == inf_pos || sx == inf_neg
			 || sy == inf_pos || sy == inf_neg)
				buildSF(false,l[0],l[1],l[2],N,dNdr);
			inv_mul_gauss(nnd,J,dNdr);
			for (j = 0; j < NDOF; j++) {
				g[j] = 0.0;
				for (k = 0; k < nnd; k++)
					g[j] += N[k]*ue[k][j];
			}
			double emod = 0.0;
			for (k = 0; k < nnd; k++)
				emod += N[k]*Ee[k];
			ematrix e(0.0);
			for (i = 0; i < NDIM; i++) {
				for (j = 0; j < NDOF; j++) {
					for (k = 0; k < nnd; k++)
						e(i,j) += dNdr[k][i]*ue[k][j];
				}
			}
			e = (e + ~e)*0.5;
			ematrix s(mat.pointstress(e)*(1+emod/mat.emod));
			ematrix pe(principle(e));
			ematrix ps(principle(s));
			d.data[4][0] = g[0];
			d.data[4][1] = g[1];
			d.data[4][2] = g[2];
			d.data[0][0] = s(0,0);
			d.data[0][1] = s(1,1);
			d.data[0][2] = s(2,2);
			d.data[1][0] = s(0,1);
			d.data[1][1] = s(0,2);
			d.data[1][2] = s(1,2);
			d.data[2][0] = ps(0,0);
			d.data[2][1] = ps(1,1);
			d.data[2][2] = ps(2,2);
			d.data[3][0] = fabs(ps(0,0) - ps(2,2))*0.5;
			d.data[3][1] = fabs(ps(0,0) - ps(1,1))*0.5;
			d.data[3][2] = fabs(ps(1,1) - ps(2,2))*0.5;
			if (d.data[3][1] < d.data[3][2])
				swap(d.data[3][1],d.data[3][2]);
			d.data[5][0] = e(0,0);
			d.data[5][1] = e(1,1);
			d.data[5][2] = e(2,2);
			d.data[6][0] = 2*e(0,1);
			d.data[6][1] = 2*e(0,2);
			d.data[6][2] = 2*e(1,2);
			d.data[7][0] = pe(0,0);
			d.data[7][1] = pe(1,1);
			d.data[7][2] = pe(2,2);
			d.data[8][0] = fabs(pe(0,0) - pe(2,2));
			d.data[8][1] = fabs(pe(0,0) - pe(1,1));
			d.data[8][2] = fabs(pe(1,1) - pe(2,2));
			if (d.data[8][1] < d.data[8][2])
				swap(d.data[8][1],d.data[8][2]);
		}
	}

private:
	// This sorts the corners into the order we expect.
	inline void coord_sort(coord3d (& cc)[8]) {
		coord3d p, m;

		do {
			m.x = cc[0].x + cc[1].x + cc[4].x + cc[5].x;
			p.x = cc[2].x + cc[3].x + cc[6].x + cc[7].x;
			m.y = cc[0].y + cc[2].y + cc[4].y + cc[6].y;
			p.y = cc[1].y + cc[3].y + cc[5].y + cc[7].y;
			m.z = cc[0].z + cc[1].z + cc[2].z + cc[3].z;
			p.z = cc[4].z + cc[5].z + cc[6].z + cc[7].z;
			if (p.x > m.x && p.y > m.y && p.z > m.z)
				break;
			if (p.x < m.x) {
				swap(cc[0],cc[2]); swap(cc[4],cc[6]);
				swap(cc[1],cc[3]); swap(cc[5],cc[7]);
				continue;
			}
			if (p.y < m.y) {
				swap(cc[0],cc[1]); swap(cc[4],cc[5]);
				swap(cc[2],cc[3]); swap(cc[6],cc[7]);
				continue;
			}
			if (p.z < m.z) {
				swap(cc[0],cc[4]); swap(cc[2],cc[6]);
				swap(cc[1],cc[5]); swap(cc[3],cc[7]);
				continue;
			}
			event_msg(EVENT_ERROR,"Bad or degenerate element shape!");
			return;
		} while(true);
	}
	// This adds nodes according to linear shape functions...
	inline void coord_add(shape_t ssx, shape_t ssy,
			const coord3d (& cc)[8]) {
		assert(nnd == 0);
		const unsigned nx = getn(ssx),    mx = nx-1;
		const unsigned ny = getn(ssy),    my = ny-1;
		const unsigned nz = unsigned(SZ), mz = nz-1;

		// Figure out our nodes using linear shape functions...
		for (unsigned i = 0; i < nz; i++) {
			for (unsigned j = 0; j < nx; j++) {
				for (unsigned k = 0; k < ny; k++) {
					coord3d p;
					p.x =  (double(cc[0].x)*(mz-i)*(mx-j)*(my-k)
					      + double(cc[1].x)*(mz-i)*(mx-j)*(   k)
					      + double(cc[2].x)*(mz-i)*(   j)*(my-k)
					      + double(cc[3].x)*(mz-i)*(   j)*(   k)
					      + double(cc[4].x)*(   i)*(mx-j)*(my-k)
					      + double(cc[5].x)*(   i)*(mx-j)*(   k)
					      + double(cc[6].x)*(   i)*(   j)*(my-k)
					      + double(cc[7].x)*(   i)*(   j)*(   k))/(mz*mx*my);
					p.y =  (double(cc[0].y)*(mz-i)*(mx-j)*(my-k)
					      + double(cc[1].y)*(mz-i)*(mx-j)*(   k)
					      + double(cc[2].y)*(mz-i)*(   j)*(my-k)
					      + double(cc[3].y)*(mz-i)*(   j)*(   k)
					      + double(cc[4].y)*(   i)*(mx-j)*(my-k)
					      + double(cc[5].y)*(   i)*(mx-j)*(   k)
					      + double(cc[6].y)*(   i)*(   j)*(my-k)
					      + double(cc[7].y)*(   i)*(   j)*(   k))/(mz*mx*my);
					p.z =  (double(cc[0].z)*(mz-i)*(mx-j)*(my-k)
					      + double(cc[1].z)*(mz-i)*(mx-j)*(   k)
					      + double(cc[2].z)*(mz-i)*(   j)*(my-k)
					      + double(cc[3].z)*(mz-i)*(   j)*(   k)
					      + double(cc[4].z)*(   i)*(mx-j)*(my-k)
					      + double(cc[5].z)*(   i)*(mx-j)*(   k)
					      + double(cc[6].z)*(   i)*(   j)*(my-k)
					      + double(cc[7].z)*(   i)*(   j)*(   k))/(mz*mx*my);
					inel[nnd++] = addnode(p);
				}
			}
		}
	}
	// Set the neighbours according to a simple pattern (no mask).
	inline void coord_neighbours() {
		const unsigned ny = getn(sy),     my = 1;
		const unsigned nx = getn(sx),     mx = my*ny;
		const unsigned nz = unsigned(SZ), mz = mx*nx;
		assert(nnd == nz*mz);

		// Loop over the nodes, setting the neighbours.
		for (unsigned i = 0; i < nz; i++) {
			for (unsigned j = 0; j < nx; j++) {
				for (unsigned k = 0; k < ny; k++) {
					node3d n = getnode(inel[mz*i+mx*j+my*k]);
					unsigned x_m, x_p, y_m, y_p, z_m, z_p;
					x_p = (j < nx-1 ? inel[mz*i+mx*j+my*k+mx] : UINT_MAX);
					x_m = (j > 0    ? inel[mz*i+mx*j+my*k-mx] : UINT_MAX);
					y_p = (k < ny-1 ? inel[mz*i+mx*j+my*k+my] : UINT_MAX);
					y_m = (k > 0    ? inel[mz*i+mx*j+my*k-my] : UINT_MAX);
					z_p = (i < nz-1 ? inel[mz*i+mx*j+my*k+mz] : UINT_MAX);
					z_m = (i > 0    ? inel[mz*i+mx*j+my*k-mz] : UINT_MAX);
					n.setneighbours(x_m,x_p,y_m,y_p,z_m,z_p);
					updatenode(n);
				}
			}
		}
	}
	// This gets the right gauss function based on the shape function.
	static inline void getgauss(const shape_t t, unsigned & n,
			const double (*& g)[2]) {
		switch (t) {
			case linear:    n = 2; g = gp_2; return;
			case quadratic: n = 3; g = gp_3; return;
			case cubic:     n = 4; g = gp_4; return;
			case absolute:  n = 4; g = gp_A; return;
			case inf_pos:
			case inf_neg:   n = 4; g = gp_4; return;
		}
	}
	// Linear shape function.
	static inline double N_2(const unsigned n, const double r) {
		switch (n) {
		case 0: return 0.5*(1.0-r);
		case 1: return 0.5*(1.0+r);
		} assert(false); return 0.0;
	}
	static inline double dN_2(const unsigned n, const double r) {
		switch (n) {
		case 0: return -0.5;
		case 1: return  0.5;
		} assert(false); return 0.0;
	}
	// Absolute value shape function.
	static inline double N_A(const unsigned n, const double r) {
		switch (n) {
		case 0: return (r < 0.0 ? -r : 0.0);
		case 1: return 1-fabs(r);
		case 2: return (r > 0.0 ?  r : 0.0);
		} assert(false); return 0.0;
	}
	static inline double dN_A(const unsigned n, const double r) {
		switch (n) {
		case 0: return (r < 0.0 ? -1.0 : r > 0.0 ?  0.0 : -0.5);
		case 1: return (r < 0.0 ?  1.0 : r > 0.0 ? -1.0 :  0.0);
		case 2: return (r < 0.0 ?  0.0 : r > 0.0 ?  1.0 :  0.5);
		} assert(false); return 0.0;
	}
	// Quadratic shape function.
	static inline double N_3(const unsigned n, const double r) {
		switch (n) {
		case 0: return 0.5*r*(r-1);
		case 1: return 1-r*r;
		case 2: return 0.5*r*(r+1);
		} assert(false); return 0.0;
	}
	static inline double dN_3(const unsigned n, const double r) {
		switch (n) {
		case 0: return r-0.5;
		case 1: return -2*r;
		case 2: return r+0.5;
		} assert(false); return 0.0;
	}
	// Cubic shape function.
	static inline double N_4(const unsigned n, const double r) {
		switch (n) {
		case 0: return 0.0625*(-1+9*r*r)*(1-1*r);
		case 1: return 0.0625*(+9-9*r*r)*(1-3*r);
		case 2: return 0.0625*(+9-9*r*r)*(1+3*r);
		case 3: return 0.0625*(-1+9*r*r)*(1+1*r);
		} assert(false); return 0.0;
	}
	static inline double dN_4(const unsigned n, const double r) {
		switch (n) {
		case 0: return 0.0625*((+18*r)*(1-1*r) - 1*(-1+9*r*r));
		case 1: return 0.0625*((-18*r)*(1-3*r) - 3*(+9-9*r*r));
		case 2: return 0.0625*((-18*r)*(1+3*r) + 3*(+9-9*r*r));
		case 3: return 0.0625*((+18*r)*(1+1*r) + 1*(-1+9*r*r));
		} assert(false); return 0.0;
	}
	// mapping functions for infinite elements.
	static inline double Mpi(const unsigned n, const double r) {
		switch (n) {
		case 0: return  -2*r/(1-r);
		case 1: return (1+r)/(1-r);
		} assert(false); return 0.0;
	}
	static inline double dMpi(const unsigned n, const double r) {
		switch (n) {
		case 0: return -2/((1-r)*(1-r));
		case 1: return  2/((1-r)*(1-r));
		} assert(false); return 0.0;
	}
	static inline double Mni(const unsigned n, const double r) {
		switch (n) {
		case 0: return (1-r)/(1+r);
		case 1: return   2*r/(1+r);
		} assert(false); return 0.0;
	}
	static inline double dMni(const unsigned n, const double r) {
		switch (n) {
		case 0: return -2/((1+r)*(1+r));
		case 1: return  2/((1+r)*(1+r));
		} assert(false); return 0.0;
	}
	// Shape functions for infinite elements are quadratic.
	static inline double Npi(const unsigned n, const double r) {
		return N_3(n,r);
	}
	static inline double dNpi(const unsigned n, const double r) {
		return dN_3(n,r);
	}
	// Shift the negative nodes to (0,+1).
	static inline double Nni(const unsigned n, const double r) {
		return N_3(n+1,r);
	}
	static inline double dNni(const unsigned n, const double r) {
		return dN_3(n+1,r);
	}
	// Return the number of nodes in each shape function.
	static inline unsigned getn(const shape_t s, const bool ismap = true) {
		switch (s) {
		case linear:    return 2;
		case absolute:  return 3;
		case quadratic: return 3;
		case cubic:     return 4;
		case inf_pos:   return (ismap ? 2 : 2);
		case inf_neg:   return (ismap ? 2 : 2);
		} assert(false); return 0;
	}
	// Build the shape/mapping function based on the type of function.
	static inline double getN(const shape_t s, const unsigned n,
			const double r, const bool ismap = true) {
		switch (s) {
		case linear:    return N_2(n,r);
		case absolute:  return N_A(n,r);
		case quadratic: return N_3(n,r);
		case cubic:     return N_4(n,r);
		case inf_pos:   return (ismap ? Mpi(n,r) : Npi(n,r));
		case inf_neg:   return (ismap ? Mni(n,r) : Nni(n,r));
		} assert(false); return 0.0;
	}
	static inline double getdN(const shape_t s, const unsigned n,
			const double r, const bool ismap = true) {
		switch (s) {
		case linear:    return dN_2(n,r);
		case absolute:  return dN_A(n,r);
		case quadratic: return dN_3(n,r);
		case cubic:     return dN_4(n,r);
		case inf_pos:   return (ismap ? dMpi(n,r) : dNpi(n,r));
		case inf_neg:   return (ismap ? dMni(n,r) : dNni(n,r));
		} assert(false); return 0.0;
	}
};

/*
 * class element_block - a finite element
 */
template<element::shape_t SZ>
class element_block : public element_base<SZ,4*int(SZ)> {
	typedef element_base<SZ,4*int(SZ)> base_t;
public:
	element_block(mesh * o, const material & m, const fset<coord3d> & c)
	  : base_t(o,element::linear,element::linear,m) {
		base_t::setup_block(c);
	}
	virtual void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const {
		assert(4*int(SZ) == this->nnd);
		base_t::buildSF_block(ismap,gx,gy,gz,N,dNdr);
	}
	virtual smatrix_elem * stiffness() const {
		assert(4*int(SZ) == this->nnd);
		return base_t::buildKe();
	}
	virtual void results(const fset<point3d> & c, fset<pavedata> & d) const {
		assert(4*int(SZ) == this->nnd);
		base_t::builddata(c,d);
	}
};

/*
 * class element_quad - a higher order finite element
 */
template<element::shape_t SZ>
class element_quad : public element_base<SZ,9*int(SZ)> {
	typedef element_base<SZ,9*int(SZ)> base_t;
public:
	element_quad(mesh * o, const material & m, const fset<coord3d> & c)
	  : base_t(o,element::quadratic,element::quadratic,m) {
		base_t::setup_block(c);
	}
	virtual void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const {
		assert(9*int(SZ) == this->nnd);
		base_t::buildSF_block(ismap,gx,gy,gz,N,dNdr);
	}
	virtual smatrix_elem * stiffness() const {
		assert(9*int(SZ) == this->nnd);
		return base_t::buildKe();
	}
	virtual void results(const fset<point3d> & c, fset<pavedata> & d) const {
		assert(9*int(SZ) == this->nnd);
		base_t::builddata(c,d);
	}
};

/*
 * class element_variable - a variable node finite element
 */
template<element::shape_t SZ>
class element_variable : public element_base<SZ,9*int(SZ)> {
	typedef element_base<SZ,9*int(SZ)> base_t;
public:
	element_variable(mesh * o, const material & m, const fset<coord3d> & c)
	  : base_t(o,element::absolute,element::absolute,m) {
		base_t::setup_mask(c,mask);
	}
	virtual void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const {
		assert(ismap == true);
		base_t::buildSF_variable(gx,gy,gz,N,dNdr,mask);
	}
	virtual smatrix_elem * stiffness() const {
		return base_t::buildKe();
	}
	virtual void results(const fset<point3d> & c, fset<pavedata> & d) const {
		base_t::builddata(c,d);
	}

protected:
	unsigned mask[5*int(SZ)];
};

/*
 * class element_adaptor - an element to adapt from quadratic to linear
 * elements.
 */
template<element::shape_t SZ>
class element_adaptor : public element_base<SZ,9*int(SZ)> {
	typedef element_base<SZ,9*int(SZ)> base_t;
public:
	element_adaptor(mesh * o, const material & m, const fset<coord3d> & c)
	  : base_t(o,element::quadratic,element::quadratic,m) {
		base_t::setup_mask(c,mask);
	}
	virtual void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const {
		assert(ismap == true);
		base_t::buildSF_adaptor(gx,gy,gz,N,dNdr,mask);
	}
	virtual smatrix_elem * stiffness() const {
		return base_t::buildKe();
	}
	virtual void results(const fset<point3d> & c, fset<pavedata> & d) const {
		base_t::builddata(c,d);
	}

protected:
	unsigned mask[5*int(SZ)];
};

/*
 * class element_infinite - a infinite element
 */
template<element::shape_t SZ>
class element_infinite : public element_base<SZ,4*int(SZ)> {
	typedef element_base<SZ,4*int(SZ)> base_t;
public:
	// For infinite elements we first have to find the eight corners, based
	// on the line or plane given.  Then we can set up the rest.
	element_infinite(mesh * o, const material & m, const fset<coord3d> & c)
	  : base_t(o,element::linear,element::linear,m) {
		fset<coord3d> cc(8);

		if (c.length() == 2) {
			// Corner infinite elements are defined by 2 points in
			// a line, which must be vertical.
			assert(c[0].x == c[1].x);
			assert(c[0].y == c[1].y);
			for (unsigned i = 0; i < 4; i++) {
				double x = c[0].x, y = c[0].y;
				x *= ((x > 0) == (i%4 < 2) ? 1 : 2);
				y *= ((y > 0) == (i%2 < 1) ? 1 : 2);
				cc[i  ] = c[0]; cc[i  ].x = x; cc[i  ].y = y;
				cc[i+4] = c[1]; cc[i+4].x = x; cc[i+4].y = y;
			}
			base_t::sx = (double(c[0].x) > 0 ?
					base_t::inf_pos : base_t::inf_neg);
			base_t::sy = (double(c[0].y) > 0 ?
					base_t::inf_pos : base_t::inf_neg);
		} else {
			// Side infinite elements are defined by 4 points in a plane.
			assert(4 == c.length());
			if (c[0].x == c[3].x) {
				assert(c[0].x == c[1].x);
				assert(c[1].x == c[2].x);
				assert(c[2].x == c[3].x);
				for (unsigned i = 0; i < 4; i++) {
					double x = c[0].x;
					x *= ((x > 0) == (i%4 < 2) ? 1 : 2);
					cc[i  ] = c[i%2  ]; cc[i  ].x = x;
					cc[i+4] = c[i%2+2]; cc[i+4].x = x;
				}
				base_t::sx = (double(c[0].x) > 0 ?
					base_t::inf_pos : base_t::inf_neg);
			} else {
				assert(c[0].y == c[1].y);
				assert(c[1].y == c[2].y);
				assert(c[2].y == c[3].y);
				for (unsigned i = 0; i < 4; i++) {
					double y = c[0].y;
					y *= ((y > 0) == (i%2 < 1) ? 1 : 2);
					cc[i  ] = c[(i<2?0:1)]; cc[i  ].y = y;
					cc[i+4] = c[(i<2?2:3)]; cc[i+4].y = y;
				}
				base_t::sy = (double(c[0].y) > 0 ?
					base_t::inf_pos : base_t::inf_neg);
			}
		}
		base_t::setup_block(cc);
	}
	virtual void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const {
		assert(4*int(SZ) == this->nnd);
		base_t::buildSF_block(ismap,gx,gy,gz,N,dNdr);
	}
	virtual smatrix_elem * stiffness() const {
		assert(4*int(SZ) == this->nnd);
		return base_t::buildKe();
	}
	virtual void results(const fset<point3d> & c, fset<pavedata> & d) const {
		assert(4*int(SZ) == this->nnd);
		base_t::builddata(c,d);
	}
};

class smatrix_diag;            // Forward declare.
class smatrix;                 // Forward declare.

/*
 * class smatrix_node - a node in a stiffness matrix
 */
class smatrix_node {
	friend class smatrix_diag;
	friend class smatrix;
	friend class mesh;

	smatrix_dof K;
	unsigned i;
	smatrix_node * col_next;
};

/*
 * class smatrix_diag - a diagonal in a stiffness matrix
 */
class smatrix_diag {
	// Add one element to the storage, fixing up as needed.
	// The new elements will always start at nodes[nnz];
	bool expand(unsigned nn) {
		if (nnz+nn > nnd) {
			// We're out of space, so allocate more.
			nnd = nnz+nn;
			smatrix_node * temp =
					static_cast<smatrix_node *>
					(realloc(nodes,nnd*sizeof(smatrix_node)));
			if (temp == 0) {
				event_msg(EVENT_ERROR,"Out of memory in smatrix_diag::insert()!");
				return false;
			}
			// If realloc moved us we need to fix up everyone's pointers.
			// This is finds the offset into the old array, and then adds
			// that to the new array.
			if (nodes != temp) {
				if (col_head != 0)
					col_head = temp + (col_head - nodes);
				for (unsigned k = 0; k < nnz; k++) {
					smatrix_node * n = &(temp[k]);
					if (n->col_next != 0)
						n->col_next = temp + (n->col_next - nodes);
				}
				nodes = temp;
			}
		}
		return true;
	}
	
	friend class smatrix;
	friend class mesh;

	smatrix_dof K;
	smatrix_node * col_head;
	smatrix_node * nodes;
	unsigned nnz, nnd;
};

/*
 * class smatrix - a symmetric positive definite stiffness matrix
 */
class smatrix {
public:
	// Create the matrix empty.
	inline explicit smatrix(const unsigned n)
	  : nnd(n), diag(0) {
		diag = static_cast<smatrix_diag *>
				(calloc(nnd,sizeof(smatrix_diag)));
		if (diag == 0)
			event_msg(EVENT_ERROR,"Out of memory in smatrix::smatrix()!");
	}
	// Copy a matrix.
	inline explicit smatrix(const smatrix & A)
	  : nnd(A.nnd), diag(0) {
		diag = static_cast<smatrix_diag *>
				(calloc(nnd,sizeof(smatrix_diag)));
		if (diag == 0) {
			event_msg(EVENT_ERROR,"Out of memory in smatrix::smatrix()!");
			return;
		}
		for (unsigned i = 0; i < nnd; i++) {
			smatrix_diag * d = &(A.diag[i]);
			// Pre-allocate the space.
			if (!diag[i].expand(d->nnz))
				return;
			append(i,i,d->K);
			smatrix_node * p = d->col_head;
			while (p) {
				append(p->i,i,p->K);
				p = p->col_next;
			}
		}
	}
	inline ~smatrix() {
		if (diag) {
			for (unsigned i = 0; i < nnd; i++) {
				// smatrix_nodes don't have a destructor...
				if (diag[i].nodes)
					free(diag[i].nodes);
			}
			free(diag);
		}
	}
	bool append(unsigned i, unsigned j, const smatrix_dof & t) {
		// If we're appending to the diagonal, just do it.
		if (i == j) {
			diag[i].K += t;
			return true;
		}
		// Here we don't know if we're above or below the diagonal,
		// so we must use MIN/MAX.  We start by trying to find the
		// node which is in a row greater than or equal to the one
		// we're adding.
		smatrix_diag * d = &diag[MAX(i,j)];
		smatrix_node * n, * o, * p = d->col_head;
		while (p != 0 && p->i < MIN(i,j))
			p = p->col_next;
		// If we found a match, add to it.
		if (p != 0 && p->i == MIN(i,j)) {
			if (j < i)
				p->K += ~t;
			else
				p->K += t;
			return true;
		}
		// We've either run off the end, or we don't have a matching
		// node.  In this case we need to add one.
		if (!d->expand(1))
			return false;
		// Now fill in the new element.
		n = &(d->nodes[(d->nnz)++]);
	  	if (j < i) {
	  		n->K = ~t;
	  		swap(i,j);
	  	} else
	  		n->K = t;
	  	n->i = i;
		// Fix up the column references.
		p = d->col_head, o = 0;
		while (p != 0 && p->i < i)
			o = p, p = p->col_next;
		assert(p == 0 || p->i > i);
		if (o == 0)
			d->col_head = n;
		else
			o->col_next = n;
		n->col_next = p;
		return true;
	}
	// This just makes things faster by sorting the node lists.
	void tidy() {
		// We need somewhere to store the switched data.
		smatrix_node temp;
		
		for (unsigned i = 0; i < nnd; i++) {
			smatrix_node * nodes = diag[i].nodes;
			// Insertion sort the column of nodes.
			for (unsigned j = 1; j < diag[i].nnz; j++) {
				for (unsigned k = j; k > 0
						&& nodes[k-1].i > nodes[k].i; k--) {
					memcpy(&temp      ,&nodes[k-1],sizeof(smatrix_node));
					memcpy(&nodes[k-1],&nodes[k]  ,sizeof(smatrix_node));
					memcpy(&nodes[k]  ,&temp      ,sizeof(smatrix_node));
				}
			}
			// Now fix up the pointers.
			if (diag[i].col_head != 0)
				diag[i].col_head = nodes;
			for (unsigned short j = 0; j < diag[i].nnz; j++)
				nodes[j].col_next = (j == diag[i].nnz-1 ? 0 : &nodes[j+1]);
		}
	}
	// Incomplete Cholesky decompostion, by block.  Look at the small
	// internal version to understand what the whole thing is doing.
	void chol() {
		smatrix_diag * d;
		smatrix_node * p;

		// We want to do a level 1 fill, so we need to insert nodes.  We
		// rely on some tricks here...  We know that we don't use nnz, only
		// col_next, while tidy with only use nnz not col_next.  So we add
		// nodes ourself without using append.  We also know that we're
		// adding zeros, so we memset().
		/*for (unsigned n = 0; n < nnd; n++) {
			d = &(diag[n]);
			p = d->col_head;
			for (unsigned r = (p ? p->i : n); r < n; r++) {
				p = d->col_head;
				while (p && p->i < r)
					p = p->col_next;
				if (p && p->i == r)
					continue;
				smatrix_node * pi = diag[r].col_head;
				smatrix_node * pj = d->col_head;
				bool wouldadd = false;
				while (pi != 0 && pj != 0) {
					if (pi->i > pj->i)
						pj = pj->col_next;
					else if (pi->i < pj->i)
						pi = pi->col_next;
					else {
						wouldadd = true;
						break;
					}
				}
				if (true || wouldadd) {
					if (!d->expand(1))
						return; // oops.
					pj = &(d->nodes[(d->nnz)++]);
					memset(pj,0,sizeof(smatrix_node));
					pj->i = r;
				}
			}
		}*/
		// Now call tidy to magically add these nodes.
		tidy();

		// Do the decomposition.
		for (unsigned n = 0; n < nnd; n++) {
			d = &(diag[n]);
			p = d->col_head;
			while (p != 0) {
				smatrix_node * pi = diag[p->i].col_head;
				smatrix_node * pj = d->col_head;
				// We need to have blocks in both columns.  So skip the
				// missing ones.
				while (pi != 0 && pj != 0) {
					if (pi->i > pj->i)
						pj = pj->col_next;
					else if (pi->i < pj->i)
						pi = pi->col_next;
					else {
						p->K -= ~(pi->K)*(pj->K);
						pi = pi->col_next;
						pj = pj->col_next;
					}
				}
				// Do the block Cholesky decomposition for an off diagonal
				// element.
				for (unsigned i = 0; i < NDOF; i++) {
					for (unsigned j = 0; j < NDOF; j++) {
						for (unsigned k = 0; k < i; k++)
							p->K(i,j) -= diag[p->i].K(k,i)*p->K(k,j);
						p->K(i,j) *= diag[p->i].K(i,i);
					}
				}
				d->K -= ~(p->K)*(p->K);
				p = p->col_next;
			}
			// Cholesky decompostion of the diagonal.
			for (unsigned i = 0; i < NDOF; i++) {
				for (unsigned k = 0; k < i; k++)
					d->K(i,i) -= d->K(k,i)*d->K(k,i);
				assert(d->K(i,i) > 0.0);
				if (d->K(i,i) < DBL_EPSILON)
					printf("OOOOPS... (%g) %i/%i %i",d->K(i,i),n,nnd,i);
				// NOTE: We store the inverse of the diagonal to avoid
				// divisions.
				d->K(i,i) = 1.0/sqrt(d->K(i,i));
				for (unsigned j = i+1; j < NDOF; j++) {
					for (unsigned k = 0; k < i; k++)
						d->K(i,j) -= d->K(k,i)*d->K(k,j);
					d->K(i,j) *= d->K(i,i);
					d->K(j,i) = 0.0;
				}
			}
		}
	}

private:
	friend class mesh;
	unsigned nnd;
	smatrix_diag * diag;
};

/*
 * class svector - a vector to go with smatrix
 */
class svector {
public:
	inline explicit svector(const unsigned n)
	  : nnd(n), V(0) {
		V = static_cast<svector_dof *>
				(calloc(nnd,sizeof(svector_dof)));
		if (V == 0)
			event_msg(EVENT_ERROR,"Out of memory in svector::svector()!");
	}
	inline ~svector() {
		if (V)
			free(V);
	}
	inline const svector_dof & operator() (const unsigned i) const {
		return V[i];
	}
	inline svector_dof & operator() (const unsigned i) {
		return V[i];
	}

private:
	unsigned nnd;
	svector_dof * V;
};

/*
 * struct mesh_bc - a holder for mesh boundary conditions and forces.
 */
struct mesh_bc_key {
	unsigned n, i;
	mesh_bc_key()
	  : n(0), i(0) {
	}
	mesh_bc_key(unsigned N, unsigned I)
	  : n(N), i(I) {
	}
	// Provide these for sorting.
	bool operator== (const mesh_bc_key & k) const {
		return (n == k.n && i == k.i);
	}
	bool operator> (const mesh_bc_key & k) const {
		return (n > k.n || (n == k.n && i > k.i));
	}
	bool operator< (const mesh_bc_key & k) const {
		return (n < k.n || (n == k.n && i < k.i));
	}
};
struct mesh_bc : public mesh_bc_key {
	double d;
	mesh_bc()
	  : mesh_bc_key(), d(0.0) {
	}
	mesh_bc(unsigned N, unsigned I, double D)
	  : mesh_bc_key(N,I), d(D) {
	}
};

/*
 * class mesh - a 3D FEM mesh
 */
class mesh : private list_owned<mesh,element> {
public:
	// Create an empty mesh.
	explicit mesh()
	  : node(), disp_bc(), f_ext() {
	}
	~mesh() {
		// The linked list will delete all the elements.
	}

	// Add an element, returning a pointer if successful.
	const element * add(const element::element_t t, const material & m,
			const fset<coord3d> & c) {
		element * e = 0;
		switch (t) {
		case element::block8:
			e = new element_block<element::linear>(this,m,c);
			break;
		case element::block12:
			e = new element_block<element::quadratic>(this,m,c);
			break;
		case element::block16:
			e = new element_block<element::cubic>(this,m,c);
			break;
		case element::block18:
			e = new element_quad<element::linear>(this,m,c);
			break;
		case element::block27:
			e = new element_quad<element::quadratic>(this,m,c);
			break;
		case element::block36:
			e = new element_quad<element::cubic>(this,m,c);
			break;
		case element::infinite8:
			e = new element_infinite<element::linear>(this,m,c);
			break;
		case element::infinite12:
			e = new element_infinite<element::quadratic>(this,m,c);
			break;
		case element::infinite16:
			e = new element_infinite<element::cubic>(this,m,c);
			break;
		case element::variable18:
			e = new element_variable<element::linear>(this,m,c);
			if (e != 0 && e->nnd == 8) {
				delete e;
				e = new element_block<element::linear>(this,m,c);
			}
			break;
		case element::variable26:
			e = new element_variable<element::quadratic>(this,m,c);
			if (e != 0 && e->nnd == 12) {
				delete e;
				e = new element_block<element::quadratic>(this,m,c);
			}
			break;
		case element::variable34:
			e = new element_variable<element::cubic>(this,m,c);
			if (e != 0 && e->nnd == 16) {
				delete e;
				e = new element_block<element::cubic>(this,m,c);
			}
			break;
		case element::adaptor18:
			e = new element_adaptor<element::linear>(this,m,c);
			if (e != 0 && e->nnd == 8) {
				delete e;
				e = new element_block<element::linear>(this,m,c);
			}
			break;
		case element::adaptor26:
			e = new element_adaptor<element::quadratic>(this,m,c);
			if (e != 0 && e->nnd == 12) {
				delete e;
				e = new element_block<element::quadratic>(this,m,c);
			}
			break;
		case element::adaptor34:
			e = new element_adaptor<element::cubic>(this,m,c);
			if (e != 0 && e->nnd == 16) {
				delete e;
				e = new element_block<element::cubic>(this,m,c);
			}
			break;
		}
		if (e == 0)
			event_msg(EVENT_ERROR,"Out of memory in mesh::add()!");
		return e;
	}
	inline unsigned addnode(const coord3d & c) {
		unsigned k = node.haskey(c);
		if (k == UINT_MAX) {
			k = node.length();
			node.add(c);
		}
		return k;
	}
	inline bool updatenode(const node3d & n) {
		unsigned k = node.haskey(n);
		if (k == UINT_MAX)
			return false;
		node[k] = n;
		return true;
	}
	inline unsigned getnodes() const {
		return node.length();
	}
	inline const node3d & getnode(const unsigned i) const {
		return node[i];
	}
	inline const node3d & getorderednode(const unsigned i) const {
		return node.getindex(i);
	}
	inline unsigned hasnode(const coord3d & p) const {
		return node.haskey(p);
	}
	inline unsigned getorderofnode(const coord3d & p) const {
		return node.getorder(node.haskey(p));
	}

	// An enum to enable handling the degrees of freedom easily.
	enum dof {
		X = 0x0001,
		Y = 0x0002,
		Z = 0x0004
	};
	// An enum to enable adding planes of BCs.
	enum bcplane {
		at    = 0x0001,
		below = 0x0002,
		above = 0x0004
	};
	// Add a constraint at a single point.
	bool add_bc(const coord3d & p, const dof f, const double d) {
		unsigned k = node.haskey(p);
		if (k == UINT_MAX)
			return false;
		if (f & mesh::X)
			disp_bc.add(mesh_bc(k,0,d));
		if (f & mesh::Y)
			disp_bc.add(mesh_bc(k,1,d));
		if (f & mesh::Z)
			disp_bc.add(mesh_bc(k,2,d));
		return true;
	}
	// Add constrints along a plane.
	bool add_bc_plane(const dof o, const bcplane p, const double c,
			const dof f, const double d) {
		assert(!((o & mesh::X) && (o & mesh::Y)));
		assert(!((o & mesh::X) && (o & mesh::Z)));
		assert(!((o & mesh::Y) && (o & mesh::Z)));
		coord3d cc(c,c,c);
		for (unsigned i = 0; i < node.length(); i++) {
			const node3d & n = node[i];
			switch (o) {
			case mesh::X:
				if ((p & mesh::at) && (n.x == cc.x))
					add_bc(n,f,d);
				if ((p & mesh::below) && (n.x < cc.x))
					add_bc(n,f,d);
				if ((p & mesh::above) && (n.x > cc.x))
					add_bc(n,f,d);
				break;
			case mesh::Y:
				if ((p & mesh::at) && (n.y == cc.y))
					add_bc(n,f,d);
				if ((p & mesh::below) && (n.y < cc.y))
					add_bc(n,f,d);
				if ((p & mesh::above) && (n.y > cc.y))
					add_bc(n,f,d);
				break;
			case mesh::Z:
				if ((p & mesh::at) && (n.z == cc.z))
					add_bc(n,f,d);
				if ((p & mesh::below) && (n.z < cc.z))
					add_bc(n,f,d);
				if ((p & mesh::above) && (n.z > cc.z))
					add_bc(n,f,d);
				break;
			default:
				return false;
			}
		}
		return true;
	}
	// Add an external force at a node.
	bool add_fext(const coord3d & p, const dof f, const double d) {
		unsigned j, k = node.haskey(p);
		if (k == UINT_MAX)
			return false;
		if (f & mesh::X) {
			mesh_bc m(k,0,d);
			if ((j = f_ext.haskey(m)) != UINT_MAX)
				m.d += f_ext[j].d;
			f_ext.add(m);
		}
		if (f & mesh::Y) {
			mesh_bc m(k,1,d);
			if ((j = f_ext.haskey(m)) != UINT_MAX)
				m.d += f_ext[j].d;
			f_ext.add(m);
		}
		if (f & mesh::Z) {
			mesh_bc m(k,2,d);
			if ((j = f_ext.haskey(m)) != UINT_MAX)
				m.d += f_ext[j].d;
			f_ext.add(m);
		}
		return true;
	}
	// Finally, solve the system.
	bool solve(const double tol) {
		unsigned i, j, n, nnd = node.length();
		printf("Solving with %i nodes!",nnd);
		smatrix K(nnd);
		svector F(nnd), U(nnd), P(nnd), W(nnd), V(nnd), R(nnd);
		smatrix_diag * d;
		smatrix_node * p;
		const element * e = first;

		timeme("\nBuilding external forces...");
		for (i = 0; i < f_ext.length(); i++) {
			const mesh_bc & f = f_ext[i];
			//F(f.n)(f.i) = f.d;
			F(node.getorder(f.n))(f.i) = f.d;
		}

		timeme("\nSetting approximation of displacement...");
		for (i = 0; i < nnd; i++) {
			//node3d & u = node[i];
			node3d & u = node.getindex(i);
			U(i)(0) = u.ux; U(i)(1) = u.uy; U(i)(2) = u.uz;
		}

		timeme("\nSetting fixed BCs...");
		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			//n = u.n;
			n = node.getorder(u.n);
			// Each bit in fixed means that DOF is fixed.
			node[u.n].fixed |= (1 << u.i);
			// The last bit means at least one DOF is non-zero.
			if (u.d != 0.0) {
				node[u.n].fixed |= (1 << NDOF);
				P(n)(u.i) = u.d;
			}
		}

		timeme("\nAssembling Stiffness Matrix...");

		// Loop over the elements building the stiffness matrix.
		while (e) {
			smatrix_elem * ke = e->stiffness();
			if (ke == 0)
				return false;
			// Fix up the element matrix to account for displacement BCs
			for (i = 0; i < ke->nnd; i++) {
				unsigned gi = e->l2g(i);
				unsigned ni = node.getorder(gi);
				// At least one of P(ni) is non-zero, so apply the
				// fixup to F for this element.
				if (node[gi].fixed & (1 << NDOF)) {
					for (j = 0; j < ke->nnd; j++) {
						unsigned nj = node.getorder(e->l2g(j));
						if (i > j)
							F(nj) -=  (*ke)(j,i)*P(ni);
						else
							F(nj) -= ~(*ke)(i,j)*P(ni);
					}
				}
			}
			for (i = 0; i < ke->nnd; i++) {
				unsigned gi = e->l2g(i);
				// If this node is free or fixed, we're done.
				if ((node[gi].fixed & ((1<<NDOF)-1)) == 0
				 || (node[gi].fixed & ((1<<NDOF)-1)) == ((1<<NDOF)-1))
					continue;
				// Zero out the columns/rows which are fixed.
				for (j = 0; j < ke->nnd; j++) {
					for (unsigned k = 0; k < NDOF; k++) {
						if ((node[gi].fixed & (1<<k)) == 0)
							continue;
						for (unsigned l = 0; l < NDOF; l++) {
							if (i >= j)
								(*ke)(j,i)(l,k) = 0.0;
							if (i <= j)
								(*ke)(i,j)(k,l) = 0.0;
						}
					}
				}
			}
			// Now assemble the matrix into K.
			for (i = 0; i < ke->nnd; i++) {
				unsigned gi = e->l2g(i);
				unsigned ni = node.getorder(gi);
				// If this node is completly fixed, don't add anything.
				if ((node[gi].fixed & ((1<<NDOF)-1)) == ((1<<NDOF)-1))
					continue;
				for (j = i; j < ke->nnd; j++) {
					unsigned gj = e->l2g(j);
					unsigned nj = node.getorder(gj);
					if ((node[gj].fixed & ((1<<NDOF)-1)) == ((1<<NDOF)-1))
						continue;
					K.append(ni,nj,(*ke)(i,j));
				}
			}
			delete ke;
			e = e->next;
		}

		timeme("\nFinalising BCs...");

		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			//n = u.n;
			n = node.getorder(u.n);
			F(n)(u.i) = U(n)(u.i) = u.d;
			K.diag[n].K(u.i,u.i) = 1.0;
		}

		timeme("\nTidying up...");

		K.tidy();

		timeme("\nCopying...");

		smatrix M(K);

		timeme("\nComputing incomplete Cholesky...");

		M.chol();

		timeme("\nBeginning CG...");

		unsigned it = 0;
		double r = 0.0, ro = 0.0;
		for (i = 0; i < nnd; i++)
			r += tmatrix_scalar<double>(~F(i)*F(i));
		double ri = r, a;
		printf("with initial residual %g\n",ri);
		
		// CG solution
		while (true) {
			// Compute first approximation on iteration 0.
			for (i = 0; it == 0 && i < nnd; i++) {
				d = &(K.diag[i]);
				R(i) = F(i) - d->K*U(i);
				p = d->col_head;
				for (unsigned k = 0; k < d->nnz; k++, p++) { //while (p) {
					   R(i) -= ~(p->K)*U(p->i);
					R(p->i) -=  (p->K)*U(i);
					//p = p->col_next;
				}
			}
			for (i = 0, r = 0.0; i < nnd; i++)
				r += tmatrix_scalar<double>(~R(i)*R(i));
			printf("CG step %i with residual %g",it,r);
			if (r < tol*ri || it >= nnd) {
				timeme("\n");
				break;
			}
			// Forward substitution of IC preconditioner.
			for (n = 0; n < nnd; n++) {
				d = &(M.diag[n]);
				svector_dof t(R(n));
				p = d->col_head;
				for (unsigned k = 0; k < d->nnz; k++, p++) { //while (p) {
					t -= ~(p->K)*V(p->i);
					//p = p->col_next;
				}
				for (i = 0; i < NDOF; i++) {
					for (j = 0; j < i; j++)
						t(i) -= t(j)*d->K(j,i);
					t(i) *= d->K(i,i);
				}
				V(n) = t;
			}
			// Backward substitution of IC preconditioner.
			for (n = nnd; n > 0; n--) {
				d = &(M.diag[n-1]);
				svector_dof t(V(n-1));
				for (i = NDOF; i > 0; i--) {
					for (j = i; j < NDOF; j++)
						t(i-1) -= t(j)*d->K(i-1,j);
					t(i-1) *= d->K(i-1,i-1);
				}
				V(n-1) = t;
				p = d->col_head;
				for (unsigned k = 0; k < d->nnz; k++, p++) { //while (p) {
					V(p->i) -= (p->K)*t;
					//p = p->col_next;
				}
			}
			// Update search direction
			for (i = 0, r = 0.0; i < nnd; i++)
				r += tmatrix_scalar<double>(~V(i)*R(i));
			if (it == 0) {
				for (i = 0; i < nnd; i++)
					P(i) = V(i);
			} else {
				for (i = 0; i < nnd; i++)
					P(i) = V(i) + (r/ro)*P(i);
			}
			// Multiply by stiffness matrix to get step.
			for (i = 0; i < nnd; i++) {
				d = &(K.diag[i]);
				W(i) = d->K*P(i);
				p = d->col_head;
				for (unsigned k = 0; k < d->nnz; k++, p++) { //while (p) {
					   W(i) += ~(p->K)*P(p->i);
					W(p->i) +=  (p->K)*P(i);
					//p = p->col_next;
				}
			}
			for (i = 0, a = 0.0; i < nnd; i++)
				a += tmatrix_scalar<double>(~W(i)*P(i));
			assert(a > 0.0);
			ro = r;
			for (i = 0; i < nnd; i++) {
				U(i) += (ro/a)*P(i);
				R(i) -= (ro/a)*W(i);
			}
			// Get real error.
			/*for (i = 0; i < nnd; i++) {
				d = &(K.diag[i]);
				W(i) = F(i) - d->K*U(i);
				p = d->col_head;
				for (unsigned k = 0; k < d->nnz; k++, p++) {
					   W(i) -= ~(p->K)*U(p->i);
					W(p->i) -=  (p->K)*U(i);
				}
			}
			for (i = 0, a = 0.0; i < nnd; i++)
				a += tmatrix_scalar<double>(~W(i)*W(i));
			printf(" (%g)",a);*/
			it++;
			timeme("\n");
		}
		printf("CG took %i steps with residual %g\n",it,r);
		printf("Setting results...");

		for (i = 0; i < nnd; i++) {
			//node3d & u = node[i];
			node3d & u = node.getindex(i);
			u.setdisp(U(i));
		}

		timeme("\n");
		return true;
	}

private:
	friend class listelement_o<mesh,element>;
	
	node_list node;
	kiset<mesh_bc_key,mesh_bc> disp_bc;
	kiset<mesh_bc_key,mesh_bc> f_ext;
};

inline mesh::dof
operator| (const mesh::dof l, const mesh::dof r)
{
	return static_cast<mesh::dof>(unsigned(l) | unsigned(r));
}
inline mesh::dof
operator& (const mesh::dof l, const mesh::dof r)
{
	return static_cast<mesh::dof>(unsigned(l) & unsigned(r));
}

inline mesh::bcplane
operator| (const mesh::bcplane l, const mesh::bcplane r)
{
	return static_cast<mesh::bcplane>(unsigned(l) | unsigned(r));
}
inline mesh::bcplane
operator& (const mesh::bcplane l, const mesh::bcplane r)
{
	return static_cast<mesh::bcplane>(unsigned(l) & unsigned(r));
}

/*
 * Add a node to the mesh's node list, returning the index.
 */
unsigned
element::addnode(const coord3d & c) const
{
	return owner->addnode(c);
}

/*
 * Add a node to the mesh's node list, returning the index.
 */
void
element::updatenode(const node3d & n) const
{
	owner->updatenode(n);
}

/*
 * Get a node from the node list, based on the index.
 */
inline const node3d &
element::getnode(const unsigned i) const
{
	return owner->getnode(i);
}

/*
 * This program is a custom 3D finite element code, intended for
 * work on my PhD.  It will hopefully form the basis for a later
 * full 3D FEM pavement modelling code, but at the moment, I need
 * to get some work done...
 */
 
int run;
bool isvar;

/*
 * The basics of a layer used for the FEM.
 */
struct femlayer {
	fixed<8> top;
	fixed<8> bot;
	fixed<8> etop;
	fixed<8> ebot;
	material mat;
	
	femlayer()
	  : top(0.0), bot(0.0), etop(0.0), ebot(0.0), mat(0.0,0.0) {
	}
	femlayer(double t, double b, double e, double v)
	  : top(t), bot(b), etop(0.0), ebot(0.0), mat(e,v) {
	}
};

sset<femlayer> layer;
double * L[6];

double
node_depth_callback_real(const coord3d & c, const mesh * FEM, const material * mat)
{
	unsigned i = 0, n;
	
	while (&(layer[i].mat) != mat)
		i++;
	double z = -c.z - layer[i].top;
	double t = layer[i].top, b = layer[i].bot, h = b-t;
	if (isvar)
		n = FEM->getorderofnode(coord3d(c.x,c.y,0.0));
	else
		n = FEM->getorderofnode(coord3d(0.0,0.0,0.0));
	if (i == 0) {
		t += L[0][n];
		b += L[1][n];
	} else if (i == 1) {
		t += L[1][n];
	} else if (i == 2) {
		b += L[2][n];
	} else if (i == 3) {
		t += L[2][n];
	}
	return -(t+(b-t)*z/h);
}

double
node_emod_callback_real(const coord3d & c, const mesh * FEM, const material * mat)
{
	unsigned i = 0, n;
	
	while (&(layer[i].mat) != mat)
		i++;
	if (isvar)
		n = FEM->getorderofnode(coord3d(c.x,c.y,0.0));
	else
		n = FEM->getorderofnode(coord3d(0.0,0.0,0.0));
	return layer[i].mat.emod + L[(i>1?i-1:i)+3][n];
}

inline double
circlearea(double x, double y, double r)
{
	x = fabs(x); y = fabs(y);
	if (hypot(x,y) >= r)
		return 0.0;
	double t = M_PI_2 - asin(x/r) - asin(y/r);
	double xr = (sqrt(r*r-y*y)-x)*y/2;
	double yr = (sqrt(r*r-x*x)-y)*x/2;
	return (r*r*t/2) - xr - yr;
}

inline double
blockarea(double x1, double x2, double y1, double y2, double r)
{
	double h1 = hypot(x1,y1), h2 = hypot(x2,y1);
	double h3 = hypot(x1,y2), h4 = hypot(x2,y2);
	if (MIN(MIN(h1,h2),MIN(h3,h4)) >= r)
		return 0.0;
	if (MAX(MAX(h1,h2),MAX(h3,h4)) <= r)
		return fabs(x1-x2)*fabs(y1-y2);
	if (x1*x2 < 0.0) {
		double y = MAX(fabs(y1),fabs(y2));
		return 2*circlearea(0.0,y,r)-circlearea(x1,y,r)
					-circlearea(x2,y,r);
	}
	if (y1*y2 < 0.0) {
		double x = MAX(fabs(x1),fabs(x2));
		return 2*circlearea(x,0.0,r) - circlearea(x,y1,r)
					- circlearea(x,y2,r);
	}
	x1 = fabs(x1); x2 = fabs(x2); if (x2 < x1) swap(x1,x2);
	y1 = fabs(y1); y2 = fabs(y2); if (y2 < y1) swap(y1,y2);
	return circlearea(x1,y1,r) - circlearea(x1,y2,r)
			- circlearea(x2,y1,r) + circlearea(x2,y2,r);
}

/*
 * A region (in the horizontal plane) which has been or is being
 * added to the mesh.  It is described by the 'stops' where there
 * are mesh points.
 */
struct region {
	cset<fixed<8> > xstop;
	cset<fixed<8> > ystop;
	
	region() {
	}
	fixed<8> xm() {
		return xstop[0];
	}
	fixed<8> xp() {
		return xstop[xstop.length()-1];
	}
	fixed<8> ym() {
		return ystop[0];
	}
	fixed<8> yp() {
		return ystop[ystop.length()-1];
	}
	bool overlaps(double x1, double y1, double x2, double y2) {
		fixed<8> fxm(x1), fym(y1), fxp(x2), fyp(y2);
		return (fxm < xp() && fxp > xm() && fym < yp() && fyp > ym());
	}
	bool overlaps(region & r) {
		return (r.xm() < xp() && r.xp() > xm()
				 && r.ym() < yp() && r.yp() > ym());
	}
};

struct region_list : sset<region> {
	double xm() {
		fixed<8> x_m = 0;
		for (unsigned i = 0; i < length(); i++)
			x_m = MIN(x_m,(*this)[i].xm());
		return x_m;
	}
	double xp() {
		fixed<8> x_p = 0;
		for (unsigned i = 0; i < length(); i++)
			x_p = MAX(x_p,(*this)[i].xp());
		return x_p;
	}
	double ym() {
		fixed<8> y_m = 0;
		for (unsigned i = 0; i < length(); i++)
			y_m = MIN(y_m,(*this)[i].ym());
		return y_m;
	}
	double yp() {
		fixed<8> y_p = 0;
		for (unsigned i = 0; i < length(); i++)
			y_p = MAX(y_p,(*this)[i].yp());
		return y_p;
	}
	bool overlaps(double x1, double y1, double x2, double y2) {
		bool rv = false;
		for (unsigned i = 0; i < length(); i++)
			rv = rv || (*this)[i].overlaps(x1,y1,x2,y2);
		return rv;
	}
	void merge() {
		for (unsigned i = 0; i < length(); i++) {
			for (unsigned j = i+1; j < length(); j++) {
				region & r1 = (*this)[i];
				region & r2 = (*this)[j];
				if (r1.overlaps(r2)
				 || ((r1.xp() == r2.xm() || r1.xm() == r2.xp())
				  && (r1.ym() == r2.ym() && r1.yp() == r2.yp()))
				 || ((r1.yp() == r2.ym() || r1.ym() == r2.yp())
				  && (r1.xm() == r2.xm() && r1.xp() == r2.xp()))) {
					r1.xstop.add(r2.xstop); r1.xstop.sort();
					r1.ystop.add(r2.ystop); r1.ystop.sort();
					remove(j--);
				}
			}
		}
	}
};

/*
 * This outputs the results for all the elements in the list
 * Obviously, the list must be ones returned by add().
 */
void
results(const fset<const element *> & e, const fset<point3d> & c,
		const char * dname)
{
	fset<pavedata> data(c.length(),c.length());
	
	char fname[128];
	sprintf(fname,"%s_%s_%05d.bin",dname,(isvar?"3d":"1d"),run);
	FILE * f = fopen(fname,"wb");
	for (unsigned i = 0; i < e.length(); i++) {
		e[i]->results(c,data);
		for (unsigned j = 0; j < data.length(); j++) {
			fwrite(&(data[j].x),sizeof(double),3,f);
			fwrite(&(data[j].data[0][0]),sizeof(double),27,f);
		}
		
	}
	fclose(f);
}

void
LEresults(const fset<const element *> & e, const fset<point3d> & c,
		const char * dname, LEsystem & le)
{
	fset<pavedata> data(c.length(),c.length());
	
	le.removepoints();
	for (unsigned i = 0; i < e.length(); i++) {
		e[i]->results(c,data);
		unsigned il = 0;
		while (&(layer[il].mat) != &(e[i]->mat))
			il++;
		for (unsigned j = 0; j < data.length(); j++) {
			data[j].z *= -1;
			le.addpoint(data[j],il);
		}
	}
	le.calculate(LEsystem::all);

	char fname[128];
	sprintf(fname,"%s_%s_%05d.bin",dname,"3d",run);
	FILE * f = fopen(fname,"wb");
	for (unsigned i = 0; i < e.length(); i++) {
		e[i]->results(c,data);
		unsigned il = 0;
		while (&(layer[il].mat) != &(e[i]->mat))
			il++;
		for (unsigned j = 0; j < data.length(); j++) {
			data[j].z *= -1;
			pavedata d(le.result(data[j],il));
			d.z *= -1; d.data[4][2] *= -1;
			swap(d.data[3][0],d.data[3][2]);
			swap(d.data[8][0],d.data[8][2]);
			d.data[1][1] *= -1; d.data[1][2] *= -1;
			d.data[6][1] *= -1; d.data[6][2] *= -1;
			fwrite(&(d.x),sizeof(double),3,f);
			fwrite(&(d.data[0][0]),sizeof(double),27,f);
		}
	}
	fclose(f);
}

int
main()
{
	timeme();

	LEsystem test;
	test.addlayer(  90.0,3000e3,0.35);
	test.addlayer( 205.0, 750e3,0.35);
	test.addlayer( 205.0, 750e3,0.35);
	test.addlayer(1500.0, 100e3,0.35);

	layer.empty();
	for (unsigned i = 0; i < test.layers(); i++) {
		const LElayer & l = test.layer(i);
		layer.add(femlayer(l.top(),l.bottom(),l.emod(),l.poissons()));
	}

	double x = 0.0, y = 0.0, z = 0.0, zmax = 0.0, rmax = 0.0;
	double dx, dy, dz, delta = 4;
	const double edge = 4096;
	unsigned step = 4;
	mesh FEM;
	fset<coord3d> coord(8);
	cset<const element *> along;
	cset<const element *> trans;
	cset<const element *> sf, ac, ab, sg;
	sset<fixed<8> > stop;
	sset<unsigned> level;
	region_list filled, filling;

	//test.addload(point2d(0.0,0.0),0.0,690.0,100.0);
	test.addload(point2d(-160.0,0.0),0.0,690.0,100.0);
	test.addload(point2d( 160.0,0.0),0.0,690.0,100.0);

	// Find the step size, which depends on our tires.
	for (unsigned i = 0; i < test.loads(); i++) {
		x += test.getload(i).x;
		y += test.getload(i).y;
		double r = test.getload(i).radius();
		if (r > rmax)
			rmax = r;
	}
	x /= test.loads(); y /= test.loads();
	assert(x == 0.0);
	assert(y == 0.0);
	while ((step-2)*delta < rmax)
		step *= 2;

	// Start with the tire grid.
	dx = delta; dy = delta;
	for (unsigned i = 0; i < test.loads(); i++) {
		double lx = test.getload(i).x;
		double ly = test.getload(i).y;
		lx = 2*dx*(lx < 0 ? floor(lx/dx/2) : ceil(lx/dx/2));
		ly = 2*dy*(ly < 0 ? floor(ly/dy/2) : ceil(ly/dy/2));
		double lr = test.getload(i).radius();
		double F = -test.getload(i).pressure();
		for (x = lx - (step-2)*dx; x < lx + (step-2)*dx; x += 2*dx) {
			for (y = ly -(step-2)*dy; y < ly + (step-2)*dy; y += 2*dy) {
				double ba = blockarea(x-lx,x+2*dx-lx,y-ly,y+2*dx-ly,lr);
				if (ba > 0.0) {
					FEM.addnode(coord3d(x     ,y     ,0.0));
					FEM.addnode(coord3d(x  +dx,y     ,0.0));
					FEM.addnode(coord3d(x+2*dx,y     ,0.0));
					FEM.addnode(coord3d(x     ,y  +dy,0.0));
					FEM.addnode(coord3d(x  +dx,y  +dy,0.0));
					FEM.addnode(coord3d(x+2*dx,y  +dy,0.0));
					FEM.addnode(coord3d(x     ,y+2*dy,0.0));
					FEM.addnode(coord3d(x  +dx,y+2*dy,0.0));
					FEM.addnode(coord3d(x+2*dx,y+2*dy,0.0));
				}
			}
		}
		// Add the tire neighbours
		for (x = lx - (step-2)*dx; x <= lx + (step-2)*dx; x += dx) {
			for (y = ly - (step-2)*dy; y <= ly + (step-2)*dy; y += dy) {
				unsigned p = FEM.hasnode(coord3d(x,y,0.0));
				if (p == UINT_MAX)
					continue;
				node3d n = FEM.getnode(p);
				unsigned x_m = FEM.hasnode(coord3d(x-dx,y,0.0));
				if (x_m == UINT_MAX && hypot(x-dx-lx,y-ly) < lr)
					x_m = FEM.hasnode(coord3d(x-2*dx,y,0.0));
				unsigned x_p = FEM.hasnode(coord3d(x+dx,y,0.0));
				if (x_p == UINT_MAX && hypot(x+dx-lx,y-ly) < lr)
					x_p = FEM.hasnode(coord3d(x+2*dx,y,0.0));
				unsigned y_m = FEM.hasnode(coord3d(x,y-dy,0.0));
				if (y_m == UINT_MAX && hypot(x-lx,y-dy-ly) < lr)
					y_m = FEM.hasnode(coord3d(x,y-2*dy,0.0));
				unsigned y_p = FEM.hasnode(coord3d(x,y+dy,0.0));
				if (y_p == UINT_MAX && hypot(x-lx,y+dy-ly) < lr)
					y_p = FEM.hasnode(coord3d(x,y+2*dy,0.0));
				n.setneighbours(x_m,x_p,y_m,y_p,UINT_MAX,UINT_MAX);
				FEM.updatenode(n);
			}
		}
		// Add the tire loads
		for (x = lx - (step-2)*dx; x <= lx + (step-2)*dx; x += dx) {
			for (y = ly - (step-2)*dy; y <= ly + (step-2)*dy; y += dy) {
				unsigned p = FEM.hasnode(coord3d(x,y,0.0));
				if (p == UINT_MAX)
					continue;
				node3d n = FEM.getnode(p);
				double fxm = 0.0, fxp = 0.0, fym = 0.0, fyp = 0.0;
				if (n.xm != UINT_MAX)
					fxm = (x - double(FEM.getnode(n.xm).x))/2;
				if (n.xp != UINT_MAX)
					fxp = (double(FEM.getnode(n.xp).x) - x)/2;
				if (n.ym != UINT_MAX)
					fym = (y - double(FEM.getnode(n.ym).y))/2;
				if (n.yp != UINT_MAX)
					fyp = (double(FEM.getnode(n.yp).y) - y)/2;
				double f = F*blockarea(x-fxm-lx,x+fxp-lx,y-fym-ly,y+fyp-ly,lr);
				if (fabs(f) > 0.0)
					FEM.add_fext(n,mesh::Z,f);
			}
		}
		region r;
		for (x = lx - step*dx; x <= lx + step*dx; x += 2*dx)
			r.xstop.add(x);
		for (y = ly - step*dy; y <= ly + step*dy; y += 2*dy)
			r.ystop.add(y);
		r.xstop.sort(); r.ystop.sort();
		filling.add(r);
	}
	filling.merge();
	printf("Loads: [%f %f] x [%f %f]\n",filling.xm(),filling.xp(),
			filling.ym(),filling.yp());
	for (unsigned i = 0; i < filling.length(); i++) {
		region & r = filling[i];
		printf(" %i: [%f %f] x [%f %f]\n",i,double(r.xm()),double(r.xp()),double(r.ym()),double(r.yp()));
	}
	// Make step power of two.
	while ((step & (step - 1)) != 0)
		step++;
	step /= 2; delta *= 2;

	// Find the equivalent layer top and bottom.
	// XXX: This relies on decreasing emod with depth.
	for (unsigned i = 0; i < layer.length(); i++) {
		double h = layer[i].bot - layer[i].top, he = h;
		layer[i].etop = (i == 0 ? 0.0 : double(layer[i-1].ebot));
		for (unsigned j = i+1; j < layer.length(); j++) {
			double v = layer[j-1].mat.v;
			double E = layer[j-1].mat.emod;
			double v1 = layer[j].mat.v;
			double E1 = layer[j].mat.emod;
			he *= pow(E/E1*(1+v1*v1)/(1+v*v),1.0/3.0);
		}
		layer[i].ebot = double(layer[i].etop) + he;
	}
	unsigned zstep = 1, zi = 0;
	// add the depth stops.
	dz = rmax/zstep;
	zmax = layer[layer.length()-1].ebot; z = layer[0].etop;
	stop.add(z); level.add(zi);
	while (z < zmax) {
		unsigned zscale = 4;
		for (; z < zscale*zstep*dz; z += dz) {
			stop.add(z+dz); level.add(zi);
		}
		dz *= 2; zi++;
	}
	for (unsigned j = 0; j < stop.length(); j++)
		stop[j] = double(stop[j])*zmax/z;
	// now scale them to the layers
	for (unsigned i = 0, t = 0, b; i < layer.length(); i++) {
		femlayer & l = layer[i];
		assert(t < stop.length());
		assert(stop[t] == l.etop);
		b = t;
		while (b < stop.length() && stop[b] < l.ebot)
			b++;
		assert(b < stop.length());
		// Find the closest stop.
		if (fabs(stop[b-1]-l.ebot) < fabs(l.ebot-stop[b]))
			b--;
		if (b == t)
			b = t+1;
		unsigned s = 1;
		//if (i == 0)
		//	s = 3; // Make sure we have three nodes in the top layer.
		while (b < stop.length() && (b-t)%s > 0)
			b++;
		if (b == stop.length()) {
			b = stop.length() - 1;
			s = b-t;
		} else
			s--;
		for (unsigned j = t+1; j < b; j++, b -= s) {
			for (unsigned k = 0; k < s; k++) {
				stop.remove(j); level.remove(j);
			}
		}
		// Now scale these stop to match the layer.
		double hl = l.bot - l.top, hs = stop[b] - stop[t];
		for (unsigned j = t; j < b; j++)
			stop[j] = double(l.bot) - double(stop[b]-stop[j])*hl/hs;
		t = b; b = stop.length() - 1;
		assert(double(stop[b]) == zmax);
		hl = zmax - double(l.ebot); hs = stop[b] - stop[t];
		for (unsigned j = t; j < b; j++)
			stop[j] = zmax - double(stop[b]-stop[j])*hl/hs;
	}
	stop[stop.length()-1] = zmax = layer[layer.length()-1].bot;

	// Now add the elements from below the tire, working outwards.
	zi = 0;
	while (zi <= level[level.length()-1] || filled.xp() < edge) {
		for (unsigned i = 0, j = 1; j < stop.length() && level[j] <= zi; j++) {
			while (i < layer.length() && stop[j] > layer[i].bot)
				i++;
			assert(i < layer.length());
			double z1 = stop[j  ];
			double z2 = stop[j-1];
			printf("%.0f\t%.0f\t%f\t%i\n",filling.xp(),filled.xp(),z1,i);
			element::element_t var = (i == 0 ?
					element::variable34 : element::variable26);
			element::element_t inf = (var == element::variable34 ?
					element::infinite16 : var == element::variable26 ?
					element::infinite12 : element::infinite8);
			for (unsigned k = 0; k < filling.length(); k++) {
				region & r = filling[k];
				for (unsigned xi = 1; xi < r.xstop.length(); xi++) {
					for (unsigned yi = 1; yi < r.ystop.length(); yi++) {
						double x1 = double(r.xstop[xi-1]);
						double x2 = double(r.xstop[xi]);
						double y1 = double(r.ystop[yi-1]);
						double y2 = double(r.ystop[yi]);
						if (filled.overlaps(x1,y1,x2,y2) && level[j] < zi)
							continue;
						coord.resize(8);
						coord[0] = coord3d(x1,y1,-z1);
						coord[1] = coord3d(x1,y2,-z1);
						coord[2] = coord3d(x2,y1,-z1);
						coord[3] = coord3d(x2,y2,-z1);
						coord[4] = coord3d(x1,y1,-z2);
						coord[5] = coord3d(x1,y2,-z2);
						coord[6] = coord3d(x2,y1,-z2);
						coord[7] = coord3d(x2,y2,-z2);
						const element * e = FEM.add(var,layer[i].mat,coord);
						if (x1 == 0.0)
							along.add(e);
						if (y1 == 0.0)
							trans.add(e);
						if (stop[j-1] == layer[0].top)
							sf.add(e);
						if (stop[j] == layer[0].bot)
							ac.add(e);
						if (stop[j] == layer[1].bot)
							ab.add(e);
						if (stop[j-1] == layer[3].top)
							sg.add(e);
						if (x1 == -edge) {
							if (y1 == -edge) {
								coord.resize(2);
								coord[0] = coord3d(x1,y1,-z1);
								coord[1] = coord3d(x1,y1,-z2);
								FEM.add(inf,layer[i].mat,coord);
							} if (y2 == edge) {
								coord.resize(2);
								coord[0] = coord3d(x1,y2,-z1);
								coord[1] = coord3d(x1,y2,-z2);
								FEM.add(inf,layer[i].mat,coord);
							}
							coord.resize(4);
							coord[0] = coord3d(x1,y1,-z1);
							coord[1] = coord3d(x1,y2,-z1);
							coord[2] = coord3d(x1,y1,-z2);
							coord[3] = coord3d(x1,y2,-z2);
							FEM.add(inf,layer[i].mat,coord);
						} else if (x2 == edge) {
							if (y1 == -edge) {
								coord.resize(2);
								coord[0] = coord3d(x2,y1,-z1);
								coord[1] = coord3d(x2,y1,-z2);
								FEM.add(inf,layer[i].mat,coord);
							} if (y2 == edge) {
								coord.resize(2);
								coord[0] = coord3d(x2,y2,-z1);
								coord[1] = coord3d(x2,y2,-z2);
								FEM.add(inf,layer[i].mat,coord);
							}
							coord.resize(4);
							coord[0] = coord3d(x2,y1,-z1);
							coord[1] = coord3d(x2,y2,-z1);
							coord[2] = coord3d(x2,y1,-z2);
							coord[3] = coord3d(x2,y2,-z2);
							FEM.add(inf,layer[i].mat,coord);
						}
						if (y1 == -edge) {
							coord.resize(4);
							coord[0] = coord3d(x1,y1,-z1);
							coord[1] = coord3d(x2,y1,-z1);
							coord[2] = coord3d(x1,y1,-z2);
							coord[3] = coord3d(x2,y1,-z2);
							FEM.add(inf,layer[i].mat,coord);
						} if (y2 == edge) {
							coord.resize(4);
							coord[0] = coord3d(x1,y2,-z1);
							coord[1] = coord3d(x2,y2,-z1);
							coord[2] = coord3d(x1,y2,-z2);
							coord[3] = coord3d(x2,y2,-z2);
							FEM.add(inf,layer[i].mat,coord);
						}
					}
				}
			}
		}
	printf("Just: [%f %f] x [%f %f]\n",filling.xm(),filling.xp(),
			filling.ym(),filling.yp());
	for (unsigned i = 0; i < filling.length(); i++) {
		region & r = filling[i];
		printf(" %i: [%f %f] x [%f %f]\n",i,double(r.xm()),double(r.xp()),double(r.ym()),double(r.yp()));
	}
		if (filling.xm() > -edge && filling.xp() < edge
		 && filling.ym() > -edge && filling.yp() < edge) {
			delta *= 2;
			for (unsigned i = 0; i < filling.length(); i++) {
				region & r = filling[i];
				assert(r.xstop.length()%2 == 1);
				assert(r.ystop.length()%2 == 1);
				for (unsigned xi = 1; xi < r.xstop.length(); xi++)
					r.xstop.remove(xi);
				for (unsigned yi = 1; yi < r.ystop.length(); yi++)
					r.ystop.remove(yi);
			}
		}
	printf("Filled: [%f %f] x [%f %f]\n",filling.xm(),filling.xp(),
			filling.ym(),filling.yp());
	for (unsigned i = 0; i < filling.length(); i++) {
		region & r = filling[i];
		printf(" %i: [%f %f] x [%f %f]\n",i,double(r.xm()),double(r.xp()),double(r.ym()),double(r.yp()));
	}
		filled = filling; dx = delta; dy = delta;
		for (unsigned i = 0; i < filling.length(); i++) {
			region & r = filling[i];
			double rx = double(r.xm() + r.xp())/2;
			double ry = double(r.ym() + r.yp())/2;
			rx = 2*dx*(rx < 0 ? floor(rx/dx/2) : ceil(rx/dx/2));
			ry = 2*dy*(ry < 0 ? floor(ry/dy/2) : ceil(ry/dy/2));
			assert(double(r.xm()) == -edge || double(r.xm()) > rx-step*dx);
			assert(double(r.xp()) ==  edge || double(r.xp()) < rx+step*dx);
			assert(double(r.ym()) == -edge || double(r.ym()) > ry-step*dy);
			assert(double(r.yp()) ==  edge || double(r.yp()) < ry+step*dy);
			for (x = r.xm(); x > MAX(rx-step*dx,-edge); x -= dx)
				r.xstop.add(MAX(x-dx,-edge));
			r.xstop.sort();
			for (x = r.xp(); x < MIN(rx+step*dx, edge); x += dx)
				r.xstop.add(MIN(x+dx, edge));
			r.xstop.sort();
			for (y = r.ym(); y > MAX(ry-step*dy,-edge); y -= dy)
				r.ystop.add(MAX(y-dy,-edge));
			r.ystop.sort();
			for (y = r.yp(); y < MIN(ry+step*dy, edge); y += dy)
				r.ystop.add(MIN(y+dy, edge));
			r.ystop.sort();
		}
		filling.merge();
	printf("Filling: [%f %f] x [%f %f]\n",filling.xm(),filling.xp(),
			filling.ym(),filling.yp());
	for (unsigned i = 0; i < filling.length(); i++) {
		region & r = filling[i];
		printf(" %i: [%f %f] x [%f %f]\n",i,double(r.xm()),double(r.xp()),double(r.ym()),double(r.yp()));
	}
		zi++;
	}
	//FEM.add_bc_plane(mesh::X,mesh::at|mesh::below,-edge,
	//		mesh::X,0.0);
	//FEM.add_bc_plane(mesh::X,mesh::at|mesh::above, edge,
	//		mesh::X,0.0);
	//FEM.add_bc_plane(mesh::Y,mesh::at|mesh::below,-edge,
	//		mesh::Y,0.0);
	//FEM.add_bc_plane(mesh::Y,mesh::at|mesh::above, edge,
	//		mesh::Y,0.0);
	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,-zmax,
			mesh::X|mesh::Y|mesh::Z,0.0);
	
	// Generate our random fields.
	unsigned np = 0;
	while (double(FEM.getorderednode(np).z) == 0.0)
		np++;
	for (unsigned l = 0; l < 6; l++) {
		L[l] = new double[np];
		/*double * C = new double[T_SIZE(np)];
		if (L[l] == 0 || C == 0) {
			event_msg(EVENT_ERROR,"Ooops, out of memory...");
			return 0;
		}

		for (unsigned i = 0; i < np; i++) {
			const node3d & p1 = FEM.getorderednode(i);
			for (unsigned j = i; j < np; j++) {
				const node3d & p2 = FEM.getorderednode(j);
				double r = fabs(hypot(p1.x-p2.x,p1.y-p2.y));
				C[T_IDX(i,j)] = (r > 6000.0 ? 0.0
					: 5.0*(1-1.5*r/6000.0+0.5*pow(r/6000.0,3)));
			}
		}
		decmp_chol_tri(np,C);

		char fname[128];
		sprintf(fname,"cor%d.tri",l);
		FILE * f = fopen(fname,"wb");
		fwrite(C,sizeof(double),T_SIZE(np),f);
		fclose(f);

		delete C;*/
	}
		
	run = 1;
	isvar = false;
	rng RNG;
	while (true) {
		//isvar = !isvar;

		for (unsigned l = 0; isvar && l < 6; l++) {
			double * A = new double[np];
			double * C = new double[T_SIZE(np)];
			if (A == 0 || C == 0) {
				event_msg(EVENT_ERROR,"Ooops, out of memory...");
				return 0;
			}
			
			char fname[128];
			sprintf(fname,"cor%d.tri",l);
			FILE * f = fopen(fname,"rb");
			fread(C,sizeof(double),T_SIZE(np),f);
			fclose(f);		

			for (unsigned i = 0; i < np; i++) {
				A[i] = RNG.stdnormal();
				L[l][i] = C[T_IDX(i,i)]*A[i];
				for (unsigned j = 0; j < i; j++)
					L[l][i] += C[T_IDX(j,i)]*A[j];
			}
			delete A;
			delete C;
		}

		for (unsigned i = 0; i < FEM.getnodes(); i++) {
			node3d n = FEM.getnode(i);
			n.ux = 0.0; n.uy = 0.0; n.uz = 0.0;
			FEM.updatenode(n);
		}
		FEM.solve(1e-33);

		fset<point3d> poly(4), face(8);
		face[0] = point3d(-1,-1,    -1);
		face[1] = point3d(-1,-1,-1/3.0);
		face[2] = point3d(-1,-1, 1/3.0);
		face[3] = point3d(-1,-1,     1);
		face[4] = point3d( 1,-1,     1);
		face[5] = point3d( 1,-1, 1/3.0);
		face[6] = point3d( 1,-1,-1/3.0);
		face[7] = point3d( 1,-1,    -1);
		results(trans,face,"face");
		LEresults(trans,face,"face",test);
		face[0] = point3d(-1,-1,    -1);
		face[1] = point3d(-1,-1,-1/3.0);
		face[2] = point3d(-1,-1, 1/3.0);
		face[3] = point3d(-1,-1,     1);
		face[4] = point3d(-1, 1,     1);
		face[5] = point3d(-1, 1, 1/3.0);
		face[6] = point3d(-1, 1,-1/3.0);
		face[7] = point3d(-1, 1,    -1);
		results(along,face,"long");
		LEresults(along,face,"long",test);
		poly[0] = point3d(-1,-1,-1);
		poly[1] = point3d(-1, 1,-1);
		poly[2] = point3d( 1, 1,-1);
		poly[3] = point3d( 1,-1,-1);
		results(ac,poly,"ac");
		LEresults(ac,poly,"ac",test);
		results(ab,poly,"ab");
		LEresults(ab,poly,"ab",test);
		poly[0] = point3d(-1,-1, 1);
		poly[1] = point3d(-1, 1, 1);
		poly[2] = point3d( 1, 1, 1);
		poly[3] = point3d( 1,-1, 1);
		results(sf,poly,"sf");
		LEresults(sf,poly,"sf",test);
		results(sg,poly,"sg");
		LEresults(sg,poly,"sg",test);

		//if (isvar)
		//	continue;
		//if (++run > 300)
			break;
	}

	/*
	unsigned i, nnd = FEM.getnodes();
	test.removepoints();
	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z; // XXX
		if (!(x == 0 || y == 0 || z == 0))
			continue;
		//if (!(x == 0))
		//	continue;
		if (x < -edge || x > edge || y < -edge || y > edge || z <= -zmax)
			continue;
		test.addpoint(point3d(x,y,-z));
	}
	test.calculate(LEsystem::all);

	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z; // XXX
		if (!(x == 0 || y == 0 || z == 0))
			continue;
		//if (!(x == 0))
		//	continue;
		if (x < -edge || x > edge || y < -edge || y > edge || z <= -zmax)
			continue;
		const pavedata & d = test.result(point3d(x,y,-z));
		double ux = n.ux;
		double uy = n.uy;
		double uz = n.uz;
		double vx =  d.result(pavedata::deflct,pavedata::xx);
		double vy =  d.result(pavedata::deflct,pavedata::yy);
		double vz = -d.result(pavedata::deflct,pavedata::zz);
		unsigned j = FEM.hasnode(n);
		double h = hypot(hypot(vx-ux,vy-uy),vz-uz);
		double v = hypot(hypot(vx,vy),vz);
		printf("Node %6i: (%+6i,%+6i,%+6i) =\t(%8.2g,%8.2g,%8.2g)\t(%8.2g,%8.2g,%8.2g)\t%8.2g\t(%8.2g)\n",j,int(x),int(y),int(z),vx,vy,vz,ux,uy,uz,h,(v == 0.0 ? 0.0 : h/v));
	}*/

	timeme("\n");
	return 0;
}

double
node_depth_callback(const coord3d & c, const mesh * FEM, const material * mat)
{
	return c.z;
}

double
node_emod_callback(const coord3d & c, const mesh * FEM, const material * mat)
{
	return mat->emod;
}

int
main_test()
{
	timeme();
	printf("Constructing Mesh...");

	material m(1000,0.2);

	const double cube[3][2] = {{-10, 10}, {-10, 10}, {-10, 0}};
	const unsigned ndiv[3] = {10, 10, 3};
	unsigned i, j, k;
	mesh FEM;
	cset<const element *> face;

	double dx = (cube[0][1]-cube[0][0])/ndiv[0];
	double dy = (cube[1][1]-cube[1][0])/ndiv[1];
	double dz = (cube[2][1]-cube[2][0])/ndiv[2];
	fset<coord3d> coord(8);

	//FEM.addnode(coord3d(cube[0][0],(cube[1][1]-cube[1][0])/2,cube[2][1]));
	//FEM.addnode(coord3d(cube[0][1],(cube[1][1]-cube[1][0])/2,cube[2][1]));
	//FEM.addnode(coord3d((cube[0][1]-cube[0][0])/2,cube[1][0],cube[2][1]));
	//FEM.addnode(coord3d((cube[1][1]-cube[1][0])/2,cube[1][1],cube[2][1]));
	//FEM.addnode(coord3d((cube[1][1]-cube[1][0])/2,(cube[1][1]-cube[1][0])/2,cube[2][1]));

	for (k = 0; k < ndiv[2]; k++) {
		double z = cube[2][1] - k*dz;
		int p = 1; //(1 << k); // 1; 
		for (j = 0; j < ndiv[1]; j += p) {
			double y = cube[1][0] + j*dy;
			for (i = 0; i < ndiv[0]; i += p) {
				double x = cube[0][0] + i*dx;
				coord[0] = coord3d(x     ,y     ,z-dz);
				coord[1] = coord3d(x     ,y+p*dy,z-dz);
				coord[2] = coord3d(x+p*dx,y     ,z-dz);
				coord[3] = coord3d(x+p*dx,y+p*dy,z-dz);
				coord[4] = coord3d(x     ,y     ,z   );
				coord[5] = coord3d(x     ,y+p*dy,z   );
				coord[6] = coord3d(x+p*dx,y     ,z   );
				coord[7] = coord3d(x+p*dx,y+p*dy,z   );
				//printf("%4.2f\t%4.2f\t%4.2f\n",x,y,z);
				element::element_t s = (k == 0 ? element::block18 :
					element::adaptor18);
				const element * e = FEM.add(s,m,coord);
				if (i == ndiv[0]/2)
					face.add(e);
			}
		}
	}
	
	timeme("\nSetting BCs and forces...");

	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,cube[2][0],mesh::Z,0.0);
	FEM.add_bc(coord3d(cube[0][0],cube[1][0],cube[2][0]),mesh::X,0.0);
	FEM.add_bc(coord3d(cube[0][0],cube[1][0],cube[2][0]),mesh::Y,0.0);
	FEM.add_bc(coord3d(cube[0][1],cube[1][0],cube[2][0]),mesh::Y,0.0);
	const double F = -0.15;
	for (i = 0; i < FEM.getnodes(); i++) {
		const node3d & n = FEM.getnode(i);
		if (double(n.z) < cube[2][1])
			continue;
		double x = n.x, y = n.y, sx = 1.0, sy = 1.0, s = 3.0;
		double fxm = 0.0, fxp = 0.0, fym = 0.0, fyp = 0.0;
		if (ROUND(2*fmod(x-cube[0][0],dx)) > 0)
			sx = 2.0;
		if (ROUND(2*fmod(y-cube[1][0],dy)) > 0)
			sy = 2.0;
		if (n.xm != UINT_MAX)
			fxm = sx*(x - double(FEM.getnode(n.xm).x))/s;
		if (n.xp != UINT_MAX)
			fxp = sx*(double(FEM.getnode(n.xp).x) - x)/s;
		if (n.ym != UINT_MAX)
			fym = sy*(y - double(FEM.getnode(n.ym).y))/s;
		if (n.yp != UINT_MAX)
			fyp = sy*(double(FEM.getnode(n.yp).y) - y)/s;
		double f = F*(fxm+fxp)*(fym+fyp);
		FEM.add_fext(n,mesh::Z,f);
	}

	//for (i = 0; i < FEM.getnodes(); i++) {
	//	node3d n = FEM.getorderednode(i);
	//	n.ux =  0.0006*(double(n.x)+10.0)/20.0;
	//	n.uy =  0.0006*(double(n.y)+10.0)/20.0;
	//	n.uz = -0.0015*(double(n.z)+10.0)/10.0;
	//	//n.ux = 0.0; n.uy = 0.0; n.uz = 0.0;
	//	FEM.updatenode(n);
	//}
	FEM.solve(1e-25);

	timeme("\nOutput...\n");

	for (i = 0; i < FEM.getnodes(); i++) {
		const node3d & n = FEM.getorderednode(i);
		double x = n.x;
		double y = n.y;
		double z = n.z;
		//if (y != 0.0 || x != 0.0)
		//	continue;
		double ux = n.ux;
		double uy = n.uy;
		double uz = n.uz;
		j = FEM.hasnode(n);
		printf("Node %3i: (%6.2f,%6.2f,%6.2f) = (%+f,%+f,%+f)\n",j,double(x),double(y),double(z),ux,uy,uz);
	}
	//const node3d & n = FEM.getnode(FEM.hasnode(coord3d(cube[0][1],cube[1][1],cube[2][1])));
	//printf("(%+f,%+f,%+f)\n",n.ux,n.uy,n.uz);

	/*fset<point3d> poly(4,4);
	poly[0] = point3d(-1,-1,-1);
	poly[1] = point3d(-1,-1,+1);
	poly[2] = point3d(-1,+1,+1);
	poly[3] = point3d(-1,+1,-1);
	results(face,poly,"test");*/

	timeme(" Done.\n");

	return 0;
}

int
main_infinity()
{
	timeme();
	
	material m(100,0.35);
	mesh FEM;

	sset<coord3d> coord(0,8);

	coord.empty();
	coord.add(coord3d(-10,-10,-10));
	coord.add(coord3d(-10,  0,-10));
	coord.add(coord3d(  0,-10,-10));
	coord.add(coord3d(  0,  0,-10));
	coord.add(coord3d(-10,-10,  0));
	coord.add(coord3d(-10,  0,  0));
	coord.add(coord3d(  0,-10,  0));
	coord.add(coord3d(  0,  0,  0));
	FEM.add(element::block8,m,coord);
	for (int i = 0; i < 8; i++)
		coord[i].x *= -1;
	FEM.add(element::block8,m,coord);
	for (int i = 0; i < 8; i++)
		coord[i].y *= -1;
	FEM.add(element::block8,m,coord);
	for (int i = 0; i < 8; i++)
		coord[i].x *= -1;
	FEM.add(element::block8,m,coord);

	coord.empty();
	coord.add(coord3d(-10,-10,-10));
	coord.add(coord3d(-10,-10,  0));
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 2; i++)
		coord[i].x *= -1;
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 2; i++)
		coord[i].y *= -1;
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 2; i++)
		coord[i].x *= -1;
	FEM.add(element::infinite8,m,coord);

	coord.empty();
	coord.add(coord3d( 10,  0,-10));
	coord.add(coord3d( 10, 10,-10));
	coord.add(coord3d( 10,  0,  0));
	coord.add(coord3d( 10, 10,  0));
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 4; i++)
		coord[i].x *= -1;
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 4; i++)
		coord[i].y *= -1;
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 4; i++)
		coord[i].x *= -1;
	FEM.add(element::infinite8,m,coord);

	coord.empty();
	coord.add(coord3d(  0, 10,-10));
	coord.add(coord3d( 10, 10,-10));
	coord.add(coord3d(  0, 10,  0));
	coord.add(coord3d( 10, 10,  0));
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 4; i++)
		coord[i].x *= -1;
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 4; i++)
		coord[i].y *= -1;
	FEM.add(element::infinite8,m,coord);
	for (int i = 0; i < 4; i++)
		coord[i].x *= -1;
	FEM.add(element::infinite8,m,coord);

	FEM.add_fext(coord3d(-10,-10,  0),mesh::Z, -2.5);
	FEM.add_fext(coord3d(-10,  0,  0),mesh::Z, -5.0);
	FEM.add_fext(coord3d(-10, 10,  0),mesh::Z, -2.5);
	FEM.add_fext(coord3d(  0,-10,  0),mesh::Z, -5.0);
	FEM.add_fext(coord3d(  0,  0,  0),mesh::Z,-10.0);
	FEM.add_fext(coord3d(  0, 10,  0),mesh::Z, -5.0);
	FEM.add_fext(coord3d( 10,-10,  0),mesh::Z, -2.5);
	FEM.add_fext(coord3d( 10,  0,  0),mesh::Z, -5.0);
	FEM.add_fext(coord3d( 10, 10,  0),mesh::Z, -2.5);

	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,-10,
			mesh::X|mesh::Y|mesh::Z,0.0);

	FEM.solve(1e-25);
	int nnd = FEM.getnodes();
	for (int i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		double x = n.x;
		double y = n.y;
		double z = n.z;
		//if (y != 0.0 || x != 0.0)
		//	continue;
		double ux = n.ux;
		double uy = n.uy;
		double uz = n.uz;
		int j = FEM.hasnode(n);
		printf("Node %i: (%i,%i,%i) =\t(%f,%f,%f)\n",j,int(x),int(y),int(z),ux,uy,uz);
	}

	timeme("\n");
	
	return 0;
}

int
main_junk()
{
	double x, y, a = 0.0, b;
	double dx, dy, delta = 4.0;

	// Start with the tire grid.
	dx = delta; dy = delta;
	for (x = -30*dx; x < 30*dx; x += dx) {
		for (y = -30*dx; y < 30*dy; y += dy) {
			a += b = blockarea(x,x+dx,y,y+dx,100);
			printf("%+f %+f %f\n",x,y,b);
		}
	}
	printf("%f %f",a,M_PI*100*100);

	return 0;
}
