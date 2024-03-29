/**************************************************************************

	FEM3D.CPP - A special 3D Finite element code for my PhD.

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

	Portions Copyright (C) 2007-2022 OpenPave.org.

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
		Gaussian integration over the point stiffness matrix, with the
		appropriate transforms.  The element stiffness matrix is then
		assembled into a global stiffness matrix.  A nodal force vector is
		assembled, and the boundary conditions are applied, and then the
		global stiffness matrix is inverted to find the nodal
		displacements.  These are then used to get the final results.

		The inversion is not done directly, but via a conjugate gradient
		iterative solver, with an incomplete Cholesky decomposition as a
		pre-conditioner.

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
#include <assert.h>
#include "event.h"
#include "autodelete.h"
#include "fixed.h"
#include "set.h"
#include "tree.h"
#include "tmatrix.h"
#include "linalg.h"
#include "rng.h"
#include "pavement.h"

#if defined(__FreeBSD__)
#include <stdlib.h>
extern "C" {
const char * _malloc_options = "ajz";
}
#endif

using namespace OP;

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
static constexpr inline double
delta() noexcept
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
			const double & l, const double & m) noexcept {
		meta_stiffness<I,J,K,L-1>::subassignR(E,l,m);
		E(J-1,L-1) = l*delta<I-1,J-1>()*delta<K-1,L-1>()
				+ m*(delta<I-1,K-1>()*delta<J-1,L-1>()
				   + delta<I-1,L-1>()*delta<J-1,K-1>());
	}
	static inline void subassign(smatrix_dof & E,
			const double & l, const double & m) noexcept {
		meta_stiffness<I,J-1,K,L>::subassign(E,l,m);
		meta_stiffness<I,J,K,L>::subassignR(E,l,m);
	}
	static inline void assignR(smatrix_dof (& E)[NDIM][NDIM],
			const double & l, const double & m) noexcept {
		meta_stiffness<I,J,K-1,L>::assignR(E,l,m);
		meta_stiffness<I,J,K,L>::subassign(E[I-1][K-1],l,m);
	}
	static inline void assign(smatrix_dof (& E)[NDIM][NDIM],
			const double & l, const double & m) noexcept {
		meta_stiffness<I-1,J,K,L>::assign(E,l,m);
		meta_stiffness<I,J,K,L>::assignR(E,l,m);
	}
	static inline void submulR(ematrix & s, const ematrix & e,
			const double & l, const double & m) noexcept {
		meta_stiffness<I,J,K,L-1>::submulR(s,e,l,m);
		s(I-1,J-1) += (l*delta<I-1,J-1>()*delta<K-1,L-1>()
				+ m*(delta<I-1,K-1>()*delta<J-1,L-1>()
				   + delta<I-1,L-1>()*delta<J-1,K-1>()))
				*e(K-1,L-1);
	}
	static inline void submul(ematrix & s, const ematrix & e,
			const double & l, const double & m) noexcept {
		meta_stiffness<I,J-1,K,L>::submul(s,e,l,m);
		meta_stiffness<I,J,K,L>::submulR(s,e,l,m);
	}
	static inline void mulR(ematrix & s, const ematrix & e,
			const double & l, const double & m) noexcept {
		meta_stiffness<I,J,K-1,L>::mulR(s,e,l,m);
		meta_stiffness<I,J,K,L>::submul(s,e,l,m);
	}
	static inline void mul(ematrix & s, const ematrix & e,
			const double & l, const double & m) noexcept {
		meta_stiffness<I-1,J,K,L>::mul(s,e,l,m);
		meta_stiffness<I,J,K,L>::mulR(s,e,l,m);
	}
};
template<unsigned I, unsigned J, unsigned K>
struct meta_stiffness<I,J,K,0>
{
	static inline void subassignR(smatrix_dof &,
			const double &, const double &) noexcept {}
	static inline void submulR(ematrix &, const ematrix &,
			const double &, const double &) noexcept {}
};
template<unsigned I, unsigned K, unsigned L>
struct meta_stiffness<I,0,K,L>
{
	static inline void subassign(smatrix_dof &,
			const double &, const double &) noexcept {}
	static inline void submul(ematrix &, const ematrix &,
			const double &, const double &) noexcept {}
};
template<unsigned I, unsigned J, unsigned L>
struct meta_stiffness<I,J,0,L>
{
	static inline void assignR(smatrix_dof (&)[NDIM][NDIM],
			const double &, const double &) noexcept {}
	static inline void mulR(ematrix &, const ematrix &,
			const double &, const double &) noexcept {}
};
template<unsigned J, unsigned K, unsigned L>
struct meta_stiffness<0,J,K,L>
{
	static inline void assign(smatrix_dof (&)[NDIM][NDIM],
			const double &, const double &) noexcept {}
	static inline void mul(ematrix &, const ematrix &,
			const double &, const double &) noexcept {}
};

/*
 * Calculate the principle stresses or strains.
 */
static ematrix
principle(const ematrix & s) noexcept
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
	inline material(const double e, const double poissons) noexcept
	  : emod(e), v(poissons) {
	}
	// Get the point stiffness tensor.  This should be made to cache
	// the result.
	inline void pointstiffness(smatrix_dof (& E)[NDIM][NDIM]) const noexcept {
		const double lambda = v*emod/(1+v)/(1-2*v);
		const double mu = emod/2/(1+v);
		meta_stiffness<>::assign(E,lambda,mu);
	}
	// Get the point stress tensor, based on the strain.
	inline ematrix pointstress(const ematrix & e) const noexcept {
		const double lambda = v*emod/(1+v)/(1-2*v);
		const double mu = emod/2/(1+v);
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

	inline coord3d() noexcept {
	}
	inline coord3d(double px, double py, double pz) noexcept
	  : x(px), y(py), z(pz) {
	}
	inline coord3d(const coord3d & p) noexcept
	  : x(p.x), y(p.y), z(p.z) {
	}
	inline coord3d(const point3d & p) noexcept
	  : x(p.x), y(p.y), z(p.z) {
	}
	inline ~coord3d () {
	}
	int compare(const coord3d & p) const noexcept {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if (z > p.z) return -1;
		if (z < p.z) return  1;
		if (x < p.x) return -1;
		if (x > p.x) return  1;
		return (y < p.y ? -1 : 1);
	}
	bool operator == (const coord3d & p) const noexcept {
		return (x == p.x && y == p.y && z == p.z ? true : false);
	}
	bool operator != (const coord3d & p) const noexcept {
		return (x != p.x || y != p.y || z != p.z ? true : false);
	}
	inline bool operator > (const coord3d & p) const noexcept {
		return (compare(p) == 1 ? true : false);
	}
	inline bool operator >= (const coord3d & p) const noexcept {
		return (compare(p) != -1 ? true : false);
	}
	inline bool operator < (const coord3d & p) const noexcept {
		return (compare(p) == -1 ? true : false);
	}
	inline bool operator <= (const coord3d & p) const noexcept {
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
const double gl_4[4][2] = {{        -1.0, 1.0/6},
						   {-sqrt(5.0)/5, 5.0/6},
						   {+sqrt(5.0)/5, 5.0/6},
						   {         1.0, 1.0/6}};

class mesh;                    // Forward declare.
class node_list;               // Forward declare.
double node_depth_callback(const coord3d &, const mesh *, const material *) noexcept;
double node_emod_callback(const coord3d &, const mesh *, const material *) noexcept;
double node_depth_callback_test(const coord3d &, const mesh *, const material *) noexcept;
double node_emod_callback_test(const coord3d &, const mesh *, const material *) noexcept;

/*
 * struct node3d
 *
 * A node in 3D, based on coord3d.  This stores the meshes neighbours,
 * so it can keep track of points for the variable node elements.
 */
struct node3d : public coord3d {
	// Neighbour in (x,y,z) plus/minus directions.
	// UINT_MAX if not there...
	//struct {
		unsigned xm, xp, ym, yp, zm, zp;
	//};
	// Final results.
	//struct {
		double ux, uy, uz;
	//};
	// Set new neighbours.  If this conflicts with what we know, we
	// complain (if we're debugging).
	void setneighbours(const unsigned x_m, const unsigned x_p,
			const unsigned y_m, const unsigned y_p,
			const unsigned z_m, const unsigned z_p) noexcept {
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
	bool red:1;                // Colour of link to parent
	unsigned fixed:31;         // Bitmask of fixed DOFs

	explicit node3d(const coord3d & c) noexcept
	  : coord3d(c), xm(UINT_MAX), xp(UINT_MAX), ym(UINT_MAX), yp(UINT_MAX),
			zm(UINT_MAX), zp(UINT_MAX), ux(0.0), uy(0.0), uz(0.0),
			order(0), left(UINT_MAX), right(UINT_MAX), red(true), fixed(0) {
	}
	// Set our results.
	void setdisp(const tmatrix<double,NDIM,1> & t) noexcept {
		ux = t(0); uy = t(1); uz = t(2);
	}
	// Placement new to support inplace init in the list.
	void * operator new (size_t, void * p) noexcept {
		return p;
	}
	void operator delete (void *, void *) noexcept {
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
	  : size(0), buffer(0), block(DFLT_BLK), root(UINT_MAX), value(nullptr) {
	}
	// Clean up.
	inline ~node_list() {
		if (value)
			free(value);
	}
	// The length. Nice for lots of things...
	inline unsigned length() const noexcept {
		return size;
	}
	// Do a key lookup, and return UINT_MAX if the key is not found.
	inline unsigned haskey(const coord3d & k) const noexcept {
		unsigned x = root;
		while (x != UINT_MAX) {
			const int cmp = k.compare(value[x]);
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
	inline node3d & operator [] (const unsigned p) const noexcept {
		return value[p];
	}
	// Allow sorted access.
	inline node3d & getindex(const unsigned i) const noexcept {
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
	inline unsigned getorder(const unsigned i) const noexcept {
		const coord3d & k = static_cast<coord3d &>(value[i]);
		unsigned x = root, order = 0;
		while (x != UINT_MAX) {
			const int cmp = k.compare(value[x]);
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
	inline void add(const coord3d & v) {
		allocate(size+1);
		append(root,v);
		value[root].red = false;
	}

protected:
	unsigned size;             // The size of the set...
	unsigned buffer;           // The allocated buffer size...
	unsigned block;            // The minimum block size.
	unsigned root;             // Root of red-black tree.
	node3d * value;            // Take a guess...

	// Make some space...
	void allocate(const unsigned s) {
		while (s > 8*block)
			block *= 8;
		const unsigned b = block*(s/block+((s%block)?1:0));
		if (b == buffer)
			return;
		node3d * temp = static_cast<node3d *>(realloc(value,
				b*sizeof(node3d)));
		if (temp == nullptr)
			throw std::bad_alloc();
		value = temp;
		buffer = b;
	}
	void append(unsigned & r, const coord3d & c) {
		// If we're UINT_MAX that means we need to make a new node...
		if (r == UINT_MAX) {
			new(&value[size]) node3d(c);
			r = size++;
			return;
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
		const int cmp = c.compare(value[r]);
		if (cmp == 0) {
			assert(false);
			return;
		} else if (cmp < 0) {
			// Re-root insert into left tree.
			append(value[r].left,c);
			value[r].order++;
		} else {
			// Or right tree.  Same issue here as above.
			append(value[r].right,c);
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
	  : K(nullptr), nnd(n) {
		// calloc() to zero the memory...
		K = static_cast<smatrix_dof *>
				(calloc(size(),sizeof(smatrix_dof)));
		if (K == nullptr)
			throw std::bad_alloc();
	}
	inline ~smatrix_elem() {
		if (K != nullptr)
			free(K);
	}
	inline smatrix_dof & operator () (unsigned i, unsigned j) const {
		assert(j >= i);
		return K[index(i,j)];
	}

private:
	friend class mesh;

	smatrix_dof * K;
	const unsigned nnd;

	// Size of the triangular matrix storage
	inline unsigned size() const noexcept {
		return nnd*(nnd+1)/2;
	}
	// Index into the triangular matrix storage
	inline unsigned index(const unsigned i, const unsigned j) const noexcept {
		return j*(j+1)/2 + i;
	}
};

/*
 * This is a specialized version of the general inv_mul_gauss function,
 * which does not pivot, and has fixed size, since we know that we never
 * pivot for the Hessian, and it is always NDIMxNDIM.  It also takes the
 * transpose of the B matrix, since the column dimension is fixed.
 */
static double
inv_mul_gauss(const unsigned nnd, double (& J)[NDIM][NDIM], double (* B)[NDIM])
{
	double det = 1.0;
	unsigned i, j, k;

	for (i = 0; i < NDIM; i++) {
		const double pvt = J[i][i];
		if (fabs(pvt) < DBL_EPSILON)
			throw std::range_error("Singular matrix in inv_mul_gauss()!");
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

	virtual unsigned l2g(const unsigned i) const noexcept = 0;
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
template<element::shape_t SZ, unsigned NND>
class element_base : public element {
public:
	element_base(mesh * o, const shape_t x, const shape_t y,const material & m)
	  : element(o,m), sx(x), sy(y) {
	}
	unsigned l2g(const unsigned i) const noexcept override {
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
		assert(mask != nullptr);
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
		for (unsigned i = 0; sx == inf_pos && i < nnd/4; i++)
			xe[i*4+2][2] = xe[i*4][2], xe[i*4+3][2] = xe[i*4+1][2];
		for (unsigned i = 0; sx == inf_neg && i < nnd/4; i++)
			xe[i*4][2] = xe[i*4+2][2], xe[i*4+1][2] = xe[i*4+3][2];
		for (unsigned i = 0; sy == inf_pos && i < nnd/4; i++)
			xe[i*4+1][2] = xe[i*4][2], xe[i*4+3][2] = xe[i*4+2][2];
		for (unsigned i = 0; sy == inf_neg && i < nnd/4; i++)
			xe[i*4][2] = xe[i*4+1][2], xe[i*4+2][2] = xe[i*4+3][2];
	}
	// Build a matrix of the emod delta values.
	void buildEe(double * Ee) const noexcept {
		for (unsigned i = 0; i < nnd; i++) {
			//const coord3d & c = getnode(inel[i]);
			//Ee[i] = node_emod_callback(c,owner,&mat);
			Ee[i] = mat.emod;
		}
		for (unsigned i = 0; sx == inf_pos && i < nnd/4; i++)
			Ee[i*4+2] = Ee[i*4], Ee[i*4+3] = Ee[i*4+1];
		for (unsigned i = 0; sx == inf_neg && i < nnd/4; i++)
			Ee[i*4] = Ee[i*4+2], Ee[i*4+1] = Ee[i*4+3];
		for (unsigned i = 0; sy == inf_pos && i < nnd/4; i++)
			Ee[i*4+1] = Ee[i*4], Ee[i*4+3] = Ee[i*4+2];
		for (unsigned i = 0; sy == inf_neg && i < nnd/4; i++)
			Ee[i*4] = Ee[i*4+1], Ee[i*4+2] = Ee[i*4+3];
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
			double * N, double (* dNdr)[NDIM]) const noexcept = 0;
	// build the shape functions for linear and infinite elements
	// in plan view.
	void buildSF_block(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const noexcept {
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
			const unsigned * mask) const noexcept {
		assert(sx == absolute);
		assert(sy == absolute);
		constexpr const unsigned nz = unsigned(SZ);
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
		// depending on the existence of the extra nodes.
		for (unsigned l = 0; l < 4*nz; l++) {
			if ((j = mask[l]) != 0) {
				N[l] -= N[j]/2;
				dNdr[l][0] -= dNdr[j][0]/2;
				dNdr[l][1] -= dNdr[j][1]/2;
				dNdr[l][2] -= dNdr[j][2]/2;
			}
			const int Jxy[4] = {+2,-1,+1,-2};
			if ((j = mask[int(l)+Jxy[l%4]]) != 0) {
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
			const unsigned * mask) const noexcept {
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
				if (mask[int(i*4+l)+Jx[l]] != 0) {
					Nx = N_3(2*(l/2),gx); dNdx = dN_3(2*(l/2),gx);
				} else {
					Nx = N_2(l/2,gx);     dNdx = dN_2(l/2,gx);
				}
				if (mask[int(i*4+l)+Jy[l]] != 0) {
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
		const double (* gx)[2] = nullptr, (* gy)[2] = nullptr,
				(* gz)[2] = nullptr;
		getgauss(sx,nx,gx);
		getgauss(sy,ny,gy);
		getgauss(SZ,nz,gz);

		smatrix_elem * K = new smatrix_elem(nnd);
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
				emod += N[k]*(log10(Ee[k])-log10(mat.emod));
			// This returns det(J);
			gw *= inv_mul_gauss(nnd,J,dNdr)*pow(10,emod);
			if (gw <= 0.0) {
				printf("\nFailed at gauss point: %d %d %d %f %f %f %f %f\n",
					gi,gj,gk,gx[gi][0],gy[gj][0],gz[gk][0],gx[gi][1]*gy[gj][1]*gz[gk][1],emod);
				printf("J residual %f %f:\n",gw,gw/gx[gi][1]*gy[gj][1]*gz[gk][1]);
				for (i = 0; i < NDIM; i++) {
					for (j = 0; j < NDIM; j++)
						printf(" %10.4g",J[i][j]);
					printf("\n");
				}
				printf("xe:\n");
				for (i = 0; i < nnd; i++) {
					for (j = 0; j < NDIM; j++)
						printf(" %10.4g",xe[i][j]);
					printf("\n");
				}
				printf("residual dNdr:\n");
				for (i = 0; i < nnd; i++) {
					for (j = 0; j < NDIM; j++)
						printf(" %10.4g",dNdr[i][j]);
					printf("\n");
				}
				buildSF(true,gx[gi][0],gy[gj][0],gz[gk][0],N,dNdr);
				printf("N:\n");
				for (i = 0; i < nnd; i++) {
					printf("%10.4g\n",N[i]);
				}
				printf("dNdr:\n");
				for (i = 0; i < nnd; i++) {
					for (j = 0; j < NDIM; j++)
						printf(" %10.4g",dNdr[i][j]);
					printf("\n");
				}
				for (i = 0; i < NDIM; i++) {
					for (j = 0; j < NDIM; j++) {
						J[i][j] = 0.0;
						for (k = 0; k < nnd; k++)
							J[i][j] += dNdr[k][i]*xe[k][j];
					}
				}
				printf("J:\n");
				for (i = 0; i < NDIM; i++) {
					for (j = 0; j < NDIM; j++)
						printf(" %10.4g",J[i][j]);
					printf("\n");
				}
				delete K;
				throw std::runtime_error("Bad determinate in element_base::buildKe()!");
			}
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
			d.deflgrad.resize(1);
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
				emod += N[k]*(log10(Ee[k])-log10(mat.emod));
			d.deflgrad[0] = pow(10,log10(mat.emod)+emod);
			ematrix e(0.0);
			for (i = 0; i < NDIM; i++) {
				for (j = 0; j < NDOF; j++) {
					for (k = 0; k < nnd; k++)
						e(i,j) += dNdr[k][i]*ue[k][j];
				}
			}
			e = (e + ~e)*0.5;
			ematrix s(mat.pointstress(e)*pow(10,emod));
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
			throw std::runtime_error("Bad or degenerate element shape!");
		} while(true);
	}
	// This adds nodes according to linear shape functions...
	inline void coord_add(shape_t ssx, shape_t ssy,
			const coord3d (& cc)[8]) {
		assert(nnd == 0);
		const double nx = getn(ssx),    mx = nx-1;
		const double ny = getn(ssy),    my = ny-1;
		const double nz = unsigned(SZ), mz = nz-1;

		// Figure out our nodes using linear shape functions...
		for (double i = 0; i < nz; i++) {
			for (double j = 0; j < nx; j++) {
				for (double k = 0; k < ny; k++) {
					double x = 0.0, y = 0.0, z = 0.0;
					for (unsigned l = 0; l < 8; l++) {
						const double s = (l   < 4 ? mz-i : i)
								*(l%4 < 2 ? mx-j : j)
								*(l%2 < 1 ? my-k : k)/(mz*mx*my);
						x += s*double(cc[l].x);
						y += s*double(cc[l].y);
						z += s*double(cc[l].z);
					}
					inel[nnd++] = addnode(coord3d(x,y,z));
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
			const double (*& g)[2]) noexcept {
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
	static inline double N_2(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return 0.5*(1.0-r);
		case 1: return 0.5*(1.0+r);
		} assert(false); return 0.0;
	}
	static inline double dN_2(const unsigned n, const double ) noexcept {
		switch (n) {
		case 0: return -0.5;
		case 1: return  0.5;
		} assert(false); return 0.0;
	}
	// Absolute value shape function.
	static inline double N_A(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return (r < 0.0 ? -r : 0.0);
		case 1: return 1-fabs(r);
		case 2: return (r > 0.0 ?  r : 0.0);
		} assert(false); return 0.0;
	}
	static inline double dN_A(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return (r < 0.0 ? -1.0 : r > 0.0 ?  0.0 : -0.5);
		case 1: return (r < 0.0 ?  1.0 : r > 0.0 ? -1.0 :  0.0);
		case 2: return (r < 0.0 ?  0.0 : r > 0.0 ?  1.0 :  0.5);
		} assert(false); return 0.0;
	}
	// Quadratic shape function.
	static inline double N_3(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return 0.5*r*(r-1);
		case 1: return 1-r*r;
		case 2: return 0.5*r*(r+1);
		} assert(false); return 0.0;
	}
	static inline double dN_3(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return r-0.5;
		case 1: return -2*r;
		case 2: return r+0.5;
		} assert(false); return 0.0;
	}
	// Cubic shape function.
	static inline double N_4(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return 0.0625*(-1+9*r*r)*(1-1*r);
		case 1: return 0.0625*(+9-9*r*r)*(1-3*r);
		case 2: return 0.0625*(+9-9*r*r)*(1+3*r);
		case 3: return 0.0625*(-1+9*r*r)*(1+1*r);
		} assert(false); return 0.0;
	}
	static inline double dN_4(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return 0.0625*((+18*r)*(1-1*r) - 1*(-1+9*r*r));
		case 1: return 0.0625*((-18*r)*(1-3*r) - 3*(+9-9*r*r));
		case 2: return 0.0625*((-18*r)*(1+3*r) + 3*(+9-9*r*r));
		case 3: return 0.0625*((+18*r)*(1+1*r) + 1*(-1+9*r*r));
		} assert(false); return 0.0;
	}
	// mapping functions for infinite elements.
	static inline double Mpi(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return  -2*r/(1-r);
		case 1: return (1+r)/(1-r);
		} assert(false); return 0.0;
	}
	static inline double dMpi(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return -2/((1-r)*(1-r));
		case 1: return  2/((1-r)*(1-r));
		} assert(false); return 0.0;
	}
	static inline double Mni(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return (1-r)/(1+r);
		case 1: return   2*r/(1+r);
		} assert(false); return 0.0;
	}
	static inline double dMni(const unsigned n, const double r) noexcept {
		switch (n) {
		case 0: return -2/((1+r)*(1+r));
		case 1: return  2/((1+r)*(1+r));
		} assert(false); return 0.0;
	}
	// Shape functions for infinite elements are quadratic.
	static inline double Npi(const unsigned n, const double r) noexcept {
		return N_3(n,r);
	}
	static inline double dNpi(const unsigned n, const double r) noexcept {
		return dN_3(n,r);
	}
	// Shift the negative nodes to (0,+1).
	static inline double Nni(const unsigned n, const double r) noexcept {
		return N_3(n+1,r);
	}
	static inline double dNni(const unsigned n, const double r) noexcept {
		return dN_3(n+1,r);
	}
	// Return the number of nodes in each shape function.
	static inline unsigned getn(const shape_t s, const bool ismap = true) noexcept {
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
			const double r, const bool ismap = true) noexcept {
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
			const double r, const bool ismap = true) noexcept {
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
	void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const noexcept override {
		assert(4*int(SZ) == this->nnd);
		base_t::buildSF_block(ismap,gx,gy,gz,N,dNdr);
	}
	smatrix_elem * stiffness() const override {
		assert(4*int(SZ) == this->nnd);
		return base_t::buildKe();
	}
	void results(const fset<point3d> & c, fset<pavedata> & d) const override {
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
	void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const noexcept override {
		assert(9*int(SZ) == this->nnd);
		base_t::buildSF_block(ismap,gx,gy,gz,N,dNdr);
	}
	smatrix_elem * stiffness() const override {
		assert(9*int(SZ) == this->nnd);
		return base_t::buildKe();
	}
	void results(const fset<point3d> & c, fset<pavedata> & d) const override {
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
	void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const noexcept override {
		assert(ismap == true);
		base_t::buildSF_variable(gx,gy,gz,N,dNdr,mask);
	}
	smatrix_elem * stiffness() const override {
		return base_t::buildKe();
	}
	void results(const fset<point3d> & c, fset<pavedata> & d) const override {
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
	void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const noexcept override {
		assert(ismap == true);
		base_t::buildSF_adaptor(gx,gy,gz,N,dNdr,mask);
	}
	smatrix_elem * stiffness() const override {
		return base_t::buildKe();
	}
	void results(const fset<point3d> & c, fset<pavedata> & d) const override {
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
	void buildSF(const bool ismap, const double & gx,
			const double & gy, const double & gz,
			double * N, double (* dNdr)[NDIM]) const noexcept override {
		assert(4*int(SZ) == this->nnd);
		base_t::buildSF_block(ismap,gx,gy,gz,N,dNdr);
	}
	smatrix_elem * stiffness() const override {
		assert(4*int(SZ) == this->nnd);
		return base_t::buildKe();
	}
	void results(const fset<point3d> & c, fset<pavedata> & d) const override {
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
	void expand(unsigned nn) {
		if (nnz+nn > nnd) {
			// We're out of space, so allocate more.
			nnd = nnz+nn;
			smatrix_node * temp =
					static_cast<smatrix_node *>
					(realloc(nodes,nnd*sizeof(smatrix_node)));
			if (temp == nullptr)
				throw std::bad_alloc();
			// If realloc moved us we need to fix up everyone's pointers.
			// This finds the offset into the old array, and then adds
			// that to the new array.
			if (nodes != temp) {
				if (col_head != nullptr)
					col_head = temp + (col_head - nodes);
				for (unsigned k = 0; k < nnz; k++) {
					smatrix_node * n = &(temp[k]);
					if (n->col_next != nullptr)
						n->col_next = temp + (n->col_next - nodes);
				}
				nodes = temp;
			}
		}
	}

	friend class smatrix;
	friend class mesh;

	smatrix_dof K;
	smatrix_node * col_head;
	smatrix_node * nodes;
	unsigned nnz, nnd;
};

// Cannot use local types in templates...
struct row_ptr {
	int col;
	row_ptr *prev;
};

/*
 * class smatrix - a symmetric positive definite stiffness matrix
 */
class smatrix {
public:
	// Create the matrix empty.
	inline explicit smatrix(const unsigned n, double d = 0.0)
	  : nnd(n), diag(nullptr) {
		diag = static_cast<smatrix_diag *>
				(calloc(nnd,sizeof(smatrix_diag)));
		if (diag == nullptr)
			throw std::bad_alloc();
		for (unsigned i = 0; d != 0.0 && i < nnd; i++) {
			for (unsigned j = 0; j < NDOF; j++)
				diag[i].K(j,j) = d;
		}
	}
	// Copy a matrix.
	inline explicit smatrix(const smatrix & A)
	  : nnd(A.nnd), diag(nullptr) {
		diag = static_cast<smatrix_diag *>
				(calloc(nnd,sizeof(smatrix_diag)));
		if (diag == nullptr)
			throw std::bad_alloc();
		for (unsigned i = 0; i < nnd; i++) {
			smatrix_diag * d = &(A.diag[i]);
			// Pre-allocate the space.
			diag[i].expand(d->nnz);
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
	void append(unsigned i, unsigned j, const smatrix_dof & t) {
		// If we're appending to the diagonal, just do it.
		if (i == j) {
			diag[i].K += t;
			return;
		}
		// Here we don't know if we're above or below the diagonal,
		// so we must use MIN/MAX.  We start by trying to find the
		// node which is in a row greater than or equal to the one
		// we're adding.
		smatrix_diag * d = &diag[MAX(i,j)];
		smatrix_node * n, * o, * p = d->col_head;
		while (p != nullptr && p->i < MIN(i,j))
			p = p->col_next;
		// If we found a match, add to it.
		if (p != nullptr && p->i == MIN(i,j)) {
			if (j < i)
				p->K += ~t;
			else
				p->K += t;
			return;
		}
		// We've either run off the end, or we don't have a matching
		// node.  In this case we need to add one.
		d->expand(1);
		// Now fill in the new element.
		n = &(d->nodes[(d->nnz)++]);
		if (j < i) {
			n->K = ~t;
			swap(i,j);
		} else
			n->K = t;
		n->i = i;
		// Fix up the column references.
		p = d->col_head, o = nullptr;
		while (p != nullptr && p->i < i)
			o = p, p = p->col_next;
		assert(p == nullptr || p->i > i);
		if (o == nullptr)
			d->col_head = n;
		else
			o->col_next = n;
		n->col_next = p;
	}
	// This just makes things faster by sorting the node lists.
	void tidy() noexcept {
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
			if (diag[i].col_head != nullptr)
				diag[i].col_head = nodes;
			for (unsigned short j = 0; j < diag[i].nnz; j++)
				nodes[j].col_next = (j == diag[i].nnz-1 ? nullptr : &nodes[j+1]);
		}
	}
	// Do a level 1 fill in on the sparse matrix for incomplete Cholesky
	// decomposition.  We need to insert nodes, which is costly, so we rely
	// on some tricks here...  We don't use nnz, only col_next, while tidy
	// will only use nnz not col_next.  So we add nodes ourself without
	// using append, and memset() them to zero, and set their row numbers.
	// Then we call tidy(), which will sort by row numbers and fix the
	// col_next pointers.
	void fillin() {
		smatrix_diag * d;
		smatrix_node * p;

		// The fill in works best if we can follow rows and columns.  Since
		// we don't normally store rows, we need to construct the row
		// linked lists here.
		unsigned nnz = 0;
		for (unsigned n = 0; n < nnd; n++)
			nnz += diag[n].nnz;
		autodelete<row_ptr> rows(new row_ptr[nnd+nnz]);
		memset(rows,0,(nnd+nnz)*sizeof(row_ptr));
		nnz = 0;
		for (unsigned n = 0; n < nnd; n++) {
			p = diag[n].col_head;
			while (p) {
				rows[nnd+nnz].col = int(n);
				rows[nnd+nnz].prev = rows[p->i].prev;
				rows[p->i].prev = &rows[nnd+nnz];
				p = p->col_next;
				nnz++;
			}
		}

		// Now that the list is constructed we do fill in from the last
		// element.
		for (unsigned n = nnd; n > 0; ) {
			d = &(diag[--n]);
			p = d->col_head;
			if (p == nullptr)
				continue; // Skip diag only
			unsigned add = 0;
			const unsigned c = p->i; // The max row for this column
			while (p) {
				// We use the first n-c columns from the diagonal to store
				// a mask for which elements to add.
				rows[n-p->i-1].col = 1;
				row_ptr * r = rows[p->i].prev->prev;
				rows[p->i].prev = r;
				while (r) {
					rows[int(n)-r->col-1].col = -1;
					r = r->prev;
				}
				p = p->col_next;
			}
			for (unsigned i = 0; i < n-c; i++) {
				if (rows[i].col == -1)
					add++;
			}
			d->expand(add);
			memset(&(d->nodes[d->nnz]),0,add*sizeof(smatrix_node));
			for (unsigned i = 0, j = add; i < n-c; i++) {
				if (rows[i].col == -1)
					d->nodes[--j+d->nnz].i = n-1-i;
				rows[i].col = 0; // cleanup
			}
			d->nnz += add;
		}
		// Now call tidy to magically add these nodes.
		tidy();
	}
	// Incomplete Cholesky decomposition, by block.  Look at the small
	// internal version to understand what the whole thing is doing.
	void incchol() {
		smatrix_diag * d;
		smatrix_node * p;

		// Do the decomposition.
		for (unsigned n = 0; n < nnd; n++) {
			d = &(diag[n]);
			p = d->col_head;
			while (p != nullptr) {
				smatrix_node * pi = diag[p->i].col_head;
				smatrix_node * pj = d->col_head;
				// We need to have blocks in both columns.  So skip the
				// missing ones.
				while (pi != nullptr && pj != nullptr) {
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
			// Cholesky decomposition of the diagonal.
			for (unsigned i = 0; i < NDOF; i++) {
				for (unsigned k = 0; k < i; k++)
					d->K(i,i) -= d->K(k,i)*d->K(k,i);
				if (d->K(i,i) < DBL_EPSILON) {
					printf("OOOOPS... (%g) %i/%i %i\n",d->K(i,i),n,nnd,i);
					throw std::range_error("Non-positive definite matrix in incchol()!");
				}
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
	// Do a complete Cholesky decomposition, by filling in all the blocks and
	// then doing an 'incomplete' Cholesky...
	void chol() {
		smatrix_diag * d;
		smatrix_node * p;

		for (unsigned n = 0; n < nnd; n++) {
			d = &(diag[n]);
			if (!(p = d->col_head))
				continue; // Diagonal only is OK.
			unsigned r = p->i;
			const unsigned a = r - n - d->nnz;
			d->expand(a);
			memset(&(d->nodes[d->nnz]),0,a*sizeof(smatrix_node));
			p = d->col_head; // The array's probably been realloced.
			for (; r < n; r++) {
				while (p && p->i < r)
					p = p->col_next;
				if (p && p->i == r)
					continue;
				d->nodes[(d->nnz)++].i = r;
			}
		}
		// Now call tidy to magically add these nodes.
		tidy();
		incchol();
	}
	unsigned nonzero() noexcept {
		smatrix_diag * d;
		smatrix_node * p;
		unsigned nnz = nnd;

		for (unsigned n = 0; n < nnd; n++) {
			d = &(diag[n]);
			if (!(p = d->col_head))
				continue; // Diagonal only is OK.
			//unsigned r = p->i;
			nnz += p->i - n;
		}
		return nnz;
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
	  : nnd(n), V(nullptr) {
		V = static_cast<svector_dof *>
				(calloc(nnd,sizeof(svector_dof)));
		if (V == nullptr)
			throw std::bad_alloc();
	}
	inline ~svector() {
		if (V != nullptr)
			free(V);
	}
	inline const svector_dof & operator () (const unsigned i) const noexcept {
		return V[i];
	}
	inline svector_dof & operator () (const unsigned i) noexcept {
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
	mesh_bc_key() noexcept
	  : n(0), i(0) {
	}
	mesh_bc_key(unsigned N, unsigned I) noexcept
	  : n(N), i(I) {
	}
	// Provide these for sorting.
	bool operator == (const mesh_bc_key & k) const noexcept {
		return (n == k.n && i == k.i);
	}
	bool operator > (const mesh_bc_key & k) const noexcept {
		return (n > k.n || (n == k.n && i > k.i));
	}
	bool operator < (const mesh_bc_key & k) const noexcept {
		return (n < k.n || (n == k.n && i < k.i));
	}
};
struct mesh_bc : public mesh_bc_key {
	double d;
	mesh_bc() noexcept
	  : mesh_bc_key(), d(0.0) {
	}
	mesh_bc(unsigned N, unsigned I, double D) noexcept
	  : mesh_bc_key(N,I), d(D) {
	}
};

/*
 * Nested dissection node renumbering, for when I figure out a multi-frontal
 * solver...
 */

/*
 * struct mesh_part - a field for breadth first walking of a node list.
 */
struct mesh_part : public fixed<8> {
	unsigned order;

	mesh_part(fixed<8> & p) noexcept
	  : fixed<8>(p), order(0) {
	}
	static void setorder(koset<fixed<8>,mesh_part> & p) {
		p.sort();
		buildtree(p,0,p.length()-1,1,0,0);
	}

private:
	static void buildtree(koset<fixed<8>,mesh_part> & p, const unsigned l,
			const unsigned r, unsigned v, unsigned o, unsigned f) {
		const unsigned i = MIN(r,l+(r-l+f%2)/2);
		p[i].order = o+(f==0?0:v/3*f);
		if (l < i && l < r)
			buildtree(p,l,i-1,v*3,p[i].order,1);
		if (r > i && l < r)
			buildtree(p,i+1,r,v*3,p[i].order,2);
	}
};
/*
 * struct mesh_part3d
 */
struct mesh_part3d {
	unsigned ox, oy, oz;

	inline mesh_part3d(mesh_part x, mesh_part y, mesh_part z) noexcept
	  : ox(x.order), oy(y.order), oz(z.order) {
	}
	inline mesh_part3d(const mesh_part3d & p) noexcept
	  : ox(p.ox), oy(p.oy), oz(p.oz) {
	}
	inline ~mesh_part3d () {
	}
	int compare(const mesh_part3d & p) const noexcept {
		unsigned xo1, xo2, yo1, yo2, zo1, zo2;
		unsigned m = 0, o1 = 0, o2 = 0;

		xo1 = ox; xo2 = p.ox;
		yo1 = oy; yo2 = p.oy;
		zo1 = oz; zo2 = p.oz;
		while (xo1 > 0 || yo1 > 0 || zo1 > 0
		    || xo2 > 0 || yo2 > 0 || zo2 > 0) {
			switch ((m++)%3) {
			case 0:
				o1 = xo1%3; xo1 /= 3;
				o2 = xo2%3; xo2 /= 3;
				break;
			case 1:
				o1 = yo1%3; yo1 /= 3;
				o2 = yo2%3; yo2 /= 3;
				break;
			case 2:
				o1 = zo1%3; zo1 /= 3;
				o2 = zo2%3; zo2 /= 3;
				break;
			}
			if (o1 != o2) {
				return (o1 < o2 ? 1 : -1); // reverse order
			}
		}
		return 0;
	}
	bool operator == (const mesh_part3d & p) const noexcept {
		return (compare(p) == 0 ? true : false);
	}
	bool operator != (const mesh_part3d & p) const noexcept {
		return (compare(p) != 0 ? true : false);
	}
	inline bool operator > (const mesh_part3d & p) const noexcept {
		return (compare(p) == 1 ? true : false);
	}
	inline bool operator >= (const mesh_part3d & p) const noexcept {
		return (compare(p) != -1 ? true : false);
	}
	inline bool operator < (const mesh_part3d & p) const noexcept {
		return (compare(p) == -1 ? true : false);
	}
	inline bool operator <= (const mesh_part3d & p) const noexcept {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * class mesh_rcm - support for Reverse Cuthill McGee ordering
 */
class mesh_rcm {
public:
	explicit inline mesh_rcm(unsigned n)
	  : nnd(n), rtail(UINT_MAX), qhead(UINT_MAX),
	    qtail(UINT_MAX), phead(0), diag(nullptr) {
		diag = static_cast<mesh_rcm_node *>
				(malloc(nnd*sizeof(mesh_rcm_node)));
		if (diag == nullptr)
			throw std::bad_alloc();
		for (unsigned i = 0; i < nnd; i++)
			new(&diag[i]) mesh_rcm_node(i-1,i+1);
		diag[0].prev = diag[nnd-1].next = UINT_MAX;
	}
	inline void append(unsigned i, unsigned j) {
		if (i == j)
			return;
		diag[i].adj.add(j);
	}
	inline void compute(const fset<unsigned> & n2o) {
		for (unsigned i = 0; i < nnd; i++)
			diag[i].adj.sort();
		for (unsigned i = phead; i != UINT_MAX; i = diag[i].next) {
			for (unsigned j = i, k = diag[j].prev; j != UINT_MAX && k != UINT_MAX
					&& diag[k].adj.length() > diag[j].adj.length(); j = k, k = diag[j].prev) {
				phead = (phead == j ? k : phead == k ? j : phead);
				if (diag[k].prev != UINT_MAX)
					diag[diag[k].prev].next = j;
				if (diag[j].next != UINT_MAX)
					diag[diag[j].next].prev = k;
				diag[k].next = diag[j].next;
				diag[j].next = k;
				diag[j].prev = diag[k].prev;
				diag[k].prev = j;
				k = j;
			}
		}
		while (phead != UINT_MAX) {
			move2r(phead);
			while (qhead != UINT_MAX)
				move2r(qhead);
		}
		for (unsigned i = 0; rtail != UINT_MAX;	i++, rtail = diag[rtail].prev)
			diag[rtail].next = i;
		for (unsigned i = 0; i < nnd; i++)
			n2o[i] = diag[n2o[i]].next;
	}

private:
	// move a node from q or p to r.
	inline void move2r(const unsigned a) {
		// remove a from the list
		if (diag[a].prev != UINT_MAX)
			diag[diag[a].prev].next = diag[a].next;
		if (diag[a].next != UINT_MAX)
			diag[diag[a].next].prev = diag[a].prev;
		// remove from q (normally qhead = a)
		if (qhead == a && qtail == a)
			qhead = qtail = UINT_MAX;
		else if (qhead == a)
			qhead = diag[a].next;
		// remove from p (only happens if q was empty)
		if (phead == a)
			phead = diag[a].next;
		// place it in r
		diag[a].prev = rtail;
		if (rtail == UINT_MAX)  // first element
			diag[a].next = phead;
		else
			diag[a].next = diag[rtail].next;
		if (diag[a].next != UINT_MAX)
			diag[diag[a].next].prev = a;
		rtail = a;
		// now move all the adjacent nodes to q.
		for (unsigned p = phead, pn; p != UINT_MAX; p = pn) {
			pn = diag[p].next;
			if (diag[a].adj.findvalue(p) != UINT_MAX)
				move2q(p);
		}
	}
	// move a node from p to q.
	inline void move2q(const unsigned a) {
		// remove a from the list
		if (diag[a].prev != UINT_MAX)
			diag[diag[a].prev].next = diag[a].next;
		if (diag[a].next != UINT_MAX)
			diag[diag[a].next].prev = diag[a].prev;
		// place it in q
		if (qhead == UINT_MAX)
			qhead = a;
		diag[a].prev = (qtail == UINT_MAX ? rtail : qtail);
		qtail = a;
		// remove from p
		if (phead == a)
			phead = diag[a].next;
		else
			diag[a].next = phead;
		// fix the list
		if (diag[a].prev != UINT_MAX)
			diag[diag[a].prev].next = a;
		if (diag[a].next != UINT_MAX)
			diag[diag[a].next].prev = a;
	}

	unsigned nnd;             // The number of nodes
	unsigned rtail;           // The tail of the final list
	unsigned qhead;           // The head of the working queue
	unsigned qtail;           // The tail of the working queue
	unsigned phead;           // The head of the unprocessed queue
	struct mesh_rcm_node {
		unsigned prev;
		unsigned next;
		cset<unsigned> adj;

		inline mesh_rcm_node(const unsigned p, const unsigned n)
		  : prev(p), next(n), adj() {
		}
		// Placement new to support in-place constructor in the list.
		void * operator new (size_t, void * p) {
			return p;
		}
		void operator delete (void *, void *) {
		}
	} * diag;
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
		element * e = nullptr;
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
			if (e->nnd == 8) {
				delete e;
				e = new element_block<element::linear>(this,m,c);
			}
			break;
		case element::variable26:
			e = new element_variable<element::quadratic>(this,m,c);
			if (e->nnd == 12) {
				delete e;
				e = new element_block<element::quadratic>(this,m,c);
			}
			break;
		case element::variable34:
			e = new element_variable<element::cubic>(this,m,c);
			if (e->nnd == 16) {
				delete e;
				e = new element_block<element::cubic>(this,m,c);
			}
			break;
		case element::adaptor18:
			e = new element_adaptor<element::linear>(this,m,c);
			if (e->nnd == 8) {
				delete e;
				e = new element_block<element::linear>(this,m,c);
			}
			break;
		case element::adaptor26:
			e = new element_adaptor<element::quadratic>(this,m,c);
			if (e->nnd == 12) {
				delete e;
				e = new element_block<element::quadratic>(this,m,c);
			}
			break;
		case element::adaptor34:
			e = new element_adaptor<element::cubic>(this,m,c);
			if (e->nnd == 16) {
				delete e;
				e = new element_block<element::cubic>(this,m,c);
			}
			break;
		}
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
		const unsigned k = node.haskey(n);
		if (k == UINT_MAX)
			return false;
		node[k] = n;
		return true;
	}
	inline unsigned getnodes() const noexcept {
		return node.length();
	}
	inline const node3d & getnode(const unsigned i) const {
		return node[i];
	}
	inline const node3d & getorderednode(const unsigned i) const {
		return node.getindex(i);
	}
	inline unsigned hasnode(const coord3d & p) const noexcept {
		return node.haskey(p);
	}
	inline unsigned getorderofnode(const coord3d & p) const noexcept {
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
		const unsigned k = node.haskey(p);
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
	// Add constraints along a plane.
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
			}
		}
		return true;
	}
	// Add an external force at a node.
	bool add_fext(const coord3d & p, const dof f, const double d) {
		unsigned j;
		const unsigned k = node.haskey(p);
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
		unsigned i, j, n;
		const unsigned nnd = node.length();
		printf("Solving with %u nodes!",nnd);
		smatrix K(nnd);
		svector F(nnd), U(nnd), P(nnd), W(nnd), V(nnd), R(nnd);
		smatrix_diag * d;
		smatrix_node * p;
		const element * e;
		//koset<fixed<8>,mesh_part> partition[3];
		//ktree_llrb<mesh_part3d,mesh_part3d> partition3d;
		fset<unsigned> n2o(nnd);
		//mesh_rcm rcm_order(nnd);

		timeme("\nOrdering nodes...");
		//for (i = 0; i < nnd; i++) {
		//	node3d & u = node[i];
		//	partition[0].add(mesh_part(u.x));
		//	partition[1].add(mesh_part(u.y));
		//	partition[2].add(mesh_part(u.z));
		//}
		//for (i = 0; i < 3; i++)
		//	mesh_part::setorder(partition[i]);
		//for (i = 0; i < nnd; i++) {
		//	node3d & u = node[i];
		//	partition3d.add(mesh_part3d(u,partition[0][u.x],
		//			partition[1][u.y],partition[2][u.z]));
		//}
		//assert(partition3d.length() == nnd);
		for (i = 0; i < nnd; i++)
			//n2o[i] = i;
			//n2o[i] = partition3d.getorder(i);
			n2o[i] = node.getorder(i);
		//e = first;
		//while (e) {
		//	for (i = 0; i < e->nnd; i++) {
		//		for (j = 0; j < e->nnd; j++)
		//			rcm_order.append(n2o[e->l2g(i)],n2o[e->l2g(j)]);
		//	}
		//	e = e->next;
		//}
		//rcm_order.compute(n2o);
		//for (i = 0; i < nnd; i++) {
		//	node3d & u = node[i];
		//	printf("%3d %3d %+04.1f %+04.1f %+04.1f\n",i,n2o[i],double(u.x),double(u.y),double(u.z));
		//}

		timeme("\nBuilding external forces...");
		for (i = 0; i < f_ext.length(); i++) {
			const mesh_bc & f = f_ext[i];
			F(n2o[f.n])(f.i) = f.d;
		}

		timeme("\nSetting approximation of displacement...");
		for (i = 0; i < nnd; i++) {
			node3d & u = node[i];
			U(n2o[i])(0) = u.ux; U(n2o[i])(1) = u.uy; U(n2o[i])(2) = u.uz;
		}

		timeme("\nSetting fixed BCs...");
		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			// Each bit in fixed means that DOF is fixed.
			node[u.n].fixed |= (1 << u.i);
			// The last bit means at least one DOF is non-zero.
			if (u.d != 0.0) {
				node[u.n].fixed |= (1 << NDOF);
				//n = u.n;
				P(n2o[u.n])(u.i) = u.d;
			}
		}

		timeme("\nAssembling Stiffness Matrix...");
		// Loop over the elements building the stiffness matrix.
		e = first;
		while (e != nullptr) {
			smatrix_elem * ke = e->stiffness();
			if (ke == nullptr)
				return false;
			// Fix up the element matrix to account for displacement BCs
			for (i = 0; i < ke->nnd; i++) {
				const unsigned gi = e->l2g(i);
				const unsigned ni = n2o[gi];
				// At least one of P(ni) is non-zero, so apply the
				// fixup to F for this element.
				if (node[gi].fixed & (1 << NDOF)) {
					for (j = 0; j < ke->nnd; j++) {
						const unsigned nj = n2o[e->l2g(j)];
						if (i > j)
							F(nj) -=  (*ke)(j,i)*P(ni);
						else
							F(nj) -= ~(*ke)(i,j)*P(ni);
					}
				}
			}
			for (i = 0; i < ke->nnd; i++) {
				const unsigned gi = e->l2g(i);
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
				const unsigned gi = e->l2g(i);
				const unsigned ni = n2o[gi];
				// If this node is completely fixed, don't add anything.
				if ((node[gi].fixed & ((1<<NDOF)-1)) == ((1<<NDOF)-1))
					continue;
				for (j = i; j < ke->nnd; j++) {
					const unsigned gj = e->l2g(j);
					if ((node[gj].fixed & ((1<<NDOF)-1)) == ((1<<NDOF)-1))
						continue;
					K.append(ni,n2o[gj],(*ke)(i,j));
				}
			}
			delete ke;
			e = e->next;
		}

		timeme("\nFinalising BCs...");
		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			n = n2o[u.n];
			F(n)(u.i) = U(n)(u.i) = u.d;
			K.diag[n].K(u.i,u.i) = 1.0;
		}

		timeme("\nTidying up...");
		K.tidy();
		timeme("\nCopying...");
		smatrix M(K);
		timeme("\nComputing fill in...");
		M.fillin();
		timeme("\nComputing incomplete Cholesky...");
		M.incchol();
		timeme("\nBeginning CG...");

		unsigned it = 0;
		double r = 0.0, ro = 0.0;
		for (i = 0; i < nnd; i++)
			r += tmatrix_scalar<double>(~F(i)*F(i));
		const double ri = r;
		double a;
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
			if (r < tol*ri || it >= nnd*NDOF) {
				timeme("\n");
				break;
			}
			// Forward substitution of IC pre-conditioner.
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
			// Backward substitution of IC pre-conditioner.
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
			node3d & u = node[i];
			u.setdisp(U(n2o[i]));
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
operator | (const mesh::dof l, const mesh::dof r) noexcept
{
	return static_cast<mesh::dof>(unsigned(l) | unsigned(r));
}
inline mesh::dof
operator & (const mesh::dof l, const mesh::dof r) noexcept
{
	return static_cast<mesh::dof>(unsigned(l) & unsigned(r));
}

inline mesh::bcplane
operator | (const mesh::bcplane l, const mesh::bcplane r) noexcept
{
	return static_cast<mesh::bcplane>(unsigned(l) | unsigned(r));
}
inline mesh::bcplane
operator & (const mesh::bcplane l, const mesh::bcplane r) noexcept
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

//#define QUAD

int run = 0;
bool isvar = false;

/*
 * The basics of a layer used for the FEM.
 */
struct femlayer {
	fixed<8> top;
	fixed<8> bot;
	fixed<8> etop;
	fixed<8> ebot;
	material mat;

	femlayer() noexcept
	  : top(0.0), bot(0.0), etop(0.0), ebot(0.0), mat(0.0,0.0) {
	}
	femlayer(double t, double b, double e, double v) noexcept
	  : top(t), bot(b), etop(0.0), ebot(0.0), mat(e,v) {
	}
};

#define NUM_VAR 3
static sset<femlayer> layer;
double * L[NUM_VAR];

double
node_depth_callback(const coord3d & c, const mesh * FEM, const material * mat) noexcept
{
	unsigned i = 0, n;

	while (&(layer[i].mat) != mat)
		i++;
	double z = -c.z - layer[i].top;
	double t = layer[i].top, b = layer[i].bot;
	const double h = b-t;
	if (isvar)
		n = FEM->getorderofnode(coord3d(c.x,c.y,0.0));
	else
		n = FEM->getorderofnode(coord3d(0.0,0.0,0.0));
	if (i == 0) {
		t += L[0][n];
		b += L[1][n];
	} else if (i == 1) {
		t += L[1][n];
		b += (L[1][n]+L[2][n])/2;
	} else if (i == 2) {
		t += (L[1][n]+L[2][n])/2;
		b += L[2][n];
	} else if (i == 3) {
		t += L[2][n];
	}
	return -t+(t-b)*z/h;
}

double
node_emod_callback(const coord3d & c, const mesh * FEM, const material * mat) noexcept
{
	unsigned i = 0, n;

	while (&(layer[i].mat) != mat)
		i++;
	if (isvar)
		n = FEM->getorderofnode(coord3d(c.x,c.y,0.0));
	else
		n = FEM->getorderofnode(coord3d(0.0,0.0,0.0));
	return L[(i>1?i-1:i)+3][n];
}

inline double
circlearea(double x, double y, double r) noexcept
{
	x = fabs(x); y = fabs(y);
	if (hypot(x,y) >= r)
		return 0.0;
	const double t = M_PI_2 - asin(x/r) - asin(y/r);
	const double xr = (sqrt(r*r-y*y)-x)*y/2;
	const double yr = (sqrt(r*r-x*x)-y)*x/2;
	return (r*r*t/2) - xr - yr;
}

inline double
blockarea(double x1, double x2, double y1, double y2, double r) noexcept
{
	const double h1 = hypot(x1,y1), h2 = hypot(x2,y1);
	const double h3 = hypot(x1,y2), h4 = hypot(x2,y2);
	if (MIN(MIN(h1,h2),MIN(h3,h4)) >= r)
		return 0.0;
	if (MAX(MAX(h1,h2),MAX(h3,h4)) <= r)
		return fabs(x1-x2)*fabs(y1-y2);
	if (x1*x2 < 0.0) {
		const double y = MAX(fabs(y1),fabs(y2));
		return 2*circlearea(0.0,y,r)-circlearea(x1,y,r)
					-circlearea(x2,y,r);
	}
	if (y1*y2 < 0.0) {
		const double x = MAX(fabs(x1),fabs(x2));
		return 2*circlearea(x,0.0,r) - circlearea(x,y1,r)
					- circlearea(x,y2,r);
	}
	x1 = fabs(x1); x2 = fabs(x2); if (x2 < x1) swap(x1,x2);
	y1 = fabs(y1); y2 = fabs(y2); if (y2 < y1) swap(y1,y2);
	return circlearea(x1,y1,r) - circlearea(x1,y2,r)
			- circlearea(x2,y1,r) + circlearea(x2,y2,r);
}

/*
 * Bessel functions from numerical recipes
 */
static double
bessi0(const double x) noexcept
{
	double ax, y;

	if ((ax = fabs(x)) < 3.75) {
		y = x/3.75;
		y *= y;
		y = 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
		    +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
	} else {
		y = 3.75/ax;
		y = 0.39894228+y*(0.1328592e-1
		    +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
		    +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
		    +y*0.392377e-2)))))));
		y *= (exp(ax)/sqrt(ax));
	}
	return y;
}

static double
bessi1(double x) noexcept
{
	double ax, y;

	if ((ax = fabs(x)) < 3.75) {
		y = x/3.75;
		y *= y;
		y = ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
			+y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
	} else {
		y = 3.75/ax;
		y = 0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
			+y*(0.163801e-2+y*(-0.1031555e-1+y*(0.2282967e-1
			+y*(-0.2895312e-1+y*(0.1787654e-1-y*0.420059e-2)))))));
		y *= (exp(ax)/sqrt(ax));
	}
	return (x < 0.0 ? -y : y);
}

static double
bessk0(const double x)
{
	double y;

	if (x <= 2.0) {
		y = x*x/4.0;
		y = (-log(x/2.0)*bessi0(x))+(-0.57721566+y*(0.42278420
		    +y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
		    +y*(0.10750e-3+y*0.74e-5))))));
	} else {
		y = 2.0/x;
		y = (exp(-x)/sqrt(x))*(1.25331414+y*(-0.7832358e-1
		    +y*(0.2189568e-1+y*(-0.1062446e-1+y*(0.587872e-2
		    +y*(-0.251540e-2+y*0.53208e-3))))));
	}
	return y;
}

static double
bessk1(const double x)
{
	double y;

	if (x <= 2.0) {
		y = x*x/4.0;
		y = (log(x/2.0)*bessi1(x))+(1.0/x)*(1.0+y*(0.15443144
			+y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
			+y*(-0.110404e-2+y*(-0.4686e-4)))))));
	} else {
		y = 2.0/x;
		y = (exp(-x)/sqrt(x))*(1.25331414+y*(0.23498619
			+y*(-0.3655620e-1+y*(0.1504268e-1+y*(-0.780353e-2
			+y*(0.325614e-2+y*(-0.68245e-3)))))));
	}
	return y;
}

static double
bessk(int k, const double x)
{
	double bk, bkm;

	k = abs(k);
	if (k == 0)
		return bessk0(x);
	if (k == 1)
		return bessk1(x);
	bkm = bessk0(x);
	bk = bessk1(x);
	for (int j = 1; j < k; j++) {
		bkm = bkm+2*j*bk/x;
		swap(bkm,bk);
	}
	return bk;
}

static constexpr double
gamma(const int n)
{
	double g = 1.0;

	for (int c = 1; c < n; c++)
		g = g * c;
	return g;
}

static double
fn_matern(const int k, const double h, const double r)
{
	double hr;

	if (h == 0.0)
		return 1.0;
	if (h > 600*r)
		return 0.0;
	hr = h/r;
	return pow(2.0,1.0-k)/gamma(k)*pow(hr,k)*bessk(k,hr);
}

static double
fn_whittle(const double h, const double r)
{
	double hr;

	if (h == 0.0)
		return 1.0;
	if (h > 600*r)
		return 0.0;
	hr = h/r;
	return (hr*bessk1(hr));
}

static double
fn_ellipse(const double t, const double a, const double b) noexcept
{
	return a*b/sqrt(pow(a*cos(t),2)+pow(b*sin(t),2));
}

/*
 * A region (in the horizontal plane) which has been or is being
 * added to the mesh.  It is described by the 'stops' where there
 * are mesh points.
 */
struct region {
	cset<fixed<8> > xstop;
	cset<fixed<8> > ystop;

	region() noexcept {
	}
	fixed<8> xm() noexcept {
		return xstop[0];
	}
	fixed<8> xp() noexcept {
		return xstop[xstop.length()-1];
	}
	fixed<8> ym() noexcept {
		return ystop[0];
	}
	fixed<8> yp() noexcept {
		return ystop[ystop.length()-1];
	}
	bool overlaps(double x1, double y1, double x2, double y2) noexcept {
		fixed<8> fxm(x1), fym(y1), fxp(x2), fyp(y2);
		return (fxm < xp() && fxp > xm() && fym < yp() && fyp > ym());
	}
	bool overlaps(region & r) noexcept {
		return (r.xm() < xp() && r.xp() > xm()
				 && r.ym() < yp() && r.yp() > ym());
	}
};

struct region_list : sset<region> {
	double xm() noexcept {
		fixed<8> x_m = 0;
		for (unsigned i = 0; i < length(); i++)
			x_m = MIN(x_m,(*this)[i].xm());
		return x_m;
	}
	double xp() noexcept {
		fixed<8> x_p = 0;
		for (unsigned i = 0; i < length(); i++)
			x_p = MAX(x_p,(*this)[i].xp());
		return x_p;
	}
	double ym() noexcept {
		fixed<8> y_m = 0;
		for (unsigned i = 0; i < length(); i++)
			y_m = MIN(y_m,(*this)[i].ym());
		return y_m;
	}
	double yp() noexcept {
		fixed<8> y_p = 0;
		for (unsigned i = 0; i < length(); i++)
			y_p = MAX(y_p,(*this)[i].yp());
		return y_p;
	}
	bool overlaps(double x1, double y1, double x2, double y2) noexcept {
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
static void
result_mesh(const fset<const element *> & e, const fset<point3d> & c)
{
	fset<pavedata> data(c.length(),c.length());

	char fname[128];
	sprintf(fname,"mesh.bin");
	FILE * f = fopen(fname,"wb");
	for (unsigned i = 0; i < e.length(); i++) {
		e[i]->results(c,data);
		unsigned il = 0;
		while (&(layer[il].mat) != &(e[i]->mat))
			il++;
		const double ild = double(il);
		for (unsigned j = 0; j < data.length(); j++) {
			fwrite(&(data[j].x),sizeof(double),3,f);
			fwrite(&ild,sizeof(double),1,f);
		}
	}
	fclose(f);
}

/*
 * This outputs the results for all the elements in the list
 * Obviously, the list must be ones returned by add().
 */
static void
results(const fset<const element *> & e, const fset<point3d> & c,
		const char * dname)
{
	fset<pavedata> data(c.length(),c.length());

	char fname[128];
	sprintf(fname,"%s_%s_%05d.bin",dname,(isvar?"3d":"1d"),run);
	FILE * f = fopen(fname,"wb");
	for (unsigned i = 0; i < e.length(); i++) {
		e[i]->results(c,data);
		unsigned il = 0;
		while (&(layer[il].mat) != &(e[i]->mat))
			il++;
		const double ild = double(il);
		for (unsigned j = 0; j < data.length(); j++) {
			fwrite(&(data[j].x),sizeof(double),3,f);
			fwrite(&ild,sizeof(double),1,f);
			fwrite(&(data[j].deflgrad[0]),sizeof(double),1,f);
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
		const double ild = double(il);
		for (unsigned j = 0; j < data.length(); j++) {
			data[j].z *= -1;
			pavedata d(le.result(data[j],il));
			d.z *= -1; d.data[4][2] *= -1;
			swap(d.data[3][0],d.data[3][2]);
			swap(d.data[8][0],d.data[8][2]);
			d.data[1][1] *= -1; d.data[1][2] *= -1;
			d.data[6][1] *= -1; d.data[6][2] *= -1;
			fwrite(&(d.x),sizeof(double),3,f);
			fwrite(&ild,sizeof(double),1,f);
			fwrite(&(layer[il].mat.emod),sizeof(double),1,f);
			fwrite(&(d.data[0][0]),sizeof(double),27,f);
		}
	}
	fclose(f);
}

#ifdef BUILD
int
main(int argc, char *argv[])
{
	timeme();

	LEsystem test;
	test.addlayer(  83.5,4874.7880e3,0.40);
	test.addlayer( 186.5, 322.1802e3,0.35);
	test.addlayer( 186.0, 322.1802e3,0.35);
	test.addlayer(1544.0, 246.4247e3,0.35);

	//test.addlayer(  83.5,9103.7370e3,0.40);
	//test.addlayer( 186.5, 684.6993e3,0.35);
	//test.addlayer( 186.0, 684.6993e3,0.35);
	//test.addlayer(1544.0, 188.9532e3,0.35);

	//test.addlayer(  83.5,2596.5850e3,0.40);
	//test.addlayer( 186.5, 684.6993e3,0.35);
	//test.addlayer( 186.0, 684.6993e3,0.35);
	//test.addlayer(1544.0, 246.4247e3,0.35);

	layer.empty();
	for (unsigned i = 0; i < test.layers(); i++) {
		const LElayer & l = test.layer(i);
		layer.add(femlayer(l.top(),l.bottom(),l.emod(),l.poissons()));
	}

	double x = 0.0, y = 0.0, z = 0.0, zmax = 0.0, rmax = 0.0;
	double dx, dy, dz;
#ifdef QUAD
	double delta = 10;
#else
	double delta = 4;
#endif
	static constexpr const double edge = 4096;
	unsigned step = 4;
	mesh FEM;
	fset<coord3d> coord(8);
	cset<const element *> along, trans, rtrans, mesh;
	cset<const element *> sf, ac, ab, sg;
	sset<fixed<8> > stop;
	sset<int> level;
	region_list filled, filling;

	//test.addload(point2d(0.0,0.0),0.0,690.0,100.0);
	test.addload(point2d(-160.0,0.0),0.0,690.0,100.0);
	test.addload(point2d( 160.0,0.0),0.0,690.0,100.0);

	// Find the step size, which depends on our tires.
	for (unsigned i = 0; i < test.loads(); i++) {
		x += test.getload(i).x;
		y += test.getload(i).y;
		const double r = test.getload(i).radius();
		if (r > rmax)
			rmax = r;
	}
	x /= test.loads(); y /= test.loads();
	assert(x == 0.0);
	assert(y == 0.0);
#ifdef QUAD
	while (step*delta <= 1.6*rmax)
		step += 2;
#else
	while ((step-2)*delta < rmax)
		step *= 2;
#endif

	// Start with the tire grid.
	dx = delta; dy = delta;
	for (unsigned i = 0; i < test.loads(); i++) {
		paveload l = test.getload(0);
		double lx = l.x, ly = l.y;
		const double lr = l.radius(), F = -l.pressure();
		lx = 2*dx*(lx < 0 ? floor(lx/dx/2) : ceil(lx/dx/2));
		ly = 2*dy*(ly < 0 ? floor(ly/dy/2) : ceil(ly/dy/2));
		l.x = lx; l.y = ly;
		test.removeload(0);
		test.addload(l);
		for (x = lx - (step-2)*dx; x < lx + (step-2)*dx; x += 2*dx) {
			for (y = ly -(step-2)*dy; y < ly + (step-2)*dy; y += 2*dy) {
				const double ba = blockarea(x-lx,x+2*dx-lx,y-ly,y+2*dx-ly,lr);
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
#ifdef QUAD
		// Add the tire loads
		for (x = lx - (step-2)*dx; x <= lx + (step-2)*dx; x += 2*dx) {
			for (y = ly - (step-2)*dy; y <= ly + (step-2)*dy; y += 2*dy) {
				double ba = blockarea(x-lx,x+2*dx-lx,y-ly,y+2*dx-ly,lr);
				if (ba == 0.0)
					continue;
				double cx = x+dx-lx, cy = y+dy-ly;
				for (unsigned pi = 0; pi < 3; pi++) {
					for (unsigned pj = 0; pj < 3; pj++) {
						unsigned p = FEM.hasnode(coord3d(x+pi*dx,y+pj*dy,0.0));
						if (p == UINT_MAX)
							continue;
						node3d n = FEM.getnode(p);
						double f = 0.0;
						for (unsigned gi = 0; gi < 4; gi++) {
							for (unsigned gj = 0; gj < 4; gj++) {
								double gx = gl_4[gi][0];
								double gy = gl_4[gj][0];
								double gr = hypot(cx+gx*dx,cy+gy*dy);
								switch (pi) {
								case 0: gx = 0.5*gx*(gx-1); break;
								case 1: gx = 1-gx*gx;       break;
								case 2: gx = 0.5*gx*(gx+1); break;
								}
								switch (pj) {
								case 0: gy = 0.5*gy*(gy-1); break;
								case 1: gy = 1-gy*gy;       break;
								case 2: gy = 0.5*gy*(gy+1); break;
								}
								if (gr <= lr)
									f += gl_4[gi][1]*gl_4[gj][1]*gx*gy;
							}
						}
						if (f == 0.0)
							continue;
						f *= F*dx*dy;
						FEM.add_fext(n,mesh::Z,f);
					}
				}
			}
		}
#else
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
				const double f = F*blockarea(x-fxm-lx,x+fxp-lx,y-fym-ly,y+fyp-ly,lr);
				if (fabs(f) > 0.0)
					FEM.add_fext(n,mesh::Z,f);
			}
		}
#endif
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
			const double v = layer[j-1].mat.v;
			const double E = layer[j-1].mat.emod;
			const double v1 = layer[j].mat.v;
			const double E1 = layer[j].mat.emod;
			he *= pow(E/E1*(1+v1*v1)/(1+v*v),1.0/3.0);
		}
		layer[i].ebot = double(layer[i].etop) + he;
	}
	const double zstep = 1;
#ifdef QUAD
	int zi = -2;
#else
	int zi = 0;
#endif
	// add the depth stops.
	dz = rmax/zstep;
	zmax = layer[layer.length()-1].ebot; z = layer[0].etop;
	stop.add(z); level.add(MAX(zi,0));
	while (z < zmax) {
		const unsigned zscale = 3;
		for (; z < zscale*zstep*dz; z += dz) {
			printf("z: %g\tzmax: %g\tdz: %g\tzi: %i\n",z,zmax,dz,zi);
			stop.add(z+dz); level.add(MAX(zi,0));
		}
		dz += rmax/zstep; zi++;
	}
	for (unsigned j = 0; j < stop.length(); j++) {
		stop[j] = double(stop[j])*zmax/z;
		if (j > 0 && level[j] > level[j-1]+1)
			level[j] = level[j-1]+1;
	}
	// now scale them to the layers
	for (unsigned i = 0, t = 0, b; i < layer.length(); i++) {
		const femlayer & l = layer[i];
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
		//  s = 3; // Make sure we have three nodes in the top layer.
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
		// Now scale these stops to match the layer.
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

	// Add 3mm thick layers to smooth stresses...
	//stop.add(1,2); level.add(1,level[0]);
#ifndef QUAD
	for (unsigned i = 0, b = 0; i < layer.length()-1; i++) {
		const femlayer & l = layer[i];
		while (b < stop.length() && stop[b] < l.bot)
			b++;
		if (i == 0 || i == 1) {
			level[b] = level[b+1];
			stop.add(b,double(stop[b])-3); level.add(b,level[b]);
		}
		if (i == 2) {
			stop.add(b+1,double(stop[b])+3); level.add(b+1,level[b-1]);
		}
	}
#endif

	// Now add the elements from below the tire, working outwards.
	zi = 0;
	while (zi <= level[level.length()-1] || filled.xp() < edge) {
		for (unsigned i = 0, j = 1; j < stop.length() && level[j] <= zi; j++) {
			while (i < layer.length() && stop[j] > layer[i].bot)
				i++;
			assert(i < layer.length());
			const double z1 = stop[j  ];
			const double z2 = stop[j-1];
			printf("%.0f\t%.0f\t%f\t%i\t%i\t%i\t%i\n",filling.xp(),filled.xp(),z1,i,zi,j,level[j]);
			element::element_t var = (i == 0 ?
					element::variable26 : element::variable26);
			if (z1-z2 < 5)
				var = element::variable18;
			const element::element_t inf = (var == element::variable34 ?
					element::infinite16 : var == element::variable26 ?
					element::infinite12 : element::infinite8);
#ifdef QUAD
			if (zi <= 0)
				var = (var == element::variable34 ?
					element::block36 : var == element::variable26 ?
					element::block27 : element::block18);
			if (zi == 1)
				var = (var == element::variable34 ?
					element::adaptor34 : var == element::variable26 ?
					element::adaptor26 : element::adaptor18);
#endif
			for (unsigned k = 0; k < filling.length(); k++) {
				const region & r = filling[k];
				for (unsigned xi = 1; xi < r.xstop.length(); xi++) {
					for (unsigned yi = 1; yi < r.ystop.length(); yi++) {
						const double x1 = double(r.xstop[xi-1]);
						const double x2 = double(r.xstop[xi]);
						const double y1 = double(r.ystop[yi-1]);
						const double y2 = double(r.ystop[yi]);
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
						mesh.add(e);
						if (x1 == 0.0)
							along.add(e);
						if (y1 == 0.0)
							trans.add(e);
						if (y2 == 0.0)
							rtrans.add(e);
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
#ifdef QUAD
			if (zi == 0)
				step *= 2;
			else {
#endif
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
#ifdef QUAD
			}
#endif
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
	//      mesh::X,0.0);
	//FEM.add_bc_plane(mesh::X,mesh::at|mesh::above, edge,
	//      mesh::X,0.0);
	//FEM.add_bc_plane(mesh::Y,mesh::at|mesh::below,-edge,
	//      mesh::Y,0.0);
	//FEM.add_bc_plane(mesh::Y,mesh::at|mesh::above, edge,
	//      mesh::Y,0.0);
	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,-zmax,
			mesh::X|mesh::Y|mesh::Z,0.0);

	// Generate our random fields.
	unsigned np = 0;
	while (double(FEM.getorderednode(np).z) == 0.0)
		np++;
	printf("np: %i\n",np);
	for (unsigned l = 0; l < NUM_VAR; l++) {
		L[l] = new double[np];

		char fname[128];
		sprintf(fname,"cor%d.tri",l);
		FILE * f = fopen(fname,"rb");
		if (f != NULL) {
			fclose(f);
			continue;
		}
		f = fopen(fname,"wb");
		double * C = new double[T_SIZE(np)];
		if (L[l] == nullptr || C == nullptr) {
			event_msg(EVENT_ERROR,"Ooops, out of memory...");
			return 0;
		}
		for (unsigned i = 0; i < np; i++) {
			const node3d & p1 = FEM.getorderednode(i);
			for (unsigned j = i; j < np; j++) {
				const node3d & p2 = FEM.getorderednode(j);
				x = p1.x-p2.x; y = p1.y-p2.y;
				const double r = fabs(hypot(x,y));
				const double t = atan2(y,x); // zero is longitudinal.
				switch (l) {
				case 0:
					C[T_IDX(i,j)] = (pow(8.886948,2)*exp(-pow(r/4420.595,2))
						+ pow(9.351780,2)*cos(2*M_PI*r/36322.982)*exp(-pow(2*M_PI*r/36322.982,2)))/(pow(8.886948,2)+pow(9.351780,2));
					break;
				case 1:
					C[T_IDX(i,j)] = fn_matern(2,r,fn_ellipse(t,1882.1698,632.9367));
					break;
				case 2:
					C[T_IDX(i,j)] = fn_whittle(r,fn_ellipse(t,3798.156,1198.276));
					break;
				default:
					assert(false);
				}
				if (i != j)
					C[T_IDX(i,j)] *= .999;
			}
		}
		decmp_chol_tri(np,C);
		fwrite(C,sizeof(double),T_SIZE(np),f);
		fclose(f);
		delete [] C;
	}

	const int startrun = (argc > 1 && argv[1] ? atoi(argv[1]) : 0);
	const int countrun = (argc > 2 && argv[2] ? atoi(argv[2]) : 250);
	run = startrun;
	isvar = false;
	rng RNG(run+clock());
	while (true) {
		isvar = !isvar;
		if (isvar)
			run++;
		if (run > startrun+countrun)
			break;

		for (unsigned l = 0; isvar && l < NUM_VAR; l++) {
			double * A = new double[np];
			double * C = new double[T_SIZE(np)];
			if (A == nullptr || C == nullptr) {
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
				switch (l) {
				case 0:
					L[l][i] *= sqrt(pow(8.886948,2)+pow(9.351780,2));
					break;
				case 1:
					// include the randomness from the previous layers.
					L[l][i] = L[l-1][i]	+ 9.0005380*L[l][i];
					break;
				case 2:
					L[l][i] = L[l-1][i]	+ 32.040745*L[l][i];
					break;
				default:
					assert(false);
				}
			}
			delete [] C;
			delete [] A;
		}

		//for (unsigned i = 0; i < np; i++) {
		//	const node3d & p1 = FEM.getorderednode(i);
		//	printf("%6.1f %6.1f:",double(p1.x),double(p1.y));
		//	for (unsigned l = 0; l < NUM_VAR; l++)
		//		printf(" %8.2f",L[l][i]);
		//	printf("\n");
		//}

		for (unsigned i = 0; i < FEM.getnodes(); i++) {
			node3d n = FEM.getnode(i);
			n.ux = 0.0; n.uy = 0.0; n.uz = 0.0;
			FEM.updatenode(n);
		}
		try {
			FEM.solve(1e-25);
		} catch (...) {
			isvar = false;
			run--;
			continue;
		}

		fset<point3d> poly(4), face(8);
		//face[0] = point3d(-1,-1,-1);
		//face[1] = point3d(-1,-1, 1);
		//face[2] = point3d(-1, 1,-1);
		//face[3] = point3d(-1, 1, 1);
		//face[4] = point3d( 1,-1,-1);
		//face[5] = point3d( 1,-1, 1);
		//face[6] = point3d( 1, 1,-1);
		//face[7] = point3d( 1, 1, 1);
		//result_mesh(mesh,face);
		face[0] = point3d(-1,-1,    -1);
		face[1] = point3d(-1,-1,-1/3.0);
		face[2] = point3d(-1,-1, 1/3.0);
		face[3] = point3d(-1,-1,     1);
		face[4] = point3d( 1,-1,     1);
		face[5] = point3d( 1,-1, 1/3.0);
		face[6] = point3d( 1,-1,-1/3.0);
		face[7] = point3d( 1,-1,    -1);
		results(trans,face,"face");
		//LEresults(trans,face,"face",test);
		face[0] = point3d(-1, 1,    -1);
		face[1] = point3d(-1, 1,-1/3.0);
		face[2] = point3d(-1, 1, 1/3.0);
		face[3] = point3d(-1, 1,     1);
		face[4] = point3d( 1, 1,     1);
		face[5] = point3d( 1, 1, 1/3.0);
		face[6] = point3d( 1, 1,-1/3.0);
		face[7] = point3d( 1, 1,    -1);
		results(rtrans,face,"rface");
		face[0] = point3d(-1,-1,    -1);
		face[1] = point3d(-1,-1,-1/3.0);
		face[2] = point3d(-1,-1, 1/3.0);
		face[3] = point3d(-1,-1,     1);
		face[4] = point3d(-1, 1,     1);
		face[5] = point3d(-1, 1, 1/3.0);
		face[6] = point3d(-1, 1,-1/3.0);
		face[7] = point3d(-1, 1,    -1);
		results(along,face,"long");
		//LEresults(along,face,"long",test);
		poly[0] = point3d(-1,-1,-1);
		poly[1] = point3d(-1, 1,-1);
		poly[2] = point3d( 1, 1,-1);
		poly[3] = point3d( 1,-1,-1);
		results(ac,poly,"ac");
		//LEresults(ac,poly,"ac",test);
		results(ab,poly,"ab");
		//LEresults(ab,poly,"ab",test);
		poly[0] = point3d(-1,-1, 1);
		poly[1] = point3d(-1, 1, 1);
		poly[2] = point3d( 1, 1, 1);
		poly[3] = point3d( 1,-1, 1);
		results(sf,poly,"sf");
		//LEresults(sf,poly,"sf",test);
		results(sg,poly,"sg");
		//LEresults(sg,poly,"sg",test);
	}

	/*
	unsigned i, nnd = FEM.getnodes();
	test.removepoints();
	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z; // XXX
		//if (!(x == 0 || y == 0 || z == 0))
		//	continue;
		//if (!(x == 0))
		//	continue;
		if (fabs(x) > edge/4 || fabs(y) > edge/4 || fabs(z) > 800.0)
			continue;
		if (x < -edge || x > edge || y < -edge || y > edge || z <= -zmax)
			continue;
		test.addpoint(point3d(x,y,-z));
	}
	test.calculate(LEsystem::all);

	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z; // XXX
		//if (!(x == 0 || y == 0 || z == 0))
		//	continue;
		//if (!(x == 0))
		//	continue;
		if (fabs(x) > edge/4 || fabs(y) > edge/4 || fabs(z) > 800.0)
			continue;
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
	}
	*/

	timeme("\n");
	return 0;
}
#endif

double
node_depth_callback_test(const coord3d & c, const mesh *, const material *) noexcept
{
	return c.z;
}

double
node_emod_callback_test(const coord3d &, const mesh *, const material * mat) noexcept
{
	return mat->emod;
}

#ifdef NOBUILD
int
main()
{
	timeme();

	LEsystem test;
	test.addlayer( 260.0,5000e3,0.35);
	test.addlayer( 240.0, 200e3,0.40);
	test.addlayer( 260.0, 200e3,0.40);
	test.addlayer(5240.0,  50e3,0.45);
	test.addload(point2d(0.0,0.0),50.0e6,700.0,0.0);

	layer.empty();
	for (unsigned i = 0; i < test.layers(); i++) {
		const LElayer & l = test.layer(i);
		layer.add(femlayer(l.top(),l.bottom(),l.emod(),l.poissons()));
	}

	double x[] = {0, 5, 50, 100, 150.8, 200, 250, 300, 400, 500, 750, 1000, 2000, 4000 };
	double y[] = {0, 5, 50, 118.4, 150, 300, 500, 1000, 2000, 4000 };
	double z[] = {0,-5,-50,-100,-200,-255,-260,-265,-300,-400,-500,-600,-700,-755,-760,-765,-800,-900,-1000,-2000,-3000,-4000,-5000,-6000 };
	unsigned nx = sizeof(x)/sizeof(double);
	unsigned ny = sizeof(y)/sizeof(double);
	unsigned nz = sizeof(z)/sizeof(double);

	mesh FEM;
	fset<coord3d> coord(8);
	cset<const element *> along;
	cset<const element *> trans;
	cset<const element *> sf, ac, ab, sg;

	for (unsigned xi = 0; xi < nx; xi++)
		x[xi] = double(fixed<8>(-x[xi]));
	for (unsigned yi = 0; yi < ny; yi++)
		y[yi] = double(fixed<8>(-y[yi]));
	for (unsigned zi = 0; zi < nz; zi++)
		z[zi] = double(fixed<8>(z[zi]));

	for (unsigned zi = 1, i = 0; zi < nz; zi++) {
		while (i < layer.length() && -double(z[zi]) > double(layer[i].bot))
			i++;
		double z1 = z[zi];
		double z2 = z[zi-1];
		for (unsigned xi = 1; xi < nx; xi++) {
			for (unsigned yi = 1; yi < ny; yi++) {
				double x1 = x[xi-1], x2 = x[xi];
				double y1 = y[yi-1], y2 = y[yi];
				coord[0] = coord3d(x1,y1,z1);
				coord[1] = coord3d(x1,y2,z1);
				coord[2] = coord3d(x2,y1,z1);
				coord[3] = coord3d(x2,y2,z1);
				coord[4] = coord3d(x1,y1,z2);
				coord[5] = coord3d(x1,y2,z2);
				coord[6] = coord3d(x2,y1,z2);
				coord[7] = coord3d(x2,y2,z2);
				const element * e = FEM.add(element::block27,
						layer[i].mat,coord);
				if (x1 == 0.0)
					along.add(e);
				if (y1 == 0.0)
					trans.add(e);
				if (-double(z2) == double(layer[0].top))
					sf.add(e);
				if (-double(z1) == double(layer[0].bot))
					ac.add(e);
				if (-double(z1) == double(layer[1].bot))
					ab.add(e);
				if (-double(z2) == double(layer[3].top))
					sg.add(e);
			}
		}
	}
	FEM.add_bc_plane(mesh::X,mesh::at|mesh::above,x[0],
			mesh::X,0.0);
	FEM.add_bc_plane(mesh::X,mesh::at|mesh::below,x[nx-1],
			mesh::X,0.0);
	FEM.add_bc_plane(mesh::Y,mesh::at|mesh::above,y[0],
			mesh::Y,0.0);
	FEM.add_bc_plane(mesh::Y,mesh::at|mesh::below,y[ny-1],
			mesh::Y,0.0);
	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,z[nz-1],
			mesh::X|mesh::Y|mesh::Z,0.0);

	// Add the tire loads
	for (unsigned xi = 1; x[xi]+150.8 > -0.1; xi++) {
		for (unsigned yi = 1; y[yi]+118.4 > -0.1; yi++) {
			double x1 = x[xi-1], y1 = y[yi-1];
			double dx = double(x[xi]-x[xi-1])/2.0, dy = double(y[yi]-y[yi-1])/2.0;
			for (unsigned pi = 0; pi < 3; pi++) {
				for (unsigned pj = 0; pj < 3; pj++) {
					unsigned p = FEM.hasnode(coord3d(x1+pi*dx,y1+pj*dy,0.0));
					assert(p != UINT_MAX);
					node3d n = FEM.getnode(p);
					double f = 0.0;
					for (unsigned gi = 0; gi < 4; gi++) {
						for (unsigned gj = 0; gj < 4; gj++) {
							double gx = gl_4[gi][0];
							double gy = gl_4[gj][0];
							switch (pi) {
							case 0: gx = 0.5*gx*(gx-1); break;
							case 1: gx = 1-gx*gx;       break;
							case 2: gx = 0.5*gx*(gx+1); break;
							}
							switch (pj) {
							case 0: gy = 0.5*gy*(gy-1); break;
							case 1: gy = 1-gy*gy;       break;
							case 2: gy = 0.5*gy*(gy+1); break;
							}
							f += gl_4[gi][1]*gl_4[gj][1]*gx*gy;
						}
					}
					if (f == 0.0)
						continue;
					f *= -700*dx*dy;
					FEM.add_fext(n,mesh::Z,f);
				}
			}
		}
	}

	for (unsigned i = 0; i < FEM.getnodes(); i++) {
		node3d n = FEM.getnode(i);
		n.ux = 0.0; n.uy = 0.0; n.uz = 0.0;
		FEM.updatenode(n);
	}
	if (!FEM.solve(1e-25))
		return 0;

	fset<point3d> poly(4), face(8);
	face[0] = point3d(-1, 1,    -1);
	face[1] = point3d(-1, 1,-1/3.0);
	face[2] = point3d(-1, 1, 1/3.0);
	face[3] = point3d(-1, 1,     1);
	face[4] = point3d( 1, 1,     1);
	face[5] = point3d( 1, 1, 1/3.0);
	face[6] = point3d( 1, 1,-1/3.0);
	face[7] = point3d( 1, 1,    -1);
	results(trans,face,"face");
	LEresults(trans,face,"face",test);
	face[0] = point3d( 1,-1,    -1);
	face[1] = point3d( 1,-1,-1/3.0);
	face[2] = point3d( 1,-1, 1/3.0);
	face[3] = point3d( 1,-1,     1);
	face[4] = point3d( 1, 1,     1);
	face[5] = point3d( 1, 1, 1/3.0);
	face[6] = point3d( 1, 1,-1/3.0);
	face[7] = point3d( 1, 1,    -1);
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

	timeme("\n");
	return 0;
}
#endif

#ifdef NOBUILD
int
main()
{
	timeme();
	printf("Constructing Mesh...");

	layer.empty();
	layer.add(femlayer(0,0,1000,0.2));

	const double cube[3][2] = {{-10, 10}, {-10, 10}, {-10, 0}};
	const unsigned ndiv[3] = {10, 10, 10};
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
		unsigned p = 1; //(1 << k); // 1;
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
				const element * e = FEM.add(s,layer[0].mat,coord);
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
	//  node3d n = FEM.getorderednode(i);
	//  n.ux =  0.0006*(double(n.x)+10.0)/20.0;
	//  n.uy =  0.0006*(double(n.y)+10.0)/20.0;
	//  n.uz = -0.0015*(double(n.z)+10.0)/10.0;
	//  //n.ux = 0.0; n.uy = 0.0; n.uz = 0.0;
	//  FEM.updatenode(n);
	//}
	FEM.solve(1e-25);

	timeme("\nOutput...\n");

	for (i = 0; i < FEM.getnodes(); i++) {
		const node3d & n = FEM.getorderednode(i);
		double x = n.x;
		double y = n.y;
		double z = n.z;
		//if (y != 0.0 || x != 0.0)
		//  continue;
		double ux = n.ux;
		double uy = n.uy;
		double uz = n.uz;
		j = FEM.hasnode(n);
		printf("Node %3i: (%6.2f,%6.2f,%6.2f) = (%+f,%+f,%+f)\n",j,double(x),double(y),double(z),ux,uy,uz);
	}
	//const node3d & n = FEM.getnode(FEM.hasnode(coord3d(cube[0][1],cube[1][1],cube[2][1])));
	//printf("(%+f,%+f,%+f)\n",n.ux,n.uy,n.uz);

	fset<point3d> poly(4);
	poly[0] = point3d(-1,-1,-1);
	poly[1] = point3d(-1,-1,+1);
	poly[2] = point3d(-1,+1,+1);
	poly[3] = point3d(-1,+1,-1);
	results(face,poly,"test");

	timeme(" Done.\n");

	return 0;
}
#endif

#ifdef NOBUILD
int
main()
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
#endif

#ifdef NOBUILD
int
main()
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
	printf("%f %f\n",a,M_PI*100*100);

	printf("\n\n");
	for (x = 0.0; x < 20.0; x+= 0.001) {
		printf("%f: %f %f %f %f %f %f %f %f %f %f\n",x,bessi0(x),bessi1(x),bessk0(x),bessk1(x),fn_matern(2,x,1),fn_whittle(x,1),pow(2.0,1-2),gamma(2),pow(x,2),bessk(2,x));
	}

	return 0;
}
#endif
