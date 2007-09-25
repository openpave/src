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

	The Original Code is OpenPave.org Core Libraries.

	The Initial Developer of the Original Code is OpenPave.org.

	Portions Copyright (C) 2007 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2007/09/07 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#define _EVENT_IMP
#define _PROGRESS_IMP
#include "event.h"
#include "fixed.h"
#include "set.h"
#include "list.h"
#include "matrix.h"
#include "pavement.h"

#define NDOF 3
#define NDIM 3

/*
 * struct material property
 *
 * A simple wrapper class around an enum of material property names,
 * to enable them to be used in a set, so we can add more properties
 * as we go.  The next set will be to work on how to make this a 3D
 * spatial estimator.
 */
struct material_property {
	enum property_t {
		null,
		emod,
		poissons,
	} property;

	inline material_property()
	  : property(null) {
	}
	inline material_property(const property_t p)
	  : property(p) {
	}
	inline bool operator== (const material_property & p) {
		return (property == p.property);
	} 
};

struct material_property_value : public material_property {
	double value;

	inline material_property_value()
	  : material_property(), value(0.0) {
	}
	inline material_property_value(
		const material_property::property_t p, const double & v)
	  : material_property(p), value(v) {
	}
	inline operator double() {
		return value;
	} 
};

class material {
public:
	inline material() {
	}
	inline void setprop(const material_property::property_t p,
		const double & v) {
		props.add(material_property_value(p,v));
	}
	inline bool hasprop(const material_property::property_t p) const {
		return props.haskey(material_property(p)) != -1;
	}
	inline double getprop(const material_property::property_t p) const {
		if (!hasprop(p))
			event_msg(EVENT_ERROR,"Access to undefined material property!");
		return double(props[material_property(p)]);
	}
private:
	ksset<material_property,material_property_value> props;
};

/*
 * A mesh coordinate point (in 2D).
 */
struct node2d {
	fixed<8> x, y;

	node2d() {
	}
	node2d(double px, double py)
	  : x(px), y(py) {
	}
	node2d(const node2d & p)
	  : x(p.x), y(p.y) {
	}
	node2d(const point2d & p)
	  : x(p.x), y(p.y) {
	}
	~node2d () {
	}
	inline double distance(const node2d & p) {
		return hypot(double(x-p.x),double(y-p.y));
	}
	int compare(const node2d & p) const {
		if (x == p.x && y == p.y)
			return 0;
		if (x < p.x)
			return -1;
		if (x > p.x)
			return 1;
		return (y < p.y ? -1 : 1);
	}
	inline bool operator == (const node2d & p) const {
		return (x == p.x && y == p.y ? true : false);
	}
	inline bool operator != (const node2d & p) const {
		return (x != p.x || y != p.y ? true : false);
	}
	inline bool operator > (const node2d & p) const {
		return (compare(p) == 1 ? true : false);
	}
	inline bool operator >= (const node2d & p) const {
		return (compare(p) != -1 ? true : false);
	}
	inline bool operator < (const node2d & p) const {
		return (compare(p) == -1 ? true : false);
	}
	inline bool operator <= (const node2d & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * A simple point in 3D space.
 */
struct node3d : public node2d {
	fixed<8> z;

	node3d() : node2d() {
	}
	node3d(double px, double py, double pz)
	  : node2d(px,py), z(pz) {
	}
	node3d(const node3d & p)
	  : node2d(p), z(p.z) {
	}
	node3d(const point3d & p)
	  : node2d(p), z(p.z) {
	}
	~node3d () {
	}
	int compare(const node3d & p) const {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if (z < p.z)
			return -1;
		if (z > p.z)
			return 1;
		return node2d::compare(static_cast<const node2d &>(p));
	}
	bool operator == (const node3d & p) const {
		return (x == p.x && y == p.y && z == p.z ? true : false);
	}
	bool operator != (const node3d & p) const {
		return (x != p.x || y != p.y || z != p.z ? true : false);
	}
	inline bool operator > (const node3d & p) const {
		return (compare(p) == 1 ? true : false);
	}
	inline bool operator >= (const node3d & p) const {
		return (compare(p) != -1 ? true : false);
	}
	inline bool operator < (const node3d & p) const {
		return (compare(p) == -1 ? true : false);
	}
	inline bool operator <= (const node3d & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * struct gauss3d - a Gauss Point in 3D
 */
struct gauss3d : public point3d {
	double gw;
	
	gauss3d()
	  : point3d() {
	}
	gauss3d(double px, double py, double pz, double w)
	  : point3d(px,py,pz), gw(w) {
	}
};

/*
 * Templated Kronecker delta
 */ 
template<int I, int J>
static inline int delta()
{
	return (I == J ? 1 : 0);
}

/*
 * struct meta_stiffness - template meta class to assign the stiffness tensor
 */
template<unsigned I, unsigned J, unsigned K, unsigned L>
struct meta_stiffness {
	static inline void assignR(tmatrix<double,3,3> & E, const double & l,
			const double & m) {
		E(J,L) = l*delta<I,J>()*delta<K,L>()
				+ m*(delta<I,K>()*delta<J,L>()
				   + delta<I,L>()*delta<J,K>());
		meta_stiffness<I,J,K,L-1>::assignR(E,l,m);
	}
	static inline void assign(tmatrix<double,3,3> & E, const double & l,
			const double & m) {
		meta_stiffness<I,J,K,L>::assignR(E,l,m);
		meta_stiffness<I,J-1,K,L>::assign(E,l,m);
	}
};
template<unsigned I, unsigned J, unsigned K>
struct meta_stiffness<I,J,K,-1> {
	static inline void assignR(tmatrix<double,3,3> &, const double &,
			const double &) {}
};
template<unsigned I, unsigned K, unsigned L>
struct meta_stiffness<I,-1,K,L> {
	static inline void assign(tmatrix<double,3,3> &, const double &,
			const double &) {}
};

/*
 * pointstiff - Determine the consitutive level stiffness
 */
#if (NDOF == 3 && NDIM == 3)
template<int I, int K>
static inline void pointstiffness(tmatrix<double,3,3> & E,
		const double lamba, const double nu)
{
	meta_stiffness<I,2,K,2>::assign(E,lamba,nu);
}
#else
#error "We're only coded for 3D!"
#endif

class mesh;

/*
 * class smatrix_elem - an element stiffness matrix
 */
class smatrix_elem {
public:
	explicit smatrix_elem(const int n, const fset<int> & in)
	  : nnd(n), inel(in), K(0) {
		K = static_cast<tmatrix<double,NDOF,NDOF> *>
				(calloc(size(),sizeof(tmatrix<double,NDOF,NDOF>)));
		if (K == 0)
			event_msg(EVENT_ERROR,"Out of memory in smatrix_elem::smatrix_elem()!");
	}
	inline ~smatrix_elem() {
		if (K)
			free(K);
	}
	inline tmatrix<double,NDOF,NDOF> & operator() (int i, int j) const {
		if (i < j)
			swap(i,j);
		return K[index(i,j)];
	}

private:
	friend class mesh;
	const int nnd;
	const fset<int> & inel;
	tmatrix<double,NDOF,NDOF> * K;

	// Size of the triangular matrix storage
	inline int size() const {
		return nnd*(nnd+1)/2;
	}
	// Index into the triangular matrix storage
	inline int index(const int i, const int j) const {
		return i*(i+1)/2 + j;
	}
};

/*
 * Special tmatrix version of this...
 */
template<unsigned N, unsigned M>
double
inv_mul_gauss(tmatrix<double,N,N> & A, tmatrix<double,N,M> & B)
{
	double det = 1.0;
	unsigned i, j, k;

	for (i = 0; i < N; i++) {
		double pvt = A(i,i);
		if (fabs(pvt) < DBL_EPSILON) {
			for (j = i+1; j < N; j++) {
				if (fabs(pvt = A(j,i)) >= DBL_EPSILON)
					break;
			}
			if (j == N) {
				event_msg(EVENT_ERROR,"Singular matrix in inv_mul_gauss()!");
				det = 0.0;
				goto abort;
			}
			for (k = 0; k < N; k++)
				swap(A(j,k),A(i,k));
			for (k = 0; k < M; k++)
				swap(B(j,k),B(i,k));
		}
		det *= pvt;
		for (k = N-1; k > i; k--) {
			double tmp = A(k,i)/pvt;
			for (j = N-1; j > i; j--)
				A(k,j) -= tmp*A(i,j);
			for (j = 0; j < M; j++)
				B(k,j) -= tmp*B(i,j);
		}
	}
	for (i = N; i > 0; ) {
		for (--i, j = N-1; j > i; j--)
			for (k = 0; k < M; k++)
				B(i,k) -= A(i,j)*B(j,k);
		for (k = 0; k < M; k++)
			B(i,k) /= A(i,i);
	}
abort:
	return det;
}

/*
 * class element - a finite element
 */
class element : public listelement_o<mesh,element> {
public:
	const enum element_t {
		block8,
		block16
	} etype;

	element(mesh * o, element * p, const element_t t, const material & m)
	  : listelement_o<mesh,element>(o,p), etype(t), mat(m) {
	}

	virtual smatrix_elem * stiffness() = 0;
	
protected:
	friend class mesh;
	const material & mat;

	inline int addnode(const node3d & c);
	inline const node3d & getnode(const int i);
};

template<int N>
class element_N : public element {
public:
	element_N(mesh * o, element * p, const element_t t, const material & m,
			const fset<node3d> & c)
	  : element(o,p,t,m), inel(8,9) {
	  	assert(N == c.length());
		switch (etype) {
		case block8:
			for (int i = 0; i < 8; i++)
				inel[i] = addnode(c[i]);
			break;
		case block16:
			break;
		}
	}
	virtual smatrix_elem * stiffness() {
		int i, j, k, l, g;
		double rx, ry, rz;
		double e = mat.getprop(material_property::emod);
		double v = mat.getprop(material_property::poissons);
		double lambda = v*e/(1+v)/(1-2*v);
		double mu = e/2/(1+v);
	
		// Build the point stiffness tensor E_abcd
		tmatrix<double,3,3> E[3][3];
		pointstiffness<0,0>(E[0][0],lambda,mu);
		pointstiffness<0,1>(E[0][1],lambda,mu);
		pointstiffness<0,2>(E[0][2],lambda,mu);
		pointstiffness<1,0>(E[1][0],lambda,mu);
		pointstiffness<1,1>(E[1][1],lambda,mu);
		pointstiffness<1,2>(E[1][2],lambda,mu);
		pointstiffness<2,0>(E[2][0],lambda,mu);
		pointstiffness<2,1>(E[2][1],lambda,mu);
		pointstiffness<2,2>(E[2][2],lambda,mu);

		// Element nodal coords
		tmatrix<double,N,NDIM> xe;
		for (i = 0; i < N; i++) {
			const node3d & p = getnode(inel[i]);
			xe[i][0] = p.x; xe[i][1] = p.y; xe[i][2] = p.z;
		}

		ksset<point3d,gauss3d> gp(0,8);
		const double gp_2[2][2] = {{-1.0/sqrt(3.0), 1.0},
		                           {+1.0/sqrt(3.0), 1.0}};
		const double gp_A[4][2] = {{-(1.0+1.0/sqrt(3.0))/2.0, 1.0/2.0},
		                           {-(1.0-1.0/sqrt(3.0))/2.0, 1.0/2.0},
		                           {+(1.0-1.0/sqrt(3.0))/2.0, 1.0/2.0},
		                           {+(1.0+1.0/sqrt(3.0))/2.0, 1.0/2.0}};
		const double gp_3[3][2] = {{-sqrt(3.0/5.0), 5.0/9.0},
		                           {             0, 8.0/9.0},
		                           {+sqrt(3.0/5.0), 5.0/9.0}};
		//const double gp_4[4][2] =
		//	{{-sqrt(525.0+70.0*sqrt(30.0))/35.0, (18.0-sqrt(30.0))/36.0},
		//	 {-sqrt(525.0-70.0*sqrt(30.0))/35.0, (18.0+sqrt(30.0))/36.0},
		//	 {+sqrt(525.0-70.0*sqrt(30.0))/35.0, (18.0+sqrt(30.0))/36.0},
		//	 {+sqrt(525.0+70.0*sqrt(30.0))/35.0, (18.0-sqrt(30.0))/36.0}};
		switch (etype) {
		case block8:
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					for (k = 0; k < 2; k++) {
						gp.add(gauss3d(gp_2[i][0],gp_2[j][0],gp_2[k][0],
								gp_2[i][1]*gp_2[j][1]*gp_2[k][1]));
					}
				}
			}
			break;
		case block16:
			for (i = 0; i < 2; i++) {
				for (j = 0; j < 2; j++) {
					for (k = 0; k < 3; k++) {
						gp.add(gauss3d(gp_2[i][0],gp_2[j][0],gp_3[k][0],
								gp_2[i][1]*gp_2[j][1]*gp_3[k][1]));
					}
				}
			}
			break;
		}
	
		smatrix_elem * _K = new smatrix_elem(N,inel);
		if (_K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}
		smatrix_elem & K = *_K;
	
		for (g = 0; g < gp.length(); g++) {
			rx = gp[g].x; ry = gp[g].y; rz = gp[g].z;
			double gw = gp[g].gw;
			tmatrix<double,NDIM,N> dHdr;
			//tmatrix<double,N,1> H;
	
			// shape functions for 8-node 3D brick:
			//   N = 1/8*(1+-x)*(1+-y)*(1+-z);
			const double Hx_8[8] = {-1, -1, +1, +1, -1, -1, +1, +1};
			const double Hy_8[8] = {-1, +1, -1, +1, -1, +1, -1, +1};
			const double Hz_8[8] = {-1, -1, -1, -1, +1, +1, +1, +1};
	
			// shape functions for 16-node 3D brick:
			const double Hx1_16[16] = {-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1};
			const double Hy1_16[16] = {-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1};
			const double Hz0_16[16] = {-1,-1,-1,-1,+9,+9,+9,+9,+9,+9,+9,+9,-1,-1,-1,-1};
			const double Hz2_16[16] = {+9,+9,+9,+9,-9,-9,-9,-9,-9,-9,-9,-9,+9,+9,+9,+9};
			const double Hz1_16[16] = {-1,-1,-1,-1,-3,-3,-3,-3,+3,+3,+3,+3,+1,+1,+1,+1};
	
			switch (etype) {
			case block8:
				assert(N == 8);
				for (l = 0; l < N; l++) {
					dHdr[0][l] = Hx_8[l]*(1+Hy_8[l]*ry)*(1+Hz_8[l]*rz)/8;
					dHdr[1][l] = Hy_8[l]*(1+Hx_8[l]*rx)*(1+Hz_8[l]*rz)/8;
					dHdr[2][l] = Hz_8[l]*(1+Hx_8[l]*rx)*(1+Hy_8[l]*ry)/8;
					//H[0][l] =(1+Hx_8[l]*rx)*(1+Hy_8[l]*ry)*(1+Hz_8[l]*rz)/8;
				}
				break;
			case block16:
				assert(N == 16);
				for (l = 0; l < N; l++) {
					dHdr[0][l] = Hx1_16[l]*(1+Hy1_16[l]*ry)
							*(Hz0_16[l]+Hz2_16[l]*rz*rz)*(1+Hz1_16[l]*rz)/64;
					dHdr[1][l] = (1+Hx1_16[l]*rx)*Hy1_16[l]
							*(Hz0_16[l]+Hz2_16[l]*rz*rz)*(1+Hz1_16[l]*rz)/64;
					dHdr[2][l] = (1+Hx1_16[l]*rx)*(1+Hy1_16[l]*ry)
							*((2*Hz2_16[l]*rz)*(1+Hz1_16[l]*rz)
						+ (Hz0_16[l]+Hz2_16[l]*rz*rz)*Hz1_16[l])/64;
					//H[l][0] = (1+Hx1_16[l]*rx)*(1+Hy1_16[l]*ry)
					//		*(Hz0_16[l]+Hz2_16[l]*rz*rz)*(1+Hz1_16[l]*rz)/64;
				}
				break;
			}
			tmatrix<double,NDIM,NDIM> J(dHdr*xe);
			// This returns det(J);
			gw *= inv_mul_gauss(J,dHdr);
			for (i = 0; i < N; i++) {
				for (j = i; j < N; j++) {
					for (k = 0; k < 3; k++)
						for (l = 0; l < 3; l++) 
							K(i,j) += E[k][l]*(dHdr[k][i]*dHdr[l][j]*gw);
				}
			}
		}
		return _K;
	}

/*
	% shape functions for 16-node 3D brick:
	% the ordering depends on the find above!!!!!!!!!!!
	Nd  = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 2; 1; 2; 4; 2; 1; 2; 1];
	N1  = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
	Nx0 = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0; 0; 0;-1;-1;-1; 0; 0; 0];
	Nx1 = [-1;-1;+1;+1;-1;-1;+1;+1;-1;-1;+1;+1;+1;+1;+1; 0; 0; 0;-1;-1;-1];
	Nxa = [ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1;-1;-1;+1;+1;+1;-1;-1;-1];
	Ny0 = [ 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 0;-1; 0; 0;-1; 0; 0;-1; 0];
	Ny1 = [-1;+1;-1;+1;-1;+1;-1;+1;-1;+1;-1;+1;+1; 0;-1;+1; 0;-1;+1; 0;-1];
	Nya = [ 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0;-1;+1;-1;-1;+1;-1;-1;+1;-1];
	Nz0 = [-1;-1;-1;-1;+9;+9;+9;+9;+9;+9;+9;+9;-1;-1;-1;-1;-1;-1;-1;-1;-1];
	Nz2 = [+9;+9;+9;+9;-9;-9;-9;-9;-9;-9;-9;-9;+9;+9;+9;+9;+9;+9;+9;+9;+9];
	Nz1 = [-1;-1;-1;-1;-3;-3;-3;-3;+3;+3;+3;+3;+1;+1;+1;+1;+1;+1;+1;+1;+1];

				dNdrx = (Nx1+Nxa*sign(rx)) ...
					  .*(Ny0+Ny1*ry+Nya*abs(ry)) ...
					  .*(Nz0+Nz2*rz^2).*(N1+Nz1*rz).*Nd/64;
				dNdry = (Nx0+Nx1*rx+Nxa*abs(rx)) ...
					  .*(Ny1+Nya*sign(ry)) ...
					  .*(Nz0+Nz2*rz^2).*(N1+Nz1*rz).*Nd/64;
				dNdrz = (Nx0+Nx1*rx+Nxa*abs(rx)) ...
					  .*(Ny0+Ny1*ry+Nya*abs(ry)) ...
					  .*((2*Nz2*rz).*(N1+Nz1*rz) ...
						 +(Nz0+Nz2*rz^2).*Nz1).*Nd/64;
				N = (Nx0+Nx1*rx+Nxa*abs(rx)) ...
				  .*(Ny0+Ny1*ry+Nya*abs(ry)) ...
				  .*(Nz0+Nz2*rz^2).*(N1+Nz1*rz).*Nd/64;
*/

protected:
	friend class mesh;
	fset<int> inel;
};

class smatrix_diag;
class smatrix;

/*
 * class smatrix_node - a node in a stifness matrix
 */
class smatrix_node {
	explicit smatrix_node(const int I, const int J,
			const tmatrix<double,NDOF,NDOF> & t, smatrix_diag * d)
	  : K(t), i(I), j(J), col_next(0), col_prev(0), col_diag(d),
	  		row_next(0) {
	}
	void *operator new(size_t, void * p) {
		return p;
	} 

	friend class smatrix_diag;
	friend class smatrix;
	friend class mesh;
	
	tmatrix<double,NDOF,NDOF> K;
	int i, j;
	smatrix_node * col_next;
	smatrix_node * col_prev;
	smatrix_diag * col_diag;
	smatrix_node * row_next;
}; 

/*
 * class smatrix_diag - a diagonal in a stiffness matrix
 */
class smatrix_diag {
	// NO CONSTRUCTOR - Make sure you zero the memory!!!
	inline ~smatrix_diag() {
		if (nodes)
			// smatrix_nodes dont have a destructor...
			free(nodes);
	}

	smatrix_node * insert(int i, int j, const tmatrix<double,NDOF,NDOF> & t,
			smatrix_diag * d) {
		int s = nnz+1, nnb = 8;
		while (s > 8*nnb)
			nnb *= 2;
		s = nnb*(s/nnb+(s%nnb?1:0));
		if (s > nnd) {
			nnd = s;
			smatrix_node * temp = static_cast<smatrix_node *>
					(realloc(nodes,nnd*sizeof(smatrix_node)));
			if (temp == 0) {
				event_msg(EVENT_ERROR,"Out of memory in smatrix_diag::insert()!");
				return 0;
			}
			if (nodes != temp) {
				// Insertion sort the row
				nodes = temp;
				int k, l;
				for (k = 1; k < nnz; k++) {
					for (l = k; l > 0
							&& nodes[l-1].j > nodes[l].j; l--) {
						memcpy(&nodes[nnz],&nodes[l-1],sizeof(smatrix_node));
						memcpy(&nodes[l-1],&nodes[l]  ,sizeof(smatrix_node));
						memcpy(&nodes[l]  ,&nodes[nnz],sizeof(smatrix_node));
					}
				}
				if (row_head)
					row_head = nodes;
				for (k = 0; k < nnz; k++) {
					nodes[k].row_next = (k < nnz-1 ? &nodes[k+1] : 0);
					if (nodes[k].col_prev)
						nodes[k].col_prev->col_next = &nodes[k];
					else
						nodes[k].col_diag->col_head = &nodes[k];
					if (nodes[k].col_next)
						nodes[k].col_next->col_prev = &nodes[k];
				}
			}
		}
		new(&nodes[nnz]) smatrix_node(i,j,t,d);
		// Fix up the row references.
		smatrix_node * p = row_head, * o = 0;
		while (p && p->j < j) {
			o = p; p = p->row_next;
		}
		if (o == 0)
			row_head = &nodes[nnz];
		else
			o->row_next = &nodes[nnz];
		if (p != 0) {
			assert(p->j > j);
			nodes[nnz].row_next = p;
		}
		// Fix up the column refereces.
		p = d->col_head; o = 0;
		while (p && p->i < i) {
			o = p, p = p->col_next;
		}
		if (o == 0)
			d->col_head = &nodes[nnz];
		else {
			o->col_next = &nodes[nnz];
			nodes[nnz].col_prev = o;
		}
		if (p != 0) {
			assert(p->i > i);
			nodes[nnz].col_next = p;
			p->col_prev = &nodes[nnz];
		}
		return &nodes[nnz++];
	}

	friend class smatrix;
	friend class mesh;

	tmatrix<double,NDOF,NDOF> K;
	int nnz, nnd;
	smatrix_node * col_head;
	smatrix_node * row_head;
	smatrix_node * nodes;
}; 

/*
 * class smatrix - a symmetric positive definite stiffness matrix
 */
class smatrix {
public:
	inline explicit smatrix(const int n)
	  : nnd(n), diag(0) {
		diag = static_cast<smatrix_diag *>
				(calloc(nnd,sizeof(smatrix_diag)));
		if (diag == 0)
			event_msg(EVENT_ERROR,"Out of memory in smatrix::smatrix()!");
	}
	inline ~smatrix() {
		if (diag) {
			for (int i = 0; i < nnd; i++)
				diag[i].~smatrix_diag();
			free(diag);
		}
	}
	
	bool append(int i, int j, const tmatrix<double,NDOF,NDOF> & t) {
		if (i == j) {
			diag[i].K += t;
			return true;
		}
		tmatrix<double,NDOF,NDOF> n;
		if (j < i) {
			swap(i,j);
			n = ~t;
		} else {
			n = t;
		}
		smatrix_node * p = diag[i].row_head;
		while (p && p->j < j)
			p = p->row_next;
		if (p == 0 || p->j > j) {
			if (!diag[i].insert(i,j,n,&diag[j]))
				return false;
		} else
			p->K += n;
		return true;
	}

private:
	friend class mesh;
	int nnd;
	smatrix_diag * diag;
};

/*
 * class svector - a vector to go with smatrix
 */
class svector {
public:
	inline explicit svector(const int n)
	  : nnd(n), V(0) {
		V = static_cast<tmatrix<double,NDOF,1> *>
				(calloc(nnd,sizeof(tmatrix<double,NDOF,1>)));
		if (V == 0)
			event_msg(EVENT_ERROR,"Out of memory in svector::svector()!");
	}
	inline ~svector() {
		if (V)
			free(V);
	}
	inline const tmatrix<double,NDOF,1> & operator() (const int i) const {
		return V[i];
	}
	inline tmatrix<double,NDOF,1> & operator() (const int i) {
		return V[i];
	}

private:
	int nnd;
	tmatrix<double,NDOF,1> * V;
};

/*
 * struct mesh_bc - a holder for mesh boundary conditions and forces.
 */
struct mesh_bc_key {
	int n, i;
	mesh_bc_key()
	  : n(0), i(0) {
	}
	mesh_bc_key(int N, int I)
	  : n(N), i(I) {
	}
	bool operator== (const mesh_bc_key & k) const {
		return (n == k.n && i == k.i);
	}  
	bool operator> (const mesh_bc_key & k) const {
		return (n > k.n || (n == k.n && i > k.i));
	}  
};
struct mesh_bc : public mesh_bc_key {
	double d;
	mesh_bc()
	  : mesh_bc_key(), d(0.0) {
	}
	mesh_bc(int N, int I, double D)
	  : mesh_bc_key(N,I), d(D) {
	}
};

/*
 * class mesh - a 3D FEM mesh
 */
class mesh : private list_owned<mesh,element> {
public:
	explicit mesh(const int b)
	  : list_owned<mesh,element>(), node(0,b), disp_bc(), f_ext() {
	}

	bool add(const element::element_t t, const material & m,
			const fset<node3d> & c) {
		assert(8 == c.length());
		element * e = new element_N<8>(this,last,t,m,c);
		if (e == 0) {
			event_msg(EVENT_ERROR,"Out of memory in mesh::add()!");
			return false;
		}
		return true;
	}
	bool add_bc(const node3d & p, const int i, const double d) {
		int k = node.haskey(p);
		assert(i >= 0 && i < NDOF && k >= 0);
		if (k == -1)
			return false;
		disp_bc.add(mesh_bc(k,i,d));
		return true;
	}
	bool add_fext(const node3d & p, const int i, const double d) {
		int k = node.haskey(p);
		assert(i >= 0 && i < NDOF && k >= 0);
		if (k == -1)
			return false;
		f_ext.add(mesh_bc(k,i,d));
		return true;
	}
	bool solve() {
		int i, j, nnd = node.length();
		printf("Solving with %i nodes!\n",nnd);
		smatrix K(nnd);
		svector F(nnd);
		svector U(nnd);
		svector P(nnd);
		svector W(nnd);
		element * e = this->first;
		smatrix_diag * d;
		smatrix_node * p;
		while (e) {
			smatrix_elem * ke = e->stiffness();
			if (ke == 0)
				return false;
			for (i = 0; i < ke->nnd; i++) {
				for (j = i; j < ke->nnd; j++) {
					K.append(ke->inel[i],ke->inel[j],(*ke)(i,j));
					//printf("Element stiffness (%i,%i):\n",i,j);
					//(*ke)(i,j).print();
				}
			}
			delete ke;
			e = e->next;
		}
		f_ext.sort();
		for (i = 0; i < f_ext.length(); i++) {
			const mesh_bc & f = f_ext[i];
			//printf("Adding force %f at node %i dof %i\n",f.d,f.n,f.i);
			F(f.n)(f.i) = f.d;
		}
		disp_bc.sort();
		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			// XXX: extract forces
			//printf("Adding bc %f at node %i dof %i\n",u.d,u.n,u.i);
			assert(u.n >= 0 && u.n < nnd && u.i >= 0 && u.i < NDOF);
			d = &(K.diag[u.n]);
			for (j = 0; j < NDOF; j++) {
				d->K(u.i,j) = 0.0; d->K(j,u.i) = 0.0;
			}
			d->K(u.i,u.i) = 1.0;
			p = d->col_head;
			while (p) {
				for (j = 0; j < NDOF; j++)
					p->K(j,u.i) = 0.0;
				p = p->col_next;
			}
			p = d->row_head;
			while (p) {
				for (j = 0; j < NDOF; j++)
					p->K(u.i,j) = 0.0;
				p = p->row_next;
			}
		}

		/*for (i = 0; i < nnd; i++) {
			d = &(K.diag[i]);
			p = d->col_head;
			while (p) {
				assert(p->j == i);
				printf("New stiffness %i,%i:\n",p->j,p->i);
				tmatrix<double,NDOF,NDOF> tc(~(p->K));
				tc.print();
				p = p->col_next;
			}
			printf("New stiffness %i,%i:\n",i,i);
			tmatrix<double,NDOF,NDOF> t(d->K);
			t.print();
			p = d->row_head;
			while (p) {
				assert(p->i == i);
				printf("New stiffness %i,%i:\n",p->i,p->j);
				tmatrix<double,NDOF,NDOF> tr(p->K);
				tr.print();
				p = p->row_next;
			}
		}*/

		// CG solution
		int it = 0;
		double r = 0.0, ro = 0.0;
		for (i = 0; i < nnd; i++)
			r += tmatrix_scalar<double>(~F(i)*F(i));
		double ri = r;
		while (r > 1e-30*ri && it < nnd*NDOF) {
			if (it == 0) {
				for (i = 0; i < nnd; i++)
					P(i) = F(i);
			} else {
				for (i = 0; i < nnd; i++)
					P(i) = F(i) + r/ro*P(i);
			}
			double a = 0.0;
			for (i = 0; i < nnd; i++) {
				d = &(K.diag[i]);
				tmatrix<double,NDOF,1> t(d->K*P(i));
				p = d->col_head;
				while (p) {
					t += ~(p->K)*P(p->i);
					p = p->col_next;
				}
				p = d->row_head;
				while (p) {
					t += (p->K)*P(p->j);
					p = p->row_next;
				}
				W(i) = t;
				a += tmatrix_scalar<double>(~P(i)*t);
			}
			if (a <= 0) {
				event_msg(EVENT_ERROR,"negative curvature!");
				return false;
			}
			ro = r; r = 0.0;
			for (i = 0; i < nnd; i++) {
				U(i) += ro/a*P(i);
				F(i) -= ro/a*W(i);
				r += tmatrix_scalar<double>(~F(i)*F(i));
			}
			it++;
			printf("CG step %i with residual %g\n",it,r);
		}
		printf("CG took %i steps with residual %g\n",it,r);
		/*for (i = 0; i < nnd; i++) {
			double x = node[i].x;
			double y = node[i].y;
			double z = node[i].z;
			double ux = U(i)(0);
			double uy = U(i)(1);
			double uz = U(i)(2);
			printf("Node %i: (%f,%f,%f) = (%f,%f,%f)\n",i,x,y,z,ux,uy,uz);
		}*/
		return true;
	}

private:
	friend class listelement_o<mesh,element>;
	friend class element;
	kiset<node3d,node3d> node;
	koset<mesh_bc_key,mesh_bc> disp_bc;
	koset<mesh_bc_key,mesh_bc> f_ext;
};

/*
 * Add a node to the mesh's node list, returning the index.
 */
int
element::addnode(const node3d & c)
{
	int k = owner->node.haskey(c);
	if (k == -1) {
		k = owner->node.length();
		owner->node.add(c);
	}
	return k;
}

/*
 * Get a node from the node list, based on the index.
 */
inline const node3d &
element::getnode(const int i)
{
	return  owner->node[i];
}

/*
 * This program is a custom 3D finite element code, intended for
 * work on my PhD.  It will hopefully form the basis for a later
 * full 3D FEM pavement modelling code, but at the moment, I need
 * to get some work done... 
 */
int
main()
{
	material m;
	m.setprop(material_property::emod,1000);
	m.setprop(material_property::poissons,0.2);

	const double domain[3][2] = {{-10, 10}, {-10, 10}, {-10, 0}};
	const int ndiv[3] = {40, 40, 20};
	int i, j, k;
	mesh FEM((ndiv[0]+1)*(ndiv[1]+1)*(ndiv[2]+1));
	
	double dx = (domain[0][1]-domain[0][0])/ndiv[0];
	double dy = (domain[1][1]-domain[1][0])/ndiv[1];
	double dz = (domain[2][1]-domain[2][0])/ndiv[2];
	fset<node3d> coord(8);

	for (i = 0; i < ndiv[0]; i++) {
		double x = domain[0][0] + i*dx;
		for (j = 0; j < ndiv[1]; j++) {
			double y = domain[1][0] + j*dy;
			for (k = 0; k < ndiv[2]; k++) {
				double z = domain[2][0] + k*dz;
				coord[0] = node3d(x   ,y   ,z   );
				coord[1] = node3d(x   ,y+dy,z   );
				coord[2] = node3d(x+dx,y   ,z   );
				coord[3] = node3d(x+dx,y+dy,z   );
				coord[4] = node3d(x   ,y   ,z+dz);
				coord[5] = node3d(x   ,y+dy,z+dz);
				coord[6] = node3d(x+dx,y   ,z+dz);
				coord[7] = node3d(x+dx,y+dy,z+dz);
				//printf("%4.2f\t%4.2f\t%4.2f\n",x,y,z);
				FEM.add(element::block8,m,coord);
			}
		}
	}
	for (i = 0; i <= ndiv[0]; i++) {
		for (j = 0; j <= ndiv[1]; j++) {
			double x = domain[0][0] + i*dx;
			double y = domain[1][0] + j*dy;
			FEM.add_bc(node3d(x,y,domain[2][0]),2,0.0);
		}
	}
	//FEM.add_bc(node3d(0.0,0.0,domain[2][0]),0,0.0);
	//FEM.add_bc(node3d(0.0,0.0,domain[2][0]),1,0.0);
	//FEM.add_bc(node3d(domain[0][0],0.0,domain[2][0]),1,0.0);
	FEM.add_bc(node3d(domain[0][0],domain[1][0],domain[2][0]),0,0.0);
	FEM.add_bc(node3d(domain[0][0],domain[1][0],domain[2][0]),1,0.0);
	FEM.add_bc(node3d(domain[0][1],domain[1][0],domain[2][0]),1,0.0);
	double F = -0.15*dx*dy;
	for (i = 0; i <= ndiv[0]; i++) {
		for (j = 0; j <= ndiv[1]; j++) {
			double x = domain[0][0] + i*dx;
			double y = domain[1][0] + j*dy;
			double f = F;
			if (x == domain[0][0] || x == domain[0][1])
				f /= 2; 
			if (y == domain[1][0] || y == domain[1][1])
				f /= 2; 
			FEM.add_fext(node3d(x,y,domain[2][1]),2,f);
		}
	}
	FEM.solve();

	return 0;
}
