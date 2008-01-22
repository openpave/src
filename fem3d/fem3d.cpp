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

	History:
		2007/09/07 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#define _EVENT_IMP
#define _PROGRESS_IMP
#if !defined(_MSC_VER) && !defined(DARWIN)
#include <time.h>
#endif
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
struct coord2d {
	fixed<7> x, y;

	coord2d() {
	}
	coord2d(double px, double py)
	  : x(px), y(py) {
	}
	coord2d(const coord2d & p)
	  : x(p.x), y(p.y) {
	}
	coord2d(const point2d & p)
	  : x(p.x), y(p.y) {
	}
	~coord2d () {
	}
	inline double distance(const coord2d & p) {
		return hypot(double(x-p.x),double(y-p.y));
	}
	int compare(const coord2d & p) const {
		if (x == p.x && y == p.y)
			return 0;
		if (x < p.x)
			return -1;
		if (x > p.x)
			return 1;
		return (y < p.y ? -1 : 1);
	}
	inline bool operator == (const coord2d & p) const {
		return (x == p.x && y == p.y ? true : false);
	}
	inline bool operator != (const coord2d & p) const {
		return (x != p.x || y != p.y ? true : false);
	}
	inline bool operator > (const coord2d & p) const {
		return (compare(p) == 1 ? true : false);
	}
	inline bool operator >= (const coord2d & p) const {
		return (compare(p) != -1 ? true : false);
	}
	inline bool operator < (const coord2d & p) const {
		return (compare(p) == -1 ? true : false);
	}
	inline bool operator <= (const coord2d & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * A simple point in 3D space.
 */
struct coord3d : public coord2d {
	fixed<7> z;

	coord3d() : coord2d() {
	}
	coord3d(double px, double py, double pz)
	  : coord2d(px,py), z(pz) {
	}
	coord3d(const coord3d & p)
	  : coord2d(p), z(p.z) {
	}
	coord3d(const point3d & p)
	  : coord2d(p), z(p.z) {
	}
	~coord3d () {
	}
	int compare(const coord3d & p) const {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if (z < p.z)
			return -1;
		if (z > p.z)
			return 1;
		return coord2d::compare(static_cast<const coord2d &>(p));
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
 * struct node3d - a Node in 3D
 */
struct node3d : public coord3d {
	int xm, xp, ym, yp, zm, zp;
	double ux, uy, uz;
	
	node3d()
	  : coord3d(), xm(-1), xp(-1), ym(-1), yp(-1), zm(-1), zp(-1),
			ux(0.0), uy(0.0), uz(0.0) {
	}
	node3d(const coord3d & c)
	  : coord3d(c), xm(-1), xp(-1), ym(-1), yp(-1), zm(-1), zp(-1),
			ux(0.0), uy(0.0), uz(0.0) {
	}
	void setneighbours(int x_m, int x_p, int y_m, int y_p, int z_m, int z_p) {
		if (x_m != -1) {
			assert(xm == -1 || xm == x_m);
			xm = x_m;
		}
		if (x_p != -1) {
			assert(xp == -1 || xp == x_p);
			xp = x_p;
		}
		if (y_m != -1) {
			assert(ym == -1 || ym == y_m);
			ym = y_m;
		}
		if (y_p != -1) {
			assert(yp == -1 || yp == y_p);
			yp = y_p;
		}
		if (z_m != -1) {
			assert(zm == -1 || zm == z_m);
			zm = z_m;
		}
		if (z_p != -1) {
			assert(zp == -1 || zp == z_p);
			zp = z_p;
		}
	}
	void setdisp(const tmatrix<double,NDIM,1> & t) {
		ux = t(0);
		uy = t(1);
		uz = t(2);
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
		block16,
		block34
	} etype;

	element(mesh * o, element * p, const element_t t,
			const material & m)
	  : listelement_o<mesh,element>(o,p), etype(t), mat(m), inel(0,8) {
	}
	virtual smatrix_elem * stiffness() const = 0;
	
protected:
	friend class mesh;
	const material & mat;
	sset<int> inel;

	inline int addnode(const coord3d & c) const;
	inline void updatenode(const node3d & n) const;
	inline const node3d & getnode(const int i) const;
};

/*
 * class element_block8 - a finite element
 */
class element_block8 : public element {
public:
	element_block8(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,block8,m) {
	  	assert(8 == c.length());
		for (int i = 0; i < 8; i++)
			inel.add(addnode(c[i]));
		for (int i = 0; i < 8; i++) {
			node3d n = getnode(inel[i]);
			int x_m = -1, x_p = -1;
			int y_m = -1, y_p = -1;
			int z_m = -1, z_p = -1;
			(i%4 < 2 ? x_p : x_m) = inel[4*(i/4)+(i%4+2)%4];
			(i%2 < 1 ? y_p : y_m) = inel[2*(i/2)+(i%2+1)%2];
			(i   < 4 ? z_p : z_m) = inel[(i+4)%8];
			n.setneighbours(x_m,x_p,y_m,y_p,z_m,z_p);
			updatenode(n);
		}
	}
	virtual smatrix_elem * stiffness() const {
		assert(8 == inel.length());
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
		tmatrix<double,8,NDIM> xe;
		for (i = 0; i < 8; i++) {
			const coord3d & p = getnode(inel[i]);
			xe(i,0) = p.x; xe(i,1) = p.y; xe(i,2) = p.z;
		}

		ksset<point3d,gauss3d> gp(0,8);
		const double gp_2[2][2] = {{-1.0/sqrt(3.0), 1.0},
		                           {+1.0/sqrt(3.0), 1.0}};
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				for (k = 0; k < 2; k++) {
					gp.add(gauss3d(gp_2[i][0],gp_2[j][0],gp_2[k][0],
							gp_2[i][1]*gp_2[j][1]*gp_2[k][1]));
				}
			}
		}
	
		smatrix_elem * _K = new smatrix_elem(8,inel);
		if (_K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}
		smatrix_elem & K = *_K;
	
		for (g = 0; g < gp.length(); g++) {
			rx = gp[g].x; ry = gp[g].y; rz = gp[g].z;
			double gw = gp[g].gw;
			tmatrix<double,NDIM,8> dHdr;
			//tmatrix<double,8,1> H;
			// shape functions for 8-node 3D brick:
			//   H = 1/8*(1+-x)*(1+-y)*(1+-z);
			const double Hx[8] = {-1, -1, +1, +1, -1, -1, +1, +1};
			const double Hy[8] = {-1, +1, -1, +1, -1, +1, -1, +1};
			const double Hz[8] = {-1, -1, -1, -1, +1, +1, +1, +1};
			for (l = 0; l < 8; l++) {
				dHdr(0,l) = Hx[l]*(1+Hy[l]*ry)*(1+Hz[l]*rz)/8;
				dHdr(1,l) = Hy[l]*(1+Hx[l]*rx)*(1+Hz[l]*rz)/8;
				dHdr(2,l) = Hz[l]*(1+Hx[l]*rx)*(1+Hy[l]*ry)/8;
				//H(0,l) =(1+Hx[l]*rx)*(1+Hy[l]*ry)*(1+Hz[l]*rz)/8;
			}
			tmatrix<double,NDIM,NDIM> J(dHdr*xe);
			// This returns det(J);
			gw *= inv_mul_gauss(J,dHdr);
			for (i = 0; i < 8; i++) {
				for (j = i; j < 8; j++) {
					for (k = 0; k < 3; k++)
						for (l = 0; l < 3; l++) 
							K(i,j) += E[k][l]*(dHdr(k,i)*dHdr(l,j)*gw);
				}
			}
		}
		return _K;
	}
};

/*
 * class element_block16 - a finite element
 */
class element_block16 : public element {
public:
	element_block16(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,block16,m) {
	  	assert(8 == c.length());
		double xb[4], yb[4], zb[4];
		double xt[4], yt[4], zt[4];
		for (int i = 0; i < 4; i++) {
			xb[i] = c[i  ].x; yb[i] = c[i  ].y; zb[i] = c[i  ].z;
			xt[i] = c[i+4].x; yt[i] = c[i+4].y; zt[i] = c[i+4].z;
		}
		for (int i = 0; i < 4; i++)
			inel.add(addnode(c[i]));
		for (int i = 0; i < 4; i++) {
			double x = xb[i]+  (xt[i]-xb[i])/3;
			double y = yb[i]+  (yt[i]-yb[i])/3;
			double z = zb[i]+  (zt[i]-zb[i])/3;
			inel.add(addnode(coord3d(x,y,z)));
		}
		for (int i = 0; i < 4; i++) {	
			double x = xb[i]+2*(xt[i]-xb[i])/3;
			double y = yb[i]+2*(yt[i]-yb[i])/3;
			double z = zb[i]+2*(zt[i]-zb[i])/3;
			inel.add(addnode(coord3d(x,y,z)));
		}
		for (int i = 0; i < 4; i++)
			inel.add(addnode(c[i+4]));
		for (int i = 0; i < 16; i++) {
			node3d n = getnode(inel[i]);
			int x_m = -1, x_p = -1;
			int y_m = -1, y_p = -1;
			int z_m = -1, z_p = -1;
			(i%4 < 2 ? x_p : x_m) = inel[4*(i/4)+(i%4+2)%4];
			(i%2 < 1 ? y_p : y_m) = inel[2*(i/2)+(i%2+1)%2];
			if (i < 12)
				z_p = inel[i+4];
			if (i >= 4)
				z_m = inel[i-4];
			n.setneighbours(x_m,x_p,y_m,y_p,z_m,z_p);
			updatenode(n);
		}
	}
	virtual smatrix_elem * stiffness() const {
		assert(16 == inel.length());
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
		tmatrix<double,16,NDIM> xe;
		for (i = 0; i < 16; i++) {
			const coord3d & p = getnode(inel[i]);
			xe(i,0) = p.x; xe(i,1) = p.y; xe(i,2) = p.z;
		}

		ksset<point3d,gauss3d> gp(0,8);
		const double gp_2[2][2] = {{-1.0/sqrt(3.0), 1.0},
		                           {+1.0/sqrt(3.0), 1.0}};
		//const double gp_3[3][2] = {{-sqrt(3.0/5.0), 5.0/9.0},
		//                           {             0, 8.0/9.0},
		//                           {+sqrt(3.0/5.0), 5.0/9.0}};
		const double gp_4[4][2] =
			{{-sqrt(525.0+70.0*sqrt(30.0))/35.0, (18.0-sqrt(30.0))/36.0},
			 {-sqrt(525.0-70.0*sqrt(30.0))/35.0, (18.0+sqrt(30.0))/36.0},
			 {+sqrt(525.0-70.0*sqrt(30.0))/35.0, (18.0+sqrt(30.0))/36.0},
			 {+sqrt(525.0+70.0*sqrt(30.0))/35.0, (18.0-sqrt(30.0))/36.0}};
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				for (k = 0; k < 4; k++) {
					gp.add(gauss3d(gp_2[i][0],gp_2[j][0],gp_4[k][0],
							gp_2[i][1]*gp_2[j][1]*gp_4[k][1]));
				}
			}
		}
	
		smatrix_elem * _K = new smatrix_elem(16,inel);
		if (_K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}
		smatrix_elem & K = *_K;
	
		for (g = 0; g < gp.length(); g++) {
			rx = gp[g].x; ry = gp[g].y; rz = gp[g].z;
			double gw = gp[g].gw;
			tmatrix<double,NDIM,16> dHdr;
			//tmatrix<double,16,1> H;
	
			// shape functions for 16-node 3D brick:
			const double Hx1[16] = {-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1};
			const double Hy1[16] = {-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1};
			const double Hz0[16] = {-1,-1,-1,-1,+9,+9,+9,+9,+9,+9,+9,+9,-1,-1,-1,-1};
			const double Hz2[16] = {+9,+9,+9,+9,-9,-9,-9,-9,-9,-9,-9,-9,+9,+9,+9,+9};
			const double Hz1[16] = {-1,-1,-1,-1,-3,-3,-3,-3,+3,+3,+3,+3,+1,+1,+1,+1};
	
			for (l = 0; l < 16; l++) {
				dHdr(0,l) = Hx1[l]*(1+Hy1[l]*ry)
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/64;
				dHdr(1,l) = (1+Hx1[l]*rx)*Hy1[l]
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/64;
				dHdr(2,l) = (1+Hx1[l]*rx)*(1+Hy1[l]*ry)
						*((2*Hz2[l]*rz)*(1+Hz1[l]*rz)
					+ (Hz0[l]+Hz2[l]*rz*rz)*Hz1[l])/64;
				//H(0,l) = (1+Hx1[l]*rx)*(1+Hy1[l]*ry)
				//		*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/64;
			}
			tmatrix<double,NDIM,NDIM> J(dHdr*xe);
			// This returns det(J);
			gw *= inv_mul_gauss(J,dHdr);
			for (i = 0; i < 16; i++) {
				for (j = i; j < 16; j++) {
					for (k = 0; k < 3; k++)
						for (l = 0; l < 3; l++) 
							K(i,j) += E[k][l]*(dHdr(k,i)*dHdr(l,j)*gw);
				}
			}
		}
		return _K;
	}
};

/*
 * class element_block34 - a finite element
 */
class element_block34 : public element {
public:
	element_block34(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,block34,m) {
	  	assert(8 == c.length());
		int i;
		double xb[4], yb[4], zb[4];
		double xt[4], yt[4], zt[4];
		for (i = 0; i < 4; i++) {
			xb[i] = c[i  ].x; yb[i] = c[i  ].y; zb[i] = c[i  ].z;
			xt[i] = c[i+4].x; yt[i] = c[i+4].y; zt[i] = c[i+4].z;
		}
		for (i = 0; i < 4; i++)
			inel.add(addnode(c[i]));
		for (i = 0; i < 4; i++) {
			double x = xb[i]+  (xt[i]-xb[i])/3;
			double y = yb[i]+  (yt[i]-yb[i])/3;
			double z = zb[i]+  (zt[i]-zb[i])/3;
			inel.add(addnode(coord3d(x,y,z)));
		}
		for (i = 0; i < 4; i++) {	
			double x = xb[i]+2*(xt[i]-xb[i])/3;
			double y = yb[i]+2*(yt[i]-yb[i])/3;
			double z = zb[i]+2*(zt[i]-zb[i])/3;
			inel.add(addnode(coord3d(x,y,z)));
		}
		for (i = 0; i < 4; i++)
			inel.add(addnode(c[i+4]));

		for (i = 0; i < 18; i++)
			mask[i] = -1;
		for (i = 0; i < 16; i++) {
			node3d n = getnode(inel[i]);
			int x_m = -1, x_p = -1;
			int y_m = -1, y_p = -1;
			int z_m = -1, z_p = -1;
			(i%4 < 2 ? x_p : x_m) = inel[4*(i/4)+(i%4+2)%4];
			(i%2 < 1 ? y_p : y_m) = inel[2*(i/2)+(i%2+1)%2];
			if (i < 12)
				z_p = inel[i+4];
			if (i >= 4)
				z_m = inel[i-4];
			switch (i%4) {
			case 0:
				if (n.xp != -1 && n.xp != x_p) {
					x_p = n.xp;
					mask[i+2] = x_p;
				}
				if (n.yp != -1 && n.yp != y_p) {
					y_p = n.yp;
					mask[i] = y_p;
					if ((i == 0 || i == 12)) {
						const node3d & mid = getnode(n.yp); 
						if (mid.xp != -1)
							mask[i == 0 ? 16 : 17] = mid.xp;
					}
				}
				break;
			case 1:
				if (n.xp != -1 && n.xp != x_p)
					x_p = n.xp;
				if (n.ym != -1 && n.ym != y_m)
					y_m = n.ym;
				break;
			case 2:
				if (n.xm != -1 && n.xm != x_m)
					x_m = n.xm;
				if (n.yp != -1 && n.yp != y_p)
					y_p = n.yp;
				break;
			case 3:
				if (n.xm != -1 && n.xm != x_m) {
					x_m = n.xm;
					mask[i-2] = x_m;
				}
				if (n.ym != -1 && n.ym != y_m) {
					y_m = n.ym;
					mask[i] = y_m;
				}
				break;
			}
			n.setneighbours(x_m,x_p,y_m,y_p,z_m,z_p);
			updatenode(n);
		}
		for (i = 0; i < 18; i++) {
			if (mask[i] != -1) {
				inel.add(mask[i]);
				mask[i] = inel.length()-1;
			} else
				mask[i] = 0;
		}
	}
	virtual smatrix_elem * stiffness() const {
		int i, j, k, l, g, nnd = inel.length();
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
		matrix_dense xe(nnd,NDIM);
		for (i = 0; i < nnd; i++) {
			const coord3d & p = getnode(inel[i]);
			xe(i,0) = p.x; xe(i,1) = p.y; xe(i,2) = p.z;
		}

		ksset<point3d,gauss3d> gp(0,8);
		const double gp_A[4][2] = {{-(1.0+1.0/sqrt(3.0))/2.0, 1.0/2.0},
		                           {-(1.0-1.0/sqrt(3.0))/2.0, 1.0/2.0},
		                           {+(1.0-1.0/sqrt(3.0))/2.0, 1.0/2.0},
		                           {+(1.0+1.0/sqrt(3.0))/2.0, 1.0/2.0}};
		//const double gp_3[3][2] = {{-sqrt(3.0/5.0), 5.0/9.0},
		//                           {             0, 8.0/9.0},
		//                           {+sqrt(3.0/5.0), 5.0/9.0}};
		const double gp_4[4][2] =
			{{-sqrt(525.0+70.0*sqrt(30.0))/35.0, (18.0-sqrt(30.0))/36.0},
			 {-sqrt(525.0-70.0*sqrt(30.0))/35.0, (18.0+sqrt(30.0))/36.0},
			 {+sqrt(525.0-70.0*sqrt(30.0))/35.0, (18.0+sqrt(30.0))/36.0},
			 {+sqrt(525.0+70.0*sqrt(30.0))/35.0, (18.0-sqrt(30.0))/36.0}};
		for (i = 0; i < 4; i++) {
			for (j = 0; j < 4; j++) {
				for (k = 0; k < 4; k++) {
					gp.add(gauss3d(gp_A[i][0],gp_A[j][0],gp_4[k][0],
							gp_A[i][1]*gp_A[j][1]*gp_4[k][1]));
				}
			}
		}
	
		smatrix_elem * _K = new smatrix_elem(nnd,inel);
		if (_K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}
		smatrix_elem & K = *_K;

		for (g = 0; g < gp.length(); g++) {
			rx = gp[g].x; ry = gp[g].y; rz = gp[g].z;
			double gw = gp[g].gw;
			matrix_dense dHdr(NDIM,nnd);
			//matrix_dense H(nnd,1);
	
			// shape functions for 16/34-node 3D brick:
			const double Sxy[16] = { 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0};
			const double Hx1[16] = {-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,+1};
			const double Hy1[16] = {-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1};
			const double Hz0[16] = {-1,-1,-1,-1,+9,+9,+9,+9,+9,+9,+9,+9,-1,-1,-1,-1};
			const double Hz2[16] = {+9,+9,+9,+9,-9,-9,-9,-9,-9,-9,-9,-9,+9,+9,+9,+9};
			const double Hz1[16] = {-1,-1,-1,-1,-3,-3,-3,-3,+3,+3,+3,+3,+1,+1,+1,+1};
	
			for (l = 0; l < 16; l++) {
				dHdr(0,l) = Hx1[l]*(1+Hy1[l]*ry)
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/64;
				dHdr(1,l) = (1+Hx1[l]*rx)*Hy1[l]
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/64;
				dHdr(2,l) = (1+Hx1[l]*rx)*(1+Hy1[l]*ry)
						*((2*Hz2[l]*rz)*(1+Hz1[l]*rz)
					+ (Hz0[l]+Hz2[l]*rz*rz)*Hz1[l])/64;
				//H(l) = (1+Hx1[l]*rx)*(1+Hy1[l]*ry)
				//		*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/64;
			}
			for (l = 0; l < 16; l++) {
				if ((i = mask[l]) == 0)
					continue;
				dHdr(0,i) = (Sxy[l] ? -SGN(rx) : Hx1[l])
				        *(1+(!Sxy[l] ? -fabs(ry) : Hy1[l]*ry))
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/32;
				dHdr(1,i) = (1+(Sxy[l] ? -fabs(rx) : Hx1[l]*rx))
				        *(!Sxy[l] ? -SGN(ry) : Hy1[l])
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/32;
				dHdr(2,i) = (1+(Sxy[l] ? -fabs(rx) : Hx1[l]*rx))
				        *(1+(!Sxy[l] ? -fabs(ry) : Hy1[l]*ry))
						*((2*Hz2[l]*rz)*(1+Hz1[l]*rz)
					+ (Hz0[l]+Hz2[l]*rz*rz)*Hz1[l])/32;
				//H(i) = (1+(Sxy[l] ? -fabs(rx) : Hx1[l]*rx))
				//      *(1+(!Sxy[l] ? -fabs(ry) : Hy1[l]*ry))
				//		*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/32;
			}
			for (l = 0; l < 16; l += 12) {
				if ((i = mask[(l == 0 ? 16 : 17)]) == 0)
					continue;
				dHdr(0,i) = -SGN(rx)*(1-fabs(ry))
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/16;
				dHdr(1,i) = -(1-fabs(rx))*SGN(ry)
						*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/16;
				dHdr(2,i) = (1-fabs(rx))*(1-fabs(ry))
						*((2*Hz2[l]*rz)*(1+Hz1[l]*rz)
					+ (Hz0[l]+Hz2[l]*rz*rz)*Hz1[l])/16;
				//H(i) = (1-fabs(rx))*(1-fabs(ry))
				//		*(Hz0[l]+Hz2[l]*rz*rz)*(1+Hz1[l]*rz)/16;
			}
			for (l = 0; l < 16; l++) {
				switch (l%4) {
				case 0: j =  2; break;
				case 1: j = -1; break;
				case 2: j =  1; break;
				case 3: j = -2; break;
				}
				for (i = 0; i < 3; i++) {
					if (mask[l])
						dHdr(i,l) -= dHdr(i,mask[l])/2;
					if (mask[l+j])
						dHdr(i,l) -= dHdr(i,mask[l+j])/2;
					if (l < 4 && mask[16])
						dHdr(i,l) += dHdr(i,mask[16])/4;
					if (l >= 12 && mask[17])
						dHdr(i,l) += dHdr(i,mask[17])/4;
				}
				//if (mask[l])
				//	H(l) -= H(mask[l])/2;
				//if (mask[l+j])
				//	H(l) -= H(mask[l+j])/2;
				//if (l < 4 && mask[16])
				//	H(l) += H(mask[16])/4;
				//if (l >= 12 && mask[17])
				//	H(l) += H(mask[17])/4;
			}
			for (l = 0; mask[16] && l < 4; l++) {
				for (i = 0; i < 3; i++)
					dHdr(i,mask[l]) -= dHdr(i,mask[16])/2;
				//H(mask[l]) -= H(mask[16])/2;
			}
			for (l = 12; mask[17] && l < 16; l++) {
				for (i = 0; i < 3; i++)
					dHdr(i,mask[l]) -= dHdr(i,mask[17])/2;
				//H(mask[l]) -= H(mask[17])/2;
			}

			matrix_dense J(dHdr*xe);
			// This returns det(J);
			gw *= inv_mul_gauss(NDIM,nnd,&J(0,0),&dHdr(0,0));
			for (i = 0; i < nnd; i++) {
				for (j = i; j < nnd; j++) {
					for (k = 0; k < 3; k++)
						for (l = 0; l < 3; l++) 
							K(i,j) += E[k][l]*(dHdr(k,i)*dHdr(l,j)*gw);
				}
			}
		}
		return _K;
	}

protected:
	int mask[18];
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
	  		row_prev(0), row_next(0) {
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
	smatrix_node * row_prev;
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
	bool insert(int i, int j, const tmatrix<double,NDOF,NDOF> & t,
			smatrix_diag * d) {
		if (nnz+1 > nnd) {
			nnd = nnz+1;
			smatrix_node * temp =
					static_cast<smatrix_node *>
					(realloc(nodes,nnd*sizeof(smatrix_node)));
			if (temp == 0) {
				event_msg(EVENT_ERROR,"Out of memory in smatrix_diag::insert()!");
				return false;
			}
			if (nodes != temp) {
				nodes = temp;
				for (int k = 1; k < nnz; k++) {
					for (int l = k; l > 0
							&& nodes[l-1].j > nodes[l].j; l--) {
						memcpy(&nodes[nnz],&nodes[l-1],sizeof(smatrix_node));
						memcpy(&nodes[l-1],&nodes[l]  ,sizeof(smatrix_node));
						memcpy(&nodes[l]  ,&nodes[nnz],sizeof(smatrix_node));
					}
				}
				if (row_head)
					row_head = nodes;
				for (int k = 0; k < nnz; k++) {
					nodes[k].row_prev = (k > 0     ? &nodes[k-1] : 0);
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
		else {
			o->row_next = &nodes[nnz];
			nodes[nnz].row_prev = o;
		}
		if (p != 0) {
			assert(p->j > j);
			nodes[nnz].row_next = p;
			p->row_prev = &nodes[nnz];
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
		nnz++;
		return true;
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
	inline explicit smatrix(const smatrix & A)
	  : nnd(A.nnd), diag(0) {
		diag = static_cast<smatrix_diag *>
				(calloc(nnd,sizeof(smatrix_diag)));
		if (diag == 0) {
			event_msg(EVENT_ERROR,"Out of memory in smatrix::smatrix()!");
			return;
		}
		for (int i = 0; i < nnd; i++) {
			smatrix_diag * d = &(A.diag[i]);
			append(i,i,d->K);
			smatrix_node * p = d->row_head;
			while (p) {
				append(p->i,p->j,p->K);
				p = p->row_next;
			}
		}
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

	bool chol() {
		int i, j;
		smatrix_diag * d;
		smatrix_node * p;
		
		for (i = 0; i < nnd; i++) {
			d = &(diag[i]);
			tmatrix<double,NDOF,NDOF> & K = d->K;
			tmatrix<double,NDOF,NDOF> t(0.0);
			p = d->col_head;
			while (p) {
				t += ~(p->K)*(p->K);
				p = p->col_next;
			}
			K(0,0) = sqrt(K(0,0)-t(0,0));
			K(0,1) = (K(0,1)-t(0,1))/K(0,0);
			K(0,2) = (K(0,2)-t(0,2))/K(0,0);
			K(1,1) = sqrt(K(1,1)-K(0,1)*K(0,1)-t(1,1));
			K(1,2) = (K(1,2)-K(0,2)*K(0,1)-t(1,2))/K(1,1);
			K(2,2) = sqrt(K(2,2)-K(0,2)*K(0,2)-K(1,2)*K(1,2)-t(2,2));
			K(1,0) = K(2,0) = K(2,1) = 0.0;
			p = d->row_head;
			while (p) {
				tmatrix<double,NDOF,NDOF> & J = p->K;
				smatrix_node * pi = d->col_head;
				smatrix_node * pj = diag[p->j].col_head;
				while (pi && pj) {
					if (pi->i > pj->i)
						pj = pj->col_next;
					else if (pi->i < pj->i)
						pi = pi->col_next;
					else {
						J -= ~(pi->K)*(pj->K);
						pi = pi->col_next;
						pj = pj->col_next;
					}
				}
				for (j = 0; j < NDOF; j++)
					J(0,j) /= K(0,0);
				for (j = 0; j < NDOF; j++)
					J(1,j) = (J(1,j)-K(0,1)*J(0,j))/K(1,1);
				for (j = 0; j < NDOF; j++)
					J(2,j) = (J(2,j)-K(0,2)*J(0,j)-K(1,2)*J(1,j))/K(2,2);
				p = p->row_next;
			}
		}
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
	explicit mesh()
	  : list_owned<mesh,element>(), node(), disp_bc(), f_ext() {
	}

	bool add(const element::element_t t, const material & m,
			const fset<coord3d> & c) {
		assert(8 == c.length());
		element * e = 0;
		switch (t) {
		case element::block8:
			e = new element_block8(this,last,m,c);
			break;
		case element::block16:
			e = new element_block16(this,last,m,c);
			break;
		case element::block34:
			e = new element_block34(this,last,m,c);
			if (e != 0 && e->inel.length() == 16) {
				delete e;
				e = new element_block16(this,last,m,c);
			}
			break;
		}
		if (e == 0) {
			event_msg(EVENT_ERROR,"Out of memory in mesh::add()!");
			return false;
		}
		return true;
	}
	inline int addnode(const coord3d & p) {
		int k = node.haskey(p);
		if (k == -1) {
			k = node.length();
			node.add(p);
		}
		return k;
	}
	inline void updatenode(const node3d & n) {
		node.replace(n);
	}
	inline int getnodes() const {
		return node.length();
	}
	inline const node3d & getnode(const int i) const {
		return node[i];
	}
	inline const node3d & getorderednode(const int i) const {
		return node.getindex(i);
	}
	inline const int hasnode(const coord3d & p) const {
		return node.haskey(p);
	}
	bool add_bc(const coord3d & p, const int i, const double d) {
		int k = node.haskey(p);
		assert(i >= 0 && i < NDOF && k >= 0);
		if (k == -1)
			return false;
		disp_bc.add(mesh_bc(k,i,d));
		return true;
	}
	bool add_fext(const coord3d & p, const int i, const double d) {
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
		svector F(nnd), U(nnd), P(nnd), W(nnd), Z(nnd);
		element * e = this->first;
		smatrix_diag * d;
		smatrix_node * p;
		while (e) {
			smatrix_elem * ke = e->stiffness();
			if (ke == 0)
				return false;
			for (i = 0; i < ke->nnd; i++) {
				for (j = i; j < ke->nnd; j++) {
					K.append(node.getorder(ke->inel[i]),
					         node.getorder(ke->inel[j]),(*ke)(i,j));
				}
			}
			delete ke;
			e = e->next;
		}
		f_ext.sort();
		for (i = 0; i < f_ext.length(); i++) {
			const mesh_bc & f = f_ext[i];
			F(node.getorder(f.n))(f.i) = f.d;
		}
		disp_bc.sort();
		for (i = 0; i < disp_bc.length(); i++) {
			if (i+2 < disp_bc.length()
			 && disp_bc[i].n == disp_bc[i+1].n
			 && disp_bc[i+1].n == disp_bc[i+2].n) {
				const mesh_bc & u0 = disp_bc[i];
				const mesh_bc & u1 = disp_bc[i+1];
				const mesh_bc & u2 = disp_bc[i+2];
				int n = node.getorder(u0.n);
				F(n)(0) = u0.d; F(n)(1) = u1.d; F(n)(2) = u2.d;
				d = &(K.diag[n]);
				d->K = tmatrix<double,NDOF,NDOF>(1.0,true);
				p = d->col_head;
				while (p) {
					if (u0.d != 0.0 || u1.d != 0.0 || u2.d != 0.0)
						F(p->i) -= p->K*F(n);
					if (p->row_prev) {
						p->row_prev->row_next = p->row_next;
					} else {
						K.diag[p->i].row_head = p->row_next;
					}
					if (p->row_next)
						p->row_next->row_prev = p->row_prev;
					p = p->col_next;
				}
				d->col_head = 0;
				p = d->row_head;
				while (p) {
					if (u0.d != 0.0 || u1.d != 0.0 || u2.d != 0.0)
						F(p->j) -= ~(p->K)*F(n);
					if (p->col_prev)
						p->col_prev->col_next = p->col_next;
					else
						p->col_diag->col_head = p->col_next;
					if (p->col_next)
						p->col_next->col_prev = p->col_prev;
					p = p->row_next;
				}
				d->row_head = 0;
				free(d->nodes);
				d->nodes = 0;
			} else {
				const mesh_bc & u = disp_bc[i];
				int n = node.getorder(u.n);
				d = &(K.diag[n]);
				for (j = 0; j < NDOF; j++) {
					if (j == u.i)
						continue;
					if (u.d != 0.0)
						F(n)(j) -= d->K(j,u.i)*u.d;
					d->K(u.i,j) = 0.0; d->K(j,u.i) = 0.0;
				}
				d->K(u.i,u.i) = 1.0;
				p = d->col_head;
				while (p) {
					for (j = 0; j < NDOF; j++) {
						if (u.d != 0.0)
							F(p->i)(j) -= p->K(j,u.i)*u.d;
						p->K(j,u.i) = 0.0;
					}
					p = p->col_next;
				}
				p = d->row_head;
				while (p) {
					for (j = 0; j < NDOF; j++) {
						if (u.d != 0.0)
							F(p->j)(j) -= p->K(u.i,j)*u.d;
						p->K(u.i,j) = 0.0;
					}
					p = p->row_next;
				}
				F(n)(u.i) = u.d;
			}
		}

		smatrix M(K);
		M.chol();

		// CG solution
		int it = 0;
		double r = 0.0, ro = 0.0;
		for (i = 0; i < nnd; i++)
			r += tmatrix_scalar<double>(~F(i)*F(i));
		double ri = r;
		while (r > 1e-30*ri && it < nnd*NDOF) {
			r = 0.0;
			for (i = 0; i < nnd; i++) {
				d = &(M.diag[i]);
				tmatrix<double,NDOF,1> t(F(i));
				p = d->col_head;
				while (p) {
					t -= ~(p->K)*Z(p->i);
					p = p->col_next;
				}
				t(0) = t(0)/d->K(0,0);
				t(1) = (t(1)-t(0)*d->K(0,1))/d->K(1,1);
				t(2) = (t(2)-t(0)*d->K(0,2)-t(1)*d->K(1,2))/d->K(2,2);
				Z(i) = t;
			}
			for (i = nnd-1; i >= 0; i--) {
				d = &(M.diag[i]);
				tmatrix<double,NDOF,1> t(Z(i));
				p = d->row_head;
				while (p) {
					t -= (p->K)*Z(p->j);
					p = p->row_next;
				}
				t(2) = t(2)/d->K(2,2);
				t(1) = (t(1)-t(2)*d->K(1,2))/d->K(1,1);
				t(0) = (t(0)-t(1)*d->K(0,1)-t(2)*d->K(0,2))/d->K(0,0);
				Z(i) = t;
				r += tmatrix_scalar<double>(~t*F(i));
			}
			if (it == 0) {
				for (i = 0; i < nnd; i++)
					P(i) = Z(i);
			} else {
				for (i = 0; i < nnd; i++)
					P(i) = Z(i) + (r/ro)*P(i);
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
				a += tmatrix_scalar<double>(~t*P(i));
			}
			if (a <= 0) {
				event_msg(EVENT_ERROR,"negative curvature!");
				return false;
			}
			ro = r; r = 0.0;
			for (i = 0; i < nnd; i++) {
				U(i) += (ro/a)*P(i);
				F(i) -= (ro/a)*W(i);
				r += tmatrix_scalar<double>(~F(i)*F(i));
			}
			it++;
			printf("CG step %i with residual %g\n",it,r);
		}
		printf("CG took %i steps with residual %g\n",it,r);
		for (i = 0; i < nnd; i++) {
			node3d n = node.getindex(i);
			n.setdisp(U(i));
			node.replace(n);
		}
		return true;
	}

private:
	friend class listelement_o<mesh,element>;
	friend class element;
	kiset<coord3d,node3d> node;
	koset<mesh_bc_key,mesh_bc> disp_bc;
	koset<mesh_bc_key,mesh_bc> f_ext;
};

/*
 * Add a node to the mesh's node list, returning the index.
 */
int
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
element::getnode(const int i) const
{
	return owner->getnode(i);
}

inline double
circlearea(double a, double b, double r)
{
	a = fabs(a); b = fabs(b);
	if (hypot(a,b) >= r)
		return 0.0;
	double t = M_PI_2 - asin(a/r) - asin(b/r);
	double ar = (sqrt(r*r-b*b)-a)*b/2;
	double br = (sqrt(r*r-a*a)-b)*a/2;
	return (r*r*t/2) - ar -br;
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
 * This program is a custom 3D finite element code, intended for
 * work on my PhD.  It will hopefully form the basis for a later
 * full 3D FEM pavement modelling code, but at the moment, I need
 * to get some work done... 
 */
int
main()
{
#if !defined(_MSC_VER) && !defined(DARWIN)
	// get starting time
	struct timespec start, stop;
	clock_gettime(CLOCK_PROF,&start);
#endif

	material m;
	m.setprop(material_property::emod,100e3); // kPa
	m.setprop(material_property::poissons,0.35);

	double x, y, z;
	double dx, dy, dz, delta = 4.0;
	mesh FEM;
	fset<coord3d> coord(8);
	
	// Start with the tyre grid.
	dx = delta; dy = delta;
	for (x = -30*dx; x < 30*dx; x += 2*dx) {
		for (y = -30*dy; y < 30*dy; y += 2*dy) {
			if (blockarea(x,x+2*dx,y,y+2*dx,100) == 0.0)
				continue;
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
	// Add the tire loads
	double F = -690;
	for (x = -30*dx; x <= 30*dx; x += dx) {
		for (y = -30*dy; y <= 30*dy; y += dy) {
			int p = FEM.hasnode(coord3d(x,y,0.0));
			if (p == -1)
				continue;
			node3d n = FEM.getnode(p);
			int x_m = FEM.hasnode(coord3d(x-dx,y,0.0));
			int x_p = FEM.hasnode(coord3d(x+dx,y,0.0));
			int y_m = FEM.hasnode(coord3d(x,y-dy,0.0));
			int y_p = FEM.hasnode(coord3d(x,y+dy,0.0));
			n.setneighbours(x_m,x_p,y_m,y_p,-1,-1);
			FEM.updatenode(n);
			double f = F*blockarea(x-dx/2,x+dx/2,y-dy/2,y+dy/2,100.0);
			if (fabs(f) > 0.0)
				FEM.add_fext(coord3d(x,y,0.0),2,f);
		}
	}

#define DOMAIN 8192

	// Now add the elements from below the tyre, working outwards.
	double xm = 0.0, xp = 0.0, ym = 0.0, yp = 0.0, zm = 0.0;
	while (zm > -4000.0) {
		dz = MIN(-30.0,zm), zm += dz;
		dx = delta*2; dy = delta*2;
		z = 0.0;
		while (z > zm) {
			dz = MIN(-30.0,z), z += dz;
			for (x = -MIN(16*dx,DOMAIN); x < MIN(16*dx,DOMAIN); x += dx) {
				for (y = -MIN(16*dy,DOMAIN); y < MIN(16*dy,DOMAIN); y += dy) {
					if (xm < xp && (x >= xm && x < xp)
					 && ym < yp && (y >= ym && y < yp)
					 && z > zm)
					 	continue;
					coord[0] = coord3d(x   ,y   ,z);
					coord[1] = coord3d(x   ,y+dy,z);
					coord[2] = coord3d(x+dx,y   ,z);
					coord[3] = coord3d(x+dx,y+dy,z);
					coord[4] = coord3d(x   ,y   ,z-dz);
					coord[5] = coord3d(x   ,y+dy,z-dz);
					coord[6] = coord3d(x+dx,y   ,z-dz);
					coord[7] = coord3d(x+dx,y+dy,z-dz);
					FEM.add(element::block34,m,coord);
					if (x == -DOMAIN) {
						//FEM.add_bc(coord[0],1,0.0);
						//FEM.add_bc(coord[1],1,0.0);
						//FEM.add_bc(coord[4],1,0.0);
						//FEM.add_bc(coord[5],1,0.0);
						FEM.add_bc(coord[0],0,0.0);
						FEM.add_bc(coord[1],0,0.0);
						FEM.add_bc(coord[4],0,0.0);
						FEM.add_bc(coord[5],0,0.0);
					}
					if (x+dx == DOMAIN) {
						//FEM.add_bc(coord[2],1,0.0);
						//FEM.add_bc(coord[3],1,0.0);
						//FEM.add_bc(coord[6],1,0.0);
						//FEM.add_bc(coord[7],1,0.0);
						FEM.add_bc(coord[2],0,0.0);
						FEM.add_bc(coord[3],0,0.0);
						FEM.add_bc(coord[6],0,0.0);
						FEM.add_bc(coord[7],0,0.0);
					}
					if (y == -DOMAIN) {
						//FEM.add_bc(coord[0],0,0.0);
						//FEM.add_bc(coord[2],0,0.0);
						//FEM.add_bc(coord[4],0,0.0);
						//FEM.add_bc(coord[6],0,0.0);
						FEM.add_bc(coord[0],1,0.0);
						FEM.add_bc(coord[2],1,0.0);
						FEM.add_bc(coord[4],1,0.0);
						FEM.add_bc(coord[6],1,0.0);
					}
					if (y+dy == DOMAIN) {
						//FEM.add_bc(coord[1],0,0.0);
						//FEM.add_bc(coord[3],0,0.0);
						//FEM.add_bc(coord[5],0,0.0);
						//FEM.add_bc(coord[7],0,0.0);
						FEM.add_bc(coord[1],1,0.0);
						FEM.add_bc(coord[3],1,0.0);
						FEM.add_bc(coord[5],1,0.0);
						FEM.add_bc(coord[7],1,0.0);
					}
					if (z < -4000.0) {
						FEM.add_bc(coord[0],0,0.0);
						FEM.add_bc(coord[1],0,0.0);
						FEM.add_bc(coord[2],0,0.0);
						FEM.add_bc(coord[3],0,0.0);
						FEM.add_bc(coord[0],1,0.0);
						FEM.add_bc(coord[1],1,0.0);
						FEM.add_bc(coord[2],1,0.0);
						FEM.add_bc(coord[3],1,0.0);
						FEM.add_bc(coord[0],2,0.0);
						FEM.add_bc(coord[1],2,0.0);
						FEM.add_bc(coord[2],2,0.0);
						FEM.add_bc(coord[3],2,0.0);
					}
				}
			}
		}
		xm = -MIN(16*dx,DOMAIN); xp = MIN(16*dx,DOMAIN);
		ym = -MIN(16*dy,DOMAIN); yp = MIN(16*dy,DOMAIN);
		delta *= 2;
	} 
	FEM.solve();

	/*int i, nnd = FEM.getnodes();
	LEsystem test;
	test.addlayer(fabs(zm),100e3,0.35);
	test.addload(point2d(0.0,0.0),0.0,690.0,100.0);
	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z;
		if (!(x == 0.0 || y == 0.0 || z == 0.0))
			continue;
		if (z <= zm)
			continue;
		test.addpoint(point3d(x,y,-z));
	}
	test.calculate(LEsystem::fast);
	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z;
		if (!(x == 0.0 || y == 0.0 || z == 0.0))
			continue;
		if (z <= zm)
			continue;
		const pavedata & d = test.result(point3d(x,y,-z));
		double ux = n.ux;
		double uy = n.uy;
		double uz = n.uz;
		double vx =  d.result(pavedata::deflct,pavedata::xx);
		double vy =  d.result(pavedata::deflct,pavedata::yy);
		double vz = -d.result(pavedata::deflct,pavedata::zz);
		int j = FEM.hasnode(n);
		double h = hypot(hypot(vx-ux,vy-uy),vz-uz);
		double v = hypot(hypot(vx,vy),vz);
		printf("Node %i: (%i,%i,%i) =\t(%4.2f,%4.2f,%4.2f)\t(%4.2f,%4.2f,%4.2f)\t%4.2f\t(%4.2f)\n",j,int(x),int(y),int(z),ux,uy,uz,vx,vy,vz,h,(v == 0.0 ? 0.0 : h/v));
	}*/

#if !defined(_MSC_VER) && !defined(DARWIN)
	// calculate run time
	clock_gettime(CLOCK_PROF,&stop);
	double run_time = (stop.tv_sec - start.tv_sec) + double(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
	fprintf(stdout,"%f\n",run_time);
#endif

	return 0;
}

int
main_test()
{
#if !defined(_MSC_VER) && !defined(DARWIN)
	// get starting time    
	struct timespec start, stop;
	clock_gettime(CLOCK_PROF,&start);
#endif

	material m;
	m.setprop(material_property::emod,1000);
	m.setprop(material_property::poissons,0.2);

	const double domain[3][2] = {{-10, 10}, {-10, 10}, {-10, 0}};
	const int ndiv[3] = {8, 8, 4};
	int i, j, k;
	mesh FEM;
	
	double dx = (domain[0][1]-domain[0][0])/ndiv[0];
	double dy = (domain[1][1]-domain[1][0])/ndiv[1];
	double dz = (domain[2][1]-domain[2][0])/ndiv[2];
	fset<coord3d> coord(8);

		for (j = 0; j < ndiv[1]; j++) {
			double y = domain[1][0] + j*dy;
	for (i = 0; i < ndiv[0]; i++) {
		double x = domain[0][0] + i*dx;
			for (k = 0; k < ndiv[2]; k++) {
				double z = domain[2][0] + k*dz;
				coord[0] = coord3d(x   ,y   ,z   );
				coord[1] = coord3d(x   ,y+dy,z   );
				coord[2] = coord3d(x+dx,y   ,z   );
				coord[3] = coord3d(x+dx,y+dy,z   );
				coord[4] = coord3d(x   ,y   ,z+dz);
				coord[5] = coord3d(x   ,y+dy,z+dz);
				coord[6] = coord3d(x+dx,y   ,z+dz);
				coord[7] = coord3d(x+dx,y+dy,z+dz);
				//printf("%4.2f\t%4.2f\t%4.2f\n",x,y,z);
				FEM.add(element::block8,m,coord);
			}
		}
	}
	for (i = 0; i <= ndiv[0]; i++) {
		for (j = 0; j <= ndiv[1]; j++) {
			double x = domain[0][0] + i*dx;
			double y = domain[1][0] + j*dy;
			FEM.add_bc(coord3d(x,y,domain[2][0]),2,1.0);
		}
	}
	//FEM.add_bc(coord3d(0.0,0.0,domain[2][0]),0,0.0);
	//FEM.add_bc(coord3d(0.0,0.0,domain[2][0]),1,0.0);
	//FEM.add_bc(coord3d(domain[0][0],0.0,domain[2][0]),1,0.0);
	FEM.add_bc(coord3d(domain[0][0],domain[1][0],domain[2][0]),0,0.0);
	FEM.add_bc(coord3d(domain[0][0],domain[1][0],domain[2][0]),1,0.0);
	FEM.add_bc(coord3d(domain[0][1],domain[1][0],domain[2][0]),1,0.0);
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
			FEM.add_fext(coord3d(x,y,domain[2][1]),2,f);
		}
	}
	FEM.solve();
	int nnd = FEM.getnodes();
	for (i = 0; i < nnd; i++) {
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
		printf("Node %i: (%i,%i,%i) =\t(%f,%f,%f)\n",j,int(x),int(y),int(z),ux,uy,uz);
	}

#if !defined(_MSC_VER) && !defined(DARWIN)
	// calculate run time
	clock_gettime(CLOCK_PROF,&stop);
	double run_time = (stop.tv_sec - start.tv_sec) + double(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
	fprintf(stdout,"%f\n",run_time);
#endif

	return 0;
}

int
main_junk()
{
	double x, y, a = 0.0, b;
	double dx, dy, delta = 4.0;
	
	// Start with the tyre grid.
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


