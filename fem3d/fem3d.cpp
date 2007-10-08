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
#include <time.h>
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

	element(mesh * o, element * p, const element_t t, const material & m)
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
			//tmatrix<double,nnd,1> H;
	
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
				//H(0,l) = (1+Hx1[l]*rx)*(1+Hy1[l]*ry)
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
				//H(0,i) = (1+(Sxy[l] ? -fabs(rx) : Hx1[l]*rx))
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
				//H(0,i) = (1-fabs(rx))*(1-fabs(ry))
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
			}
			for (l = 0; mask[16] && l < 4; l++) {
				for (i = 0; i < 3; i++)
					dHdr(i,mask[l]) -= dHdr(i,mask[16])/2;
			}
			for (l = 12; mask[17] && l < 16; l++) {
				for (i = 0; i < 3; i++)
					dHdr(i,mask[l]) -= dHdr(i,mask[17])/2;
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
	bool insert(int i, int j, const tmatrix<double,NDOF,NDOF> & t,
			smatrix_diag * d) {
		int s = nnz+1, nnb = 1;
		//while (s > 8*nnb)
		//	nnb *= 2;
		s = nnb*(s/nnb+(s%nnb?1:0));
		if (s > nnd) {
			nnd = s;
			smatrix_node * temp = static_cast<smatrix_node *>
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
				printf("Ooops...\n");
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
	inline const node3d & getnode(const int i) {
		return node[i];
	}
	inline const int hasnode(const coord3d & p) {
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
					//K.append(ke->inel[i],ke->inel[j],(*ke)(i,j));
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
			printf("Adding force %f at node %i dof %i\n",f.d,f.n,f.i);
			F(node.getorder(f.n))(f.i) = f.d;
			//F(f.n)(f.i) = f.d;
		}
		disp_bc.sort();
		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			// XXX: extract forces
			printf("Adding bc %f at node %i dof %i\n",u.d,u.n,u.i);
			d = &(K.diag[node.getorder(u.n)]);
			//d = &(K.diag[u.n]);
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

		smatrix M(K);
		M.chol();
		/*for (i = 0; i < nnd; i++) {
			d = &(M.diag[i]);
			printf("Chol %i,%i:\n",i,i);
			tmatrix<double,NDOF,NDOF> t(d->K);
			t.print();
			p = d->row_head;
			while (p) {
				assert(p->i == i);
				printf("Chol %i,%i:\n",p->i,p->j);
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
			r = 0.0;
			for (i = 0; i < nnd; i++) {
				//d = &(K.diag[i]);
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
			//for (i = 0; i < nnd; i++) {
			//	d = &(K.diag[i]);
			//	for (j = 0; j < NDOF; j++)
			//		Z(i)(j) *= d->K(j,j);
			//}
			for (i = nnd-1; i >= 0; i--) {
				//d = &(K.diag[i]);
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
					//P(i) = F(i);
			} else {
				for (i = 0; i < nnd; i++)
					P(i) = Z(i) + (r/ro)*P(i);
					//P(i) = F(i) + (r/ro)*P(i);
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
			//node3d n = node[i];
			node3d n = node.getindex(i);
			n.setdisp(U(i));
			node.replace(n);
		}
		for (i = 0; i < nnd; i++) {
			const node3d & n = node.getindex(i);
			double x = n.x;
			double y = n.y;
			double z = n.z;
			double ux = n.ux;
			double uy = n.uy;
			double uz = n.uz;
			j = node.haskey(n);
			printf("Node %i: (%f,%f,%f) = (%f,%f,%f)\n",j,x,y,z,ux,uy,uz);
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

/*
 * This program is a custom 3D finite element code, intended for
 * work on my PhD.  It will hopefully form the basis for a later
 * full 3D FEM pavement modelling code, but at the moment, I need
 * to get some work done... 
 */
int
main()
{
    // get starting time    
    struct timespec start, stop;
    clock_gettime(CLOCK_PROF,&start);

	material m;
	m.setprop(material_property::emod,1000);
	m.setprop(material_property::poissons,0.2);

	double x, y, z;
	double dx, dy, dz;
	mesh FEM;
	fset<coord3d> coord(8);
	
	// Start with the tyre grid.
	dx = 4.0; dy = 4.0;
	for (x = -128.0; x <= 128.0; x += dx) {
		for (y = -128.0; y <= 128.0; y += dy) {
			FEM.addnode(coord3d(x,y,0.0));
		}
	} 
	double F = -1*dx*dy;
	for (x = -128.0; x <= 128.0; x += dx) {
		for (y = -128.0; y <= 128.0; y += dy) {
			node3d n = FEM.getnode(FEM.hasnode(coord3d(x,y,0.0)));
			int x_m = FEM.hasnode(coord3d(x-dx,y,0.0));
			int x_p = FEM.hasnode(coord3d(x+dx,y,0.0));
			int y_m = FEM.hasnode(coord3d(x,y-dy,0.0));
			int y_p = FEM.hasnode(coord3d(x,y+dy,0.0));
			n.setneighbours(x_m,x_p,y_m,y_p,-1,-1);
			if (hypot(x,y) > 100)
				continue;
			FEM.add_fext(coord3d(x,y,0.0),2,F);
		}
	}
	// Now add the elements below the tyre.
	dx = 8.0; dy = 8.0;
	for (x = -128.0; x < 128.0; x += dx) {
		for (y = -128.0; y < 128.0; y += dy) {
			z = 0.0; dz = -30.0;
			coord[0] = coord3d(x   ,y   ,z   );
			coord[1] = coord3d(x   ,y+dy,z   );
			coord[2] = coord3d(x+dx,y   ,z   );
			coord[3] = coord3d(x+dx,y+dy,z   );
			coord[4] = coord3d(x   ,y   ,z+dz);
			coord[5] = coord3d(x   ,y+dy,z+dz);
			coord[6] = coord3d(x+dx,y   ,z+dz);
			coord[7] = coord3d(x+dx,y+dy,z+dz);
			FEM.add(element::block34,m,coord);
		}
	} 
	for (x = -128.0; x <= 128.0; x += dx) {
		for (y = -128.0; y <= 128.0; y += dy) {
			FEM.add_bc(coord3d(x,y,-30.0),2,0.0);
		}
	}
	FEM.add_bc(coord3d(   0.0,0.0,-30.0),0,0.0);
	FEM.add_bc(coord3d(   0.0,0.0,-30.0),1,0.0);
	FEM.add_bc(coord3d(-128.0,0.0,-30.0),1,0.0);
	FEM.solve();

    // calculate run time
    clock_gettime(CLOCK_PROF,&stop);
    double run_time = (stop.tv_sec - start.tv_sec) + double(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
    fprintf(stdout,"%f\n",run_time);

	return 0;
}

int
main_test()
{
    // get starting time    
    struct timespec start, stop;
    clock_gettime(CLOCK_PROF,&start);

	material m;
	m.setprop(material_property::emod,1000);
	m.setprop(material_property::poissons,0.2);

	const double domain[3][2] = {{-10, 10}, {-10, 10}, {-10, 0}};
	const int ndiv[3] = {40, 40, 20};
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
			FEM.add_bc(coord3d(x,y,domain[2][0]),2,0.0);
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

    // calculate run time
    clock_gettime(CLOCK_PROF,&stop);
    double run_time = (stop.tv_sec - start.tv_sec) + double(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
    fprintf(stdout,"%f\n",run_time);

	return 0;
}
