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
		thesis work, which requires and fast linear elastic 3D FEM
		analysis.  The code is fast.  It does most of the maths using
		expression templates, and to aid that uses small 3x3 matrices to
		store the DOF level data.  This works because the math for the
		FEM is really a tensor based math, with the matrix and vectors
		being a simplification.  In essence, this code solves it as two
		nested matrices, to make a forth order tensor.  This also has an
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
 * Templated Kronecker delta
 */
template<int I, int J>
static inline int delta()
{
	return (I == J ? 1 : 0);
}

/*
 * struct meta_stiffness - template meta class to assign the stiffness
 * tensor based on lambda and mu.  DO NOT PLAY WITH THIS.
 */
template<unsigned I, unsigned J, unsigned K, unsigned L>
struct meta_stiffness
{
	static inline void assignR(tmatrix<double,NDOF,NDOF> & E,
			const double & l, const double & m) {
		E(J,L) = l*delta<I,J>()*delta<K,L>()
				+ m*(delta<I,K>()*delta<J,L>()
				   + delta<I,L>()*delta<J,K>());
		meta_stiffness<I,J,K,L-1>::assignR(E,l,m);
	}
	static inline void assign(tmatrix<double,NDOF,NDOF> & E,
			const double & l, const double & m) {
		meta_stiffness<I,J,K,L>::assignR(E,l,m);
		meta_stiffness<I,J-1,K,L>::assign(E,l,m);
	}
};
template<unsigned I, unsigned J, unsigned K>
struct meta_stiffness<I,J,K,-1>
{
	static inline void assignR(tmatrix<double,NDOF,NDOF> &,
			const double &, const double &) {}
};
template<unsigned I, unsigned K, unsigned L>
struct meta_stiffness<I,-1,K,L>
{
	static inline void assign(tmatrix<double,NDOF,NDOF> &,
			const double &, const double &) {}
};

/*
 * struct material_property
 *
 * A simple wrapper class around an enum of material property names,
 * to enable them to be used in a set, so we can add more properties
 * as we go.  The next step will be to work on how to make this a 3D
 * spatial estimator.
 */
struct material_property {
	enum property_t {
		emod,
		poissons,
	} property;

	inline material_property(const property_t p)
	  : property(p) {
	}
	inline bool operator== (const material_property & p) {
		return (property == p.property);
	}
};

/*
 * struct material_property_value.
 *
 * A simple class to store the properties value (as a double).
 */
struct material_property_value
  : public material_property
{
	double value;

	inline material_property_value(
		const material_property::property_t p, const double & v)
	  : material_property(p), value(v) {
	}
	inline operator double() {
		return value;
	}
};

/*
 * class material
 *
 * This class represents a material.  It will eventually grow to
 * support a wide variety of material types.  At the moment it
 * really is overkill...
 */
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

	// Get the point stiffness matrix.  This should be made to cache
	// the result.
	inline void pointstiffness(
			tmatrix<double,NDOF,NDOF> (& E)[NDIM][NDIM]) const {
		double e = getprop(material_property::emod);
		double v = getprop(material_property::poissons);
		double lambda = v*e/(1+v)/(1-2*v);
		double mu = e/2/(1+v);

		meta_stiffness<0,2,0,2>::assign(E[0][0],lambda,mu);
		meta_stiffness<0,2,1,2>::assign(E[0][1],lambda,mu);
		meta_stiffness<0,2,2,2>::assign(E[0][2],lambda,mu);
		meta_stiffness<1,2,0,2>::assign(E[1][0],lambda,mu);
		meta_stiffness<1,2,1,2>::assign(E[1][1],lambda,mu);
		meta_stiffness<1,2,2,2>::assign(E[1][2],lambda,mu);
		meta_stiffness<2,2,0,2>::assign(E[2][0],lambda,mu);
		meta_stiffness<2,2,1,2>::assign(E[2][1],lambda,mu);
		meta_stiffness<2,2,2,2>::assign(E[2][2],lambda,mu);
	}

private:
	ksset<material_property,material_property_value> props;
};

/*
 * struct coord3d
 *
 * A mesh point in 3D space.  This is not a point3d for to reasons:
 * 1. It sorts differently.
 * 2. We use a fixed point type to make sure the nodes line up.
 */
struct coord3d {
	fixed<8> x, y, z;

	coord3d() {
	}
	coord3d(double px, double py, double pz)
	  : x(px), y(py), z(pz) {
	}
	coord3d(const coord3d & p)
	  : x(p.x), y(p.y), z(p.z) {
	}
	coord3d(const point3d & p)
	  : x(p.x), y(p.y), z(p.z) {
	}
	~coord3d () {
	}
	int compare(const coord3d & p) const {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if (z < p.z) return -1;
		if (z > p.z) return  1;
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
 * struct gauss3d
 *
 * A Gauss Point in 3D.  This uses point3d to get full double accuracy. 
 * This should grow to support a material property list for local
 * properties.
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

/*
 * struct node3d
 *
 * A node in 3D, based on coord3d.  This stores the meshes neighbours,
 * so it can keep track of points for the variable node elements.
 */
struct node3d : public coord3d {
	union {
		struct {
			int xm, xp, ym, yp, zm, zp;
		};
		struct {
			double ux, uy, uz;
		};
	};

	node3d()
	  : coord3d(), xm(-1), xp(-1), ym(-1), yp(-1), zm(-1), zp(-1) {
	}
	node3d(const coord3d & c)
	  : coord3d(c), xm(-1), xp(-1), ym(-1), yp(-1), zm(-1), zp(-1) {
	}
	void setneighbours(const int x_m, const int x_p, const int y_m,
			const int y_p, const int z_m, const int z_p) {
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

class mesh;                            // Forward declare.

/*
 * class smatrix_elem - an element stiffness matrix
 */
class smatrix_elem {
public:
	explicit smatrix_elem(const int n)
	  : nnd(n), K(0) {
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
		assert(j >= i);
		return K[index(i,j)];
	}

private:
	friend class mesh;
	const int nnd;
	tmatrix<double,NDOF,NDOF> * K;

	// Size of the triangular matrix storage
	inline int size() const {
		return nnd*(nnd+1)/2;
	}
	// Index into the triangular matrix storage
	inline int index(const int i, const int j) const {
		return j*(j+1)/2 + i;
	}
};

const double H2x[4] = {-1,-1,+1,+1};
const double H2y[4] = {-1,+1,-1,+1};
const double V2z[2] = {-1,+1};
const double V4z0[4] = {-1,+9,+9,-1};
const double V4z2[4] = {+9,-9,-9,+9};
const double V4z1[4] = {-1,-3,+3,+1};

/*
 * class element - a finite element
 */
class element : public listelement_o<mesh,element> {
public:
	enum element_t {
		block8,
		block16,
		block34,
		infinite16
	};
	enum shape_t {
		linear,
		quadratic,
		cubic,
		absolute,
		inf_pos,
		inf_neg
	};

	element(mesh * o, element * p, const shape_t x, const shape_t y,
			const shape_t z, const material & m)
	  : listelement_o<mesh,element>(o,p), sx(x), sy(y), sz(z), mat(m),
	  		inel(0,8) {
	}
	virtual smatrix_elem * stiffness() const = 0;

protected:
	friend class mesh;

	shape_t sx, sy, sz;
	const material & mat;
	sset<int> inel;

	inline int addnode(const coord3d & c) const;
	inline void updatenode(const node3d & n) const;
	inline const node3d & getnode(const int i) const;
	inline void setup(const fset<coord3d> & c, int * mask = 0) {
		assert(8 == c.length());
		const int nz = getnz(sz);
		for (int i = 0; i < nz; i++) {
			for (int j = 0; j < 4; j++) {
				double x = double(c[j].x)+i*double(c[j+4].x-c[j].x)/(nz-1);
				double y = double(c[j].y)+i*double(c[j+4].y-c[j].y)/(nz-1);
				double z = double(c[j].z)+i*double(c[j+4].z-c[j].z)/(nz-1);
				inel.add(addnode(coord3d(x,y,z)));
			}
		}
		for (int i = 0; mask && i < 4*nz+2; i++)
			mask[i] = -1;
		for (int i = 0; i < 4*nz; i++) {
			node3d n = getnode(inel[i]);
			int x_m = -1, x_p = -1;
			int y_m = -1, y_p = -1;
			int z_m = -1, z_p = -1;
			(i%4 < 2 ? x_p : x_m) = inel[4*(i/4)+(i%4+2)%4];
			(i%2 < 1 ? y_p : y_m) = inel[2*(i/2)+(i%2+1)%2];
			if (i < 4*(nz-1))
				z_p = inel[i+4];
			if (i >= 4)
				z_m = inel[i-4];
			if (mask) {
				switch (i%4) {
				case 0:
					if (n.xp != -1 && n.xp != x_p) {
						mask[i+2] = x_p = n.xp;
					}
					if (n.yp != -1 && n.yp != y_p) {
						mask[i] = y_p = n.yp;
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
						mask[i-2] = x_m = n.xm;
					}
					if (n.ym != -1 && n.ym != y_m) {
						mask[i] = y_m = n.ym;
					}
					break;
				}
			}
			n.setneighbours(x_m,x_p,y_m,y_p,z_m,z_p);
			updatenode(n);
		}
		for (int i = 0; mask && i < 4*nz+2; i++) {
			if (mask[i] != -1) {
				inel.add(mask[i]);
				mask[i] = inel.length()-1;
			} else
				mask[i] = 0;
		}
	}
	void buildxe(matrix_dense & xe) const {
		assert(xe.rows() == inel.length());
		for (int i = 0; i < inel.length(); i++) {
			const coord3d & p = getnode(inel[i]);
			xe(i,0) = p.x; xe(i,1) = p.y; xe(i,2) = p.z;
		}
	}
	void buildgp(ksset<point3d,gauss3d> & gp) const {
		int nx = 0, ny = 0, nz = 0;
		const double (* gx)[2] = 0, (* gy)[2] = 0, (* gz)[2] = 0;
		
		getgauss(sx,nx,gx);
		getgauss(sy,ny,gy);
		getgauss(sz,nz,gz);
		for (int i = 0; i < nx; i++) {
			for (int j = 0; j < ny; j++) {
				for (int k = 0; k < nz; k++) {
					gp.add(gauss3d(gx[i][0],gy[j][0],gz[k][0],
							gx[i][1]*gy[j][1]*gz[k][1]));
				}
			}
		}
	}
	void builddNdr(const point3d & gp,
			matrix_dense & dNdr, const int * mask = 0) const {
		int nz = getnz(sz);
		double rx = gp.x, ry = gp.y, rz = gp.z;
		double Nx = 0, Ny = 0, Nz = 0, dNdx = 0, dNdy = 0, dNdz = 0;

		assert(dNdr.cols() == inel.length());
		assert(mask != 0 || dNdr.cols() == 4*nz);

		for (int i = 0, j; i < nz; i++) {
			switch (sz) {
			case linear:
				Nz = (1+V2z[i]*rz)/2; dNdz = V2z[i]/2;
				break;
			case cubic:
				Nz = (V4z0[i]+V4z2[i]*rz*rz)*(1+V4z1[i]*rz)/16;
				dNdz = ((2*V4z2[i]*rz)*(1+V4z1[i]*rz)
					+ (V4z0[i]+V4z2[i]*rz*rz)*V4z1[i])/16;
				break;
			case quadratic:
			case absolute:
			case inf_pos:
			case inf_neg:
				assert(false);
				break;
			}
			for (int l = 0; l < 4; l++) {
				Nx = (1+H2x[l]*rx)/2; dNdx = H2x[l]/2;
				Ny = (1+H2y[l]*ry)/2; dNdy = H2y[l]/2;
				dNdr(0,i*4+l) = dNdx*Ny*Nz;
				dNdr(1,i*4+l) = Nx*dNdy*Nz;
				dNdr(2,i*4+l) = Nx*Ny*dNdz;
			}
			if (!mask)
				continue;
			for (int l = 0; l < 4; l++) {
				if ((j = mask[i*4+l]) == 0)
					continue;
				const int Sxy[4] = { 0, 1, 1, 0};
				Nx = (1+( Sxy[l] ? -fabs(rx) : H2x[l]*rx));
				dNdx = ( Sxy[l] ? -SGN(rx) : H2x[l]);
				Ny = (1+(!Sxy[l] ? -fabs(ry) : H2y[l]*ry));
				dNdy = (!Sxy[l] ? -SGN(ry) : H2y[l]);
				dNdr(0,j) = dNdx*Ny*Nz/2;
				dNdr(1,j) = Nx*dNdy*Nz/2;
				dNdr(2,j) = Nx*Ny*dNdz/2;
			}
			if (i > 0 && i < nz-1)
				continue;
			if ((j = mask[(i == 0 ? 4*nz : 4*nz+1)]) == 0)
				continue;
			Nx = 1-fabs(rx); dNdx = -SGN(rx);
			Ny = 1-fabs(ry); dNdy = -SGN(ry);
			dNdr(0,j) = dNdx*Ny*Nz;
			dNdr(1,j) = Nx*dNdy*Nz;
			dNdr(2,j) = Nx*Ny*dNdz;
		}
		if (!mask)
			return;
		for (int l = 0; l < 4*nz; l++) {
			if (mask[l] != 0) {
				dNdr(0,l) -= dNdr(0,mask[l])/2;
				dNdr(1,l) -= dNdr(1,mask[l])/2;
				dNdr(2,l) -= dNdr(2,mask[l])/2;
			}
			const int Jxy[4] = {+2,-1,+1,-2};
			if (mask[l+Jxy[l%4]] != 0) {
				dNdr(0,l) -= dNdr(0,mask[l+Jxy[l%4]])/2;
				dNdr(1,l) -= dNdr(1,mask[l+Jxy[l%4]])/2;
				dNdr(2,l) -= dNdr(2,mask[l+Jxy[l%4]])/2;
			}
			if (l < 4 && mask[4*nz] != 0) {
				dNdr(0,l) += dNdr(0,mask[4*nz])/4;
				dNdr(1,l) += dNdr(1,mask[4*nz])/4;
				dNdr(2,l) += dNdr(2,mask[4*nz])/4;
			}
			if (l >= 12 && mask[4*nz+1] != 0) {
				dNdr(0,l) += dNdr(0,mask[4*nz+1])/4;
				dNdr(1,l) += dNdr(1,mask[4*nz+1])/4;
				dNdr(2,l) += dNdr(2,mask[4*nz+1])/4;
			}
		}
		for (int l = 0; mask[4*nz] != 0 && l < 4; l++) {
			dNdr(0,mask[l]) -= dNdr(0,mask[4*nz])/2;
			dNdr(1,mask[l]) -= dNdr(1,mask[4*nz])/2;
			dNdr(2,mask[l]) -= dNdr(2,mask[4*nz])/2;
		}
		for (int l = 4*(nz-1); mask[4*nz+1] != 0 && l < 4*nz; l++) {
			dNdr(0,mask[l]) -= dNdr(0,mask[4*nz+1])/2;
			dNdr(1,mask[l]) -= dNdr(1,mask[4*nz+1])/2;
			dNdr(2,mask[l]) -= dNdr(2,mask[4*nz+1])/2;
		}
	}
	smatrix_elem * buildKe(const int * mask = 0) const {
		int i, j, k, l, g, nnd = inel.length();

		// Build the point stiffness tensor E_abcd
		tmatrix<double,NDOF,NDOF> E[NDIM][NDIM];
		mat.pointstiffness(E);

		// Element nodal coords
		matrix_dense xe(nnd,NDIM);
		buildxe(xe);

		ksset<point3d,gauss3d> gp(0,8);
		buildgp(gp);

		smatrix_elem * K = new smatrix_elem(nnd);
		if (K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}

		matrix_dense dNdr(NDIM,nnd);
		matrix_dense J(NDIM,NDIM);
		for (g = 0; g < gp.length(); g++) {
			double gw = gp[g].gw;
			builddNdr(gp[g],dNdr,mask);
			J = dNdr*xe;
			// This returns det(J);
			gw *= inv_mul_gauss(NDIM,nnd,&J(0,0),&dNdr(0,0));
			for (i = 0; i < nnd; i++) {
				for (j = i; j < nnd; j++) {
					for (k = 0; k < NDIM; k++)
						for (l = 0; l < NDIM; l++)
							(*K)(i,j) += E[k][l]*(dNdr(k,i)*dNdr(l,j)*gw);
				}
			}
		}
		return K;
	}

private:
	static const int getnz(const shape_t z) {
		switch (z) {
		case linear:    return 2;
		case cubic:     return 4;
		case quadratic:
		case absolute:
		case inf_pos:
		case inf_neg:   break;
		}
		assert(false);  return 0;
	}
	static void getgauss(const shape_t t, int & n,
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
};

/*
 * class element_block8 - a finite element
 */
class element_block8 : public element {
public:
	element_block8(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,linear,linear,linear,m) {
		setup(c);
	}
	virtual smatrix_elem * stiffness() const {
		assert(8 == inel.length());
		return buildKe();
	}
};

/*
 * class element_block16 - a finite element
 */
class element_block16 : public element {
public:
	element_block16(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,linear,linear,cubic,m) {
		setup(c);
	}
	virtual smatrix_elem * stiffness() const {
		assert(16 == inel.length());
		return buildKe();
	}
};

/*
 * class element_block34 - a finite element
 */
class element_block34 : public element {
public:
	element_block34(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,absolute,absolute,cubic,m) {
		setup(c,mask);
	}
	virtual smatrix_elem * stiffness() const {
		return buildKe(mask);
	}

protected:
	int mask[18];
};

/*
 * class element_infinite16 - a infinite element
 */
class element_infinite16 : public element {
public:
	element_infinite16(mesh * o, element * p, const material & m,
			const fset<coord3d> & c)
	  : element(o,p,linear,linear,cubic,m) {
		fset<coord3d> cc(8);
		if (c.length() == 2) {
			// Corner infinite elements are defined by 2 points in
			// a line, which must be vertical.
			assert(c[0].x == c[1].x);
			assert(c[0].y == c[1].y);
			for (int i = 0; i < 4; i++) {
				double x = c[0].x, y = c[0].y;
				x *= ((x > 0) == (i%4 < 2) ? 1 : 2);
				y *= ((y > 0) == (i%2 < 1) ? 1 : 2);
				cc[i  ].x = x; cc[i  ].y = y; cc[i  ].z = c[0].z;
				cc[i+4].x = x; cc[i+4].y = y; cc[i+4].z = c[1].z;
			}
			sx = (double(c[0].x) > 0 ? inf_pos : inf_neg);
			sy = (double(c[0].y) > 0 ? inf_pos : inf_neg);
		} else {
			// Side infinite elements are defined by 4 points in a plane.
			assert(4 == c.length());
			if (c[0].x == c[3].x) {
				assert(c[0].x == c[1].x);
				assert(c[1].x == c[2].x);
				assert(c[2].x == c[3].x);
				for (int i = 0; i < 4; i++) {
					double x = c[0].x;
					x *= ((x > 0) == (i%4 < 2) ? 1 : 2);
					cc[i  ] = c[i%2  ]; cc[i  ].x = x;
					cc[i+4] = c[i%2+2]; cc[i+4].x = x;
				}
				sx = (double(c[0].x) > 0 ? inf_pos : inf_neg);
			} else {
				assert(c[0].y == c[1].y);
				assert(c[1].y == c[2].y);
				assert(c[2].y == c[3].y);
				for (int i = 0; i < 4; i++) {
					double y = c[0].y;
					y *= ((y > 0) == (i%2 < 1) ? 1 : 2);
					cc[i  ] = c[(i<2?0:1)]; cc[i  ].y = y;
					cc[i+4] = c[(i<2?2:3)]; cc[i+4].y = y;
				}
				sy = (double(c[0].y) > 0 ? inf_pos : inf_neg);
			}
		}
		setup(cc);
	}
	virtual smatrix_elem * stiffness() const {
		assert(16 == inel.length());
		int i, j, k, l, g;
		double rx, ry, rz;

		// Build the point stiffness tensor E_abcd
		tmatrix<double,NDOF,NDOF> E[NDIM][NDIM];
		mat.pointstiffness(E);

		matrix_dense xe(16,NDIM);
		buildxe(xe);

		ksset<point3d,gauss3d> gp(0,64);
		buildgp(gp);

		smatrix_elem * K = new smatrix_elem(16);
		if (K == 0) {
			event_msg(EVENT_ERROR,"Out of memory in element::stiffness()!");
			return 0;
		}

		// mapping functions for infinite elements.
		const double Mx0[16] = { 0, 0, 1, 1};
		const double Mx1[16] = {-2,-2,+1,+1};
		const double My0[16] = { 0, 1, 0, 1};
		const double My1[16] = {-2,+1,-2,+1};

		matrix_dense dNdr(NDIM,16);
		matrix_dense J(NDIM,NDIM);
		for (g = 0; g < gp.length(); g++) {
			rx = gp[g].x; ry = gp[g].y; rz = gp[g].z;
			double gw = gp[g].gw;
			for (l = 0; l < 16; l++) {
				double Nx = 0.0, dNdx = 0.0, Ny = 0.0, dNdy = 0.0;
				switch (sx) {
				case linear:
				case absolute:
					Nx = (1+H2x[l%4]*rx);
					dNdx = H2x[l%4];
					break;
				case inf_pos: {
						double rxi = 1.0/(1-rx);
						Nx = (Mx0[l%4]+Mx1[l%4]*rx)*rxi;
						dNdx = H2x[l%4]*2*rxi*rxi;
					} break;
				case inf_neg: {
						double rxi = 1.0/(1-rx);
						Nx = (Mx0[(l+2)%4]+Mx1[(l+2)%4]*rx)*rxi;
						dNdx = H2x[(l+2)%4]*2*rxi*rxi;
					} break;
				case quadratic:
				case cubic:
					assert(false);
					break;
				}
				switch (sy) {
				case linear:
				case absolute:
					Ny = (1+H2y[l%4]*ry);
					dNdy = H2y[l%4];
					break;
				case inf_pos: {
						double ryi = 1.0/(1-ry);
						Ny = (My0[l%4]+My1[l%4]*ry)*ryi;
						dNdy = H2y[l%4]*2*ryi*ryi;
					} break;
				case inf_neg: {
						double ryi = 1.0/(1-ry);
						Ny = (My0[(l+1)%4]+My1[(l+1)%4]*ry)*ryi;
						dNdy = H2y[(l+1)%4]*2*ryi*ryi;
					} break;
				case quadratic:
				case cubic:
					assert(false);
					break;
				}
				double Nz = (V4z0[l/4]+V4z2[l/4]*rz*rz)*(1+V4z1[l/4]*rz)/16;
				double dNdz = ((2*V4z2[l/4]*rz)*(1+V4z1[l/4]*rz)
					+ (V4z0[l/4]+V4z2[l/4]*rz*rz)*V4z1[l/4])/16;
				dNdr(0,l) = dNdx*Ny*Nz;
				dNdr(1,l) = Nx*dNdy*Nz;
				dNdr(2,l) = Nx*Ny*dNdz;
			}
			J = dNdr*xe;

			for (l = 0; l < 16; l++) {
				double Nx = 0.0, dNdx = 0.0, Ny = 0.0, dNdy = 0.0;
				switch (sx) {
				case linear:
				case absolute:
					Nx = (1+H2x[l%4]*rx);
					dNdx = H2x[l%4];
					break;
				case inf_pos:
					Nx = (rx*rx+(1-Mx0[l%4])*rx-Mx0[l%4])/(-Mx1[l%4]);
					dNdx = (2*rx+(1-Mx0[l%4]))/(-Mx1[l%4]);
					break;
				case inf_neg:
					Nx = (rx*rx+(1-Mx0[(l+2)%4])*rx-Mx0[(l+2)%4])/(-Mx1[(l+2)%4]);
					dNdx = (2*rx+(1-Mx0[(l+2)%4]))/(-Mx1[(l+2)%4]);
					break;
				case quadratic:
				case cubic:
					assert(false);
					break;
				}
				switch (sy) {
				case linear:
				case absolute:
					Ny = (1+H2y[l%4]*ry);
					dNdy = H2y[l%4];
					break;
				case inf_pos:
					Ny = (ry*ry+(1-My0[l%4])*ry-My0[l%4])/(-My1[l%4]);
					dNdy = (2*ry+(1-My0[l%4]))/(-My1[l%4]);
					break;
				case inf_neg:
					Ny = (ry*ry+(1-My0[(l+1)%4])*ry-My0[(l+1)%4])/(-My1[(l+1)%4]);
					dNdy = (2*ry+(1-My0[(l+1)%4]))/(-My1[(l+1)%4]);
					break;
				case quadratic:
				case cubic:
					assert(false);
					break;
				}
				double Nz = (V4z0[l/4]+V4z2[l/4]*rz*rz)*(1+V4z1[l/4]*rz)/16;
				double dNdz = ((2*V4z2[l/4]*rz)*(1+V4z1[l/4]*rz)
					+ (V4z0[l/4]+V4z2[l/4]*rz*rz)*V4z1[l/4])/16;
				dNdr(0,l) = dNdx*Ny*Nz;
				dNdr(1,l) = Nx*dNdy*Nz;
				dNdr(2,l) = Nx*Ny*dNdz;
			}
			// This returns det(J);
			// The absolute value is to take care of all of negative
			// definite Jacobians above.  This is the quickest way to fix up
			// all of the shape functions.
			gw *= fabs(inv_mul_gauss(NDIM,16,&J(0,0),&dNdr(0,0)));
			for (i = 0; i < 16; i++) {
				for (j = i; j < 16; j++) {
					for (k = 0; k < NDIM; k++)
						for (l = 0; l < NDIM; l++)
							(*K)(i,j) += E[k][l]*(dNdr(k,i)*dNdr(l,j)*gw);
				}
			}
		}
		return K;
	}
};

class smatrix_diag;
class smatrix;

/*
 * class smatrix_node - a node in a stifness matrix
 */
class smatrix_node {
	explicit smatrix_node(const int I, const int J,
			const tmatrix<double,NDOF,NDOF> & t)
	  : K(t), i(I), j(J), col_next(0), col_prev(0), row_next(0) {
	}
 	explicit smatrix_node(const int I, const int J,
			const tmatrix_expr<double,tmatrix_expr_trans<double,
				tmatrix<double,NDOF,NDOF> >,NDOF,NDOF> & t)
	  : K(t), i(J), j(I), col_next(0), col_prev(0), row_next(0) {
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
	bool insert(int i, int j, const tmatrix<double,NDOF,NDOF> & t) {
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
				if (row_head)
					row_head += (temp-nodes);
				for (int k = 0; k < nnz; k++) {
					if (temp[k].row_next)
						temp[k].row_next += (temp-nodes);
					if (temp[k].col_prev)
						temp[k].col_prev->col_next += (temp-nodes);
					else {
						smatrix_diag * d = this
								- temp[k].i + temp[k].j;
						d->col_head += (temp-nodes);
					}
					if (temp[k].col_next)
						temp[k].col_next->col_prev += (temp-nodes);
				}
				nodes = temp;
			}
		}
		if (j < i) {
			new(&nodes[nnz]) smatrix_node(i,j,~t);
			swap(i,j);
		} else
			new(&nodes[nnz]) smatrix_node(i,j,t);
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
		smatrix_diag * d = this - i + j;
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
	bool isfixed;
	bool nonzero;
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
			// Pre-allocate the space.
			diag[i].nodes = static_cast<smatrix_node *>
					(malloc(d->nnz*sizeof(smatrix_node)));
			if (diag[i].nodes == 0) {
				event_msg(EVENT_ERROR,"Out of memory in smatrix::smatrix()!");
				return;
			}
			diag[i].nnd = d->nnz;
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
		if (j < i) {
			smatrix_node * p = diag[j].row_head;
			while (p && p->j < i)
				p = p->row_next;
			if (p == 0 || p->j > i) {
				if (!diag[j].insert(i,j,t))
					return false;
			} else
				p->K += ~t;
		} else {
			smatrix_node * p = diag[i].row_head;
			while (p && p->j < j)
				p = p->row_next;
			if (p == 0 || p->j > j) {
				if (!diag[i].insert(i,j,t))
					return false;
			} else
				p->K += t;
		}
		return true;
	}
	void tidy() {
		void * temp = malloc(sizeof(smatrix_node));
		if (temp == 0) {
			event_msg(EVENT_ERROR,"Out of memory in smatrix::tidy()!");
			return;
		}
		for (int i = 0; i < nnd; i++) {
			smatrix_node * nodes = diag[i].nodes;
			for (int j = 1; j < diag[i].nnz; j++) {
				for (int k = j; k > 0
						&& nodes[k-1].j > nodes[k].j; k--) {
					memcpy(temp       ,&nodes[k-1],sizeof(smatrix_node));
					memcpy(&nodes[k-1],&nodes[k]  ,sizeof(smatrix_node));
					memcpy(&nodes[k]  ,temp       ,sizeof(smatrix_node));
				}
			}
			if (diag[i].row_head)
				diag[i].row_head = nodes;
			for (int j = 0; j < diag[i].nnz; j++) {
				nodes[j].row_next = (j < diag[i].nnz-1 ? &nodes[j+1] : 0);
				if (nodes[j].col_prev)
					nodes[j].col_prev->col_next = &nodes[j];
				else
					diag[nodes[j].j].col_head = &nodes[j];
				if (nodes[j].col_next)
					nodes[j].col_next->col_prev = &nodes[j];
			}
		}
		free(temp);
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
			assert(K(0,0) >= t(0,0));
			K(0,0) = sqrt(K(0,0)-t(0,0));
			K(0,1) = (K(0,1)-t(0,1))/K(0,0);
			K(0,2) = (K(0,2)-t(0,2))/K(0,0);
			assert(K(1,1)-K(0,1)*K(0,1) >= t(1,1));
			K(1,1) = sqrt(K(1,1)-K(0,1)*K(0,1)-t(1,1));
			K(1,2) = (K(1,2)-K(0,2)*K(0,1)-t(1,2))/K(1,1);
			assert(K(2,2)-K(0,2)*K(0,2)-K(1,2)*K(1,2) >= t(2,2));
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
		element * e = 0;
		switch (t) {
		case element::block8:
			e = new element_block8(this,last,m,c);
			break;
		case element::block16:
			e = new element_block16(this,last,m,c);
			break;
		case element::infinite16:
			e = new element_infinite16(this,last,m,c);
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
	bool add_bc(const coord3d & p, const dof f, const double d) {
		int k = node.haskey(p);
		assert(k >= 0);
		if (k == -1)
			return false;
		if (f & mesh::X)
			disp_bc.add(mesh_bc(k,0,d));
		if (f & mesh::Y)
			disp_bc.add(mesh_bc(k,1,d));
		if (f & mesh::Z)
			disp_bc.add(mesh_bc(k,2,d));
		return true;
	}
	bool add_bc_plane(const dof o, const bcplane p, const fixed<8> c,
			const dof f, const double d) {
		assert(!((o & mesh::X) && (o & mesh::Y)));
		assert(!((o & mesh::X) && (o & mesh::Z)));
		assert(!((o & mesh::Y) && (o & mesh::Z)));
		for (int i = 0; i < node.length(); i++) {
			const node3d & n = node[i];
			switch (o) {
			case mesh::X:
				if ((p & mesh::at) && (n.x == c))
					add_bc(n,f,d);
				if ((p & mesh::below) && (n.x < c))
					add_bc(n,f,d);
				if ((p & mesh::above) && (n.x > c))
					add_bc(n,f,d);
				break;
			case mesh::Y:
				if ((p & mesh::at) && (n.y == c))
					add_bc(n,f,d);
				if ((p & mesh::below) && (n.y < c))
					add_bc(n,f,d);
				if ((p & mesh::above) && (n.y > c))
					add_bc(n,f,d);
				break;
			case mesh::Z:
				if ((p & mesh::at) && (n.z == c))
					add_bc(n,f,d);
				if ((p & mesh::below) && (n.z < c))
					add_bc(n,f,d);
				if ((p & mesh::above) && (n.z > c))
					add_bc(n,f,d);
				break;
			default:
				return false;
			}
		}
		return true;
	}
	bool add_fext(const coord3d & p, const dof f, const double d) {
		int k = node.haskey(p);
		assert(k >= 0);
		if (k == -1)
			return false;
		if (f & mesh::X)
			f_ext.add(mesh_bc(k,0,d));
		if (f & mesh::Y)
			f_ext.add(mesh_bc(k,1,d));
		if (f & mesh::Z)
			f_ext.add(mesh_bc(k,2,d));
		return true;
	}
	bool solve(const double tol) {
		int i, j, nnd = node.length();
		printf("Solving with %i nodes!\n",nnd);
		smatrix K(nnd);
		svector F(nnd), U(nnd), P(nnd), W(nnd), V(nnd);
		element * e = this->first;
		smatrix_diag * d;
		smatrix_node * p;

		// Start by building the force vector.
		f_ext.sort();
		for (i = 0; i < f_ext.length(); i++) {
			const mesh_bc & f = f_ext[i];
			//F(f.n)(f.i) = f.d;
			F(node.getorder(f.n))(f.i) = f.d;
		}

		// Now check for totally fixed nodes.  We're going to skip these in
		// the later steps...
		disp_bc.sort();
		for (i = 0; i < disp_bc.length(); i++) {
			if (i+2 < disp_bc.length()
			 && disp_bc[i].n == disp_bc[i+1].n
			 && disp_bc[i+1].n == disp_bc[i+2].n) {
				const mesh_bc & u0 = disp_bc[i];
				const mesh_bc & u1 = disp_bc[i+1];
				const mesh_bc & u2 = disp_bc[i+2];
				//int n = u0.n;
				int n = node.getorder(u0.n);
				F(n)(0) = u0.d; F(n)(1) = u1.d; F(n)(2) = u2.d;
				d = &(K.diag[n]);
				d->isfixed = true;
				if (u0.d != 0.0 || u1.d != 0.0 || u2.d != 0.0)
					d->nonzero = true;
				d->K = tmatrix<double,NDOF,NDOF>(1.0,true);
				i += 2; // skip the two BC's we just checked.
			}
		}
		while (e) {
			smatrix_elem * ke = e->stiffness();
			if (ke == 0)
				return false;
			for (i = 0; i < ke->nnd; i++) {
				//int gi = e->inel[i];
				int gi = node.getorder(e->inel[i]);
				if (K.diag[gi].isfixed) {
					for (j = 0; K.diag[gi].nonzero && j < ke->nnd; j++) {
						//int gj = e->inel[j];
						int gj = node.getorder(e->inel[j]);
						if (i == j)
							continue;
						if ((j < i) == (gj < gi))
							F(gj) -= ~(*ke)(i,j)*F(gi);
						else
							F(gj) -=  (*ke)(j,i)*F(gi);
					}
				} else {
					for (j = i; j < ke->nnd; j++) {
						//int gj = ke->inel[j];
						int gj = node.getorder(e->inel[j]);
						K.append(gi,gj,(*ke)(i,j));
					}
				}
			}
			delete ke;
			e = e->next;
		}
		for (i = 0; i < disp_bc.length(); i++) {
			const mesh_bc & u = disp_bc[i];
			//int n = u.n;
			int n = node.getorder(u.n);
			d = &(K.diag[n]);
			if (d->isfixed) {
				i += 2; // we know the next two can be skipped.
				continue;
			}
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

		K.tidy();
		smatrix M(K);
		M.chol();

		// CG solution
		int it = 0;
		double r = 0.0, ro = 0.0;
		for (i = 0; i < nnd; i++)
			r += tmatrix_scalar<double>(~F(i)*F(i));
		double ri = r;
		while (r > tol*ri && it < nnd*NDOF) {
			r = 0.0;
			for (i = 0; i < nnd; i++) {
				d = &(M.diag[i]);
				tmatrix<double,NDOF,1> t(F(i));
				p = d->col_head;
				while (p) {
					t -= ~(p->K)*V(p->i);
					p = p->col_next;
				}
				t(0) = t(0)/d->K(0,0);
				t(1) = (t(1)-t(0)*d->K(0,1))/d->K(1,1);
				t(2) = (t(2)-t(0)*d->K(0,2)-t(1)*d->K(1,2))/d->K(2,2);
				V(i) = t;
			}
			for (i = nnd-1; i >= 0; i--) {
				d = &(M.diag[i]);
				tmatrix<double,NDOF,1> t(V(i));
				p = d->row_head;
				while (p) {
					t -= (p->K)*V(p->j);
					p = p->row_next;
				}
				t(2) = t(2)/d->K(2,2);
				t(1) = (t(1)-t(2)*d->K(1,2))/d->K(1,1);
				t(0) = (t(0)-t(1)*d->K(0,1)-t(2)*d->K(0,2))/d->K(0,0);
				V(i) = t;
				r += tmatrix_scalar<double>(~t*F(i));
			}
			if (it == 0) {
				for (i = 0; i < nnd; i++)
					P(i) = V(i);
			} else {
				for (i = 0; i < nnd; i++)
					P(i) = V(i) + (r/ro)*P(i);
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
		return true;
	}

private:
	friend class listelement_o<mesh,element>;
	friend class element;
	kiset<coord3d,node3d> node;
	koset<mesh_bc_key,mesh_bc> disp_bc;
	koset<mesh_bc_key,mesh_bc> f_ext;
};

inline mesh::dof
operator| (const mesh::dof l, const mesh::dof r)
{
	return static_cast<mesh::dof>(int(l) | int(r));
}
inline mesh::dof
operator& (const mesh::dof l, const mesh::dof r)
{
	return static_cast<mesh::dof>(int(l) & int(r));
}

inline mesh::bcplane
operator| (const mesh::bcplane l, const mesh::bcplane r)
{
	return static_cast<mesh::bcplane>(int(l) | int(r));
}
inline mesh::bcplane
operator& (const mesh::bcplane l, const mesh::bcplane r)
{
	return static_cast<mesh::bcplane>(int(l) & int(r));
}

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

	int x, y, z;
	int dx, dy, dz, delta = 4, step = 1;
	mesh FEM;
	fset<coord3d> coord(8);

	// Find the step size, which depends on our tyre.
	while ((step-2)*delta < 100)
		step *= 2;
	// Start with the tyre grid.
	dx = delta; dy = delta;
	for (x = -(step-2)*dx; x < (step-2)*dx; x += 2*dx) {
		for (y = -(step-2)*dy; y < (step-2)*dy; y += 2*dy) {
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
	for (x = -(step-2)*dx; x <= (step-2)*dx; x += dx) {
		for (y = -(step-2)*dy; y <= (step-2)*dy; y += dy) {
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
				FEM.add_fext(coord3d(x,y,0.0),mesh::Z,f);
		}
	}
	step /= 2;

#define DOMAIN (4096+2048)

	// Now add the elements from below the tyre, working outwards.
	int xm = 0, xp = 0, ym = 0, yp = 0, zp = 0, zm = 0;
	while (zp > -4000 || xp < DOMAIN) {
		zp += MIN(-30,zp);
		dx = delta*2; dy = delta*2;
		z = 0;
		while (z > zp) {
			dz = MIN(-30,z);
			z += dz;
			//printf("%i\t%i\t%i\t%i\t%i\n",MIN(16*dx,DOMAIN),xp,zm,z,zp);
			for (x = -MIN(16*dx,DOMAIN); x < MIN(16*dx,DOMAIN); x += dx) {
				for (y = -MIN(16*dy,DOMAIN); y < MIN(16*dy,DOMAIN); y += dy) {
					if (xm < xp && (x >= xm && x < xp)
					 && ym < yp && (y >= ym && y < yp)
					 && z >= zm)
						continue;
					coord.resize(8);
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
						if (y == -DOMAIN) {
							coord.resize(2);
							coord[0] = coord3d(x,y,z);
							coord[1] = coord3d(x,y,z-dz);
							FEM.add(element::infinite16,m,coord);
						} else if (y+dy == DOMAIN) {
							coord.resize(2);
							coord[0] = coord3d(x,y+dy,z);
							coord[1] = coord3d(x,y+dy,z-dz);
							FEM.add(element::infinite16,m,coord);
						}
						if (y+dy <= DOMAIN) {
							coord.resize(4);
							coord[0] = coord3d(x,y   ,z);
							coord[1] = coord3d(x,y+dy,z);
							coord[2] = coord3d(x,y   ,z-dz);
							coord[3] = coord3d(x,y+dy,z-dz);
							FEM.add(element::infinite16,m,coord);
						}
					} else if (x+dx == DOMAIN) {
						if (y == -DOMAIN) {
							coord.resize(2);
							coord[0] = coord3d(x+dx,y,z);
							coord[1] = coord3d(x+dx,y,z-dz);
							FEM.add(element::infinite16,m,coord);
						} else if (y+dy == DOMAIN) {
							coord.resize(2);
							coord[0] = coord3d(x+dx,y+dy,z);
							coord[1] = coord3d(x+dx,y+dy,z-dz);
							FEM.add(element::infinite16,m,coord);
						}
						if (y+dy <= DOMAIN) {
							coord.resize(4);
							coord[0] = coord3d(x+dx,y   ,z);
							coord[1] = coord3d(x+dx,y+dy,z);
							coord[2] = coord3d(x+dx,y   ,z-dz);
							coord[3] = coord3d(x+dx,y+dy,z-dz);
							FEM.add(element::infinite16,m,coord);
						}
					}
					if (x+dx <= DOMAIN) {
						if (y == -DOMAIN) {
							coord.resize(4);
							coord[0] = coord3d(x   ,y,z);
							coord[1] = coord3d(x+dx,y,z);
							coord[2] = coord3d(x   ,y,z-dz);
							coord[3] = coord3d(x+dx,y,z-dz);
							FEM.add(element::infinite16,m,coord);
						}
						if (y+dy == DOMAIN) {
							coord.resize(4);
							coord[0] = coord3d(x   ,y+dy,z);
							coord[1] = coord3d(x+dx,y+dy,z);
							coord[2] = coord3d(x   ,y+dy,z-dz);
							coord[3] = coord3d(x+dx,y+dy,z-dz);
							FEM.add(element::infinite16,m,coord);
						}
					}
				}
			}
		}
		xm = -MIN(step*dx,DOMAIN); xp = MIN(step*dx,DOMAIN);
		ym = -MIN(step*dy,DOMAIN); yp = MIN(step*dy,DOMAIN);
		if (xm > -DOMAIN && xp < DOMAIN
		 && ym > -DOMAIN && xp < DOMAIN)
			delta *= 2;
		zm = z;
	}
	//FEM.add_bc_plane(mesh::X,mesh::at|mesh::below,-DOMAIN,
	//		mesh::X,0.0);
	//FEM.add_bc_plane(mesh::X,mesh::at|mesh::above, DOMAIN,
	//		mesh::X,0.0);
	//FEM.add_bc_plane(mesh::Y,mesh::at|mesh::below,-DOMAIN,
	//		mesh::Y,0.0);
	//FEM.add_bc_plane(mesh::Y,mesh::at|mesh::above, DOMAIN,
	//		mesh::Y,0.0);
	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,zm,
			mesh::X|mesh::Y|mesh::Z,0.0);
	FEM.solve(1e-20);

	int i, nnd = FEM.getnodes();
	LEsystem test;
	test.addlayer(-zm,100e3,0.35);
	test.addload(point2d(0.0,0.0),0.0,690.0,100.0);
	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z;
		if (!(x == 0 || y == 0 || z == 0))
			continue;
		if (x < xm || x > xp || y < ym || y > yp || z <= zm)
			continue;
		test.addpoint(point3d(x,y,-z));
	}
	test.calculate(LEsystem::fast);
	for (i = 0; i < nnd; i++) {
		const node3d & n = FEM.getorderednode(i);
		x = n.x; y = n.y; z = n.z;
		if (!(x == 0 || y == 0 || z == 0))
			continue;
		if (x < xm || x > xp || y < ym || y > yp || z <= zm)
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
		printf("Node %6i: (%+6i,%+6i,%+6i) =\t(%8.2g,%8.2g,%8.2g)\t(%8.2g,%8.2g,%8.2g)\t%8.2g\t(%8.2g)\n",j,int(x),int(y),int(z),vx,vy,vz,ux,uy,uz,h,(v == 0.0 ? 0.0 : h/v));
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

	for (k = 0; k < ndiv[2]; k++) {
		double z = domain[2][1] - k*dz;
		int p = (1 << k); // 1; 
		for (j = 0; j < ndiv[1]; j += p) {
			double y = domain[1][0] + j*dy;
			for (i = 0; i < ndiv[0]; i += p) {
				double x = domain[0][0] + i*dx;
				coord[0] = coord3d(x     ,y     ,z-dz);
				coord[1] = coord3d(x     ,y+p*dy,z-dz);
				coord[2] = coord3d(x+p*dx,y     ,z-dz);
				coord[3] = coord3d(x+p*dx,y+p*dy,z-dz);
				coord[4] = coord3d(x     ,y     ,z   );
				coord[5] = coord3d(x     ,y+p*dy,z   );
				coord[6] = coord3d(x+p*dx,y     ,z   );
				coord[7] = coord3d(x+p*dx,y+p*dy,z   );
				//printf("%4.2f\t%4.2f\t%4.2f\n",x,y,z);
				//FEM.add(element::block16,m,coord);
				FEM.add(element::block34,m,coord);
			}
		}
	}
	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,domain[2][0],mesh::Z,0.0);
	FEM.add_bc(coord3d(domain[0][0],domain[1][0],domain[2][0]),mesh::X,0.0);
	FEM.add_bc(coord3d(domain[0][0],domain[1][0],domain[2][0]),mesh::Y,0.0);
	FEM.add_bc(coord3d(domain[0][1],domain[1][0],domain[2][0]),mesh::Y,0.0);
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
			FEM.add_fext(coord3d(x,y,domain[2][1]),mesh::Z,f);
		}
	}
	FEM.solve(1e-12);
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
main_infinite()
{
#if !defined(_MSC_VER) && !defined(DARWIN)
	// get starting time
	struct timespec start, stop;
	clock_gettime(CLOCK_PROF,&start);
#endif

	material m;
	m.setprop(material_property::emod,100); // kPa
	m.setprop(material_property::poissons,0.35);
	mesh FEM;

	sset<coord3d> coord(0,8);

	coord.empty();
	coord.add(coord3d(-10,-10,-10));
	coord.add(coord3d(-10, 10,-10));
	coord.add(coord3d( 10,-10,-10));
	coord.add(coord3d( 10, 10,-10));
	coord.add(coord3d(-10,-10,  0));
	coord.add(coord3d(-10, 10,  0));
	coord.add(coord3d( 10,-10,  0));
	coord.add(coord3d( 10, 10,  0));
	FEM.add(element::block16,m,coord);

	coord.empty();
	coord.add(coord3d(-10,-10,-10));
	coord.add(coord3d(-10,-10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d(-10, 10,-10));
	coord.add(coord3d(-10, 10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d(-10,-10,-10));
	coord.add(coord3d(-10, 10,-10));
	coord.add(coord3d(-10,-10,  0));
	coord.add(coord3d(-10, 10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d( 10,-10,-10));
	coord.add(coord3d( 10,-10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d( 10, 10,-10));
	coord.add(coord3d( 10, 10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d( 10,-10,-10));
	coord.add(coord3d( 10, 10,-10));
	coord.add(coord3d( 10,-10,  0));
	coord.add(coord3d( 10, 10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d(-10,-10,-10));
	coord.add(coord3d( 10,-10,-10));
	coord.add(coord3d(-10,-10,  0));
	coord.add(coord3d( 10,-10,  0));
	FEM.add(element::infinite16,m,coord);
	coord.empty();
	coord.add(coord3d(-10, 10,-10));
	coord.add(coord3d( 10, 10,-10));
	coord.add(coord3d(-10, 10,  0));
	coord.add(coord3d( 10, 10,  0));
	FEM.add(element::infinite16,m,coord);

	FEM.add_fext(coord3d(-10,-10,  0),mesh::Z,-10);
	FEM.add_fext(coord3d(-10, 10,  0),mesh::Z,-10);
	FEM.add_fext(coord3d( 10,-10,  0),mesh::Z,-10);
	FEM.add_fext(coord3d( 10, 10,  0),mesh::Z,-10);

	FEM.add_bc_plane(mesh::Z,mesh::at|mesh::below,-10,
			mesh::X|mesh::Y|mesh::Z,0.0);

	FEM.solve(1e-30);
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
