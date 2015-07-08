/**************************************************************************

	PAVEMENT.H - Interface for the pavement class.

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
		This file implements various pavement related classes.  Basically,
		this is an N layered elastic structure, with N loads and
		evaluation points.

		The calculations are based on a C++ conversion of ELSYM5M, with
		some minor changes.  It should produce the same results.  The
		units for the various properties are mentioned below.

		At the moment this class is limited to just doing the layered
		elastic calculations.

	Future:
		This class needs to be designed more carefully to accommodate
		structural and material models other than layered elastic.
		Obviously this will need another calculation engine...

	History:
		2002/01/23 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __PAVEMENT_H
#define __PAVEMENT_H

#include <memory.h>
#if defined(_MSC_VER)
#include <stddef.h>
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#else
#include <stdint.h>
#endif
#include "mathplus.h"
#include "event.h"
#include "set.h"
#include "list.h"

class LEsystem;
class LEsystem_cache;
class LEbackcalc;

/*
 * struct point2d - A simple point in 2D space.
 */
struct point2d {
	double x, y;

	point2d() {
	}
	point2d(double px, double py)
	  : x(px), y(py) {
	}
	point2d(const point2d & p)
	  : x(p.x), y(p.y) {
	}
	double r() const {
		return hypot(x,y);
	}
	double theta() const {
		return atan2(x,y);
	}
	double distance(const point2d & p) const {
		return hypot(x-p.x,y-p.y);
	}
	int compare(const point2d & p) const  {
		if (x == p.x && y == p.y)
			return 0;
		if (r() < p.r())
			return -1;
		if (r() > p.r())
			return 1;
		if (theta() < p.theta())
			return -1;
		else
			return 1;
	}
	bool operator == (const point2d & p) const {
		return (x == p.x && y == p.y ? true : false);
	}
	bool operator != (const point2d & p) const {
		return (x != p.x || y != p.y ? true : false);
	}
	bool operator > (const point2d & p) const {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const point2d & p) const {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const point2d & p) const {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const point2d & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * struct point3d - A simple point in 3D space.
 */
struct point3d : public point2d {
	// WARNING: Z is positive down.
	double z;

	point3d()
	  : point2d() {
	}
	point3d(double px, double py, double pz)
	  : point2d(px,py), z(pz) {
	}
	point3d(const point3d & p)
	  : point2d(p), z(p.z) {
	}
	int compare(const point3d & p) const {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if ((x*x+y*y) < (p.x*p.x+p.y*p.y))
			return -1;
		if ((x*x+y*y+z*z) > (p.x*p.x+p.y*p.y+p.z*p.z))
			return 1;
		return point2d::compare(static_cast<const point2d &>(p));
	}
	bool operator == (const point3d & p) const {
		return (x == p.x && y == p.y && z == p.z ? true : false);
	}
	bool operator != (const point3d & p) const {
		return (x != p.x || y != p.y || z != p.z ? true : false);
	}
	bool operator > (const point3d & p) const {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const point3d & p) const {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const point3d & p) const {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const point3d & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * class LElayer - A linear elastic pavement layer.
 *
 * This implements the basics structure of a pavement layer.
 */
class LElayer : public listelement_o<LEsystem,LElayer> {
public:
	double bottom() const {
		return (h == 0.0 ? 0.0 : top() + h);
	}
	double top() const {
		return (prev == 0 ? 0.0 : prev->bottom());
	}
	double thickness() const {
		return h;
	}
	double thickness(double th) {
		return h = th;
	}
	double emod() const {
		return E;
	}
	double emod(double te) {
		return E = te;
	}
	double poissons() const {
		return v;
	}
	double poissons(double tv) {
		return v = tv;
	}
	double slip() const {
		return s;
	}
	double slip(double ts) {
		return s = ts;
	}
	LElayer(LEsystem * o, LElayer * p)
	  : listelement_o<LEsystem,LElayer>(o,p), h(0.0), E(0.0), v(0.0),
			s(1.0) {
	}
	LElayer(LEsystem * o, LElayer * p, double th, double te, double tv,
			double ts = 1.0)
	  : listelement_o<LEsystem,LElayer>(o,p), h(th), E(te), v(tv),
			s(ts) {
	}
	LElayer(LEsystem * o, LElayer * p, const LElayer & l)
	  : listelement_o<LEsystem,LElayer>(o,p), h(l.h), E(l.E), v(l.v),
			s(l.s) {
	}

private:
	friend class LEsystem;

	double h;                           // The thickness of the layer.
	double E;                           // The elastic modulus.
	double v;                           // Poisson's ratio.
	double s;                           // Layer slip (on bottom boundary).
};

/*
 * struct paveload - A load on the pavement structure.
 */
class paveload : public point2d {
public:
	paveload()
	  : point2d(0.0,0.0), f(0.0), p(0.0) {
	}
	paveload(const point2d & l, double lf, double lp, double lr = 0.0)
	  : point2d(l), f(lf), p(lp)  {
		if (lf == 0.0) {
			f = M_PI*lr*lr*p;
		} else if (lp == 0.0) {
			p = f/(M_PI*lr*lr);
		}
	}
	paveload(const paveload & pl)
	  : point2d(pl), f(pl.f), p(pl.p) {
	}
	double force() const {
		return f;
	}
	double pressure() const {
		return p;
	}
	double radius() const {
		return (p == 0.0 ? 0.0 : sqrt(f/(p*M_PI)));
	}
	double area() const {
		return (p == 0.0 ? 0.0 : f/p);
	}

private:
	double f;                          // Force.
	double p;                          // Pressure.
};

/*
 * struct pavepoint - A point in a pavement (3D and layer)
 */
struct pavepoint : public point3d {
	unsigned il;

	pavepoint()
	  : point3d(0.0,0.0,0.0), il(UINT_MAX) {
	}
	pavepoint(double px, double py, double pz, unsigned l = UINT_MAX)
	  : point3d(px,py,pz), il(l) {
	}
	pavepoint(const point3d & p, unsigned l = UINT_MAX)
	  : point3d(p), il(l) {
	}
	pavepoint(const pavepoint & p)
	  : point3d(p), il(p.il) {
	}
	pavepoint & operator= (const pavepoint & p) {
		x = p.x; y = p.y; z = p.z;
		if (il == UINT_MAX || p.il != UINT_MAX)
			il = p.il;
		return *this;
	}
	int compare(const pavepoint & p) const {
		if (x == p.x && y == p.y && z == p.z
		  && (il == p.il || il == UINT_MAX || p.il == UINT_MAX))
			return 0;
		if (il != UINT_MAX && p.il != UINT_MAX && il < p.il)
			return -1;
		if (il != UINT_MAX && p.il != UINT_MAX && il > p.il)
			return 1;
		return point3d::compare(static_cast<const point3d &>(p));
	}
	bool operator == (const pavepoint & p) const {
		return (x == p.x && y == p.y && z == p.z
				 && (il == p.il || il == UINT_MAX || p.il == UINT_MAX)
				 ? true : false);
	}
	bool operator != (const pavepoint & p) const {
		return (x != p.x || y != p.y || z != p.z
				 || (il != p.il && il != UINT_MAX && p.il != UINT_MAX)
				 ? true : false);
	}
	bool operator > (const pavepoint & p) const {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const pavepoint & p) const {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const pavepoint & p) const {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const pavepoint & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * class pavedata - An evaluation point within the pavement structure.
 */
struct pavedata : public pavepoint {
	// Return the results based on a more rational system...
	enum type {deflct, stress, strain};
	enum direction {xx, yy, zz, xy, xz, yz, p1, p2, p3, s1, s2, s3};
	double result(type t, direction d) const {
		switch (t) {
		case stress:
			switch (d) {
			case xx: case yy: case zz:
				return data[0][d-xx];
			case xy: case xz: case yz:
				return data[1][d-xy];
			case p1: case p2: case p3:
				return data[2][d-p1];
			case s1: case s2: case s3:
				return data[3][d-s1];
			}
		case deflct:
			switch (d) {
			case xx: case yy: case zz:
				return data[4][d-xx];
			case xy: case xz: case yz:
			case p1: case p2: case p3:
			case s1: case s2: case s3:
				return 0.0;
			}
		case strain:
			switch (d) {
			case xx: case yy: case zz:
				return data[5][d-xx];
			case xy: case xz: case yz:
				return data[6][d-xy];
			case p1: case p2: case p3:
				return data[7][d-p1];
			case s1: case s2: case s3:
				return data[8][d-s1];
			}
		default:
			return 0.0;
		}
	}

	pavedata()
	  : pavepoint(), deflgrad(0) {
	}
	pavedata(const point3d & p, unsigned l = UINT_MAX)
	  : pavepoint(p,l), deflgrad(0) {
	}
	pavedata(const pavedata & pd)
	  : pavepoint(pd), deflgrad(pd.deflgrad) {
		memcpy(data,pd.data,sizeof(data));
	}
	~pavedata() {
	}

	double data[9][3];
	fset<double> deflgrad;
	void principle(double v, double E);
};

/*
 * class LEsystem - A layered elastic pavement system.
 */
class LEsystem : private list_owned<LEsystem, LElayer> {
public:
	LEsystem()
	  : list_owned<LEsystem,LElayer>(), data(), load(),
	    cache(0), cache_state(cachestate::reset) {
	    cache_reset();
	}
	LEsystem(LEsystem & p)
	  : list_owned<LEsystem,LElayer>(p), data(p.data), load(p.load),
	    cache(0), cache_state(cachestate::reset) {
	    cache_reset();
	}
	~LEsystem() {
		cache_free();
	}
	void addlayer(const double h, const double e, const double v,
	              const double s = 1.0, const unsigned p = UINT_MAX);
	void removelayer(const unsigned l);
	void removelayers();
	inline unsigned layers() const {
		return length();
	}
	LElayer & layer(const unsigned l) const;

	void addload(const point2d & l, double f, double p, double r = 0);
	void addload(const paveload & l);
	void removeload(const unsigned i);
	void removeloads();
	inline unsigned loads() const {
		return load.length();
	}
	const paveload & getload(const unsigned i) {
		return load[i];
	}

	void addpoint(const point3d & p, unsigned l = UINT_MAX);
	void addgrid(const unsigned nx, const double * xp,
				 const unsigned ny, const double * yp,
				 const unsigned nz, const double * zp);
	void removepoint(const point3d & p, unsigned l = UINT_MAX);
	void removepoints();
	inline unsigned results() const {
		return data.length();
	}
	const pavedata & result(const point3d & p, unsigned l = UINT_MAX) const {
		return data[pavepoint(p,l)];
	}
	const pavedata & result(const unsigned i) const {
		return data[i];
	}

	bool check();
	enum resulttype {
		all       = 0x0000,
		fast      = 0x0001,
		dirty     = 0x0002,
		odemark   = 0x0003,
		fastnum   = 0x0004,
		accurate  = 0x00FF,
		mask      = 0x00FF,
		disp      = 0x0100,
		grad      = 0x0200,
		dispgrad  = 0x0300,
		fastdisp  = 0x0101,
		dirtydisp = 0x0102,
		fastgrad  = 0x0301
	};
	bool calc_accurate();
	bool calculate(resulttype result = all, const double * Q = 0);
	bool calc_odemark();
	bool calc_fastnum();

protected:
	void * cache_alloc(size_t len);
	void cache_reset();
	void cache_free();
	enum cachestate {
	    reset    = 0x0000,
		geom     = 0x0001,
		pressure = 0x0002,
		emod     = 0x0004,
		checked  = 0x0008
	} cache_state;
	friend cachestate & operator&= (cachestate & l, const cachestate r);
	friend cachestate operator~ (const cachestate l);

private:
	friend class listelement_o<LEsystem, LElayer>;
	friend class list_owned<LEsystem, LElayer>;
	friend class LEbackcalc;

	ksset<pavepoint,pavedata> data;
	sset<paveload> load;
	
	LEsystem_cache * cache;
};

inline LEsystem::cachestate &
operator&= (LEsystem::cachestate & l, const LEsystem::cachestate r)
{
	return l = static_cast<LEsystem::cachestate>(unsigned(l) & unsigned(r));
}
inline LEsystem::cachestate
operator~ (const LEsystem::cachestate l)
{
	return static_cast<LEsystem::cachestate>(~unsigned(l));
}

/*
 * class defldata - A measured point within the pavement structure.
 */
struct defldata : public point3d {
	defldata()
	  : point3d() {
	}
	defldata(const point3d & p, const double m)
	  : point3d(p) {
		measured = m;
		calculated = 0.0;
	}
	defldata(const defldata & dd)
	  : point3d(dd) {
		measured = dd.measured;
		calculated = dd.calculated;
	}
	~defldata() {
	}
	double measured;
	double calculated;
};

/*
 * class LEbackcalc - A layered elastic backcalculation.
 */
class LEbackcalc : public LEsystem {
public:
	inline void adddefl(const point3d & p, double d) {
		defl.add(defldata(p,d));
	}
	inline void adddefl(const defldata & d) {
		defl.add(d);
	}
	inline unsigned deflections() {
		return defl.length();
	}
	const defldata & getdefl(const unsigned i) {
		return defl[i];
	}
	void removedeflections() {
		defl.empty();
	}
	void setup(double p, double n, double t, unsigned m) {
		precision = MAX(0.0,p);
		noise = MAX(0.0,n);
		tolerance = MAX(1e-6,t);
		maxsteps = MAX(3,m);
	}
	bool backcalc();

	LEbackcalc() {
		setup(0.0,0.0,1e-6,5);
	}
	LEbackcalc(LEbackcalc & b)
	  : LEsystem(b), defl(b.defl) {
		setup(b.precision,b.noise,b.tolerance,b.maxsteps);
	}
	~LEbackcalc() {
	}
private:
	sset<defldata> defl;
	double precision;
	double noise;
	double tolerance;
	unsigned maxsteps;

	enum calctype {slow, fast, reuse};
	bool seed(unsigned nl, double * P);
	double deflgrad(unsigned nl, double * P, double * Q,
			calctype cl = slow);
	double gaussnewton(unsigned nl, double * P, calctype cl = slow);
	double kalman(unsigned nl, double * P);
	double bowlerror(unsigned nl = 0, const double * P = 0,
			const double s = 0.0, const double * D = 0);
	double brent(unsigned nl, double * P, double * D);
	double conjgrad(unsigned nl, double * P);
	double swarm(unsigned nl, double * P);
};

#endif // PAVEMENT_H
