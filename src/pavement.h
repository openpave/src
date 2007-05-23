/**************************************************************************

	PAVEMENT.H - Interface for the pavement class.

	$OpenPave$

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
		2002/01/23 - Created by Jeremy Lea <jlea@csir.co.za>

**************************************************************************/

#ifndef __PAVEMENT_H
#define __PAVEMENT_H

#include "config.h"
#include <memory.h>
#include "mathplus.h"
#include "event.h"
#include "set.h"
#include "list.h"
#include "fixed.h"

struct LElayer;
struct paveload;
struct pavedata;
class LEsystem;
class LEbackcalc;

/*
 * A simple point in 3D space.
 */
struct point2d {
	double x, y;

	point2d() {
	};
	point2d(double px, double py) {
		x = px;
		y = py;
	};
	point2d(const point2d & p) {
		x = p.x;
		y = p.y;
	};
	~point2d () {
	};
	inline double r() const {
		return sqrt(x*x+y*y);
	};
	inline double theta() const {
		return (y >= 0.0 ? acos(x/r()) : M_2PI - acos(x/r()));
	};
	inline double distance(const point2d & p) {
		return sqrt((x-p.x)*(x-p.x)+(y-p.y)*(y-p.y));
	};
	int compare(const point2d & p) const {
		if (x == p.x && y == p.y)
			return 0;
		if (fabs(x) <= fabs(p.x) && fabs(y) <= fabs(p.y))
			return -1;
		if (fabs(x) >= fabs(p.x) && fabs(y) >= fabs(p.y))
			return 1;
		if ((x*x+y*y) < (p.x*p.x+p.y*p.y))
			return -1;
		if ((x*x+y*y) > (p.x*p.x+p.y*p.y))
			return 1;
		if (theta() < p.theta())
			return -1;
		else
			return 1;
	};
	inline bool operator == (const point2d & p) const {
		return (x == p.x && y == p.y ? true : false);
	};
	inline bool operator != (const point2d & p) const {
		return (x != p.x || y != p.y ? true : false);
	};
	inline bool operator > (const point2d & p) const {
		return (compare(p) == 1 ? true : false);
	};
	inline bool operator >= (const point2d & p) const {
		return (compare(p) != -1 ? true : false);
	};
	inline bool operator < (const point2d & p) const {
		return (compare(p) == -1 ? true : false);
	};
	inline bool operator <= (const point2d & p) const {
		return (compare(p) != 1 ? true : false);
	};
};

/*
 * A simple point in 3D space.
 */
struct point3d : public point2d {
	// WARNING: Z is positive down.
	double z;

	point3d() : point2d() {
	};
	point3d(double px, double py, double pz) : point2d(px,py) {
		z = pz;
	};
	point3d(const point3d & p) : point2d(p) {
		z = p.z;
	};
	~point3d () {
	};
	int compare(const point3d & p) const {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if ((x*x+y*y) < (p.x*p.x+p.y*p.y))
			return -1;
		if ((x*x+y*y+z*z) > (p.x*p.x+p.y*p.y+p.z*p.z))
			return 1;
		return point2d::compare((const point2d &)p);
	};
	bool operator == (const point3d & p) const {
		return (x == p.x && y == p.y && z == p.z ? true : false);
	};
	bool operator != (const point3d & p) const {
		return (x != p.x || y != p.y || z != p.z ? true : false);
	};
	inline bool operator > (const point3d & p) const {
		return (compare(p) == 1 ? true : false);
	};
	inline bool operator >= (const point3d & p) const {
		return (compare(p) != -1 ? true : false);
	};
	inline bool operator < (const point3d & p) const {
		return (compare(p) == -1 ? true : false);
	};
	inline bool operator <= (const point3d & p) const {
		return (compare(p) != 1 ? true : false);
	};
};

/*
 * struct LElayer - A linear elastic pavement layer.
 *
 * This implements the basics structure of a pavement layer.
 */
struct LElayer : public listelement_o<LEsystem,LElayer> {
public:
	double bottom() const {
		return (h == 0.0 ? 0.0 : top() + h);
	};
	double top() const {
		return (prev == 0 ? 0.0 : prev->bottom());
	};
	double thickness() const {
		return h;
	};
	double thickness(const double th) {
		return h = th;
	};
	double emod() const {
		return E;
	};
	double emod(const double te) {
		return E = te;
	};
	double poissons() const {
		return v;
	};
	double poissons(const double tv) {
		return v = tv;
	};
	double slip() const {
		return s;
	};
	double slip(const double ts) {
		return s = ts;
	};
	LElayer(LEsystem * o, LElayer * p)
		: listelement_o<LEsystem,LElayer>(o,p) {
		h = E = v = 0.0;
		s = 1.0;
	};
	LElayer(LEsystem * o, LElayer * p,
			const double th, const double te, const double tv, const double ts = 1.0)
		: listelement_o<LEsystem,LElayer>(o,p) {
		h = th;
		E = te;
		v = tv;
		s = ts;
	};
	LElayer(LEsystem * o, LElayer * p, const LElayer & pl)
		: listelement_o<LEsystem,LElayer>(o,p) {
		h = pl.h;
		E = pl.E;
		v = pl.v;
		s = pl.s;
	};
	virtual ~LElayer() {
	};
private:
	friend class LEsystem;
	double h;							// The thickness of the layer.
	double E;							// The elastic modulus.
	double v;							// Poisson's ratio.
	double s;							// Layer slip (on bottom boundary).
};

/*
 * class paveload - A load on the pavement structure.
 */
struct paveload : point2d {
	double inline force() const {
		return f;
	};
	double inline pressure() const {
		return p;
	};
	double inline radius() const {
		return (p == 0.0 ? 0.0 : sqrt(f/p/M_PI));
	};
	double inline area() const {
		return f/p;
	};
	paveload() : point2d() {
		x = y = f = p = 0.0;
	};
	paveload(const point2d & l, double lf, double lp, double lr = 0)
		: point2d(l) {
		if (lr == 0.0) {
			f = lf;
			p = lp;
		} else if (lf == 0.0) {
			p = lp;
			f = M_PI*lr*lr*lp;
		} else {
			f = lf;
			p = f/M_PI/lr/lr;
		}
	};
	paveload(const paveload & pl) : point2d(pl) {
		f = pl.f;
		p = pl.p;
	};
	virtual ~paveload() {
	};
private:
	friend class LEsystem;
	double f;
	double p;
};

/*
 * class pavedata - An evaluation point within the pavement structure.
 */
struct pavedata : point3d {
	enum type {deflct, stress, strain};
	enum direction {xx, yy, zz, xy, xz, yz, p1, p2, p3, s1, s2, s3};
	double result(type t, direction d) const;

	pavedata() : point3d() {
	};
	pavedata(const point3d & p) : point3d(p) {
		memset(data,0,sizeof(data));
	};
	pavedata(const pavedata & pd) : point3d(pd), deflgrad(pd.deflgrad) {
		memcpy(data, pd.data, sizeof(data));
	};
	virtual ~pavedata() {
	};
//private:
	double data[9][3];
	sset<double> deflgrad;
	friend class LEsystem;
	friend class LEbackcalc;
	void principle(double v, double E);
};

/*
 * class LEsystem - A layered elastic pavement system.
 */
class LEsystem : private list_owned<LEsystem, LElayer> {
public:
	bool addlayer(const double h, const double e, const double v,
	              const double s = 1.0, const int p = -1);
	bool removelayer(const int l);
	bool removelayers();
	bool addload(const point2d & l, double f, double p, double r = 0) {
		return load.add(paveload(l,f,p,r));
	};
	bool addload(const paveload & l) {
		return load.add(l);
	};
	bool removeload(const int i);
	bool removeloads();
	bool addpoint(const point3d & p);
	bool addgrid(const int nx, const double *xp,
				 const int ny, const double *yp,
				 const int nz, const double *zp);
	bool removepoints();
	bool removepoint(const point3d & p);

	bool check();
	enum resulttype {all, fast, disp, fastdisp, dispgrad, fastgrad};
	bool accurate();
	bool calculate(resulttype result = all, double * Q = 0);
	bool odemark();
	bool fastnum();

	const pavedata & result(const point3d & p) const {
		return data[p];
	};
	const int inline layers() const {
		return length();
	};
	LElayer & layer(const int l);

	LEsystem()
		: list_owned<LEsystem,LElayer>(), data(), load() {
		callcount = 0;
	};
	LEsystem::LEsystem(LEsystem & p)
		: list_owned<LEsystem,LElayer>(p), data(p.data), load(p.load) {
		callcount = 0;
	};
	~LEsystem() {
	};
//private:
	int callcount;
	ksset<point3d,pavedata> data;
	sset<paveload> load;
	friend class listelement_o<LEsystem, LElayer>;
	friend class list_owned<LEsystem, LElayer>;
	friend class LEbackcalc;
};

/*
 * class defldata - A measured point within the pavement structure.
 */
struct defldata : point3d {
	defldata() : point3d() {
	};
	defldata(const point3d & p, const double m) : point3d(p) {
		measured = m;
		calculated = 0.0;
	};
	defldata(const defldata & dd) : point3d(dd) {
		measured = dd.measured;
		calculated = dd.calculated;
	};
	virtual ~defldata() {
	};
//private:
	double measured;
	double calculated;
	friend class LEbackcalc;
};

/*
 * class LEbackcalc - A layered elastic backcalculation.
 */
class LEbackcalc : public LEsystem {
public:
	bool adddefl(const point3d & p, double d) {
		return defl.add(defldata(p,d));
	};
	bool adddefl(const defldata & d) {
		return defl.add(d);
	};
	bool removedeflections() {
		return defl.empty();
	};
	void setup(double p, double n, double t, int m) {
		precision = MAX(0.0,p);
		noise = MAX(0.0,n);
		tolerance = MAX(1e-6,t);
		maxsteps = MAX(3,m);
	};
	bool backcalc();

	LEbackcalc() {
		setup(0.0,0.0,1e-6,5);
	};
	LEbackcalc::LEbackcalc(LEbackcalc & b)
		: LEsystem(b), defl(b.defl) {
		setup(b.precision,b.noise,b.tolerance,b.maxsteps);
	};
	~LEbackcalc() {
	};
//private:
	sset<defldata> defl;
	double precision;
	double noise;
	double tolerance;
	int	maxsteps;

	enum calctype {slow, fast, reuse};
	bool seed(int nl, double * P);
	double deflgrad(int nl, double * P, double * Q, calctype cl = slow);
	double gaussnewton(int nl, double * P, calctype cl = slow);
	double kalman(int nl, double * P);
	double bowlerror(int nl = 0, double * P = 0,
					 double s = 0.0, double * D = 0);
	double brent(int nl, double * P, double *D);
	double conjgrad(int nl, double * P);
};

#endif // PAVEMENT_H
