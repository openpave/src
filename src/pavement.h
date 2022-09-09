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

	Portions Copyright (C) 2006-2022 OpenPave.org.

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

#pragma once
#ifndef __PAVEMENT_H
#define __PAVEMENT_H

#include <memory.h>
#include <stdexcept>
#include "../include/cache.h"
#include "../include/mathplus.h"
#include "../include/set.h"
#include "../include/list.h"

namespace OP {

class LEsystem;
class LEbackcalc;

/*
 * struct point2d - A simple point in 2D space.
 */
struct point2d {
	double x, y;

	point2d() noexcept {
	}
	point2d(double px, double py) noexcept
	  : x(px), y(py) {
	}
	point2d(const point2d &) noexcept = default;
	point2d & operator = (const point2d &) noexcept = default;
	double r() const noexcept {
		return hypot(x,y);
	}
	double theta() const noexcept {
		return atan2(x,y);
	}
	double distance(const point2d & p) const noexcept {
		return hypot(x-p.x,y-p.y);
	}
	int compare(const point2d & p) const noexcept {
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
	bool operator == (const point2d & p) const noexcept {
		return (x == p.x && y == p.y ? true : false);
	}
	bool operator != (const point2d & p) const noexcept {
		return (x != p.x || y != p.y ? true : false);
	}
	bool operator > (const point2d & p) const noexcept {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const point2d & p) const noexcept {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const point2d & p) const noexcept {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const point2d & p) const noexcept {
		return (compare(p) != 1 ? true : false);
	}
};

/*
 * struct point3d - A simple point in 3D space.
 */
struct point3d : public point2d {
	// WARNING: Z is positive down.
	double z;

	point3d() noexcept
	  : point2d() {
	}
	point3d(double px, double py, double pz) noexcept
	  : point2d(px,py), z(pz) {
	}
	point3d(const point3d &) noexcept = default;
	point3d & operator = (const point3d &) noexcept = default;
	int compare(const point3d & p) const noexcept {
		if (x == p.x && y == p.y && z == p.z)
			return 0;
		if ((x*x+y*y) < (p.x*p.x+p.y*p.y))
			return -1;
		if ((x*x+y*y+z*z) > (p.x*p.x+p.y*p.y+p.z*p.z))
			return 1;
		return point2d::compare(static_cast<const point2d &>(p));
	}
	bool operator == (const point3d & p) const noexcept {
		return (x == p.x && y == p.y && z == p.z ? true : false);
	}
	bool operator != (const point3d & p) const noexcept {
		return (x != p.x || y != p.y || z != p.z ? true : false);
	}
	bool operator > (const point3d & p) const noexcept {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const point3d & p) const noexcept {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const point3d & p) const noexcept {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const point3d & p) const noexcept {
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
		return (prev == nullptr ? 0.0 : prev->bottom());
	}
	double thickness() const noexcept {
		return h;
	}
	double emod() const noexcept {
		return E;
	}
	double poissons() const noexcept {
		return v;
	}
	double slip() const noexcept {
		return s;
	}
	double thickness(double th) noexcept;
	double emod(double te) noexcept;
	double poissons(double tv) noexcept;
	double slip(double ts) noexcept;
	LElayer(LEsystem * o, LElayer * p, LElayer * n = nullptr)
	  : listelement_o<LEsystem,LElayer>(o,p,n), h(0.0), E(0.0), v(0.0),
			s(1.0) {
	}
	LElayer(LEsystem * o, LElayer * p, LElayer * n, double th, double te,
			double tv, double ts = 1.0)
	  : listelement_o<LEsystem,LElayer>(o,p,n), h(th), E(te), v(tv),
			s(ts) {
	}
	LElayer(LEsystem * o, LElayer * p, LElayer * n, const LElayer & l)
	  : listelement_o<LEsystem,LElayer>(o,p,n), h(l.h), E(l.E), v(l.v),
			s(l.s) {
	}

private:
	friend class LEsystem;

	double h;                  // The thickness of the layer.
	double E;                  // The elastic modulus.
	double v;                  // Poisson's ratio.
	double s;                  // Layer slip (on bottom boundary).
};

/*
 * struct paveload - A load on the pavement structure.
 */
class paveload : public point2d {
public:
	paveload() noexcept
	  : point2d(0.0,0.0), f(0.0), r(0.0) {
	}
	paveload(const point2d & l, double lf, double lp, double lr = 0.0) noexcept
	  : point2d(l), f(lf), r(lr)  {
		if (lf == 0.0) {
			f = M_PI*lr*lr*lp;
		} else if (lr == 0.0) {
			r = sqrt(lf/(lp*M_PI));
		}
	}
	paveload(const paveload & pl) noexcept
	  : point2d(pl), f(pl.f), r(pl.r) {
	}
	double force() const noexcept {
		return f;
	}
	double pressure() const noexcept {
		return (r == 0.0 ? 0.0 : f/(M_PI*r*r));
	}
	double radius() const noexcept {
		return r;
	}
	double area() const noexcept {
		return M_PI*r*r;
	}

private:
	double f;                  // Force
	double r;                  // Load radius
};

/*
 * struct pavepoint - A point in a pavement (3D and layer)
 */
struct pavepoint : public point3d {
	unsigned il;

	pavepoint() noexcept
	  : point3d(0.0,0.0,0.0), il(UINT_MAX) {
	}
	pavepoint(double px, double py, double pz, unsigned l = UINT_MAX) noexcept
	  : point3d(px,py,pz), il(l) {
	}
	pavepoint(const point3d & p, unsigned l = UINT_MAX) noexcept
	  : point3d(p), il(l) {
	}
	pavepoint(const pavepoint & p) noexcept
	  : point3d(p), il(p.il) {
	}
	pavepoint & operator = (const pavepoint & p) noexcept {
		x = p.x; y = p.y; z = p.z;
		if (il == UINT_MAX || p.il != UINT_MAX)
			il = p.il;
		return *this;
	}
	int compare(const pavepoint & p) const noexcept {
		if (x == p.x && y == p.y && z == p.z
		  && (il == p.il || il == UINT_MAX || p.il == UINT_MAX))
			return 0;
		if (il != UINT_MAX && p.il != UINT_MAX && il < p.il)
			return -1;
		if (il != UINT_MAX && p.il != UINT_MAX && il > p.il)
			return 1;
		return point3d::compare(static_cast<const point3d &>(p));
	}
	bool operator == (const pavepoint & p) const noexcept {
		return (x == p.x && y == p.y && z == p.z
				 && (il == p.il || il == UINT_MAX || p.il == UINT_MAX)
				 ? true : false);
	}
	bool operator != (const pavepoint & p) const noexcept {
		return (x != p.x || y != p.y || z != p.z
				 || (il != p.il && il != UINT_MAX && p.il != UINT_MAX)
				 ? true : false);
	}
	bool operator > (const pavepoint & p) const noexcept {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const pavepoint & p) const noexcept {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const pavepoint & p) const noexcept {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const pavepoint & p) const noexcept {
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
	double result(type t, direction d) const noexcept {
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
				return 0.0; // XXX: throw?
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
#if defined(_MSC_VER)
		default:
			return 0.0;
#endif
		}
	}

	pavedata()
	  : pavepoint(), deflgrad(0) {
	}
	pavedata(const point3d & p, unsigned l = UINT_MAX)
	  : pavepoint(p,l), deflgrad(0) {
	}
	pavedata(const pavepoint & p)
	  : pavepoint(p), deflgrad(0) {
	}
	pavedata(const pavedata & pd)
	  : pavepoint(pd), deflgrad(pd.deflgrad) {
		memcpy(data,pd.data,sizeof(data));
	}
	pavedata & operator = (const pavedata & pd) {
		pavepoint::operator=(pd);
		deflgrad = pd.deflgrad;
		memcpy(data,pd.data,sizeof(data));
		return *this;
	}
	~pavedata() {
	}

	double data[9][3];
	fset<double> deflgrad;
	void principle(double v, double E) noexcept;
};

/*
 * class LEsystem_cache - A wrapper for the slab_cache.
 */
template<typename S>
class LEsystem_cache {
protected:
	LEsystem_cache() noexcept
	  : cache_state(S::empty), cache(nullptr) {
	}
	LEsystem_cache(const LEsystem_cache&) = delete;
	LEsystem_cache& operator = (const LEsystem_cache&) = delete;
	LEsystem_cache(LEsystem_cache&&) noexcept
	  : cache_state(S::empty), cache(nullptr) {
	}
	~LEsystem_cache() {
		cache_free();
	}

	template<typename T> T* cache_alloc(unsigned count) {
		if (cache == nullptr)
			cache = slab_cache::create();
		return static_cast<T*>(cache->allocate(count * sizeof(T)));
	}
	void cache_reset() noexcept {
		if (cache != nullptr)
			cache->reset();
	}
	void cache_free() noexcept {
		cached_state(S::empty);
		slab_cache::destroy(cache);
		cache = nullptr;
	}
	void cached_state(S s) noexcept {
		cache_state = MIN(cache_state, s);
	}

	S cache_state;
	slab_cache* cache;
};

// declare externally to avoid ordering issues
enum class LEsystem_cachestate {
	empty, emod, all
};

/*
 * class LEsystem - A layered elastic pavement system.
 */
class LEsystem
  : private list_owned<LEsystem, LElayer>,
	private LEsystem_cache<LEsystem_cachestate>
{
public:
	LEsystem()
	  : list_owned<LEsystem,LElayer>(),
		LEsystem_cache<LEsystem_cachestate>(),
		points(), data(), lg(), clg(0) {
	}
	LEsystem(const LEsystem &) = delete;
	LEsystem & operator = (const LEsystem &) = delete;
	LEsystem(LEsystem && s)
	  : list_owned<LEsystem,LElayer>(std::move(s)),
		LEsystem_cache<LEsystem_cachestate>(std::move(s)),
		points(std::move(s.points)), data(std::move(s.data)),
		lg(std::move(s.lg)), clg(s.clg) {
	}
	~LEsystem() {
		cache_free();
	}
	LElayer & addlayer(const double h, const double e, const double v,
	                   const double s = 1.0, const unsigned p = UINT_MAX);
	void removelayer(const unsigned l) noexcept;
	void removelayers() noexcept;
	unsigned layers() const noexcept {
		return length();
	}
	LElayer & layer(const unsigned l) const noexcept;

	unsigned defaultgroup() const noexcept {
		return clg;
	}
	unsigned defaultgroup(unsigned g) noexcept {
		return clg = g;
	}
	unsigned groups() const noexcept {
		return lg.length();
	}
	void addload(const point2d & l, double f, double p, double r = 0) {
		addload(clg,l,f,p,r);
	}
	void addload(const paveload & l) {
		addload(clg,l);
	}
	void removeload(const unsigned i) {
		removeload(clg,i);
	}
	void removeloads() {
		removeloads(clg);
	}
	unsigned loads() const {
		if (lg.length() == 0)
			return 0;
		return loads(clg);
	}
	const paveload & getload(const unsigned i) {
		return getload(clg,i);
	}
	void addload(const unsigned g, const point2d & l, double f, double p,
			double r = 0) {
		addload(g,paveload(l,f,p,r));
	}
	void addload(const unsigned g, const paveload & l);
	void removeload(const unsigned g, const unsigned i);
	void removeloads(const unsigned g);
	unsigned loads(const unsigned g) const
	{
		if (!lg.inbounds(g))
			throw std::out_of_range("Load group out of range!");
		return lg[g].length();
	}
	const paveload & getload(const unsigned g, const unsigned i)
	{
		if (!lg.inbounds(g))
			throw std::out_of_range("Load group out of range!");
		return lg[g][i];
	}

	void addpoint(const point3d & p, unsigned l = UINT_MAX);
	void addgrid(const unsigned nx, const double * xp,
				 const unsigned ny, const double * yp,
				 const unsigned nz, const double * zp);
	void removepoint(const point3d & p, unsigned l = UINT_MAX);
	void removepoints() noexcept;
	unsigned results() const noexcept {
		if (cache_state < cachestate::all)
			return 0;
		return points.length();
	}
	const pavedata & result(const point3d & p, unsigned l = UINT_MAX) const {
		return result(clg,p,l);
	}
	const pavedata & result(const unsigned i) const {
		return result(clg,i);
	}
	const pavedata & result(const unsigned g, const point3d & p,
		unsigned l = UINT_MAX) const {
		if (cache_state < cachestate::all)
			throw std::runtime_error("No results available!");
		return result(g,points.haskey(pavepoint(p,l)));
	}
	const pavedata & result(const unsigned g, const unsigned i) const {
		if (cache_state < cachestate::all)
			throw std::runtime_error("No results available!");
		return data[g][i];
	}

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
		fastgrad  = 0x0301,
		failure   = 0xFFFF
	};
	bool calc_accurate();
	bool calculate(resulttype result = all, const double * Q = nullptr);
	bool calc_odemark();
	bool calc_fastnum();

protected:
	bool check();

private:
	friend class listelement_o<LEsystem,LElayer>;
	friend class list_owned<LEsystem,LElayer>;
	friend class LElayer;
	friend class LEbackcalc;

	ksset<pavepoint,pavepoint> points;
	sset<sset<pavedata>> data;
	sset<sset<paveload>> lg;
	unsigned clg;

	resulttype cache_res = failure;
	using cachestate = LEsystem_cachestate;
	using LEsystem_cache<LEsystem_cachestate>::cache_alloc;
	using LEsystem_cache<LEsystem_cachestate>::cache_reset;
	using LEsystem_cache<LEsystem_cachestate>::cache_free;
	using LEsystem_cache<LEsystem_cachestate>::cached_state;
	using LEsystem_cache<LEsystem_cachestate>::cache_state;
};

/*
 * class defldata - A measured point within the pavement structure.
 */
struct defldata : public point3d {
	defldata() noexcept
	  : point3d() {
	}
	defldata(const point3d & p, const double m) noexcept
	  : point3d(p) {
		measured = m;
		calculated = 0.0;
	}
	defldata(const defldata & dd) noexcept
	  : point3d(dd) {
		measured = dd.measured;
		calculated = dd.calculated;
	}
	~defldata() {
	}
	double measured;
	double calculated;
};

enum class LEbackcalc_cachestate {
	empty, deflgrad, kalman, conjgrad, gausssnewton, swarm
};

/*
 * class LEbackcalc - A layered elastic backcalculation.
 */
class LEbackcalc
  : private LEsystem_cache<LEbackcalc_cachestate>,
	public LEsystem
{
public:
	void adddefl(const point3d & p, double d) {
		defl.add(defldata(p,d));
	}
	void adddefl(const defldata & d) {
		defl.add(d);
	}
	unsigned deflections() const noexcept {
		return defl.length();
	}
	const defldata & getdefl(const unsigned i) const noexcept {
		return defl[i];
	}
	void removedeflections() noexcept {
		defl.empty();
	}
	void setup(double p, double n, double t, unsigned m) noexcept {
		precision = std::max(0.0,p);
		noise = std::max(0.0,n);
		tolerance = std::max(1e-6,t);
		maxsteps = std::max(3u,m);
	}
	bool backcalc();

	LEbackcalc() noexcept {
		setup(0.0,0.0,1e-6,5);
	}
	LEbackcalc(const LEbackcalc &) = delete;
	LEbackcalc & operator = (const LEbackcalc &) = delete;
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
	double bowlerror(unsigned nl = 0, const double * P = nullptr,
			const double s = 0.0, const double * D = nullptr);
	double brent(unsigned nl, double * P, double * D);
	double conjgrad(unsigned nl, double * P);
	double swarm(unsigned nl, double * P);

	using cachestate = LEbackcalc_cachestate;
	using LEsystem_cache<LEbackcalc_cachestate>::cache_alloc;
	using LEsystem_cache<LEbackcalc_cachestate>::cache_reset;
	using LEsystem_cache<LEbackcalc_cachestate>::cache_free;
	using LEsystem_cache<LEbackcalc_cachestate>::cached_state;
	using LEsystem_cache<LEbackcalc_cachestate>::cache_state;
};

} // namespace OP

#endif // PAVEMENT_H
