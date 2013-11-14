/*************************************************************************

	PAVEMENT.CPP - Implementation for the pavement class.

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

	See PAVEMENT.H.

	History:
		2002/01/23 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#include "pavement.h"
#include "matrix.h"
#include <time.h>
#include <stdio.h>

/*
 * Calculate the principle stresses, and the strains.
 */
void
pavedata::principle(double v, double E)
{
	unsigned i, j;
	double t1, t2;

	for (i = 0; i < 3; i++) {
		data[2][i] = data[0][i];
		data[3][i] = data[1][i];
	}
	// Find the eigen values of the stress tensor.  This is heavily
	// optimised for a symmetric 3x3 matrix, to use as few FP registers
	// as possible.
	t1 = data[3][0]*data[3][0] + data[3][1]*data[3][1]
			+ data[3][2]*data[3][2];
	while (t1 > 1e-30) {
		if (data[3][1] != 0.0) {
			t1 = data[2][2]-data[2][0];
			t2 = (t1 < 0.0 ? -2 : 2)*data[3][1];
			t1 = t2/(fabs(t1) + hypot(t1,t2));
			data[2][0] -= t1*data[3][1];
			data[2][2] += t1*data[3][1];
			t2 = 1.0/hypot(t1,1.0);
			data[3][2] *= t2;
			t2 = data[3][0] *= t2;
			data[3][0] -= t1*data[3][2];
			data[3][2] += t1*t2;
			data[3][1] = 0.0;
		}
		if (data[3][0] != 0.0) {
			t1 = data[2][1]-data[2][0];
			t2 = (t1 < 0.0 ? -2 : 2)*data[3][0];
			t1 = t2/(fabs(t1) + hypot(t1,t2));
			data[2][0] -= t1*data[3][0];
			data[2][1] += t1*data[3][0];
			data[3][2] /= hypot(t1,1.0);
			data[3][1] = -t1*data[3][2];
			data[3][0] = 0.0;
		}
		if (data[3][2] != 0.0) {
			t1 = data[2][2]-data[2][1];
			t2 = (t1 < 0.0 ? -2 : 2)*data[3][2];
			t1 = t2/(fabs(t1) + hypot(t1,t2));
			data[2][1] -= t1*data[3][2];
			data[2][2] += t1*data[3][2];
			data[3][1] /= hypot(t1,1.0);
			data[3][0] = -t1*data[3][1];
			data[3][2] = 0.0;
		}
		t1 = data[3][0]*data[3][0] + data[3][1]*data[3][1];
	}
	data[3][0] = fabs(data[2][0]-data[2][2])*0.5;
	data[3][1] = fabs(data[2][0]-data[2][1])*0.5;
	data[3][2] = fabs(data[2][1]-data[2][2])*0.5;
	for (i = 0; i < 2; i++) {
		for (j = i + 1; j < 3; j++) {
			if (data[2][i] > data[2][j])
				swap(data[2][i],data[2][j]);
			if (data[3][i] > data[3][j])
				swap(data[3][i],data[3][j]);
		}
	}

	t1 = 2*(v+1)/E; t2 = 1.0/E;
	data[5][0] = (data[0][0]-v*(data[0][1]+data[0][2]))*t2;
	data[5][1] = (data[0][1]-v*(data[0][2]+data[0][0]))*t2;
	data[5][2] = (data[0][2]-v*(data[0][0]+data[0][1]))*t2;
	data[6][0] = t1*data[1][0];
	data[6][1] = t1*data[1][1];
	data[6][2] = t1*data[1][2];

	data[7][0] = (data[2][0]-v*(data[2][1]+data[2][2]))*t2;
	data[7][1] = (data[2][1]-v*(data[2][2]+data[2][0]))*t2;
	data[7][2] = (data[2][2]-v*(data[2][0]+data[2][1]))*t2;
	data[8][0] = t1*data[3][0];
	data[8][1] = t1*data[3][1];
	data[8][2] = t1*data[3][2];
	for (i = 0; i < 2; i++) {
		for (j = i + 1; j < 3; j++) {
			if (data[7][i] > data[7][j])
				swap(data[7][i],data[7][j]);
			if (data[8][i] > data[8][j])
				swap(data[8][i],data[8][j]);
		}
	}
}

void
LEsystem::addlayer(double h, double e, const double v, const double s,
		const unsigned p)
{
	new LElayer(this,(p == UINT_MAX ? last : &layer(p)),h,e,v,s);
}

bool
LEsystem::removelayer(const unsigned l)
{
	LElayer * pl = first;
	unsigned i = 0;

	while (i++ < l && pl != 0)
		pl = pl->next;
	if (pl != 0) {
		delete pl;
		return true;
	}
	return false;
}

bool
LEsystem::addgrid(const unsigned nx, const double * xp,
                  const unsigned ny, const double * yp,
                  const unsigned nz, const double * zp)
{
	bool rv = true;
	for (unsigned ix = 0; ix < nx; ix++) {
		for (unsigned iy = 0; iy < ny; iy++) {
			for (unsigned iz = 0; iz < nz; iz++)
				rv = rv && data.add(pavedata(point3d(xp[ix],yp[iy],zp[iz])));
		}
	}
	return rv;
}

LElayer &
LEsystem::layer(const unsigned l) const
{
	LElayer * pl = first;
	unsigned i = 0;

	while (i++ < l && pl->next != 0)
		pl = pl->next;
	return *pl;
}

/*
 * This checks the structure of the pavement to ensure that is is
 * good...
 */
bool
LEsystem::check()
{
	unsigned il, nl = layers();
	const LElayer * pl;
	bool rv = true;

	if (nl == 0) {
		event_msg(EVENT_WARN,
			"Cannot calculate a pavement without any layers!");
		rv = false;
	}
	if (load.length() == 0) {
		event_msg(EVENT_WARN,
			"Cannot calculate a pavement without any loads!");
		rv = false;
	}
	if (data.length() == 0) {
		event_msg(EVENT_WARN,
			"Cannot calculate a pavement without any evaluation points!");
		rv = false;
	}
	for (il = 1, pl = first; pl != 0; il++, pl = pl->next) {
		if (pl->emod() <= 0.0) {
			event_msg(EVENT_WARN,
				"Error: Elastic modulus of layer %d must be greater"
				" than zero not %f!", il, pl->emod());
			rv = false;
		}
		if (pl->poissons() <= 0.0 || pl->poissons() > 0.5) {
			event_msg(EVENT_WARN,
				"Error: Poisson's ratio of layer %d must be between"
				" zero and one half not %f!", il, pl->poissons());
			rv = false;
		}
		if (pl->thickness() < 0.0) {
			event_msg(EVENT_WARN,
				"Error: Layer %d cannot have negative thickness!", il);
			rv = false;
		}
		if (pl->slip() < 0.0 || pl->slip() > 1.0 ) {
			event_msg(EVENT_WARN,
				"Error: Layer %d has an invalid bonding coefficent!", il);
			rv = false;
		}
		if (pl->thickness() == 0.0 && pl->next != 0) {
			event_msg(EVENT_WARN,
				"Error: Layer %d cannot have zero thickness!", il);
			rv = false;
		}
		if (pl->thickness() == 0.0 && pl->next == 0 && pl->slip() != 1.0) {
			event_msg(EVENT_WARN,
				"Error: Infinite layer %d cannot have imperfect bonding!", il);
			rv = false;
		}
	}
	if (last && last->bottom() > 0.0 && last->poissons() == 0.75) {
 		event_msg(EVENT_WARN,
			"Error: Last layer cannot have a Poisson's ratio of 0.75!");
		rv = false;
	}
	if (rv == false)
		return rv;
	for (unsigned ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		// Clear the data.
		memset(d.data,0,sizeof(d.data));
		// Set the il variables correctly.
		double z = first->top(), h = 0.0;
		if (d.z < z) {
 			event_msg(EVENT_WARN,
				"Error: evaluation point %d (%f,%f,%f) not within pavement!",
				ixy+1,d.x,d.y,d.z);
			rv = false;
			continue;
		}
		for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
			h = pl->thickness();
			if (d.il != UINT_MAX && d.il == il) {
				break;
			} else if (d.il == UINT_MAX
					&& d.z >= z && (h == 0.0 || d.z < (z+h))) {
				d.il = il;
				break;
			}
			z += h;
		}
		assert(pl == 0 || d.il != UINT_MAX);
		// Points must be within the layer or on the boundary.
		if (pl == 0 || d.z < z || (h > 0.0 && d.z > (z+h))) {
 			event_msg(EVENT_WARN,
				"Error: evaluation point %d (%f,%f,%f) not within layer %d!",
				ixy+1,d.x,d.y,d.z,d.il+1);
			rv = false;
		}
	}
	return rv;
}

/*
 * As simple struct to store the axisymetric data.  These are
 * stored in the ax set, indexed by the radius and depth.
 */
struct axialdata {
	bool   active;
	double rse;
	double tse;
	double vse;
	double sse;
	double rdp;
	double vdp;
	double rse2;
	double tse2;
	double vse2;
	double sse2;
	double rdp2;
	double vdp2;
	
	// This function takes the results from the xxx2 variables, and adds
	// them to the main varaibles, reseting them, and then deactivating
	// the point if none changed by more than epsilon.
	void accumulate(LEsystem::resulttype res, const double & eps) {
		bool isactive = false;
		
		if (!active)
			return;
		if (!(res & LEsystem::disp)) {
			vse += vse2;
			isactive = isactive || (fabs(vse2) > eps*fabs(vse));
			rse += rse2;
			isactive = isactive || (fabs(rse2) > eps*fabs(rse));
			tse += tse2;
			isactive = isactive || (fabs(tse2) > eps*fabs(tse));
			sse += sse2;
			isactive = isactive || (fabs(sse2) > eps*fabs(sse));
			rdp += rdp2;
			isactive = isactive || (fabs(rdp2) > eps*fabs(rdp));
			vse2 = 0.0; rse2 = 0.0; tse2 = 0.0;
			sse2 = 0.0; rdp2 = 0.0;
		}
		vdp += vdp2;
		isactive = isactive || (fabs(vdp2) > eps*fabs(vdp));
		vdp2 = 0.0;
		if (!isactive && eps > 0.0)
			active = false;
	} 
	// This function does the final step of the axial data calculations,
	// which is to multiply each result by a constant and fixup the radial
	// and tangential shear stresses.
	void finalize(LEsystem::resulttype res, const double & r,
			const double & a, const double & v, const double & E) {
		double t1 = a*(v+1)/E, t2;
 
		// Always make sure we accumulate the results.
		accumulate(res,0.0);
		vdp *= t1;
		if (res & LEsystem::disp)
			return;
		if (r > 0.0) {
			t2 = rdp/r; rse -= t2; tse += t2;
		}
		rdp *= t1; rse *= a; tse *= a; vse *= a; sse *= a;
	}
	// This function takes the data from this axial data point and
	// adds it to the general data point, accounting for rotation.
	void addtodata(LEsystem::resulttype res, pavedata * d,
			const paveload & l, const double & r,
			unsigned gl = UINT_MAX) const {
		double p = l.pressure();
		
		assert(fabs(r-l.distance(*d)) < DBL_MIN);
		if (gl != UINT_MAX) {
			d->deflgrad[gl] += p*vdp;
			return;
		}
		d->data[4][2] += p*vdp;
		if (res & LEsystem::disp)
			return;
		d->data[0][2] += p*vse;
		if (r == 0.0) {
			d->data[0][0] += p*(rse+tse)/2;
			d->data[0][1] += p*(rse+tse)/2;
			d->data[1][1] += p*sse;
		} else {
			double sint = (d->y-l.y)/r;
			double cost = (d->x-l.x)/r;
			if (cost == 0.0) {
				d->data[0][0] += p*tse;
				d->data[0][1] += p*rse;
				d->data[1][2] += p*SGN(sint)*sse;
				d->data[4][1] += p*SGN(sint)*rdp;
			} else if (sint == 0.0) {
				d->data[0][0] += p*rse;
				d->data[0][1] += p*tse;
				d->data[1][1] += p*SGN(cost)*sse;
				d->data[4][0] += p*SGN(cost)*rdp;
			} else {
				d->data[0][0] += p*(cost*cost*rse+sint*sint*tse);
				d->data[0][1] += p*(cost*cost*tse+sint*sint*rse);
				d->data[1][0] += p*cost*sint*(rse-tse);
				d->data[1][1] += p*cost*sse;
				d->data[1][2] += p*sint*sse;
				d->data[4][0] += p*cost*rdp;
				d->data[4][1] += p*sint*rdp;
			}
		}
	}
};

/*
 * Store the depth and layer.
 */
struct zpoint {
	double z;
	unsigned il;
	
	zpoint() {
	}
	zpoint(double d, unsigned l)
	  : z(d), il(l) {
	}
	int compare(const zpoint & p) const {
		if (z == p.z && il == p.il)
			return 0;
		if (z < p.z)
			return -1;
		if (z > p.z)
			return 1;
		if (il < p.il)
			return -1;
		else
			return 1;
	}
	bool operator == (const zpoint & p) const {
		return (z == p.z && il == p.il ? true : false);
	}
	bool operator != (const zpoint & p) const {
		return (z != p.z || il != p.il ? true : false);
	}
	bool operator > (const zpoint & p) const {
		return (compare(p) == 1 ? true : false);
	}
	bool operator >= (const zpoint & p) const {
		return (compare(p) != -1 ? true : false);
	}
	bool operator < (const zpoint & p) const {
		return (compare(p) == -1 ? true : false);
	}
	bool operator <= (const zpoint & p) const {
		return (compare(p) != 1 ? true : false);
	}
};

#define GRADSTEP 1e-8

#define NBZ		8192						// Number of zeros in bessels.
static double j0r[NBZ+1], j1r[NBZ+1], j0p[NBZ+1], j0m[NBZ+1], j1max;
#define j0pj1(x)	(j0(x)+j1(x))
#define j0mj1(x)	(j0(x)-j1(x))
#define j1d(x)		(j0(x)-j1(x)/x)

#define NGQP	16						// Max number of Gauss points
static double gu[NGQP+1][NGQP];
static double gf[NGQP+1][NGQP];

//#define NLQP	16						// Max number of Lobatto points
//static double lu[NLQP+1][NLQP];
//static double lf[NLQP+1][NLQP];

/*
 * This initialises the arrays above to max accuracy, and
 * find the zeros of the various bessel functions...
 */
static void
initarrays()
{
	unsigned ib;
	static bool done = false;

	if (done)
		return;
	memset(gu,0,sizeof(double)*(NGQP+1)*NGQP);
	memset(gf,0,sizeof(double)*(NGQP+1)*NGQP);
	for (ib = 1; ib <= NGQP; ib++) {
		for (unsigned j, i = 0; i < (ib+1)/2; i++) {
			double p1, p2, p3, pp, z = cos(M_PI*(i+0.75)/(ib+0.5));
			do {
				p1 = 1.0, p2 = 0.0;
				for (j = 0; j < ib; j++)
					p3 = p2, p2 = p1, p1 = z*p2+j*(z*p2-p3)/(j+1);
				pp = ib*(z*p1-p2)/(z*z-1.0), z -= p1/pp;
			} while (fabs(p1/pp) > DBL_EPSILON);
			gu[ib][i] = -z, gu[ib][ib-1-i] = z;
			gf[ib][i] =     gf[ib][ib-1-i] = 2.0/(1.0-z*z)/pp/pp;
		}
	}
	//memset(lu,0,sizeof(double)*(NLQP+1)*NLQP);
	//memset(lf,0,sizeof(double)*(NLQP+1)*NLQP);
	//for (ib = 2; ib <= NLQP; ib++) {
	//	for (unsigned j, i = 0; i < (ib+1)/2; i++) {
	//		double p1, p2, p3, z = cos(M_PI*i/(ib-1));
	//		do {
	//			p2 = 1.0, p1 = z;
	//			for (j = 1; j < (ib-1); j++)
	//				p3 = p2, p2 = p1, p1 = z*p2+j*(z*p2-p3)/(j+1);
	//			z += (p2/p1-z)/ib;
	//		} while (fabs(p2/p1-z)/ib > DBL_EPSILON);
	//		lu[ib][i] = -z, lu[ib][ib-1-i] = z;
	//		lf[ib][i] =     lf[ib][ib-1-i] = 2.0/(ib*(ib-1)*p1*p1);
	//	}
	//}
	memset(j0r,0,sizeof(double)*(NBZ+1));
	memset(j1r,0,sizeof(double)*(NBZ+1));
	memset(j0p,0,sizeof(double)*(NBZ+1));
	memset(j0m,0,sizeof(double)*(NBZ+1));
	for (ib = 0, j1max = 0.0; ib <= NBZ; ib++) {
		double x1, x2, y1, y2, xi;
		x1 = (ib == 0 ? M_PI_2 : j0r[ib-1] + M_PI), x2 = x1 + 0.1;
		while (fabs(j0(x2)) > 1e-30 && x1 != x2) {
			y1 = j0(x1), y2 = j0(x2);
			xi = x2 - y2*(x2-x1)/(y2-y1);
			x1 = x2, x2 = xi;
		}
		j0r[ib] = x2;
		x1 = x2 + M_PI_2, x2 = x1 + 0.1;
		while (fabs(j1(x2)) > 1e-30 && x1 != x2) {
			y1 = j1(x1), y2 = j1(x2);
			xi = x2 - y2*(x2-x1)/(y2-y1);
			x1 = x2, x2 = xi;
		}
		j1r[ib] = x2;
		if (fabs(j1(x2)) > j1max)
			j1max = fabs(j1(x2));
		x1 = j0r[ib] - M_PI_4, x2 = x1 + 0.1;
		while (fabs(j0mj1(x2)) > 1e-30 && x1 != x2) {
			y1 = j0mj1(x1), y2 = j0mj1(x2);
			xi = x2 - y2*(x2-x1)/(y2-y1);
			x1 = x2, x2 = xi;
		}
		j0m[ib] = x2;
		x1 = j0r[ib] + M_PI_4, x2 = x1 + 0.1;
		while (fabs(j0pj1(x2)) > 1e-30 && x1 != x2) {
			y1 = j0pj1(x1), y2 = j0pj1(x2);
			xi = x2 - y2*(x2-x1)/(y2-y1);
			x1 = x2, x2 = xi;
		}
		j0p[ib] = x2;
		x1 = j0r[ib], x2 = x1 + 0.1;
		while (fabs(j1d(x2)) > 1e-30 && x1 != x2) {
			y1 = j1d(x1), y2 = j1d(x2);
			xi = x2 - y2*(x2-x1)/(y2-y1);
			x1 = x2, x2 = xi;
		}
		j0p[ib] *= j0r[ib]/x2, j0m[ib] *= j0r[ib]/x2;
	}
	done = true;
}

// These are the asymptotic expansions of the integrals of J1(ma)*J0(mr)
// and J1(ma)*J0(mr).  They are not scaled, since we only use them to find
// roots, and the scale actually hurts us for large m values.  They
// should not be used if r+-=a.  Also, the search should be started very
// close to their true zero, since they are not stable.
#define j1j0(m,a,r) (sin((a+r)*m)/(a+r)+cos((a-r)*m)/(a-r))
#define j1j1(m,a,r) (sin((a-r)*m)/(a-r)+cos((a+r)*m)/(a+r))

static double
refine_m0(double ma, double mb, double a, double r)
{
	double m0, ya = j1j0(ma,a,r), yb = j1j0(mb,a,r), y0;
	do {
		m0 = ma - (ma-mb)*ya/(ya-yb), y0 = j1j0(m0,a,r);
		if (ma == m0 || mb == m0)
			break;
		(yb*y0 < 0 ? ma : mb) = m0, (ma == m0 ? ya : yb) = y0;
	} while (fabs(y0) > DBL_EPSILON);
	return m0;
}

static double
refine_m1(double ma, double mb, double a, double r)
{
	double m1, ya = j1j1(ma,a,r), yb = j1j1(mb,a,r), y1;
	do {
		m1 = ma - ya*(ma-mb)/(ya-yb), y1 = j1j1(m1,a,r);
		if (ma == m1 || mb == m1)
			break;
		(yb*y1 < 0 ? ma : mb) = m1, (ma == m1 ? ya : yb) = y1;
	} while (fabs(y1) > DBL_EPSILON);
	return m1;
}

/*
 * For the functions J1(ma)*J0(mr) and J1(ma)*J1(mr) find the correct
 * points to stop integrating so that we are very close to the m
 * values where the integrals are exact.
 */
static void
stoppingpoints(const unsigned nbz, const double a, const double r,
               double * m0, double * m1)
{
	unsigned ib;
	double ra = r/a, r1 = fabs(ra-1.0);

	// At r=0 or r=a just use the last root.
	*m0 = j0r[nbz]/a, *m1 = (ra == 0.0 ? 0.0 : j0r[nbz]/a);
	if (ra > 0.0 && ra <= 0.5) {
		// Find the last root of J1(mr) & J0(mr).
		for (ib = nbz; ib > 0 && j1r[ib] > j1r[nbz]*ra; ib--)
			;
		*m0 = j1r[ib]/ra;
		for (ib = nbz; ib > 0 && j0r[ib] > j1r[nbz]*ra; ib--)
			;
		*m1 = j0r[ib]/ra;
	} else if (ra > 0.5 && ra != 1.0 && ra < 2.0) {
		// Between 0.5 and 2.0 use the J0(m|r-1|)-J1(m|r-1|) and
		// J0(m|r-1|)+J1(m|r-1|) approximations.
		for (ib = nbz; ib > 0 && j0m[ib] >
				r1*(ra < 1.0 ? j1r[nbz] : j0r[nbz]/ra); ib--)
			;
		*m0 = j0m[ib]/r1;
		for (ib = nbz; ib > 0 && j0p[ib] >
				r1*(ra < 1.0 ? j1r[nbz] : j1r[nbz]/ra); ib--)
			;
		*m1 = j0p[ib]/r1;
	} else if (ra >= 2.0) {
		// Find the roots of J0(m).
		for (ib = nbz; ib > 0 && j0r[ib] > j0r[nbz]/ra; ib--)
			;
		*m0 = j0r[ib];
		for (ib = nbz; ib > 0 && j0r[ib] > j1r[nbz]/ra; ib--)
			;
		*m1 = j0r[ib];
	}
	if (ra > 0.0 && ra < 1.0) {
		// Now find the closest roots of J0(ma).
		for (ib = nbz; ib > 0 && j1r[ib] >= *m0; ib--)
			;
		*m0 = (ib == nbz ? j0r[nbz]/a
				: refine_m0(j1r[ib]/a,j1r[ib+1]/a,a,r));
		for (ib = nbz; ib > 0 && j1r[ib] >= *m1; ib--)
			;
		*m1 = (ib == nbz ? j0r[nbz]/a
				: refine_m1(j1r[ib]/a,j1r[ib+1]/a,a,r));
	} else if (ra > 1.0) {
		for (ib = nbz; ib > 0 && j0r[ib]/ra >= *m0; ib--)
			;
		*m0 = (ib == nbz ? j1r[nbz]/r
				: refine_m0(j0r[ib]/r,j0r[ib+1]/r,a,r));
		for (ib = nbz; ib > 0 && j1r[ib]/ra >= *m1; ib--)
			;
		*m1 = (ib == nbz ? j0r[nbz]/r
				: refine_m1(j1r[ib]/r,j1r[ib+1]/r,a,r));
	}
}

/*
 * This build the ABCD matrix, based on the structure.
 */
static void
buildabcd(const double m, const unsigned nl, const double * h,
          const double * v, const double * E,
          const double * f, double (* __restrict R)[4][2],
          double (* __restrict ABCD)[4])
{
	unsigned i, j, il;
	double B1[2][4], X[4][4], F[4][4], D[4][4];
	double CDi[4] = {0, 0, 0, 0};
	double mi = 1.0/m;

	memset(ABCD,0,sizeof(double)*nl*4);
	if (m <= 0.0)
		return;
	// We start at the last layer...
	il = nl-1;
	memset(&R[il][0][0],0,8*sizeof(double));
	if (h[il] > 0.0) {
		double z = h[il];
		double s = f[il];
		double v1 = v[il];
		double t1 = 2*m*E[il]*(1-s)*(v1-1)/(v1+1);
		R[il][0][0] =   s*(1-4*v1-2*m*z) - t1;
		R[il][0][1] = 4*s*(2*v1-1)*mi - 2*s*m*z*z - 2*z*t1;
		R[il][1][0] = 2*s*m;
		R[il][1][1] =   s*(1-4*v1+2*m*z) + t1;
		t1 = 1.0/(t1+s*(3-4*v1));
		for (i = 0; i < 2; i++)
			for (j = 0; j < 2; j++)
				R[il][i][j] *= exp(-2*m*h[il])*t1;
	}
	R[il][2][0] = 1.0; R[il][3][1] = 1.0;
	// Now we work back up, building the 4x2 R matrices...
	for ( ; il > 0; il--) {
		double z = h[il-1];
		double v1 = v[il-1];
		double v2 = v[il];
		double K = ((1+v2)*E[il-1])/((1+v1)*E[il]);
		double t1 = m*z, t2 = 2*t1;
		X[0][0] = X[2][2] = (4*v1-3-K);
		X[0][1] = ((1+8*v1*v2-t2)*(1-K) + (4*t1*v1-6*v2)
			- K*(4*t1*v2-6*v1))*mi;
		X[2][3] = ((1+8*v1*v2+t2)*(K-1) + (4*t1*v1+6*v2)
			- K*(4*t1*v2+6*v1))*mi;
		X[1][0] = X[3][2] = 0.0;
		X[1][1] = X[3][3] = (K*(4*v2-3)-1);
		X[0][2] = ( t2+4*v1-1)*(1-K);
		X[2][0] = (-t2+4*v1-1)*(1-K);
		X[1][3] = (-t2+4*v2-1)*(1-K);
		X[3][1] = ( t2+4*v2-1)*(1-K);
		X[0][3] = ((1+2*(t1+2*v1)*(t1-2*v2))*(1-K) + 2*v2
			- 2*K*v1)*mi;
		X[2][1] = ((1+2*(t1-2*v1)*(t1+2*v2))*(K-1) - 2*v2
			+ 2*K*v1)*mi;
		X[1][2] = 2*m*(K-1);
		X[3][0] = 2*m*(1-K);
		if (f[il-1] != 1.0) {
			double s = E[il-1]*(1-f[il-1])/(1+v1);
			double r1 = (t1+2*v1-1)*s;
			double r2 = (t1-2*v1+1)*s;
			double r3 =  t1+2*v1;
			double r4 =  t1-2*v1;
			F[0][0] = r1;   F[0][1] = r1*r3*mi;
			F[0][2] = r1;   F[0][3] = r1*r4*mi;
			F[1][0] = -s*m; F[1][1] = -s*r3;
			F[1][2] = -s*m; F[1][3] = -s*r4;
			F[2][0] = r2;   F[2][1] = r2*r3*mi;
			F[2][2] = r2;   F[2][3] = r2*r4*mi;
			F[3][0] = -s*m; F[3][1] = -s*r3;
			F[3][2] = -s*m; F[3][3] = -s*r4;
		}
		for (i = 0; i < 2; i++) {
			for (j = 0; j < 2; j++) {
				D[i][j] = D[i+2][j+2] = 1.0;
				D[i][j+2] = exp(-t2);
				D[i+2][j] = (t2<MAX_EXP?exp(t2):DBL_MAX);
			}
		}
		if (f[il-1] == 1.0) {
			double t3 = 1.0/(4*v1-4);
			for (i = 0; i < 4; i++) {
				R[il-1][i][0] = R[il-1][i][1] = 0.0;
				for (j = 0; j < 4; j++) {
					R[il-1][i][0] += D[i][j]*(X[i][j]*R[il][j][0])*t3;
					R[il-1][i][1] += D[i][j]*(X[i][j]*R[il][j][1])*t3;
				}
			}
		} else {
			double t3 = 1.0/(4*v1-4)/f[il-1];
			for (i = 0; i < 4; i++) {
				R[il-1][i][0] = R[il-1][i][1] = 0.0;
				for (j = 0; j < 4; j++) {
					R[il-1][i][0] += D[i][j]*(f[il-1]*X[i][j]*R[il][j][0]
						+ F[i][j]*R[il][j][0])*t3;
					R[il-1][i][1] += D[i][j]*(f[il-1]*X[i][j]*R[il][j][1]
						+ F[i][j]*R[il][j][1])*t3;
				}
			}
		}
	}
	// Fill in B1...
	B1[0][0] = B1[1][0] = B1[0][2] = m, B1[1][2] = -m;
	B1[0][1] = 2*v[0], B1[0][3] = -2*v[0];
	B1[1][1] = B1[1][3] = 2*v[0]-1;
	// Solve for Cinf and Dinf
	// Abuse CDi as a 2x2 matrix and temp storage
	for (i = 0; i < 4; i++) {
		CDi[0] += B1[0][i]*R[0][i][0];
		CDi[1] += B1[0][i]*R[0][i][1];
		CDi[2] += B1[1][i]*R[0][i][0];
		CDi[3] += B1[1][i]*R[0][i][1];
	}
	CDi[3] = CDi[0]*CDi[3]-CDi[1]*CDi[2];
	CDi[2] = -(CDi[1]/(m*m))/CDi[3];
	CDi[3] =  (CDi[0]/(m*m))/CDi[3];
	// With Cinf and Dinf, go back and solve for ABCD.
	for (il = 0; il < nl; il++) {
		ABCD[il][0] = R[il][0][0]*CDi[2]+R[il][0][1]*CDi[3];
		ABCD[il][1] = R[il][1][0]*CDi[2]+R[il][1][1]*CDi[3];
		ABCD[il][2] = R[il][2][0]*CDi[2]+R[il][2][1]*CDi[3];
		ABCD[il][3] = R[il][3][0]*CDi[2]+R[il][3][1]*CDi[3];
		//printf("%g\t%g\t%g\t%0.16g\t%0.16g\t%0.16g\t%0.16g\n",m,h[il],E[il],ABCD[il][0],ABCD[il][1],ABCD[il][2],ABCD[il][3]);
	}
}

/*
 * This function builds the common variables for the various integrations.
 */
static void
buildT(const double m, const double z, const double (&ABCD)[4],
		double (&T)[4])
{
	T[1] = exp(-m*z);
	T[0] = m*(ABCD[2] + ABCD[3]*z)*T[1];
	T[1] *= ABCD[3];
	T[3] = (m*z < MAX_EXP ? exp(m*z) : DBL_MAX);
	T[2] = m*(ABCD[0] + ABCD[1]*z)*T[3];
	T[3] *= ABCD[1];
}

/*
 * This is the 'correct' code, according to exactly how the
 * math should be done.  As a result, it is slow...  It is only
 * really here for checking the optimised code below, but it
 * is exposed so that people can use it if they really want.
 */
bool
LEsystem::calc_accurate()
{
	unsigned ixy, ild, ib, igp, il;
	const LElayer * pl;
	
	initarrays();
	if (!check())
		return false;
	unsigned nl = layers();

	// The integration constants, per layer.
	double (* R)[4][2] = new double[nl][4][2];
	double (* ABCD)[4] = new double[nl][4];
	double * h = new double[nl];
	double * f = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	// We allocate these as big as we ever make them,
	// so we never have to worry about the add()'s failing.
	cset<double> bm0(0, 2*NBZ+2);
	cset<double> bm1(0, 2*NBZ+2);
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
		h[il] = pl->bottom();
		f[il] = MAX(0.0,pl->slip());
		v[il] = pl->poissons();
		E[il] = pl->emod();
	}

	// We loop through all of the evaluation points
	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		il = d.il;
		// Now loop through the list of loads...
		for (ild = 0; ild < load.length(); ild++) {
			double a = load[ild].radius();
			double r = load[ild].distance(d);
			double m0, m1;
			axialdata s;
			memset(&s,0,sizeof(axialdata));
			s.active = true;

			// Now gerenate a list of integration intervals, then sort them.
			bm0.empty(), bm1.empty();
			bm0.add(0.0), bm1.add(0.0);
			stoppingpoints(NBZ,a,r,&m0,&m1);
			// Account for big z's.
			for (ib = 1; ib <= 5; ib++) {
				if (ib*7*a/d.z < j0r[0])
					bm0.add(ib*7*a/d.z);
				if (ib*7*a/d.z < j1r[0])
					bm1.add(ib*7*a/d.z);
			}
			for (ib = 0; r > 0 && ib <= NBZ && j0r[ib]/r < m0; ib++)
				bm0.add(j0r[ib]/r);
			for (ib = 0; r > 0 && ib <= NBZ && j1r[ib]/r < m1; ib++)
				bm1.add(j1r[ib]/r);
			for (ib = 0; ib <= NBZ && j1r[ib]/a < m0; ib++)
				bm0.add(j1r[ib]/a);
			for (ib = 0; ib <= NBZ && j1r[ib]/a < m1; ib++)
				bm1.add(j1r[ib]/a);
			bm0.add(m0), bm1.add(m1);
			bm0.sort(), bm1.sort();

			// We loop through all of our roots and gauss points.
			// In this version we do this in two parts, one for the
			// J1(ma)*J0(mr) integrals and one for the J1(ma)*J1(mr).
			for (ib = bm0.length()-1; ib > 0; ib--) {
				for (igp = 0; igp < NGQP; igp++) {
					// Calculate the gauss point and weight.
					double m = (bm0[ib]+bm0[ib-1])/2
							 + gu[NGQP][igp]*(bm0[ib]-bm0[ib-1])/2;
					double w = gf[NGQP][igp]*(bm0[ib]-bm0[ib-1])/2;
					// First build a new ABCD matrix.
					buildabcd(m,nl,h,v,E,f,R,ABCD);
					double t = m*j1(m*a)*w*j0(m*r), T[4];
					buildT(m,d.z,ABCD[il],T);
					s.vse2 += t*m*((1-2*v[il])*(T[1]+T[3])+(T[0]-T[2]));
					s.rse2 += t*m*((1+2*v[il])*(T[1]+T[3])-(T[0]-T[2]));
					s.tse2 += t*m*(2*v[il])*(T[1]+T[3]);
					s.vdp2 += t*((2-4*v[il])*(T[3]-T[1])-(T[0]+T[2]));
				}
				s.accumulate(accurate,0.0);
			}
			for (ib = 1; ib < bm1.length(); ib++) {
				for (igp = 0; igp < NGQP; igp++) {
					// Calculate the gauss point and weight.
					double m = (bm1[ib]+bm1[ib-1])/2
							 + gu[NGQP][igp]*(bm1[ib]-bm1[ib-1])/2;
					double w = gf[NGQP][igp]*(bm1[ib]-bm1[ib-1])/2;
					// First build a new ABCD matrix.
					buildabcd(m,nl,h,v,E,f,R,ABCD);
					double t = m*j1(m*a)*j1(m*r)*w, T[4];
					buildT(m,d.z,ABCD[il],T);
					s.sse2 += t*m*((2*v[il])*(T[3]-T[1]) + (T[0]+T[2]));
					s.rdp2 += t*((T[1]+T[3])-(T[0]-T[2]));
				}
				s.accumulate(accurate,0.0);
			}
			s.finalize(accurate,r,a,v[il],E[il]);
			// Now add this load's axial data into the point's data
			s.addtodata(accurate,&d,load[ild],r);
		}
		// Calculate the derived results (principal stresses and strains).
		d.principle(v[il],E[il]);
	}
	delete [] R;
	delete [] ABCD;
	delete [] h;
	delete [] f;
	delete [] v;
	delete [] E;
	return true;
}

/*
 * This used to be ELSYM5M, now it's a NxNxN layered elastic code...
 */
bool
LEsystem::calculate(resulttype res, const double * Q)
{
	unsigned ixy, nr, ir, nz, iz, ild, ia, ib, igp, il;
	const LElayer * pl;
	double x1, x2;
	bool interpolate = false;
	
	initarrays();
	if (!check())
		return false;
	unsigned ngqp = NGQP, nbz = NBZ, gl = UINT_MAX, nl = layers();
	if ((res & mask) == dirty) {
		ngqp = MIN(NGQP,8);
		nbz = MIN(NBZ,64);
		interpolate = true;
	} else if ((res & mask) == fast) {
		ngqp = MIN(NGQP,8);
		nbz = MIN(NBZ,64);
		//interpolate = true;
	} else {
		ngqp = MIN(NGQP,12);
		nbz = MIN(NBZ,256);
	}
	const double eps = ((res & mask) == dirty ? 1e-6 : 
			((res & mask) == fast ? 1e-8 : 0.0));

	// The integration constants, per layer.
	double (* R)[4][2] = new double[nl][4][2];
	double (* ABCD)[4] = new double[nl][4];
	double (* iT)[2][4] = 0; // interpolated T's.
	// Local variables, so we don't have to walk the list.
	double * h = new double[nl];
	double * f = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	// Some place to store our data...
	cset<double> a, r, bm;
	cset<zpoint> z;
	fset<double> m0(data.length()), m1(data.length()); 
	fset<axialdata> ax(data.length()); 
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
		h[il] = pl->bottom();
		f[il] = MAX(0.0,pl->slip());
		v[il] = pl->poissons();
		E[il] = pl->emod();
	}

	// Collect and sort the z positions.
	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		// If we're collecting displacement gradient results
		// resize the array as needed then zero it.
		if (res & grad) {
			d.deflgrad.resize(nl);
			memset(&(d.deflgrad[0]),0,nl*sizeof(double));
		}
		if (!z.add(zpoint(d.z,d.il)))
			goto abort;
	}
	z.sort();
	nz = z.length();

	if (interpolate)
		iT = new double[nz][2][4];

	// Gerenate a list of load radii, then sort them and map from loads.
	for (ild = 0; ild < load.length(); ild++) {
		if (!a.add(load[ild].radius()))
			goto abort;
	}
	a.sort();

	// Now loop through the list of load radii, calculating only for
	// the applicable loads... (load[ild].radius() == a[ia])
	for (ia = 0; ia < a.length(); ia++) {
		// Gerenate a list of radii, then sort them.
		r.empty();
		for (ild = 0; ild < load.length(); ild++) {
			if (fabs(load[ild].radius()-a[ia]) > DBL_MIN)
				continue;
			for (ixy = 0; ixy < data.length(); ixy++) {
				if (!r.add(load[ild].distance(data[ixy])))
					goto abort;
			}
		}
		r.sort();
		nr = r.length();
		m0.resize(nr);
		m1.resize(nr);

		// Now gerenate a list of integration intervals, then sort them.
		bm.empty();
		if (!bm.add(0.0))
			goto abort;
		for (ib = 0; ib < (nbz+1); ib++) {
			if (!bm.add(j1r[ib]/a[ia]))
				goto abort;
		}

		// The correct stopping points for each radius.
		for (ir = 0; ir < nr; ir++) {
			stoppingpoints(nbz,a[ia],r[ir],&m0[ir],&m1[ir]);
			if (!bm.add(m0[ir]) || !bm.add(m1[ir]))
				goto abort;
		}
		bm.sort();

		// Account for big r's by adding extra integration intervals...
		x1 = 0.0, x2 = 0.0;
		for (ir = nr; ir > 0 && r[ir-1] > a[ia]*MAX(4,ngqp-6); ir--) {
			for (unsigned i = MAX(4,ngqp-6);
					i <= nbz; i += MAX(4,ngqp-6)) {
				if ((x1 = j1r[i]/r[ir-1]) < x2)
					continue;
				for (ib = 1; ib < bm.length() && bm[ib] < x1; ib++)
					;
				if (ib == bm.length() ||
						MIN(x1-bm[ib-1],bm[ib]-x1)*r[ir-1]
							 < MAX(4,ngqp-6)*M_PI_4)
					continue;
				if (!bm.add(ib,x2 = x1))
					goto abort;
			}
		}
		// Account for big z's.  We drop approximately three orders of
		// magnitude for exp(-7).  Add 5 intervals, so we drop 15 orders
		// of magnitude.
		x1 = 0.0, x2 = 0.0;
		for (iz = nz; iz > 0 && z[iz-1].z > 0.0; iz--) {
			for (unsigned i = 1; i <= 5
					&& i*7*a[ia] < j0r[0]*z[iz-1].z; i++) {
				if ((x1 = i*7*a[ia]/z[iz-1].z) < x2)
					continue;
				for (ib = 1; ib < bm.length() && bm[ib] < x1; ib++)
					;
				if (ib == bm.length() ||
						MIN(x1-bm[ib-1],bm[ib]-x1)*z[iz-1].z < 5*a[ia])
					continue;
				if (!bm.add(ib,x2 = x1))
					goto abort;
			}
		}
		bm.sort();
		
gradloop:
		// And finally, somewhere to stick the radial data...
		ax.resize(nr*nz);
		memset(&ax[0],0,sizeof(axialdata)*nz*nr);
		// Compute the active set.
		for (ixy = 0; ixy < data.length(); ixy++) {
			pavedata & d = data[ixy];
			iz = z.findvalue(zpoint(d.z,d.il));
			for (ild = 0; ild < load.length(); ild++) {
				if (fabs(load[ild].radius()-a[ia]) > DBL_MIN)
					continue;
				ir = r.findvalue(load[ild].distance(d));
				ax[iz*nr+ir].active = true;
			}
		}
		
		// Now that we know the radii, get down to work.
		// We loop through all of our roots and gauss points.
		for (ib = 1; ib < bm.length(); ib++) {
			bool alldone = true;
			bool firstpanel = (bm[ib]*a[ia] <= j1r[1]);
			unsigned agqp = (firstpanel && (res & mask) != dirty ? NGQP : ngqp);
			if (interpolate) {
				if (!firstpanel) {
					for (iz = 0; iz < nz; iz++)
						memcpy(iT[iz][0],iT[iz][1],4*sizeof(double));
				}
				buildabcd(bm[ib],nl,h,v,E,f,R,ABCD);
				for (iz = 0; iz < nz; iz++)
					buildT(bm[ib],z[iz].z,ABCD[z[iz].il],iT[iz][1]);
			}
			for (igp = 0; igp < agqp; igp++) {
				// Calculate the gauss point and weight.
				double dm = (bm[ib]-bm[ib-1])/2;
				double m = (bm[ib]+bm[ib-1])/2 + gu[agqp][igp]*dm;
				double w = gf[agqp][igp]*dm;
				w *= m*j1(m*a[ia]);

				// First build a new ABCD matrix.
				if (!interpolate || firstpanel)
					buildabcd(m,nl,h,v,E,f,R,ABCD);
				// Now calculate the integrals.
				for (iz = 0; iz < nz; iz++) {
					double T[4];
					il = z[iz].il;
					const double tv = 2*v[il];
					if (!interpolate || firstpanel) {
						buildT(m,z[iz].z,ABCD[il],T);
					} else {
						for (unsigned i = 0; i < 4; i++)
							T[i] = iT[iz][1][i] -
									(iT[iz][1][i]-iT[iz][0][i])
										*(bm[ib]-m)/(2*dm);
					}
					for (ir = 0; ir < nr; ir++) {
						axialdata & s = ax[iz*nr+ir];
						if (!s.active)
							continue;
						alldone = false;
						if (m < m0[ir]) {
							double t = w*j0(m*r[ir]);
							s.vdp2 += t*(2*(1-tv)*(T[3]-T[1])-(T[0]+T[2]));
							if (!(res & disp)) {
								double t1 = T[1]+T[3], t2 = T[0]-T[2];
								t *= m;
								s.vse2 += t*((1-tv)*t1+t2);
								s.rse2 += t*((1+tv)*t1-t2);
								s.tse2 += t*tv*t1;
							}
						}
						if (m < m1[ir] && !(res & disp)) {
							double t = w*j1(m*r[ir]);
							s.sse2 += t*m*(tv*(T[3]-T[1])+(T[0]+T[2]));
							s.rdp2 += t*((T[1]+T[3])-(T[0]-T[2]));
						}
					}
				}
			}
			if (alldone)
				break;
			// Don't accumulate for non-full J1(m*a) panels.
			if (fabs(j1(bm[ib]*a[ia])) > j1max)
				continue;
			for (iz = 0; iz < nz; iz++) {
				for (ir = 0; ir < nr; ir++)
					ax[iz*nr+ir].accumulate(res,eps);
			}
		}

		// Finalise the calculations, now that we have done the integration.
		for (iz = 0; iz < nz; iz++) {
			il = z[iz].il;
			for (ir = 0; ir < nr; ir++)
				ax[iz*nr+ir].finalize(res,r[ir],a[ia],v[il],E[il]);
		}

		// After doing everything in radial coords, translate to cartesian.
		for (ixy = 0; ixy < data.length(); ixy++) {
			pavedata & d = data[ixy];
			iz = z.findvalue(zpoint(d.z,d.il));
			for (ild = 0; ild < load.length(); ild++) {
				if (fabs(load[ild].radius()-a[ia]) > DBL_MIN)
					continue;
				ir = r.findvalue(load[ild].distance(d));
				ax[iz*nr+ir].addtodata(res,&d,load[ild],r[ir],gl);
			}
		}
		
		// Take care of the deflection gradient calculation.
		if (res & grad) {
			res = LEsystem::resulttype(res | disp);
			gl = (gl == UINT_MAX ? 0 : gl+1);
			if (gl != nl) {
				// Deflection gradients are calculated in a log(E) space.
				for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
					E[il] = pl->emod();
					if (Q != 0)
						E[il] = pow(10,log10(E[il])+Q[il*nl+gl]*GRADSTEP);
					else if (il == gl)
						E[il] = pow(10,log10(E[il])+GRADSTEP);
				}
				goto gradloop;
			}
			gl = UINT_MAX;
		}
	}
	// After everything, loop through the answers, and calculate the
	// derived results (principal stresses and strains).
	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		if (!(res & disp))
			d.principle(v[d.il],E[d.il]);
		if (res & grad) {
			if (Q == 0) {
				for (il = 0; il < nl; il++)
					d.deflgrad[il] = (d.deflgrad[il]-d.data[4][2])/GRADSTEP;
			} else {
				for (il = 0; il < nl; il++)
					h[il] = (d.deflgrad[il]-d.data[4][2])/GRADSTEP;
				for (gl = 0; gl < nl; gl++) {
					d.deflgrad[gl] = 0.0;
					for (il = 0; il < nl; il++)
						d.deflgrad[gl] += h[il]*Q[gl*nl+il];
				}
			}
		}
	}
abort:
	delete [] R;
	delete [] ABCD;
	delete [] iT;
	delete [] h;
	delete [] f;
	delete [] v;
	delete [] E;
	return true;
}

/*
 * Boussinesq's equations for a point and circular load.
 */
static inline double
boussinesq_vse(double z, double r, double a, double v, double E,
		double R, double A)
{
	if (r > 0.0)
		return -3*a*a*pow(z,3)/(2*pow(R,5));
	else
		return pow(z,3)/pow(A,3)-1.0;
}

static inline double
boussinesq_rse(double z, double r, double a, double v, double E,
		double R, double A)
{
	if (r > 0.0)
		return -a*a*(3*z*r*r/pow(R,4)-(1-2*v)/(R+z))/(2*R);
	else
		// Note: tse is zero, so we'll divide this by 2 later.
		return 2*(1+v)*z/A-pow(z,3)/pow(A,3)-1-2*v;
}

static inline double
boussinesq_tse(double z, double r, double a, double v, double E,
		double R)
{
	if (r > 0.0)
		return -a*a*(1-2*v)*(-z/(R*R) + 1.0/(R+z))/(2*R);
	else
		return 0.0;
}

static inline double
boussinesq_sse(double z, double r, double a, double v, double E,
		double R)
{
	if (r > 0.0)
		return -3*a*a*z*z*r/(2*pow(R,5));
	else
		return 0.0;
}

static inline double
boussinesq_rdp(double z, double r, double a, double v, double E,
		double R)
{
	if (r > 0.0)
		return a*a*((1+v)/(2*R*E))*(z*r/(R*R)-(1-2*v)*r/(R+z));
	else
		return 0.0;
}

static inline double
boussinesq_vdp(double z, double r, double a, double v, double E,
		double R, double A)
{
	if (r > 0.0)
		return a*a*((1+v)/(2*R*E))*(2*(1-v) + z*z/(R*R));
	else
		return ((1+v)/E)*(a*a/A+(1-2*v)*(A-z));
}

/*
 * This is Odemark's method of equivalent thicknesses, with Bosussinesq's
 * solution for a point load on a linear elastic half-space.  This can be
 * used for very fast approximations.
 */
bool
LEsystem::calc_odemark()
{
	unsigned ixy, ild, il;
	LElayer * pl;

	if (!check())
		return false;
	unsigned nl = layers();
	double * h = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++)
		h[il] = pl->bottom(), v[il] = pl->poissons(), E[il] = pl->emod();

	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		for (ild = 0; ild < load.length(); ild++) {
			double a = load[ild].radius();
			double r = load[ild].distance(d);
			double z = 0.0, he = 0.0, hc = 1.0, R, A;
			axialdata s;
			memset(&s,0,sizeof(axialdata));
			for (il = 0; il < nl; il++) {
				if (il > 0) {
					z = h[il-1]-(il > 1 ? h[il-2] : 0.0);
					he = (he+z)*pow(E[il-1]/E[il]*
						(1+v[il]*v[il])/(1+v[il-1]*v[il-1]),1.0/3.0);
					if (nl > 2 && il > 1) // Per's fudge factors...
						hc = 0.8;
					else if (nl > 2)
						hc = 1.0;
					else
						hc = 0.9;
				}
				if (il < d.il)
					continue;
				if (il > d.il) {
					// We're below, but we still need the vertical
					// deflection contributions from this layer...
					z = hc*he;
					R = hypot(r,z); A = (r > 0.0 ? 0.0 : hypot(z,a));
					s.vdp += boussinesq_vdp(z,r,a,v[il],E[il],R,A);
					if (h[il] == 0.0)
						continue;
					z = hc*he+h[il]-h[il-1];
					R = hypot(r,z); A = (r > 0.0 ? 0.0 : hypot(z,a));
					s.vdp -= boussinesq_vdp(z,r,a,v[il],E[il],R,A);
					continue;
				}
				z = hc*he+d.z-(il > 0 ? h[il-1] : 0.0);
				R = hypot(r,z); A = (r > 0.0 ? 0.0 : hypot(z,a));
				double tv = v[d.il];
				double tE = E[d.il];
				s.vse = boussinesq_vse(z,r,a,tv,tE,R,A);
				s.rse = boussinesq_rse(z,r,a,tv,tE,R,A);
				s.tse = boussinesq_tse(z,r,a,tv,tE,R);
				s.sse = boussinesq_sse(z,r,a,tv,tE,R);
				s.rdp = boussinesq_rdp(z,r,a,tv,tE,R);
				s.vdp = boussinesq_vdp(z,r,a,tv,tE,R,A);
				if (h[il] == 0.0)
					continue;
				z = hc*he+h[il]-(il > 0 ? h[il-1] : 0.0);
				R = hypot(r,z); A = (r > 0.0 ? 0.0 : hypot(z,a));
				s.vdp -= boussinesq_vdp(z,r,a,tv,tE,R,A);
			}
			// After doing everything in radial coords, translate
			// to cartesian.
			s.addtodata(odemark,&d,load[ild],r);
		}
		d.principle(v[d.il],E[d.il]);
	}
	delete [] h;
	delete [] v;
	delete [] E;
	return true;
}

/*
 * The following equations are derived from the Boussinesq solution,
 * integrated and simplified as much as possible.  However, it cannot
 * completely integrated because there is a problem at points where
 * the load is directly over the evaluation point, so...
 *
 * An implementation of adaptive, recursive Newton-Cotes integration.
 * Based on the Matlab implementation, but covered in a lot of books...
 */
#define LEVMAX	12

static double
quad8_vdp(double r, double z, double a, double v, double A = 0.0,
          double B = M_PI, double Q = 10.0)
{
	// The magic Newton-Cotes weights
	const double w[9] = {3956, 23552, -3712, 41984, -18160, 41984,
						-3712, 23552, 3956};
	const double dw = 14175;
	static unsigned level = 0;
	static double tol = 1e-6;
	register double h, t, Q1 = 0.0, Q2 = 0.0;
	register unsigned i;

	if (a == 0.0)
		return 0.0;
	if (r == 0.0)
		return (2*(v*v-1)*(z-hypot(z,a))
						+(v+1)*(z-z*z/hypot(z,a)));
	level++;
	h = (B-A)/16.0;
	for (i = 0; i < 9; i++) {
		t = cos(A+i*h)*r;
		if (z == 0.0)
			Q1 += h*w[i]/dw*(t == r ? 2*(v*v-1)*(r-sqrt(a*a-2*a*t+r*r)) :
				2*(v*v-1)*(t*(log(r-t)-log(a-t+sqrt(a*a-2*a*t+r*r)))
					-sqrt(a*a-2*a*t+r*r)+r)
			);
		else
			Q1 += h*w[i]/dw*(t*t == r*r+z*z ? 2*(v*v-1)*a :
				2*(v*v-1)*
					(t*(log(hypot(z,r)-t)
					   -log(a-t+sqrt(z*z-2*a*t+a*a+r*r)))
					 -sqrt(z*z-2*a*t+a*a+r*r)+hypot(z,r))
				 +(v+1)*z*z*((z*z+r*r-t*a)/sqrt(z*z-2*a*t+a*a+r*r)
					 -hypot(z,r))/(t*t-z*z-r*r)
			);
		t = cos(A+(i+8)*h)*r;
		if (z == 0.0)
			Q2 += h*w[i]/dw*(t == r ? 2*(v*v-1)*(r-sqrt(a*a-2*a*t+r*r)) :
				2*(v*v-1)*(t*(log(r-t)-log(a-t+sqrt(a*a-2*a*t+r*r)))
					-sqrt(a*a-2*a*t+r*r)+r)
			);
		else
			Q2 += h*w[i]/dw*(t*t == r*r+z*z ? 2*(v*v-1)*a :
				2*(v*v-1)*
					(t*(log(hypot(z,r)-t)
					   -log(a-t+sqrt(z*z-2*a*t+a*a+r*r)))
					 -sqrt(z*z-2*a*t+a*a+r*r)+hypot(z,r))
				 +(v+1)*z*z*((z*z+r*r-t*a)/sqrt(z*z-2*a*t+a*a+r*r)
					 -hypot(z,r))/(t*t-z*z-r*r)
			);
	}
	// This is the adaptive recursive bit.  We only recurse if we
	// can improve...
	if (fabs(Q1+Q2-Q) > tol*fabs(Q1+Q2) && level <= LEVMAX) {
		tol = tol/2;
		Q1 = quad8_vdp(r,z,a,v,A,(A+B)/2,Q1);
		Q2 = quad8_vdp(r,z,a,v,(A+B)/2,B,Q2);
		tol = tol*2;
	}
	level--;
	if (level == 0)
		return (Q1 + Q2)/M_PI;
	else
		return (Q1 + Q2);
}

static double
quad8_vse(double r, double z, double a, double A = 0.0,
          double B = M_PI, double Q = 10.0)
{
	// The magic Newton-Cotes weights
	const double w[9] = {3956, 23552, -3712, 41984, -18160, 41984,
						-3712, 23552, 3956};
	const double dw = 14175;
	static unsigned level = 0;
	static double tol = 1e-6;
	register double h, t, t1, Q1 = 0.0, Q2 = 0.0;
	register unsigned i;

	if (a == 0.0)
		return 0.0;
	if (z == 0.0)
		return (r == a ? -0.5 : r < a ? -1.0 : 0.0);
	if (r == 0.0)
		return (pow(z/hypot(z,a),3)-1);

	level++;
	h = (B-A)/16.0;
	for (i = 0; i < 9; i++) {
		t = cos(A+i*h);
		t1 = t*t*r*r + r*r + z*z;
		Q1 += h*w[i]/dw*(
			z*z*z*(((r*r+z*z-3*r*a*t)*t1
					- 2*t*r*a*(a*a-3*a*t*r)) / pow(sqrt(a*a-2*r*a*t+r*r+z*z),3)
				-t1/hypot(r,z)) / pow((1-t*t)*r*r+z*z,2));
		t = cos(A+(i+8)*h);
		t1 = t*t*r*r + r*r + z*z;
		Q2 += h*w[i]/dw*(
			z*z*z*(((r*r+z*z-3*r*a*t)*t1
					- 2*t*r*a*(a*a-3*a*t*r)) / pow(sqrt(a*a-2*r*a*t+r*r+z*z),3)
				-t1/hypot(r,z)) / pow((1-t*t)*r*r+z*z,2));
	};
	// This is the adaptive recursive bit.  We only recurse if we
	// can improve...
	if (fabs(Q1+Q2-Q) > tol*fabs(Q1+Q2) && level <= LEVMAX) {
		tol = tol/2;
		Q1 = quad8_vse(r,z,a,A,(A+B)/2,Q1);
		Q2 = quad8_vse(r,z,a,(A+B)/2,B,Q2);
		tol = tol*2;
	}
	level--;
	if (level == 0)
		return (Q1 + Q2)/M_PI;
	else
		return (Q1 + Q2);
}

#undef LEVMAX

/*
 * This is Odemark's method of equivalent thicknesses, with Boussinesq's
 * solution for a circular load on a linear elastic half-space.  This can
 * be used for very fast approximations of vertical deflections under a
 * circular load.
 *
 * XXX: This could be expanded to all of the outputs, if some work was
 * done to derive more integration functions above.  But these are long...
 */
bool
LEsystem::calc_fastnum()
{
	unsigned ixy, ild, il;
	LElayer * pl;

	if (!check())
		return false;
	unsigned nl = layers();
	double * h = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++)
		h[il] = pl->bottom(), v[il] = pl->poissons(), E[il] = pl->emod();

	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		for (ild = 0; ild < load.length(); ild++) {
			double a = load[ild].radius();
			double z, r = load[ild].distance(d);
			double he = 0.0, hc = 1.0;
			double vdp = 0.0, vse = 0.0;
			for (il = 0; il < nl; il++) {
				if (il > 0) {
					z = h[il-1]-(il > 1 ? h[il-2] : 0.0);
					he = (he+z)*pow(E[il-1]/E[il]*
						(1+v[il]*v[il])/(1+v[il-1]*v[il-1]),1.0/3.0);
					if (nl > 2 && il > 1) // Per's fudge factors...
						hc = 0.8;
					else if (nl > 2)
						hc = 1.0;
					else
						hc = 0.9;
				}
				if (h[il] != 0.0 && d.z >= h[il])
					continue;
				if (il > 0 && h[il-1] > d.z) {
					// We're below, but we still need the vertical
					// deflection contributions from this layer...
					z = hc*he;
					vdp += quad8_vdp(r,z,a,v[il])/E[il];
					if (h[il] != 0.0) {
						z = hc*he+h[il]-h[il-1];
						vdp -= quad8_vdp(r,z,a,v[il])/E[il];
					}
					continue;
				}
				z = hc*he+d.z-(il > 0 ? h[il-1] : 0.0);
				vse = quad8_vse(r,z,a);
				vdp += quad8_vdp(r,z,a,v[il])/E[il];
				if (h[il] != 0.0) {
					z = hc*he+h[il]-(il > 0 ? h[il-1] : 0.0);
					vdp -= quad8_vdp(r,z,a,v[il])/E[il];
				}
			}
			d.data[0][2] += load[ild].pressure()*vse;
			d.data[4][2] += load[ild].pressure()*vdp;
		}
	}
	delete [] h;
	delete [] v;
	delete [] E;
	return true;
}

/*
 * The overall backcalculation routine.
 */
bool
LEbackcalc::backcalc()
{
	unsigned i, j, nl = layers();
	unsigned steps = 0;
	//double derr = DBL_MAX, oerr = 0.0;
	double astep, tstep = 0.0, step = 1.0;
	bool seeded = true;
	//bool badstep = false;
	calctype speed = (precision >= 1e-4 ? fast : slow);
	ksset<pavepoint,pavedata> orig(data);

	if (nl == 0) {
		event_msg(EVENT_ERROR,"LEbackcalc::backcalc() called without layers!");
		return false;
	}
	if (defl.length() == 0) {
		event_msg(EVENT_ERROR,"LEbackcalc::backcalc() called without deflections!");
		return false;
	}
	double * P = new double[nl];
	double * T = new double[nl];
	double * O = new double[nl];
	double * Q = new double[nl*nl];

	// Remove all existing evaluation points, after saving...
	removepoints();
	// Add all of the measure points as evaluation points
	for (i = 0; i < defl.length(); i++) {
		defldata & d = defl[i];
		addpoint(d);
	}
	// We work in a log(E) space to avoid having to constrain the
	// optimisation.  Throughout P[] is the current point.
	for (i = 0; i < nl; i++) {
		if (layer(i).emod() <= 0.0)
			seeded = false;
		else
			O[i] = T[i] = P[i] = log10(layer(i).emod());
		for (j = 0; j < nl; j++)
			Q[i*nl+j] = (i == j ? 1.0 : 0.0);
	}
	// If we have starting values, try to work out if they're good.
	// If not, then we resort to reseeding the problem. We do this by
	// using the kalman method, since this is the best test of if we
	// can get the error to improve.
	if (seeded) {
		calculate(LEsystem::dispgrad);
		step = kalman(nl,P);
		if (step < tolerance)
			return true;
	}
	if (!seeded || step > 0.5) {
		//badstep =
		seed(nl,P);
		memcpy(T,P,sizeof(double)*nl);
		memcpy(O,P,sizeof(double)*nl);
	}
	//swarm(nl,P);
	// Now settle into our true minimisation algorithm.
	// First we start by using the deflection gradient in
	// emod space, since it converges faster...
	while (true) {
		if (++steps > maxsteps*nl)
			break;
		memcpy(T,P,sizeof(double)*nl);
		astep = step, tstep += step = deflgrad(nl,P,Q,speed);
		//printf("E%s ",speed == fast ? "*" : speed == slow ? "!" : "@");
		//printf("step = %10.8f ",step);
		if (astep < 0.1 && step > astep) {
			//printf(" Bailing...\n");
			memcpy(P,T,sizeof(double)*nl);
			break;
		}
		//for (i = 0, astep = 0.0; i < nl; i++)
		//	astep += (O[i]-P[i])*(O[i]-P[i]);
		//astep = sqrt(astep);
		//printf("tstep = %10.8f (%10.8f) (%6.4f) %d\n",tstep,astep,tstep/astep,ns);
		//for (i = 0; i < nl; i++)
		//	printf(" %8.5g",P[i]);
		//printf("\n");
		if (step < tolerance)
			break;
		//if (tstep/astep > 3.0 && P[0] > P[1] && !badstep) {
		//	// We're wandering around...  This normally happens if the
		//	// surface stiffness is low.
		//	// XXX: We should grid search this...
		//	P[0] = 1.0, badstep = true;
		//	//printf("Bad step correction!\n");
		//}
		//if (ns >= 2 && !badstep) {
		//	// We're not converging...  This nornally happens for stiff
		//	// surfaces on weak subgrade.
		//	astep = bowlerror(nl,P);
		//	//printf("%f %f %f\n",P[0],P[1],astep);
		//	for (i = 0; i < nl-2; i++) {
		//		for (j = 0; j < 10; j++) {
		//			P[i] += 0.1;
		//			step = bowlerror(nl,P);
		//			//printf("%f %f %f (%f)\n",P[0],P[1],step,step/astep);
		//			if (step/astep < 1.0)
		//				break;
		//			if (step/astep > 100.0) {
		//				P[i] -= 0.1;
		//				break;
		//			}
		//		}
		//		if (step/astep < 1.0)
		//			break;
		//	}
		//	//for (i = 0; i < nl; i++)
		//	//	printf(" %8.5g",P[i]);
		//	//printf("\n");
		//	badstep = true, step = 1.0;
		//	//printf("Slow step correction!\n");
		//}
		//if (false && ns >= 6)
		//	break;
	}
	// Now change to the Kalman filter to converge in deflection space
	while (true) {
		if (++steps > maxsteps*nl)
			break;
		for (i = 0; i < nl; i++)
			layer(i).emod(pow(10,P[i]));
		calculate((speed == fast ? LEsystem::fastgrad : LEsystem::dispgrad));
		//oerr = derr;
		//derr = (precision > 0.0
		//	? ROUND(bowlerror()/precision)*precision : bowlerror());
		//printf("K%s ",speed == fast ? "*" : speed == slow ? "!" : "@");
		//printf("err = %g ",derr);
		//if (oerr + FLT_EPSILON < derr && ++ns > 2) { // Ooops, our error increased...
		//	printf(" Bailing...\n");
		//	memcpy(P,T,sizeof(double)*nl);
		//	break;
		//}
		memcpy(T,P,sizeof(double)*nl);
		step = kalman(nl,P);
		//printf("step = %g\n",step);
		if (step < tolerance) // We're not moving any more...
			break;
	}
	// If we still haven't converged, resort to brute force...
	//if (step > tolerance) {
	//	step = conjgrad(nl,P);
	//	printf("G! ");
	//	printf("step = %g\n",step);
	//}
	// Now put the modulli back and calculate the final deflections.
	for (i = 0; i < nl; i++)
		layer(i).emod(pow(10,P[i]));
	calculate(LEsystem::disp);
	for (i = 0; i < defl.length(); i++) {
		defldata & d = defl[i];
		d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
	}
	// Restore original evaluation points...
	data = orig;
	delete [] Q;
	delete [] O;
	delete [] T;
	delete [] P;
	return (step <= tolerance);
}

/*
 * E mod seeding algorithm.
 */
bool
LEbackcalc::seed(unsigned nl, double * P)
{
	double E = 0.0, v = layer(nl-1).poissons();
	unsigned i, j;
	bool negdefl = false;

	// Calculate the average surface modulus based on each deflection.
	for (i = 0; i < defl.length(); i++) {
		defldata & d = defl[i];
		if (d.measured <= 0.0) {
			negdefl = true;
		} else {
			for (j = 0; j < load.length(); j++) {
				E += load[j].pressure()*quad8_vdp(load[j].distance(d),
					d.z,load[j].radius(),v)/MAX(d.measured,1e-4);
			}
		}
	}
	E = log10(MIN(MAX(1e3,E/i),1e6));
	for (i = 0; i < nl; i++)
		P[i] = (negdefl ? 2.0 : E);
	return negdefl;
}

/*
 * Deflection gradient approach.
 *
 * This works by computing the gradient in deflection at each
 * deflection point, with respect to a change in each layer modulus.
 * It then computes the best step to make the calculated delfection
 * equal to the measured deflection.
 *
 * The does not minimise the error, but rather finds the modulus values
 * which minimises the distance between the 'zero' moduli for each
 * deflection point.  The minimisation is performed using a conjugate
 * gradient algorithm, which can solve the quadratic problem exactly.
 *
 * Only one step is performed, and the algorithm should be called
 * repeatedly to improve the approximation.
 *
 * See the ICAP 2006 paper for details.
 */
double
LEbackcalc::deflgrad(unsigned nl, double * P, double * Q,
                     calctype cl)
{
	double step = 0.0, dgg = 0.0, gg = 0.0, dd = 0.0;
	unsigned i, j, k, dl = defl.length();

	// Initial setup.
	double * PG = new double[dl*nl];
	double * MD = new double[dl];
	double * CD = new double[dl];
	double * DG = new double[dl];
	double * D = new double[nl];
	double * G = new double[nl];
	double * W = new double[nl];
	memset(W,0,sizeof(double)*nl);
	for (i = 0; i < nl; i++)
		layer(i).emod(pow(10,P[i]));
	switch (cl) {
	case slow:  calculate(LEsystem::dispgrad,Q); break;
	case fast:  calculate(LEsystem::fastgrad,Q); break;
	case reuse: break;
	}
	for (j = 0; j < dl; j++) {
		defldata & d = defl[j];
		d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
		CD[j] = d.calculated;
		MD[j] = d.measured;
		for (i = 0, DG[j] = 0.0; i < nl; i++) {
			PG[j*nl+i] = result(d).deflgrad[i];
			DG[j] += PG[j*nl+i]*PG[j*nl+i];
		}
		gg += 1.0/sqrt(DG[j]);
	}
	for (j = 0; j < dl; j++)
		DG[j] *= (dl*dl)/(gg*gg);
	for (k = 0; k < (nl+1); k++) {
		if (dgg != 0.0) {
			// Calculate our step, which we know exactly, because
			// we're quadratic.
			double q = 0.0, p = 0.0, r;
			for (i = 0, r = 0.0; i < nl; i++)
				G[i] = D[i], r += D[i]*D[i];
			for (i = 0, r = sqrt(r); i < nl; i++)
				D[i] /= r;
			for (j = 0; j < dl; j++) {
				for (i = 0, r = 0.0; i < nl; i++)
					r += PG[j*nl+i]*D[i];
				q += (CD[j]-MD[j])*r/DG[j], p += r*r/DG[j];
			}
			// Apply the step and update all of the variables.
			for (i = 0, r = 0.0; i < nl; i++) {
				if (Q != 0)
					Q[i*nl+(k-1)] = D[i];
				D[i] *= -q/p;
				r += pow(W[i]+D[i],2);
				for (j = 0; j < dl; j++)
					CD[j] += PG[j*nl+i]*D[i];
			}
			if (r > 1.0) {
				for (i = 0, r = sqrt((1.0-dd)/r); i < nl; i++)
					W[i] += r*D[i];
				break;
			}
			for (i = 0, dd = r; i < nl; i++)
				W[i] += D[i];
		}
		// Find our gradient vector.
		for (i = 0, gg = dgg, dgg = 0.0; i < nl; i++) {
			for (j = 0, D[i] = 0.0; j < dl; j++)
				D[i] -= 2*(CD[j]-MD[j])*PG[j*nl+i]/DG[j];
			dgg += D[i]*D[i];
		}
		// We're pushing the limits of numerical stability
		if (dgg <= precision)
			break;
		for (i = 0; gg != 0.0 && i < nl; i++)
			D[i] += dgg*G[i]/gg;
	}
	for (i = 0, step = 0.0; i < nl; i++)
		P[i] += W[i], step += W[i]*W[i], P[i] = MIN(MAX(1,P[i]),9);
	step = sqrt(step);
	// Orthonormalize Q.
	if (Q != 0)
		orth_gs(nl,Q);
	delete [] W;
	delete [] G;
	delete [] D;
	delete [] DG;
	delete [] CD;
	delete [] MD;
	delete [] PG;
	return step;
}

/*
 * Gauss-Newton approach.
 */
double
LEbackcalc::gaussnewton(unsigned nl, double * P, calctype cl)
{
	unsigned i, j, k, dl = defl.length();
	double step = 0.0;

	// Initial setup.
	double * H = new double[dl*nl];
	double * S = new double[nl*nl];
	double * W = new double[nl];
	double * Y = new double[nl];
	memset(Y,0,sizeof(double)*nl);
	for (i = 0; i < nl; i++)
		layer(i).emod(pow(10,P[i]));
	switch (cl) {
	case slow:  calculate(LEsystem::dispgrad); break;
	case fast:  calculate(LEsystem::fastgrad); break;
	case reuse: break;
	}
	for (j = 0; j < dl; j++) {
		defldata & d = defl[j];
		d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
		double e = d.measured-d.calculated;
		for (i = 0; i < nl; i++) {
			H[j*nl+i] = result(d).deflgrad[i];
			Y[i] += e*H[j*nl+i];
		}
	}
	for (i = 0; i < nl; i++) {
		for (j = 0; j < nl; j++) {
			for (k = 0, S[i*nl+j] = 0.0; k < dl; k++)
				S[i*nl+j] += H[k*nl+i]*H[k*nl+j];
		}
	}
	equ_eig(nl,S,Y,W); // try
	for (i = 0, step = 0.0; i < nl; i++)
		step += W[i]*W[i], P[i] += W[i];
	step = sqrt(step);
	delete [] Y;
	delete [] W;
	delete [] S;
	delete [] H;
	return step;
}

/*
 * This is an alternative algorithm based on the Extended Kalman
 * Filter, which is used in signal processing.  Original idea by
 * Rongzong Wu <wurong@berkeley.edu>.  Filter implementation based
 * on Wikipedia article.
 */
double
LEbackcalc::kalman(unsigned nl, double * P)
{
	unsigned i, j, k, dl = defl.length();
	double step = 0.0;

	// Initial setup.
	double * H = new double[dl*nl];
	double * S = new double[dl*dl];
	double * K = new double[nl*dl];
	double * Y = new double[dl];
	double * W = new double[nl];
	for (j = 0; j < dl; j++) {
		defldata & d = defl[j];
		d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
		Y[j] = d.measured-d.calculated;
		for (i = 0; i < nl; i++)
			H[j*nl+i] = result(d).deflgrad[i];
	}
	for (i = 0; i < dl; i++) {
		for (j = 0; j < dl; j++) {
			for (k = 0, S[i*dl+j] = 0.0; k < nl; k++)
				S[i*dl+j] += H[i*nl+k]*H[j*nl+k];
		}
	}
	for (i = 0; i < dl; i++)
		S[i*dl+i] += noise*noise + (precision*precision)/3.0;
	inv_eig(dl,S); // try
	for (i = 0; i < nl; i++) {
		for (j = 0; j < dl; j++) {
			for (k = 0, K[i*dl+j] = 0.0; k < dl; k++)
				K[i*dl+j] += H[k*nl+i]*S[k*dl+j];
		}
	}
	for (i = 0, step = 0.0; i < nl; i++) {
		for (j = 0, W[i] = 0.0; j < dl; j++)
			W[i] += K[i*dl+j]*Y[j];
		step += W[i]*W[i];
	}
	step = sqrt(step);
	//if (step > 1.0) {
	//	step = step*brent(nl,P,W);
	//} else {
		for (i = 0; i < nl; i++)
			P[i] = MIN(MAX(1,P[i]+W[i]),9);
	//}
	delete [] W;
	delete [] Y;
	delete [] K;
	delete [] S;
	delete [] H;
	return step;
}

/*
 * Bowl error function.  This works out mean deviation of the difference
 * between the caclulated and measured deflections.  Call with P=0 to
 * avoid it rerunning the calcs.  Call with s and D to get it to do
 * prediction at distance s along the vector D.
 */
double
LEbackcalc::bowlerror(unsigned nl, const double * P, const double s,
                      const double * D)
{
	unsigned i;
	double err;
	
	for (i = 0; P != 0 && i < nl; i++)
		layer(i).emod(pow(10,P[i]+(D != 0 ? s*D[i] : 0.0)));
	if (P != 0)
		calculate(LEsystem::disp);
	for (i = 0, err = 0.0; i < defl.length(); i++) {
		defldata & d = defl[i];
		if (P != 0)
			d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
		err += pow(d.measured-d.calculated,2);
	}
	return sqrt(err/defl.length());
}

/*
 * Brent's method for 1D function minimisation, with initial bracketing.
 * This is based around the code from "Numerical Recipies in C",
 * although the algorithms are in much more general use.  The function
 * takes a direction vector D, which defines the search direction.
 * We work in the log space of elastic modulus, so that we don't have to
 * constrain the minimisation.
 *
 * We start with one value, and construct an interval around it.  Then
 * we use that interval to determine if we have bracketed the minimum.
 * Once we have bracketed the minimum we close in on the final minimum.
 */
double
LEbackcalc::brent(unsigned nl, double * P, double * D)
{
	const double GC = (3.0-sqrt(5.0))/2.0;
	const double GR = 1/(1-GC);
	const double TOL = 1e-6;

	double A, B, C, X, Y, Z, a, b, c, x, y, z;
	double r, q;
	unsigned i;
	
	b = 0.0, B = bowlerror(nl,P,b,D);
	// To start with we only have one modulus.  We need three,
	// so dream up two more...  We think that D is a descent
	// direction, so take a step in the opposite direction to get a
	// left bracket.
	a = -0.05, A = bowlerror(nl,P,a,D);
	if (B > A) {
		// D is not a descent direction!  Go the other way...
		for (i = 0; i < nl; i++)
			D[i] = -D[i];
		c = -a, C = A, a = -c/GR, A = bowlerror(nl,P,a,D);
	} else {
		c = -a*GR, C = bowlerror(nl,P,c,D);
	}
	// Now search until we bracket a minimum...
	// We limit the step size to one order of magnitude.
	while (A >= B && B >= C ) { // && c < 3.0
		r = (b-a)*(B-C), q = (b-c)*(B-A);
		x = (r != q ? b-((b-a)*r-(b-c)*q)/(2*(q-r)) : 1.0);
		// x is the turning point of a parabola through A, B, C.
		if (x > b && x < c) {
			X = bowlerror(nl,P,x,D);
			if (X < C) {
				// A > B > C && B > X < C so we've turned between
				// b and c.  Swap c & x, then let them be shifted so we
				// terminate normally.
				swap(x,c), swap(C,X);
			} else if (X > B) {
				// A > B > C && A > B < X so we've turned between
				// A and X.  Replace c with x and skip the final shift.
				c = x, C = X;
				continue;
			}
		} else if (x > c) { // && x < 3.0
			X = bowlerror(nl,P,x,D);
			if (X < C) {
				// The function keeps going down, so take another step.
				b = c, c = x, x = c+(c-b)*GR;
				B = C, C = X, X = bowlerror(nl,P,x,D);
			}
		} else {
			// The parabola didn't work out, so take a golden section step.
			x = c+(c-b)*GR, X = bowlerror(nl,P,x,D);
		}
		a = b, b = c, c = x, A = B, B = C, C = X;
	}
	// Now that we know we have a minimum, hunt it down...
	// Since we've already evaluated the function at the end points,
	// use this to initialise Y and Z.
	if (C < A)
		z = a, y = c, Z = A, Y = C;
	else
		z = c, y = a, Z = C, Y = A;
	// Continue tightening the bracket until we're as close as possible...
	while ((c-a) > TOL*(1.0+b)) {
		// Make a parabolic approximatation of the minimum.
		r = (b-z)*(B-Y), q = (b-y)*(B-Z);
		x = (r != q ? b-((b-z)*r-(b-y)*q)/(2*(r-q)) : 0.0);
		// if we failed or are outside the centre of the two intervals
		// then just take a conservative golden section step.
		if (x == 0.0 || x < (a+b)/2 || x > (b+c)/2)
			x = b+GC*(c-b < b-a ? a-b : c-b);
		X = bowlerror(nl,P,x,D);
		if (X <= B) {
			// We found a new best point.  Update everything.
			(x < b ? c : a) = b;
			z = y, y = b, b = x, Z = Y, Y = B, B = X;
		} else {
			// Our new point wasn't better, but update what we can.
			(x < b ? a : c) = x;
			if (X <= Y || x == y)
				z = y, y = x, Z = Y, Y = X;
			else if (X <= Z || x == z || y == z)
				z = x, Z = x;
		}
	}
	for (i = 0; i < nl; i++)
		P[i] = P[i]+b*D[i];
	return b;
}

/*
 * This is a straight conjugate gradient minimization used to converge on
 * a final best solution.  This is only used if the other methods fail
 * completely.  The routine is entirely self contained, since it works
 * better if it can keep track of its last direction.  It is slow.
 */
double
LEbackcalc::conjgrad(unsigned nl, double * P)
{
	double step = 1.0, dgg = 0.0, gg = 0.0;
	double derr = DBL_MAX, oerr = 0.0;
	unsigned i, j, k, dl = defl.length();

	// Initial setup.
	double * D = new double[nl];
	double * G = new double[nl];
	double * W = new double[nl];
	memcpy(W,P,sizeof(double)*nl);
	memset(D,0,sizeof(double)*nl);
	for (k = 0; k <= maxsteps*nl; k++) {
		if (dgg != 0.0)
			step = dgg*brent(nl,W,D);
		for (i = 0; i < nl; i++)
			layer(i).emod(pow(10,W[i]));
		calculate(LEsystem::dispgrad);
		gg = dgg, dgg = 0.0;
		// Find our gradient vector.
		for (i = 0; i < nl; i++) {
			G[i] = D[i], D[i] = 0.0;
			for (j = 0; j < dl; j++) {
				defldata & d = defl[j];
				d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
				D[i] -= 2*(d.calculated-d.measured)*result(d).deflgrad[i];
			}
			dgg += D[i]*D[i];
		}
		oerr = derr, derr = bowlerror();
		if (dgg == 0.0 || log10(oerr)-log10(derr) < 0.001 || step < tolerance) 
			break;
		for (i = 0; gg != 0.0 && i < nl; i++)
			D[i] += dgg*G[i]/gg;
	}
	for (i = 0, gg = 0.0; i < nl; i++)
		gg += (P[i]-W[i])*(P[i]-W[i]), P[i] = MIN(MAX(1,W[i]),9);
	delete [] W;
	delete [] G;
	delete [] D;
	return sqrt(gg);
}

/*
 * This is a new particle swarm optimisation based solver...
 */
#define PARTICLES 100
double
LEbackcalc::swarm(unsigned nl, double * P)
{
	double gg = 0.0, dgg = 0.0, derr = DBL_MAX, werr = -DBL_MAX;
	//double r1 = 0.0;
	double r2 = 0.0;
	double sderr = sqrt(noise*noise + (precision*precision)/3.0);
	const double w = 0.0;
	//const double c1 = 0.2;
	//const double c2 = 0.3;
	const double G = 0.1;
	unsigned p, i, iter = 0;
	unsigned best = 0;

	double * X = new double[PARTICLES*nl];
	double * V = new double[PARTICLES*nl];
	double * B = new double[PARTICLES*nl];
	double * E = new double[2*PARTICLES];
	double * R1 = new double[nl*nl];
	double * R2 = new double[nl*nl];
	for (p = 0; p < PARTICLES; p++) {
		for (i = 0; i < nl; i++) {
			X[p*nl+i] = RAND(MAX(2,P[i]-6.0),MIN(P[i]+4.0,8));
			V[p*nl+i] = 0.0;
		}
		E[2*p] = DBL_MAX;
	}
	while (iter++ < 10000) {
		for (p = 0; p < PARTICLES; p++) {
			double err = 0.0;
			for (i = 0; i < nl; i++)
				layer(i).emod(pow(10,X[p*nl+i]));
			calculate(LEsystem::disp);
			for (i = 0; i < defl.length(); i++) {
				defldata & d = defl[i];
				d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
				err += pow(d.measured-d.calculated,2);
			}
			err = sqrt(err/defl.length());
			printf("%4i:",p);
			for (i = 0; i < nl; i++)
				printf(" %+4.2f",X[p*nl+i]);
			for (i = 0; i < nl; i++)
				printf(" %+4.2f",V[p*nl+i]);
			printf(" (%g)\n",log10(err));
			E[2*p+1] = log10(err);
			if (log10(err) < E[2*p]) {
				E[2*p] = log10(err);
				for (i = 0; i < nl; i++)
					B[p*nl+i] = X[p*nl+i];
			}
			if (err < derr)
				best = p, derr = err;
			if (err > werr)
				werr = err;
		}
		printf("%i: (%g)\n",best,derr);
		if (derr < MAX(1e-6,sderr))
			break;
		for (p = 0, dgg = 0.0; p < PARTICLES; p++) {
			/*
			for (i = 0; i < nl; i++) {
				R1[i*nl+i] = R2[i*nl+i] = 1.0;
				for (j = i+1; j < nl; j++) {
					r1 = RAND(-0.5,0.5); r2 = RAND(-0.5,0.5);
					R1[i*nl+j] = R1[j*nl+i] = 10*M_PI/180*(r1+r2);
					r1 = RAND(-0.5,0.5); r2 = RAND(-0.5,0.5);
					R2[i*nl+j] = R2[j*nl+i] = 10*M_PI/180*(r1+r2);
				}
			}
			orth_gs(nl,R1); orth_gs(nl,R2);
			r1 = 0.0; r2 = 0.0;
			for (i = 0; i < nl; i++) {
				r1 += pow(B[p*nl+i]-X[p*nl+i],2);
				r2 += pow(B[best*nl+i]-X[p*nl+i],2);
			}
			r1 = (r1 > 0.0 ? RAND(0.0,c1)/MAX(sqrt(r1),1.0) : 0.0);
			r2 = (r2 > 0.0 ? RAND(0.0,c2)/MAX(sqrt(r2),1.0) : 0.0);
			for (i = 0; i < nl; i++) {
				V[p*nl+i] *= w;
				for (j = 0; j < nl; j++)
					V[p*nl+i] += r1*R1[i*nl+j]*(B[p*nl+j]-X[p*nl+j]);
				for (j = 0; j < nl; j++)
					V[p*nl+i] += r2*R2[i*nl+j]*(B[best*nl+j]-X[p*nl+j]);
			}
			for (i = 0, gg = 0.0; i < nl; i++) {
				double x = MIN(MAX(1,X[p*nl+i]+V[p*nl+i]),9);
				V[p*nl+i] = x-X[p*nl+i], X[p*nl+i] = x;
				gg += pow(V[p*nl+i],2);
			}
			*/
			for (i = 0; i < nl; i++)
				V[p*nl+i] *= w;
			double m2 = log10(werr) - E[2*p+1] + 1.0;
			m2 = pow(m2*3/4/M_PI,1.0/3.0);
			for (unsigned q = 0; q < PARTICLES; q++) {
				if (E[2*q] >= E[2*p+1])
					continue;
				double m1 = log10(werr) - E[2*q] + 1.0;
				m1 = pow(m1*3/4/M_PI,1.0/3.0);
				printf("%g %g\n",m1,m2);
				for (i = 0, gg = 0.0; i < nl; i++)
					gg += pow(B[q*nl+i]-X[p*nl+i],2);
				r2 = sqrt(gg); gg = MAX(0.1,sqrt(gg));
				//gg = G*(E[2*p]-E[2*q])/pow(gg,3);
				gg = G*(E[2*p+1]-E[2*q])/pow(gg,2);
				for (i = 0; i < nl; i++)
					V[p*nl+i] += gg*(B[q*nl+i]-X[p*nl+i])/r2;
			}
			double r = 0.0;
			for (i = 0; i < nl; i++)
				r += pow(V[p*nl+i],2);
			r = MAX(1.0,sqrt(r));
			for (i = 0, gg = 0.0; i < nl; i++) {
				double x = MIN(MAX(1,X[p*nl+i]+V[p*nl+i]/r),9);
				V[p*nl+i] = x-X[p*nl+i], X[p*nl+i] = x;
				gg += pow(V[p*nl+i],2);
			}
			dgg += sqrt(gg)/PARTICLES;
		}
		printf("dgg = %g\n",dgg);
		if (dgg < MAX(1e-8,tolerance))
			break;
	}
	for (i = 0, gg = 0.0; i < nl; i++) {
		gg += pow(P[i]-B[best*nl+i],2);
		P[i] = MIN(MAX(1,B[best*nl+i]),9);
	}
	delete [] X;
	delete [] V;
	delete [] B;
	delete [] E;
	delete [] R1;
	delete [] R2;
	return sqrt(gg);
}
