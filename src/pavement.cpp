/*************************************************************************

	PAVEMENT.CPP - Implementation for the pavement class.

	$OpenPave$

	See PAVEMENT.H.

	History:
		2002/01/23 - Created by Jeremy Lea <jlea@csir.co.za>

*************************************************************************/

#include "config.h"
#include "pavement.h"
#include "matrix.h"
#include <time.h>
#include <stdio.h>

/*
 * Return the results based on a more rational system...
 */
double
pavedata::result(type t, direction d) const
{
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
		default:
			return 0.0;
		}
	case deflct:
		switch (d) {
		case xx: case yy: case zz:
			return data[4][d-xx];
		default:
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
		default:
			return 0.0;
		}
	default:
		return 0.0;
	}
}

/*
 * Calculate the principle stresses, and the strains.
 */
void
pavedata::principle(double v, double E)
{
	int k1, k2;
    double t1, t2;

	for (k1 = 0; k1 < 3; k1++) {
		data[2][k1] = data[0][k1];
		data[3][k1] = data[1][k1];
	}
	/*
	 * Find the eigen values of the stress tensor.  This is heavily
	 * optimised for a symmetric 3x3 matrix, to use as few FP registers
	 * as possible.
	 */
	t1 = data[3][0]*data[3][0] + data[3][1]*data[3][1]
			+ data[3][2]*data[3][2];
	while (t1 > 1e-30) {
		if (data[3][1] != 0.0) {
			t1 = data[2][2]-data[2][0];
			t2 = (t1 < 0.0 ? -2 : 2)*data[3][1];
			t1 = t2/(fabs(t1) + pythag(t1,t2));
			data[2][0] -= t1*data[3][1];
			data[2][2] += t1*data[3][1];
			t2 = pythag(t1,1.0);
			data[3][2] /= t2;
			t2 = data[3][0] /= t2;
			data[3][0] -= t1*data[3][2];
			data[3][2] += t1*t2;
			data[3][1] = 0.0;
		}
		if (data[3][0] != 0.0) {
			t1 = data[2][1]-data[2][0];
			t2 = (t1 < 0.0 ? -2 : 2)*data[3][0];
			t1 = t2/(fabs(t1) + pythag(t1,t2));
			data[2][0] -= t1*data[3][0];
			data[2][1] += t1*data[3][0];
			data[3][2] /= pythag(t1,1.0);
			data[3][1] = -t1*data[3][2];
			data[3][0] = 0.0;
		}
		if (data[3][2] != 0.0) {
			t1 = data[2][2]-data[2][1];
			t2 = (t1 < 0.0 ? -2 : 2)*data[3][2];
			t1 = t2/(fabs(t1) + pythag(t1,t2));
			data[2][1] -= t1*data[3][2];
			data[2][2] += t1*data[3][2];
			data[3][1] /= pythag(t1,1.0);
			data[3][0] = -t1*data[3][1];
			data[3][2] = 0.0;
		}
		t1 = data[3][0]*data[3][0] + data[3][1]*data[3][1];
	}
	data[3][0] = (data[2][0]-data[2][2])/2;
	data[3][1] = (data[2][0]-data[2][1])/2;
	data[3][2] = (data[2][1]-data[2][2])/2;
	for (k1 = 0; k1 < 2; k1++) {
		for (k2 = k1 + 1; k2 < 3; k2++) {
			if (data[2][k1] > data[2][k2]) {
				t1 = data[2][k1];
				data[2][k1] = data[2][k2];
				data[2][k2] = t1;
			}
			if (data[3][k1] > data[3][k2]) {
				t1 = data[3][k1];
				data[3][k1] = data[3][k2];
				data[3][k2] = t1;
			}
		}
	}

	t1 = 2*(v+1)/E;
	data[5][0] = (data[0][0]-v*(data[0][1]+data[0][2]))/E;
	data[5][1] = (data[0][1]-v*(data[0][2]+data[0][0]))/E;
	data[5][2] = (data[0][2]-v*(data[0][0]+data[0][1]))/E;
	data[6][0] = t1*data[1][0];
	data[6][1] = t1*data[1][1];
	data[6][2] = t1*data[1][2];

	data[7][0] = (data[2][0]-v*(data[2][1]+data[2][2]))/E;
	data[7][1] = (data[2][1]-v*(data[2][2]+data[2][0]))/E;
	data[7][2] = (data[2][2]-v*(data[2][0]+data[2][1]))/E;
	data[8][0] = t1*data[3][0];
	data[8][1] = t1*data[3][1];
	data[8][2] = t1*data[3][2];
}

bool
LEsystem::addlayer(double h, double e, const double v, const double s, const int p)
{
	LElayer * pl = new LElayer(this,(p < 0 ? last : &layer(p)),h,e,v,s);
	if (pl == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEsystem::addlayer()!");
		return false;
	}
	return true;
}

bool
LEsystem::removelayers()
{
	empty();
	return isempty();
}

bool
LEsystem::removelayer(const int l)
{
	LElayer * pl = first;
	int i = 0;

	if (l < 0)
		return false;
	while (i++ < l && pl != 0)
		pl = pl->next;
	if (pl != 0) {
		delete pl;
		return true;
	}
	return false;
};

bool
LEsystem::removeloads()
{
	return load.empty();
};

bool
LEsystem::removeload(const int i)
{
	return load.remove(i);
};

bool
LEsystem::addpoint(const point3d & p)
{
	return data.add(pavedata(p));
};
	
bool
LEsystem::addgrid(const int nx, const double *xp,
                  const int ny, const double *yp,
                  const int nz, const double *zp)
{
	pavedata * pd = new pavedata[nx*ny*nz];
	memset(pd,0,nx*ny*nz*sizeof(pavedata));
	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				pavedata & pde = pd[ix*ny*nz+iy*nz+iz];
				pde.x = xp[ix], pde.y = yp[iy], pde.z = zp[iz];
			}
		}
	}
	return data.add(pd,nx*ny*nz);
};

bool
LEsystem::removepoints()
{
	return data.empty();
};

bool
LEsystem::removepoint(const point3d & p)
{
	return data.remove(p);
};

LElayer &
LEsystem::layer(const int l)
{
	LElayer * pl = first;
	int i = 0;
	while (i++ < l && pl->next != 0)
		pl = pl->next;
	return *pl;
};

/*
 * This checks the structure of the pavement to ensure that is is
 * good...
 */
bool
LEsystem::check()
{
    int il, ixy;
	const LElayer * pl;

	if (layers() <= 0) {
		event_msg(EVENT_WARN,
			"Cannot calculate a pavement without any layers!");
		return false;
	}
	if (load.length() <= 0) {
		event_msg(EVENT_WARN,
			"Cannot calculate a pavement without any loads!");
		return false;
	}
	if (data.length() <= 0) {
		event_msg(EVENT_WARN,
			"Cannot calculate a pavement without any evaluation points!");
		return false;
	}
	for (il = 1, pl = first; pl != 0; il++, pl = pl->next) {
		if (pl->emod() <= 0.0) {
			event_msg(EVENT_WARN,
				"Error: Elastic modulus of layer %d must be greater"
				" than zero not %f!", il, pl->emod());
			return false;
		}
		if (pl->poissons() <= 0.0 || pl->poissons() > 0.5) {
			event_msg(EVENT_WARN,
				"Error: Poisson's ratio of layer %d must be between"
				" zero and one half not %f!", il, pl->poissons());
			return false;
		}
		if (pl->thickness() < 0.0) {
			event_msg(EVENT_WARN,
				"Error: Layer %d cannot have negative thickness!", il);
			return false;
		}
		if (pl->slip() < 0.0 || pl->slip() > 1.0 ) {
			event_msg(EVENT_WARN,
				"Error: Layer %d has an invalid bonding coefficent!", il);
			return false;
		}
		if (pl->thickness() == 0.0 && pl->next != 0) {
			event_msg(EVENT_WARN,
				"Error: Layer %d cannot have zero thickness!", il);
			return false;
		}
		if (pl->thickness() == 0.0 && pl->next == 0 && pl->slip() != 1.0) {
			event_msg(EVENT_WARN,
				"Error: Infinite layer %d cannot have imperfect bonding!", il);
			return false;
		}
		if (pl->next != 0 && pl->slip() < 0.0) {
			event_msg(EVENT_WARN,
				"Error: Layer %d's bonding coefficent will be clamped to 0.0!", il);
		}
	}
	if (last->bottom() > 0.0  && last->poissons() == 0.75) {
 		event_msg(EVENT_WARN,
			"Error: Last layer cannot have a Poisson's ratio of 0.75!");
		return false;
	}
	for (ixy = 0; ixy < data.length(); ixy++) {
		// Ignore evaluation points above the surface
		// or below the rigid interface.
		if (data[ixy].z < 0.0 || (last->bottom() > 0.0
		 && last->bottom() <= data[ixy].z)) {
 			event_msg(EVENT_WARN,
				"Error: evaluation point %d (%f,%f,%f) not within the pavement!",
				ixy,data[ixy].x,data[ixy].y,data[ixy].z);
			return false;
		}
	}
	return true;
};

/*
 * As simple struct to store the axisymetric data.  These are
 * stored in the ax set, indexed by the radius and depth.
 */
struct axialdata {
	double rse;
	double tse;
	double vse;
	double sse;
	double rdp;
	double vdp;
};
#define GRADSTEP 1e-8

#define NBZ		8192						// Number of zeros in bessels.
static double j0r[NBZ+1], j1r[NBZ+1], j0p[NBZ+1], j0m[NBZ+1];
#define j0pj1(x)	(j0(x)+j1(x))
#define j0mj1(x)	(j0(x)-j1(x))
#define j1d(x)		(j0(x)-j1(x)/x)

#define NGQP	16						// Max number of Gauss points
static double gu[NGQP+1][NGQP];
static double gf[NGQP+1][NGQP];

#define NLQP	16						// Max number of Lobatto points
static double lu[NLQP+1][NLQP];
static double lf[NLQP+1][NLQP];

/*
 * This initialises the arrays above to max accuracy, and
 * find the zeros of the various bessel functions...
 */
static void
initarrays() {
    int ib;
	static bool done = false;

	if (done)
		return;
	memset(gu,0,sizeof(double)*(NGQP+1)*NGQP);
	memset(gf,0,sizeof(double)*(NGQP+1)*NGQP);
	for (ib = 1; ib <= NGQP; ib++) {
		for (int j, i = 0; i < (ib+1)/2; i++) {
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
	memset(lu,0,sizeof(double)*(NLQP+1)*NLQP);
	memset(lf,0,sizeof(double)*(NLQP+1)*NLQP);
	for (ib = 2; ib <= NLQP; ib++) {
		for (int j, i = 0; i < (ib+1)/2; i++) {
			double p1, p2, p3, z = cos(M_PI*i/(ib-1));
			do {
				p2 = 1.0, p1 = z;
				for (j = 1; j < (ib-1); j++)
					p3 = p2, p2 = p1, p1 = z*p2+j*(z*p2-p3)/(j+1);
				z += (p2/p1-z)/ib;
			} while (fabs(p2/p1-z)/ib > DBL_EPSILON);
			lu[ib][i] = -z, lu[ib][ib-1-i] = z;
			lf[ib][i] =     lf[ib][ib-1-i] = 2.0/(ib*(ib-1)*p1*p1);
		}
	}
	memset(j0r,0,sizeof(double)*(NBZ+1));
	memset(j1r,0,sizeof(double)*(NBZ+1));
	memset(j0p,0,sizeof(double)*(NBZ+1));
	memset(j0m,0,sizeof(double)*(NBZ+1));
	for (ib = 0; ib <= NBZ; ib++) {
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
		};
		j1r[ib] = x2;
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
stoppingpoints(const int nbz, const double a, const double r,
               double * m0, double * m1)
{
	int ib;
	double ra = r/a, r1 = fabs(ra-1);

	// At r=0 or r=a just use the last root.
	*m0 = j0r[nbz]/a, *m1 = (ra == 0.0 ? 0.0 : j0r[nbz]/a);
	if (ra > 0.0 && ra <= 0.5) {
		// Find the last root of J1(mr) & J0(mr).
		for (ib = nbz; ib > 0 && j1r[ib]/ra > j1r[nbz]; ib--)
			*m0 = j1r[ib-1]/ra;
		for (ib = nbz; ib > 0 && j0r[ib]/ra > j1r[nbz]; ib--)
			*m1 = j0r[ib-1]/ra;
	} else if (ra > 0.5 && ra != 1.0 && ra < 2.0) {
		// Between 0.5 and 2.0 use the J0(m|r-1|)-J1(m|r-1|) and
		// J0(m|r-1|)+J1(m|r-1|) approximations.
		for (ib = nbz; ib > 0 && j0m[ib]/r1 >
					(ra < 1.0 ? j1r[nbz] : j0r[nbz]/ra); ib--)
			*m0 = j0m[ib-1]/r1;
		for (ib = nbz; ib > 0 && j0p[ib]/r1 >
					(ra < 1.0 ? j1r[nbz] : j1r[nbz]/ra); ib--)
			*m1 = j0p[ib-1]/r1;
	} else if (ra >= 2.0) {
		// Find the roots of J0(m).
		for (ib = nbz; ib > 0 && j0r[ib] > j0r[nbz]/ra; ib--)
			*m0 = j0r[ib-1];
		for (ib = nbz; ib > 0 && j0r[ib] > j1r[nbz]/ra; ib--)
			*m1 = j0r[ib-1];
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
 * This build the ABCD matrix, based on the structure, with slip.  This has
 * to be done by solving the entire system, since the piecewise solution
 * used above is singular.
 */
static void
buildabcd_full(const double m, const int nl, const double * h,
			   const double * v, const double * E, const double * f,
               double (* ABCD)[4])
{
	int k1, k2, il = (nl-1);

	memset(ABCD,0,(nl*4)*sizeof(double));
	if (m <= 0.0)
		return;
	// First allocate space for the arrays...
	double * A = new double[(nl*4)*(nl*4)];
	double * B = &ABCD[0][0];
	if (A == 0) {
		event_msg(EVENT_ERROR,"Out of memory in buildabcd_slip()!");
		goto abort;
	}
	memset(A,0,(nl*4)*(nl*4)*sizeof(double));
	memset(B,0,(nl*4)*sizeof(double));
	// We start with the last layer...
	if (h[il] > 0.0) {
		double h1 = h[il];
		double v1 = v[il];
		double s1 = f[il];
		double t1 = (m*E[il]*(s1-1))/(1+v1);
		A[(il*4+2-2)*(nl*4)+il*4+0]  = m*t1 + m*s1;
		A[(il*4+2-2)*(nl*4)+il*4+1]  = (m*h1+2*v1)*t1 + s1*(1+m*h1);
		A[(il*4+2-2)*(nl*4)+il*4+2]  = m*t1 - m*s1;
		A[(il*4+2-2)*(nl*4)+il*4+2] *= exp(-2*m*h1);
		A[(il*4+2-2)*(nl*4)+il*4+3]  = (m*h1-2*v1)*t1 + s1*(1-m*h1);
		A[(il*4+2-2)*(nl*4)+il*4+3] *= exp(-2*m*h1);
		A[(il*4+3-2)*(nl*4)+il*4+0]  = m;
		A[(il*4+3-2)*(nl*4)+il*4+1]  = m*h1+2*(2*v1-1);
		A[(il*4+3-2)*(nl*4)+il*4+2]  = m;
		A[(il*4+3-2)*(nl*4)+il*4+2] *= exp(-2*m*h1);
		A[(il*4+3-2)*(nl*4)+il*4+3]  = m*h1-2*(2*v1-1);
		A[(il*4+3-2)*(nl*4)+il*4+3] *= exp(-2*m*h1);
	} else {
		A[(il*4+2-2)*(nl*4)+il*4+0] = 1.0;
		A[(il*4+3-2)*(nl*4)+il*4+1] = 1.0;
	}
	// Now we work back up, building the two 4x4 matrices...
	for ( ; il > 0; il--) {
		double h1 = h[il-1];
		double s1 = f[il-1];
		double v1 = v[il-1];
		double v2 = v[il];
		double K = (1+v2)/(1+v1)*E[il-1]/E[il];
		double R = E[il-1]*m*(s1-1)/(1+v1);
		double t1 = m*h1;
		A[(il*4-4)*(nl*4)+il*4-4] =  m;
		A[(il*4-4)*(nl*4)+il*4-3] =  m*h1+2*v1-1;
		A[(il*4-4)*(nl*4)+il*4-2] = -m;
		A[(il*4-4)*(nl*4)+il*4-1] = -m*h1+2*v1-1;
		A[(il*4-4)*(nl*4)+il*4+0] =  m;
		A[(il*4-4)*(nl*4)+il*4+1] =  m*h1+2*v2-1;
		A[(il*4-4)*(nl*4)+il*4+2] = -m;
		A[(il*4-4)*(nl*4)+il*4+3] = -m*h1+2*v2-1;
		A[(il*4-3)*(nl*4)+il*4-4] =  m;
		A[(il*4-3)*(nl*4)+il*4-3] =  m*h1+2*v1;
		A[(il*4-3)*(nl*4)+il*4-2] =  m;
		A[(il*4-3)*(nl*4)+il*4-1] =  m*h1-2*v1;
		A[(il*4-3)*(nl*4)+il*4+0] =  m;
		A[(il*4-3)*(nl*4)+il*4+1] =  m*h1+2*v2;
		A[(il*4-3)*(nl*4)+il*4+2] =  m;
		A[(il*4-3)*(nl*4)+il*4+3] =  m*h1-2*v2;
		A[(il*4-2)*(nl*4)+il*4-4] =  m*(R+s1);
		A[(il*4-2)*(nl*4)+il*4-3] =  R*(m*h1+2*v1)+s1*(1+m*h1);
		A[(il*4-2)*(nl*4)+il*4-2] =  m*(R-s1);
		A[(il*4-2)*(nl*4)+il*4-1] =  R*(m*h1-2*v1)+s1*(1-m*h1);
		A[(il*4-2)*(nl*4)+il*4+0] =  m*K*s1;
		A[(il*4-2)*(nl*4)+il*4+1] = (1+m*h1)*K*s1;
		A[(il*4-2)*(nl*4)+il*4+2] = -m*K*s1;
		A[(il*4-2)*(nl*4)+il*4+3] = (1-m*h1)*K*s1;
		A[(il*4-1)*(nl*4)+il*4-4] =  m;
		A[(il*4-1)*(nl*4)+il*4-3] =  m*h1+4*v1-2;
		A[(il*4-1)*(nl*4)+il*4-2] =  m;
		A[(il*4-1)*(nl*4)+il*4-1] =  m*h1-4*v1+2;
		A[(il*4-1)*(nl*4)+il*4+0] =  m*K;
		A[(il*4-1)*(nl*4)+il*4+1] = (m*h1+4*v2-2)*K;
		A[(il*4-1)*(nl*4)+il*4+2] =  m*K;
		A[(il*4-1)*(nl*4)+il*4+3] = (m*h1-4*v2+2)*K;
		for (k1 = 0; k1 < 2; k1++) {
			for (k2 = 0; k2 < 2; k2++) {
				A[(il*4-4+k1)*(nl*4)+il*4-2+k2] *=  exp(-2*t1);
				A[(il*4-2+k1)*(nl*4)+il*4-4+k2] *=  (2*t1<MAX_EXP?exp(2*t1):DBL_MAX);
				A[(il*4-4+k1)*(nl*4)+il*4+0+k2] *= -1.0;
				A[(il*4-4+k1)*(nl*4)+il*4+2+k2] *= -exp(-2*t1);
				A[(il*4-2+k1)*(nl*4)+il*4+0+k2] *= -(2*t1<MAX_EXP?exp(2*t1):DBL_MAX);
				A[(il*4-2+k1)*(nl*4)+il*4+2+k2] *= -1.0;
			}
		}
	}
	// Then we fill in the first final values...
	A[(nl*4-2)*(nl*4)+0] =  m;
	A[(nl*4-2)*(nl*4)+1] =  2*v[0];
	A[(nl*4-2)*(nl*4)+2] =  m;
	A[(nl*4-2)*(nl*4)+3] = -2*v[0];
	A[(nl*4-1)*(nl*4)+0] =  m;
	A[(nl*4-1)*(nl*4)+1] =  2*v[0]-1;
	A[(nl*4-1)*(nl*4)+2] = -m;
	A[(nl*4-1)*(nl*4)+3] =  2*v[0]-1;
	B[(nl*4-1)] = 1/m/m;
	//for (k1 = 0; k1 < nl*4; k1++) {
	//	for (k2 = 0; k2 < nl*4; k2++) {
	//		printf("%g",A[k1*(nl*4)+k2]);
	//		if (k2 == nl*4-1)
	//			printf("\n");
	//		else
	//			printf("\t");
	//	}
	//}
	//printf("\n");

	// Solve the system using Gaussian elimination with full pivoting.
	// this is slow, but stable - even LU decomposition is not stable enough
	// to ensure good results...
	for (il = 0; il < nl*4; il++) {
		double pvt = A[il*(nl*4)+il];
		if (fabs(pvt) < DBL_EPSILON) {
			for (k1 = il+1; k1 < nl*4; k1++) {
				if (fabs(pvt = A[k1*(nl*4)+il]) >= DBL_EPSILON)
					break;
			}
			if (k1 == nl*4)
				break; // XXX
			for (k2 = 0; k2 < nl*4; k2++)
				swap(A[k1*(nl*4)+k2],A[il*(nl*4)+k2]);
			swap(B[k1],B[il]);
		}
		for (k2 = nl*4-1; k2 > il; k2--) {
			double tmp = A[k2*(nl*4)+il]/pvt;
			for (k1 = nl*4-1; k1 > il; k1--)
				A[k2*(nl*4)+k1] -= tmp*A[il*(nl*4)+k1];
			B[k2] -= tmp*B[il];
		}
	}
	for (il = nl*4-1; il >= 0; il--) {
		for (k1 = nl*4 - 1; k1 > il; k1--)
			B[il] -= A[il*(nl*4)+k1]*B[k1];
		B[il] /= A[il*(nl*4)+il];
	}
	for (il = 0; il < nl; il++) {
		printf("%g\t%g\t%g\t%0.16g\t%0.16g\t%0.16g\t%0.16g\n",m,h[il],E[il],ABCD[il][0],ABCD[il][1],ABCD[il][2],ABCD[il][3]);
	}
abort:
	delete [] A;
}

/*
 * This build the ABCD matrix, based on the structure.
 */
static void
buildabcd(const double m, const int nl, const double * h, const double * v,
          const double * E, const double * f, double (* R)[4][2],
          double (* ABCD)[4])
{
	int k1, k2, il;
	double B1[2][4], X[4][4], F[4][4], D[4][4];
	double CDi[4] = {0, 0, 0, 0};

	memset(ABCD,0,sizeof(double)*nl*4);
	if (m <= 0.0)
		return;
	//if (slip) {
	//	buildabcd_full(m,nl,h,v,E,f,ABCD);
	//	return;
	//}
	// We start at the last layer...
	il = nl-1;
	if (h[il] > 0.0) {
		double z = h[il];
		double s = f[il];
		double v1 = v[il];
		double t1 = 2*m*E[il]*(1-s)*(v1-1)/(v1+1);
		R[il][0][0] =   s*(1-4*v1-2*m*z) - t1;
		R[il][0][1] = 4*s*(2*v1-1)/m - 2*s*m*z*z - 2*z*t1;
		R[il][1][0] = 2*s*m;
		R[il][1][1] =   s*(1-4*v1+2*m*z) + t1;
		t1 += s*(3-4*v1);
		for (k1 = 0; k1 < 2; k1++)
			for (k2 = 0; k2 < 2; k2++)
				R[il][k1][k2] *= exp(-2*m*h[il])/t1;
	} else {
		memset(&R[il][0][0],0,4*sizeof(double));
	}
	R[il][2][0] = 1.0, R[il][2][1] = 0.0;
	R[il][3][0] = 0.0, R[il][3][1] = 1.0;
	// Now we work back up, building the 4x2 R matrices...
	for ( ; il > 0; il--) {
		double z = h[il-1];
		double v1 = v[il-1];
		double v2 = v[il];
		double K = (1+v2)/(1+v1)*E[il-1]/E[il];
		double t1 = m*z, t2 = 2*t1;
		X[0][0] = X[2][2] = (4*v1-3-K);
		X[0][1] = ((1+8*v1*v2-t2)*(1-K) + (4*t1*v1-6*v2)
			- K*(4*t1*v2-6*v1))/m;
		X[2][3] = ((1+8*v1*v2+t2)*(K-1) + (4*t1*v1+6*v2)
			- K*(4*t1*v2+6*v1))/m;
		X[1][0] = X[3][2] = 0.0;
		X[1][1] = X[3][3] = (K*(4*v2-3)-1);
		X[0][2] = ( t2+4*v1-1)*(1-K);
		X[2][0] = (-t2+4*v1-1)*(1-K);
		X[1][3] = (-t2+4*v2-1)*(1-K);
		X[3][1] = ( t2+4*v2-1)*(1-K);
		X[0][3] = ((1+2*(t1+2*v1)*(t1-2*v2))*(1-K) + 2*v2
			- 2*K*v1)/m;
		X[2][1] = ((1+2*(t1-2*v1)*(t1+2*v2))*(K-1) - 2*v2
			+ 2*K*v1)/m;
		X[1][2] = 2*m*(K-1);
		X[3][0] = 2*m*(1-K);
		if (f[il-1] != 1.0) {
			double s = E[il-1]*(1-f[il-1])/(1+v1);
			double r1 = (t1+2*v1-1)*s;
			double r2 = (t1-2*v1+1)*s;
			double r3 =  t1+2*v1;
			double r4 =  t1-2*v1;
			F[0][0] = r1;   F[0][1] = r1*r3/m;
			F[0][2] = r1;   F[0][3] = r1*r4/m;
			F[1][0] = -s*m; F[1][1] = -s*r3;
			F[1][2] = -s*m; F[1][3] = -s*r4;
			F[2][0] = r2;   F[2][1] = r2*r3/m;
			F[2][2] = r2;   F[2][3] = r2*r4/m;
			F[3][0] = -s*m; F[3][1] = -s*r3;
			F[3][2] = -s*m; F[3][3] = -s*r4;
		} else
			memset(F,0,4*4*sizeof(double));
		for (k1 = 0; k1 < 2; k1++) {
			for (k2 = 0; k2 < 2; k2++) {
				D[k1][k2] = D[k1+2][k2+2] = 1.0;
				D[k1][k2+2] = exp(-t2);
				D[k1+2][k2] = (t2<MAX_EXP?exp(t2):DBL_MAX);
			}
		}
		double t3 = 4*v1-4;
		for (k1 = 0; k1 < 4; k1++) {
			R[il-1][k1][0] = R[il-1][k1][1] = 0.0;
			for (k2 = 0; k2 < 4; k2++) {
				if (f[il-1] == 1.0) {
					R[il-1][k1][0] += D[k1][k2]*(X[k1][k2]*R[il][k2][0])/t3;
					R[il-1][k1][1] += D[k1][k2]*(X[k1][k2]*R[il][k2][1])/t3;
				} else {
					R[il-1][k1][0] += D[k1][k2]*(f[il-1]*X[k1][k2]*R[il][k2][0]/m+F[k1][k2]*R[il][k2][0])/t3/f[il-1];
					R[il-1][k1][1] += D[k1][k2]*(f[il-1]*X[k1][k2]*R[il][k2][1]/m+F[k1][k2]*R[il][k2][1])/t3/f[il-1];
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
	for (k1 = 0; k1 < 4; k1++) {
		CDi[0] += B1[0][k1]*R[0][k1][0];
		CDi[1] += B1[0][k1]*R[0][k1][1];
		CDi[2] += B1[1][k1]*R[0][k1][0];
		CDi[3] += B1[1][k1]*R[0][k1][1];
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
 * This is the 'correct' code, according to exactly how the
 * math should be done.  As a result, it is slow...  It is only
 * really here for checking the optimised code below, but it
 * is exposed so that people can use it if they really want.
 */
bool
LEsystem::accurate()
{
    int ixy, ild, ib, igp, il;
	const LElayer * pl;
	bool rv = true;
	
	initarrays();
	if (!check())
		return false;
	int nl = layers();

	// The integration constants, per layer.
	double (* R)[4][2] = new double[nl][4][2];
	double (* ABCD)[4] = new double[nl][4];
	double * h = new double[nl];
	double * f = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	// We allocate these as big as we ever make them,
	// so we never have to worry about the add()'s failing.
	cset<double> * _bm0 = new cset<double>(0, 2*NBZ+2);
	cset<double> & bm0 = *_bm0;
	cset<double> * _bm1 = new cset<double>(0, 2*NBZ+2);
	cset<double> & bm1 = *_bm1;
	if (R == 0 || ABCD == 0 || h == 0 || f == 0 || v == 0 || E == 0
	 || _bm0 == 0 || _bm1 == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEsystem::calculate()!");
		rv =  false;
		goto abort;
	}
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
		h[il] = pl->bottom();
		f[il] = MAX(0.0,pl->slip());
		v[il] = pl->poissons();
		E[il] = pl->emod();
	}

	// We loop through all of the evaluation points
	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		// Zero the data...
		memset(d.data,0,sizeof(d.data));
		// Ignore evaluation points above the surface or below
		// the rigid interface.
		if (d.z < 0.0 || (last->bottom() > 0.0 && last->bottom() <= d.z))
			continue;
		// Find our layer...
		for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
			if (pl->top() <= d.z && (h[il] == 0.0 || d.z < h[il]))
				break;
		}
		// Now loop through the list of loads...
		for (ild = 0; ild < load.length(); ild++) {
			double a = load[ild].radius();
			double p = load[ild].pressure();
			double r = load[ild].distance(d);
			axialdata s;
			double m0, m1;
			memset(&s,0,sizeof(axialdata));

			// Now gerenate a list of integration intervals, then sort them.
			bm0.empty(), bm1.empty();
			bm0.add(0.0), bm1.add(0.0);
			stoppingpoints(NBZ,a,r,&m0,&m1);
			for (ib = 0; r > 0 && j0r[ib]/r < m0 && ib <= NBZ; ib++)
				bm0.add(j0r[ib]/r);
			for (ib = 0; r > 0 && j1r[ib]/r < m1 && ib <= NBZ; ib++)
				bm1.add(j1r[ib]/r);
			for (ib = 0; j1r[ib]/a < m0 && ib <= NBZ; ib++)
				bm0.add(j1r[ib]/a);
			for (ib = 0; j1r[ib]/a < m1 && ib <= NBZ; ib++)
				bm1.add(j1r[ib]/a);
			bm0.add(m0), bm1.add(m1);
			bm0.sort(), bm1.sort();
			// Now remove every second zero.  We allow the sort
			// to remove these, since it does things in one pass.
			for (ib = 1; bm0[ib] < m0 && ib < bm0.length()-1; ib += 2)
				bm0[ib] = bm0[ib-1];
			for (ib = 1; bm1[ib] < m1 && ib < bm1.length()-1; ib += 2)
				bm1[ib] = bm1[ib-1];
			bm0.sort(), bm1.sort();

			// We loop through all of our roots and gauss points.
			// In this version we do this in two parts, one for the
			// J1(ma)*J0(mr) integrals and one for the J1(ma)*J1(mr).
			// This is so we can use Guass-Lebatto, with exact end
			// points, which gives us a higher accuracy integral.
			// XXX: Except that Guass-Lebatto is not accurate...
			for (ib = bm0.length()-1; ib > 0; ib--) {
				double vse = 0.0, rse = 0.0, tse = 0.0, vdp = 0.0;
				for (igp = 0; igp < NGQP; igp++) {
					//if (ib >= 2 && igp == NLQP-1 && bm0[ib] != m0)
					//	continue;
					// Calculate the gauss point and weight.
					double m = (bm0[ib]+bm0[ib-1])/2
							 + gu[NGQP][igp]*(bm0[ib]-bm0[ib-1])/2;
					double w = gf[NGQP][igp]*(bm0[ib]-bm0[ib-1])/2;
					// First build a new ABCD matrix.
					buildabcd(m,nl,h,v,E,f,R,ABCD);
					double t1, t3, t4, t5, t6;
					t1 = m*j1(m*a)*w*j0(m*r);
					t3 = m*(ABCD[il][2] + ABCD[il][3]*d.z)*exp(-m*d.z);
					t4 = ABCD[il][3]*exp(-m*d.z);
					t5 = m*(ABCD[il][0] + ABCD[il][1]*d.z)*
									 (m*d.z<MAX_EXP?exp(m*d.z):DBL_MAX);
					t6 = ABCD[il][1]*(m*d.z<MAX_EXP?exp(m*d.z):DBL_MAX);
					vse += t1*m*((1-2*v[il])*(t4+t6)+(t3-t5));
					rse += t1*m*((1+2*v[il])*(t4+t6)-(t3-t5));
					tse += t1*m*(2*v[il])*(t4+t6);
					vdp += t1*((2-4*v[il])*(t6-t4)-(t3+t5));
				}
				s.vse += vse, s.rse += rse, s.tse += tse, s.vdp += vdp;
			}
			for (ib = 1; ib < bm1.length(); ib++) {
				double sse = 0.0, rdp = 0.0; 
				for (igp = 0; igp < NGQP; igp++) {
					// Calculate the gauss point and weight.
					double m = (bm1[ib]+bm1[ib-1])/2
							 + gu[NGQP][igp]*(bm1[ib]-bm1[ib-1])/2;
					double w = gf[NGQP][igp]*(bm1[ib]-bm1[ib-1])/2;
					// First build a new ABCD matrix.
					buildabcd(m,nl,h,v,E,f,R,ABCD);
					double t1, t3, t4, t5, t6;
					t1 = m*j1(m*a)*j1(m*r)*w;
					t3 = m*(ABCD[il][2] + ABCD[il][3]*d.z)*exp(-m*d.z);
					t4 = ABCD[il][3]*exp(-m*d.z);
					t5 = m*(ABCD[il][0] + ABCD[il][1]*d.z)*
									 (m*d.z<MAX_EXP?exp(m*d.z):DBL_MAX);
					t6 = ABCD[il][1]*(m*d.z<MAX_EXP?exp(m*d.z):DBL_MAX);
					sse += t1*m*((2*v[il])*(t6-t4) + (t3+t5));
					rdp += t1*((t4+t6)-(t3-t5));
				}
				s.sse += sse, s.rdp += rdp;
			}
			if (r > 0.0) {
				s.rse -= s.rdp/r;
				s.tse += s.rdp/r;
			}
			s.rse *= a; s.tse *= a;
			s.vse *= a; s.sse *= a;
			s.rdp *= a*(v[il]+1)/E[il];
			s.vdp *= a*(v[il]+1)/E[il];
			// Now add this load's axial data into the point's data
			if (r == 0.0) {
				d.data[0][0] += p*(s.rse+s.tse)/2;
				d.data[0][1] += p*(s.rse+s.tse)/2;
				d.data[0][2] += p*s.vse;
				d.data[1][1] += p*s.sse;
			} else {
				double cost, sint;
				sint = (d.y-load[ild].y)/r;
				cost = (d.x-load[ild].x)/r;
				d.data[0][0] += p*(cost*cost*s.rse+sint*sint*s.tse);
				d.data[0][1] += p*(cost*cost*s.tse+sint*sint*s.rse);
				d.data[0][2] += p*s.vse;
				d.data[1][0] += p*cost*sint*(s.rse-s.tse);
				d.data[1][1] += p*cost*s.sse;
				d.data[1][2] += p*sint*s.sse;
				d.data[4][0] += p*cost*s.rdp;
				d.data[4][1] += p*sint*s.rdp;
			}
			d.data[4][2] += p*s.vdp;
		}
		// Calculate the derived results (principal stresses and strains).
		d.principle(v[il],E[il]);
    }
abort:
	delete [] R;
	delete [] ABCD;
	delete [] h;
	delete [] f;
	delete [] v;
	delete [] E;
	delete _bm0;
	delete _bm1;
	return rv;
};

/*
 * This used to be ELSYM5M, now it's a NxNxN layered elastic code...
 */
bool
LEsystem::calculate(resulttype result, double * Q)
{
    int ixy, ir, iz, ild, ia, ib, igp, il;
	const LElayer * pl;
	double x1, x2;
	bool rv = false;
	
	initarrays();
	if (!check())
		return false;
	int ngqp = NGQP, nbz = NBZ, gl = -1, nl = layers();
	if (result == fast || result == fastdisp || result == fastgrad) {
		ngqp = MIN(MAX(8,NGQP),8);
		nbz = MIN(NBZ,128);
		result = (result == fast ? all :
						result == fastdisp ? disp : dispgrad);
	} else {
		ngqp = MIN(MAX(8,NGQP),12);
		nbz = MIN(NBZ,512);
	}
	callcount++;
	if (result == dispgrad)
		callcount += nl;

	// The integration constants, per layer.
	double (* R)[4][2] = new double[nl][4][2];
	double (* ABCD)[4] = new double[nl][4];
	// Local variables, so we don't have to walk the list.
	double * h = new double[nl];
	double * f = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	// Some place to store our data...
	cset<double> * _z = new cset<double>(data.length());
	cset<double> * _a = new cset<double>(load.length());
	cset<double> & z = *_z;
	cset<double> & a = *_a;
	int * zl = new int[data.length()];
	cset<double> * _r = new cset<double>(0,data.length());
	cset<double> * _bm = new cset<double>(0,4*nbz);
	cset<double> & r = *_r;
	cset<double> & bm = *_bm;
	// And finally, somewhere to stick the radial data...
	axialdata * ax = new axialdata[data.length()*data.length()];
	double * m0 = new double[data.length()];
	double * m1 = new double[data.length()];
	if (R == 0 || ABCD == 0 || h == 0 || f == 0 || v == 0 || E == 0
	 || _z == 0 || _a == 0 || zl == 0 || _r == 0 || _bm == 0
	 || ax == 0 || m0 == 0 || m1 == 0)
		goto abort;
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
		h[il] = pl->bottom();
		f[il] = MAX(0.0,pl->slip());
		v[il] = pl->poissons();
		E[il] = pl->emod();
	}

	// Collect and sort the z positions.
	for (ixy = 0; ixy < data.length(); ixy++) {
		// Zero the data...
		memset(data[ixy].data,0,sizeof(data[ixy].data));
		// If we're collecting displacement gradient results
		// resize the array as needed then zero it.
		// We don't use memset since it's not a real array.
		if (result == dispgrad) {
			if (data[ixy].deflgrad.length() != nl) {
				data[ixy].deflgrad.empty();
				data[ixy].deflgrad.add(1,h,nl);
			}
			for (il = 0; il < nl; il++)
				data[ixy].deflgrad[il] = 0.0;
		}
		z[ixy] = data[ixy].z;
	}
	z.sort();
	// Map z values to layers.
	for (iz = 0; iz < z.length(); iz++) {
		for (pl = first, il = 0; pl != 0; pl = pl->next, il++) {
			if (pl->top() <= z[iz] && (h[il] == 0.0 || z[iz] < h[il]))
				zl[iz] = il;
		}
	}

	// Gerenate a list of load radii, then sort them and map from loads.
	for (ild = 0; ild < load.length(); ild++)
		a[ild] = load[ild].radius();
	a.sort();

	// Now loop through the list of load radii, calculating only for
	// the applicable loads... (load[ild].radius() == a[ia])
	for (ia = 0; ia < a.length(); ia++) {
		// Gerenate a list of radii, then sort them.
		r.empty();
		for (ild = 0; ild < load.length(); ild++) {
			if (load[ild].radius() != a[ia])
				continue;
			for (ixy = 0; ixy < data.length(); ixy++) {
				if (!r.add(load[ild].distance(data[ixy])))
					goto abort;
			}
		}
		r.sort();

		// Now gerenate a list of integration intervals, then sort them.
		bm.empty();
		if (!bm.add(0.0))
			goto abort;
		for (ib = 0; ib <= nbz; ib++) {
			if (!bm.add(j0r[ib]/a[ia]))
				goto abort;
		}

		// The correct stopping points for each radius.
		for (ir = 0; ir < r.length(); ir++) {
			stoppingpoints(nbz,a[ia],r[ir],&m0[ir],&m1[ir]);
			if (!bm.add(m0[ir]) || !bm.add(m1[ir]))
				goto abort;
		}
		bm.sort();

		// Account for big r's by adding extra integration intervals...
		x1 = 0.0, x2 = 0.0;
		for (ir = r.length()-1; ir >= 0 && r[ir]/a[ia] > ngqp-6; ir--) {
			for (int k1 = ngqp-6; k1 <= nbz; k1 += ngqp-6) {
				if ((x1 = j1r[k1]/r[ir]) < x2)
					continue;
				for (ib = 1; ib < bm.length() && bm[ib] < x1; ib++)
					;
				if (MIN(x1-bm[ib-1],bm[ib]-x1) < (ngqp-6)*M_PI_4/r[ir])
					continue;
				if (!bm.add(ib+1,x2 = x1))
					goto abort;
			}
		}
		// Account for big z's.  We drop approximately three orders of
		// magnitude for exp(-7).  Add 5 intervals, so we drop 15 orders
		// of magnitude.
		x1 = 0.0, x2 = 0.0;
		for (iz = z.length()-1; iz >= 0; iz--) {
			for (int k1 = 1; k1 <= 5 && k1*7*a[ia] < j0r[0]*z[iz]; k1++) {
				if ((x1 = k1*7*a[ia]/z[iz]) < x2)
					continue;
				for (ib = 1; ib < bm.length() && bm[ib] < x1; ib++)
					;
				if (MIN(x1-bm[ib-1],bm[ib]-x1) < 5*a[ia]/z[iz])
					continue;
				if (!bm.add(ib+1,x2 = x1))
					goto abort;
			}
		}
		bm.sort();

gradloop:
		memset(ax,0,sizeof(axialdata)*r.length()*z.length());
		// Now that we know the radii, get down to work.
		// We loop through all of our roots and gauss points.
		for (ib = 1; ib < bm.length(); ib++) {
			for (igp = 0; igp < ngqp; igp++) {
				// Calculate the gauss point and weight.
				double m = (bm[ib]+bm[ib-1])/2
						 + gu[ngqp][igp]*(bm[ib]-bm[ib-1])/2;
				double w = gf[ngqp][igp]*(bm[ib]-bm[ib-1])/2;

				// First build a new ABCD matrix.
				buildabcd(m,nl,h,v,E,f,R,ABCD);
				// Now calculate the integrals.
				for (ir = 0; ir < r.length(); ir++) {
					// Break if we are wasting time...
					if (m >= m0[ir] && m >= m1[ir])
						continue;
					double t1, t2, t3, t4, t5, t6;
					t2 = m*j1(m*a[ia])*w;
					t1 = t2*j0(m*r[ir]);
					t2 *= j1(m*r[ir]);
					for (iz = 0; iz < z.length(); iz++) {
						axialdata & s = ax[ir*z.length()+iz];
						double & tz = z[iz];
						il = zl[iz];
						t3 = m*(ABCD[il][2] + ABCD[il][3]*tz)*exp(-m*tz);
						t4 = ABCD[il][3]*exp(-m*tz);
						t5 = m*(ABCD[il][0] + ABCD[il][1]*tz)*
										 (m*tz<MAX_EXP?exp(m*tz):DBL_MAX);
						t6 = ABCD[il][1]*(m*tz<MAX_EXP?exp(m*tz):DBL_MAX);
						if (m < m0[ir]) {
							if (result == all) {
								s.vse += t1*m*((1-2*v[il])*(t4+t6)+(t3-t5));
								s.rse += t1*m*((1+2*v[il])*(t4+t6)-(t3-t5));
								s.tse += t1*m*(2*v[il])*(t4+t6);
							}
							s.vdp += t1*((2-4*v[il])*(t6-t4)-(t3+t5));
						}
						if (m < m1[ir] && result == all) {
							s.sse += t2*m*((2*v[il])*(t6-t4)+(t3+t5));
							s.rdp += t2*((t4+t6)-(t3-t5));
						}
					}
				}
			}
		}

		// Finalise the calculations, now that we have done the integration.
		for (ir = 0; ir < r.length(); ir++) {
			for (iz = 0; iz < z.length(); iz++) {
				axialdata & s = ax[ir*z.length()+iz];
				il = zl[iz];
				if (result == all) {
					if (r[ir] > 0.0) {
						s.rse -= s.rdp/r[ir];
						s.tse += s.rdp/r[ir];
					}
					s.rse *= a[ia];
					s.tse *= a[ia];
					s.vse *= a[ia];
					s.sse *= a[ia];
					s.rdp *= a[ia]*(v[il]+1)/E[il];
				}
				s.vdp *= a[ia]*(v[il]+1)/E[il];
			}
		}

		// After doing everything in radial coords, translate to cartesian.
		for (ixy = 0; ixy < data.length(); ixy++) {
			pavedata & d = data[ixy];
			for (iz = 0; z[iz] != d.z; iz++)
				;
			for (ild = 0; ild < load.length(); ++ild) {
				if (load[ild].radius() != a[ia])
					continue;
				for (ir = 0; r[ir] != load[ild].distance(d); ir++)
					;
				axialdata & s = ax[ir*z.length()+iz];
				double p = load[ild].pressure();
				if (r[ir] == 0.0 && result == all) {
					d.data[0][0] += p*(s.rse+s.tse)/2;
					d.data[0][1] += p*(s.rse+s.tse)/2;
					d.data[0][2] += p*s.vse;
					d.data[1][1] += p*s.sse;
				} else if (result == all) {
					double cost, sint;
					sint = (d.y-load[ild].y)/r[ir];
					cost = (d.x-load[ild].x)/r[ir];
					d.data[0][0] += p*(cost*cost*s.rse+sint*sint*s.tse);
					d.data[0][1] += p*(cost*cost*s.tse+sint*sint*s.rse);
					d.data[0][2] += p*s.vse;
					d.data[1][0] += p*cost*sint*(s.rse-s.tse);
					d.data[1][1] += p*cost*s.sse;
					d.data[1][2] += p*sint*s.sse;
					d.data[4][0] += p*cost*s.rdp;
					d.data[4][1] += p*sint*s.rdp;
				}
				if (gl != -1)
					d.deflgrad[gl] += p*s.vdp;
				else
					d.data[4][2] += p*s.vdp;
			}
		}
		if (result == dispgrad) {
			if (++gl != nl) {
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
			gl = -1;
		}
	}
	// After everything, loop through the answers, and calculate the
	// derived results (principal stresses and strains).
	for (ixy = 0; ixy < data.length(); ixy++) {
		if (result == all) {
			for (iz = 0; z[iz] != data[ixy].z; iz++)
				;
			il = zl[iz];
			data[ixy].principle(v[il],E[il]);
		}
		if (result == dispgrad) {
			if (Q == 0) {
				for (il = 0; il < nl; il++)
					data[ixy].deflgrad[il] = (data[ixy].deflgrad[il]
						-data[ixy].data[4][2])/GRADSTEP;
			} else {
				for (il = 0; il < nl; il++)
					h[il] = (data[ixy].deflgrad[il]-data[ixy].data[4][2])
								/GRADSTEP;
				for (gl = 0; gl < nl; gl++) {
					data[ixy].deflgrad[gl] = 0.0;
					for (il = 0; il < nl; il++)
						data[ixy].deflgrad[gl] += h[il]*Q[gl*nl+il];
				}
			}
		}
    }
	rv = true;

abort:
	if (rv == false)
		event_msg(EVENT_ERROR,"Out of memory in LEsystem::calculate()!");
	delete [] R;
	delete [] ABCD;
	delete [] h;
	delete [] f;
	delete [] v;
	delete [] E;
	delete _z;
	delete _a;
	delete [] zl;
	delete _r;
	delete _bm;
	delete [] ax;
	delete [] m0;
	delete [] m1;
	return rv;
};

/*
 * This is Odemark's method of equivalent thicknesses, with Bosussinesq's
 * solution for a point load on a linear elastic half-space.  This can be
 * used for very fast approximations.
 */
bool
LEsystem::odemark()
{
	int ixy, ild, il;
	LElayer * pl;
	bool rv = true;
	double de = 0.0;

	if (!check())
		return false;
	int nl = layers();
	double * h = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	if (h == 0 || v == 0 || E == 0) {
		rv = false;
		goto abort;
	}
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++)
		h[il] = pl->bottom(), v[il] = pl->poissons(), E[il] = pl->emod();

	for (ixy = 0; ixy < data.length(); ixy++)
		memset(data[ixy].data,0,sizeof(data[ixy].data));
	for (pl = first; pl != 0; pl = pl->next) {
		double v = pl->poissons();
		double E = pl->emod();
		double he = pl->thickness();
		if (pl->next != 0) {
			double v1 = pl->next->poissons();
			double E1 = pl->next->emod();
			he *= pow(E/E1*(1+v1*v1)/(1+v*v),1.0/3.0);
		}
		for (ild = 0; ild < load.length(); ild++) {
			double p = load[ild].force();
			double a = load[ild].radius();
			for (ixy = 0; ixy < data.length(); ixy++) {
				pavedata & d = data[ixy];
				double z, r, R;
				if (pl->bottom() != 0.0 && d.z >= pl->bottom())
					continue;
				r = load[ild].distance(d);
				if (pl->top() > d.z) {
					// We're below, but we still need the vertical
					// deflection contributions from this layer...
					z = de;
					R = pythag(r,z);
					if (R > 0.0)
						d.data[4][2] += p*((1+v)/(M_2PI*R*E))*
													(2*(1-v) + z*z/R/R);
					else
						d.data[4][2] += p*2*(1-v*v)/E/M_PI/a;
					if (pl->bottom() != 0.0) {
						z = de+pl->thickness();
						R = pythag(r,z);
						if (R > 0.0)
							d.data[4][2] -= p*((1+v)/(M_2PI*R*E))*
													(2*(1-v) + z*z/R/R);
						else
							d.data[4][2] -= p*2*(1-v*v)/E/M_PI/a;
					}
					continue;
				}
				z = de+(d.z-pl->top())*(pl->bottom() != 0.0 ?
											he/pl->thickness() : 1.0);
				R = pythag(r,z);
				axialdata s;
				if (R > 0.0) {
					s.vse = -3*z*z*z/pow(R,5)/M_2PI;
					s.rse = -(3*z*r*r/pow(R,5)-(1-2*v)/R/(R+z))/M_2PI;
					s.tse = -(1-2*v)*(-z/pow(R,3) + 1.0/R/(R+z))/M_2PI;
					s.sse = -3*z*z*r/pow(R,5)/M_2PI;
					s.rdp = ((1+v)/(M_2PI*R*E))*(z*r/R/R-(1-2*v)*r/(R+z));
					s.vdp = ((1+v)/(M_2PI*R*E))*(2*(1-v) + z*z/R/R);
				} else {
					s.vse = -1.0/(M_PI*a*a);
					s.rse = s.tse = -(1+2*v)/2.0;
					s.sse = s.rdp = 0.0;
					s.vdp = 2*(1-v*v)/(E*M_PI*a);
				}
				if (pl->bottom() != 0.0) {
					z = de+he;
					R = pythag(r,z);
					if (R > 0.0)
						s.vdp -= ((1+v)/(M_2PI*R*E))*(2*(1-v) + z*z/R/R);
					else
						s.vdp -= 2*(1-v*v)/(E*M_PI*a);
				}

				// After doing everything in radial coords, translate
				// to cartesian.
				if (r == 0.0) {
					d.data[0][0] += p*(s.rse+s.tse)/2;
					d.data[0][1] += p*(s.rse+s.tse)/2;
					d.data[1][1] += p*s.sse;
					d.data[4][0] += p*s.rdp;
				} else {
					double sint = (d.y-load[ild].y)/r;
					double cost = (d.x-load[ild].x)/r;
					d.data[0][0] += p*(cost*cost*s.rse+sint*sint*s.tse);
					d.data[0][1] += p*(cost*cost*s.tse+sint*sint*s.rse);
					d.data[1][0] += p*cost*sint*(s.rse-s.tse);
					d.data[1][1] += p*cost*s.sse;
					d.data[1][2] += p*sint*s.sse;
					d.data[4][0] += p*cost*s.rdp;
					d.data[4][1] += p*sint*s.rdp;
				}
				d.data[0][2] += p*s.vse;
				d.data[4][2] += p*s.vdp;
			}
		}
		de += he;
	}
	// After everything, loop through the answers, and calculate the
	// derived results (principal stresses and strains).
	for (ixy = 0; ixy < data.length(); ixy++) {
		for (pl = first; pl != 0; pl = pl->next) {
			if (pl->top() <= data[ixy].z
			  && (pl->bottom() == 0.0 || data[ixy].z < pl->bottom()))
				break;
		}
		data[ixy].principle(pl->poissons(),pl->emod());
    }
abort:
	if (rv == false)
		event_msg(EVENT_ERROR,"Out of memory in LEsystem::odemark()!");
	delete [] h;
	delete [] v;
	delete [] E;
	return rv;
};

/*
 * The deflection equation is derived from the Boissenq solution,
 * integrated and simplified as much as possible.  However, it cannot
 * completely integrated because there is a problem at points where
 * the load is directly over the evaluation point, so...
 *
 * An implementation of adaptive, recursive Newton-Cotes integration.
 * Based on the Matlab implementation, but covered in a lot of books...
 */
#define LEVMAX	12
static double
quad8_vdp(double r, double z, double s, double v, double a = 0.0,
          double b = M_PI, double Q = 10.0)
{
	// The magic Newton-Cotes weights
	const double w[9] = {3956, 23552, -3712, 41984, -18160, 41984,
						-3712, 23552, 3956};
	const double dw = 14175;
	static int level = -1;
	static double tol = 1e-6;
	register double h, t, Q1 = 0.0, Q2 = 0.0;
	register int i;

	if (s == 0.0)
		return 0.0;
	if (r == 0.0)
		return (2*(v*v-1)*(z-pythag(z,s))
						+(v+1)*(z-z*z/pythag(z,s)));
	level++;
	h = (b-a)/16.0;
	for (i = 0; i < 9; i++) {
		t = cos(a+i*h)*r;
		if (z == 0.0)
			Q1 += h*w[i]/dw*(t == r ? 2*(v*v-1)*(r-sqrt(s*s-2*s*t+r*r)) :
				2*(v*v-1)*(t*(log(r-t)-log(s-t+sqrt(s*s-2*s*t+r*r)))
					-sqrt(s*s-2*s*t+r*r)+r)
			);
		else
			Q1 += h*w[i]/dw*(t*t == r*r+z*z ? 2*(v*v-1)*s :
				2*(v*v-1)*
					(t*(log(pythag(z,r)-t)
					   -log(s-t+sqrt(z*z-2*s*t+s*s+r*r)))
					 -sqrt(z*z-2*s*t+s*s+r*r)+pythag(z,r))
				 +(v+1)*z*z*((z*z+r*r-t*s)/sqrt(z*z-2*s*t+s*s+r*r)
					 -pythag(z,r))/(t*t-z*z-r*r)
			);
		t = cos(a+(i+8)*h)*r;
		if (z == 0.0)
			Q2 += h*w[i]/dw*(t == r ? 2*(v*v-1)*(r-sqrt(s*s-2*s*t+r*r)) :
				2*(v*v-1)*(t*(log(r-t)-log(s-t+sqrt(s*s-2*s*t+r*r)))
					-sqrt(s*s-2*s*t+r*r)+r)
			);
		else
			Q2 += h*w[i]/dw*(t*t == r*r+z*z ? 2*(v*v-1)*s :
				2*(v*v-1)*
					(t*(log(pythag(z,r)-t)
					   -log(s-t+sqrt(z*z-2*s*t+s*s+r*r)))
					 -sqrt(z*z-2*s*t+s*s+r*r)+pythag(z,r))
				 +(v+1)*z*z*((z*z+r*r-t*s)/sqrt(z*z-2*s*t+s*s+r*r)
					 -pythag(z,r))/(t*t-z*z-r*r)
			);
	};
	// This is the adaptive recursive bit.  We only recurse if we
	// can improve...
	if (fabs(Q1+Q2-Q) > tol*fabs(Q1+Q2) && level <= LEVMAX) {
		tol = tol/2;
		Q1 = quad8_vdp(r,z,s,v,a,(a+b)/2,Q1);
		Q2 = quad8_vdp(r,z,s,v,(a+b)/2,b,Q2);
		tol = tol*2;
	}
	level--;
	if (level == -1)
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
LEsystem::fastnum()
{
	int ixy, ild, il;
	LElayer * pl;
	bool rv = true;

	if (!check())
		return false;
	int nl = layers();
	double * h = new double[nl];
	double * v = new double[nl];
	double * E = new double[nl];
	if (h == 0 || v == 0 || E == 0) {
		rv = false;
		goto abort;
	}
	for (pl = first, il = 0; pl != 0; pl = pl->next, il++)
		h[il] = pl->bottom(), v[il] = pl->poissons(), E[il] = pl->emod();

	for (ixy = 0; ixy < data.length(); ixy++) {
		pavedata & d = data[ixy];
		d.data[4][2] = 0.0;
		for (ild = 0; ild < load.length(); ild++) {
			double a = load[ild].radius();
			double z, r = load[ild].distance(d);
			double he = 0.0, hc = 1.0, vdp = 0.0;
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
				vdp += quad8_vdp(r,z,a,v[il])/E[il];
				if (h[il] != 0.0) {
					z = hc*he+h[il]-(il > 0 ? h[il-1] : 0.0);
					vdp -= quad8_vdp(r,z,a,v[il])/E[il];
				}
			}
			d.data[4][2] += load[ild].pressure()*vdp;
		}
	}
abort:
	if (rv == false)
		event_msg(EVENT_ERROR,"Out of memory in LEsystem::fastnum()!");
	delete [] h;
	delete [] v;
	delete [] E;
	return rv;
};

/*
 * The overall backcalculation routine.
 */
bool
LEbackcalc::backcalc()
{
	int i, j, ns = 0, nl = layers();
	double derr = DBL_MAX, oerr = 0.0;
	double astep, tstep = 0.0, step = 1.0;
	bool seeded = true, badstep = false;
	calctype speed = (precision >= 1e-4 ? fast : slow);
	ksset<point3d,pavedata> orig(data);

	double * P = new double[nl];
	double * T = new double[nl];
	double * O = new double[nl];
	double * Q = new double[nl*nl];
	if (P == 0 || T == 0 || O == 0 || Q == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEbackcalc::backcalc()!");
		goto abort;
	}

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
		badstep = seed(nl,P);
		memcpy(T,P,sizeof(double)*nl);
		memcpy(O,P,sizeof(double)*nl);
	}
	// Now settle into our true minimisation algorithm.
	// First we start by using the deflection gradient in
	// emod space, since it converges faster...
	for (j = 0; j < maxsteps*nl; j++) {
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
	for (ns = 0; j < maxsteps*nl; j++) {
		for (i = 0; i < nl; i++)
			layer(i).emod(pow(10,P[i]));
		calculate((speed == fast ? LEsystem::fastgrad : LEsystem::dispgrad));
		oerr = derr, derr = (precision > 0.0
			? ROUND(bowlerror()/precision)*precision : bowlerror());
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
	for (j = 0; j < defl.length(); j++) {
		defldata & d = defl[j];
		d.calculated = result(d).result(pavedata::deflct, pavedata::zz);
	}
	// Restore original evaluation points...
	data = orig;
abort:
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
LEbackcalc::seed(int nl, double * P)
{
	double E = 0.0, v = layer(nl-1).poissons();
	int i, j;
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
 * The algorithm can also be called for a straight Gauss-Newton method
 * or for a differential deflection method, which aids convergence for
 * some bowls.
 *
 * See the ICAP 2006 paper for details.
 */
double
LEbackcalc::deflgrad(int nl, double * P, double * Q, calctype cl) {
	double step = 0.0, dgg = 0.0, gg = 0.0, dd = 0.0;
	int i, j, k, dl = defl.length();

	// Initial setup.
	double * PG = new double[dl*nl];
	double * MD = new double[dl];
	double * CD = new double[dl];
	double * DG = new double[dl];
	double * D = new double[nl];
	double * G = new double[nl];
	double * W = new double[nl];
	if (PG == 0 || MD == 0 || CD == 0 || DG == 0
	  || D == 0 || G == 0 || W == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEbackcalc::deflgrad()!");
		goto abort;
	}
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
	for (k = 0; k <= nl; k++) {
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
abort:
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
LEbackcalc::gaussnewton(int nl, double * P, calctype cl)
{
	int i, j, k, dl = defl.length();
	double step = 0.0;

	// Initial setup.
	double * H = new double[dl*nl];
	double * S = new double[nl*nl];
	double * W = new double[nl];
	double * Y = new double[nl];
	if (H == 0 || S == 0 || W == 0 || Y == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEbackcalc::kalman()!");
		goto abort;
	}
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
	equ_eig(nl,S,Y,W);
	for (i = 0, step = 0.0; i < nl; i++)
		step += W[i]*W[i], P[i] += W[i];
	step = sqrt(step);
abort:
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
double LEbackcalc::kalman(int nl, double * P) {
	int i, j, k, dl = defl.length();
	double step = 0.0;

	// Initial setup.
	double * H = new double[dl*nl];
	double * S = new double[dl*dl];
	double * K = new double[nl*dl];
	double * Y = new double[dl];
	double * W = new double[nl];
	if (H == 0 || S == 0 || K == 0 || Y == 0 || W == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEbackcalc::kalman()!");
		goto abort;
	}
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
	inv_eig(dl,S);
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
	if (step > 1.0) {
		step = step*brent(nl,P,W);
	} else {
		for (i = 0; i < nl; i++)
			P[i] = MIN(MAX(1,P[i]+W[i]),9);
	}
abort:
	if (W != 0)
		delete [] W;
	if (Y != 0)
		delete [] Y;
	if (K != 0)
		delete [] K;
	if (S != 0)
		delete [] S;
	if (H != 0)
		delete [] H;
	return step;
}

/*
 * Bowl error function.  This works out mean deviation of the difference
 * between the caclulated and measured deflections.  Call with P=0 to
 * avoid it rerunning the calcs.  Call with s and D to get it to do
 * prediction at distance s along the vector D.
 */
double LEbackcalc::bowlerror(int nl, double * P, double s, double * D) {
	int i;
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
double LEbackcalc::brent(int nl, double * P, double * D) {
	const double GC = (3.0-sqrt(5.0))/2.0;
	const double GR = 1/(1-GC);
	const double TOL = 1e-6;

	double A, B, C, X, Y, Z, a, b, c, x, y, z;
	double r, q;
	int i;
	
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
double LEbackcalc::conjgrad(int nl, double * P) {
	double step = 1.0, dgg = 0.0, gg = 0.0;
	double derr = DBL_MAX, oerr = 0.0;
	int i, j, k, dl = defl.length();

	// Initial setup.
	double * D = new double[nl];
	double * G = new double[nl];
	double * W = new double[nl];
	if (D == 0 || G == 0 || W == 0) {
		event_msg(EVENT_ERROR,"Out of memory in LEbackcalc::conjgrad()!");
		goto abort;
	}
	memcpy(W,P,sizeof(double)*nl);
	memset(D,0,sizeof(double)*nl);
	for (k = 0; k <= nl*maxsteps; k++) {
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
abort:
	if (W != 0)
		delete [] W;
	if (G != 0)
		delete [] G;
	if (D != 0)
		delete [] D;
	return sqrt(gg);
}
