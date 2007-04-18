/**************************************************************************

	RELIABILIY.CPP - Implementation for the reliability class.

	$OpenPave$

	See RELIABILITY.H

	History:
		2002/01/10 - Created by Jeremy Lea <jlea@csir.co.za>

**************************************************************************/

#include <stdlib.h>
#include "..\include\mathplus.h"
#include "reliability.h"

/*
 * A normally distributed random number generator.  We avoid
 * the uniform rv's being 0.0 since this will result in infinte
 * values, and double count the 0 == 2pi.   
 */
double stdnormal_rnd() {
	static int i = 1;
	static double u[2] = {0.0, 0.0};
	register double r[2];

	if (i == 1) {
		r[0] = sqrt(-2*log((double)(rand()+1)/(double)(RAND_MAX+1)));
		r[1] = M_2PI*(double)(rand()+1)/(double)(RAND_MAX+1);
		u[0] = r[0]*sin(r[1]);
		u[1] = r[0]*cos(r[1]);
		i = 0;
	} else {
		i = 1;
	}

	return u[i];
};

/*
 * The standard normal PDF, for one random variable.
 */
double stdnormal_pdf(const double u)
{
	return exp(-u*u/2)/M_SQRT2PI;
};

/*
 * An implementation of adaptive, recursive Newton-Cotes integration.
 * Based on the Matlab implementation, but covered in a lot of books...
 *
 * This only does integration over the standard normal PDF.  It's just
 * here to check the error function approximations.
 */
#define LEVMAX	10
double quad8_stdnormal_pdf(const double a, const double b, const double Q = 1.0)
{
	/* The magic Newton-Cotes weights */
	const int w[9] = {3956, 23552, -3712, 41984, -18160, 41984, -3712, 23552, 3956};
	const int dw = 14175;
	static int level = -1;
	static double tol = 1e-30;
	register double h, Q1 = 0.0, Q2 = 0.0;
	register int i;

	level++;
	h = (b-a)/16.0;
	for (i = 0; i < 9; i++) {
		Q1 += h*w[i]*stdnormal_pdf(a+i*h)/dw;
		Q2 += h*w[i]*stdnormal_pdf(a+(i+8)*h)/dw;
	};
	/* This is the adaptive recursive bit.  We only recurse if we can improve... */
	if (fabs(Q1+Q2-Q) > tol*fabs(Q1+Q2) && level <= LEVMAX) {
		tol = tol/2;
		Q1 = quad8_stdnormal_pdf(a,(a+b)/2,Q1);
		Q2 = quad8_stdnormal_pdf((a+b)/2,b,Q2);
		tol = tol*2;
	}
	level--;
	return Q1 + Q2;
}

/*
 * The standard normal CDF, for one random variable.
 *
 *   Author:		W. J. Cody
 *   URL:			http://www.netlib.org/specfun/erf
 *
 * This is the erfc() routine only, adapted by the
 * transform stdnormal_cdf(u)=(erfc(-u/sqrt(2))/2;
 */
double stdnormal_cdf(const double u)
{
	const double a[5] = {
		1.161110663653770e-002,3.951404679838207e-001,2.846603853776254e+001,
		1.887426188426510e+002,3.209377589138469e+003
	};
	const double b[5] = {
		1.767766952966369e-001,8.344316438579620e+000,1.725514762600375e+002,
		1.813893686502485e+003,8.044716608901563e+003
	};
	const double c[9] = {
		1.07655767737201917e-8,3.98941512088134664e-1,8.88314979438837682e00,
		9.35066561321778664e01,5.97270276394800248e02,2.49453758529037214e03,
		6.84819045053628270e03,1.16026514376473497e04,9.84271483838397762e03
	};
	const double d[9] = {
		1.00000000000000000e00,2.22666880443281165e01,2.35387901782624994e02,
		1.51937759940755473e03,6.48555829826676108e03,1.86155716408850967e04,
		3.49009527211459790e04,3.89120032860932697e04,1.96854296768599925e04
	};
	const double p[6] = {
		2.30734417649401738e-2,2.15898534057957003e-1,1.27401161160247384e-1,
		2.22352778706498104e-2,1.42161919322789337e-3,2.91128749511687932e-5
	};
	const double q[6] = {
		1.00000000000000000e00,1.28426009614491110e00,4.68238212480865112e-1,
		6.59881378689285564e-2,3.78239633202758245e-3,7.29751555083966178e-5
	};
	register double y, z;

	if (_isnan(u))
#if _MSC_VER > 1200
		return _Nan._Double;
#else
		return _Nan._D;
#endif
	if (!_finite(u))
		return (u < 0 ? 0.0 : 1.0);
	y = fabs(u);
	if (y <= 0.662912607362388) {
		// evaluate for |u| <= sqrt(2)*0.46875
		y = u*u;
		y = u*((((a[0]*y+a[1])*y+a[2])*y+a[3])*y+a[4])
		     /((((b[0]*y+b[1])*y+b[2])*y+b[3])*y+b[4]);
		return 0.5+y;
	}
	if (y <= 5.65685424949238) {
		// evaluate for sqrt(2)*0.46875 <= |u| <= sqrt(2)*4.0
		y = ((((((((c[0]*y+c[1])*y+c[2])*y+c[3])*y+c[4])*y+c[5])*y+c[6])*y+c[7])*y+c[8])
		   /((((((((d[0]*y+d[1])*y+d[2])*y+d[3])*y+d[4])*y+d[5])*y+d[6])*y+d[7])*y+d[8]);
	} else {
		// evaluate for |u| > sqrt(2)*4.0
		z = 1.0/y;
		y = z*z;
		y *= (((((p[0]*y+p[1])*y+p[2])*y+p[3])*y+p[4])*y+p[5])
		    /(((((q[0]*y+q[1])*y+q[2])*y+q[3])*y+q[4])*y+q[5]);
		y = z*(M_1_SQRT2PI-y);
	}
	z = floor(u*16)/16;
	y *= exp(-0.5*z*z)*exp(-0.5*(u-z)*(u+z));
	return (u < 0.0 ? y : 1-y);
};

/*
 * The inverse standard normal distribution.
 *
 *   Author:      Peter J. Acklam <jacklam@math.uio.no>
 *   URL:         http://www.math.uio.no/~jacklam
 *
 * This function is based on the Matlab code from the address above,
 * translated to C, and adapted for our purposes.
 */
double stdnormal_inv(double p)
{
	const double a[6] = {
		-3.969683028665376e+01,  2.209460984245205e+02,
		-2.759285104469687e+02,  1.383577518672690e+02,
		-3.066479806614716e+01,  2.506628277459239e+00
	};
	const double b[5] = {
		-5.447609879822406e+01,  1.615858368580409e+02,
		-1.556989798598866e+02,  6.680131188771972e+01,
		-1.328068155288572e+01
	};
	const double c[6] = {
		-7.784894002430293e-03, -3.223964580411365e-01,
		-2.400758277161838e+00, -2.549732539343734e+00,
		 4.374664141464968e+00,  2.938163982698783e+00
	};
	const double d[4] = {
		 7.784695709041462e-03,  3.224671290700398e-01,
		 2.445134137142996e+00,  3.754408661907416e+00
	};

	register double q, t, u;

	if (_isnan(p) || p > 1.0 || p < 0.0)
#if _MSC_VER > 1200
		return _Nan._Double;
#else
		return _Nan._D;
#endif
	if (p == 0.0)
#if _MSC_VER > 1200
		return -_Inf._Double;
#else
		return -_Inf._D;
#endif
	if (p == 1.0)
#if _MSC_VER > 1200
		return _Inf._Double;
#else
		return _Inf._D;
#endif
	q = MIN(p,1-p);
	if (q > 0.02425) {
		// Rational approximation for central region.
		u = q-0.5;
		t = u*u;
		u = u*(((((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4])*t+a[5])
			 /(((((b[0]*t+b[1])*t+b[2])*t+b[3])*t+b[4])*t+1);
	} else {
		// Rational approximation for tail region.
		t = sqrt(-2*log(q));
		u = (((((c[0]*t+c[1])*t+c[2])*t+c[3])*t+c[4])*t+c[5])
			/((((d[0]*t+d[1])*t+d[2])*t+d[3])*t+1);
	}
	// The relative error of the approximation has absolute value less
	// than 1.15e-9.  One iteration of Halley's rational method (third
	// order) gives full machine precision...
	t = stdnormal_cdf(u)-q;				// error
	t *= M_SQRT2PI*exp(u*u/2);			// f(u)/df(u)
	u -= t/(1+u*t/2);					// Halley's method
	return (p > 0.5 ? -u : u);
};

/*
 *	Simple constructor.  We always have an owner, and optionally a previous list member.
 */
randomvar::randomvar(reliability * o, randomvar * p) {
	owner = o;
	prev = p, next = 0;
	if (prev != 0) {
		if (prev->next != 0) {
			next = prev->next;
			next->prev = this;
		}
		prev->next = this;
	}

	stdnormal = 0.0;
	ustdnormal = 0.0;
	for (int i = 0; i < RV_DIST_PARAM; i++)
		d[i] = 0.0;
};

/*
 *	A simple destructor which does a bit of list management...
 */
randomvar::~randomvar() {
	if (prev != 0) {
		prev->next = next;
	}
	if (next != 0) {
		next->prev = prev;
	}
};

/*
 * An implementation of a normally distributed random variable.
 */
struct rv_normal : public randomvar
{
public:
	virtual distribution type() {
		return randomvar::NORMAL;
	};
	virtual double x() {
		return d[0] + d[1]*z();
	};
	virtual double mean() {
		return d[0];
	};
	virtual double stddev() {
		return d[1];
	};

protected:
	friend reliability;

	rv_normal(reliability * o, randomvar * p, double m, double s)
		: randomvar(o,p) {
		param(0, m);
		param(1, s);
	}
	virtual ~rv_normal() {
	};
};

/*
 * An implementation of a natural log normally distributed random variable.
 */
struct rv_lognormal : public randomvar
{
public:
	virtual distribution type() {
		return randomvar::LOGNORMAL;
	};
	virtual double x() {
		return exp(d[0] + d[1]*z());
	};
	virtual double mean() {
		return exp(d[0] + d[1]*d[1]/2);
	};
	virtual double stddev() {
		return mean()*sqrt(exp(d[1]*d[1])-1);
	};

protected:
	friend reliability;

	rv_lognormal(reliability * o, randomvar * p, double m, double s)
		: randomvar(o,p) {
		param(0, sqrt(log(1+s*s/m/m)));
		param(1, log(m) - d[1]*d[1]/2);
	}
	virtual ~rv_lognormal() {
	};
};

/*
 *	Simple constructor.  We always have an owner, and optionally a previous list member.
 */
gfunction::gfunction(reliability * o, gfunction * p) {
	owner = o;
	prev = p, next = 0;
	if (prev != 0) {
		if (prev->next != 0) {
			next = prev->next;
			next->prev = this;
		}
		prev->next = this;
	}
};

/*
 *	A simple destructor which does a bit of list management...
 */
gfunction::~gfunction() {
	if (prev != 0) {
		prev->next = next;
	}
	if (next != 0) {
		next->prev = prev;
	}
};

/*
 *	A simple construcutor for a reliability problem.
 */
reliability::reliability() {
	rv_head = 0;
	gf_head = 0;
};

/*
 *	Simple destuctor.  delete's all of our lists.
 */
reliability::~reliability() {
	if (rv_head != 0) {
		while (rv_head->next != 0) {
			delete rv_head->next;
		}
		delete rv_head;
	}
	rv_head = 0;
	if (gf_head != 0) {
		while (gf_head->next != 0) {
			delete gf_head->next;
		}
		delete gf_head;
	}
	gf_head = 0;
};

/*
 * Creates a new random variable of the type specified.
 */
randomvar *
reliability::NewRV(randomvar::distribution t, double m, double s) {
	randomvar * rv;

	rv = rv_head;
	while (rv != 0)
		rv = rv->next;
	switch (t) {
	case (randomvar::NORMAL):
		rv = new rv_normal(this, rv, m, s);
		break;
	case (randomvar::LOGNORMAL):
		rv = new rv_lognormal(this, rv, m, s);
		break;
	}
	if (!rv) {
		// XXX: Error.
		return 0;
	}
	if (rv_head == 0)
		rv_head = rv;

	return rv;
};

