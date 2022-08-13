/**************************************************************************

	STATISTICS.CPP - Implementation for basic statistical functions.

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

	See STATISTICS.H

	History:
		2002/01/10 - Created by Jeremy Lea <reg@openpave.org>
		2016/02/02 - Seperated out of RELIABILITY.H

**************************************************************************/

#include <cstdlib>
#include <cmath>
#include "mathplus.h"
#include "statistics.h"
#include "randomvar.h"

namespace OP {

/*
 * A normally distributed random number generator.  We avoid
 * the uniform rv's being 0.0 since this will result in infinte
 * values, and double count the 0 == 2pi.
 */
double stdnormal_rnd() {
	static int i = 1;
	static double u[2] = {0.0, 0.0};
	double r[2];

	if (i == 1) {
		r[0] = sqrt(-2*log((double(std::rand())+1)/(double(RAND_MAX)+1)));
		r[1] = M_2PI*(double(std::rand())+1)/(double(RAND_MAX)+1);
		u[0] = r[0]*sin(r[1]);
		u[1] = r[0]*cos(r[1]);
		i = 0;
	} else {
		i = 1;
	}

	return u[i];
}

/*
 * The standard normal PDF, for one random variable.
 */
double
stdnormal_pdf(double u)
{
	return exp(-u*u/2)/M_SQRT2PI;
}

/*
 * An implementation of adaptive, recursive Newton-Cotes integration.
 * Based on the Matlab implementation, but covered in a lot of books...
 *
 * This only does integration over the standard normal PDF.  It's just
 * here to check the error function approximations.
 */
#define LEVMAX	10
double
quad8_stdnormal_pdf(double a, double b, double Q = 1.0)
{
	// The magic Newton-Cotes weights
	const int w[9] = {3956, 23552, -3712, 41984, -18160, 41984, -3712, 23552, 3956};
	const int dw = 14175;
	static int level = -1;
	static double tol = 1e-30;
	double h, Q1 = 0.0, Q2 = 0.0;
	int i;

	level++;
	h = (b-a)/16.0;
	for (i = 0; i < 9; i++) {
		Q1 += h*w[i]*stdnormal_pdf(a+i*h)/dw;
		Q2 += h*w[i]*stdnormal_pdf(a+(i+8)*h)/dw;
	}
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
double
stdnormal_cdf(double u)
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
	double y, z;

	if (std::isnan(u))
		return NAN;
	if (!std::isfinite(u))
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
}

/*
 * The inverse standard normal distribution.
 *
 *   Author:      Peter J. Acklam <jacklam@math.uio.no>
 *   URL:         http://www.math.uio.no/~jacklam
 *
 * This function is based on the Matlab code from the address above,
 * translated to C, and adapted for our purposes.
 */
double
stdnormal_inv(double p)
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

	double q, t, u;

	if (std::isnan(p) || p > 1.0 || p < 0.0)
		return NAN;
	if (p == 0.0)
		return -INFINITY;
	if (p == 1.0)
		return INFINITY;
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
}

// Creates a new random variable of the type specified.
random *
random::make_rv(distribution d, float m, float s, house * h) {
	switch (d) {
	case distribution::normal:
		return new rv_normal(m,s,h);
	case distribution::lognormal:
		return new rv_lognormal(m,s,h);
	default:
		throw std::runtime_error("Trying to create an unknown rv type!");
	}
}

random *
random::make_rv(const random & r, house * h) {
	switch (r.type()) {
	case distribution::normal:
		return new rv_normal(r,h);
	case distribution::lognormal:
		return new rv_lognormal(r,h);
	default:
		throw std::runtime_error("Trying to create an unknown rv type!");
	}
}

} // namespace OP
