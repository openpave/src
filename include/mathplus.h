/**************************************************************************

	MATHPLUS.H - Some additional mathematical constants and macros.

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

	Portions Copyright (C) 2006 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header includes all of the normal fun maths headers, and
		defines a few other constants not in the standard headers.  It
		also defines a few simple inline helper functions, and smooths
		out the differences between compilers...

	History:
		Ages ago   - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __MATHPLUS_H
#define __MATHPLUS_H

#include <float.h>
#include <limits.h>
#include <math.h>
#include <ymath.h>
#include <stdlib.h>

#define	M_E			2.71828182845904523536		// e
#define M_1_E       .367879441171442321596		// 1/e
#define	M_LOG2E		1.44269504088896340736		// log 2e
#define	M_LOG10E	.434294481903251827651		// log 10e
#define	M_LN2		.693147180559945309417		// log e2
#define	M_LN10		2.30258509299404568402		// log e10
#define	M_PI		3.14159265358979323846		// pi
#define M_2PI       6.28318530717958647693		// 2*pi
#define	M_PI_2		1.57079632679489661923		// pi/2
#define	M_PI_4		.785398163397448309616		// pi/4
#define	M_1_PI		.318309886183790671538		// 1/pi
#define	M_2_PI		.636619772367581343076		// 2/pi
#define M_PISQUARE  9.86960440108935861883		// pi*pi
#define M_SQRTPI    1.77245385090551602730		// sqrt(pi)
#define M_SQRT2PI   2.50662827463100050242		// sqrt(2*pi)
#define M_1_SQRTPI  .564189583547756286948		// 1/sqrt(pi)
#define	M_2_SQRTPI	1.12837916709551257390		// 2/sqrt(pi)
#define M_180_PI    57.2957795130823208768		// 180/pi
#define M_PI_180    1.74532925199432957692e-2 	// pi/180
#define M_LNPI      1.14472988584940017414		// log pi
#define	M_SQRT2		1.41421356237309504880		// sqrt(2)
#define	M_SQRT1_2	.707106781186547524401		// 1/sqrt(2)
#define M_EULER     .577215664901532860607		// euler's constant

#define	MAX_EXP		log(DBL_MAX)				// Max/Min value for exp()

#define ABS(x)		((x)>=0.0?(x):-(x))
#define MAX(x,y)	((x)>=(y)?(x):(y))
#define MIN(x,y)	((x)<=(y)?(x):(y))
#define SGN(x)		((x)>0.0?1:((x)<0.0?-1:0))
#define FLOOR(x)	((int)((x)+(1-(int)(x)))-(1-(int)(x)))
#define CEIL(x)		((2+(int)(x))-(int)((2+(int)(x))-(x)))
#define ROUND(x)	((x)>=0?(int)((x)+0.5):-((int)(0.5-(x))))

#define RAND(a,b)	((a)+((b)-(a))*(double)(rand()+1)/(double)(RAND_MAX+2))

// Make some protable names for the bessel functions.
#define	j0(x)	_j0(x)
#define j1(x)	_j1(x)

/*
 * Simple template functions for swapping things. One day
 * someone will make these part of the language...
 */
template <class T>
void inline swap(T &a, T &b) {
	T temp = b;
	b = a;
	a = temp;
};

/*
 * Computes sqrt(a^2+b^2) without destructive underflow or overflow
 */
inline double pythag(const double a, const double b) {
	double aa = fabs(a), ab = fabs(b);
	return (aa > ab ? aa*sqrt(1.0+(ab/aa)*(ab/aa)) :
		ab == 0.0 ? 0.0 : ab*sqrt(1.0+(aa/ab)*(aa/ab)));
}

/*
 * Computes n choose k
 */
inline int choose(const int n, const int k) {
	unsigned int i, j = n-k+1, c = j;
	if (k < 0 || k > n)
		return 0;
	if (k == 0 || k == n)
		return 1;
	for (i = 2; i <= k; i++)
		c = c*++j/i;
	return c;
}

#endif // MATHPLUS_H
