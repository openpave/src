/**************************************************************************

	FFT.H - Templated Fast Fourier Transform functions

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
		This header file implements some basic templated FFT functions,
		which allow one to perform forward and inverse FFTs, and
		convolutions.  They are templated to use a fixed size, so one must
		know the size of FFT you are planning on doing at compile time.

	Design:
		These are based on DJBFFT, which is in the public domain.  The code
		has been converted to use templates, so it is more generic, while
		still being about as fast.

	History:
		2008/09/12 - Initial check-in.

**************************************************************************/

#ifndef __FFT_H
#define __FFT_H

#include <mathplus.h>

namespace OP {

typedef struct {
	double re;
	double im;
} complex;

template<unsigned A, unsigned B, unsigned N = 34, unsigned M = 2>
struct SIN {
	static inline double series() {
		return 1-(A*M_PI/B)*(A*M_PI/B)/M/(M+1)
		       *SIN<A,B,N,M+2>::series();
	}
	static inline double sin() {
		return (A*M_PI/B)*SIN<A,B>::series();
	}
};
template<unsigned A, unsigned B, unsigned N>
struct SIN<A,B,N,N> {
	static inline double series() { return 1.0f; }
};

inline void
TRANSFORM(complex & a0, complex & a1, complex & a2, complex & a3,
		const double & wr,const double & wi)
{
	double t1, t2, t3, t4, t5, t6, t7, t8;
	t1 = a0.re - a2.re;
	t2 = a0.im - a2.im;
	t3 = a1.re - a3.re;
	t4 = a1.im - a3.im;
	a0.re += a2.re;
	a0.im += a2.im;
	a1.re += a3.re;
	a1.im += a3.im;
	t5 = t1 - t4;
	t6 = t2 + t3;
	t7 = t1 + t4;
	t8 = t2 - t3;
	a2.re = t5 * wr - t6 * wi;
	a2.im = t6 * wr + t5 * wi;
	a3.re = t7 * wr + t8 * wi;
	a3.im = t8 * wr - t7 * wi;
}

template<unsigned N>
inline void
fftc(complex * a)
{
	const unsigned n = N/8;
	TRANSFORM(a[0],a[2*n],a[4*n],a[6*n],1.0,0.0);
	TRANSFORM(a[n],a[3*n],a[5*n],a[7*n],M_SQRT1_2,M_SQRT1_2);
	double wt = -SIN<1,N>::sin(), wpr = -2.0*wt*wt, wpi = -SIN<2,N>::sin();
	double wr = 1.0, wi = 0.0;
	for (unsigned i = 1; i < n; i++) {
		wt = wr; wr += wr*wpr - wi*wpi; wi += wi*wpr + wt*wpi;
		TRANSFORM(a[    i],a[2*n+i],a[4*n+i],a[6*n+i],wr,-wi);
		TRANSFORM(a[2*n-i],a[4*n-i],a[6*n-i],a[8*n-i],-wi,wr);
	}
	fftc<N/2>(a);
	fftc<N/4>(a+  N/2);
	fftc<N/4>(a+3*N/4);
}
template<>
inline void
fftc<2>(complex * a)
{
	double t0 = a[0].re - a[1].re, t1 = a[0].im - a[1].im;
	a[0].re += a[1].re; a[0].im += a[1].im;
	a[1].re  = t0;      a[1].im  = t1;
}
template<>
inline void
fftc<4>(complex * a)
{
	TRANSFORM(a[0],a[1],a[2],a[3],1.0,0.0);
	fftc<2>(a);
}

inline void
R(double & a0, double & a1, double & b0, double & b1,
		const double & wr, const double & wi)
{
	double t1, t2, t3, t4;
	t1 = a0 - a1;
	t2 = b0 - b1;
	t3 = a0 + a1;
	t4 = b0 + b1;
	b0 = t1 * wr - t2 * wi;
	b1 = t2 * wr + t1 * wi;
	a0 = t3;
	a1 = t4;
}

template<unsigned N>
inline void
fftr(double * a)
{
	const unsigned n = N/8;
	R(a[  0],a[    1],a[4*n],a[4*n+1],1.0,0.0);
	R(a[2*n],a[2*n+1],a[6*n],a[6*n+1],M_SQRT1_2,M_SQRT1_2);
	double wt = -SIN<1,N>::sin(), wpr = -2.0*wt*wt, wpi = -SIN<2,N>::sin();
	double wr = 1.0, wi = 0.0;
	for (unsigned i = 1; i < n; i++) {
		wt = wr; wr += wr*wpr - wi*wpi; wi += wi*wpr + wt*wpi;
		R(a[    2*i],a[    1+2*i],a[4*n+2*i],a[4*n+1+2*i],wr,-wi);
		R(a[4*n-2*i],a[4*n+1-2*i],a[8*n-2*i],a[8*n+1-2*i],-wi,wr);
	}
	fftr<N/2>(a);
	fftc<N/4>(reinterpret_cast<complex *>(a+N/2));
}
template<>
inline void
fftr<2>(double * a)
{
	double t1 = a[0] + a[1], t2 = a[0] - a[1];
	a[0] = t1; a[1] = t2;
}
template<>
inline void
fftr<4>(double * a)
{
	R(a[0],a[1],a[2],a[3],1.0,0.0);
	fftr<2>(a);
}

inline void
UNTRANSFORM(complex & a0, complex & a1, complex & a2, complex & a3,
		const double & wr, const double & wi)
{
	double t1, t2, t3, t4, t5, t6, t7, t8;
	t5 = a2.re * wr + a2.im * wi;
	t6 = a2.im * wr - a2.re * wi;
	t7 = a3.re * wr - a3.im * wi;
	t8 = a3.im * wr + a3.re * wi;
	t1 = t5 + t7;
	t2 = t6 + t8;
	t3 = t6 - t8;
	t4 = t7 - t5;
	a2.re = a0.re - t1;
	a2.im = a0.im - t2;
	a3.re = a1.re - t3;
	a3.im = a1.im - t4;
	a0.re += t1;
	a0.im += t2;
	a1.re += t3;
	a1.im += t4;
}

template<unsigned N>
inline void
fftc_un(complex * a)
{
	fftc_un<N/2>(a);
	fftc_un<N/4>(a+  N/2);
	fftc_un<N/4>(a+3*N/4);
	const unsigned n = N/8;
	UNTRANSFORM(a[0],a[2*n],a[4*n],a[6*n],1.0,0.0);
	UNTRANSFORM(a[n],a[3*n],a[5*n],a[7*n],M_SQRT1_2,M_SQRT1_2);
	double wt = -SIN<1,N>::sin(), wpr = -2.0*wt*wt, wpi = -SIN<2,N>::sin();
	double wr = 1.0, wi = 0.0;
	for (unsigned i = 1; i < n; i++) {
		wt = wr; wr += wr*wpr - wi*wpi; wi += wi*wpr + wt*wpi;
		UNTRANSFORM(a[    i],a[2*n+i],a[4*n+i],a[6*n+i],wr,-wi);
		UNTRANSFORM(a[2*n-i],a[4*n-i],a[6*n-i],a[8*n-i],-wi,wr);
	}
}
template<>
inline void
fftc_un<2>(complex * a)
{
	double t0 = a[0].re - a[1].re, t1 = a[0].im - a[1].im;
	a[0].re += a[1].re; a[0].im += a[1].im;
	a[1].re  = t0;      a[1].im  = t1;
}
template<>
inline void
fftc_un<4>(complex *a)
{
	fftc_un<2>(a);
	UNTRANSFORM(a[0],a[1],a[2],a[3],1.0,0.0);
}

inline void
V(double & a0, double & a1, double & b0, double & b1,
		const double & wr, const double & wi)
{
	double t1, t2, t3, t4, t5, t6;
	t5 = b0 * wr + b1 * wi;
	t6 = b1 * wr - b0 * wi;
	t1 = a0 + t5;
	t2 = a0 - t5;
	t3 = a1 + t6;
	t4 = a1 - t6;
	a0 = t1;
	a1 = t2;
	b0 = t3;
	b1 = t4;
}

template<unsigned N>
inline void
fftr_un(double * a)
{
	fftc_un<N/4>(reinterpret_cast<complex *>(a + N/2));
	fftr_un<N/2>(a);
	const unsigned n = N/8;
	V(a[  0],a[    1],a[4*n],a[4*n+1],1.0,0.0);
	V(a[2*n],a[2*n+1],a[6*n],a[6*n+1],M_SQRT1_2,M_SQRT1_2);
	double wt = -SIN<1,N>::sin(), wpr = -2.0*wt*wt, wpi = -SIN<2,N>::sin();
	double wr = 1.0, wi = 0.0;
	for (unsigned i = 1; i < n; i++) {
		wt = wr; wr += wr*wpr - wi*wpi; wi += wi*wpr + wt*wpi;
		V(a[    2*i],a[    1+2*i],a[4*n+2*i],a[4*n+1+2*i],wr,-wi);
		V(a[4*n-2*i],a[4*n+1-2*i],a[8*n-2*i],a[8*n+1-2*i],-wi,wr);
	}
}
template<>
inline void
fftr_un<2>(double * a)
{
	double t1 = a[0] + a[1], t2 = a[0] - a[1];
	a[0] = t1; a[1] = t2;
}
template<>
inline void
fftr_un<4>(double * a)
{
	fftr_un<2>(a);
	V(a[0],a[1],a[2],a[3],1.0,0.0);
}

template<unsigned N>
inline void
fftc_scale(complex * a)
{
	const double n = 1.0/N;
	for (unsigned i = 0; i < N; i++)
		a[i].re *= n, a[i].im *= n;
}

/* n even, n >= 2 */
template<unsigned N>
inline void
fftr_scale(double * a)
{
	a[0] /= N; a[1] /= N;
	const double n = 2.0/N;
	for (unsigned i = 2; i < N; i++)
		a[i] *= n;
}

template<unsigned N>
inline void
fftc_mul(complex * a, complex * b)
{
	for (unsigned i = 0; i < N; i++) {
		double t1 = a[i].re, t2 = a[i].im;
		a[i].re = t1 * b[i].re - t2 * b[i].im;
		a[i].im = t1 * b[i].im + t2 * b[i].re;
	}
}

template<unsigned N>
inline void
fftr_mul(double * a, double * b)
{
	double t0 = a[0]*b[0], t1 = a[1]*b[1];
	fftc_mul<N/2>(reinterpret_cast<complex *>(a),
	              reinterpret_cast<complex *>(b));
	a[0] = t0; a[1] = t1;
}
template<>
inline void
fftr_mul<2>(double * a, double * b)
{
	a[0] *= b[0];
	a[1] *= b[1];
}

} // namespace OP

#endif // FFT_H
