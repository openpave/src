/**************************************************************************

	LINALG.H - Simple linear algebra (matrix and vector) operations

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
		This header contains prototypes for several basic linear algebra
		opertations, implemented in C style functions.  Most are based on
		reference implementations from online sources.

	Design:
		These provide the basic building blocks for doing linear algebra.

	History:
		1993       - Created by Jeremy Lea <reg@openpave.org>
		2007/09/04 - Complete redesign (again).
		2016/02/01 - Split from matrix.h

**************************************************************************/

#ifndef __LINALG_H
#define __LINALG_H

/*
 * The minimum condition number for SVD and eigenvalue decompositions.
 * NOTE: This is recorded here to show the default.  It cannot be changed
 * without recompiling.
 */
#define EIG_TOL 10e-18
/*
 * The minimum error in the result of an equals operation to trigger
 * a refinement, and the maximum number of refinements.
 */
#define ERR_TOL 10e-6
#define ITER_MAX 0

/*
 * Returns the size of the triangular matrix storage array
 */
#define T_SIZE(n)    ((n)*((n)+1)/2)
/*
 * Returns the index into the triangular matrix storage array
 * This is in column major format, so loops over i are efficient.
 */
#define T_IDX(i,j)   ((j)*((j)+1)/2+(i))

/*
 * Returns the size of the special banded matrix storage array
 */
#define B_SIZE(n,w)    (((w)+1)*(n)-(w)*((w)+1)/2)
/*
 * Returns the index into the special banded matrix storage array
 * This is in column major format, so loops over i are efficient.
 */
#define B_IDX(n,w,i,j) ((j) <= (w) ? (j)*((j)+1)/2+(i) : \
                                     ((j)+1)*(w)+(i)-(w)*((w)+1)/2)

void transpose(const unsigned m, const unsigned n, double * A);
void orth_gs(const unsigned n, double * Q);
void equ_gauss(const unsigned n, const double * A, const double * b,
	double * x);
double inv_mul_gauss(const unsigned n, const unsigned m, double * A, double * B);
void decmp_lu(const unsigned n, double * A, unsigned * idx, int & d);
void bksub_lu(const unsigned n, const double * A, const unsigned * idx,
	double * b, const unsigned m = 1, const unsigned c = 0);
void equ_lu(const unsigned n, const double * A, const double * b, double * x);
double inv_mul_lu(const unsigned n, const unsigned m, double * A, double * B);
void inv_lu(const unsigned n, double * A);
void decmp_chol(const unsigned n, double * A);
void decmp_chol_tri(const unsigned n, double * A);
void decmp_chol(const unsigned n, const unsigned w, double * A);
void bksub_chol(const unsigned n, const double * A, double * b,
	const unsigned m = 1, const unsigned c = 0);
void bksub_chol(const unsigned n, const unsigned w, const double * A, double * b,
	const unsigned m = 1, const unsigned c = 0);
void equ_chol(const unsigned n, const double * A, const double * b, double * x);
void equ_chol(const unsigned n, const unsigned w, const double * A, const double * b,
	double * x);
void inv_chol(const unsigned n, double * A);
void decmp_ldl(const unsigned n, double * A);
void bksub_ldl(const unsigned n, const double * A, double * b,
	const unsigned m = 1, const unsigned c = 0);
void equ_ldl(const unsigned n, const double * A, const double * b, double * x);
void decmp_svd(const unsigned m, const unsigned n, double * A, double * W, double * V);
void bksub_svd(const unsigned m, const unsigned n, const double * U, const double * W,
	const double * V, double * b, const unsigned p = 1, const unsigned c = 0);
void equ_svd(const unsigned n, const double * A, const double * b, double * x);
void inv_svd(const unsigned n, double * A);
void orth_svd(const unsigned n, double * Q);
bool decmp_qr(const unsigned n, double * A, double * s, double * d);
void bksub_qr(const unsigned n, const double * A, const double * s,
	const double * d, double * b, const unsigned m = 1, const unsigned c = 0);
void tridiag_hh(const unsigned n, double * A, double * d, double * e);
void eig_tri_ql(const unsigned n, double * d, double * e, double * A);
void eig_ql(const unsigned n, double * A, double * d, bool sorted = true);
void bksub_eig(const unsigned n, const double * Q, const double * d, double * b,
	const unsigned p = 1, const unsigned c = 0);
void equ_eig(const unsigned n, const double * A, const double * b, double * x);
void inv_eig(const unsigned n, double * A);

#endif // LINALG_H
