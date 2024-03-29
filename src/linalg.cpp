/*************************************************************************

	LINALG.CPP - Implementation for linear algebra operations.

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

	See LINALG.H.

	History:
		2006/07/19 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#include <memory>
#include <stdexcept>
#include <cstring>
#include "autodelete.h"
#include "mathplus.h"
#include "linalg.h"

namespace OP {

/*
 * In-place transpose a rectangular matrix.
 */
void
transpose(unsigned m, unsigned n, double * A) noexcept
{
	for (unsigned s = 0, j = 0, i = 0; s < m * n; j = ++s, i = 0) {
		do {
			i++;
			j = (j%m)*n+j/m;
		} while (j > s);
		if (j < s || i == 1)
			continue;
		const double tmp = A[j = s];
		do {
			i = (j%m)*n+j/m;
			A[j] = (i == s ? tmp : A[i]);
			j = i;
		} while (j > s);
	}
}

/*
 * Orthonormalize the nxn matrix Q, using the Gramm-Schmidt algorithm.
 */
void
orth_gs(unsigned n, double * Q) noexcept
{
	double r;
	unsigned i, j, k;

	for (k = 0; k < n; k++) {
		for (i = 0, r = 0.0; i < n; i++)
			r += pow(Q[i*n+k],2);
		if ((r = sqrt(r)) == 0.0)
			break;
		for (i = 0; i < n; i++)
			Q[i*n+k] /= r;
		for (j = k+1; j < n; j++) {
			for (i = 0, r = 0.0; i < n; i++)
				r += Q[i*n+k]*Q[i*n+j];
			for (i = 0; i < n; i++)
				Q[i*n+j] -= r*Q[i*n+k];
		}
	}
}

/*
 * This function returns the solution of Ax = b
 *
 * This function employs Gaussian elimination with full pivoting.
 */
void
equ_gauss(unsigned n, const double * A, const double * b,
	      double * x)
{
	unsigned i, j, k;

	if (n == 0)
		return;
	autodelete<double> a(new double[n*n]);
	// avoid destroying A, B by copying them to a, x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	for (i = 0; i < n; i++) {
		double pvt = 0.0;
		for (j = i; i < n && j < n; j++) {
			if (fabs(pvt = a[j*n+i]) >= DBL_EPSILON)
				break;
		}
		if (j == n)
			throw std::range_error("Singular matrix in equ_gauss()!");
		if (j != i) {
			for (k = 0; k < n; k++)
				swap(a[j*n+k],a[i*n+k]);
			swap(x[j],x[i]);
		}
		for (k = i+1; i < n && k < n; k++) {
			double tmp = a[k*n+i]/pvt;
			for (j = i+1; i < n && j < n; j++)
				a[k*n+j] -= tmp*a[i*n+j];
			x[k] -= tmp*x[i];
		}
	}
	for (i = n; i > 0; i--) {
		for (j = i; i < n && j < n; j++)
			x[i-1] -= a[(i-1)*n+j]*x[j];
		x[i-1] /= a[(i-1)*n+(i-1)];
	}
}

/*
 * This function returns the nxm matrix X = A^-1*B.  Both A and B are destroyed...
 * The result is returned in B, and the determinant in the return value.
 *
 * This function employs Gaussian elimination with full pivoting.
 */
double
inv_mul_gauss(unsigned n, unsigned m, double * A, double * B)
{
	double det = 1.0;
	unsigned i, j, k;

	for (i = 0; i < n; i++) {
		double pvt = 0.0;
		for (j = i; j < n; j++) {
			if (fabs(pvt = A[j*n+i]) >= DBL_EPSILON)
				break;
		}
		if (j == n)
			throw std::range_error("Singular matrix in inv_mul_gauss()!");
		if (j != i) {
			for (k = 0; k < n; k++)
				swap(A[j*n+k],A[i*n+k]);
			for (k = 0; k < m; k++)
				swap(B[j*m+k],B[i*m+k]);
		}
		det *= pvt;
		for (k = i+1; i < n && k < n; k++) {
			const double tmp = A[k*n+i]/pvt;
			for (j = i+1; i < n && j < n; j++)
				A[k*n+j] -= tmp*A[i*n+j];
			for (j = 0; j < m; j++)
				B[k*m+j] -= tmp*B[i*m+j];
		}
	}
	for (i = n; i > 0; i--) {
		for (j = i; j < n; j++)
			for (k = 0; k < m; k++)
				B[(i-1)*m+k] -= A[(i-1)*n+j]*B[j*m+k];
		for (k = 0; k < m; k++)
			B[(i-1)*m+k] /= A[(i-1)*n+(i-1)];
	}
	return det;
}

/*
 * Perform LU decompostion of an nxn matrix A.
 *
 * The algorithm comes from the net, based on the implementation in
 * Numerical Recipes, with my own improvements to the numerical accuracy.
 *
 * idx holds the row interchange index, and d is negative if an odd number
 * of interchanges have been performed (this is for the determinant calc).
 */
void
decmp_lu(unsigned n, double * A, unsigned * idx, int & d)
{
	unsigned i, j, k;

	if (n == 0)
		return;
	autodelete<double> work(new double[n]);
	d = 1;
	// compute the LU decomposition of a row permutation of matrix a;
	// the permutation itself is saved in idx[]
	for (i = 0; i < n; i++) {
		double tmp, max = 0.0;
		for (j = 0; j < n; j++)
			if ((tmp = fabs(A[i*n+j])) > max)
				max = tmp;
		if (max < DBL_EPSILON)
			throw std::range_error("Singular matrix A in decmp_lu()!");
		work[i] = max;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			for (k = 0; k < j; k++)
				A[j*n+i] -= A[j*n+k]*A[k*n+i];
		}
		double tmp = 0.0, max = 0.0;
		for (j = i; j < n; j++) {
			for (k = 0; k < i; k++)
				A[j*n+i] -= A[j*n+k]*A[k*n+i];
			if ((tmp = fabs(A[j*n+i])/work[j]) >= max) {
				max = tmp; idx[i] = j;
			}
		}
		if (i != idx[i]) {
			for (j = 0; j < n; j++)
				swap(A[idx[i]*n+j],A[i*n+j]);
			d *= -1;
			swap(work[idx[i]],work[i]);
		}
		if (A[i*n+i] == 0.0) {
			throw std::range_error("Singular matrix A in decmp_lu()!");
		}
		for (j = i+1; j < n; j++)
			A[j*n+i] /= A[i*n+i];
	}
}

/*
 * This function performs forward/back substitution of an LU decomposed
 * matrix, to solve a particular system.
 *
 * The algorithm comes from the net, based on the implementation in
 * Numerical Recipes, with my own improvements to the numerical accuracy.
 *
 * m is the number of columns in b, and c is the column index.
 */
void
bksub_lu(unsigned n, const double * A, const unsigned * idx,
	     double * b, unsigned m, unsigned c) noexcept
{
	unsigned i, j, k;

	for (i = 0, k = 0; i < n; i++) {
		j = idx[i];
		double sum = b[j*m+c];
		b[j*m+c] = b[i*m+c];
		if (k > 0)
			for (j = k-1; k <= i && j < i; j++)
				sum -= A[i*n+j]*b[j*m+c];
		else if (sum != 0.0)
			k = i+1;
		b[i*m+c] = sum;
	}
	for (i = n; i > 0; i--) {
		for (j = i; i < n && j < n; j++)
			b[(i-1)*m+c] -= A[(i-1)*n+j]*b[j*m+c];
		b[(i-1)*m+c] /= A[(i-1)*n+(i-1)];
	}
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LU decomposition followed by forward/back
 * substitution, with a single step refinement.
 */
void
equ_lu(unsigned n, const double * A, const double * b, double * x)
{
	unsigned i, j;
	int d;

	if (n == 0)
		return;
	autodelete<unsigned> idx(new unsigned[n]);
	autodelete<double> a(new double[n*n]);
	autodelete<double> r(new double[n]);
	// avoid destroying A, B by copying them to a, x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	decmp_lu(n,a,idx,d);
	bksub_lu(n,a,idx,x);
	for (i = 0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
	}
	bksub_lu(n,a,idx,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];
}

/*
 * This function returns the nxm matrix X = A^-1*B.  Both A and B are destoryed...
 * The result is returned in B, and the determinant in the return value.
 */
double
inv_mul_lu(unsigned n, unsigned m, double * A, double * B)
{
	double det = 0.0;
	unsigned i;
	int d;

	if (n == 0)
		return 0.0;
	autodelete<unsigned> idx(new unsigned[n]);
	decmp_lu(n,A,idx,d);
	for (i = 0, det = d; i < n; i++)
		det *= A[i*n+i];
	for (i = 0; i < m; i++)
		bksub_lu(n,A,idx,B,m,i);
	return det;
}

/*
 * Matrix inverse of the real nxn matrix A using LU decomposition.
 */
void
inv_lu(unsigned n, double * A)
{
	unsigned i;
	int d;

	if (n == 0)
		return;
	autodelete<unsigned> idx(new unsigned[n]);
	autodelete<double> a(new double[n*n]);
	memcpy(a,A,sizeof(double)*n*n);
	memset(A,0,sizeof(double)*n*n);
	for (i = 0; i < n; i++)
		A[i*n+i] = 1.0;
	decmp_lu(n,a,idx,d);
	for (i = 0; i < n; i++)
		bksub_lu(n,a,idx,A,n,i);
}

/*
 * Cholesky Decomposition of symmetic postive definite nxn matrix A.
 *
 * Returns L in the lower triangle of A.
 */
void
decmp_chol(unsigned n, double * A)
{
	unsigned i, j, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			sum = A[j*n+i];
			for (k = 0; k < i; k++)
				sum -= A[i*n+k]*A[j*n+k];
			if (i == j) {
				if (sum <= 0.0)
					throw std::runtime_error("Non-positive definite matrix in decmp_chol()!");
				A[i*n+i] = 1.0/sqrt(sum);
			} else
				A[j*n+i] = sum*A[i*n+i];
		}
	}
}

/*
 * Cholesky Decomposition of symmetic postive definite triangular matrix A.
 *
 * Returns U in the trianglar matrix.
 */
void
decmp_chol_tri(unsigned n, double * A)
{
	unsigned i, j, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			sum = A[T_IDX(i,j)];
			for (k = 0; k < i; k++)
				sum -= A[T_IDX(k,i)]*A[T_IDX(k,j)];
			if (i == j) {
				if (sum <= 0.0)
					throw std::runtime_error("Non-positive definite matrix in decmp_chol_tri()!");
				A[T_IDX(i,i)] = sqrt(sum);
			} else
				A[T_IDX(i,j)] = sum/A[T_IDX(i,i)];
		}
	}
}

/*
 * Cholesky Decomposition of postive definite nxn matrix A.  The matrix is banded,
 * with bandwidth 2*w+1, and only the upper triangle is stored.
 *
 * Returns U in the upper triangle of A.
 */
void
decmp_chol(unsigned n, unsigned w, double * A)
{
	unsigned i, j, k;

	for (i = 0; i < n; i++) {
		for (j = i; j <= i+w && j < n; j++) {
			double sum = A[B_IDX(n,w,i,j)];
			for (k = (w>j?0:j-w); k < i; k++)
				sum -= A[B_IDX(n,w,k,i)]*A[B_IDX(n,w,k,j)];
			if (i == j) {
				if (sum <= 0)
					throw std::runtime_error("Non-positive definite matrix in decmp_chol()!");
				A[B_IDX(n,w,i,i)] = 1.0/sqrt(sum);
			} else
				A[B_IDX(n,w,i,j)] = sum*A[B_IDX(n,w,i,i)];
		}
	}
}

/*
 * This function performs forward/back substitution of a Cholesky decomposed
 * matrix, to solve a particular system.
 *
 * m is the number of columns in b, and c is the column index.
 */
void
bksub_chol(unsigned n, const double * A,
	       double * b, unsigned m, unsigned c) noexcept
{
	unsigned i, k;

	for (i = 0; i < n; i++) {
		for (k = 0; k < i; k++)
			b[i*m+c] -= A[i*n+k]*b[k*m+c];
		b[i*m+c] *= A[i*n+i];
	}
	for (i = n; i > 0; i--) {
		for (k = i; k < n; k++)
			b[(i-1)*m+c] -= A[k*n+(i-1)]*b[k*m+c];
		b[(i-1)*m+c] *= A[(i-1)*n+(i-1)];
	}
}

/*
 * This function performs forward/back substitution of a Cholesky decomposed
 * matrix, to solve a particular system.
 */
void
bksub_chol(unsigned n, unsigned w, const double * A,
	       double * b, unsigned m, unsigned c) noexcept
{
	unsigned i, k;

	for (i = 0; i < n; i++) {
		for (k = (w>i?0:i-w); k < i; k++)
			b[i*m+c] -= A[B_IDX(n,w,k,i)]*b[k*m+c];
		b[i*m+c] *= A[B_IDX(n,w,i,i)];
	}
	for (i = n; i > 0; i--) {
		for (k = i; k < i+w && k < n; k++)
			b[(i-1)*m+c] -= A[B_IDX(n,w,i-1,k)]*b[k*m+c];
		b[(i-1)*m+c] *= A[B_IDX(n,w,i-1,i-1)];
	}
}

/*
 * Solve of Ax = b by Cholesky decomposition followed by forward/back
 * substitution.
 */
void
equ_chol(unsigned n, const double * A, const double * b, double * x)
{
	//unsigned i, j;

	if (n == 0)
		return;
	autodelete<double> a(new double[n*n]);
	//autodelete<double> r(new double[n]);
	// avoid destroying A and B by copying them to a and x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	decmp_chol(n,a);
	bksub_chol(n,a,x);
	/*
	for (i = 0, dot = 0.0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
		dot += r[i]*r[i];
	}
	if (dot > ERR_TOL)
		return;
	bksub_chol(n,a,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];
	*/
}

/*
 * This function returns the solution of Ax = b, where A is banded.
 *
 * The function employs Cholesky decomposition followed by forward/back
 * substitution, with a single step refinement.
 */
void
equ_chol(unsigned n, unsigned w, const double * A,
	     const double * b, double * x)

{
	unsigned i, j, iter = 0;
	double dot, c1, y1, t1, c, y, t;

	if (n == 0)
		return;
	autodelete<double> a(new double[B_SIZE(n,w)]);
	autodelete<double> r(new double[n]);
	// avoid destroying A and B by copying them to a and x resp.
	memcpy(a,A,sizeof(double)*B_SIZE(n,w));
	memcpy(x,b,sizeof(double)*n);
	decmp_chol(n,w,a);
	bksub_chol(n,w,a,x);
	while (1) {
		dot = 0.0; c1 = 0.0;
		for (i = 0; i < n; i++) {
			r[i] = -b[i];
			for (j = (w>i?0:i-w), c = 0.0; j < i; j++) {
				y = fma(A[B_IDX(n,w,j,i)],x[j],-c); t = r[i] + y;
				c = t - r[i] - y; r[i] = t;
			}
			for (j = i; j <= i+w && j < n; j++) {
				y = fma(A[B_IDX(n,w,i,j)],x[j],-c); t = r[i] + y;
				c = t - r[i] - y; r[i] = t;
			}
			y1 = fma(r[i],r[i],-c1); t1 = dot + y1;
			c1 = t1 - dot - y1; dot = t1;
		}
		if (++iter > ITER_MAX || sqrt(dot) <= ERR_TOL)
			break;
		bksub_chol(n,w,a,r);
		for (i = 0; i < n; i++)
			x[i] -= r[i];
	}
}

/*
 * Matrix inverse of the real positive definite nxn matrix A using
 * Cholesky decomposition.
 */
void
inv_chol(unsigned n, double * A)
{
	unsigned i, j, k;
	double sum;

	decmp_chol(n,A);
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n; j++) {
			for (k = i, sum = 0.0; i < j && k < j; k++)
				sum -= A[j*n+k]*A[k*n+i];
			A[j*n+i] = sum*A[j*n+j];
		}
	}
	for (i = 0; i < n; i++) {
		A[i*n+i] *= A[i*n+i];
		for (j = i+1; j < n; j++)
			A[i*n+i] += A[j*n+i]*A[j*n+i];
		for (j = i+1; j < n; j++) {
			A[i*n+j] = 0.0;
			for (k = j; k < n; k++)
				A[i*n+j] += A[k*n+i]*A[k*n+j];
			A[j*n+i] = A[i*n+j];
		}
	}
}

/*
 * LDL^T (Square root free Cholesky) Decomposition of postive definite
 * nxn matrix A.
 *
 * Returns L in the lower triangle of A (all diagonal elements are 1),
 * with D in the diagonal.
 */
void
decmp_ldl(unsigned n, double * A)
{
	unsigned i, j, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (k = 0, sum = 0.0; k < i; k++)
			sum += A[k*n+k]*(A[i*n+k]*A[i*n+k]);
		A[i*n+i] -= sum;
		if (fabs(A[i*n+i]) < DBL_EPSILON)
			throw std::runtime_error("Indefinite matrix in decmp_ldl()!");
		for (j = i+1; j < n; j++) {
			for (k = 0, sum = 0.0; k < i; k++) {
				sum += A[k*n+k]*(A[i*n+k]*A[j*n+k]);
			}
			A[j*n+i] = (A[j*n+i]-sum)/A[i*n+i];
		}
	}
}

/*
 * This function performs forward/back substitution of a LDL^T decomposed
 * matrix, to solve a particular system.
 *
 * m is the number of columns in b, and c is the column index.
 */
void
bksub_ldl(unsigned n, const double * A,
	      double * b, unsigned m, unsigned c) noexcept
{
	unsigned i, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (k = 0, sum = 0.0; k < i; k++)
			sum += A[i*n+k]*b[k*m+c];
		b[i*m+c] -= sum;
	}
	for (i = n; i > 0; i--) {
		b[(i-1)*m+c] /= A[(i-1)*n+(i-1)];
		for (k = i, sum = 0.0; k < n; k++)
			sum += A[k*n+(i-1)]*b[k*m+c];
		b[(i-1)*m+c] -= sum;
	}
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LDL^T decomposition followed by forward/back
 * substitution, (with a single step refinement).
 */
void
equ_ldl(unsigned n, const double * A, const double * b, double * x)
{
	unsigned i, j;

	if (n == 0)
		return;
	autodelete<double> a(new double[n*n]);
	autodelete<double> r(new double[n]);
	// avoid destroying A and B by copying them to a and x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	decmp_ldl(n,a);
	bksub_ldl(n,a,x);
	for (i = 0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
	}
	bksub_ldl(n,a,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];
}

/*
 * Given a mxn matrix A, this routine computes its singular value
 * decomposition, A = U * W * V'. The matrix U replaces A on output.
 * The diagonal matrix of singular values, W, is output as a vector W.
 * The matrix V (not the transpose of V) is output as V.
 */
void
decmp_svd(unsigned m, unsigned n, double * A,
	      double * W, double * V)
{
	double F, G = 0.0, H;
	double C, S, X, Y;
	double scale = 0.0;
	double anorm = 0.0;
	unsigned iter, i, j, k, q;

	if (m == 0 || n == 0)
		return;
	autodelete<double> rv1(new double[n]);
	// Householder reduction to bidiagonal form.
	for (i = 0; i < n; i++) {
		rv1[i] = scale*G;
		for (k = i, scale = 0.0; i < m && k < m; k++)
			scale += fabs(A[k*n+i]);
		if (scale != 0.0) {
			for (k = i, S = 0.0; i < m && k < m; k++) {
				A[k*n+i] /= scale;
				S += A[k*n+i]*A[k*n+i];
			}
			G = (A[i*n+i] > 0.0 ? -1 : 1)*sqrt(S);
			H = A[i*n+i]*G - S;
			A[i*n+i] -= G;
			for (j = i+1; i < n && j < n; j++) {
				for (k = i, S = 0.0; i < n && k < m; k++)
					S += A[k*n+i]*A[k*n+j];
				for (k = i; i < m && k < m; k++)
					A[k*n+j] += S*A[k*n+i]/H;
			}
			for (k = i; i < m && k < m; k++)
				A[k*n+i] *= scale;
			W[i] = scale*G;
		}
		for (k = i+1, scale = 0.0; i < n && i < m && k < n; k++)
			scale += fabs(A[i*n+k]);
		if (scale != 0.0) {
			for (k = i+1, S = 0.0; i < n && k < n; k++) {
				A[i*n+k] /= scale;
				S += A[i*n+k]*A[i*n+k];
			}
			G = (A[i*n+i+1] > 0.0 ? -1 : 1)*sqrt(S);
			H = A[i*n+i+1]*G - S;
			A[i*n+i+1] -= G;
			for (k = i+1; i < n && k < n; k++)
				rv1[k] = A[i*n+k] / H;
			for (j = i+1; i < m && j < m; j++ ) {
				for (k = i+1, S = 0.0; i < n && k < n; k++)
					S += A[j*n+k]*A[i*n+k];
				for (k = i+1; i < n && k < n; k++ )
					A[j*n+k] += S*rv1[k];
			}
			for (k = i+1; i < n && k < n; k++)
				A[i*n+k] *= scale;
		}
		const double tmp = fabs(W[i]) + fabs(rv1[i]);
		if (tmp > anorm)
			anorm = tmp;
	}
	// Accumulation of right-hand transformations.
	for (i = n; i > 0; i--) {
		if (rv1[i] != 0.0 ) {
			for (j = i; i < n && j < n; j++)
				V[j*n+(i-1)] = (A[(i-1)*n+j]/A[(i-1)*n+i])/rv1[i];
			for (j = i; i < n && j < n; j++) {
				for (k = i, S = 0.0; i < n && k < n; k++)
					S += A[(i-1)*n+k]*V[k*n+j];
				for (k = i; i < n && k < n; k++)
					V[k*n+j] += S*V[k*n+(i-1)];
			}
		}
		for (j = i; i < n && j < n; j++)
			V[(i-1)*n+j] = V[j*n+(i-1)] = 0.0;
		V[(i-1)*n+(i-1)] = 1.0;
	}
	// Accumulation of left-hand transformations.
	for (i = std::min(m,n); i > 0; i--) {
		for (j = i; i < n && j < n; j++)
			A[(i-1)*n+j] = 0.0;
		if (W[i-1] != 0.0) {
			for (j = i; i < n && j < n; j++) {
				for (k = i, S = 0.0; i < m && k < m; k++)
					S += A[k*n+(i-1)]*A[k*n+j];
				const double f = (S/A[(i-1)*n+(i-1)])/W[i-1];
				for (k = i-1; k < m; k++)
					A[k*n+j] += f*A[k*n+(i-1)];
			}
			for (j = i-1; j < m; j++)
				A[j*n+(i-1)] /= W[i-1];
		} else {
			for (j = i-1; j < m; j++)
				A[j*n+(i-1)] = 0.0;
		}
		A[(i-1)*n+(i-1)] += 1.0;
	}
	// Diagonalization of the bidiagonal form, looping over
	// singular values.
	for (i = n; i > 0; ) {
		i--; iter = 0;
		while (iter++ < 30) {
			// Test for splitting. Note that rv1[0] is always zero.
			bool flag = true;
			for (q = i+1; q > 0; ) {
				if ((fabs(rv1[--q]) + anorm) == anorm) {
					flag = false;
					break;
				} else if ((fabs(W[q-1]) + anorm) == anorm)
					break;
			}
			// Cancellation of rv1[q], if q > 0.
			S = 1.0;
			for (k = q; flag && q <= i && k <= i; k++) {
				F = S*rv1[k];
				if ((fabs(F) + anorm) != anorm) {
					G = W[k];
					W[k] = hypot(F,G);
					C = G/W[k]; S = -F/W[k];
					for (j = 0; j < m; j++) {
						G = A[j*n+q-1];
						A[j*n+q-1] = A[j*n+k]*S + G*C;
						A[j*n+k]   = A[j*n+k]*C - G*S;
					}
				}
			}
			if (q == i) {
				// Singular value is made nonnegative.
				if (W[i] < 0.0) {
					W[i] *= -1;
					for (j = 0; j < n; j++)
						V[j*n+i] *= -1;
				}
				break; // Convergence.
			}
			F = ((W[i-1]-W[i])*(W[i-1]+W[i])
					+ (rv1[i-1]-rv1[i])*(rv1[i-1]+rv1[i]))
					/ (2.0*rv1[i]*W[i-1]);
			G = (F < 0.0 ? -1 : 1)*hypot(F,1.0);
			F = ((W[q]-W[i])*(W[q]+W[i])
					+ rv1[i]*((W[i-1]/(F+G))-rv1[i])) / W[q];
			// Next QR transformation.
			X = W[q]; C = 1.0; S = 1.0;
			for (j = q; q < i && j < i; j++) {
				H = S*rv1[j+1], G = C*rv1[j+1];
				rv1[j] = hypot(F,H);
				C = F/rv1[j], S = H/rv1[j];
				F = G*S+X*C,  G = G*C-X*S;
				H = W[j+1]*S, Y = W[j+1]*C;
				for (k = 0; k < n; k++) {
					X = V[k*n+j];
					V[k*n+j]   = V[k*n+j+1]*S+X*C;
					V[k*n+j+1] = V[k*n+j+1]*C-X*S;
				}
				W[j] = hypot(F,H);
				if (W[j] != 0.0)
					C = F/W[j], S = H/W[j];
				F = S*Y+C*G, X = C*Y-S*G;
				for (k = 0; k < m; k++) {
					Y = A[k*n+j];
					A[k*n+j]   = A[k*n+j+1]*S + Y*C;
					A[k*n+j+1] = A[k*n+j+1]*C - Y*S;
				}
			}
			rv1[q] = 0.0;
			rv1[i] = F; W[i] = X;
		}
	}
	return;
}

void
bksub_svd(unsigned m, unsigned n, const double * U,
	      const double * W, const double * V,
	      double * b, unsigned p, unsigned c)
{
	unsigned i, j;

	if (m == 0 || n == 0)
		return;
	autodelete<double> tmp(new double[n]);
	for (j = 0; j < n; j++) {
		tmp[j] = 0.0;
		if (W[j] != 0.0) {
			for (i = 0; i < m; i++)
				tmp[j] += U[i*n+j]*b[i*p+c];
			tmp[j] /= W[j];
		}
	}
	for (j = 0; j < n; j++) {
		b[j*p+c] = 0.0;
		for (i = 0; i < n; i++)
			b[j*p+c] += V[j*n+i]*tmp[i];
	}
	return;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs SVD decomposition followed by forward/back
 * substitution, with refinement.
 */
void
equ_svd(unsigned n, const double * A, const double * b, double * x)
{
	double max;
	unsigned i, j;

	if (n == 0)
		return;
	autodelete<double> U(new double[n*n]);
	autodelete<double> W(new double[n]);
	autodelete<double> V(new double[n*n]);
	autodelete<double> r(new double[n]);
	// avoid destroying A, B by copying them to U, x resp.
	memcpy(U,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	decmp_svd(n,n,U,W,V);
	for (i = 0, max = 0.0; i < n; i++)
		if (W[i] > max)
			max = W[i];
	for (i = 0, max *= EIG_TOL; i < n; i++)
		if (W[i] < max)
			W[i] = 0.0;
	bksub_svd(n,n,U,W,V,x);
	for (i = 0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
	}
	bksub_svd(n,n,U,W,V,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];
}

/*
 * Matrix inverse of the real nxn matrix A using SVD decomposition.
 */
void
inv_svd(unsigned n, double * A)
{
	double max;
	unsigned i, j, k;

	if (n == 0)
		return;
	autodelete<double> U(new double[n*n]);
	autodelete<double> W(new double[n]);
	autodelete<double> V(new double[n*n]);
	memcpy(U,A,sizeof(double)*n*n);
	decmp_svd(n,n,U,W,V);
	for (i = 0, max = 0.0; i < n; i++)
		if (W[i] > max)
			max = W[i];
	for (i = 0, max *= EIG_TOL; i < n; i++)
		if (W[i] < max)
			W[i] = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0, A[i*n+j] = 0.0; k < n; k++) {
				if (W[k] != 0.0)
					A[i*n+j] += V[i*n+k]*U[j*n+k]/W[k];
			}
		}
	}
}

/*
 * Orthonormalize the nxn matrix Q, using the SVD decomposition.
 */
void
orth_svd(unsigned n, double * Q)
{
	if (n == 0)
		return;
	autodelete<double> W(new double[n]);
	autodelete<double> V(new double[n*n]);
	decmp_svd(n,n,Q,W,V);
}

/*
 * Performs a QR decomposition of the nxn matrix A.
 *
 * XXX: this routine is using a stupid access scheme.
 */
bool
decmp_qr(unsigned n, double * A, double * s, double * d) noexcept
{
	bool rv = true;
	unsigned i, j, k;

	if (n == 0)
		return true;
	for (k = 0; n > 1 && k < n-1; k++) {
		double sum, tmp, scale = 0.0;
		for (i = k; i < n; i++) {
			if ((tmp = fabs(A[i*n+k])) > scale)
				scale = tmp;
		}
		if (scale == 0.0) {
			rv = false;
			s[k] = d[k] = 0.0;
		} else {
			for (i = k, sum = 0.0; i < n; i++) {
				A[i*n+k] /= scale;
				sum += A[i*n+k]*A[i*n+k];
			}
			sum = (A[k*n+k] < 0.0 ? -1 : 1)*sqrt(sum);
			A[k*n+k] += sum;
			s[k] = sum*A[k*n+k];
			d[k] = -scale*sum;
			for (j = k+1; j < n; j++) {
				for (i = k, sum = 0.0; i < n; i++)
					sum += A[i*n+k]*A[i*n+j];
				sum /= s[k];
				for (i = k; i < n; i++)
					A[i*n+j] -= sum*A[i*n+k];
			}
		}
	}
	d[n-1] = A[n*n-1];
	if (d[n-1] == 0.0)
		rv = false;
	return rv;
}

/*
 * Performs a backsubstitution for a QR decomposed matrix A, solving
 * Ax = b.
 */
void
bksub_qr(unsigned n, const double * A,
	     const double * s, const double * d,
	     double * b, unsigned m, unsigned c) noexcept
{
	unsigned i, j;
	double sum;

	for (j = 0; j < n; j++) {
		for (i = j, sum = 0.0; j < n && i < n; i++)
			sum += A[i*n+j]*b[i*m+c];
		sum /= s[j];
		for (i = j; j < n && i < n; i++)
			b[i*m+c] -= sum*A[i*n+j];
	}
	for (i = n; i > 0; i--) {
		for (j = i, sum = 0.0; i < n && j < n; j++)
			sum += A[(i-1)*n+j]*b[j*m+c];
		b[(i-1)*m+c] = (b[(i-1)*m+c]-sum)/d[i-1];
	}
}

/*
 * Householder reduction of a nxn matrix A to tridiagonal form.
 */
void
tridiag_hh(unsigned n, double * A, double * d, double * e) noexcept
{
	unsigned i, j, k;
	double scale, f;

	if (n == 0)
		return;
	for (i = n-1; i > 0; i--)  {
		for (j = 0, d[i] = 0.0, scale = 0.0; j < i-1; j++)
			scale += fabs(A[i*n+j]);
		if (scale == 0.0) {
			e[i] = A[i*n+i-1];
			continue;
		}
		scale += fabs(A[i*n+i-1]);
		for (j = 0; j < i; j++) {
			A[i*n+j] /= scale;
			d[i] += A[i*n+j]*A[i*n+j];
		}
		f = (A[i*n+i-1] >= 0 ? -1 : 1)*sqrt(d[i]);
		e[i] = scale*f;
		d[i] -= A[i*n+i-1]*f;
		A[i*n+i-1] -= f;
		for (j = 0, f = 0.0; j < i; j++) {
			A[j*n+i] = A[i*n+j]/d[i];
			for (k = 0, e[j] = 0.0; k < j+1; k++)
				e[j] += A[j*n+k]*A[i*n+k];
			for (k = j+1; j < i && k < i; k++)
				e[j] += A[k*n+j]*A[i*n+k];
			e[j] /= d[i];
			f += e[j]*A[i*n+j];
		}
		for (j = 0; j < i; j++) {
			e[j] -= f/(2*d[i])*A[i*n+j];
			for (k = 0; k < j+1; k++)
				A[j*n+k] -= (A[i*n+j]*e[k] + e[j]*A[i*n+k]);
		}
	}
	d[0] = e[0] = 0.0;
	for (i = 0; i < n; i++) {
		if (d[i] != 0.0) {
			for (j = 0; j < i; j++) {
				for (k = 0, f = 0.0; k < i; k++)
					f += A[i*n+k]*A[k*n+j];
				for (k = 0; k < i; k++)
					A[k*n+j] -= f*A[k*n+i];
			}
		}
		d[i] = A[i*n+i], A[i*n+i] = 1.0;
		for (j = 0; j < i; j++)
			A[j*n+i] = A[i*n+j] = 0.0;
	}
}

/*
 * This subroutine implements the QL algorithm with
 * implicit shifts to determine the eigenvalues and
 * eigenvectors of a real symmetric tridiagonal matrix.
 */
void
eig_tri_ql(unsigned n, double * d, double * e, double * A) noexcept
{
	unsigned i, j, k, m;
	double b,c,f,g,p,r,s;

	if (n == 0)
		return;
	for (i = 0; i < n; i++) {
		while (1) {
			for (m = i; m < n-1; m++) {
				g = fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m+1])+g == g)
					break;
			}
			if (m == i)
				break;
			g = (d[i+1]-d[i])/(2.0*e[i+1]);
			r = (g > 0.0 ? 1 : -1)*sqrt(g*g+1.0);
			g = d[m] - d[i] + e[i+1]/(g+r);
			s = c = 1.0; p = 0.0;
			for (j = m; j > i; j--) {
				f = s*e[j], b = c*e[j];
				if ((e[j+1<n?j+1:0] = r = hypot(f,g)) == 0.0) {
					d[j] -= p;
					break;
				}
				c = g/r, s = f/r;
				d[j] -= p;
				r = (d[j-1]-d[j])*s+2.0*c*b;
				d[j] += p = s*r;
				g = c*r-b;
				for (k = 0; k < n; k++) {
					r = A[k*n+j];
					A[k*n+j  ] = s*A[k*n+j-1] + c*r;
					A[k*n+j-1] = c*A[k*n+j-1] - s*r;
				}
			}
			e[m+1<n?m+1:0] = 0.0;
			if (j <= i)
				d[i] -= p, e[i+1] = g;
		}
	}
}

/*
 * Obtain the eigenvalues and eigenvectors of a real symmetric matrix by
 * QL transform.
 */
void
eig_ql(unsigned n, double * A, double * d, bool sorted)
{
	double t;
	unsigned i, j, k;

	if (n == 0)
		return;
	autodelete<double> e(new double[n]);
	tridiag_hh(n,A,d,e);
	eig_tri_ql(n,d,e,A);
	// Sort the eigenvalues into descending order
	for (i = 0; sorted && i < n-1; i++) {
		for (j = i+1, t = d[k = i]; j < n; j++) {
			if (d[j] >= t)
				t = d[k = j];
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = t;
			for (j = 0; j < n; j++) {
				t=A[j*n+i];
				A[j*n+i]=A[j*n+k];
				A[j*n+k]=t;
			}
		}
	}
}

void
bksub_eig(unsigned n, const double * Q, const double * d,
	      double * b, unsigned p, unsigned c)
{
	unsigned i, j;

	if (n == 0)
		return;
	autodelete<double> tmp(new double[n]);
	memset(tmp,0,sizeof(double)*n);
	for (j = 0; j < n; j++) {
		if (d[j] != 0.0) {
			for (i = 0; i < n; i++)
				tmp[j] += Q[i*n+j]*b[i*p+c];
			tmp[j] /= d[j];
		}
	}
	for (j = 0; j < n; j++) {
		for (i = 0, b[j] = 0.0; i < n; i++)
			b[j] += Q[j*n+i]*tmp[i];
	}
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs an Eigenvalue decomposition followed by forward/back
 * substitution, with refinement.
 */
void
equ_eig(unsigned n, const double * A, const double * b, double * x)
{
	double max;
	unsigned i, j;

	if (n == 0)
		return;
	autodelete<double> Q(new double[n*n]);
	autodelete<double> d(new double[n]);
	autodelete<double> e(new double[n]);
	autodelete<double> r(new double[n]);
	memcpy(Q,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	tridiag_hh(n,Q,d,e);
	eig_tri_ql(n,d,e,Q);
	for (i = 0, max = 0.0; i < n; i++) {
		if (fabs(d[i]) > max)
			max = fabs(d[i]);
	}
	for (i = 0, max *= EIG_TOL; i < n; i++) {
		if (fabs(d[i]) < max)
			d[i] = 0.0;
	}
	bksub_eig(n,Q,d,x);
	for (i = 0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
	}
	bksub_eig(n,Q,d,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];
}

/*
 * Matrix inverse of the real symmetric nxn matrix A using eigenvalue decomposition.
 */
void
inv_eig(unsigned n, double * A)
{
	double max;
	unsigned i, j, k;

	if (n == 0)
		return;
	autodelete<double> Q(new double[n*n]);
	autodelete<double> d(new double[n]);
	autodelete<double> e(new double[n]);
	memcpy(Q,A,sizeof(double)*n*n);
	tridiag_hh(n,Q,d,e);
	eig_tri_ql(n,d,e,Q);
	for (i = 0, max = 0.0; i < n; i++) {
		if (fabs(d[i]) > max)
			max = fabs(d[i]);
	}
	for (i = 0, max *= EIG_TOL; i < n; i++) {
		if (fabs(d[i]) < max)
			d[i] = 0.0;
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0, A[i*n+j] = 0.0; k < n; k++) {
				if (d[k] != 0.0)
					A[i*n+j] += Q[i*n+k]*Q[j*n+k]/d[k];
			}
		}
	}
}

} // namespace OP
