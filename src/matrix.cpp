/*************************************************************************

	MATRIX.CPP - Implementation for matrix operations.

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

	See MATRIX.H.

	History:
		2006/07/19 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#include "matrix.h"
#include <stdio.h>

/*
 * Orthonormalize the nxn matrix Q, using the Gramm-Schmidt algorithm.
 */
void
orth_gs(const int n, double * __restrict Q)
{
	double r;
	int i, j, k;

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
bool
equ_gauss(const int n, const double * __restrict A, const double * __restrict b,
          double * __restrict x)
{
	bool rv = true;
	int i, j, k;

	double * __restrict a = new double[n*n];
	if (a == 0) {
		event_msg(EVENT_ERROR,"Out of memory in equ_gauss()!");
		rv = false;
		goto abort;
	}
	// avoid destroying A, B by copying them to a, x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	for (i = 0; i < n; i++) {
		double pvt = a[i*n+i];
		if (fabs(pvt) < DBL_EPSILON) {
			for (j = i+1; j < n; j++) {
				if (fabs(pvt = a[j*n+i]) >= DBL_EPSILON)
					break;
			}
			if (j == n) {
				event_msg(EVENT_ERROR,"Singular matrix in equ_gauss()!");
				rv = false;
				goto abort;
			}
			for (k = 0; k < n; k++)
				swap(a[j*n+k],a[i*n+k]);
			swap(x[j],x[i]);
		}
		for (k = n-1; k > i; k--) {
			double tmp = a[k*n+i]/pvt;
			for (j = n-1; j > i; j--)
				a[k*n+j] -= tmp*a[i*n+j];
			x[k] -= tmp*x[i];
		}
	}
	for (i = n-1; i >= 0; i--) {
		for (j = n-1; j > i; j--)
			x[i] -= a[i*n+j]*x[j];
		x[i] /= a[i*n+i];
	}
  abort:
	delete [] a;
	return rv;
}

/*
 * This function returns the nxm matrix X = A^-1*B.  Both A and B are destoryed...
 * The result is returned in B, and the determinant in the return value.
 *
 * This function employs Gaussian elimination with full pivoting.
 */
double
inv_mul_gauss(const int n, const int m, double * __restrict A,
              double * __restrict B)
{
	double det = 1.0;
	int i, j, k;

	for (i = 0; i < n; i++) {
		double pvt = A[i*n+i];
		if (fabs(pvt) < DBL_EPSILON) {
			for (j = i+1; j < n; j++) {
				if (fabs(pvt = A[j*n+i]) >= DBL_EPSILON)
					break;
			}
			if (j == n) {
				event_msg(EVENT_ERROR,"Singular matrix in inv_mul_gauss()!");
				det = 0.0;
				goto abort;
			}
			for (k = 0; k < n; k++)
				swap(A[j*n+k],A[i*n+k]);
			for (k = 0; k < m; k++)
				swap(B[j*m+k],B[i*m+k]);
		}
		det *= pvt;
		for (k = n-1; k > i; k--) {
			double tmp = A[k*n+i]/pvt;
			for (j = n-1; j > i; j--)
				A[k*n+j] -= tmp*A[i*n+j];
			for (j = 0; j < m; j++)
				B[k*m+j] -= tmp*B[i*m+j];
		}
	}
	for (i = n-1; i >= 0; i--) {
		for (j = n-1; j > i; j--)
			for (k = 0; k < m; k++)
				B[i*m+k] -= A[i*n+j]*B[j*m+k];
		for (k = 0; k < m; k++)
			B[i*m+k] /= A[i*n+i];
	}
abort:
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
bool
decmp_lu(const int n, double * __restrict A, int * __restrict idx,
         int & __restrict d)
{
	bool rv = true;
	int i, j, k;

	double * __restrict work = new double[n];
	if (work == 0) {
		event_msg(EVENT_ERROR,"Out of memory in decmp_lu()!");
		rv = false;
		goto abort;
	}

	d = 1;
	// compute the LU decomposition of a row permutation of matrix a;
	// the permutation itself is saved in idx[]
	for (i = 0; i < n; i++) {
		double tmp, max = 0.0;
		for (j = 0; j < n; j++)
			if ((tmp = fabs(A[i*n+j])) > max)
				max = tmp;
		if (max < DBL_EPSILON) {
			event_msg(EVENT_WARN,"Singular matrix A in decmp_lu()!");
			rv = false;
			goto abort;
		}
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
			event_msg(EVENT_WARN,"Singular matrix A in decmp_lu()!");
			rv = false;
			goto abort;
		}
		for (j = i+1; j < n; j++)
			A[j*n+i] /= A[i*n+i];
	}
abort:
	delete [] work;
	return rv;
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
bksub_lu(const int n, const double * __restrict A, const int * __restrict idx,
         double * __restrict b, const int m = 1, const int c = 0)
{
	int i, j, k;

	for (i = 0, k = 0; i < n; i++) {
		j = idx[i];
		double sum = b[j*m+c];
		b[j*m+c] = b[i*m+c];
		if (k != 0)
			for (j = k-1; j < i; j++)
				sum -= A[i*n+j]*b[j*m+c];
		else if (sum != 0.0)
			k = i+1;
		b[i*m+c] = sum;
	}
	for (i = n-1; i >= 0; i--) {
		double sum = b[i*m+c];
		for (j = i+1; j < n; j++)
			sum -= A[i*n+j]*b[j*m+c];
		b[i*m+c] = sum/A[i*n+i];
	}
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LU decomposition followed by forward/back
 * substitution, with a single step refinement.
 */
bool
equ_lu(const int n, const double * __restrict A, const double * __restrict b,
       double * __restrict x, const double tol)
{
	bool rv = true;
	int i, j, d;

	int * __restrict idx = new int[n];
	double * __restrict a = new double[n*n];
	double * __restrict r = new double[n];
	if (idx == 0 || a == 0 || r == 0) {
		event_msg(EVENT_ERROR,"Out of memory in equ_lu()!");
		rv = false;
		goto abort;
	}
	// avoid destroying A, B by copying them to a, x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	if (!decmp_lu(n,a,idx,d)) {
		memset(x,0,sizeof(double)*n);
		rv = false;
		goto abort;
	}
	bksub_lu(n,a,idx,x);
	for (i = 0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
	}
	bksub_lu(n,a,idx,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];

  abort:
	delete [] r;
	delete [] a;
	delete [] idx;
	return rv;
}

/*
 * This function returns the nxm matrix X = A^-1*B.  Both A and B are destoryed...
 * The result is returned in B, and the determinant in the return value.
 */
double
inv_mul_lu(const int n, const int m, double * __restrict A,
           double * __restrict B)
{
	double det = 0.0;
	int i, d;

	int * __restrict idx = new int[n];
	if (idx == 0) {
		event_msg(EVENT_ERROR,"Out of memory in inv_mul_lu()!");
		goto abort;
	}
	if (!decmp_lu(n,A,idx,d)) {
		memset(B,0,sizeof(double)*n*m);
		goto abort;
	}
	for (i = 0, det = d; i < n; i++)
		det *= A[i*n+i];
	for (i = 0; i < m; i++)
		bksub_lu(n,A,idx,B,m,i);
abort:
	delete [] idx;
	return det;
}

/*
 * Matrix inverse of the real nxn matrix A using LU decomposition.  
 */
bool
inv_lu(const int n, double * __restrict A)
{
	bool rv = true;
	int i, d;

	int * __restrict idx = new int[n];
	double * __restrict a = new double[n*n];
	if (idx == 0 || a == 0) {
		event_msg(EVENT_ERROR,"Out of memory in inv_lu()!");
		rv = false;
		goto abort;
	}
	memcpy(a,A,sizeof(double)*n*n);
	memset(A,0,sizeof(double)*n*n);
	for (i = 0; i < n; i++)
		A[i*n+i] = 1.0;
	if (!decmp_lu(n,a,idx,d)) {
		rv = false;
		goto abort;
	}
	for (i = 0; i < n; i++)
		bksub_lu(n,a,idx,A,n,i);
abort:
	delete [] a;
	delete [] idx;
	return rv;
}

/*
 * Cholesky Decomposition of postive definite nxn matrix A.
 *
 * Returns L in the lower triangle of A.
 */
bool
decmp_chol(const int n, double * __restrict A)
{
	int i, j, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			sum = A[i*n+j];
			for (k = i-1; k >= 0; k--)
                sum -= A[i*n+k]*A[j*n+k];
			if (i == j) {
                if (sum <= 0.0) {
					event_msg(EVENT_WARN,"Non-positive definite matrix in decmp_chol(%f)!",sum);
					return false;
				}
                A[i*n+i] = sqrt(sum);
			} else
                A[j*n+i] = sum/A[i*n+i];
		}
	}
	return true;
}

/*
 * Cholesky Decomposition of postive definite nxn matrix A.  The matrix is banded,
 * with bandwidth 2*w+1, and only the upper triangle is stored.
 *
 * Returns U in the upper triangle of A.
 */
bool
decmp_chol(const int n, const int w, double * __restrict A)
{
	int i, j, k;
	double sum, c, y, t;
	
    for (i = 0; i < n; i++) {
        for (j = i; j <= i+w && j < n; j++) {
			sum = A[B_IDX(n,w,i,j)];
			for (k = i-1, c = 0.0; k >= j-w && k >= 0; k--) {
				y = -fma(A[B_IDX(n,w,k,i)],A[B_IDX(n,w,k,j)],c);
				t = sum + y; c = (t - sum) - y; sum = t;
			}
			if (i == j) {
				if (sum <= 0) {
					event_msg(EVENT_WARN,"Non-positive definite matrix in decmp_chol(%f)!",double(sum));
					return false;
				}
				A[B_IDX(n,w,i,i)] = sqrt(sum);
			} else
				A[B_IDX(n,w,i,j)] = sum/A[B_IDX(n,w,i,i)];
		}
	}
	return true;
}

/*
 * This function performs forward/back substitution of a Cholesky decomposed
 * matrix, to solve a particular system.
 *
 * m is the number of columns in b, and c is the column index.
 */
void
bksub_chol(const int n, const double * __restrict A,
           double * __restrict b, const int m = 1, const int c = 0)
{
	int i, k;

	for (i = 0; i < n; i++) {
		for (k = i-1; k >= 0; k--)
			b[i*m+c] -= A[i*n+k]*b[k*m+c];
		b[i*m+c] /= A[i*n+i];
	}
	for (i = n-1; i >= 0; i--) {
		for (k = i+1; k < n; k++)
			b[i*m+c] -= A[k*n+i]*b[k*m+c];
		b[i*m+c] /= A[i*n+i];
	}
}

/*
 * This function performs forward/back substitution of a Cholesky decomposed
 * matrix, to solve a particular system.
 */
void
bksub_chol(const int n, const int w, const double * __restrict A,
           double * __restrict b, const int m = 1, const int c = 0)
{
	int i, k;
	double s, y, t;
	
    for (i = 0; i < n; i++) {
		for (k = i-1, s = 0.0; k >= i-w && k >= 0; k--) {
			y = -fma(A[B_IDX(n,w,k,i)],b[k*m+c],s); t = b[i*m+c] + y;
			s = t - b[i*m+c] - y; b[i*m+c] = t;
		}
		b[i*m+c] /= A[B_IDX(n,w,i,i)];
	}
	for (i = n-1; i >= 0; i--) {
		for (k = i+1, s = 0.0; k <= i+w && k < n; k++) {
			y = -fma(A[B_IDX(n,w,i,k)],b[k*m+c],s); t = b[i*m+c] + y;
			s = t - b[i*m+c] - y; b[i*m+c] = t;
		}
		b[i*m+c] /= A[B_IDX(n,w,i,i)];
	}
}

/*
 * Solve of Ax = b by Cholesky decomposition followed by forward/back
 * substitution.
 */
bool
equ_chol(const int n, const double * __restrict A,
         const double * __restrict b, double * __restrict x, const double tol)
{
	bool rv = true;
	//int i, j;

	double * __restrict a = new double[n*n];
	//double * r = new double[n];
	if (a == 0) { // || r == 0
		event_msg(EVENT_ERROR,"Out of memory in equ_chol()!");
		rv = false;
		goto abort;
	}
	// avoid destroying A and B by copying them to a and x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	if (!decmp_chol(n,a)) {
		memset(x,0,sizeof(double)*n);
		rv = false;
		goto abort;
	}
	bksub_chol(n,a,x);
	/*for (i = 0, dot = 0.0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
		dot += r[i]*r[i];
	}
	printf(" %g",dot);
	if (dot > ERR_TOL)
		return rv;
	bksub_chol(n,a,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];*/
abort:
	//delete [] r;
	delete [] a;
	return rv;
}

/*
 * This function returns the solution of Ax = b, where A is banded.
 *
 * The function employs Cholesky decomposition followed by forward/back
 * substitution, with a single step refinement.
 */
bool
equ_chol(const int n, const int w, const double * __restrict A,
         const double * __restrict b, double * __restrict x, const double tol)
{
	bool rv = true;
	int i, j, iter = 0;
	double dot, c1, y1, t1, c, y, t;

	double * __restrict a = new double[B_SIZE(n,w)];
	double * __restrict r = new double[n];
	if (a == 0 || r == 0) {
		event_msg(EVENT_ERROR,"Out of memory in equ_chol()!");
		rv = false;
		goto abort;
	}
	// avoid destroying A and B by copying them to a and x resp.
	memcpy(a,A,sizeof(double)*B_SIZE(n,w));
	memcpy(x,b,sizeof(double)*n);
	if (!decmp_chol(n,w,a)) {
		memset(x,0,sizeof(double)*n);
		rv = false;
		goto abort;
	}
	bksub_chol(n,w,a,x);
	while (1) {
		dot = 0.0; c1 = 0.0;
		for (i = 0; i < n; i++) {
			r[i] = -b[i];
			for (j = MAX(i-w,0), c = 0.0; j < i; j++) {
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
		if (++iter > ITER_MAX || sqrt(dot) <= tol)
			break;
		bksub_chol(n,w,a,r);
		for (i = 0; i < n; i++)
			x[i] -= r[i];
	}
abort:
	delete [] r;
	delete [] a;
	return rv;
}

/*
 * Matrix inverse of the real positive definite nxn matrix A using
 * Cholesky decomposition.  
 */
bool
inv_chol(int n, double * __restrict A)
{
	int i, j, k;
	double sum;

	if (!decmp_chol(n,A))
		return false;
	for (i = 0; i < n; i++) {
		for (j = i+1; j < n; j++) {
			for (k = i, sum = 0.0; k < j; k++)
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
	return true;
}

/*
 * LDL^T (Square root free Cholesky) Decomposition of postive definite
 * nxn matrix A.
 *
 * Returns L in the upper triangle of A (all diagonal elements are 1),
 * with D in the diagonal.
 */
bool
decmp_ldl(const int n, double * __restrict A)
{
	int i, j, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (k = i-1, sum = 0.0; k >= 0; k--) {
			sum += A[k*n+k]*(A[i*n+k]*A[i*n+k]);
		}
		A[i*n+i] -= sum;
        if (fabs(A[i*n+i]) < DBL_EPSILON) {
			event_msg(EVENT_WARN,"Indefinite matrix in decmp_ldl(%f)!",A[i*n+i]);
			return false;
		}
		for (j = i+1; j < n; j++) {
			for (k = i-1, sum = 0.0; k >= 0; k--) {
				sum += A[k*n+k]*(A[i*n+k]*A[j*n+k]);
			}
			A[j*n+i] = (A[j*n+i]-sum)/A[i*n+i];
		}
	}
	return true;
}

/*
 * This function performs forward/back substitution of a LDL^T decomposed
 * matrix, to solve a particular system.
 *
 * m is the number of columns in b, and c is the column index.
 */
void
bksub_ldl(const int n, const double * __restrict A,
          double * __restrict b, const int m = 1, const int c = 0)
{
	int i, k;
	double sum;

	for (i = 0; i < n; i++) {
		for (k = i-1, sum = 0.0; k >= 0; k--)
			sum += A[i*n+k]*b[k*m+c];
		b[i*m+c] -= sum;
	}
	for (i = n-1; i >= 0; i--) {
		b[i*m+c] /= A[i*n+i];
		for (k = i+1, sum = 0.0; k < n; k++)
			sum += A[k*n+i]*b[k*m+c];
		b[i*m+c] -= sum;
	}
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs LDL^T decomposition followed by forward/back
 * substitution, (with a single step refinement).
 */
bool
equ_ldl(const int n, const double * __restrict A,
        const double * __restrict b, double * __restrict x, const double tol)
{
	bool rv = true;
	int i, j;

	double * __restrict a = new double[n*n];
	double * __restrict r = new double[n];
	if (a == 0 || r == 0) { 
		event_msg(EVENT_ERROR,"Out of memory in equ_chol()!");
		rv = false;
		goto abort;
	}
	// avoid destroying A and B by copying them to a and x resp.
	memcpy(a,A,sizeof(double)*n*n);
	memcpy(x,b,sizeof(double)*n);
	if (!decmp_ldl(n,a)) {
		memset(x,0,sizeof(double)*n);
		rv = false;
		goto abort;
	}
	bksub_ldl(n,a,x);
	for (i = 0; i < n; i++) {
		for (j = 0, r[i] = -b[i]; j < n; j++)
			r[i] += A[i*n+j]*x[j];
	}
	bksub_ldl(n,a,r);
	for (i = 0; i < n; i++)
		x[i] -= r[i];
  abort:
	delete [] r;
	delete [] a;
	return rv;
}

/*
 * Given a mxn matrix A, this routine computes its singular value
 * decomposition, A = U * W * V'. The matrix U replaces A on output.
 * The diagonal matrix of singular values, W, is output as a vector W.
 * The matrix V (not the transpose of V) is output as V.
 */
void
decmp_svd(const int m, const int n, double * __restrict A,
          double * __restrict W, double * __restrict V)
{
	double F, G = 0.0, H;
	double C, S, X, Y;
	double scale = 0.0;
	double anorm = 0.0;
	int iter, i, j, k, q;
	bool flag;
	
	double * __restrict rv1 = new double[n];
	if (rv1 == 0) {
		event_msg(EVENT_ERROR,"Out of memory in decmp_svd()!");
		goto abort;
	}
	// Householder reduction to bidiagonal form.
	for (i = 0; i < n; i++) {
		rv1[i] = scale*G;
		if (i < m) {
			for (k = i, scale = 0.0; k < m; k++)
				scale += fabs(A[k*n+i]);
			if (scale != 0.0) {
				for (k = i, S = 0.0; k < m; k++) {
					A[k*n+i] /= scale;
					S += A[k*n+i]*A[k*n+i];
				}
				G = (A[i*n+i] > 0.0 ? -1 : 1)*sqrt(S);
				H = A[i*n+i]*G - S;
				A[i*n+i] -= G;
				for (j = i+1; j < n; j++) {
					for (k = i, S = 0.0; k < m; k++)
						S += A[k*n+i]*A[k*n+j];
					for (k = i; k < m; k++)
						A[k*n+j] += S*A[k*n+i]/H;
				}
				for (k = i; k < m; k++)
					A[k*n+i] *= scale;
				W[i] = scale*G;
			}
			for (k = i+1, scale = 0.0; k < n; k++)
				scale += fabs(A[i*n+k]);
			if (scale != 0.0) {
				for (k = i+1, S = 0.0; k < n; k++) {
					A[i*n+k] /= scale;
					S += A[i*n+k]*A[i*n+k];
				}
				G = (A[i*n+i+1] > 0.0 ? -1 : 1)*sqrt(S);
				H = A[i*n+i+1]*G - S;
				A[i*n+i+1] -= G;
				for (k = i+1; k < n; k++)
					rv1[k] = A[i*n+k] / H;
				for (j = i+1; j < m; j++ ) {
					for (k = i+1, S = 0.0; k < n; k++)
						S += A[j*n+k]*A[i*n+k];
					for (k = i+1; k < n; k++ )
						A[j*n+k] += S*rv1[k];
				}
				for (k = i+1; k < n; k++)
					A[i*n+k] *= scale;
			}
		}
		double tmp = fabs(W[i]) + fabs(rv1[i]);
		if (tmp > anorm)
			anorm = tmp;
	}
	// Accumulation of right-hand transformations.
	V[n*n-1] = 1.0;
	for (i = n-2; i >= 0; i--) {
		if (rv1[i+1] != 0.0 ) {
			for (j = i+1; j < n; j++)
				V[j*n+i] = (A[i*n+j]/A[i*n+i+1])/rv1[i+1];
			for (j = i+1; j < n; j++) {
				for (k = i+1, S = 0.0; k < n; k++)
					S += A[i*n+k]*V[k*n+j];
				for (k = i+1; k < n; k++)
					V[k*n+j] += S*V[k*n+i];
			}
		}
		for (j = i+1; j < n; j++)
			V[i*n+j] = V[j*n+i] = 0.0;
		V[i*n+i] = 1.0;
	}
	// Accumulation of left-hand transformations.
	for (i = MIN(m-1,n-1); i >= 0; i--) {
		for (j = i+1; j < n; j++)
			A[i*n+j] = 0.0;
		if (W[i] != 0.0) {
			for (j = i+1; j < n; j++) {
				for (k = i+1, S = 0.0; k < m; k++)
					S += A[k*n+i]*A[k*n+j];
				double f = (S/A[i*n+i])/W[i];
				for (k = i; k < m; k++)
					A[k*n+j] += f*A[k*n+i];
			}
			for (j = i; j < m; j++)
				A[j*n+i] /= W[i];
		} else {
			for (j = i; j < m; j++)
				A[j*n+i] = 0.0;
		}
		A[i*n+i] += 1.0;
	}
	// Diagonalization of the bidiagonal form, looping over
	// singular values.
	for (i = n-1; i >= 0; i--) {
		for (iter = 0; iter < 30; iter++) {
			// Test for splitting. Note that rv1[0] is always zero.
			flag = true;
			for (q = i; q >= 0; q--) {
				if ((fabs(rv1[q]) + anorm) == anorm) {
					flag = false;
					break;
				} else if ((fabs(W[q-1]) + anorm) == anorm)
					break;
			}
			// Cancellation of rv1[q], if q > 0.
			if (flag) {
				C = 0.0; S = 1.0;
				for (k = q; k <= i; k++) {
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
			for (j = q; q <= i-1 && j <= i-1; j++) {
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
abort:
	delete [] rv1;
}

void
bksub_svd(const int m, const int n, const double * __restrict U,
          const double * __restrict W, const double * __restrict V,
          double * __restrict b, const int p = 1, const int c = 0)
{
	int i, j;

	double * __restrict tmp = new double[n];
	if (tmp == 0) {
		event_msg(EVENT_ERROR,"Out of memory in bksub_svd()!");
		goto abort;
	}
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

abort:
	delete [] tmp;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs SVD decomposition followed by forward/back
 * substitution, with refinement.
 */
void
equ_svd(const int n, const double * __restrict A, const double * __restrict b,
        double * __restrict x, const double tol)
{
	double max;
	int i, j;

	double * __restrict U = new double[n*n];
	double * __restrict W = new double[n];
	double * __restrict V = new double[n*n];
	double * __restrict r = new double[n];
	if (U == 0 || W == 0 || V == 0 || r == 0) {
		event_msg(EVENT_ERROR,"Out of memory in equ_svd()!");
		goto abort;
	}
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
abort:
	delete [] r;
	delete [] U;
	delete [] W;
	delete [] V;
}

/*
 * Matrix inverse of the real nxn matrix A using SVD decomposition.
 */
void
inv_svd(const int n, double * __restrict A)
{
	double max;
	int i, j, k;

	double * __restrict U = new double[n*n];
	double * __restrict W = new double[n];
	double * __restrict V = new double[n*n];
	if (U == 0 || W == 0 || V == 0) {
		event_msg(EVENT_ERROR,"Out of memory in inv_svd()!");
		goto abort;
	}
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
abort:
	delete [] U;
	delete [] W;
	delete [] V;
}

/*
 * Orthonormalize the nxn matrix Q, using the SVD decomposition.
 */
void
orth_svd(const int n, double * __restrict Q)
{
	double * __restrict W = new double[n];
	double * __restrict V = new double[n*n];
	if (W == 0 || V == 0) {
		event_msg(EVENT_ERROR,"Out of memory in orth_svd()!");
		goto abort;
	}
	decmp_svd(n,n,Q,W,V);
abort:
	delete [] W;
	delete [] V;
}

/*
 * Performs a QR decomposition of the nxn matrix A.
 *
 * XXX: this routine is using a stupid access scheme.
 */
bool
decmp_qr(const int n, double * __restrict A, double * __restrict s,
         double * __restrict d)
{
	bool rv = true;
	int i, j, k;
	
	for (k = 0; k < n-1; k++) {
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
bksbp_qr(const int n, const double * __restrict A,
         const double * __restrict s, const double * __restrict d, 
         double * __restrict b, const int m = 1, const int c = 0)
{
	int i, j;
	double sum;
	for (j = 0; j < n-1; j++) {
		for (i = j, sum = 0.0; i < n; i++)
			sum += A[i*n+j]*b[i*m+c];
		sum /= s[j];
		for (i = j; i < n; i++)
			b[i*m+c] -= sum*A[i*n+j];
	}
	for (i = n-1; i >= 0; i--) {
		for (j = i+1, sum = 0.0; j < n; j++)
			sum += A[i*n+j]*b[j*m+c];
		b[i*m+c] = (b[i*m+c]-sum)/d[i];
	}
}

/*
 * Householder reduction of a nxn matrix A to tridiagonal form.
 */
void
tridiag_hh(const int n, double * __restrict A, double * __restrict d,
           double * __restrict e)
{
	int i, j, k;
	double scale, f;
	
	for (i = n-1; i >= 1; i--)  {
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
			for (k = 0, e[j] = 0.0; k <= j; k++)
				e[j] += A[j*n+k]*A[i*n+k];
			for (k = j+1; k < i; k++)
				e[j] += A[k*n+j]*A[i*n+k];
			e[j] /= d[i];
			f += e[j]*A[i*n+j];
		}
		for (j = 0; j < i; j++) {
			e[j] -= f/(2*d[i])*A[i*n+j];
			for (k = 0; k <= j; k++)
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
eig_tri_ql(const int n, double * __restrict d, double * __restrict e,
           double * __restrict A)
{
	int i, j, k, m;
	double b,c,f,g,p,r,s;

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
			for (j = m-1; m-1 >= i && j >= i; j--) {
				f = s*e[j+1], b = c*e[j+1];
				if ((e[j+2<n?j+2:0] = r = hypot(f,g)) == 0.0) {
					d[j+1] -= p;
					break;
				}
				c = g/r, s = f/r;
				d[j+1] -= p;
				r = (d[j]-d[j+1])*s+2.0*c*b;
				d[j+1] += p = s*r;
				g = c*r-b;
				for (k = 0; k < n; k++) {
					r = A[k*n+j+1];
					A[k*n+j+1] = s*A[k*n+j] + c*r;
					A[k*n+j]   = c*A[k*n+j] - s*r;
				}
			}
			e[m+1<n?m+1:0] = 0.0;
			if (j < i)
				d[i] -= p, e[i+1] = g;
		}
	}
}

/*
 * Obtain the eigenvalues and eigenvectors of a real symmetric matrix by
 * QL transform.
 */
void
eig_ql(const int n, double * __restrict A, double * __restrict d, bool sorted)
{
	double t;
	int i, j, k;

	double * __restrict e = new double[n];
	if (e == 0) {
		event_msg(EVENT_ERROR,"Out of memory in inv_eig()!");
		goto abort;
	}
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
abort:
	delete [] e;
}

void
bksub_eig(const int n, const double * __restrict Q, const double * __restrict d,
          double * __restrict b, const int p = 1, const int c = 0)
{
	int i, j;

	double * __restrict tmp = new double[n];
	if (tmp == 0) {
		event_msg(EVENT_ERROR,"Out of memory in bksub_eig()!");
		goto abort;
	}
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
abort:
	delete [] tmp;
}

/*
 * This function returns the solution of Ax = b
 *
 * The function employs an Eigenvalue decomposition followed by forward/back
 * substitution, with refinement.
 */
void
equ_eig(const int n, const double * __restrict A, const double * __restrict b,
        double * __restrict x, const double tol)
{
	double max;
	int i, j;

	double * __restrict Q = new double[n*n];
	double * __restrict d = new double[n];
	double * __restrict e = new double[n];
	double * __restrict r = new double[n];
	if (Q == 0 || d == 0 || e == 0 || r == 0) {
		event_msg(EVENT_ERROR,"Out of memory in inv_eig()!");
		goto abort;
	}
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
abort:
	delete [] r;
	delete [] e;
	delete [] d;
	delete [] Q;
}

/*
 * Matrix inverse of the real symmetric nxn matrix A using eigenvalue decomposition.
 */
void
inv_eig(const int n, double * __restrict A)
{
	double max;
	int i, j, k;

	double * __restrict Q = new double[n*n];
	double * __restrict d = new double[n];
	double * __restrict e = new double[n];
	if (Q == 0 || d == 0 || e == 0) {
		event_msg(EVENT_ERROR,"Out of memory in inv_eig()!");
		goto abort;
	}
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
abort:
	delete [] e;
	delete [] d;
	delete [] Q;
}

