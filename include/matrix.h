/**************************************************************************

	MATRIX.H - Simple matrix and vector operations

	$OpenPave$

	Purpose:
		This header implements simple C++ vector and matrix structures.

	Design:
		The idea is to behave like a vector or a matrix...  Funny that.

	History:
		1993       - Created by Jeremy Lea <jlea@csir.co.za>
		2002/01/30 - Updated for use in mePADS.
		           - Removed exception handling, and made template
				     class.
		2006/07/16 - Removed templates.  Good idea, but useless...

**************************************************************************/

#ifndef __MATRIX_H
#define __MATRIX_H

#include "../include/mathplus.h"
#include "../include/event.h"
#include <memory.h>

/*
 * class vector - A simple vector.
 *
 * A vector of length N.  Has no concept of row or column vectors.
 */
struct vector {
	// Various constuctors.
	inline vector() {
		N = 0, data = 0;
	};
	inline vector(const int n) {
		resize(n);
	};
	inline vector(const int n, const double d) {
		resize(n);
		for (int i = 0; i < N; i++)
			data[i] = d;
	};
	inline vector(const int n, const double * v) {
		resize(n);
		memcpy(data,v,sizeof(double)*N);
	};
	inline vector(const vector & v) {
		resize(v.N);
		memcpy(data,v.data,sizeof(double)*N);
	};

	// Allow copying.
	inline vector & operator = (const vector & v) {
		if (N != v.N)
			resize(v.N);
		memcpy(data,v.data,sizeof(double)*N);
		return *this;
	};
	inline int n() const {
		return N;
	};
	inline double & operator [] (const int p) const {
		return data[p];
	};
	inline void zero() {
		memset(data,0,sizeof(double)*N);
	}
	inline ~vector() {
		if (data != 0)
			delete [] data;
		N = 0, data = 0;
	};

protected:
	int	N;							// The size.
	double * data;					// The data.
	inline void resize(const int n) {
		if (data != 0)
			delete [] data;
		N = n;
		data = new double[N];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct vector!");
	};
};

/*
 * These are some basic math operations for vectors, which try to be
 * as quick as possible.  They are here so they can be inlined.
 */
// Negation.
inline void neg(vector & v) {
	for (int i = 0; i < v.n(); i++)
		v[i] *= -1;
};
// Inversion.
inline void inv(vector & v) {
	for (int i = 0; i < v.n(); i++)
		v[i] = 1.0/v[i];
};
// Addition.
inline void add(vector & v, const double d) {
	for (int i = 0; d != 0.0 && i < v.n(); i++)
		v[i] += d;
};
// Subraction.
inline void sub(vector & v, const double d) {
	for (int i = 0; d != 0.0 && i < v.n(); i++)
		v[i] -= d;
};
// Multiplication.
inline void mul(vector & v, const double d) {
	for (int i = 0; i < v.n(); i++)
		v[i] *= d;
};
// Division.
inline void div(vector & v, const double d) {
	for (int i = 0; i < v.n(); i++)
		v[i] /= d;
};
// Add two vectors.
inline void add(vector & v, const vector & w) {
	for (int i = 0; i < v.n(); i++)
		v[i] += w[i];
};
// Subtraction.
inline void sub(vector & v, const vector & w) {
	for (int i = 0; i < v.n(); i++)
		v[i] -= w[i];
};
// Multiplication.   Maths buffs are turning in their graves...
inline void mul(vector & v, const vector & w) {
	for (int i = 0; i < v.n(); i++)
		v[i] *= w[i];
};
// Division.
inline void div(vector & v, const vector & w) {
	for (int i = 0; i < v.n(); i++)
		v[i] /= w[i];
};
// Dot product.
inline double dot(const vector & a, const vector & b) {
	double d = 0.0;
	for (int i = 0; i < a.n(); i++)
		d += a[i]*b[i];
	return d;
};

/*
 * struct matrix - A simple MxN matrix.
 */
struct matrix {
public:
	// A few constructors, for various uses...
	inline matrix() {
		M = N = 0, data = 0;
	};
	inline matrix(const int m, const int n) {
		resize(m,n);
	};
	inline matrix(const int m, const int n, const double d) {
		resize(m,n);
		for (int i = 0; i < M*N; i++)
			data[i] = d;
	};
	inline matrix(const int m, const int n, const double * v) {
		resize(m,n);
		memcpy(data,v,sizeof(double)*M*N);
	};
	inline matrix(const matrix & A) {
		resize(A.M,A.N);
		memcpy(data,A.data,sizeof(double)*M*N);
	};
	inline ~matrix() {
		if (data != 0)
			delete [] data;
		M = N = 0, data = 0;
	};

	// Assignment operator...
	matrix & operator = (const matrix & m) {
		if (M != m.M || N != m.N)
			resize(m.M,m.N);
		memcpy(data,m.data,sizeof(double)*M*N);
		return *this;
	};
	inline int m() const {
		return M;
	};
	inline int n() const {
		return N;
	};
	inline double * operator [] (const int r) const {
		return &data[r*N];
	};

protected:
	int M;						// The rows
	int N;						// The cols
	double * data;				// The data
	void resize(const int m, const int n) {
		if (data != 0)
			delete [] data;
		M = m, N = n;
		data = new double[M*N];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct matrix!");
	};
};

/*
 * These are some basic math operations for matrices, which try to be
 * as quick as possible.  They are here so they can be inlined.
 */
// Negation.
inline void neg(matrix & A) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] *= -1;
};
// Inversion.
inline void inv(matrix & A) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] = 1.0/A[i][j];
};
// Addition.
inline void add(matrix & A, const double d) {
	for (int i = 0; d != 0.0 && i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] += d;
};
// Subraction.
inline void sub(matrix & A, const double d) {
	for (int i = 0; d != 0.0 && i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] -= d;
};
// Multiplication.
inline void mul(matrix & A, const double d) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] *= d;
};
// Division.
inline void div(matrix & A, const double d) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] /= d;
};
// Add two matrices.
inline void add(matrix & A, const matrix & B) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] += B[i][j];
};
// Subtraction...
inline void sub(matrix & A, const matrix & B) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] -= B[i][j];
};
// Multiplication...
inline void mul(matrix & A, const matrix & B) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] *= B[i][j];
};
// Division...
inline void div(matrix & A, const matrix & B) {
	for (int i = 0; i < A.m(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] /= B[i][j];
};

/*
 * struct sqrmat - A square NxN matrix.
 */
struct sqrmat {
public:
	// A few constructors, for various uses...
	inline sqrmat() {
		N = 0, data = 0;
	};
	inline sqrmat(const int n) {
		resize(n);
	};
	inline sqrmat(const int n, const double d, const bool eye = false) {
		resize(n);
		if (eye) {
			memset(data,0,sizeof(double)*N*N);
			for (int i = 0; i < N; i++)
				data[i*N+i] = d;
		} else {
			for (int i = 0; i < N*N; i++)
				data[i] = d;
		}
	};
	inline sqrmat(const int n, const double * v) {
		resize(n);
		memcpy(data,v,sizeof(double)*N*N);
	};
	inline sqrmat(const sqrmat & A) {
		resize(A.N);
		memcpy(data,A.data,sizeof(double)*N*N);
	};
	inline ~sqrmat() {
		if (data != 0)
			delete [] data;
		N = 0, data = 0;
	};

	// Assignment operator...
	sqrmat & operator = (const sqrmat & m) {
		if (N != m.N)
			resize(m.N);
		memcpy(data,m.data,sizeof(double)*N*N);
		return *this;
	};
	inline int n() const {
		return N;
	};
	inline double * operator [] (const int r) const {
		return &data[r*N];
	};

protected:
	int N;						// The size
	double * data;				// The data
	void resize(const int n) {
		if (data != 0)
			delete [] data;
		N = n;
		data = new double[N*N];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct sqrmat!");
	};
};

/*
 * These are some basic math operations for matrices, which try to be
 * as quick as possible.  They are here so they can be inlined.
 */
// Negation.
inline void neg(sqrmat & A) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] *= -1;
};
// Inversion.
inline void inv(sqrmat & A) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] = 1.0/A[i][j];
};
// Addition.
inline void add(sqrmat & A, const double d) {
	for (int i = 0; d != 0.0 && i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] += d;
};
// Subraction.
inline void sub(sqrmat & A, const double d) {
	for (int i = 0; d != 0.0 && i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] -= d;
};
// Multiplication.
inline void mul(sqrmat & A, const double d) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] *= d;
};
// Division.
inline void div(sqrmat & A, const double d) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] /= d;
};
// Add two matrices.
inline void add(sqrmat & A, const sqrmat & B) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] += B[i][j];
};
// Subtraction...
inline void sub(sqrmat & A, const sqrmat & B) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] -= B[i][j];
};
// Multiplication...
inline void mul(sqrmat & A, const sqrmat & B) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] *= B[i][j];
};
// Division...
inline void div(sqrmat & A, const sqrmat & B) {
	for (int i = 0; i < A.n(); i++)
		for (int j = 0; j < A.n(); j++)
			A[i][j] /= B[i][j];
};

/*
 * struct dlgmat - A diagonal NxN matrix.
 */
struct dlgmat {
public:
	// A few constructors, for various uses...
	inline dlgmat() {
		N = 0, data = 0;
	};
	inline dlgmat(const int n) {
		resize(n);
	};
	inline dlgmat(const int n, const double d) {
		resize(n);
		for (int i = 0; i < N; i++)
				data[i] = d;
	};
	inline dlgmat(const int n, const double * v) {
		resize(n);
		memcpy(data,v,sizeof(double)*N);
	};
	inline dlgmat(const dlgmat & A) {
		resize(A.N);
		memcpy(data,A.data,sizeof(double)*N);
	};
	inline ~dlgmat() {
		if (data != 0)
			delete [] data;
		N = 0, data = 0;
	};

	// Assignment operator...
	dlgmat & operator = (const dlgmat & m) {
		if (N != m.N)
			resize(m.N);
		memcpy(data,m.data,sizeof(double)*N);
		return *this;
	};
	inline int n() const {
		return N;
	};
	inline double & operator [] (const int r) const {
		return data[r];
	};

protected:
	int N;						// The size
	double * data;				// The data
	void resize(const int n) {
		if (data != 0)
			delete [] data;
		N = n;
		data = new double[N];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct dlgmat!");
	};
};

/*
 * These are some basic math operations for matrices, which try to be
 * as quick as possible.  They are here so they can be inlined.
 */
// Negation.
inline void neg(dlgmat & A) {
	for (int i = 0; i < A.n(); i++)
		A[i] *= -1;
};
// Inversion.
inline void inv(dlgmat & A) {
	for (int i = 0; i < A.n(); i++)
		A[i] = 1.0/A[i];
};
// Addition.
inline void add(dlgmat & A, const double d) {
	for (int i = 0; d != 0.0 && i < A.n(); i++)
		A[i] += d;
};
// Subraction.
inline void sub(dlgmat & A, const double d) {
	for (int i = 0; d != 0.0 && i < A.n(); i++)
		A[i] -= d;
};
// Multiplication.
inline void mul(dlgmat & A, const double d) {
	for (int i = 0; i < A.n(); i++)
		A[i] *= d;
};
// Division.
inline void div(dlgmat & A, const double d) {
	for (int i = 0; i < A.n(); i++)
		A[i] /= d;
};
// Add two matrices.
inline void add(dlgmat & A, const dlgmat & B) {
	for (int i = 0; i < A.n(); i++)
		A[i] += B[i];
};
// Subtraction...
inline void sub(dlgmat & A, const dlgmat & B) {
	for (int i = 0; i < A.n(); i++)
		A[i] -= B[i];
};
// Multiplication...
inline void mul(dlgmat & A, const dlgmat & B) {
	for (int i = 0; i < A.n(); i++)
		A[i] *= B[i];
};
// Division...
inline void div(dlgmat & A, const dlgmat & B) {
	for (int i = 0; i < A.n(); i++)
		A[i] /= B[i];
};

/*
 * Returns the size of the special banded matrix storage array
 */
#define B_SIZE(n,m)		((m+1)*n-m*(m+1)/2)
/*
 * Returns the index into the special banded matrix storage array
 * This is in column major format, so loops over i are efficent.
 */
#define B_IDX(n,m,i,j)		(j <= m ? j*(j+1)/2+i : (j+1)*m+i-m*(m+1)/2)

void orth_gs(const int n, double * Q);
bool equ_lu(const int n, const double * A, const double * b, double * x);
bool inv_lu(const int n, double * A);
bool equ_chol(const int n, const double * A, const double * b, double * x);
bool inv_chol(const int n, double * A);
bool equ_udu(const int n, const double * A, const double * b, double * x);
void equ_svd(const int n, const double * A, const double * b, double * x);
void inv_svd(const int n, double * A);
void orth_svd(const int n, double * Q);
void eig_ql(const int n, double * A, double * d);
void equ_eig(const int n, const double * A, const double * b, double * x);
void inv_eig(const int n, double * A);

#endif // __MATRIX_H

