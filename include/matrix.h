/**************************************************************************

	MATRIX.H - Simple matrix and vector operations

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
		This header implements simple C++ matrix structures and
		mathmematics.  The idea is to create classes which are as easy to
		use as MATLAB, but without the overhead of a class based matrix
		approach.  There are a number of low level functions which can also
		be called directly.

	Design:
		The design of these classes is based on a COM type model.  There is
		a matrix class which provides the working interface for consumers,
		but does not provide storage.  The storage is provided by a number
		of container classes, which are tailored to special storage
		requirements, like banded or diagonal matrices.  Then there are
		operators, which are collected into equations which are evaluated
		lazily.  This allows the system to minimise the use of temporaries,
		and to perform optimisations which are not possible with a direct
		approach.  This also allows for reference counting and copy on write
		syntax, so that storage is shared.
		
		The idea is to provide basic classes, and then to tie into
		BLAS/LAPACK etc.

	History:
		1993       - Created by Jeremy Lea <reg@openpave.org>
		2007/09/04 - Complete redesign (again).

**************************************************************************/

#ifndef __MATRIX_H
#define __MATRIX_H

#include "mathplus.h"
#include "event.h"
#include <memory.h>

class matrix_storage_ptr;

/*
 * class matrix_storage - Storage for a matrix.
 *
 * This interface class supports all of the details necessary for
 * storing and retrieving a matrix.
 */
class matrix_storage {
protected:
	explicit matrix_storage() : count(0), flags(is_t(0)) {
	}
	virtual ~matrix_storage() {
		if (count != 0)
			event_msg(EVENT_ERROR,"Reference counting failed!");
	}

public:
	virtual int rows() const = 0;
	virtual int cols() const = 0;
	virtual double operator() (const int i, const int j) const = 0;

	bool iszero() const {
		return (flags & zero);
	}
	bool iseye() const {
		return (flags & eye);
	}
	bool isscalar() const {
		return (flags & scalar);
	}

private:
	friend class matrix_storage_ptr;
	
	virtual int addref() {
		return ++count;
	}
	virtual int refcount() const {
		return count;
	}
	virtual int release() {
		return --count;
	}
	int count;
	enum is_t {
		zero   = 0x0001,
		eye    = 0x0002,
		diag   = 0x0004,
		sym    = 0x0008,
		banded = 0x0010,
		posdef = 0x0020,
		negdef = 0x0040,
		semi   = 0x0080,
		scalar = 0x0100
	} flags;
};

/*
 * class matrix_storage_ptr - Smart pointer for matrix storage.
 *
 * This class supports all of the details necessary for
 * storing and retrieving a matrix.
 */
class matrix_storage_ptr {
public:
	explicit matrix_storage_ptr(matrix_storage * p = 0)
		: ptr(p) {
		addref();
	}
	virtual ~matrix_storage_ptr() {
		release();
		ptr = 0;
	}
	virtual int addref() const {
		if (ptr == 0)
			return 0;
		else
			return ptr->addref();
	}
	virtual int refcount() const {
		return ptr->refcount();
	}
	virtual int release() {
		int count = ptr->release();
		if (count == 0) {
			delete ptr;
			ptr = 0;
		}
		return count;
	}
    // assignment
    matrix_storage_ptr & operator= (const matrix_storage_ptr & p) {
        if (this != &p) {
			matrix_storage * t = p.ptr;
            t->addref();
            release();
            ptr = t;
        }
        return *this;
    }
    // access the value to which the pointer refers
    matrix_storage & operator* () const {
        return *ptr;
    }
    matrix_storage * operator-> () const {
        return ptr;
    }
	
private:
	matrix_storage * ptr;
};

/*
 * struct matrix - A simple MxN matrix.
 */
class matrix_full : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_full(const int _m, const int _n)
	  : matrix_storage() {
		resize(_m,_n);
	}
	inline matrix_full(const int _m, const int _n, const double d)
	  : matrix_storage() {
		resize(_m,_n);
		for (int i = 0; i < M*N; i++)
			data[i] = d;
	}
	inline matrix_full(const int _m, const int _n, const double * v)
	  : matrix_storage() {
		resize(_m,_n);
		memcpy(data,v,sizeof(double)*M*N);
	}
	inline matrix_full(const matrix_full & A)
	  : matrix_storage() {
		resize(A.M,A.N);
		memcpy(data,A.data,sizeof(double)*M*N);
	}
	inline ~matrix_full() {
		if (data != 0)
			delete [] data;
		M = N = 0, data = 0;
	}

	// Assignment operator...
	matrix_full & operator = (const matrix_full & _m) {
		if (M != _m.M || N != _m.N)
			resize(_m.M,_m.N);
		memcpy(data,_m.data,sizeof(double)*M*N);
		return *this;
	}
	virtual int rows() const {
		return M;
	}
	virtual int cols() const {
		return N;
	}
	virtual double operator () (const int i, const int j) const {
		return data[i*N+j];
	}

protected:
	int M;						// The rows
	int N;						// The cols
	double * data;				// The data
	void resize(const int _m, const int _n) {
		if (data != 0)
			delete [] data;
		M = _m, N = _n;
		data = new double[M*N];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct matrix!");
	}
};

class matrix_operator : public matrix_storage {
public:
	inline matrix_operator()
	  : matrix_storage(), op(null), op1(0), op2(0) {
	}
	inline ~matrix_operator() {
	}

	// Assignment operator...
	matrix_operator & operator = (const matrix_operator & _m) {
		return *this;
	}
	inline int rows() const {
		switch (op) {
		case ref:
		case unref:
		case neg:
			return op1->rows();
		case trans:
			return op1->cols();
		case add:
		case sub:
			if (op1->isscalar())
				return op2->rows();
			else
				return op1->rows();
		case mul:
			return op1->rows();
		case inv:
			return op1->rows();						
		case null:
		default:
			return 0;
		}
	}
	inline int cols() const {
		switch (op) {
		case ref:
		case unref:
		case neg:
			return op1->cols();
		case trans:
			return op1->rows();
		case add:
		case sub:
			if (op1->isscalar())
				return op2->cols();
			else
				return op1->cols();
		case mul:
			return op2->cols();
		case inv:
			return op1->cols();						
		case null:
		default:
			return 0;
		}
	}
	virtual double operator () (const int i, const int j) const {
		return 0;
	}

protected:
	enum op_t {
		null, ref, unref, neg, trans, add, sub, mul, inv
	} op;
	matrix_storage_ptr op1;
	matrix_storage_ptr op2; 
};

/*
 * class matrix - Matrix interface class.
 */
class matrix {
public:
	explicit matrix()
		: data(0) {
	}
	explicit matrix(const int n, const int m)
		: data(0) {
	}

private:
	matrix_storage_ptr data;
};

/*
 * The minimum condition number for SVD and eigenvalue decompositions.
 * NOTE: This is recorded here to show the default.  It cannot be changed
 * without recompiling.
 */
#define EIG_TOL 10e-18
/*
 * The miminimum error in the result of an equals operation to trigger
 * a refinement, and the maximum number of refinements.
 */
#define ERR_TOL 10e-6
#define ITER_MAX 0

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
bool equ_gauss(const int n, const double * A, const double * b, double * x);
bool equ_lu(const int n, const double * A, const double * b, double * x, const double tol = ERR_TOL);
bool inv_lu(const int n, double * A);
bool equ_chol(const int n, const double * A, const double * b, double * x, const double tol = ERR_TOL);
bool equ_chol(const int n, const int w, const double * A, const double * b, double * x, const double tol = ERR_TOL);
bool inv_chol(const int n, double * A);
bool equ_ldl(const int n, const double * A, const double * b, double * x, const double tol = ERR_TOL);
void equ_svd(const int n, const double * A, const double * b, double * x, const double tol = ERR_TOL);
void inv_svd(const int n, double * A);
void orth_svd(const int n, double * Q);
void eig_ql(const int n, double * A, double * d);
void equ_eig(const int n, const double * A, const double * b, double * x, const double tol = ERR_TOL);
void inv_eig(const int n, double * A);

#endif // __MATRIX_H
