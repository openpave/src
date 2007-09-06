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
#include <assert.h>

class matrix_storage_ptr;
class matrix;

/*
 * class matrix_storage - Storage for a matrix.
 *
 * This interface class supports all of the details necessary for
 * storing and retrieving a matrix.
 */
class matrix_storage {
public:
	enum is_t {
		zero   = 0x0001,
		eye    = 0x0002,
		diag   = 0x0004,
		sym    = 0x0008,
		banded = 0x0010,
		posdef = 0x0020,
		negdef = 0x0040,
		semi   = 0x0080,
		scalar = 0x0100,
		op     = 0x0200,
		error  = 0x8000
	};
	inline is_t getflags() const {
		return flags;
	};
	inline is_t setflags(const is_t f) {
		return flags = f;
	};
	virtual int rows() const = 0;
	virtual int cols() const = 0;
	virtual double operator() (const int i, const int j) const = 0;

protected:
	explicit matrix_storage(const is_t f)
	  : count(0), flags(f) {
		flags = f;
	}
	virtual ~matrix_storage() {
		assert(count == 0);
	}

private:
	friend class matrix_storage_ptr;
	
	inline int addref() {
		return ++count;
	}
	inline int refcount() const {
		return count;
	}
	inline int release() {
		return --count;
	}
	int count;
	is_t flags;
};

/*
 * class matrix_storage_ptr - Smart pointer for matrix storage.
 *
 * This class supports all of the details necessary for
 * storing and retrieving a matrix.
 */
class matrix_storage_ptr {
public:
	inline matrix_storage_ptr(matrix_storage * p = 0)
	  : ptr(p) {
		addref();
	}
	inline matrix_storage_ptr(const matrix_storage_ptr & p)
	  : ptr(p.ptr) {
		addref();
	}
	inline ~matrix_storage_ptr() {
		release();
		ptr = 0;
	}
	inline int addref() const {
		if (ptr == 0)
			return 0;
		return ptr->addref();
	}
	inline int refcount() const {
		if (ptr == 0)
			return 0;
		return ptr->refcount();
	}
	inline int release() {
		if (ptr == 0)
			return 0;
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
			if (t != 0)
				t->addref();
			release();
			ptr = t;
		}
		return *this;
	}
	// access the value to which the pointer refers
	matrix_storage & operator* () const {
		if (ptr == 0)
			event_msg(EVENT_ERROR,"Null pointer deference in class matrix!");
		return *ptr;
	}
	matrix_storage * operator-> () const {
		if (ptr == 0)
			event_msg(EVENT_ERROR,"Null pointer deference in class matrix!");
		return ptr;
	}
	
private:
	matrix_storage * ptr;
};

/*
 * class matrix_zero - A MxN zero matrix.
 */
class matrix_zero : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_zero(const int m_, const int n_)
	  : matrix_storage(zero) {
		m = m_, n = n_;
	}
	virtual ~matrix_zero() {
	}
	virtual int rows() const {
		return m;
	}
	virtual int cols() const {
		return n;
	}
	virtual double operator () (const int i, const int j) const {
		return 0.0;
	}

protected:
	int m;						// The rows
	int n;						// The cols
	void resize(const int m_, const int n_) {
		m = m_, n = n_;
	}
};

/*
 * class matrix_eye - A MxN zero matrix.
 */
class matrix_eye : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_eye(const int n_)
	  : matrix_storage(eye) {
		n = n_;
	}
	virtual ~matrix_eye() {
	}
	virtual int rows() const {
		return n;
	}
	virtual int cols() const {
		return n;
	}
	virtual double operator () (const int i, const int j) const {
		return (i == j ? 1.0 : 0.0);
	}

protected:
	int n;						// The rows and cols
	void resize(const int n_) {
		n = n_;
	}
};

/*
 * class matrix_dense - A dense MxN matrix.
 */
class matrix_dense : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_dense(const int m_, const int n_)
	  : matrix_storage(is_t(0)) {
		resize(m_,n_);
	}
	inline matrix_dense(const int m_, const int n_, const double d)
	  : matrix_storage(is_t(0)) {
		resize(m_,n_);
		for (int i = 0; i < m*n; i++)
			data[i] = d;
	}
	inline matrix_dense(const int m_, const int n_, const double * v)
	  : matrix_storage(is_t(0)) {
		resize(m_,n_);
		memcpy(data,v,sizeof(double)*m*n);
	}
	inline matrix_dense(const matrix_dense & A)
	  : matrix_storage(A.getflags()) {
		resize(A.m,A.n);
		memcpy(data,A.data,sizeof(double)*m*n);
	}
	virtual ~matrix_dense() {
		if (data != 0)
			delete [] data;
		m = n = 0, data = 0;
	}

	// Assignment operator...
	matrix_dense & operator = (const matrix_dense & m_) {
		if (m != m_.m || n != m_.n)
			resize(m_.m,m_.n);
		memcpy(data,m_.data,sizeof(double)*m*n);
		return *this;
	}
	virtual int rows() const {
		return m;
	}
	virtual int cols() const {
		return n;
	}
	virtual double operator () (const int i, const int j) const {
		return data[i*n+j];
	}

protected:
	int m;						// The rows
	int n;						// The cols
	double * data;				// The data
	void resize(const int m_, const int n_) {
		if (data != 0)
			delete [] data;
		m = m_, n = n_;
		data = new double[m*n];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct matrix!");
	}
};

class matrix_operator : public matrix_storage {
public:
	enum op_t {
		ref, unref, neg, trans, add, sub, mul, inv
	};

public:
	inline matrix_operator(const op_t op_,
		matrix_storage * op1_ = 0, matrix_storage * op2_ = 0)
	  : matrix_storage(is_t(0)), op(op_), op1(op1_), op2(op2_) {
	}
	inline matrix_operator(const op_t op_,
		const matrix_storage_ptr & op1_)
	  : matrix_storage(is_t(0)), op(op_), op1(op1_), op2(0) {
	}
	inline matrix_operator(const op_t op_,
		const matrix_storage_ptr & op1_,
		const matrix_storage_ptr & op2_)
	  : matrix_storage(is_t(0)), op(op_), op1(op1_), op2(op2_) {
	}
	virtual ~matrix_operator() {
	}

	// Assignment operator...
	matrix_operator & operator = (const matrix_operator & m_) {
		return *this;
	}
	virtual int rows() const {
		switch (op) {
		case ref:
		case unref:
		case neg:
			return op1->rows();
		case trans:
			return op1->cols();
		case add:
		case sub:
			if (op1->getflags() & scalar)
				return op2->rows();
			else
				return op1->rows();
		case mul:
			return op1->rows();
		case inv:
			return op1->rows();						
		default:
			return 0; // XXX
		}
	}
	virtual int cols() const {
		switch (op) {
		case ref:
		case unref:
		case neg:
			return op1->cols();
		case trans:
			return op1->rows();
		case add:
		case sub:
			if (op1->getflags() & scalar)
				return op2->cols();
			else
				return op1->cols();
		case mul:
			return op2->cols();
		case inv:
			return op1->cols();						
		default:
			return 0; // XXX
		}
	}
	virtual double operator () (const int i, const int j) const {
		switch (op) {
		case ref:
		case unref:
			return (*op1)(i,j);
		case neg:
			return -(*op1)(i,j);
		case trans:
			return (*op1)(j,i);
		case add:
		case sub:
		case mul:
		case inv:
			// We use the const cast here to get around the compiler.  We
			// need this function to be const, so it is called in the right
			// way.  We're doing delayed evaluation, so we're not really
			// changing the value of the object.
			const_cast<matrix_operator *>(this)->evaluate();
			return (*op1)(i,j);
		default:
			return 0; // XXX
		}
	}

protected:
	op_t op;
	matrix_storage_ptr op1;
	matrix_storage_ptr op2; 

	inline void evaluate() {
		switch (op) {
		case ref:
		case unref:
			return;
		case neg:
			return;
		case trans:
			return;
		case add:
			return;
		case sub:
			return;
		case mul:
			return;
		case inv:
			return;
		default:
			return; // XXX
		}
	}
};

/*
 * class matrix - Matrix interface class.
 */
class matrix {
public:
	inline matrix()
		: data(0) {
	}
	inline matrix(const int m, const int n)
		: data(0) {
		matrix_storage * d = new matrix_zero(m,n);
		if (d == 0)
			event_msg(EVENT_ERROR,"Out of memory for class matrix!");
		else
			data = matrix_storage_ptr(d);
	}
	inline matrix(const matrix & m)
		: data(m.data) {
	}
	inline matrix(matrix_storage * d)
		: data(matrix_storage_ptr(d)) {
	}

	inline int rows() const {
		return data->rows();
	}
	inline int cols() const {
		return data->cols();
	}
	inline double operator() (const int i, const int j) const {
		return (*data)(i,j);
	}

	inline bool iszero() const {
		return (data->getflags() & matrix_storage::zero);
	}
	inline bool iseye() const {
		return (data->getflags() & matrix_storage::eye);
	}
	inline bool isscalar() const {
		return (data->getflags() & matrix_storage::scalar);
	}

	// Assignment operator
	inline matrix & operator = (const matrix & m) {
		data = m.data;
		return *this;
	}
	// Comparison operators...
	inline bool operator == (const matrix & a) const {
		int m = rows(), n = cols(), i, j;
		if (m != a.rows() || n != a.cols())
			return false;
		for (i = 0; i < m; i++) {
			for (j = 0; j < n; j++) {
				if ((*data)(i,j) != a(i,j))
					return false;
			}
		} 
		return true;
	}

private:
	friend matrix operator! (matrix &);
	friend matrix operator- (const matrix &);

	matrix_storage_ptr data;
};

// The magic unref operator...
// Yes I know.  You're not supposed to overload operators with different
// syntax to their normal C/C++ meaning.  But this has no meaning
// for a matrix, so we can put it to a very good use.
inline matrix operator! (matrix & b) {
	matrix_storage * d = new matrix_operator(matrix_operator::unref,b.data);
	if (d == 0) {
		event_msg(EVENT_ERROR,"Out of memory for class matrix!");
		return matrix();
	}
	b.data.release();
	return matrix(d);
}

// Some math operators...
inline matrix operator- (const matrix & b) {
	matrix_storage * d = new matrix_operator(matrix_operator::neg,b.data);
	if (d == 0) {
		event_msg(EVENT_ERROR,"Out of memory for class matrix!");
		return matrix();
	} else
		return matrix(d);
}

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
