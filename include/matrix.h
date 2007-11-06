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
#include "tmatrix.h"
#include <memory.h>
#include <assert.h>

class matrix_storage_ptr;
class matrix_dense;
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
		none   = 0x0000,
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
	virtual unsigned rows() const = 0;
	virtual unsigned cols() const = 0;
	virtual double operator() (const unsigned i,
		const unsigned j) const = 0;
	virtual matrix_dense getdense() const = 0;

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
		assert(count > 0);
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
	inline matrix_storage_ptr(matrix_storage * restrict p = 0)
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
		if (count == 0)
			delete ptr;
		ptr = 0;
		return count;
	}
	// assignment
	matrix_storage_ptr & operator= (const matrix_storage_ptr & p) {
		if (this != &p) {
			matrix_storage * restrict t = p.ptr;
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
	matrix_storage * restrict ptr;
};

/*
 * class matrix_dense - A dense MxN matrix.
 */
class matrix_dense : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_dense(const unsigned m, const unsigned n)
	  : matrix_storage(none), data(0) {
		resize(m,n);
	}
	inline matrix_dense(const unsigned m, const unsigned n,
		const double d, bool mkdiag = true)
	  : matrix_storage(none), data(0) {
		resize(m,n);
		if (mkdiag) {
			memset(data,0,M*N*sizeof(double));
			for (unsigned i = 0; d != 0.0 && i < MIN(M,N); i++)
				data[i*N+i] = d;
		} else {
			if (d == 0.0)
				memset(data,0,M*N*sizeof(double));
			else {
				for (unsigned i = 0; i < M*N; i++)
					data[i] = d;
			}
		}
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				data[i*N+j] = (!mkdiag || i == j ? d : 0.0);
	}
	inline matrix_dense(const unsigned m, const unsigned n,
		const double * restrict v)
	  : matrix_storage(none), data(0) {
		resize(m,n);
		memcpy(data,v,M*N*sizeof(double));
	}
	inline matrix_dense(const matrix_dense & A)
	  : matrix_storage(A.getflags()), data(0) {
		resize(A.M,A.N);
		memcpy(data,A.data,M*N*sizeof(double));
	}
	virtual ~matrix_dense() {
		if (data != 0)
			delete [] data;
		M = N = 0, data = 0;
	}

	// Assignment operator...
	inline matrix_dense & operator= (const matrix_dense & m) {
		if (M != m.M || N != m.N)
			resize(m.M,m.N);
		memcpy(data,m.data,M*N*sizeof(double));
		return *this;
	}
	inline matrix_dense & operator+= (const double & d) {
		for (unsigned i = 0; d != 0.0 && i < M*N; i++)
			data[i] += d;
		return *this;
	}
	inline matrix_dense & operator+= (const matrix_dense & m) {
		assert(rows() == m.rows() && cols() == m.cols());
		for (unsigned i = 0; i < M*N; i++)
			data[i] += m.data[i];
		return *this;
	}
	inline matrix_dense & operator-= (const matrix_dense & m) {
		assert(rows() == m.rows() && cols() == m.cols());
		for (unsigned i = 0; i < M*N; i++)
			data[i] -= m.data[i];
		return *this;
	}
	inline matrix_dense & operator*= (const double & d) {
		if (d == 0.0)
			memset(data,0,M*N*sizeof(double));
		else {
			for (unsigned i = 0; i < M*N; i++)
				data[i] *= d;
		}
		return *this;
	}
	virtual unsigned rows() const {
		return M;
	}
	virtual unsigned cols() const {
		return N;
	}
	virtual double operator() (const unsigned i, const unsigned j) const {
		return data[i*N+j];
	}
	inline double & operator() (const unsigned i, const unsigned j) {
		return data[i*N+j];
	}
	virtual double operator() (const unsigned i) const {
		return data[i];
	}
	inline double & operator() (const unsigned i) {
		return data[i];
	}
	virtual matrix_dense getdense() const {
		return matrix_dense(*this);
	} 

protected:
	unsigned M;					// The rows
	unsigned N;					// The cols
	double * restrict data;		// The data
	void resize(const unsigned m, const unsigned n) {
		if (data != 0)
			delete [] data;
		M = m, N = n;
		data = new double[M*N];
		if (data == 0)
			event_msg(EVENT_ERROR,"Out of memory for struct matrix!");
	}
};

inline matrix_dense
operator~ (const matrix_dense & m) {
	matrix_dense t(m.cols(),m.rows());
	for (unsigned i = 0; i < t.rows(); i++)
		for (unsigned j = 0; j < t.cols(); j++)
			t(i,j) = m(j,i);
	return t;
}

inline matrix_dense
operator- (const matrix_dense & m) {
	matrix_dense t(m.rows(),m.cols());
	for (unsigned i = 0; i < t.rows(); i++)
		for (unsigned j = 0; j < t.cols(); j++)
			t(i,j) = -m(i,j);
	return t;
}

inline matrix_dense
operator+ (const matrix_dense & m, const double & d) {
	matrix_dense t(m);
	return t += d;
}

inline matrix_dense
operator+ (const double & d, const matrix_dense & m) {
	return m + d;
}

inline matrix_dense
operator- (const matrix_dense & m, const double & d) {
	return m + (-d);
}

inline matrix_dense
operator- (const double & d, const matrix_dense & m) {
	matrix_dense t(m.rows(),m.cols());
	for (unsigned i = 0; i < t.rows(); i++)
		for (unsigned j = 0; j < t.cols(); j++)
			t(i,j) = d - m(i,j);
	return t;
}

inline matrix_dense
operator* (const matrix_dense & m, const double & d) {
	matrix_dense t(m);
	return t *= d;
}

inline matrix_dense
operator* (const double & d, const matrix_dense & m) {
	return m*d;
}

inline matrix_dense
operator/ (const matrix_dense & m, const double & d) {
	return m*(1.0/d);
}

inline matrix_dense
operator/ (const double & d, const matrix_dense & m) {
	matrix_dense t(m.rows(),m.cols());
	for (unsigned i = 0; i < t.rows(); i++)
		for (unsigned j = 0; j < t.cols(); j++)
			t(i,j) = d/m(i,j);
	return t;
}

inline matrix_dense
operator+ (const matrix_dense & a, const matrix_dense & b) {
	matrix_dense t(a);
	return t += b;
}

inline matrix_dense
operator- (const matrix_dense & a, const matrix_dense & b) {
	matrix_dense t(a);
	return t -= b;
}

inline matrix_dense
operator* (const matrix_dense & a, const matrix_dense & b) {
	assert(a.cols() == b.rows());
	matrix_dense t(a.rows(),b.cols());
	for (unsigned i = 0; i < a.rows(); i++) {
		for (unsigned j = 0; j < b.cols(); j++) {
			t(i,j) = 0.0;
			for (unsigned k = 0; k < a.cols(); k++)
				t(i,j) += a(i,k)*b(k,j);
		}
	}
	return t;
}

/*
 * class matrix_zero - A MxN zero matrix.
 */
class matrix_zero : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_zero(const unsigned m, const unsigned n)
	  : matrix_storage(zero) {
		M = m, N = n;
	}
	virtual ~matrix_zero() {
	}
	virtual unsigned rows() const {
		return M;
	}
	virtual unsigned cols() const {
		return N;
	}
	virtual double operator() (const unsigned i, const unsigned j) const {
		return 0.0;
	}
	virtual matrix_dense getdense() const {
		return matrix_dense(M,N,0.0,false);
	}

protected:
	unsigned M;						// The rows
	unsigned N;						// The cols
	void resize(const unsigned m, const unsigned n) {
		M = m, N = n;
	}
};

/*
 * class matrix_eye - A MxN zero matrix.
 */
class matrix_eye : public matrix_storage {
public:
	// A few constructors, for various uses...
	inline matrix_eye(const unsigned n)
	  : matrix_storage(eye) {
		N = n;
	}
	virtual ~matrix_eye() {
	}
	virtual unsigned rows() const {
		return N;
	}
	virtual unsigned cols() const {
		return N;
	}
	virtual double operator() (const unsigned i, const unsigned j) const {
		return (i == j ? 1.0 : 0.0);
	}
	virtual matrix_dense getdense() const {
		return matrix_dense(N,N,1.0,true);
	}

protected:
	unsigned N;						// The rows and cols
	void resize(const unsigned n) {
		N = n;
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
	  : matrix_storage(none), op(op_), op1(op1_), op2(op2_) {
	}
	inline matrix_operator(const op_t op_,
		const matrix_storage_ptr & op1_)
	  : matrix_storage(none), op(op_), op1(op1_), op2(0) {
	}
	inline matrix_operator(const op_t op_,
		const matrix_storage_ptr & op1_,
		const matrix_storage_ptr & op2_)
	  : matrix_storage(none), op(op_), op1(op1_), op2(op2_) {
	}
	virtual ~matrix_operator() {
	}

	// Assignment operator...
	matrix_operator & operator= (const matrix_operator & m) {
		// XXX
		return *this;
	}
	virtual unsigned rows() const {
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
	virtual unsigned cols() const {
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
	virtual double operator() (const unsigned i, const unsigned j) const {
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
	virtual matrix_dense getdense() const {
		const_cast<matrix_operator *>(this)->evaluate();
		matrix_dense t(rows(),cols());
		for (unsigned i = 0; i < rows(); i++)
			for (unsigned j = 0; j < cols(); j++) 
				t(i,j) = (*op1)(i,j);
		return t;
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
	inline matrix(const unsigned m, const unsigned n)
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

	inline unsigned rows() const {
		return data->rows();
	}
	inline unsigned cols() const {
		return data->cols();
	}
	inline double operator() (const unsigned i, const unsigned j) const {
		return (*data)(i,j);
	}

	inline bool iszero() const {
		return ((data->getflags() & matrix_storage::zero) != 0);
	}
	inline bool iseye() const {
		return ((data->getflags() & matrix_storage::eye) != 0);
	}
	inline bool isscalar() const {
		return ((data->getflags() & matrix_storage::scalar) != 0);
	}

	// Assignment operator
	inline matrix & operator= (const matrix & m) {
		data = m.data;
		return *this;
	}
	// Comparison operators...
	inline bool operator== (const matrix & a) const {
		unsigned m = rows(), n = cols(), i, j;
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
	// Comparison operators...
	inline bool operator!= (const matrix & a) const {
		return !(*this == a);
	}

private:
	friend matrix operator! (matrix &);
	friend matrix operator~ (const matrix &);
	friend matrix operator- (const matrix &);

	matrix_storage_ptr data;
};

// The magic unref operator...
// Yes I know.  You're not supposed to overload operators with different
// syntax to their normal C/C++ meaning.  But this has no meaning
// for a matrix, so we can put it to a very good use.
inline matrix
operator! (matrix & b) {
	matrix_storage * d = new matrix_operator(matrix_operator::unref,b.data);
	if (d == 0) {
		event_msg(EVENT_ERROR,"Out of memory for class matrix!");
		return matrix();
	}
	b.data.release();
	return matrix(d);
}

// Matrix transpose...
inline matrix
operator~ (const matrix & b) {
	matrix_storage * d = new matrix_operator(matrix_operator::trans,b.data);
	if (d == 0) {
		event_msg(EVENT_ERROR,"Out of memory for class matrix!");
		return matrix();
	} else
		return matrix(d);
}

// Some math operators...
inline matrix
operator- (const matrix & b) {
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
#define B_SIZE(n,w)		((w+1)*n-w*(w+1)/2)
/*
 * Returns the index into the special banded matrix storage array
 * This is in column major format, so loops over i are efficent.
 */
#define B_IDX(n,w,i,j)		(j <= w ? j*(j+1)/2+i : (j+1)*w+i-w*(w+1)/2)

void orth_gs(const int n, double * restrict Q);
bool equ_gauss(const int n, const double * restrict A,
	const double * restrict b, double * restrict x);
double inv_mul_gauss(const int n, const int m, double * restrict A,
	double * restrict B);
bool equ_lu(const int n, const double * restrict A,
	const double * restrict b, double * restrict x,
	const double tol = ERR_TOL);
bool inv_lu(const int n, double * restrict A);
double inv_mul_lu(const int n, const int m, double * restrict A,
	double * restrict B);
bool equ_chol(const int n, const double * restrict A,
	const double * restrict b, double * restrict x,
	const double tol = ERR_TOL);
bool equ_chol(const int n, const int w, const double * restrict A,
	const double * restrict b, double * restrict x,
	const double tol = ERR_TOL);
bool inv_chol(const int n, double * restrict A);
bool equ_ldl(const int n, const double * restrict A,
	const double * restrict b, double * restrict x,
	const double tol = ERR_TOL);
void equ_svd(const int n, const double * restrict A,
	const double * restrict b, double * restrict x,
	const double tol = ERR_TOL);
void inv_svd(const int n, double * restrict A);
void orth_svd(const int n, double * restrict Q);
void eig_ql(const int n, double * restrict A,
	double * restrict d, bool sorted = true);
void equ_eig(const int n, const double * restrict A,
	const double * restrict b, double * restrict x,
	const double tol = ERR_TOL);
void inv_eig(const int n, double * restrict A);

#endif // MATRIX_H
