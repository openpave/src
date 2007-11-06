/**************************************************************************

	TMATRIX.H - Templated matrix and vector operations

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
		This header file implements some basic templated matrix operators,
		designed for use with small vectors and matrices.  These are stored
		on the stack, and use template meta programs to do their math.  You
		should not use these classes for more than a 10x10 matrix or so.

	Design:
		The design of these classes follows a number of online references
		for C++ template meta programs and expression templates.  There are
		a number of existing classes which do all of this, and better than
		here, but they do not integrate well into the OpenPave.org world.
		
		These are quite difficult to understand.  Please don't play without
		extensive testing.  Also, this class will not compile on VC6.

	History:
		2007/09/04 - Complete redesign (again).

**************************************************************************/

#ifndef __TMATRIX_H
#define __TMATRIX_H

#include <stdio.h>

/*
 * struct meta_assign - template meta class to assign one matrix to another
 */
template<unsigned I, unsigned J>
struct meta_assign {
	template<class A, class E>
	static inline void assignR(A & lhs, const E & rhs) {
	    lhs(I,J) = rhs(I,J);
	    meta_assign<I,J-1>::assignR(lhs,rhs);
	}
	template<class A, class E>
	static inline void assign(A & lhs, const E & rhs) {
	    meta_assign<I,J>::assignR(lhs,rhs);
	    meta_assign<I-1,J>::assign(lhs,rhs);
	}
	template<class A, class E>
	static inline void assignR_add(A & lhs, const E & rhs) {
	    lhs(I,J) += rhs(I,J);
	    meta_assign<I,J-1>::assignR_add(lhs,rhs);
	}
	template<class A, class E>
	static inline void assign_add(A & lhs, const E & rhs) {
	    meta_assign<I,J>::assignR_add(lhs,rhs);
	    meta_assign<I-1,J>::assign_add(lhs,rhs);
	}
	template<class A, class E>
	static inline void assignR_sub(A & lhs, const E & rhs) {
	    lhs(I,J) -= rhs(I,J);
	    meta_assign<I,J-1>::assignR_sub(lhs,rhs);
	}
	template<class A, class E>
	static inline void assign_sub(A & lhs, const E & rhs) {
	    meta_assign<I,J>::assignR_sub(lhs,rhs);
	    meta_assign<I-1,J>::assign_sub(lhs,rhs);
	}
	template<class A, typename T>
	static inline void assignR_d(A & lhs, const T & d) {
	    lhs(I,J) = d;
	    meta_assign<I,J-1>::assignR_d(lhs,d);
	}
	template<class A, typename T>
	static inline void assign_d(A & lhs, const T & d) {
	    meta_assign<I,J>::assignR_d(lhs,d);
	    meta_assign<I-1,J>::assign_d(lhs,d);
	}
	template<class A, typename T>
	static inline void assignR_add_d(A & lhs, const T & d) {
	    lhs(I,J) += d;
	    meta_assign<I,J-1>::assignR_add_d(lhs,d);
	}
	template<class A, typename T>
	static inline void assign_add_d(A & lhs, const T & d) {
	    meta_assign<I,J>::assignR_add_d(lhs,d);
	    meta_assign<I-1,J>::assign_add_d(lhs,d);
	}
	template<class A, typename T>
	static inline void assignR_mul_d(A & lhs, const T & d) {
	    lhs(I,J) *= d;
	    meta_assign<I,J-1>::assignR_mul_d(lhs,d);
	}
	template<class A, typename T>
	static inline void assign_mul_d(A & lhs, const T & d) {
	    meta_assign<I,J>::assignR_mul_d(lhs,d);
	    meta_assign<I-1,J>::assign_mul_d(lhs,d);
	}
};
template<unsigned I> struct meta_assign<I,-1> {
	template<class A, class E>
	static inline void assignR(A &, const E &) {}
	template<class A, class E>
	static inline void assignR_add(A &, const E &) {}
	template<class A, class E>
	static inline void assignR_sub(A &, const E &) {}
	template<class A, typename T>
	static inline void assignR_d(A &, const T & d) {}
	template<class A, typename T>
	static inline void assignR_add_d(A &, const T & d) {}
	template<class A, typename T>
	static inline void assignR_mul_d(A &, const T & d) {}
};
template<unsigned J> struct meta_assign<-1,J> {
	template<class A, class E>
	static inline void assign(A &, const E &) {}
	template<class A, class E>
	static inline void assign_add(A &, const E &) {}
	template<class A, class E>
	static inline void assign_sub(A &, const E &) {}
	template<class A, typename T>
	static inline void assign_d(A &, const T & d) {}
	template<class A, typename T>
	static inline void assign_add_d(A &, const T & d) {}
	template<class A, typename T>
	static inline void assign_mul_d(A &, const T & d) {}
};

/*
 * struct tmatrix_expr - Encapsulates a matrix as an expression
 */
template<typename T, class E, unsigned M, unsigned N>
struct tmatrix_expr {
	inline explicit tmatrix_expr(const E & e) : m_e(e) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j);
	}
	inline T operator() (const int i) const {
		return m_e(i/N,i%N);
	}
private:
	const E m_e;
};

/*
 * class tmatrix - A simple MxN matrix.
 */
template<typename T, unsigned M, unsigned N>
class tmatrix {
public:
	// A few constructors, for various uses...
	inline explicit tmatrix() {
	}
	inline explicit tmatrix(const T & d, const bool eye = false) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				data[i][j] = (!eye || i == j ? d : 0.0);
	}
	inline explicit tmatrix(const T * __restrict v) {
		memcpy(data,v,M*N*sizeof(T));
	}
	inline tmatrix(const tmatrix & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
	}
	template<class E>
	inline explicit tmatrix(const tmatrix_expr<T,E,M,N> & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
	}
	~tmatrix() {
	}

	inline tmatrix & operator= (const tmatrix & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
		return *this;
	}
	inline tmatrix & operator= (const T & d) {
		meta_assign<M-1,N-1>::assign_d(*this,d);
		return *this;
	}
	template<class E>
	inline tmatrix & operator= (const tmatrix_expr<T,E,M,N> & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
		return *this;
	}
	inline tmatrix & operator+= (const tmatrix & m) {
		meta_assign<M-1,N-1>::assign_add(*this,m);
		return *this;
	}
	template<class E>
	inline tmatrix & operator+= (const tmatrix_expr<T,E,M,N> & m) {
		meta_assign<M-1,N-1>::assign_add(*this,m);
		return *this;
	}
	inline tmatrix & operator+= (const T & d) {
		meta_assign<M-1,N-1>::assign_add_d(*this,d);
		return *this;
	}
	inline tmatrix & operator-= (const tmatrix & m) {
		meta_assign<M-1,N-1>::assign_sub(*this,m);
		return *this;
	}
	template<class E>
	inline tmatrix & operator-= (const tmatrix_expr<T,E,M,N> & m) {
		meta_assign<M-1,N-1>::assign_sub(*this,m);
		return *this;
	}
	inline tmatrix & operator-= (const T & d) {
		meta_assign<M-1,N-1>::assign_add_d(*this,-d);
		return *this;
	}
	inline tmatrix & operator*= (const T & d) {
		meta_assign<M-1,N-1>::assign_mul_d(*this,d);
		return *this;
	}
	inline tmatrix & operator/= (const T & d) {
		meta_assign<M-1,N-1>::assign_mul_d(*this,1.0/d);
		return *this;
	}
	inline int rows() const {
		return M;
	}
	inline int cols() const {
		return N;
	}
	inline T operator() (const int r, const int c) const {
		return data[r][c];
	}
	inline T operator() (const int i) const {
		return data[i/N][i%N];
	}
	inline T & operator() (const int r, const int c) {
		return data[r][c];
	}
	inline T & operator() (const int i) {
		return data[i/N][i%N];
	}
	inline T * operator[] (const int r) {
		return data[r];
	}

	void print() const {
		printf("[\n");
		for (unsigned i = 0; i < M; ++i) {
			printf("\t[");
			for (unsigned j = 0; j < N; ++j)
				printf("\t%+6.4f", static_cast<double>(
						this->operator()(i,j)));
			printf("]\n");
		}
		printf("]\n");
	}

private:
	T data[M][N];
};

/*
 * class tmatrix_scalar - A 1x1 matrix.
 */
template<typename T>
class tmatrix_scalar {
public:
	inline tmatrix_scalar(const tmatrix<T,1,1> & m) {
		m_d = m(0,0);
	}
	template<class E>
	inline tmatrix_scalar(const tmatrix_expr<T,E,1,1> & m) {
		m_d = m(0,0);
	}
	~tmatrix_scalar() {
	}
	inline operator T () const {
		return m_d;
	}
private:
	T m_d;
};

/*
 * Support for matrix transpose
 */
template<typename T, class E>
struct tmatrix_expr_trans {
	explicit tmatrix_expr_trans(const E & e) : m_e(e) { }
	inline T operator() (const unsigned i, const unsigned j) const {
  		return m_e(j,i);
  	}
private:
	const E m_e;
};
template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_trans<T,tmatrix<T,M,N> >,N,M>
operator~ (const tmatrix<T,M,N> & m) {  
	typedef tmatrix_expr_trans<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,N,M>(expr_t(m));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_trans<T,tmatrix_expr<T,E,M,N> >,N,M>
operator~ (const tmatrix_expr<T,E,M,N> & m) {  
	typedef tmatrix_expr_trans<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,N,M>(expr_t(m));
}

/*
 * Support for matrix negate
 */
template<typename T, class E>
struct tmatrix_expr_neg {
	explicit tmatrix_expr_neg(const E & e) : m_e(e) { }
	inline T operator() (const unsigned i, const unsigned j) const {
  		return -m_e(i,j);
  	}
private:
	const E m_e;
};
template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_neg<T,tmatrix<T,M,N> >,M,N>
operator- (const tmatrix<T,M,N> & m) {  
	typedef tmatrix_expr_neg<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(m));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_neg<T,tmatrix_expr<T,E,M,N> >,M,N>
operator- (const tmatrix_expr<T,E,M,N> & m) {  
	typedef tmatrix_expr_neg<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(m));
}

/*
 * Add a value to a matrix
 */
template<typename T, class E>
struct tmatrix_expr_add_d {
	explicit tmatrix_expr_add_d(const E & e, const T & d)
	  : m_e(e), m_d(d) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j) + m_d;
	}
private:
	const E m_e;
	const T m_d;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add_d<T,tmatrix<T,M,N> >,M,N> 
operator+ (const tmatrix<T,M,N> & e, const T & d) {
	typedef tmatrix_expr_add_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator+ (const tmatrix_expr<T,E,M,N> & e, const T & d) {
	typedef tmatrix_expr_add_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add_d<T,tmatrix<T,M,N> >,M,N> 
operator+ (const T & d, const tmatrix<T,M,N> & e) {
	typedef tmatrix_expr_add_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator+ (const T & d, const tmatrix_expr<T,E,M,N> & e) {
	typedef tmatrix_expr_add_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add_d<T,tmatrix<T,M,N> >,M,N> 
operator- (const tmatrix<T,M,N> & e, const T & d) {
	typedef tmatrix_expr_add_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,-d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator- (const tmatrix_expr<T,E,M,N> & e, const T & d) {
	typedef tmatrix_expr_add_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,-d));
}

/*
 * Subtract a matrix from a value
 */
template<typename T, class E>
struct tmatrix_expr_sub_d {
	explicit tmatrix_expr_sub_d(const E & e, const T & d)
	  : m_e(e), m_d(d) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_d - m_e(i,j);
	}
private:
	const E m_e;
	const T m_d;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_sub_d<T,tmatrix<T,M,N> >,M,N> 
operator- (const T & d, const tmatrix<T,M,N> & e) {
	typedef tmatrix_expr_sub_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_sub_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator- (const T & d, const tmatrix_expr<T,E,M,N> & e) {
	typedef tmatrix_expr_sub_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}

/*
 * Multiply a matrix by a value
 */
template<typename T, class E>
struct tmatrix_expr_mul_d {
	explicit tmatrix_expr_mul_d(const E & e, const T & d)
	  : m_e(e), m_d(d) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_d*m_e(i,j);
	}
private:
	const E m_e;
	const T m_d;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_mul_d<T,tmatrix<T,M,N> >,M,N> 
operator* (const tmatrix<T,M,N> & e, const T & d) {
	typedef tmatrix_expr_mul_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_mul_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator* (const tmatrix_expr<T,E,M,N> & e, const T & d) {
	typedef tmatrix_expr_mul_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_mul_d<T,tmatrix<T,M,N> >,M,N> 
operator* (const T & d, const tmatrix<T,M,N> & e) {
	typedef tmatrix_expr_mul_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_mul_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator* (const T & d, const tmatrix_expr<T,E,M,N> & e) {
	typedef tmatrix_expr_mul_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}

/*
 * Divide a matrix by a value
 */
template<typename T, class E>
struct tmatrix_expr_div_d {
	explicit tmatrix_expr_div_d(const E & e, const T & d)
	  : m_e(e), m_d(d) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j)/m_d;
	}
private:
	const E m_e;
	const T m_d;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_div_d<T,tmatrix<T,M,N> >,M,N> 
operator/ (const tmatrix<T,M,N> & e, const T & d) {
	typedef tmatrix_expr_div_d<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_mul_d<T,tmatrix_expr<T,E,M,N> >,M,N>
operator/ (const tmatrix_expr<T,E,M,N> & e, const T & d) {
	typedef tmatrix_expr_div_d<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}

/*
 * Divide a value by a matrix elementwise
 */
template<typename T, class E>
struct tmatrix_expr_d_div {
	explicit tmatrix_expr_d_div(const E & e, const T & d)
	  : m_e(e), m_d(d) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_d/m_e(i,j);
	}
private:
	const E m_e;
	const T m_d;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_d_div<T,tmatrix<T,M,N> >,M,N> 
operator/ (const T & d, const tmatrix<T,M,N> & e) {
	typedef tmatrix_expr_d_div<T,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_d_div<T,tmatrix_expr<T,E,M,N> >,M,N>
operator/ (const T & d, const tmatrix_expr<T,E,M,N> & e) {
	typedef tmatrix_expr_d_div<T,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e,d));
}

/*
 * Matrix addition
 */
template<typename T, class E1, class E2>
struct tmatrix_expr_add {
	explicit tmatrix_expr_add(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_e1(i,j) + m_e2(i,j);
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add<T,tmatrix<T,M,N>,tmatrix<T,M,N> >,M,N> 
operator+ (const tmatrix<T,M,N> & e1, const tmatrix<T,M,N> & e2) {
	typedef tmatrix_expr_add<T,tmatrix<T,M,N>,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add<T,tmatrix_expr<T,E,M,N>,tmatrix<T,M,N> >,M,N>
operator+ (const tmatrix_expr<T,E,M,N> & e1, const tmatrix<T,M,N> & e2) {
	typedef tmatrix_expr_add<T,tmatrix_expr<T,E,M,N>,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add<T,tmatrix<T,M,N>,tmatrix_expr<T,E,M,N> >,M,N>
operator+ (const tmatrix<T,M,N> & e1, const tmatrix_expr<T,E,M,N> & e2) {
	typedef tmatrix_expr_add<T,tmatrix<T,M,N>,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E1, class E2, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_add<T,tmatrix_expr<T,E1,M,N>,tmatrix_expr<T,E2,M,N> >,M,N>
operator+ (const tmatrix_expr<T,E1,M,N> & e1, const tmatrix_expr<T,E2,M,N> & e2) {
	typedef tmatrix_expr_add<T,tmatrix_expr<T,E1,M,N>,tmatrix_expr<T,E2,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}

/*
 * Matrix subtraction
 */
template<typename T, class E1, class E2>
struct tmatrix_expr_sub {
	explicit tmatrix_expr_sub(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return m_e1(i,j) - m_e2(i,j);
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};

template<typename T, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_sub<T,tmatrix<T,M,N>,tmatrix<T,M,N> >,M,N> 
operator- (const tmatrix<T,M,N> & e1, const tmatrix<T,M,N> & e2) {
	typedef tmatrix_expr_sub<T,tmatrix<T,M,N>,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_sub<T,tmatrix_expr<T,E,M,N>,tmatrix<T,M,N> >,M,N>
operator- (const tmatrix_expr<T,E,M,N> & e1, const tmatrix<T,M,N> & e2) {
	typedef tmatrix_expr_sub<T,tmatrix_expr<T,E,M,N>,tmatrix<T,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_sub<T,tmatrix<T,M,N>,tmatrix_expr<T,E,M,N> >,M,N>
operator- (const tmatrix<T,M,N> & e1, const tmatrix_expr<T,E,M,N> & e2) {
	typedef tmatrix_expr_sub<T,tmatrix<T,M,N>,tmatrix_expr<T,E,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E1, class E2, unsigned M, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_sub<T,tmatrix_expr<T,E1,M,N>,tmatrix_expr<T,E2,M,N> >,M,N>
operator- (const tmatrix_expr<T,E1,M,N> & e1, const tmatrix_expr<T,E2,M,N> & e2) {
	typedef tmatrix_expr_sub<T,tmatrix_expr<T,E1,M,N>,tmatrix_expr<T,E2,M,N> > expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}

/*
 * Matrix product support
 */
template<typename T, unsigned K>
struct meta_multiply {
	template<class E1, class E2>
	static inline T prod(const E1 & e1, const E2 & e2,
			const unsigned i, const unsigned j) {
		return e1(i,K)*e2(K,j) + meta_multiply<T,K-1>::prod(e1,e2,i,j);
	}
};
template<typename T>
struct meta_multiply<T,0> {
	template<class E1, class E2>
	static inline T prod(const E1 & e1, const E2 & e2,
			const unsigned i, const unsigned j) {
		return e1(i,0)*e2(0,j);
	}
};

template<typename T, class E1, class E2, unsigned K>
struct tmatrix_expr_prod {
	explicit tmatrix_expr_prod(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline T operator() (const unsigned i, const unsigned j) const {
		return meta_multiply<T,K-1>::prod(m_e1,m_e2,i,j);
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};

template<typename T, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_prod<T,tmatrix<T,M,K>,tmatrix<T,K,N>,K>,M,N> 
operator* (const tmatrix<T,M,K> & e1, const tmatrix<T,K,N> & e2) {
	typedef tmatrix_expr_prod<T,tmatrix<T,M,K>,tmatrix<T,K,N>,K> expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_prod<T,tmatrix_expr<T,E,M,K>,tmatrix<T,K,N>,K>,M,N>
operator* (const tmatrix_expr<T,E,M,K> & e1, const tmatrix<T,K,N> & e2) {
	typedef tmatrix_expr_prod<T,tmatrix_expr<T,E,M,K>,tmatrix<T,K,N>,K> expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_prod<T,tmatrix<T,M,K>,tmatrix_expr<T,E,K,N>,K>,M,N>
operator* (const tmatrix<T,M,K> & e1, const tmatrix_expr<T,E,K,N> & e2) {
	typedef tmatrix_expr_prod<T,tmatrix<T,M,K>,tmatrix_expr<T,E,K,N>,K> expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}
template<typename T, class E1, class E2, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_prod<T,tmatrix_expr<T,E1,M,K>,tmatrix_expr<T,E2,K,N>,K>,M,N>
operator* (const tmatrix_expr<T,E1,M,K> & e1, const tmatrix_expr<T,E2,K,N> & e2) {
	typedef tmatrix_expr_prod<T,tmatrix_expr<T,E1,M,K>,tmatrix_expr<T,E2,K,N>,K> expr_t;
	return tmatrix_expr<T,expr_t,M,N>(expr_t(e1,e2));
}

/*
 * Support for matrix diag
 */
template<typename T, class E>
struct tmatrix_expr_diag {
	explicit tmatrix_expr_diag(const E & e) : m_e(e) { }
	inline T operator() (const unsigned i, const unsigned) const {
  		return m_e(i,i);
  	}
private:
	const E m_e;
};
template<typename T, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_diag<T,tmatrix<T,N,N> >,N,1> 
diag(const tmatrix<T,N,N> & e) {
	typedef tmatrix_expr_diag<T,tmatrix<T,N,N> > expr_t;
	return tmatrix_expr<T,expr_t,N,1>(expr_t(e));
}
template<typename T, class E, unsigned N>
inline tmatrix_expr<T,tmatrix_expr_diag<T,tmatrix_expr<T,E,N,N> >,N,1>
diag(const tmatrix_expr<T,E,N,N> & e) {
	typedef tmatrix_expr_diag<T,tmatrix_expr<T,E,N,N> > expr_t;
	return tmatrix_expr<T,expr_t,N,1>(expr_t(e));
}

#endif // TMATRIX_H
