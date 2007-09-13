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
	template<class A, class E>
	static inline void assignR_add_d(A & lhs, const double & d) {
	    lhs(I,J) += d;
	    meta_assign<I,J-1>::assignR_add_d(lhs,d);
	}
	template<class A, class E>
	static inline void assign_add_d(A & lhs, const double & d) {
	    meta_assign<I,J>::assignR_add_d(lhs,d);
	    meta_assign<I-1,J>::assign_add_d(lhs,d);
	}
	template<class A, class E>
	static inline void assignR_mul_d(A & lhs, const double & d) {
	    lhs(I,J) *= d;
	    meta_assign<I,J-1>::assignR_mul_d(lhs,d);
	}
	template<class A, class E>
	static inline void assign_mul_d(A & lhs, const double & d) {
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
	template<class A, class E>
	static inline void assignR_add_d(A &, const E &) {}
	template<class A, class E>
	static inline void assignR_mul_d(A &, const E &) {}
};
template<unsigned J> struct meta_assign<-1,J> {
	template<class A, class E>
	static inline void assign(A &, const E &) {}
	template<class A, class E>
	static inline void assign_add(A &, const E &) {}
	template<class A, class E>
	static inline void assign_sub(A &, const E &) {}
	template<class A, class E>
	static inline void assign_add_d(A &, const E &) {}
	template<class A, class E>
	static inline void assign_mul_d(A &, const E &) {}
};

/*
 * struct tmatrix_expr - Encapsulates a matrix as an expression
 */
template<class E, unsigned M, unsigned N>
struct tmatrix_expr {
	inline explicit tmatrix_expr(const E & e) : m_e(e) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j);
	}
private:
	const E m_e;
};

/*
 * class tmatrix - A simple MxN matrix.
 */
template <unsigned M, unsigned N>
class tmatrix {
public:
	// A few constructors, for various uses...
	inline explicit tmatrix() {
	}
	inline explicit tmatrix(double d, bool eye = false) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				data[i][j] = (!eye || i == j ? d : 0.0);
	}
	inline explicit tmatrix(const double * v) {
		memcpy(data,v,M*N*sizeof(double));
	}
	inline explicit tmatrix(const tmatrix<M,N> & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
	}
	template<class E>
	inline explicit tmatrix(const tmatrix_expr<E,M,N> & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
	}
	~tmatrix() {
	}

	inline tmatrix<M,N> & operator= (const tmatrix<M,N> & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
		return *this;
	}
	template<class E>
	inline tmatrix<M,N> & operator= (const tmatrix_expr<E,M,N> & m) {
		meta_assign<M-1,N-1>::assign(*this,m);
		return *this;
	}
	inline tmatrix<M,N> & operator+= (const tmatrix<M,N> & m) {
		meta_assign<M-1,N-1>::assign_add(*this,m);
		return *this;
	}
	template<class E>
	inline tmatrix<M,N> & operator+= (const tmatrix_expr<E,M,N> & m) {
		meta_assign<M-1,N-1>::assign_add(*this,m);
		return *this;
	}
	inline tmatrix<M,N> & operator+= (const double & d) {
		meta_assign<M-1,N-1>::assign_add_d(*this,d);
		return *this;
	}
	inline tmatrix<M,N> & operator-= (const tmatrix<M,N> & m) {
		meta_assign<M-1,N-1>::assign_sub(*this,m);
		return *this;
	}
	template<class E>
	inline tmatrix<M,N> & operator-= (const tmatrix_expr<E,M,N> & m) {
		meta_assign<M-1,N-1>::assign_sub(*this,m);
		return *this;
	}
	inline tmatrix<M,N> & operator-= (const double & d) {
		meta_assign<M-1,N-1>::assign_add_d(*this,-d);
		return *this;
	}
	inline tmatrix<M,N> & operator*= (const double & d) {
		meta_assign<M-1,N-1>::assign_mul_d(*this,d);
		return *this;
	}
	inline tmatrix<M,N> & operator/= (const double & d) {
		meta_assign<M-1,N-1>::assign_mul_d(*this,1.0/d);
		return *this;
	}
	inline int rows() const {
		return M;
	}
	inline int cols() const {
		return N;
	}
	inline double operator() (const int r, const int c) const {
		return data[r][c];
	}
	inline double & operator() (const int r, const int c) {
		return data[r][c];
	}
	inline double * operator[] (const int r) {
		return data[r];
	}

	void print() const {
		printf("[\n");
		for (unsigned i = 0; i < M; ++i) {
			printf("\t[");
			for (unsigned j = 0; j < N; ++j)
				printf("\t%+6.4f", this->operator()(i,j));
			printf("]\n");
		}
		printf("]\n");
	}

private:
	double data[M][N];
};

/*
 * Support for matrix transpose
 */
template<class E>
struct tmatrix_expr_trans {
	explicit tmatrix_expr_trans(const E & e) : m_e(e) { }
	inline double operator() (const unsigned i, const unsigned j) const {
  		return m_e(j,i);
  	}
private:
	const E m_e;
};
template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_trans<tmatrix<M,N> >,N,M>
operator~ (const tmatrix<M,N> & m) {  
	typedef tmatrix_expr_trans<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,N,M>(expr_t(m));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_trans<tmatrix_expr<E,M,N> >,N,M>
operator~ (const tmatrix_expr<E,M,N> & m) {  
	typedef tmatrix_expr_trans<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,N,M>(expr_t(m));
}

/*
 * Support for matrix negate
 */
template<class E>
struct tmatrix_expr_neg {
	explicit tmatrix_expr_neg(const E & e) : m_e(e) { }
	inline double operator() (const unsigned i, const unsigned j) const {
  		return -m_e(i,j);
  	}
private:
	const E m_e;
};
template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_neg<tmatrix<M,N> >,M,N>
operator- (const tmatrix<M,N> & m) {  
	typedef tmatrix_expr_neg<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(m));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_neg<tmatrix_expr<E,M,N> >,M,N>
operator- (const tmatrix_expr<E,M,N> & m) {  
	typedef tmatrix_expr_neg<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(m));
}

/*
 * Add a double to a matrix
 */
template<class E>
struct tmatrix_expr_add_d {
	explicit tmatrix_expr_add_d(const E & e, const double & d)
	  : m_e(e), m_d(d) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j) + m_d;
	}
private:
	const E m_e;
	const double m_d;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add_d<tmatrix<M,N> >,M,N> 
operator+ (const tmatrix<M,N> & e, const double & d) {
	typedef tmatrix_expr_add_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add_d<tmatrix_expr<E,M,N> >,M,N>
operator+ (const tmatrix_expr<E,M,N> & e, const double & d) {
	typedef tmatrix_expr_add_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add_d<tmatrix<M,N> >,M,N> 
operator+ (const double & d, const tmatrix<M,N> & e) {
	typedef tmatrix_expr_add_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add_d<tmatrix_expr<E,M,N> >,M,N>
operator+ (const double & d, const tmatrix_expr<E,M,N> & e) {
	typedef tmatrix_expr_add_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add_d<tmatrix<M,N> >,M,N> 
operator- (const tmatrix<M,N> & e, const double & d) {
	typedef tmatrix_expr_add_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,-d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add_d<tmatrix_expr<E,M,N> >,M,N>
operator- (const tmatrix_expr<E,M,N> & e, const double & d) {
	typedef tmatrix_expr_add_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,-d));
}

/*
 * Subtract a matrix from a double
 */
template<class E>
struct tmatrix_expr_sub_d {
	explicit tmatrix_expr_sub_d(const E & e, const double & d)
	  : m_e(e), m_d(d) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_d - m_e(i,j);
	}
private:
	const E m_e;
	const double m_d;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_sub_d<tmatrix<M,N> >,M,N> 
operator- (const double & d, const tmatrix<M,N> & e) {
	typedef tmatrix_expr_sub_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_sub_d<tmatrix_expr<E,M,N> >,M,N>
operator- (const double & d, const tmatrix_expr<E,M,N> & e) {
	typedef tmatrix_expr_sub_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}

/*
 * Multiply a matrix by a double
 */
template<class E>
struct tmatrix_expr_mul_d {
	explicit tmatrix_expr_mul_d(const E & e, const double & d)
	  : m_e(e), m_d(d) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_d*m_e(i,j);
	}
private:
	const E m_e;
	const double m_d;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_mul_d<tmatrix<M,N> >,M,N> 
operator* (const tmatrix<M,N> & e, const double & d) {
	typedef tmatrix_expr_mul_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_mul_d<tmatrix_expr<E,M,N> >,M,N>
operator* (const tmatrix_expr<E,M,N> & e, const double & d) {
	typedef tmatrix_expr_mul_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_mul_d<tmatrix<M,N> >,M,N> 
operator* (const double & d, const tmatrix<M,N> & e) {
	typedef tmatrix_expr_mul_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_mul_d<tmatrix_expr<E,M,N> >,M,N>
operator* (const double & d, const tmatrix_expr<E,M,N> & e) {
	typedef tmatrix_expr_mul_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}

/*
 * Divide a matrix by a double
 */
template<class E>
struct tmatrix_expr_div_d {
	explicit tmatrix_expr_div_d(const E & e, const double & d)
	  : m_e(e), m_d(d) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j)/m_d;
	}
private:
	const E m_e;
	const double m_d;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_div_d<tmatrix<M,N> >,M,N> 
operator/ (const tmatrix<M,N> & e, const double & d) {
	typedef tmatrix_expr_div_d<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_mul_d<tmatrix_expr<E,M,N> >,M,N>
operator/ (const tmatrix_expr<E,M,N> & e, const double & d) {
	typedef tmatrix_expr_div_d<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}

/*
 * Divide a double by a matrix elementwise
 */
template<class E>
struct tmatrix_expr_d_div {
	explicit tmatrix_expr_d_div(const E & e, const double & d)
	  : m_e(e), m_d(d) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_d/m_e(i,j);
	}
private:
	const E m_e;
	const double m_d;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_d_div<tmatrix<M,N> >,M,N> 
operator/ (const double & d, const tmatrix<M,N> & e) {
	typedef tmatrix_expr_d_div<tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_d_div<tmatrix_expr<E,M,N> >,M,N>
operator/ (const double & d, const tmatrix_expr<E,M,N> & e) {
	typedef tmatrix_expr_d_div<tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e,d));
}

/*
 * Matrix addition
 */
template<class E1, class E2>
struct tmatrix_expr_add {
	explicit tmatrix_expr_add(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e1(i,j) + m_e2(i,j);
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add<tmatrix<M,N>,tmatrix<M,N> >,M,N> 
operator+ (const tmatrix<M,N> & e1, const tmatrix<M,N> & e2) {
	typedef tmatrix_expr_add<tmatrix<M,N>,tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add<tmatrix_expr<E,M,N>,tmatrix<M,N> >,M,N>
operator+ (const tmatrix_expr<E,M,N> & e1, const tmatrix<M,N> & e2) {
	typedef tmatrix_expr_add<tmatrix_expr<E,M,N>,tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add<tmatrix<M,N>,tmatrix_expr<E,M,N> >,M,N>
operator+ (const tmatrix<M,N> & e1, const tmatrix_expr<E,M,N> & e2) {
	typedef tmatrix_expr_add<tmatrix<M,N>,tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E1, class E2, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_add<tmatrix_expr<E1,M,N>,tmatrix_expr<E2,M,N> >,M,N>
operator+ (const tmatrix_expr<E1,M,N> & e1, const tmatrix_expr<E2,M,N> & e2) {
	typedef tmatrix_expr_add<tmatrix_expr<E1,M,N>,tmatrix_expr<E2,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}

/*
 * Matrix subtraction
 */
template<class E1, class E2>
struct tmatrix_expr_sub {
	explicit tmatrix_expr_sub(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e1(i,j) - m_e2(i,j);
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};

template<unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_sub<tmatrix<M,N>,tmatrix<M,N> >,M,N> 
operator- (const tmatrix<M,N> & e1, const tmatrix<M,N> & e2) {
	typedef tmatrix_expr_sub<tmatrix<M,N>,tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_sub<tmatrix_expr<E,M,N>,tmatrix<M,N> >,M,N>
operator- (const tmatrix_expr<E,M,N> & e1, const tmatrix<M,N> & e2) {
	typedef tmatrix_expr_sub<tmatrix_expr<E,M,N>,tmatrix<M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_sub<tmatrix<M,N>,tmatrix_expr<E,M,N> >,M,N>
operator- (const tmatrix<M,N> & e1, const tmatrix_expr<E,M,N> & e2) {
	typedef tmatrix_expr_sub<tmatrix<M,N>,tmatrix_expr<E,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E1, class E2, unsigned M, unsigned N>
inline tmatrix_expr<tmatrix_expr_sub<tmatrix_expr<E1,M,N>,tmatrix_expr<E2,M,N> >,M,N>
operator- (const tmatrix_expr<E1,M,N> & e1, const tmatrix_expr<E2,M,N> & e2) {
	typedef tmatrix_expr_sub<tmatrix_expr<E1,M,N>,tmatrix_expr<E2,M,N> > expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}

/*
 * Matrix product support
 */
template<unsigned K>
struct meta_multiply {
	template<class E1, class E2>
	static inline double prod(const E1 & e1, const E2 & e2,
			const unsigned i, const unsigned j) {
		return e1(i,K)*e2(K,j) + meta_multiply<K-1>::prod(e1,e2,i,j);
	}
};
template<>
struct meta_multiply<0> {
	template<class E1, class E2>
	static inline double prod(const E1 & e1, const E2 & e2,
			const unsigned i, const unsigned j) {
		return e1(i,0)*e2(0,j);
	}
};

template<class E1, class E2, unsigned K>
struct tmatrix_expr_prod {
	explicit tmatrix_expr_prod(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return meta_multiply<K-1>::prod(m_e1,m_e2,i,j);
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};

template<unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<tmatrix_expr_prod<tmatrix<M,K>,tmatrix<K,N>,K>,M,N> 
operator* (const tmatrix<M,K> & e1, const tmatrix<K,N> & e2) {
	typedef tmatrix_expr_prod<tmatrix<M,K>,tmatrix<K,N>,K> expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<tmatrix_expr_prod<tmatrix_expr<E,M,K>,tmatrix<K,N>,K>,M,N>
operator* (const tmatrix_expr<E,M,K> & e1, const tmatrix<K,N> & e2) {
	typedef tmatrix_expr_prod<tmatrix_expr<E,M,K>,tmatrix<K,N>,K> expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<tmatrix_expr_prod<tmatrix<M,K>,tmatrix_expr<E,K,N>,K>,M,N>
operator* (const tmatrix<M,K> & e1, const tmatrix_expr<E,K,N> & e2) {
	typedef tmatrix_expr_prod<tmatrix<M,K>,tmatrix_expr<E,K,N>,K> expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}
template<class E1, class E2, unsigned M, unsigned K, unsigned N>
inline tmatrix_expr<tmatrix_expr_prod<tmatrix_expr<E1,M,K>,tmatrix_expr<E2,K,N>,K>,M,N>
operator* (const tmatrix_expr<E1,M,K> & e1, const tmatrix_expr<E2,K,N> & e2) {
	typedef tmatrix_expr_prod<tmatrix_expr<E1,M,K>,tmatrix_expr<E2,K,N>,K> expr_t;
	return tmatrix_expr<expr_t,M,N>(expr_t(e1,e2));
}

#endif // __TMATRIX_H