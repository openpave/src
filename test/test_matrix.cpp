/**************************************************************************

	TEST_MATRIX.CPP - A test harness for matrix.h.

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

	History:
		2002/11/05 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "mathplus.h"
#include "tmatrix.h"
#include "matrix.h"
#include "linalg.h"

using namespace OP;

static void
print_matlab(const char * c, unsigned m, unsigned n, double * A)
{
	printf("%s = [ ",c);
	for (unsigned i = 0; i < m; i++) {
		for (unsigned j = 0; j < n; j++)
			printf("%.60e%s ...\n",A[i*n+j],(j == n-1 ? ";" : ","));
	}
	printf("];\n");
}

static void
init_matrix_u(unsigned n, unsigned w, double * A)
{
	memset(A,0,sizeof(double)*n*n);
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = i; j <= i+w && j < n; j++)
			A[i*n+j] = RAND(-1.0,1.0);
	}
}

static void
init_matrix_pd(unsigned n, double * A, double * B)
{
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			B[i*n+j] = 0.0;
			for (unsigned k = 0; k < n; k++)
				B[i*n+j] += A[i*n+k]*A[j*n+k];
		}
	}
	memcpy(A,B,sizeof(double)*n*n);
}

static void
init_vector(unsigned n, double * b)
{
	for (unsigned i = 0; i < n; i++)
		b[i] = RAND(0.0,1.0);
}

static double
residual(unsigned n, double * B, double * x, double * b)
{
	unsigned i, j;
	double dot, s;
	double c1, y1, t1, c2, y2, t2;

	for (i = 0, dot = 0.0, c1 = 0.0; i < n; i++) {
		for (j = 0, s = -b[i], c2 = 0.0; j < n; j++) {
			y2 = fma(B[i*n+j],x[j],-c2); t2 = s + y2;
			c2 = t2 - s - y2; s = t2;
		}
		y1 = fma(s,s,-c1); t1 = dot + y1;
		c1 = t1 - dot - y1; dot = t1;
	}
	dot = sqrt(dot);
	if (dot > n*n*1e-8) {
		printf("\n%g\n",dot);
		print_matlab("A",n,n,B);
		print_matlab("b",n,1,b);
		print_matlab("x",n,1,x);
	}
	return dot;
}

static double
identity(unsigned n, double * A, double * B)
{
	unsigned i, j, k;
	double res = 0.0;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			double I = -(i==j?1.0:0.0);
			for (k = 0; k < n; k++)
				I += A[i*n+k]*B[k*n+j];
			if (fabs(I) > res)
				res = fabs(I);
		}
	}
	if (res > n*n*1e-8) {
		printf("\n%g\n",res);
		print_matlab("A",n,n,A);
		print_matlab("B",n,n,B);
	}
	return res;
}

#define N 10
#define M 5
#define W 3

static void
test1()
{
	int i, j, iter = 30;
	double * A = new double[N*N];
	double * B = new double[N*N];
	double * b = new double[N];
	double * x = new double[N];

	while (iter-- > 0) {
		init_matrix_u(N,W,A);
		init_matrix_pd(N,A,B);
		init_vector(N,b);
		printf("Repeat: %d\n",30-iter);
		fflush(NULL);

		equ_gauss(N,A,b,x);
		printf("gauss: %9.5g\n",residual(N,B,x,b));
		equ_lu(N,A,b,x);
		printf("lu:    %9.5g\n",residual(N,B,x,b));
		equ_chol(N,A,b,x);
		printf("chol:  %9.5g\n",residual(N,B,x,b));
		equ_ldl(N,A,b,x);
		printf("ldl:   %9.5g\n",residual(N,B,x,b));
		equ_svd(N,A,b,x);
		printf("svd:   %9.5g\n",residual(N,B,x,b));
		equ_eig(N,A,b,x);
		printf("eig:   %9.5g\n",residual(N,B,x,b));

		for (i = 0; i < N; i++) {
			for (j = i; j <= i+W && j < N; j++)
				A[B_IDX(N,W,i,j)] = B[i*N+j];
		}
		equ_chol(N,W,A,b,x);
		printf("bchol  %9.5g\n",residual(N,B,x,b));

		printf("\n");
	}
	delete [] x;
	delete [] b;
	delete [] B;
	delete [] A;
}

static void
test2()
{
	unsigned iter = 30;
	double * A = new double[N*N];
	double * B = new double[N*N];

	while (iter-- > 0) {
		init_matrix_u(N,N-1,A);
		init_matrix_pd(N,A,B);
		printf("Repeat: %d\n",30-iter);
		fflush(NULL);

		inv_lu(N,A);
		printf("lu:    %9.5g\n",identity(N,A,B));
		memcpy(A,B,N*N*sizeof(double));
		inv_chol(N,A);
		printf("chol:  %9.5g\n",identity(N,A,B));
		memcpy(A,B,N*N*sizeof(double));
		inv_svd(N,A);
		printf("svd:   %9.5g\n",identity(N,A,B));
		memcpy(A,B,N*N*sizeof(double));
		inv_eig(N,A);
		printf("eig:   %9.5g\n",identity(N,A,B));

		printf("\n");
	}
	delete [] A;
	delete [] B;
}

static void
test3()
{
	unsigned iter = 30;
	double * A = new double[N*N];
	double * B = new double[N*N];
	double * I = new double[N*N];

	while (iter-- > 0) {
		init_matrix_u(N,N-1,A);
		init_matrix_pd(N,A,B);
		printf("Repeat: %d\n",30-iter);
		fflush(NULL);

		memset(I,0,N*N*sizeof(double));
		for (unsigned i = 0; i < N; i++)
			I[i*N+i] = 1.0;
		inv_mul_gauss(N,N,A,I);
		printf("gauss:  %9.5g\n",identity(N,I,B));
		memcpy(A,B,N*N*sizeof(double));
		memset(I,0,N*N*sizeof(double));
		for (unsigned i = 0; i < N; i++)
			I[i*N+i] = 1.0;
		inv_mul_lu(N,N,A,I);
		printf("lu:     %9.5g\n",identity(N,I,B));

		printf("\n");
	}
	delete [] A;
	delete [] B;
	delete [] I;
}

static void
test4()
{
	int i;
	double * A = new double[N*N];

	for (i = 0; i < N*N; i++)
		A[i] = RAND(0.0,1.0);

		matrix_dense *s = new matrix_dense(N,N,A);
	matrix a(s);
	matrix b;

	b = -a;
	matrix c(~a);

	printf("%g\t%g\t%g\n",a(1,1),b(1,1),c(1,1));
	printf("%g\t%g\t%g\n",a(1,5),b(1,5),c(5,1));

	delete [] A;
}

static void
test5()
{
	int i, j;
	double * A = new double[N*M];
	double * B = new double[M*N];

	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++)
			B[j*M+i] = A[i*N+j] = RAND(0.0,1.0);
	}
	print_matlab("A",M,N,A);
	transpose(M,N,A);
	print_matlab("At",N,M,A);
	print_matlab("B",N,M,B);
	for (i = 0; i < M; i++) {
		for (j = 0; j < N; j++) {
			if (B[j*M+i] != A[j*M+i])
				printf("Transpose doesn't work!");
		}
	}
	delete [] A;
	delete [] B;
}

static void
test_tmatrix()
{
	tmatrix<double,3,2> B;
	tmatrix<double,3,3> D;

	B[0][0] = -0.05; B[0][1] =  0;
	B[1][0] =  0;    B[1][1] =  0.05;
	B[2][0] =  0.05; B[2][1] = -0.05;

	D[0][0] = 2000; D[0][1] = 1000; D[0][2] = 0;
	D[1][0] = 1000; D[1][1] = 2000; D[1][2] = 0;
	D[2][0] = 0;    D[2][1] = 0;    D[2][2] = 500;

	printf("B = "); B.print();
	printf("D = "); D.print();

	tmatrix<double,2,2> BT;
	BT = ~B*B;
	printf(" BT = "); BT.print();

	tmatrix<double,2,2> K;
	K = ~B*D*B;
	printf(" K = "); K.print();
	K = ~B*D*B + ~B*B;
	printf(" K = "); K.print();
	K -= BT;
	printf(" K = "); K.print();

	tmatrix<double,2,2> I;
	I[0][0] = 1; I[0][1] = 0;
	I[1][0] = 0; I[1][1] = 1;
	printf(" I = "); I.print();
	I = I + 1.0;
	printf(" I = "); I.print();
	I = I * 2.0;
	printf(" I = "); I.print();
	I = I / 2.0;
	printf(" I = "); I.print();
	I = I - 1.0;
	printf(" I = "); I.print();
	I = I * 0.5;
	printf(" I = "); I.print();
	I = 1.0 / I;
	printf(" I = "); I.print();
}

int
main()
{
	printf("Test 1:\n");
	test1();
	printf("Test 2:\n");
	test2();
	printf("Test 3:\n");
	test3();
	printf("Test 4:\n");
	test4();
	printf("Test 5:\n");
	test5();
	printf("Test tmatrix:\n");
	test_tmatrix();
	return 0;
}

#undef N
#undef M
#undef W

#if 0

/*
 * struct tensor_indexX - tensor indices
 *
 * We use a struct here rather than a raw char, because these create
 * distinct types, so the compiler sees the various functions as distinct
 * function signatures.  The 2, 3 and 4 valued ones are just for convenience.
 */
template<char I>
struct tensor_index { explicit tensor_index() {} };
template<char I, char J>
struct tensor_index2 { explicit tensor_index2() {} };
template<char I, char J, char K>
struct tensor_index3 { explicit tensor_index3() {} };
template<char I, char J, char K, char L>
struct tensor_index4 { explicit tensor_index4() {} };

/*
 * struct tensor2_expr - an encapsulated expression.
 *
 * This struct holds the result of an expression.  struct E is the result,
 * and then we have two dimensions and two implicit indices.  The struct E
 * is copied in the constructor, so that temporaries are not removed by
 * the compiler.
 */
template<class E, unsigned M, unsigned N, char I, char J>
struct tensor2_expr {
	inline explicit tensor2_expr(const E & e) : m_e(e) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j);
	}
private:
	const E m_e;
};

// Forward declare
template<unsigned M, unsigned N> struct tensor2;

/*
 * struct tensor2_ref - an encapsulated tensor2.
 *
 * This struct holds a reference to a tensor2, but adds the two implict
 * indices, so that it can be used in expressions.
 */
template<unsigned M, unsigned N, char I, char J>
struct tensor2_ref {
	inline explicit tensor2_ref(tensor2<M,N> & e) : m_e(e) {}
	inline double & operator() (const unsigned i, const unsigned j) {
		return m_e(i,j);
	}
	inline double operator() (const unsigned i, const unsigned j) const {
		return m_e(i,j);
	}
	tensor2_ref & operator= (const tensor2_ref & d) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				m_e(i+1,j+1) = d(i+1,j+1);
		return *this;
	}
	template<class E1>
	tensor2_ref & operator= (const tensor2_expr<E1,M,N,I,J> & d) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				m_e(i+1,j+1) = d(i+1,j+1);
		return *this;
	}
private:
	tensor2<M,N> & m_e;
};

/*
 * struct tensor2 - a second order tensor of size M by N.
 *
 * This is the real tensor.
 */
template<unsigned M, unsigned N>
struct tensor2 {
	explicit tensor2() {
	}
	explicit tensor2(double d) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				data[i][j] = d;
	}
	explicit tensor2(const tensor2 & d) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				data[i][j] = d.data[i][j];
	}
	template<char I, char J>
	const tensor2_expr<tensor2<M,N>,M,N,I,J> operator() (
			tensor_index2<I,J>) const {
		return tensor2_expr<tensor2<M,N>,M,N,I,J>(*this);
	}
	template<char I, char J>
	tensor2_ref<M,N,I,J> operator() (tensor_index2<I,J>) {
		return tensor2_ref<M,N,I,J>(*this);
	}
	inline double & operator() (const unsigned i, const unsigned j) {
		return data[i-1][j-1];
	}
	inline double operator() (const unsigned i, const unsigned j) const {
		return data[i-1][j-1];
	}
private:
	double data[M][N];
};

/*
 * struct tensor4_expr - an encapsulated expression.
 *
 * This struct holds the result of an expression.  struct E is the result,
 * and then we have four dimensions and four implicit indices.  The struct E
 * is copied in the constructor, so that temporaries are not removed by
 * the compiler.
 */
template<class E, unsigned M, unsigned N, unsigned P, unsigned Q,
		char I, char J, char K, char L>
struct tensor4_expr {
	inline explicit tensor4_expr(const E & e) : m_e(e) {}
	inline double operator() (const unsigned i, const unsigned j,
			const unsigned k, const unsigned l) const {
		return m_e(i,j,k,l);
	}
private:
	const E m_e;
};

// Forward declare
template<unsigned M, unsigned N, unsigned P, unsigned Q> struct tensor4;

/*
 * struct tensor4_ref - an encapsulated tensor4.
 *
 * This struct holds a reference to a tensor4, but adds the four implicit
 * indices, so that it can be used in expressions.
 */
template<unsigned M, unsigned N, unsigned P, unsigned Q,
		char I, char J, char K, char L>
struct tensor4_ref {
	inline explicit tensor4_ref(tensor4<M,N,P,Q> & e) : m_e(e) {}
	inline double & operator() (const unsigned i, const unsigned j,
			const unsigned k, const unsigned l) {
		return m_e(i,j,k,l);
	}
	inline double operator() (const unsigned i, const unsigned j,
			const unsigned k, const unsigned l) const {
		return m_e(i,j,k,l);
	}
private:
	tensor4<M,N,P,Q> & m_e;
};

/*
 * struct tensor4_trans_XXXX - support for general fourth order transpose
 */
template<unsigned M, unsigned N, unsigned P, unsigned Q,
		char I, char J, char K, char L>
struct tensor4_trans_ikjl {
	explicit tensor4_trans_ikjl(const tensor4<M,P,N,Q> & e) : m_e(e) { }
	inline double operator() (const unsigned i, const unsigned j,
			const unsigned k, const unsigned l) const {
		return m_e(i,k,j,l);
	}
private:
	const tensor4<M,P,N,Q> & m_e;
};

template<unsigned M, unsigned N, unsigned P, unsigned Q>
struct tensor4 {
public:
	explicit tensor4() {
	}
	explicit tensor4(double d) {
		for (unsigned i = 0; i < M; i++)
			for (unsigned j = 0; j < N; j++)
				for (unsigned k = 0; k < P; k++)
					for (unsigned l = 0; l < Q; l++)
						data[i][j][k][l] = d;
	}
	template<char I, char J, char K, char L>
	tensor4_ref<M,N,P,Q,I,J,K,L>
	operator() (tensor_index4<I,J,K,L>) {
		return tensor4_ref<M,N,P,Q,I,J,K,L>(*this);
	}
	template<char I, char J, char K, char L>
	const tensor4_expr<tensor4_trans_ikjl<M,P,N,Q,I,K,J,L>,M,P,N,Q,I,K,J,L>
	operator() (tensor_index4<I,J,K,L>, tensor_index4<I,K,J,L>) {
		return tensor4_expr<tensor4_trans_ikjl<M,P,N,Q,I,K,J,L>,
				M,P,N,Q,I,K,J,L>(tensor4_trans_ikjl<M,P,N,Q,I,K,J,L>(*this));
	}
	inline double & operator() (unsigned i, unsigned j, unsigned k, unsigned l) {
		return data[i-1][j-1][k-1][l-1];
	}
	inline double operator() (unsigned i, unsigned j, unsigned k, unsigned l) const {
		return data[i-1][j-1][k-1][l-1];
	}
private:
	double data[M][N][P][Q];
};

template<class E1, class E2, unsigned P, unsigned Q>
struct tensor_expr_prod_ijkl_kl {
	explicit tensor_expr_prod_ijkl_kl(const E1 & e1, const E2 & e2)
	  : m_e1(e1), m_e2(e2) {}
	inline double operator() (const unsigned i, const unsigned j) const {
		double m = 0.0;
		for (unsigned k = 0; k < P; k++)
			for (unsigned l = 0; k < Q; k++)
				m += m_e1(i,j,k,l)*m_e2(k,l);
		return m;
	}
private:
	const E1 m_e1;
	const E2 m_e2;
};
template<unsigned M, unsigned N, unsigned P, unsigned Q,
		char I, char J, char K, char L>
inline tensor2_expr<tensor_expr_prod_ijkl_kl<
	tensor4_ref<M,N,P,Q,I,J,K,L>,
	tensor2_ref<P,Q,K,L>,P,Q>,M,N,I,J>
operator* (const tensor4_ref<M,N,P,Q,I,J,K,L> & e1,
		const tensor2_ref<P,Q,K,L> & e2) {
	typedef tensor_expr_prod_ijkl_kl<
		tensor4_ref<M,N,P,Q,I,J,K,L>,
		tensor2_ref<P,Q,K,L>,P,Q> expr_t;
	return tensor2_expr<expr_t,M,N,I,J>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned N, unsigned P, unsigned Q,
		char I, char J, char K, char L>
inline tensor2_expr<tensor_expr_prod_ijkl_kl<
	tensor4_expr<E,M,N,P,Q,I,J,K,L>,
	tensor2_ref<P,Q,K,L>,P,Q>,M,N,I,J>
operator* (const tensor4_expr<E,M,N,P,Q,I,J,K,L> & e1,
		const tensor2_ref<P,Q,K,L> & e2) {
	typedef tensor_expr_prod_ijkl_kl<
		tensor4_expr<E,M,N,P,Q,I,J,K,L>,
		tensor2_ref<P,Q,K,L>,P,Q> expr_t;
	return tensor2_expr<expr_t,M,N,I,J>(expr_t(e1,e2));
}
template<class E, unsigned M, unsigned N, unsigned P, unsigned Q,
		char I, char J, char K, char L>
inline tensor2_expr<tensor_expr_prod_ijkl_kl<
	tensor4_ref<M,N,P,Q,I,J,K,L>,
	tensor2_expr<E,P,Q,K,L>,P,Q>,M,N,I,J>
operator* (const tensor4_ref<M,N,P,Q,I,J,K,L> & e1,
		const tensor2_expr<E,P,Q,K,L> & e2) {
	typedef tensor_expr_prod_ijkl_kl<
		tensor4_ref<M,N,P,Q,I,J,K,L>,
		tensor2_expr<E,P,Q,K,L>,P,Q> expr_t;
	return tensor2_expr<expr_t,M,N,I,J>(expr_t(e1,e2));
}

int main()
{
	tensor_index2<'i','j'> ij;
	tensor_index2<'k','l'> kl;
	tensor_index2<'i','k'> ik;
	tensor_index2<'j','l'> jl;
	tensor_index4<'i','j','k','l'> ijkl;
	tensor_index4<'i','k','j','l'> ikjl;
	tensor_index2<'p','q'> pq;
	const tensor2<3,3> e;
	tensor2<3,3> b;
	tensor2<3,3> s;
	tensor2<4,4> a;
	tensor4<3,3,3,3> E;
	tensor4<4,3,4,3> K;

	s(ij) = E(ijkl)*e(kl);
	tensor2<3,3> w(s);
	//w(pq) = E(ijkl)*e(kl); // Compile error...
	a(ik) = K(ijkl,ikjl)*b(jl);
}

#endif
