/**************************************************************************

	TEST_MATRIX.CPP - A test harness for matrix.h.

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

	The Initial Developer of the Original Software is Jeremy Lea.

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2002/11/05 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "mathplus.h"
#include "matrix.h"

void print_matlab(const char * c, unsigned m, unsigned n, double * A) {
	printf("%s = [ ",c);
	for (unsigned i = 0; i < m; i++) {
		for (unsigned j = 0; j < n; j++)
			printf("%.60e%s ...\n",A[i*n+j],(j == n-1 ? ";" : ","));
	}
	printf("];\n");
}

void init_matrix_u(unsigned n, unsigned w, double * A) {
	memset(A,0,sizeof(double)*n*n);
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = i; j <= i+w && j < n; j++)
			A[i*n+j] = RAND(-1.0,1.0);
	}
}

void init_matrix_pd(unsigned n, double * A, double * B) {
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < n; j++) {
			B[i*n+j] = 0.0;
			for (unsigned k = 0; k < n; k++)
				B[i*n+j] += A[i*n+k]*A[j*n+k];
		}
	}
	memcpy(A,B,sizeof(double)*n*n);
}

void init_vector(unsigned n, double * b) {
	for (unsigned i = 0; i < n; i++)
		b[i] = RAND(0.0,1.0);
}

double residual(unsigned n, double * B, double * x, double * b) {
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
#if defined(DEBUG)
	if (dot > n*n*1e-8) {
		printf("\n%g\n",dot);
		print_matlab("A",n,n,B);
		print_matlab("b",n,1,b);
		print_matlab("x",n,1,x);
	}
#endif
	return dot;
}

double identity(unsigned n, double * A, double * B) {
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
#if defined(DEBUG)
	if (res > n*n*1e-8) {
		printf("\n%g\n",res);
		print_matlab("A",n,n,A);
		print_matlab("B",n,n,B);
	}
#endif
	return res;
}

#define N 10
#define W 3

void
test_matrix1()
{
	int i, j, iter = 30;
	double * A = new double[N*N];
	double * B = new double[N*N];
	double * b = new double[N];
	double * x = new double[N];
	if (A == 0 || B == 0 || b == 0 || x == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}
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
abort:
	delete [] x;
	delete [] b;
	delete [] B;
	delete [] A;
}

void
test_matrix2()
{
	unsigned iter = 30;
	double * A = new double[N*N];
	double * B = new double[N*N];
	if (A == 0 || B == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}
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
abort:
	delete [] A;
	delete [] B;
}

void
test_matrix3()
{
	unsigned iter = 30;
	double * A = new double[N*N];
	double * B = new double[N*N];
	double * I = new double[N*N];
	if (A == 0 || B == 0 || I == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}
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
abort:
	delete [] A;
	delete [] B;
	delete [] I;
}

void
test_matrix4()
{
	int i;

	double * A = new double[N*N];
	if (A == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		exit(1);
	}
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

void
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
	printf("D = ");	D.print();

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

#ifdef NOBUILD
int
main()
{
	printf("Test 1:\n");
	test_matrix1();
	printf("Test 2:\n");
	test_matrix2();
	printf("Test 3:\n");
	test_matrix3();
	printf("Test 4:\n");
	test_matrix4();
	printf("Test tmatrix:\n");
	test_tmatrix();
	return 0;
}
#endif
