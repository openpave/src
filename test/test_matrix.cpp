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

	The Original Code is OpenPave.org Core Libraries.

	The Initial Developer of the Original Code is OpenPave.org.

	Portions Copyright (C) 2006 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2002/11/05 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "matrix.h"

#ifdef NOBUILD
#define n 3
#define m 2
int
main()
{
	int i, j, k, r = 0;
	double s, dot = 0.0;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * b = new double[n];
	double * x = new double[n];
	if (A == 0 || B == 0 || b == 0 || x == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}

again:
	printf(".");
	memset(A,0,sizeof(double)*n*n);
	for (i = 0; i < n; i++) {
		for (j = i; j <= i+m && j < n; j++)
			A[i*n+j] = RAND(-1.0,1.0);
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			B[i*n+j] = 0.0;
			for (k = 0; k < n; k++)
				B[i*n+j] += A[i*n+k]*A[j*n+k];
		}
	}
	for (i = 0; i < n; i++)
		b[i] = RAND(0.0,1.0);
	printf("(%d",r++);
	fflush(NULL);
	
	memcpy(A,B,sizeof(double)*n*n);
	//equ_gauss(n,A,b,x);
	equ_lu(n,A,b,x);
	//equ_chol(n,A,b,x);
	//equ_ldl(n,A,b,x);
	//equ_svd(n,A,b,x);
	//equ_eig(n,A,b,x);

	//for (i = 0; i < n; i++) {
	//	for (j = i; j <= i+m && j < n; j++)
	//		A[B_IDX(n,m,i,j)] = B[i*n+j];
	//}
	//equ_chol(n,m,A,b,x,n*n*10e-12);
	
	printf(")");
	double c1, y1, t1, c2, y2, t2;
	for (i = 0, dot = 0.0, c1 = 0.0; i < n; i++) {
		for (j = 0, s = -b[i], c2 = 0.0; j < n; j++) {
			y2 = fma(B[i*n+j],x[j],-c2); t2 = s + y2;
			c2 = t2 - s - y2; s = t2;
		}
		y1 = fma(s,s,-c1); t1 = dot + y1;
		c1 = t1 - dot - y1; dot = t1;
	}
	if (sqrt(dot) > n*n*1e-8) {
		printf("\n%g\nA = [ ",sqrt(dot));
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.60e%s ...\n",B[i*n+j],(j == n-1 ? ";" : ","));
		}
		printf("];\n");
		printf(" b = [");
		for (i = 0; i < n; i++) {
			printf("%.60e; ...\n",b[i]);
		}
		printf("];\n x = [");
		for (i = 0; i < n; i++) {
			printf("%.60e; ...\n",x[i]);
		}
		printf("];\n");
		exit(1);
	}
	goto again;

abort:
	delete [] x;
	delete [] b;
	delete [] B;
	delete [] A;
}
#endif

#ifdef NOBUILD
int
main()
{
	int n, m, i, j, k;
	int * A;

	for (n = 0; n < 1000; n++) {
		for (m = 0; m < n; m++) {
			A = new int[B_SIZE(n,m)];
			memset(A,0,sizeof(int)*B_SIZE(n,m));
			for (i = 0; i < n; i++) {
				for (j = i; j >= i-m && j >= 0; j--) {
					k = B_IDX(n,m,j,i);
					if (k < 0 || k >= B_SIZE(n,m)) {
						printf("%d\t%d\t%d\t%d\t%d\t%d\n",n,m,i,j,B_SIZE(n,m),k);
						exit(1);
					}
					printf("%d\t%d\t%d\t%d\t%d\n",n,m,i,j,k);
					if (A[k] == 1) {
						printf("OVERLAP\n");
						exit(1);
					}
					A[k] = 1;
				}
			}
			for (i = 0; i < B_SIZE(n,m); i++) {
				if (A[i] != 1) {
					printf("MISSED\n");
					exit(1);
				}
			}
			delete [] A;
		}
	}
}
#endif

#ifdef NOBUILD
#define n 10
int
main()
{
	bool rv = false;
	int i, j, k, iter = 0;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * I = new double[n*n];
	if (A == 0 || B == 0 || I == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}

again:
	printf(".");
	for (i = 0; i < n*n; i++)
		B[i] = RAND(0.0,1.0), I[i] = A[i] = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				A[i*n+j] += B[i*n+k]*B[j*n+k];
		}
	}
	for (i = 0; i < n*n; i++)
		B[i] = A[i];
	printf("(");
	//inv_lu(n,A);
	//inv_chol(n,A);
	//inv_svd(n,A);
	inv_eig(n,A);
	printf("%d)",++iter);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				I[i*n+j] += A[i*n+k]*B[k*n+j];
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//printf("%f\t",I[i*n+j]-(i==j?1.0:0.0));
			if (fabs(I[i*n+j]-(i==j?1.0:0.0)) > 1e-6)
				rv = true;
		}
		//printf("\n");
	}
	if (rv)
		exit(1);
	goto again;

abort:
	if (A != 0)
		delete [] A;
	if (B != 0)
		delete [] B;
	if (I != 0)
		delete [] I;

}
#endif

#ifdef BUILD
#define n 3
int
main()
{
	bool rv = false;
	int i, j, k, iter = 0;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * I = new double[n*n];
	double * R = new double[n*n];
	if (A == 0 || B == 0 || I == 0 || R == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}
again:
	printf(".");
	for (i = 0; i < n*n; i++)
		R[i] = RAND(0.0,1.0), B[i] = A[i] = I[i] = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				A[i*n+j] += R[i*n+k]*R[j*n+k];
		}
	}
	for (i = 0; i < n*n; i++)
		R[i] = A[i];
	for (i = 0; i < n; i++)
		B[i*n+i] = 1.0;
	//printf("\nA = [ ");
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n; j++)
	//		printf("%.60e%s",A[i*n+j],(j == n-1 ? "; ...\n" : ", "));
	//}
	//printf("];\n");
	//printf("\nB = [ ");
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n; j++)
	//		printf("%.60e%s",B[i*n+j],(j == n-1 ? "; ...\n" : ", "));
	//}
	//printf("];\n");
	printf("(");
	//inv_mul_gauss(n,n,A,B);
	inv_mul_lu(n,n,A,B);
	printf("%d)",++iter);
	//printf("\nR = [ ");
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n; j++)
	//		printf("%.60e%s",B[i*n+j],(j == n-1 ? "; ...\n" : ", "));
	//}
	//printf("];\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				I[i*n+j] += B[i*n+k]*R[k*n+j];
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//printf("%f\t",I[i*n+j]-(i==j?1.0:0.0));
			if (fabs(I[i*n+j]-(i==j?1.0:0.0)) > 1e-6)
				rv = true;
		}
		//printf("\n");
	}
	if (rv) {
		printf("\nA = [ ");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.60e%s",R[i*n+j],(j == n-1 ? "; ...\n" : ", "));
		}
		printf("];\n");
		printf("\nI = [ ");
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.60e%s",I[i*n+j],(j == n-1 ? "; ...\n" : ", "));
		}
		printf("];\n");
		exit(1);
	}
	goto again;

abort:
	delete [] A;
	delete [] B;
	delete [] I;
	delete [] R;
}
#endif

#ifdef NOBUILD
#define n 5
int
main()
{
	int i, j, k;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * I = new double[n*n];
	double * d = new double[n];
	double * e = new double[n];
	if (A == 0 || d == 0 || e == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}

again:
	//printf(".");
	for (i = 0; i < n*n; i++)
		B[i] = RAND(0.0,1.0), I[i] = A[i] = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				A[i*n+j] += B[i*n+k]*B[j*n+k];
		}
	}
	for (i = 0; i < n*n; i++)
		B[i] = A[i];
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%8.4f%s",A[i*n+j],(j==n-1?"\n":"\t"));
	}
	printf("\n");
	//printf("(");
	tridiag_hh(n,A,d,e);
	eig_tri_ql(n,d,e,A);
	//printf(")");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%8.4f%s",A[i*n+j],(j==n-1?"\n":"\t"));
	}
	printf("\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				printf("%8.4f%s",d[i],(j==n-1?"\n":"\t"));
			else if (i == j+1)
				printf("%8.4f%s",e[i],(j==n-1?"\n":"\t"));
			else if (i == j-1)
				printf("%8.4f%s",e[i],(j==n-1?"\n":"\t"));
			else
				printf("%8.4f%s",0.0,(j==n-1?"\n":"\t"));
		}
	}
	//goto again;
abort:
	if (e != 0)
		delete [] e;
	if (d != 0)
		delete [] d;
	if (B != 0)
		delete [] B;
	if (A != 0)
		delete [] A;
	if (I != 0)
		delete [] I;
}
#endif

#ifdef NOBUILD
#define n 5
int
main()
{
	int i, j;

	double * A = new double[n*n];
	if (A == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		exit(1);
	}
	for (i = 0; i < n*n; i++)
		A[i] = RAND(0.0,1.0);

	matrix_dense *s = new matrix_dense(n,n,A);
	matrix a(s);
	matrix b;
	
	b = -a;
	matrix c(~a);

	printf("%g\t%g\t%g\n",a(1,1),b(1,1),c(1,1));
	printf("%g\t%g\t%g\n",a(1,5),b(1,5),c(5,1));

	delete [] A;
}
#endif

#ifdef NOBUILD
int
main()
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
	I = I + 1;
	printf(" I = "); I.print();
	I = I * 2;
	printf(" I = "); I.print();
	I = I / 2;
	printf(" I = "); I.print();
	I = I - 1;
	printf(" I = "); I.print();
	I = I * 0.5;
	printf(" I = "); I.print();
	I = 1.0 / I;
	printf(" I = "); I.print();
}
#endif
