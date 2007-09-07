/**************************************************************************

	FEM3D.CPP - A special 3D Finite element code for my PhD.

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

	Portions Copyright (C) 2007 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2007/09/07 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#define _EVENT_IMP
#define _PROGRESS_IMP
#include "pavement.h"
#include "matrix.h"
#include <stdlib.h>
#include <time.h>

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
