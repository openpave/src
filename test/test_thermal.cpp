/**************************************************************************

	TEST_THERMAL.CPP - A test harness for OpenPave.org code.

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

	Portions Copyright (C) 2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2008/11/21 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "thermal.h"
#include "mathplus.h"
#include <stdio.h>

using namespace OP;

#define NT 800
#define TB 15.0
#define ALPHA 2000.0

static double
solution(double t, double z)
{
	double d = 1.0/24.0;
	double y = 1.0/24.0/365.0;
	double D = sqrt(2*ALPHA/d);
	double Y = sqrt(2*ALPHA/y);
	return TB+10.0*exp(-z/Y)*sin(t*y-z/Y)+5.0*exp(-z/D)*sin(t*d-z/D);
}

int
main()
{
	int t;
	double tt, tb = TB;
	double h, D = ALPHA;
	double nd[NT+1], nt[NT+1];
	double pd = 50.0, pt, at;

	for (t = 0; t < NT+1; t++) {
		nd[t] = (2000.0/NT)*t;
		nt[t] = solution(0,nd[t]);
	}
	h = nd[NT];
	FEMthermal * system = new FEMthermal(1,&h,&D,NT+1,nd,nt,1,1);

	for (t = 1; t < 24*365*30; t++) {
		tt = solution(t,0);
		tb = solution(t,h);
		system->step(tt,tb);
		system->interpolate(1,&pd,&pt);
		at = solution(t,pd);
		printf("%10.6g\t%10.6g\t%10.6g\n",pt,at,pt-at);
	}
	delete system;
}
