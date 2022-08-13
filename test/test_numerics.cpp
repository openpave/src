/**************************************************************************

	TEST_NUMERICS.CPP - A test harness for numerics.h.

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

	Portions Copyright (C) 2017-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2017/05/09 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mathplus.h"
#include "numerics.h"

using namespace OP;

static void
test1()
{
	std::function<double(double)> f =
		[](double x) { return x*x-2*x; };

	double m = fmin_brent(f,0,2);
	printf("fmin_brent: %g should be %g\n",m,1.0);
	m = fmin_brent([](double x) { return x*x+2*x; });
	printf("fmin_brent: %g should be %g\n",m,-1.0);
	m = fmin_brent([](double x) { return std::exp(x)-4*x-(4-4*std::log(4)); });
	printf("fmin_brent: %g should be %g\n",m,1.38629437556774);
	m = fmin_brent([](double x) { return std::exp(x)-20*x+90; });
	printf("fmin_brent: %g should be %g\n",m,2.99573192325833);
}

static void
test2()
{
	double m;
	std::function<double(double)> f =
		[](double x) { return x*x-2*x; };

	m = fzero(f,-2,1);
	printf("fzero: %g should be %g\n",m,0.0);
	m = fzero(f,-1);
	printf("fzero: %g should be %g\n",m,0.0);
	m = fzero(f,-1e-18);
	printf("fzero: %g should be %g\n",m,0.0);
	m = fzero(f, 1e-18);
	printf("fzero: %g should be %g\n",m,0.0);
}

static void
test3()
{
	double m;

	m = fzero([](double x) { return std::exp(x)+x-2; },-5,4);
	printf("fzero: %g should be %g\n",m,0.44285440100239);
	m = fzero([](double x) { return std::exp(6*x-std::pow(x,4)-1)-1; },1,3);
	printf("fzero: %g should be %g\n",m,1.75777201824726);
	m = fzero([](double x) { return std::log(6*x-std::pow(x,4)); },0.1,1);
	printf("fzero: %g should be %g\n",m,0.16679566609859);
	m = fzero([](double x) { return std::sqrt(x)-4; },0.5);
	printf("fzero: %g should be %g\n",m,16.0000000000000);
	m = fzero([](double x) { return std::sqrt(x)-4; },0.5,40);
	printf("fzero: %g should be %g\n",m,16.0000000000000);
	m = fzero([](double x) { return std::sqrt(x)-4; },777,0,INFINITY);
	printf("fzero: %g should be %g\n",m,16.0000000000000);
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
	return 0;
}
