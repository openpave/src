/**************************************************************************

	TEST_VARIANT.CPP - A test harness for variant.h.

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

	Portions Copyright (C) 2006-2016 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2016/02/22 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "event.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "variant.h"

using namespace OP;

double f1() { return 0.4; };
int f2() { return 5; };

void
test1()
{
	printf("Test 1:\n");
	variant<int,double> v1(1.0);
	variant<int,double> v2(2);
	//variant<int,double> v3(nullptr);

	double d = v1;
	int i = v2;
	//double e = v2;
	//double * p = v1;
	printf("v1:%f, v2: %i\n",d,i);
	v1 = 0.5;
	v2 = 3;
	printf("v1:%f, v2: %i\n",(double)(v1),(int)(v2));
}

void
test2()
{
	printf("Test 2:\n");
	vunctor<int,double> v1(1.0);
	vunctor<int,double> v2(2);

	double d = v1;
	int i = v2;
	printf("v1:%f, v2: %i\n",d,i);
	v1 = 0.5;
	v2 = 3;
	printf("v1:%f, v2: %i\n",(double)(v1),(int)(v2));
	v1 = []() -> double { return 0.3; };
	v2 = []() -> int { return 4; };
	printf("v1:%f, v2: %i\n",(double)(v1),(int)(v2));
	v1 = f1;
	v2 = f2;
	printf("v1:%f, v2: %i\n",(double)(v1),(int)(v2));
	std::function<double()> s1(f1);
	std::function<int()> s2(f2);
	v1 = s1;
	v2 = s2;
	printf("v1:%f, v2: %i\n",(double)(v1),(int)(v2));
}

double t1(const double & d) { return d == 0.3; };
int t2(const int & i) { return i == 2; };

void
test3()
{
	printf("Test 3:\n");
	variant<int,double> v1(1.0);
	variant<int,double> v2(2);
	validator<int,double> f1([](const double & d) -> bool {
		return d == 0.3;
	});
	validator<int,double> f2([](const int & i) -> bool {
		return i == 2;
	});

	bool b1 = f1.validate(v1);
	bool b2 = f2.validate((int)(v2));
	printf("v1:%s, v2: %s\n",b1 ? "T" : "F",b2 ? "T" : "F");
}

int
main()
{
	test1();
	test2();
	test3();
	return 0;
}
