/**************************************************************************

	TEST_TABLE.CPP - A test harness for axis.h and table.h.

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

	Portions Copyright (C) 2016 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2016/03/05 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include <cstdio>
#include "mathplus.h"
#include "axis.h"
#include "table.h"

using namespace OP;

struct keyA {
	int i;
	keyA() : i(0) {}
	explicit keyA(int s) : i(s) {}
	int compare (const keyA & k) const { printf("."); return (i==k.i ? 0 : SGN(i-k.i)); }
	bool operator== (const keyA & k) const { return (i==k.i); }
	bool operator<  (const keyA & k) const { return (i< k.i); }
	bool operator<= (const keyA & k) const { return (i<=k.i); }
	bool operator>  (const keyA & k) const { return (i> k.i); }
	bool operator>= (const keyA & k) const { return (i>=k.i); }
};

static void
test_axis1()
{
	const int d1[12] = {1,2,3,4,5,6,7,8,9,10,11,12};
	const int d2[4] = {1,2,3,4};
	
	axis<keyA> ax;
	axis<keyA,keyA> ax2(ax);
	for (unsigned i = 0; i < 12; i++)
		ax.add(keyA(d1[i]));
	for (unsigned i = 0; i < 4; i++)
		ax2.add(keyA(d2[i]),keyA(1));
	ax2.add(keyA(1),keyA(12));
	keyA a, b;
	a = ax[1];
	std::tie(b,a) = ax2[1];
}

static void
test_axis2()
{
	const int d[4] = {1,2,3,4};

	axis<int> ax;
	axis<int,int> ax2(ax);
	axis<int,int,int> ax3(ax2);
	for (unsigned i = 0; i < 4; i++)
		ax.add(d[i]);
	for (unsigned i = 0; i < 4; i++)
		ax2.add(d[i],1);
	ax2.add(1,4);
	ax3.add(1,1,1);
	int a, b, c;
	a = ax[1];
	std::tie(b,a) = ax2[1];
	std::tie(c,b,a) = ax3[0];
}

static void
test_table1()
{
	const int d1[4] = {1,2,4,5};
	
	axis<keyA> ax;
	for (unsigned i = 0; i < 4; i++)
		ax.add(keyA(d1[i]));
	table<double,axis<keyA>> tbl(ax);
	ax.add(keyA(3));
	for (unsigned i = 0; i < 5; i++)
		tbl[i] = (double)(i);
	axis<unsigned,keyA> ax2(ax);
	for (unsigned i = 0; i < 4; i++)
		ax2.add(i,keyA(2));
	table<double,axis<keyA>,axis<unsigned,keyA>> tbl2(ax,ax2);
	(void)tbl2(1,1);
	(void)tbl2(keyA(2),ax2[2]);
}

int
main()
{
	std::printf("Test Axis 1:\n");
	test_axis1();
	std::printf("Test Axis 2:\n");
	test_axis2();
	std::printf("Test Table 1:\n");
	test_table1();
	return 0;
}
