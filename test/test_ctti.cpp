/**************************************************************************

	TEST_CTTI.CPP - A test harness for conststr.h and ctti.h.

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

	Portions Copyright (C) 2016-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2017/09/06 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "ctti.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace OP;

void
test1()
{
	printf("Test 1:\n");
	conststr v1("test");

	printf("test: %s\n",static_cast<const char *>(v1));
}

void
test2()
{
	printf("Test 2:\n");
	OP::type_info t1(type_id<int>());
	std::string v1{t1.name().operator std::string()};
	std::size_t h1 = t1.hash_code();

	printf("type: %s, hash: %zu\n",v1.c_str(),h1);
}

int
main()
{
	test1();
	test2();
	return 0;
}
