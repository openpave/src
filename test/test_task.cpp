/**************************************************************************

	TEST_TASK.CPP - A test harness for task.h.

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
		2017/12/20 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "event.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "task.h"

using namespace OP;

void
test1()
{
	printf("Test 1:\n");
	task_queue<> t;
	
	t.enqueue([]() {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		return;
	});
}

void
test2()
{
	printf("Test 2:\n");
	task_queue<int> t;
	
	t.enqueue([]() -> int {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		return 1;
	});
	t.drain();
	t.enqueue([]() -> int {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		return 2;
	}, [](int r) {
		printf("Returned %i.\n",r);
	});
}

int
main()
{
	test1();
	test2();
	return 0;
}
