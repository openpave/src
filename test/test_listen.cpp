/**************************************************************************

	TEST_LISTEN.CPP - A test harness for listen.h.

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
		2016/02/22 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "event.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "listen.h"

using namespace OP;

enum class test_event { test };

class test_dispatcher : public dispatcher<message<test_event>> {
public:
	test_dispatcher() {
		printf("Constructing test dispatcher\n");
	}
	~test_dispatcher() {
		printf("Deleting test dispatcher\n");
	}
	void event() {
		dispatch(test_event::test);
	}
};

class test_listener : public listener {
public:
	test_listener() {
		printf("Constructing listener\n");
	}
	void test_listen(test_dispatcher & d) {
		printf("Attaching listener\n");
		listen(d,message<test_event>(
			[this](test_event e){ this->onevent(e); }));
	}
	~test_listener() {
		printf("Deleting listener\n");
	}
	void onevent(test_event) {
		printf("Got event!\n");
	}
};

void
test1()
{
	printf("Test 1:\n");
	test_dispatcher dispatcher;
	test_listener listener;
	listener.test_listen(dispatcher);
	printf("Sending event...\n");
	dispatcher.event();
}

void
test2()
{
	printf("Test 2:\n");
	test_dispatcher * dispatcher1 = new test_dispatcher();
	test_dispatcher * dispatcher2 = new test_dispatcher();
	test_listener * listener1 = new test_listener();
	listener1->test_listen(*dispatcher1);
	test_listener * listener2 = new test_listener();
	listener2->test_listen(*dispatcher2);
	listener2->test_listen(*dispatcher1);
	printf("Sending event...\n");
	dispatcher1->event();
	delete listener1;
	delete listener2;
	delete dispatcher1;
}

int
main()
{
	test1();
	test2();
	return 0;
}
