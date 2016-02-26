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

	Portions Copyright (C) 2006-2016 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2016/02/22 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "event.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "listen.h"

enum class test_event { test };

class test_source : public dispatcher<message<test_event>> {
public:
	test_source() {
		printf("Constructing source\n");
	}
	~test_source() {
		printf("Deleting source\n");
	}
	void event() {
		dispatch(test_event::test);
	}
};

class test_sink : public listener {
public:
	test_sink(test_source & d) {
		printf("Constructing sink\n");
		listen(d,message<test_event>(
			[this](test_event e){ this->onevent(e); }));
	}
	~test_sink() {
		printf("Deleting sink\n");
	}
	void onevent(test_event) {
		printf("Got event!\n");
	}
};

int
main()
{
	printf("Test 1:\n");
	test_source source;
	test_sink sink(source);
	printf("Sending event...\n");
	source.event();
	return 0;
}
