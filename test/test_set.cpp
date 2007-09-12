/**************************************************************************

	TEST_SET.CPP - A test harness for set.h.

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

#include "set.h"
#include <stdio.h>

#ifdef BUILD

struct key {
	int i;
	key() : i(0) {}
};
struct pair : public key {
	double d;
	pair() : key(), d(0.0) {}
};
struct value {
	double d;
	value() : d(0.0) {}
};

template<class T>
void fset_print(const fset<T> & t) {
	printf("t = {");
	for (int i = 0; i < t.length(); i++)
		printf("%i%s",t[i],(i == t.length()-1 ? "" : ","));
	printf("}\n");
}

int
main()
{
	const int d1[7] = {1,2,2,3,4,4,5};
	const int d2[4] = {1,1,1,1};
	cset<int> t(7,10,d1);
	fset_print(t);
	t.empty();
	fset_print(t);
	t.add(7,d2);
	fset_print(t);
	t.sort();
	fset_print(t);
	t.add(2);
	t.add(2);
	fset_print(t);
	t.sort();
	fset_print(t);
	return 0;
}
#endif
