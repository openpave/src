/**************************************************************************

	TEST_RNG.CPP - A test harness for rng.h.

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

	The Initial Developer of the Original Software is Jeremy Lea.

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2008/09/12 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "event.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "rng.h"

#ifdef NOBUILD
int
main(int argc, char * argv[])
{
	int i, cnt, seed;
	double x, y, pi;
	const int NUM = 1000000;

	if (argc >= 2)
		seed = strtol(argv[1], NULL, 10);
	else
		seed = 12345;
	cnt = 0;
	rng RNG(seed);
	for (i = 0; i < NUM; i++) {
		x = RNG.c0o1();
		y = RNG.c0o1();
		if (x*x + y*y < 1.0)
			cnt++;
	}
	pi = double(cnt) / NUM * 4;
	printf("%f\n", pi);
	return 0;
}
#endif
