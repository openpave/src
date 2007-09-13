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

#include "mathplus.h"
#include "set.h"
#include <stdio.h>

#ifdef NOBUILD
struct key {
	int i;
	key() : i(0) {}
	key(int s) : i(s) {}
	bool operator== (const key & k) const { return (i==k.i); }  
	bool operator<= (const key & k) const { return (i<=k.i); }  
	bool operator> (const key & k) const { return (i>k.i); }  
};
struct pair : public key {
	double d;
	pair() : key(), d(0.0) {}
	pair(int s, double v) : key(s), d(v) {}
};
struct value {
	double d;
	value() : d(0.0) {}
	value(double v) : d(v) {}
};

template<class T>
void fset_print(const fset<T> & t) {
	printf("t = {");
	for (int i = 0; i < t.length(); i++)
		printf("%i%s",t[i],(i == t.length()-1 ? "" : ","));
	printf("}\n");
}
template<class K, class V>
void kfset_print(const kfset<K,V> & t) {
	printf("k = {");
	for (int i = 0; i < t.length(); i++) {
		K & k = t[i];
		printf("%i:%4.2f%s",k.i,t[k].d,(i == t.length()-1 ? "" : ","));
	}
	printf("}\n");
}
template<class K, class V>
void afset_print(const afset<K,V> & t) {
	printf("a = {");
	for (int i = 0; i < t.length(); i++) {
		K & k = t.getkey(i);
		printf("%i:%4.2f%s",k.i,t[k].d,(i == t.length()-1 ? "" : ","));
	}
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
	t.add(d2,4);
	fset_print(t);
	t.sort();
	fset_print(t);
	t.add(2);
	t.add(2);
	t.add(3);
	fset_print(t);
	t.sort();
	t.remove();
	t.remove(0);
	fset_print(t);
	t.add(0,d1,7);
	fset_print(t);
	t.remove();
	t.remove(0);
	t.sort();
	fset_print(t);

	koset<key,pair> k(7);
	kfset_print(k);
	k.add(pair(2,1.5));
	k.add(pair(1,0.5));
	kfset_print(k);
	k.sort();
	kfset_print(k);
	k.add(pair(1,1.5));
	kfset_print(k);
	k.sort();
	kfset_print(k);

	aoset<key,value> a(7);
	afset_print(a);
	a.add(key(2),value(1.5));
	a.add(key(1),value(0.5));
	afset_print(a);
	a.sort();
	afset_print(a);
	a.add(key(1),value(1.5));
	afset_print(a);
	a.sort();
	afset_print(a);

	return 0;
}
#endif

#ifdef BUILD
#define N 30000
int
main()
{
	int i, j;
	oset<double> t(0,N);
again:
	t.empty();
	for (i = 0; i < N; i++)
		t.add(RAND(-N,N));
	//for (j = 0; j < N; j++)
	//	printf("%4.2f\n",t[j]);
	t.sort();
	printf(".");
	fflush(NULL);
	for (i = 1; i < N; i++) {
		if (t[i-1] > t[i]) {
			printf("Failed! (%i)\n",i);
			for (j = MAX(i-2,0); j <= MIN(i+2,N-1); j++)
				printf("%4.2f\n",t[j]);
			return 1;
		}
	}
	goto again;
	return 0;
}
#endif