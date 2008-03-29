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

	The Initial Developer of the Original Software is Jeremy Lea.

	Portions Copyright (C) 2006-2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	History:
		2002/11/05 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "mathplus.h"
#include "set.h"
#include "tree.h"
#include <stdio.h>

#define N 10

struct key {
	int i;
	key() : i(0) {}
	key(int s) : i(s) {}
	bool operator== (const key & k) const { return (i==k.i); }  
	bool operator<  (const key & k) const { return (i< k.i); }  
	bool operator<= (const key & k) const { return (i<=k.i); }  
	bool operator>  (const key & k) const { return (i> k.i); }  
	bool operator>= (const key & k) const { return (i>=k.i); }  
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
	for (unsigned i = 0; i < t.length(); i++)
		printf("%i%s",t[i],(i == t.length()-1 ? "" : ","));
	printf("}\n");
}

void
test_set1a()
{
	const int d1[14] = {1,2,2,3,4,4,5,6,7,8,9,10,11,12};
	const int d2[4] = {1,1,1,1};
	cset<int> t(14,10,d1);
	fset_print(t);
	t.empty();
	t.resize(10);
	t.add(d1,14);
	t.add(d2,4);
	fset_print(t);
	t.sort();
	fset_print(t);
	t.empty();
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
}

template<class K, class V>
void kfset_print(const kfset<K,V> & t) {
	printf("k = {");
	for (unsigned i = 0; i < t.length(); i++) {
		K & k = t[i];
		printf("%i:%4.2f%s",k.i,t[k].d,(i == t.length()-1 ? "" : ","));
	}
	printf("}\n");
}

void
test_set1b()
{
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
}

template<class K, class V>
void kiset_print(const kiset<K,V> & t) {
	printf("i = {");
	for (unsigned i = 0; i < t.length(); i++) {
		const K & k = t.getindex(i);
		printf("%i:%i %4.2f%s",i,k.i,t[k].d,(i == t.length()-1 ? "" : ","));
	}
	printf("}\n");
}

void
test_set1c()
{
	kiset<key,pair> i(7);
	kiset_print(i);
	i.add(pair(2,1.5));
	i.add(pair(1,0.5));
	kiset_print(i);
	i.add(pair(1,1.5));
	kiset_print(i);
}

template<class K, class V>
void afset_print(const afset<K,V> & t) {
	printf("a = {");
	for (unsigned i = 0; i < t.length(); i++) {
		K & k = t.getkey(i);
		printf("%i:%4.2f%s",k.i,t[k].d,(i == t.length()-1 ? "" : ","));
	}
	printf("}\n");
}

void
test_set1d()
{
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
}

void
test_set2()
{
	unsigned i, j, l = 10;
	oset<double> t(0,N);
	while (l-- > 0) {
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
				for (j = (i<2?0:i-2); j < MIN(i+3,N); j++)
					printf("%4.2f\n",t[j]);
			}
		}
	}
}

void
test_set3()
{
	unsigned i, j, l = 10;
	iset<double> t(0,N);
	while (l-- > 0) {
		t.empty();
		for (i = 0; i < N; i++)
			t.add(RAND(-N,N));
		//for (j = 0; j < N; j++)
		//	printf("%4.2f\n",t[j]);
		printf(".");
		fflush(NULL);
		for (i = 1; i < N; i++) {
			if (t.getindex(i-1) > t.getindex(i)) {
				printf("Failed! (%i)\n",i);
				for (j = (i<2?0:i-2); j < MIN(i+3,N); j++)
					printf("%4.2f\n",t[j]);
			}
		}
	}
}

void
test_set4()
{
	unsigned i, l = 100;
	kiset<key,pair> t(0,N);

	while (l-- > 0) {
		t.empty();
		for (i = 0; i < N; i++) {
			int s = int(floor(RAND(0,N)));
			double v = double(i);
			//printf("Add %i:%i %4.2f\n",i,s,v);
			t.add(pair(s,v));
			//kiset_print(t);
		}
		//for (j = 0; j < N; j++)
		//	printf("%4.2f\n",t[j]);
		printf(".");
		fflush(NULL);
		for (i = 1; i < t.length(); i++) {
			if (t.getindex(i-1).i > t.getindex(i).i) {
				printf("Failed! (%i)\n",i);
				kiset_print(t);
			}
		}
	}
}

#ifdef NOBUILD
int
main()
{
	printf("Test 1a:\n");
	test_set1a();
	printf("\nTest 1b:\n");
	test_set1b();
	printf("\nTest 1c:\n");
	test_set1c();
	printf("\nTest 1d:\n");
	test_set1d();
	printf("\nTest 2:\n");
	test_set2();
	printf("\nTest 3:\n");
	test_set3();
	printf("\nTest 4:\n");
	test_set4();
	printf("\n");
	return 0;
}
#endif

#ifdef BUILD
int main()
{
	BST<key,value> bst;
	
	unsigned i, l = 1;
	oset<double> t(0,N);
	while (l-- > 0) {
		printf(".");
		fflush(NULL);
		for (i = 0; i < N; i++) {
			int s = int(floor(RAND(0,N)));
			double v = double(i);
			printf("Inserting %d: %f\n",s,v);
			bst.insert(key(s),value(v));
			bst.print();
		}
		//bst.print();
		//for (i = 0; i < N; i++) {
		//	bst.remove(key(i));
		//}
		//bst.print();
	}
	printf("\n");
	return 0;
}
#endif
