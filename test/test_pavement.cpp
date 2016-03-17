/**************************************************************************

	TEST_PAVEMENT.CPP - A test harness for LEsystem.

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

#define _EVENT_IMP
#define _PROGRESS_IMP
#include "event.h"
#include "pavement.h"
#include <stdio.h>

static void
test0()
{
	LEsystem Pavement;

	Pavement.addlayer(0.0,1.0,0.5);
	Pavement.addload(point2d(0,0),0,1.0,1.0);
	Pavement.addpoint(point3d(0,0,0));
	Pavement.calculate();
	const pavedata & d1 = Pavement.result(point3d(0,0,0));
	printf("First result: %f\n",d1.result(pavedata::deflct, pavedata::zz));
	Pavement.calculate();
	const pavedata & d2 = Pavement.result(point3d(0,0,0));
	printf("Second result: %f\n",d2.result(pavedata::deflct, pavedata::zz));
	Pavement.layer(0).emod(1.0);
	Pavement.calculate();
	const pavedata & d3 = Pavement.result(point3d(0,0,0));
	printf("Third result: %f\n",d3.result(pavedata::deflct, pavedata::zz));
	Pavement.removelayer(0);
	Pavement.addlayer(0.0,1.0,0.5);
	Pavement.calculate();
	const pavedata & d4 = Pavement.result(point3d(0,0,0));
	printf("Fourth result: %f\n",d4.result(pavedata::deflct, pavedata::zz));
}

static void
test1()
{
	LEsystem Pavement;
	unsigned i;

	//Pavement.addlayer(1000.0,1000.0,0.45);
	//for (i =0; i< 100; i++)
	//	Pavement.addlayer(10.0,1000.0,0.45);
	//Pavement.addlayer(1000.0, 500.0,0.35);
	//for (i =0; i< 100; i++)
	//	Pavement.addlayer(10.0, 500.0,0.35);
	Pavement.addlayer(   0.0, 200.0,0.5);
	//Pavement.addload(point2d(0,0),0,700.0,134.867);
	Pavement.addload(point2d(0,0),0,1.0,1.0);
	//Pavement.addpoint(point3d(0,0,0));
	//Pavement.addpoint(point3d(0,0,20000));
	for (i=0; i<=20; i++)
		Pavement.addpoint(point3d(i/20.0,0,0));
	for (i=1; i<=20; i++)
		Pavement.addpoint(point3d(1+i/20.0,0,0));
	for (i=1; i<99; i++)
		Pavement.addpoint(point3d(2+i,0,0));
	for (i=0; i<=20; i++)
		Pavement.addpoint(point3d(0,0,i*0.1));
	for (i=0; i<=18; i++)
		Pavement.addpoint(point3d(0,0,2+i));
	for (i=0; i<=18; i++)
		Pavement.addpoint(point3d(0,0,20+i*20));
	//Pavement.addpoint(point3d(0,0,500));
	//Pavement.addpoint(point3d(200,0,0));
	Pavement.calculate();
	const pavedata & d = Pavement.result(point3d(0,0,0));
	printf("%f\n",d.result(pavedata::deflct, pavedata::zz));
}

static void
test2()
{
	LEsystem Pavement;
	unsigned i;

	Pavement.addload(point2d(0,0),0,1.0,1.0);
	Pavement.addpoint(point3d(0,0,0));
	for (i = 0; i < 100; i++)
		Pavement.addlayer(1.0,1.0*pow(i+0.5,2*0.1),0.5);
	Pavement.addlayer(0.0,1.0*pow(200.0,2*0.1),0.5);
	Pavement.calculate();
	const pavedata & d = Pavement.result(point3d(0,0,0));
	printf("Result: %f\n",d.result(pavedata::deflct, pavedata::zz));
}

/*
 * This test adds a bunch of points at various depths below the center
 * of circular load on an infinite layer and compares the results of the
 * two methods (accurate and normal) against the exact solution.
 */
static void
test3()
{
	cset<point3d> p(0,1000);
	LEsystem Best;
	LEsystem Slow;
	LEsystem Fast;
	unsigned i;

	Best.addlayer(0.0,1.0,0.5);
	Slow.addlayer(0.0,1.0,0.5);
	Fast.addlayer(0.0,1.0,0.5);

	Best.addload(point2d(0.0,0.0),0,1.0,1.0);
	Slow.addload(point2d(0.0,0.0),0,1.0,1.0);
	Fast.addload(point2d(0.0,0.0),0,1.0,1.0);

	for (i = 0; i < 100; i++)
		p.add(point3d(0,0,double(i)/100.0));
	for (i = 10; i < 100; i++)
		p.add(point3d(0,0,double(i)/10.0));
	for (i = 10; i < 100; i++)
		p.add(point3d(0,0,double(i)));
	for (i = 10; i <= 100; i++)
		p.add(point3d(0,0,double(i)*10.0));
	p.sort();
	for (i = 0; i < p.length(); i++) {
		Best.addpoint(p[i]);
		Slow.addpoint(p[i]);
		Fast.addpoint(p[i]);
	}
	Best.calc_accurate();
	Slow.calculate(LEsystem::all);
	Fast.calc_fastnum();

	for (i = 0; i < p.length(); i++) {
		const pavedata & f = Fast.result(p[i]);
		const pavedata & b = Best.result(p[i]);
		const pavedata & s = Slow.result(p[i]);
		printf("%6.4g",p[i].z);
		printf("\t%+0.9e",f.result(pavedata::deflct,pavedata::zz));
		printf("\t%+0.9e",b.result(pavedata::deflct,pavedata::zz));
		printf("\t%+0.9e",s.result(pavedata::deflct,pavedata::zz));
		printf("\n");
	}
}

/*
 * This test adds a lot of points along the surface of a single
 * layer pavement, and compares the three basic solutions for this
 * problem (accurate, normal, numerical).
 */
static void
test4()
{
	cset<point3d> p(0,100);
	LEsystem Best;
	LEsystem Slow;
	LEsystem Fast;
	unsigned i;

	for (i = 0; i < 9; i++) {
		Best.addlayer(0.1,1.0,0.5);
		Slow.addlayer(0.1,1.0,0.5);
		Fast.addlayer(0.1,1.0,0.5);
	}
	Best.addlayer(0.0,1.0,0.5);
	Slow.addlayer(0.0,1.0,0.5);
	Fast.addlayer(0.0,1.0,0.5);

	Best.addload(point2d(0.0,0.0),0,1.0,1.0);
	Slow.addload(point2d(0.0,0.0),0,1.0,1.0);
	Fast.addload(point2d(0.0,0.0),0,1.0,1.0);

	p.add(point3d(0,0,0));
	p.add(point3d(1,0,0));
	for (i = 0; i < 50; i++) {
		double r = 1+pow(10,-3.0*(25-i)/25);
		p.add(point3d(r,0,0));
		p.add(point3d(1.0/r,0,0));
	}
	p.add(point3d(0,0,1000.0));
	p.add(point3d(1,0,1000.0));
	for (i = 0; i < 50; i++) {
		double r = 1+pow(10,-3.0*(25-i)/25);
		p.add(point3d(r,0,1000.0));
		p.add(point3d(1.0/r,0,1000.0));
	}
	p.sort();
	for (i = 0; i < p.length(); i++) {
		Best.addpoint(p[i]);
		Slow.addpoint(p[i]);
		Fast.addpoint(p[i]);
	}
	//Best.calc_accurate();
	Best.calculate(LEsystem::all);
	//Slow.calculate(LEsystem::all);
	Slow.calculate(LEsystem::fast);
	Fast.calculate(LEsystem::dirty);
	//Fast.calc_fastnum();

	for (i = 0; i < p.length(); i++) {
		const pavedata & b = Best.result(p[i]);
		const pavedata & s = Slow.result(p[i]);
		const pavedata & f = Fast.result(p[i]);
		printf("%8.4g",p[i].x);
		//double l = (p[i].x < 1 ? -1 : p[i].x > 1 ? 0 : -0.5);
		//printf("\t%+0.12e",b.result(pavedata::stress,pavedata::zz)-l);
		//printf("\t%+0.12e",s.result(pavedata::stress,pavedata::zz)-l);
		//printf("\t%+0.12e",f.result(pavedata::stress,pavedata::zz)-l);
		//double e = b.result(pavedata::stress,pavedata::zz)
		//       -s.result(pavedata::stress,pavedata::zz);
		//printf("\t%+0.8e",e);
		printf("\t%+0.8e",b.result(pavedata::deflct,pavedata::zz));
		printf("\t%+0.8e",s.result(pavedata::deflct,pavedata::zz));
		printf("\t%+0.8e",f.result(pavedata::deflct,pavedata::zz));
		printf("\n");
	}
}

static void
test5()
{
	cset<point3d> p(0,100);
	LEsystem pave;
	unsigned i, j;

	for (i = 0; i < 9; i++)
		pave.addlayer(0.1,1.0,0.5);
	pave.addlayer(0.0,1.0,0.5);

	pave.addload(point2d(0.0,0.0),0,1.0,1.0);

	for (j = 0; j < 20; j++) {
		double z = j*0.05;
		p.add(point3d(0,0,z));
		p.add(point3d(1,0,z));
		for (i = 0; i < 15; i++) {
			double r = 1+pow(10,-3.0*(25-i)/25);
			p.add(point3d(r,0,z));
			p.add(point3d(1.0/r,0,z));
		}
	}
	p.sort();
	for (i = 0; i < p.length(); i++)
		pave.addpoint(p[i]);
	pave.calculate(LEsystem::fast);

	for (i = 0; i < p.length(); i++) {
		const pavedata & f = pave.result(p[i]);
		printf("%8.4g",p[i].x);
		printf("\t%+0.8e",f.result(pavedata::deflct,pavedata::zz));
		printf("\n");
	}
}

static void
test6()
{
	//LEsystem Best;
	LEsystem Slow;
	//LEsystem Fast;
	int i;

	Slow.addlayer(  50.0,3000.0*1e3,0.44);
	Slow.addlayer( 150.0, 600.0*1e3,0.35);
	Slow.addlayer( 250.0, 250.0*1e3,0.35);
	Slow.addlayer(   0.0,  90.0*1e3,0.35);

	Slow.addload(point2d(   0.0,  0.0),20*1e6,720.0);
	Slow.addload(point2d(   0.0,350.0),20*1e6,720.0);

	for (i = -898; i < (1024-898); i++)
		Slow.addpoint(point3d(7.4*i,175.0,0.0));
	Slow.calculate(LEsystem::all);
	for (i = -898; i < (1024-898); i++) {
		const pavedata & s = Slow.result(point3d(7.4*i,175.0,0.0));
		printf("%10.4e\t%+0.8e\n",7.4*i,s.result(pavedata::deflct,pavedata::zz));
	}
}

static void
test7()
{
	double t = 0;

	double E1i = 5800.0;
	double E2i =  150.0;
	double E3i =  150.0;
	double E4i =   90.0;
	do {
		LEsystem Pavement;
		unsigned i, j;

		double E1 = E1i*(0.05+0.95*exp(-t/10.0));
		double E2 = E2i*(1+0.40*((E1-E1i)/E1i));
		double E3 = E3i*(1+0.35*((E1-E1i)/E1i));
		double E4 = E4i*(1+0.30*((E1-E1i)/E1i));
		//double E2 = E2i;
		//double E3 = E3i;
		//double E4 = E4i;

		Pavement.addlayer( 150.0,E1*1e3,0.35);
		Pavement.addlayer( 274.0,E2*1e3,0.35);
		Pavement.addlayer( 215.0,E3*1e3,0.35);
		Pavement.addlayer(   0.0,E4*1e3,0.35);

		Pavement.addload(point2d(   0.0,  0.0),20*1e6,720.0);
		Pavement.addload(point2d(   0.0,350.0),20*1e6,720.0);

		for (i = 0; i <= 100; i++) {
			for (j = 0; j <= 200; j++)
				Pavement.addpoint(point3d(0.0,i*5.0+175.0,j*5.0));
		}
		Pavement.calculate();
		//Pavement.odemark();
		for (i = 0; i <= 100; i++) {
			for (j = 0; j <= 200; j++) {
				const pavedata & d = Pavement.result(point3d(0.0,i*5.0+175.0,j*5.0));
				printf("%6.2f %6.2f %6.2f %6.2f ", E1, E2, E3, E4);
				printf("%7.2f %7.2f %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g %10.6g\n", d.y, d.z,
				d.result(pavedata::stress, pavedata::xx),
				d.result(pavedata::stress, pavedata::yy),
				d.result(pavedata::stress, pavedata::zz),
				d.result(pavedata::stress, pavedata::p1),
				d.result(pavedata::stress, pavedata::p3),
				d.result(pavedata::stress, pavedata::s1),
				d.result(pavedata::stress, pavedata::s3),
				d.result(pavedata::deflct, pavedata::yy),
				d.result(pavedata::deflct, pavedata::zz));
			}
		}
		printf("\n");
	} while (++t < 20);
}

static void
test8()
{
	unsigned i;
	double t;

	LEsystem Test;
	LEsystem Fast;

	//srand((unsigned)time(NULL));
//redo:
	Test.removelayers();
	Fast.removelayers();
	for (i = 0; i < 4; i++) {
		t = 25.0+RAND(0,50.0*(i+1));
		Test.addlayer( t,1000.0,0.5);
		Fast.addlayer( t,1000.0,0.5);
		//printf("%0.4g\n",t);
	}
	Test.addlayer( 0.0,1000.0,0.35);
	Fast.addlayer( 0.0,1000.0,0.35);
	//printf("0.0\n");

	Test.removeloads();
	Fast.removeloads();
	Test.addload(point2d(0.0,  0.0),20*1e6,690);
	Test.addload(point2d(0.0,305.0),20*1e6,520);
	Fast.addload(point2d(0.0,  0.0),20*1e6,690);
	Fast.addload(point2d(0.0,305.0),20*1e6,520);

	double T[5];
	for (i = 0; i < Test.layers(); i++) {
		T[i] = pow(10,RAND(4.0,7.0));
		while (i > 0 && (T[i] > 20*T[i-1] || T[i] < T[i-1]/20.0))
			T[i] = pow(10,RAND(4.0,7.0));
		Test.layer(i).emod(T[i]);
		Fast.layer(i).emod(T[i]);
		//printf("%0.4e\t%g\n",log10(T[i]),T[i]);
	}

	Test.removepoints();
	Fast.removepoints();
	for (i = 0; i < 10; i++) {
		//Test.addpoint(point3d(   26.5*(i-100),175.0,  0.0));
		Test.addpoint(point3d(   0.0,200.0, i*25.0));
		Fast.addpoint(point3d(   0.0,200.0, i*25.0));
		Test.addpoint(point3d(   0.0,95.0, i*25.0));
		Fast.addpoint(point3d(   0.0,95.0, i*25.0));
	}
	Fast.addpoint(point3d(0.0,0.0,0.0));
	//Test.addpoint(point3d(0.0,0.0,t));
	Fast.addpoint(point3d(   0.0,175.0,  0.0));
	Fast.addpoint(point3d( 100.0,175.0,  0.0));
	Fast.addpoint(point3d( 200.0,175.0,  0.0));
	Fast.addpoint(point3d( 400.0,175.0,  0.0));
	Fast.addpoint(point3d( 600.0,175.0,  0.0));
	Fast.addpoint(point3d( 900.0,175.0,  0.0));
	Fast.addpoint(point3d(1200.0,175.0,  0.0));
	Fast.addpoint(point3d(   0.0,175.0,  0.0));
	Fast.addpoint(point3d(   0.0,175.0, 50.0));
	Fast.addpoint(point3d(   0.0,175.0,150.0));
	Fast.addpoint(point3d(   0.0,175.0,300.0));
	Fast.addpoint(point3d(   0.0,175.0,500.0));
	Fast.addpoint(point3d(   0.0,  0.0,  0.0));
	Fast.addpoint(point3d(   0.0,  0.0, 50.0));
	Fast.addpoint(point3d(   0.0,  0.0,150.0));
	Fast.addpoint(point3d(   0.0,  0.0,300.0));
	Fast.addpoint(point3d(   0.0,  0.0,500.0));

	//printf("\n");
	Test.calculate();
	//Fast.calc_odemark();
	Fast.calculate(LEsystem::fast);
	for (i = 0; i < Test.results(); i++) {
		const point3d & p = Test.result(i);
		const pavedata & d = Test.result(p);
		printf("%7.2f\t%7.2f\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\n", d.y, d.z,
		d.result(pavedata::stress, pavedata::xx),
		d.result(pavedata::stress, pavedata::yy),
		d.result(pavedata::stress, pavedata::zz),
		d.result(pavedata::stress, pavedata::p1),
		d.result(pavedata::stress, pavedata::p3),
		d.result(pavedata::stress, pavedata::s1),
		d.result(pavedata::stress, pavedata::s3),
		d.result(pavedata::deflct, pavedata::yy),
		d.result(pavedata::deflct, pavedata::zz));
		const pavedata & f = Fast.result(p);
		printf("%7.2f\t%7.2f\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\t%10.6g\n", d.y, d.z,
		f.result(pavedata::stress, pavedata::xx),
		f.result(pavedata::stress, pavedata::yy),
		f.result(pavedata::stress, pavedata::zz),
		f.result(pavedata::stress, pavedata::p1),
		f.result(pavedata::stress, pavedata::p3),
		f.result(pavedata::stress, pavedata::s1),
		f.result(pavedata::stress, pavedata::s3),
		f.result(pavedata::deflct, pavedata::yy),
		f.result(pavedata::deflct, pavedata::zz));
	}
	printf("\n");
	//goto redo;
}

static void
test9()
{
	LEsystem Pavement;

	Pavement.addlayer(10.0,1.0,0.5,0.0001);
	Pavement.addlayer(0.0,1.0,0.5);
	Pavement.addload(point2d(0.0,0.0),1.0,0.0,1.0);
	Pavement.addpoint(point3d(0.0,0.0,0.0));
	Pavement.addpoint(point3d(10.0,0.0,9.99999999999999));
	Pavement.addpoint(point3d(10.0,0.0,10.0));
	Pavement.calculate();
	const pavedata & d0 = Pavement.result(point3d(0.0,0.0,0.0));
	printf("%0.16f\n",d0.result(pavedata::deflct,pavedata::zz));
	printf("%0.16f\n",d0.result(pavedata::deflct,pavedata::yy));
	printf("%0.16f\n",d0.result(pavedata::deflct,pavedata::xx));
	printf("%0.16f\n",d0.result(pavedata::stress,pavedata::zz));
	printf("%0.16f\n",d0.result(pavedata::stress,pavedata::yy));
	printf("%0.16f\n",d0.result(pavedata::stress,pavedata::xx));
	printf("%0.16f\n",d0.result(pavedata::stress,pavedata::xy));
	printf("%0.16f\n",d0.result(pavedata::stress,pavedata::xz));
	printf("%0.16f\n",d0.result(pavedata::stress,pavedata::yz));
	printf("\n");
	const pavedata & d1 = Pavement.result(point3d(10.0,0.0,9.99999999999999));
	printf("%0.16f\n",d1.result(pavedata::deflct,pavedata::zz));
	printf("%0.16f\n",d1.result(pavedata::deflct,pavedata::yy));
	printf("%0.16f\n",d1.result(pavedata::deflct,pavedata::xx));
	printf("%0.16f\n",d1.result(pavedata::stress,pavedata::zz));
	printf("%0.16f\n",d1.result(pavedata::stress,pavedata::yy));
	printf("%0.16f\n",d1.result(pavedata::stress,pavedata::xx));
	printf("%0.16f\n",d1.result(pavedata::stress,pavedata::xy));
	printf("%0.16f\n",d1.result(pavedata::stress,pavedata::xz));
	printf("%0.16f\n",d1.result(pavedata::stress,pavedata::yz));
	printf("\n");
	const pavedata & d2 = Pavement.result(point3d(10.0,0.0,10.0));
	printf("%0.16f\n",d2.result(pavedata::deflct,pavedata::zz));
	printf("%0.16f\n",d2.result(pavedata::deflct,pavedata::yy));
	printf("%0.16f\n",d2.result(pavedata::deflct,pavedata::xx));
	printf("%0.16f\n",d2.result(pavedata::stress,pavedata::zz));
	printf("%0.16f\n",d2.result(pavedata::stress,pavedata::yy));
	printf("%0.16f\n",d2.result(pavedata::stress,pavedata::xx));
	printf("%0.16f\n",d2.result(pavedata::stress,pavedata::xy));
	printf("%0.16f\n",d2.result(pavedata::stress,pavedata::xz));
	printf("%0.16f\n",d2.result(pavedata::stress,pavedata::yz));
}

static void
test10()
{
	unsigned l, i;
	double T[10], h[10], v[10];
	sset<point3d> P;

	LEsystem Base;
	LEsystem Test;

	//srand((unsigned)time(NULL));
//redo:
	l = unsigned(ceil(RAND(1,10)));
	Base.removelayers();
	Test.removelayers();
	for (i = 0; i < l; i++) {
		do {
			T[i] = pow(10,RAND((i<2?5.0:4.0),7.0));
		} while (i > 0 && (T[i] > 2.0*T[i-1] || 5.0*T[i] < T[i-1]));
		h[i] = RAND(25.0+i*50.0,75.0+i*100.0);
		if (i == l-1 && int(ceil(RAND(0.0,5.0))) < 5.0)
			h[i] = 0.0;
		v[i] = RAND(0.15,(l < 2 ? 0.45 : 0.4));
		Base.addlayer(h[i],T[i],v[i]);
		Test.addlayer(h[i],T[i],v[i]);
		//printf("Layer %d: %8.4f %12.6f %6.4f\n",i+1,h[i],T[i]/1000,v[i]);
	}

	Base.removeloads();
	Test.removeloads();
	Base.addload(point2d(0.0,0.0),20*1e6,690);
	Test.addload(point2d(0.0,0.0),20*1e6,690);

	P.empty();
	Base.removepoints();
	Test.removepoints();
	for (i = 0; i < 1; i++) {
		P.add(        point3d(   0.0, 0.0, i*25.0));
		//P.add(        point3d(   0.0,95.0, i*25.0));
		Base.addpoint(point3d(   0.0, 0.0, i*25.0));
		//Base.addpoint(point3d(   0.0,95.0, i*25.0));
		Test.addpoint(point3d(   0.0, 0.0, i*25.0));
		//Test.addpoint(point3d(   0.0,95.0, i*25.0));
	}

	Base.calc_accurate();
	Test.calculate();
	//Test.calculate(LEsystem::fast);
	for (i = 0; i < P.length(); i++) {
		point3d & p = P[i];
		const pavedata & f = Base.result(p);
		const pavedata & d = Test.result(p);
		double fz = f.result(pavedata::deflct, pavedata::zz);
		double dz = d.result(pavedata::deflct, pavedata::zz);
		double err = fabs(200.0*(fz-dz)/(fz+dz));
		if (err > 1.0) {
			printf("\n");
			printf("%7.2f\t%7.2f\t%10.6g\t", f.y, f.z, fz);
			printf("%7.2f\t%7.2f\t%10.6g\t", d.y, d.z, dz);
			printf("%10.6g\n",err);
			for (unsigned j = 0; j < l; j++)
				printf("Layer %d: %8.4f %12.6f %6.4f\n",j+1,h[j],T[j]/1000,v[j]);
			printf("\n");
			exit(1);
		} else
			printf(".");
	}
	fflush(NULL);
	//goto redo;
}

/*
 * MB Road under a single wheel.
 */
static void
test11()
{
	LEsystem Pavement;
	cset<point3d> p;

	Pavement.addlayer(  90.0,3000e3,0.35);
	Pavement.addlayer( 410.0, 750e3,0.35);
	Pavement.addlayer(4500.0, 100e3,0.35);

	Pavement.addload(point2d(0.0,0.0),0.0,690.0,100.0);

	//for (double x = 0.0; x <= 8192.0; x += 8.0) {
	//	for (double z = 0; z < 5000.0; z += 5.0) {
	//		p.add(point3d(x,0.0,z));
	//	}
	//}
	for (double x = -1024.0; x <= 1024.0; x += 8.0) {
		for (double y = -1024.0; y <= 1024.0; y += 8.0) {
			p.add(point3d(x,y,90.0-1e-6));
		}
	}
	//p.sort();
	for (unsigned i = 0; i < p.length(); i++)
		Pavement.addpoint(p[i]);

	Pavement.calculate();

	for (unsigned i = 0; i < p.length(); i++) {
		const pavedata & d = Pavement.result(p[i]);
		printf("%7.1f\t%7.1f\t%7.1f\t",d.x,d.y,d.z);
		printf("%0.16f\t",d.result(pavedata::deflct,pavedata::xx));
		printf("%0.16f\t",d.result(pavedata::deflct,pavedata::yy));
		printf("%0.16f\t",d.result(pavedata::deflct,pavedata::zz));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::xx));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::yy));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::zz));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::xy));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::xz));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::yz));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::p1));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::p2));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::p3));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::s1));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::s2));
		printf("%0.16f\t",d.result(pavedata::strain,pavedata::s3));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::xx));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::yy));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::zz));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::xy));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::xz));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::yz));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::p1));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::p2));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::p3));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::s1));
		printf("%0.16f\t",d.result(pavedata::stress,pavedata::s2));
		printf("%0.16f\n",d.result(pavedata::stress,pavedata::s3));
	}
}

int
main()
{
	printf("Test 0:\n");
	test0();
	printf("Test 1:\n");
	test1();
	printf("Test 2:\n");
	test2();
	printf("Test 3:\n");
	test3();
	printf("Test 4:\n");
	test4();
	printf("Test 5:\n");
	test5();
	printf("Test 6:\n");
	test6();
	printf("Test 7:\n");
	test7();
	printf("Test 8:\n");
	test8();
	printf("Test 9:\n");
	test9();
	printf("Test 10:\n");
	test10();
	printf("Test 11:\n");
	test11();
}
