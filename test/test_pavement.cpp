/**************************************************************************

	TEST_PAVEMENT.CPP - A test harness for LEsystem.

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

#include "pavement.h"
#include <stdio.h>

#ifdef NOBUILD
int
main()
{
	LEsystem Pavement;
	int i;

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
	cout << d.result(pavedata::deflct, pavedata::zz) << endl;
}
#endif

#ifdef NOBUILD
int
main()
{
	LEsystem Pavement;
	int i;

	Pavement.addload(point2d(0,0),0,1.0,1.0);
	Pavement.addpoint(point3d(0,0,0));
	for (i = 0; i < 100; i++)
		Pavement.addlayer(1.0,1.0*pow(i+0.5,2*0.1),0.5);
	Pavement.addlayer(0.0,1.0*pow(200.0,2*0.1),0.5);
	Pavement.calculate();
	const pavedata & d = Pavement.result(point3d(0,0,0));
	printf("Result: %f\n",d.result(pavedata::deflct, pavedata::zz));
}
#endif

#ifdef NOBUILD
int
main()
{
	cset<point3d> * _p = new cset<point3d>(0,100);
	cset<point3d> & p = *_p;
	LEsystem Best;
	LEsystem Slow;
	LEsystem Fast;
	int i;

	//Best.addlayer(3.0,1600.0*1e3,0.5);
	Best.addlayer(0.0,1.0,0.5);
	//Slow.addlayer(3.0,1600.0*1e3,0.5);
	Slow.addlayer(0.0,1.0,0.5);
	//Fast.addlayer(3.0,1600.0*1e3,0.5);
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
		const pavedata & b = Best.result(p[i]);
		const pavedata & s = Slow.result(p[i]);
		const pavedata & f = Fast.result(p[i]);
		printf("%10.4e",p[i].x);
		//double l = (p[i].x < 1 ? -1 : p[i].x > 1 ? 0 : -0.5);
		//printf("\t%+0.12e",b.result(pavedata::stress,pavedata::zz)-l);
		//printf("\t%+0.12e",s.result(pavedata::stress,pavedata::zz)-l);
		//printf("\t%+0.12e",f.result(pavedata::stress,pavedata::zz)-l);
		//double e = b.result(pavedata::stress,pavedata::zz)
		//	      -s.result(pavedata::stress,pavedata::zz);
		//printf("\t%+0.8e",e);
		printf("\t%+0.8e",b.result(pavedata::deflct,pavedata::zz));
		printf("\t%+0.8e",s.result(pavedata::deflct,pavedata::zz));
		printf("\t%+0.8e",f.result(pavedata::deflct,pavedata::zz));
		printf("\n");
	}
	delete _p;
}
#endif

#ifdef NOBUILD
int
main()
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
#endif

#ifdef BUILD
int
main()
{
	double t = 0;
	
	double E1i = 5800.0;
	double E2i =  150.0;
	double E3i =  150.0;
	double E4i =   90.0;
	do {
		LEsystem Pavement;
		int i, j;

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
#endif

#ifdef NOBUILD
int
main(int argc, char* argv[])
{
	int i;
	double t;

	LEsystem Test;
	LEsystem Fast;

	//srand((unsigned)time(NULL));
redo:
	Test.removelayers();
	Fast.removelayers();
	for (i = 0; i < 4; i++) {
		t = 25.0+50.0*(i+1)*((double)(rand())/(double)(RAND_MAX));
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
		T[i] = pow(10,4.0+3.0*(double)(rand())/(double)(RAND_MAX));
		while (i > 0 && (T[i] > 20*T[i-1] || T[i] < T[i-1]/20.0))
			T[i] = pow(10,4.0+3.0*(double)(rand())/(double)(RAND_MAX));
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
	for (i = 0; i < Test.data.length(); i++) {
		point3d & p = Test.data[i];
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
	goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
double quad8_vdp(double r, double z, double s, double v, double a = 0.0, double b = M_PI, double Q = 10.0);
int
main()
{
	int i;

	FILE * fp = fopen("real.dat","w");
	for (i = 0; i < 1<<16; i++) {
		double d = quad8_vdp(0.01*i,0.0,1.0,0.5);
		fprintf(fp,"%0.16f\t%0.16f\n",0.01*i,d);
	}
	fclose(fp);
}
#endif

#ifdef NOBUILD
int
main()
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
#endif
