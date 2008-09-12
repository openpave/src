/**************************************************************************

	TEST_BACKCALC.CPP - A test harness for LEbackcalc.

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

#include "pavement.h"
#include "reliability.h"
#include <stdlib.h>
#include <stdio.h>

#ifdef NOBUILD
int
main()
{
	unsigned i, l;
	LEbackcalc Bowl;
	cset<point3d> p;
	p.add(point3d(   0.0,0.0,  0.0));
	p.add(point3d( 200.0,0.0,  0.0));
	p.add(point3d( 300.0,0.0,  0.0));
	p.add(point3d( 460.0,0.0,  0.0));
	p.add(point3d( 610.0,0.0,  0.0));
	p.add(point3d( 910.0,0.0,  0.0));
	p.add(point3d(1525.0,0.0,  0.0));

	//srand((unsigned)time(NULL));
	double T[5], h[5], v[5], mz[7], t, err;

redo:
	l = unsigned(floor(RAND(2,6)));
	//Bowl.setup(0.0,0.0,1e-6,5);
	Bowl.setup(0.0001,0.0005,1e-4,5);

	Bowl.removelayers();
	for (i = 0; i < l; i++) {
		do {
			T[i] = pow(10,RAND((i<2?5.0:4.0),7.0));
		} while (i > 0 && (T[i] > 2.0*T[i-1] || 5.0*T[i] < T[i-1]));
		h[i] = RAND(25.0+i*50.0,75.0+i*100.0);
		if (i == l-1 && RAND(0.0,5.0) < 5.0)
			h[i] = 0.0;
		v[i] = RAND(0.3,(l < 2 ? 0.5 : 0.4));
		Bowl.addlayer(h[i],T[i],v[i]);
		printf("Layer %d: %8.4f %12.6f %6.4f\n",i+1,h[i],T[i]/1000,v[i]);
	}
	Bowl.removeloads();
	Bowl.addload(point2d(0.0,0.0),40*1e6,0.0,150.0);

	Bowl.removepoints();
	for (i = 0; i < p.length(); i++)
		Bowl.addpoint(p[i]);

	Bowl.calculate(LEsystem::disp);
	for (i = 0; i < p.length(); i++) {
		mz[i] = Bowl.result(p[i]).result(pavedata::deflct,pavedata::zz);
		printf("Measured defl. at position %d = %f\n",i+1,mz[i]);
		if (mz[i] < 0.0) {
			printf("Not a nice pavement! Bailing...\n");
			goto redo;
		}
	}

	Bowl.removedeflections();
	printf("M! ");
	for (i = 0; i < Bowl.layers(); i++)
		printf(" %0.4g",log10(Bowl.layer(i).emod()));
	for (i = 0; i < p.length(); i++) {
		t = Bowl.result(p[i]).result(pavedata::deflct,pavedata::zz);
		t += RAND(2,5)*0.0001*stdnormal_rnd();
		t = ROUND(t/0.0001)*0.0001;
		Bowl.adddefl(p[i],t);
		printf(" %0.4g",t);
	}
	printf("\n");
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(0.0);

	bool rv = Bowl.backcalc();
	printf("\n%s! (with %i layers)\n",(rv?"DONE":"FAILED"),l);

	for (i = 0; i < Bowl.layers(); i++) {
		t = log10(Bowl.layer(i).emod()/T[i]);
		printf("Layer %d/%d:\t%6.1f\t%0.2f\t%6.1f\t",i+1-l+5,l,h[i],v[i],T[i]/1000);
		printf("%6.1f\t",Bowl.layer(i).emod()/1000);
		if (rv) {
			if (fabs(t) > log10(2.0))
				printf("W\t");
			else
				printf("S\t");
		} else {
			if (fabs(t) > log10(2.0))
				printf("P\t");
			else
				printf("F\t");
		}
		printf("%g\n",(Bowl.layer(i).emod()/T[i]));
	}
	for (i = 0, err = 0.0; i < p.length(); i++) {
		const defldata & d = Bowl.getdefl(i);
		printf("Calculated vs. measured defl. at position %d = %f vs %f (%f)\n",
			i+1,d.calculated,d.measured,mz[i]);
		err += pow(d.measured-d.calculated,2);
	}
	err = sqrt(err/p.length());
	printf("\nRMSE in backcalculation = %f\n\n\n",err);
	goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
int
main()
{
	unsigned i;

double Ebad[17][5] = {
	{   16.322063,   39.174653,   29.353888,  230.667895,  649.966568}, 
	{   11.712960,   28.542286,   13.832604,   10.659610,  111.108490}, 
	{   13.264302,   33.157726,   35.584088,  116.874717,  222.735240},
	{   15.510243,   27.700517,   26.590475,   81.123447,  190.402050},
	{   28.488183,   47.349302,  188.564558, 2562.375905,  355.890892},
	{   38.212162,  759.059618,  939.373980, 3361.048713, 1695.809708},
	{ 7364.663840,  780.808133,  132.661665,   24.824413,   12.119889},
	{ 1329.603156,  710.290320,   44.946861,   15.770719,   32.118840}, 
	{ 5032.730142,  777.195194,  161.157751,   57.690076,   25.681444},
	{  762.106108, 1335.220981,   71.892747,   60.339601,  167.073355}, 
	{  876.798643, 5526.552481,  817.357172,  197.307914,   11.565737}, 
	{  666.900116, 8980.627448, 1851.246358,  233.998421, 1364.243244}, 
	{ 2559.676400,  201.936867,   10.316924,  104.343091,  134.180512}, 
	{ 4058.125273, 1923.655609,  177.831065,   15.968101,   31.230790}, 
	{ 4232.885833, 2459.175446,  153.045334,  279.037000,   14.379712}, 
	{ 1240.509997, 2480.000309,  203.603987,  806.573371,   76.893795}, 
	{  567.332552,  113.212841,  304.356563,   17.424377,   28.011755},
};

	LEbackcalc Bowl;
	cset<point3d> p;
	p.add(point3d(   0.0,175.0,0.0));
	p.add(point3d( 100.0,175.0,0.0));
	p.add(point3d( 200.0,175.0,0.0));
	p.add(point3d( 400.0,175.0,0.0));
	p.add(point3d( 600.0,175.0,0.0));
	p.add(point3d( 900.0,175.0,0.0));
	p.add(point3d(1200.0,175.0,0.0));
	//for (i = 0; i < 256; i++)
	//	p.add(point3d(   26.5*(i-100),175.0,  0.0));
	//p.add(point3d(   0.0,175.0,  0.0));
	//p.add(point3d(   0.0,175.0, 50.0));
	//p.add(point3d(   0.0,175.0,150.0));
	//p.add(point3d(   0.0,175.0,300.0));
	//p.add(point3d(   0.0,175.0,500.0));
	//p.add(point3d(   0.0,  0.0,  0.0));
	//p.add(point3d(   0.0,  0.0, 50.0));
	//p.add(point3d(   0.0,  0.0,150.0));
	//p.add(point3d(   0.0,  0.0,300.0));
	//p.add(point3d(   0.0,  0.0,500.0));

	//srand((unsigned)time(NULL));
	Bowl.addlayer( 50.0,800000.0,0.5);
	Bowl.addlayer(100.0,400000.0,0.5);
	Bowl.addlayer(150.0,200000.0,0.5);
	Bowl.addlayer(200.0,100000.0,0.5);
	Bowl.addlayer(  0.0, 50000.0,0.5);

	Bowl.addload(point2d(0.0,  0.0),20*1e6,520);
	Bowl.addload(point2d(0.0,350.0),20*1e6,520);

	double T[5], mz, err;
	unsigned bad = 0;

redo:
	for (i = 0; i < Bowl.layers(); i++) {
		//do {
		//	T[i] = pow(10,RAND(4.0,7.0));
		//} while (i > 0 && (T[i] > 20.0*T[i-1] || T[i] < T[i-1]/20.0));
		//printf("Layer %d: E=%0.0f\n",i+1,T[i]/1e3);
		//Bowl.layer(i).emod(T[i]);
		T[i] = Ebad[bad][i]*1000;
		printf("Layer %d: E=%0.0f\n",i+1,T[i]/1e3);
		Bowl.layer(i).emod(T[i]);
	}
	
	Bowl.removepoints();
	for (i = 0; i < p.length(); i++)
		Bowl.addpoint(p[i]);

	Bowl.calculate(LEsystem::disp);
	//Bowl.odemark(true);
	//Bowl.fastnum();

	Bowl.removedeflections();
	printf("M! ");
	for (i = 0; i < Bowl.layers(); i++)
		printf(" %0.4e",log10(Bowl.layer(i).emod()));
	for (i = 0; i < p.length(); i++) {
		mz = Bowl.result(p[i]).result(pavedata::deflct,pavedata::zz);
		//mz += 0.01*stdnormal_rnd();
		Bowl.adddefl(p[i],mz);
		printf(" %0.4e",mz);
	}
	printf("\n");
	Bowl.removepoints();
	for (i = 0; i < p.length(); i++) {
		const defldata & d = Bowl.getdefl(i);
		printf("Measured defl. at position %d = %f\n",i+1,d.measured);
	}
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(0.0);

	bool rv = Bowl.backcalc();
	printf("\nDONE!\n\n");
	
	for (i = 0, err = 0.0; i < Bowl.layers(); i++) {
		err += (Bowl.layer(i).emod()/T[i])-1;
		printf("E-modulus for layer %d = %f (%f) or %4.2f (%4.2f)\n",i+1,Bowl.layer(i).emod()/1000,T[i]/1000,log10(Bowl.layer(i).emod()),log10(T[i]));
	}
	if (rv) {
		if (err > 0.1)
			printf("\nBUT I FAILED!\n\n");
	} else {
		printf("\nFAILED TO CONVERGE!\n\n");
	}
	for (i = 0, err = 0.0; i < p.length(); i++) {
		const defldata & d = Bowl.getdefl(i);
		printf("Calculated vs. measured defl. at position %d = %f vs %f\n",
			i+1,d.calculated,d.measured);
		err += pow(d.measured-d.calculated,2);
	}
	err = sqrt(err/p.length());
	printf("\nRMSE in backcalculation = %f\n\n\n",err);
	if (++bad < 17)
		goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
int
main()
{
	unsigned line = 0;
	unsigned i;
	double err;

// PengCheng's Redbluff data
/*
#define LINES 63
#define DEFLS 7
double data[LINES][DEFLS+1] = {
 {0.1574, 0.1341, 0.1182, 0.0917, 0.0675, 0.0360, 0.0147, 520},
 {0.1592, 0.1304, 0.1132, 0.0884, 0.0678, 0.0410, 0.0150, 513},
 {0.1789, 0.1503, 0.1344, 0.1110, 0.0951, 0.0666, 0.0332, 522},
 {0.1892, 0.1654, 0.1576, 0.1422, 0.1280, 0.1019, 0.0558, 505},
 {0.2765, 0.2118, 0.1786, 0.1440, 0.1187, 0.0811, 0.0297, 505},
 {0.4457, 0.2948, 0.2179, 0.1307, 0.0835, 0.0415, 0.0215, 518},
 {0.2954, 0.2190, 0.1807, 0.1259, 0.0886, 0.0470, 0.0198, 519},
 {0.2843, 0.2245, 0.1716, 0.1261, 0.0910, 0.0466, 0.0193, 516},
 {0.3042, 0.2308, 0.1877, 0.1310, 0.0908, 0.0443, 0.0167, 506},
 {0.3257, 0.2753, 0.2170, 0.1055, 0.0794, 0.0460, 0.0210, 503},
 {0.1811, 0.1551, 0.1285, 0.0999, 0.0808, 0.0504, 0.0221, 491},
 {0.2876, 0.2262, 0.1845, 0.1365, 0.0985, 0.0570, 0.0231, 508},
 {0.2667, 0.2205, 0.1933, 0.1506, 0.1150, 0.0602, 0.0299, 506},
 {0.2200, 0.1814, 0.1615, 0.1317, 0.1059, 0.0712, 0.0421, 498},
 {0.2556, 0.2048, 0.1721, 0.1279, 0.1004, 0.0623, 0.0382, 497},
 {0.2551, 0.2091, 0.1787, 0.1310, 0.0983, 0.0575, 0.0314, 526},
 {0.2359, 0.1844, 0.1522, 0.1076, 0.0775, 0.0427, 0.0209, 507},
 {0.1786, 0.1564, 0.1427, 0.1151, 0.0901, 0.0521, 0.0215, 504},
 {0.2400, 0.2025, 0.1777, 0.1393, 0.1115, 0.0707, 0.0294, 492},
 {0.2135, 0.1754, 0.1488, 0.1109, 0.0864, 0.0524, 0.0251, 492},
 {0.1986, 0.1684, 0.1440, 0.1076, 0.0812, 0.0462, 0.0234, 488},
 {0.1984, 0.1689, 0.1491, 0.1158, 0.0871, 0.0472, 0.0197, 687},
 {0.1972, 0.1612, 0.1401, 0.1096, 0.0853, 0.0512, 0.0203, 681},
 {0.2244, 0.1888, 0.1690, 0.1402, 0.1197, 0.0852, 0.0420, 692},
 {0.2428, 0.2117, 0.2003, 0.1818, 0.1636, 0.1301, 0.0714, 672},
 {0.3397, 0.2596, 0.2214, 0.1813, 0.1512, 0.1064, 0.0415, 680},
 {0.5274, 0.3555, 0.2658, 0.1630, 0.1067, 0.0554, 0.0305, 692},
 {0.3572, 0.2669, 0.2216, 0.1554, 0.1114, 0.0615, 0.0254, 699},
 {0.3472, 0.2762, 0.2135, 0.1588, 0.1164, 0.0611, 0.0257, 692},
 {0.3743, 0.2853, 0.2339, 0.1658, 0.1161, 0.0591, 0.0212, 685},
 {0.4023, 0.3399, 0.2710, 0.1389, 0.1046, 0.0619, 0.0282, 679},
 {0.2367, 0.2030, 0.1676, 0.1315, 0.1060, 0.0675, 0.0295, 679},
 {0.3547, 0.2791, 0.2284, 0.1714, 0.1245, 0.0738, 0.0323, 676},
 {0.3328, 0.2748, 0.2431, 0.1913, 0.1440, 0.0792, 0.0405, 682},
 {0.2849, 0.2348, 0.2100, 0.1721, 0.1381, 0.0943, 0.0560, 668},
 {0.3170, 0.2552, 0.2160, 0.1619, 0.1276, 0.0816, 0.0492, 675},
 {0.3165, 0.2599, 0.2227, 0.1661, 0.1261, 0.0760, 0.0387, 708},
 {0.2938, 0.2299, 0.1897, 0.1357, 0.0993, 0.0554, 0.0283, 675},
 {0.2291, 0.2002, 0.1824, 0.1477, 0.1172, 0.0686, 0.0292, 675},
 {0.3116, 0.2631, 0.2302, 0.1825, 0.1468, 0.0944, 0.0380, 665},
 {0.2734, 0.2259, 0.1934, 0.1451, 0.1136, 0.0719, 0.0338, 669},
 {0.2560, 0.2168, 0.1855, 0.1401, 0.1063, 0.0619, 0.0291, 669},
 {0.2844, 0.2419, 0.2138, 0.1661, 0.1249, 0.0699, 0.0310, 1082},
 {0.2828, 0.2297, 0.1986, 0.1552, 0.1225, 0.0753, 0.0317, 1080},
 {0.3238, 0.2736, 0.2445, 0.2039, 0.1756, 0.1264, 0.0649, 1095},
 {0.3557, 0.3106, 0.2936, 0.2657, 0.2384, 0.1900, 0.1044, 1062},
 {0.4703, 0.3643, 0.3139, 0.2600, 0.2195, 0.1583, 0.0645, 1069},
 {0.7118, 0.4867, 0.3686, 0.2316, 0.1538, 0.0855, 0.0475, 1088},
 {0.4864, 0.3643, 0.3064, 0.2195, 0.1598, 0.0905, 0.0410, 1088},
 {0.4763, 0.3813, 0.3002, 0.2256, 0.1671, 0.0923, 0.0405, 1089},
 {0.5070, 0.3881, 0.3210, 0.2300, 0.1641, 0.0863, 0.0342, 1061},
 {0.5478, 0.4607, 0.3728, 0.2047, 0.1563, 0.0936, 0.0432, 1053},
 {0.3369, 0.2897, 0.2405, 0.1890, 0.1532, 0.0987, 0.0459, 1039},
 {0.4936, 0.3900, 0.3207, 0.2431, 0.1800, 0.1107, 0.0505, 1065},
 {0.4650, 0.3854, 0.3415, 0.2706, 0.2084, 0.1183, 0.0610, 1061},
 {0.4112, 0.3400, 0.3049, 0.2511, 0.2028, 0.1402, 0.0854, 1051},
 {0.4502, 0.3632, 0.3102, 0.2350, 0.1857, 0.1225, 0.0760, 1061},
 {0.4459, 0.3657, 0.3142, 0.2367, 0.1823, 0.1133, 0.0614, 1107},
 {0.4154, 0.3258, 0.2700, 0.1966, 0.1441, 0.0849, 0.0437, 1079},
 {0.3324, 0.2888, 0.2636, 0.2134, 0.1713, 0.1015, 0.0444, 1055},
 {0.4545, 0.3864, 0.3399, 0.2716, 0.2210, 0.1448, 0.0597, 1031},
 {0.3966, 0.3290, 0.2830, 0.2157, 0.1703, 0.1102, 0.0534, 1057},
 {0.3669, 0.3099, 0.2662, 0.2040, 0.1567, 0.0944, 0.0460, 1050},
};
*/

// CTMETWG Data set 1
#define LINES 63
#define DEFLS 8
double data[LINES][DEFLS+1] = {
{0.1293, 0.1204, 0.1123, 0.1019, 0.0904, 0.0734, 0.0612, 0.0498, 616.8}, 
{0.1278, 0.1201, 0.1074, 0.1003, 0.0904, 0.0732, 0.0622, 0.0498, 612.4}, 
{0.1295, 0.1191, 0.1092, 0.0988, 0.0899, 0.0739, 0.0602, 0.0500, 632.0}, 
{0.1019, 0.0973, 0.0899, 0.0838, 0.0808, 0.0699, 0.0605, 0.0538, 615.6}, 
{0.0975, 0.0914, 0.0879, 0.0792, 0.0803, 0.0627, 0.0579, 0.0483, 606.1}, 
{0.0993, 0.0899, 0.0846, 0.0798, 0.0759, 0.0655, 0.0577, 0.0460, 615.6}, 
{0.1417, 0.1207, 0.1168, 0.1087, 0.0988, 0.0861, 0.0668, 0.0554, 613.7}, 
{0.1420, 0.1240, 0.1184, 0.1090, 0.1001, 0.0848, 0.0671, 0.0559, 615.6}, 
{0.1455, 0.1255, 0.1168, 0.1087, 0.0998, 0.0828, 0.0635, 0.0541, 624.4}, 
{0.2670, 0.2352, 0.2154, 0.1864, 0.1671, 0.1313, 0.0940, 0.0737, 569.0}, 
{0.2680, 0.2403, 0.2149, 0.1915, 0.1689, 0.1300, 0.1064, 0.0752, 589.1}, 
{0.2629, 0.2316, 0.2118, 0.1859, 0.1651, 0.1273, 0.0914, 0.0742, 612.4}, 
{0.1834, 0.1615, 0.1461, 0.1293, 0.1207, 0.0914, 0.0734, 0.0612, 616.8}, 
{0.1859, 0.1623, 0.1506, 0.1298, 0.1207, 0.0942, 0.0780, 0.0688, 618.1}, 
{0.1854, 0.1636, 0.1483, 0.1308, 0.1196, 0.0935, 0.0739, 0.0620, 632.0}, 
{0.1910, 0.1742, 0.1613, 0.1392, 0.1318, 0.1107, 0.0925, 0.0798, 626.3}, 
{0.1908, 0.1737, 0.1626, 0.1417, 0.1280, 0.1102, 0.0892, 0.0826, 618.1}, 
{0.1819, 0.1628, 0.1532, 0.1351, 0.1270, 0.1057, 0.0897, 0.0747, 596.7}, 
{0.0986, 0.0871, 0.0795, 0.0813, 0.0838, 0.0744, 0.0665, 0.0549, 638.3}, 
{0.0978, 0.0866, 0.0843, 0.0818, 0.0815, 0.0742, 0.0660, 0.0526, 641.4}, 
{0.0942, 0.0848, 0.0782, 0.0782, 0.0805, 0.0724, 0.0658, 0.0569, 610.5}, 
{0.0935, 0.0663, 0.0622, 0.0592, 0.0561, 0.0508, 0.0432, 0.0386, 624.4}, 
{0.0927, 0.0673, 0.0607, 0.0584, 0.0551, 0.0470, 0.0437, 0.0404, 630.7}, 
{0.0912, 0.0665, 0.0625, 0.0602, 0.0554, 0.0485, 0.0439, 0.0378, 644.6}, 
{0.1128, 0.0958, 0.0881, 0.0815, 0.0775, 0.0663, 0.0564, 0.0498, 635.1}, 
{0.1125, 0.0963, 0.0897, 0.0836, 0.0780, 0.0676, 0.0561, 0.0490, 598.6}, 
{0.1097, 0.0942, 0.0892, 0.0808, 0.0757, 0.0645, 0.0561, 0.0480, 577.1}, 
{0.1857, 0.1488, 0.1339, 0.1168, 0.1128, 0.0902, 0.0752, 0.0635, 618.1}, 
{0.1859, 0.1567, 0.1346, 0.1181, 0.1087, 0.0902, 0.0765, 0.0630, 635.1}, 
{0.1814, 0.1483, 0.1326, 0.1166, 0.1064, 0.0904, 0.0757, 0.0638, 613.7}, 
{0.1694, 0.1494, 0.1359, 0.1212, 0.1105, 0.0871, 0.0719, 0.0582, 620.0}, 
{0.1669, 0.1466, 0.1331, 0.1219, 0.1113, 0.0886, 0.0729, 0.0569, 620.0}, 
{0.1651, 0.1445, 0.1318, 0.1196, 0.1092, 0.0869, 0.0709, 0.0574, 603.0}, 
{0.1722, 0.1544, 0.1415, 0.1270, 0.1199, 0.0996, 0.0836, 0.0704, 613.7}, 
{0.1732, 0.1516, 0.1420, 0.1303, 0.1219, 0.1008, 0.0833, 0.0701, 615.6}, 
{0.1725, 0.1501, 0.1412, 0.1295, 0.1219, 0.1029, 0.0848, 0.0711, 606.1}, 
{0.1864, 0.1689, 0.1580, 0.1499, 0.1392, 0.1189, 0.1003, 0.0838, 618.1}, 
{0.1892, 0.1702, 0.1600, 0.1506, 0.1435, 0.1207, 0.1006, 0.0843, 633.9}, 
{0.1842, 0.1659, 0.1552, 0.1461, 0.1359, 0.1161, 0.0975, 0.0818, 613.7}, 
{0.2370, 0.2154, 0.2002, 0.1811, 0.1631, 0.1323, 0.1085, 0.0889, 650.9}, 
{0.2344, 0.2139, 0.1979, 0.1781, 0.1613, 0.1303, 0.1059, 0.0856, 621.3}, 
{0.2344, 0.2129, 0.1976, 0.1788, 0.1628, 0.1316, 0.1074, 0.0876, 615.6}, 
{0.2174, 0.1938, 0.1783, 0.1565, 0.1392, 0.1095, 0.0833, 0.0724, 613.7}, 
{0.2164, 0.1938, 0.1770, 0.1562, 0.1384, 0.1087, 0.0828, 0.0683, 610.5}, 
{0.2129, 0.1923, 0.1753, 0.1539, 0.1364, 0.1077, 0.0826, 0.0658, 609.3}, 
{0.2880, 0.2438, 0.2174, 0.1897, 0.1646, 0.1316, 0.1019, 0.0798, 610.5}, 
{0.2875, 0.2421, 0.2177, 0.1910, 0.1707, 0.1331, 0.1034, 0.0813, 591.0}, 
{0.2868, 0.2441, 0.2177, 0.1908, 0.1692, 0.1318, 0.1016, 0.0808, 610.5}, 
{0.2522, 0.2062, 0.1788, 0.1593, 0.1382, 0.1052, 0.0820, 0.0673, 612.4}, 
{0.2471, 0.2024, 0.1770, 0.1549, 0.1374, 0.1064, 0.0833, 0.0683, 610.5}, 
{0.2454, 0.2009, 0.1778, 0.1532, 0.1369, 0.1054, 0.0828, 0.0676, 606.1}, 
{0.2151, 0.1808, 0.1588, 0.1425, 0.1262, 0.0986, 0.0782, 0.0625, 591.0}, 
{0.2129, 0.1788, 0.1598, 0.1415, 0.1265, 0.1003, 0.0780, 0.0627, 603.0}, 
{0.2101, 0.1765, 0.1567, 0.1402, 0.1247, 0.0988, 0.0782, 0.0612, 606.1}, 
{0.2484, 0.2057, 0.1806, 0.1605, 0.1430, 0.1118, 0.0889, 0.0681, 607.4}, 
{0.2471, 0.2019, 0.1796, 0.1595, 0.1407, 0.1115, 0.0859, 0.0673, 598.6}, 
{0.2477, 0.2050, 0.1791, 0.1585, 0.1420, 0.1118, 0.0869, 0.0726, 621.3}, 
{0.2205, 0.1913, 0.1758, 0.1567, 0.1402, 0.1097, 0.0853, 0.0640, 629.4}, 
{0.2146, 0.1887, 0.1715, 0.1557, 0.1374, 0.1080, 0.0836, 0.0650, 603.0}, 
{0.2156, 0.1900, 0.1709, 0.1575, 0.1359, 0.1105, 0.0866, 0.0630, 637.0}, 
{0.2550, 0.2134, 0.1902, 0.1626, 0.1430, 0.1102, 0.0800, 0.0658, 610.5}, 
{0.2527, 0.2103, 0.1857, 0.1608, 0.1412, 0.1087, 0.0818, 0.0643, 615.6}, 
{0.2553, 0.2149, 0.1875, 0.1626, 0.1433, 0.1090, 0.0836, 0.0655, 603.0}, 
};

	LEbackcalc Bowl;
	Bowl.setup(0.0001,0.0005,1e-4,6);

	// PengCheng's Redbluff data
	//Bowl.addlayer(177.8,0.0,0.35);
	//Bowl.addlayer(  0.0,0.0,0.35);

	// CTMETWG Data Set 1
	Bowl.addlayer(152.0,0.0,0.35);
	Bowl.addlayer(213.0,0.0,0.35);
	Bowl.addlayer(213.0,0.0,0.35);
	Bowl.addlayer(  0.0,0.0,0.35);

	cset<point3d> p;
	// PengCheng's Redbluff data
	//p.add(point3d(   0.0,0.0,  0.0));
	//p.add(point3d( 200.0,0.0,  0.0));
	//p.add(point3d( 300.0,0.0,  0.0));
	//p.add(point3d( 460.0,0.0,  0.0));
	//p.add(point3d( 610.0,0.0,  0.0));
	//p.add(point3d( 910.0,0.0,  0.0));
	//p.add(point3d(1525.0,0.0,  0.0));

	// CTMETWG Data Set 1
	p.add(point3d(   0.0,0.0,  0.0));
	p.add(point3d( 203.0,0.0,  0.0));
	p.add(point3d( 305.0,0.0,  0.0));
	p.add(point3d( 457.0,0.0,  0.0));
	p.add(point3d( 610.0,0.0,  0.0));
	p.add(point3d( 914.0,0.0,  0.0));
	p.add(point3d(1219.0,0.0,  0.0));
	p.add(point3d(1525.0,0.0,  0.0));

redo:
	Bowl.removeloads();
	Bowl.addload(point2d(0.0,  0.0),0.0,data[line][DEFLS],150.0);
	//for (double E1 = 4.0; E1 <= 7.0; E1 += 0.01) {
	//	for (double E2 = 4.0; E2 <= 6.0; E2 += 0.01) {
	//		Bowl.layer(0).emod(pow(10,E1));
	//		Bowl.layer(1).emod(pow(10,E2));
	//		Bowl.calculate(LEsystem::disp);
	//		printf("%0.4g\t%0.4g",E1,E2);
	//		double mz;
	//		for (i = 0; i < DEFLS; i++) {
	//			mz = Bowl.getdefl(i).result(pavedata::deflct,pavedata::zz);
	//			printf("\t%0.16e",mz);
	//		}
	//		printf("\n");
	//	}
	//}
	//exit(1);

	Bowl.removedeflections();
	for (i = 0; i < DEFLS; i++)
		Bowl.adddefl(p[i],data[line][i]);
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(0.0);

	bool rv = Bowl.backcalc();
	printf("Line %2i: %s",line+1,rv ? "DONE!   " : "FAILED! ");
	for (i = 0; i < Bowl.layers(); i++)
		printf(" %7.1f",Bowl.layer(i).emod()/1000);
	for (i = 0, err = 0.0; i < p.length(); i++) {
		const defldata & d = Bowl.getdefl(i);
		printf(" %6.1f",d.calculated*1000);
		err += pow(d.measured-d.calculated,2);
	}
	err = sqrt(err/p.length());
	printf("(err=%5.1f)\n",err*1000);
	while (++line < LINES)
		goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
int
main()
{
	unsigned i;
	LEbackcalc Bowl;

	Bowl.setup(0.0001,0.0,1e-4,5);

	Bowl.addlayer( 40.0,0.0,0.45);
	Bowl.addlayer(143.0,0.0,0.35);
	Bowl.addlayer(272.0,0.0,0.35);
	Bowl.addlayer(  0.0,0.0,0.35);

	Bowl.addload(point2d(0.0,0.0),0.0,999.9,150.0);
	
	Bowl.adddefl(point3d(   0.0,0.0,0.0), 0.2939);
	Bowl.adddefl(point3d( 203.0,0.0,0.0), 0.1831);
	Bowl.adddefl(point3d( 305.0,0.0,0.0), 0.1631);
	Bowl.adddefl(point3d( 457.0,0.0,0.0), 0.1463);
	Bowl.adddefl(point3d( 610.0,0.0,0.0), 0.1265);
	Bowl.adddefl(point3d( 915.0,0.0,0.0), 0.0914);
	Bowl.adddefl(point3d(1219.0,0.0,0.0), 0.0676);
	Bowl.adddefl(point3d(1524.0,0.0,0.0), 0.0516);

	for (i = 0; i < Bowl.deflections(); i++) {
		printf("Measured defl. at position %d = %f\n",i+1,Bowl.getdefl(i).measured);
	}

	if (Bowl.backcalc()) {
		printf("\nDONE!\n\n");
	} else {
		printf("\nFAILED TO CONVERGE!\n\n");
	}
	for (i = 0; i < Bowl.layers(); i++) {
		printf("E-modulus for layer %d = %f\n",i+1,Bowl.layer(i).emod()/1000);
	}
	double err = 0.0;
	for (i = 0; i < Bowl.deflections(); i++) {
		const defldata & d = Bowl.getdefl(i);
		printf("Calculated vs. measured defl. at position %d = %f vs %f\n",
			i+1,d.calculated,d.measured);
		err += pow(d.measured-d.calculated,2);
	}
	err = sqrt(err/i);
	printf("\nAbsolute error in backcalculation = %f\n",err);
	return 0;
}
#endif

#ifdef NOBUILD
// This is the special program for the paper...
int
main()
{
	unsigned i;
	double mz;
	LEbackcalc Bowl;
	static double T[2]= {800000.0, 50000.0};

	//srand((unsigned)time(NULL));
	Bowl.setup(0.0,0.0,1e-12,100);

	Bowl.addlayer(75.0,T[0],0.45);
	Bowl.addlayer( 0.0,T[1],0.35);

	Bowl.addload(point2d(0.0,  0.0),20*1e6,520);
	Bowl.addload(point2d(0.0,350.0),20*1e6,520);

	Bowl.addpoint(point3d(   0.0,175.0,  0.0));
	Bowl.addpoint(point3d( 100.0,175.0,  0.0));
	Bowl.addpoint(point3d( 200.0,175.0,  0.0));
	Bowl.addpoint(point3d( 400.0,175.0,  0.0));
	Bowl.addpoint(point3d( 600.0,175.0,  0.0));
	Bowl.addpoint(point3d( 900.0,175.0,  0.0));
	Bowl.addpoint(point3d(1200.0,175.0,  0.0));

	Bowl.calculate(LEsystem::disp);

	//printf("%8.6g\t%8.6g\t",log10(T[0]),log10(T[1]));
	for (i = 0; i < Bowl.results(); i++) {
		const point3d & p = Bowl.result(i);
		mz = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
		Bowl.adddefl(p,mz);
		//printf("%0.16e\t",mz);
	}
	//printf("%0.4e\n",0.0);

	//for (double E1 = 3.0; E1 <= 8.0; E1 += (E1<3.999||E1>7.999?0.01:0.005)) {
	//	for (double E2 = 3.0; E2 <= 6.0; E2 += (E2<3.999||E2>4.999?0.01:0.002)) {
	//		Bowl.layer(0).emod(pow(10,E1));
	//		Bowl.layer(1).emod(pow(10,E2));
	//		Bowl.calculate(LEsystem::disp);
	//		printf("%0.4g\t%0.4g",E1,E2);
	//		double err = 0.0;
	//		for (i = 0; i < Bowl.deflections(); i++) {
	//			point3d & p = Bowl.getdefl(i);
	//			mz = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
	//			printf("\t%0.16e",mz);
	//			err += pow(mz-Bowl.defl[i].measured,2);
	//		}
	//		printf("\t%0.16e\n",sqrt(err/i));
	//	}
	//}
	static double S1[2] = {3, 3};
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(pow(10,S1[i]));
	Bowl.backcalc();
	static double S2[2] = {8, 6};
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(pow(10,S2[i]));
	Bowl.backcalc();
	static double S3[2] = {0, 0};
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(pow(10,S3[i]));
	Bowl.backcalc();
	
	static double X[2] = {7.2260901743325290, 4.2287893718141945};
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(pow(10,X[i]));
	//Bowl.calculate(LEsystem::dispgrad);

	//printf("%8.6g\t%8.6g\t",log10(T[0]),log10(T[1]));
	//for (i = 0; i < Bowl.deflections(); i++) {
	//	point3d & p = Bowl.getdefl(i);
	//	mz = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
	//	//Bowl.adddefl(p,mz);
	//	printf("%0.16e\n",mz);
	//	for (unsigned j = 0; j < Bowl.layers(); j++) {
	//		mz = Bowl.result(p).deflgrad[j];
	//		printf("%0.16e\t",mz);
	//	}
	//	printf("\n");
	//}
	//printf("%0.4e\n",0.0);

	Bowl.backcalc();
	//XXX: Should have a conjgrad here...
	//XXX: Should have a gauss-newton here...
	
	return 0;
}
#endif

#ifdef NOBUILD
// This is the second special program for the paper...
int
main()
{
	unsigned i;
	double mz[7], t;
	LEbackcalc Bowl;
	static double T[2]= {800000.0, 50000.0};

	Bowl.setup(0.0001,0.0005,1e-6,5);

	Bowl.addlayer(75.0,T[0],0.45);
	Bowl.addlayer( 0.0,T[1],0.35);

	Bowl.addload(point2d(0.0,  0.0),20*1e6,520);
	Bowl.addload(point2d(0.0,350.0),20*1e6,520);

	Bowl.addpoint(point3d(   0.0,175.0,  0.0));
	Bowl.addpoint(point3d( 100.0,175.0,  0.0));
	Bowl.addpoint(point3d( 200.0,175.0,  0.0));
	Bowl.addpoint(point3d( 400.0,175.0,  0.0));
	Bowl.addpoint(point3d( 600.0,175.0,  0.0));
	Bowl.addpoint(point3d( 900.0,175.0,  0.0));
	Bowl.addpoint(point3d(1200.0,175.0,  0.0));

	Bowl.calculate(LEsystem::disp);
	for (i = 0; i < Bowl.results(); i++) {
		mz[i] = Bowl.result(i).result(pavedata::deflct,pavedata::zz);
	}

redo:
	Bowl.removedeflections();
	for (i = 0; i < Bowl.results(); i++) {
		const point3d & p = Bowl.result(i);
		t = mz[i] + RAND(2,6)*0.0001*stdnormal_rnd();
		t = round(t/0.0001)*0.0001;
		Bowl.adddefl(p,t);
	}

	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(T[i]);
	Bowl.backcalc();
	
	printf("%12.3f\t",Bowl.layer(0).emod()/1000);
	printf("%12.3f\n",Bowl.layer(1).emod()/1000);
	goto redo;
	return 0;
}
#endif
