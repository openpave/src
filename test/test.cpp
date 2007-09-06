/**************************************************************************

	TEST.CPP - A test harness for OpenPave.org code.

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

#define _EVENT_IMP
#define _PROGRESS_IMP
#include "pavement.h"
#include "traffic.h"
#include "reliability.h"
#include "matrix.h"
#include <stdlib.h>
#include <time.h>

#ifdef NOBUILD
int
main()
{
	WIMsurvey WIM(str2time("000101",NULL),str2time("040101",NULL));

	WIM.ProcessRSADir("C:/Work/SANRAL LEF/NEURAL/2030","C:/Work/SANRAL LEF/NEURAL/2030.WIM");
}
#endif

#ifdef NOBUILD
int
main(int argc, char* argv[])
{
	int i, l;
	LEbackcalc Bowl;

	//srand((unsigned)time(NULL));
	double T[5];
	double h[5], v[5], mz[7], t, f;

redo:
	l = (int)floor(RAND(2,6));
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
	Bowl.addpoint(point3d(   0.0,0.0,  0.0));
	Bowl.addpoint(point3d( 200.0,0.0,  0.0));
	Bowl.addpoint(point3d( 300.0,0.0,  0.0));
	Bowl.addpoint(point3d( 460.0,0.0,  0.0));
	Bowl.addpoint(point3d( 610.0,0.0,  0.0));
	Bowl.addpoint(point3d( 910.0,0.0,  0.0));
	Bowl.addpoint(point3d(1525.0,0.0,  0.0));

	Bowl.calculate(LEsystem::disp);
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
		mz[i] = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
		printf("Measured defl. at position %d = %f\n",i+1,mz[i]);
		if (mz[i] < 0.0) {
			printf("Not a nice pavement! Bailing...\n");
			goto redo;
		}
	}

	Bowl.callcount = 0;
	Bowl.removedeflections();
	printf("M! ");
	for (i = 0; i < Bowl.layers(); i++)
		printf(" %0.4g",log10(Bowl.layer(i).emod()));
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
		f = t = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
		t += RAND(2,5)*0.0001*random_normal();
		t = ROUND(t/0.0001)*0.0001;
		Bowl.adddefl(p,t);
		printf(" %0.4g",t);
	}
	printf("\n");
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(0.0);

	bool rv = Bowl.backcalc();
	printf("\n%s! (with %i calls for %i layers)\n",(rv?"DONE":"FAILED"),Bowl.callcount,l);

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
	printf("\nAbsolute error in backcalculation = %f\n",Bowl.bowlerror(i,0,0,0));
	for (i = 0; i < Bowl.defl.length(); i++) {
		printf("Calculated vs. measured defl. at position %d = %f vs %f (%f)\n",
			i+1,Bowl.defl[i].calculated,Bowl.defl[i].measured,mz[i]);
	}
	goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
int
main(int argc, char* argv[])
{
	int i;

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
}

	LEbackcalc Bowl;

	srand((unsigned)time(NULL));
	Bowl.addlayer( 50.0,800000.0,0.5);
	Bowl.addlayer(100.0,400000.0,0.5);
	Bowl.addlayer(150.0,200000.0,0.5);
	Bowl.addlayer(200.0,100000.0,0.5);
	Bowl.addlayer(  0.0, 50000.0,0.5);

	Bowl.addload(point2d(0.0,  0.0),20*1e6,520);
	Bowl.addload(point2d(0.0,350.0),20*1e6,520);

	double T[5];
	double mz;
	int bad = 10;

redo:
	for (i = 0; i < Bowl.layers(); i++) {
		T[i] = Bowl.layer(i).emod();
		//do {
		//	T[i] = pow(10,4.0+3.0*(double)(rand())/(double)(RAND_MAX));
		//} while (i > 0 && (T[i] > 20.0*T[i-1] || T[i] < T[i-1]/20.0));
		//printf("Layer %d: E=%0.0f\n",i+1,T[i]/1e3);
		//Bowl.layer(i).emod(T[i]);
		T[i] = Ebad[bad][i]*1000;
		printf("Layer %d: E=%0.0f\n",i+1,T[i]/1e3);
		Bowl.layer(i).emod(T[i]);
	}

	//for (i = 0; i < 256; i++)
	//	Bowl.addpoint(point3d(   26.5*(i-100),175.0,  0.0));
	Bowl.addpoint(point3d(   0.0,175.0,  0.0));
	Bowl.addpoint(point3d( 100.0,175.0,  0.0));
	Bowl.addpoint(point3d( 200.0,175.0,  0.0));
	Bowl.addpoint(point3d( 400.0,175.0,  0.0));
	Bowl.addpoint(point3d( 600.0,175.0,  0.0));
	Bowl.addpoint(point3d( 900.0,175.0,  0.0));
	Bowl.addpoint(point3d(1200.0,175.0,  0.0));
	//Bowl.addpoint(point3d(   0.0,175.0,  0.0));
	//Bowl.addpoint(point3d(   0.0,175.0, 50.0));
	//Bowl.addpoint(point3d(   0.0,175.0,150.0));
	//Bowl.addpoint(point3d(   0.0,175.0,300.0));
	//Bowl.addpoint(point3d(   0.0,175.0,500.0));
	//Bowl.addpoint(point3d(   0.0,  0.0,  0.0));
	//Bowl.addpoint(point3d(   0.0,  0.0, 50.0));
	//Bowl.addpoint(point3d(   0.0,  0.0,150.0));
	//Bowl.addpoint(point3d(   0.0,  0.0,300.0));
	//Bowl.addpoint(point3d(   0.0,  0.0,500.0));

	Bowl.calculate(LEsystem::disp);
	//Bowl.odemark(true);
	//Bowl.fastnum();

	Bowl.removedeflections();
	printf("M! ");
	for (i = 0; i < Bowl.layers(); i++)
		printf(" %0.4e",log10(Bowl.layer(i).emod()));
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
		mz = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
		//mz += 0.01*random_normal();
		Bowl.adddefl(p,mz);
		printf(" %0.4e",mz);
	}
	printf("\n");
	Bowl.removepoints();
	for (i = 0; i < Bowl.defl.length(); i++) {
		printf("Measured defl. at position %d = %f\n",i+1,Bowl.defl[i].measured);
	}
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(0.0);

	bool rv = Bowl.backcalc();
	printf("\nDONE!\n\n");
	
	for (i = 0, mz = 0.0; i < Bowl.layers(); i++) {
		mz += (Bowl.layer(i).emod()/T[i])-1;
		printf("E-modulus for layer %d = %f (%f)\n",i+1,Bowl.layer(i).emod()/1000,T[i]/1000);
	}
	if (rv) {
		if (mz > 0.1)
			printf("\nBUT I FAILED!\n\n");
	} else {
		printf("\nFAILED TO CONVERGE!\n\n");
	}
	printf("\nAbsolute error in backcalculation = %f\n",Bowl.bowlerror(i,0,0,0));
	for (i = 0; i < Bowl.defl.length(); i++) {
		printf("Calculated vs. measured defl. at position %d = %f vs %f\n",
			i+1,Bowl.defl[i].calculated,Bowl.defl[i].measured);
	}
	if (++bad < 17)
		goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
int
main(int argc, char* argv[])
{
	int line = 0;
	int i;

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


redo:
	// PengCheng's Redbluff data
	//Bowl.addpoint(point3d(   0.0,0.0,  0.0));
	//Bowl.addpoint(point3d( 200.0,0.0,  0.0));
	//Bowl.addpoint(point3d( 300.0,0.0,  0.0));
	//Bowl.addpoint(point3d( 460.0,0.0,  0.0));
	//Bowl.addpoint(point3d( 610.0,0.0,  0.0));
	//Bowl.addpoint(point3d( 910.0,0.0,  0.0));
	//Bowl.addpoint(point3d(1525.0,0.0,  0.0));

	// CTMETWG Data Set 1
	Bowl.addpoint(point3d(   0.0,0.0,  0.0));
	Bowl.addpoint(point3d( 203.0,0.0,  0.0));
	Bowl.addpoint(point3d( 305.0,0.0,  0.0));
	Bowl.addpoint(point3d( 457.0,0.0,  0.0));
	Bowl.addpoint(point3d( 610.0,0.0,  0.0));
	Bowl.addpoint(point3d( 914.0,0.0,  0.0));
	Bowl.addpoint(point3d(1219.0,0.0,  0.0));
	Bowl.addpoint(point3d(1525.0,0.0,  0.0));

	Bowl.removeloads();
	Bowl.addload(point2d(0.0,  0.0),0.0,data[line][DEFLS],150.0);
	//for (double E1 = 4.0; E1 <= 7.0; E1 += 0.01) {
	//	for (double E2 = 4.0; E2 <= 6.0; E2 += 0.01) {
	//		Bowl.layer(0).emod(pow(10,E1));
	//		Bowl.layer(1).emod(pow(10,E2));
	//		Bowl.calculate(LEsystem::disp);
	//		printf("%0.4g\t%0.4g",E1,E2);
	//		double mz;
	//		for (i = 0; i < Bowl.data.length(); i++) {
	//			mz = Bowl.data[i].result(pavedata::deflct,pavedata::zz);
	//			printf("\t%0.16e",mz);
	//		}
	//		printf("\n");
	//	}
	//}
	//exit(1);

	Bowl.removedeflections();
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
		Bowl.adddefl(p,data[line][i]);
	}
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(0.0);

	bool rv = Bowl.backcalc();
	printf("Line %2i: ", line+1);
	if (rv) {
		printf("DONE!   ");
	} else {
		printf("FAILED! ");
	}
	printf("(err=%5.1f) ",Bowl.bowlerror(i,0,0,0)*1000);
	for (i = 0; i < Bowl.layers(); i++)
		printf(" %7.1f",Bowl.layer(i).emod()/1000);
	for (i = 0; i < Bowl.defl.length(); i++)
		printf(" %6.1f",Bowl.defl[i].calculated*1000);
	printf("\n");
	while (++line < LINES)
		goto redo;
	return 0;
}
#endif

#ifdef NOBUILD
int
main(int argc, char* argv[])
{
	int i;
	LEbackcalc Bowl;

	Bowl.setup(0.0001,1e-4,5);

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

	for (i = 0; i < Bowl.defl.length(); i++) {
		printf("Measured defl. at position %d = %f\n",i+1,Bowl.defl[i].measured);
	}

	if (Bowl.backcalc()) {
		printf("\nDONE!\n\n");
	} else {
		printf("\nFAILED TO CONVERGE!\n\n");
	}
	printf("\nAbsolute error in backcalculation = %f\n",Bowl.bowlerror(i,0,0,0));
	for (i = 0; i < Bowl.layers(); i++) {
		printf("E-modulus for layer %d = %f\n",i+1,Bowl.layer(i).emod()/1000);
	}
	for (i = 0; i < Bowl.defl.length(); i++) {
		printf("Calculated vs. measured defl. at position %d = %f vs %f\n",
			i+1,Bowl.defl[i].calculated,Bowl.defl[i].measured);
	}
	return 0;
}
#endif

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
	Best.accurate();
	Slow.calculate(LEsystem::all);
	Fast.fastnum();

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

#ifdef NOBUILD
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
// This is the special program for the paper...
int
main()
{
	int i;
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
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
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
	//		for (i = 0; i < Bowl.defl.length(); i++) {
	//			point3d & p = Bowl.defl[i];
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
	Bowl.seed(Bowl.layers(),S3);
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(pow(10,S3[i]));
	Bowl.backcalc();
	
	static double X[2] = {7.2260901743325290, 4.2287893718141945};
	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(pow(10,X[i]));
	//Bowl.calculate(LEsystem::dispgrad);

	//printf("%8.6g\t%8.6g\t",log10(T[0]),log10(T[1]));
	//for (i = 0; i < Bowl.data.length(); i++) {
	//	point3d & p = Bowl.data[i];
	//	mz = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
	//	//Bowl.adddefl(p,mz);
	//	printf("%0.16e\n",mz);
	//	for (int j = 0; j < Bowl.layers(); j++) {
	//		mz = Bowl.result(p).deflgrad[j];
	//		printf("%0.16e\t",mz);
	//	}
	//	printf("\n");
	//}
	//printf("%0.4e\n",0.0);

	Bowl.backcalc();
	//XXX: Should have a conjgrad here...
	//XXX: Should have a gauss-newton here...
	
	return;
}
#endif

#ifdef NOBUILD
// This is the second special program for the paper...
int
main()
{
	int i;
	double mz[7], E[2], t;
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
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
		mz[i] = Bowl.result(p).result(pavedata::deflct,pavedata::zz);
	}

redo:
	Bowl.removedeflections();
	for (i = 0; i < Bowl.data.length(); i++) {
		point3d & p = Bowl.data[i];
		t = mz[i] + RAND(2,6)*0.0001*random_normal();
		t = ROUND(t/0.0001)*0.0001;
		Bowl.adddefl(p,t);
	}

	for (i = 0; i < Bowl.layers(); i++)
		Bowl.layer(i).emod(T[i]);
	Bowl.backcalc();
	
	printf("%12.3f\t",Bowl.layer(0).emod()/1000);
	printf("%12.3f\n",Bowl.layer(1).emod()/1000);
	goto redo;
	return;
}
#endif

#ifdef NOBUILD
/* Bessel J_0(x) function in double precision */
double dbesj0(double x)
{
    int k;
    double w, t, y, v, theta;
    static double a[8] = {
        -2.3655394e-12, 4.70889868e-10, 
        -6.78167892231e-8, 6.7816840038636e-6, 
        -4.340277777716935e-4, 0.0156249999999992397, 
        -0.2499999999999999638, 0.9999999999999999997
    };
    static double b[65] = {
        6.26681117e-11, -2.2270614428e-9, 
        6.62981656302e-8, -1.6268486502196e-6, 
        3.21978384111685e-5, -5.00523773331583e-4, 
        0.0059060313537449816, -0.0505265323740109701, 
        0.2936432097610503985, -1.0482565081091638637, 
        1.9181123286040428113, -1.13191994752217001, 
        -0.1965480952704682, 
        4.57457332e-11, -1.5814772025e-9, 
        4.55487446311e-8, -1.0735201286233e-6, 
        2.02015179970014e-5, -2.942392368203808e-4, 
        0.0031801987726150648, -0.0239875209742846362, 
        0.1141447698973777641, -0.2766726722823530233, 
        0.1088620480970941648, 0.5136514645381999197, 
        -0.2100594022073706033, 
        3.31366618e-11, -1.1119090229e-9, 
        3.08823040363e-8, -6.956602653104e-7, 
        1.23499947481762e-5, -1.66295194539618e-4, 
        0.0016048663165678412, -0.0100785479932760966, 
        0.0328996815223415274, -0.0056168761733860688, 
        -0.2341096400274429386, 0.2551729256776404262, 
        0.2288438186148935667, 
        2.38007203e-11, -7.731046439e-10, 
        2.06237001152e-8, -4.412291442285e-7, 
        7.3107766249655e-6, -8.91749801028666e-5, 
        7.34165451384135e-4, -0.0033303085445352071, 
        0.0015425853045205717, 0.0521100583113136379, 
        -0.1334447768979217815, -0.1401330292364750968, 
        0.2685616168804818919, 
        1.6935595e-11, -5.308092192e-10, 
        1.35323005576e-8, -2.726650587978e-7, 
        4.151324014176e-6, -4.43353052220157e-5, 
        2.815740758993879e-4, -4.393235121629007e-4, 
        -0.0067573531105799347, 0.0369141914660130814, 
        0.0081673361942996237, -0.257338128589888186, 
        0.0459580257102978932
    };
    static double c[70] = {
        -3.009451757e-11, -1.4958003844e-10, 
        5.06854544776e-9, 1.863564222012e-8, 
        -6.0304249068078e-7, -1.47686259937403e-6, 
        4.714331342682714e-5, 6.286305481740818e-5, 
        -0.00214137170594124344, -8.9157336676889788e-4, 
        0.04508258728666024989, -0.00490362805828762224, 
        -0.27312196367405374426, 0.04193925184293450356, 
        -7.1245356e-12, -4.1170814825e-10, 
        1.38012624364e-9, 5.704447670683e-8, 
        -1.9026363528842e-7, -5.33925032409729e-6, 
        1.736064885538091e-5, 3.0692619152608375e-4, 
        -9.2598938200644367e-4, -0.00917934265960017663, 
        0.02287952522866389076, 0.10545197546252853195, 
        -0.16126443075752985095, -0.19392874768742235538, 
        2.128344556e-11, -3.1053910272e-10, 
        -3.34979293158e-9, 4.50723289505e-8, 
        3.6437959146427e-7, -4.46421436266678e-6, 
        -2.523429344576552e-5, 2.7519882931758163e-4, 
        9.7185076358599358e-4, -0.00898326746345390692, 
        -0.01665959196063987584, 0.11456933464891967814, 
        0.07885001422733148815, -0.23664819446234712621, 
        3.035295055e-11, 5.486066835e-11, 
        -5.01026824811e-9, -5.0124684786e-9, 
        5.8012340163034e-7, 1.6788922416169e-7, 
        -4.373270270147275e-5, 1.183898532719802e-5, 
        0.00189863342862291449, -0.0011375924956163613, 
        -0.03846797195329871681, 0.02389746880951420335, 
        0.22837862066532347461, -0.06765394811166522844, 
        1.279875977e-11, 3.5925958103e-10, 
        -2.28037105967e-9, -4.852770517176e-8, 
        2.8696428000189e-7, 4.40131125178642e-6, 
        -2.366617753349105e-5, -2.4412456252884129e-4, 
        0.00113028178539430542, 0.0070847051391978908, 
        -0.02526914792327618386, -0.08006137953480093426, 
        0.16548380461475971846, 0.14688405470042110229
    };
    static double d[52] = {
        1.059601355592185731e-14, -2.71150591218550377e-13, 
        8.6514809056201638e-12, -4.6264028554286627e-10, 
        5.0815403835647104e-8, -1.76722552048141208e-5, 
        0.16286750396763997378, 2.949651820598278873e-13, 
        -8.818215611676125741e-12, 3.571119876162253451e-10, 
        -2.63192412099371706e-8, 4.709502795656698909e-6, 
        -0.005208333333333283282, 
        7.18344107717531977e-15, -2.51623725588410308e-13, 
        8.6017784918920604e-12, -4.6256876614290359e-10, 
        5.0815343220437937e-8, -1.7672255176494197e-5, 
        0.16286750396763433767, 2.2327570859680094777e-13, 
        -8.464594853517051292e-12, 3.563766464349055183e-10, 
        -2.631843986737892965e-8, 4.70950234228865941e-6, 
        -0.0052083333332278466225, 
        5.15413392842889366e-15, -2.27740238380640162e-13, 
        8.4827767197609014e-12, -4.6224753682737618e-10, 
        5.0814848128929134e-8, -1.7672254763876748e-5, 
        0.16286750396748926663, 1.7316195320192170887e-13, 
        -7.971122772293919646e-12, 3.544039469911895749e-10, 
        -2.631443902081701081e-8, 4.709498228695400603e-6, 
        -0.005208333331514365361, 
        3.84653681453798517e-15, -2.04464520778789011e-13, 
        8.3089298605177838e-12, -4.6155016158412096e-10, 
        5.081326369646665e-8, -1.76722528311426167e-5, 
        0.1628675039665006593, 1.3797879972460878797e-13, 
        -7.448089381011684812e-12, 3.51273379710695978e-10, 
        -2.630500895563592722e-8, 4.709483934775839193e-6, 
        -0.0052083333227940760113
    };

    w = fabs(x);
    if (w < 1) {
        t = w * w;
        y = ((((((a[0] * t + a[1]) * t + 
            a[2]) * t + a[3]) * t + a[4]) * t + 
            a[5]) * t + a[6]) * t + a[7];
    } else if (w < 8.5) {
        t = w * w * 0.0625;
        k = (int) t;
        t -= k + 0.5;
        k *= 13;
        y = (((((((((((b[k] * t + b[k + 1]) * t + 
            b[k + 2]) * t + b[k + 3]) * t + b[k + 4]) * t + 
            b[k + 5]) * t + b[k + 6]) * t + b[k + 7]) * t + 
            b[k + 8]) * t + b[k + 9]) * t + b[k + 10]) * t + 
            b[k + 11]) * t + b[k + 12];
    } else if (w < 12.5) {
        k = (int) w;
        t = w - (k + 0.5);
        k = 14 * (k - 8);
        y = ((((((((((((c[k] * t + c[k + 1]) * t + 
            c[k + 2]) * t + c[k + 3]) * t + c[k + 4]) * t + 
            c[k + 5]) * t + c[k + 6]) * t + c[k + 7]) * t + 
            c[k + 8]) * t + c[k + 9]) * t + c[k + 10]) * t + 
            c[k + 11]) * t + c[k + 12]) * t + c[k + 13];
    } else {
        v = 24 / w;
        t = v * v;
        k = 13 * ((int) t);
        y = ((((((d[k] * t + d[k + 1]) * t + 
            d[k + 2]) * t + d[k + 3]) * t + d[k + 4]) * t + 
            d[k + 5]) * t + d[k + 6]) * sqrt(v);
        theta = (((((d[k + 7] * t + d[k + 8]) * t + 
            d[k + 9]) * t + d[k + 10]) * t + d[k + 11]) * t + 
            d[k + 12]) * v - 0.78539816339744830962;
        y *= cos(w + theta);
    }
    return y;
}

/* Bessel J_1(x) function in double precision */
double dbesj1(double x)
{
    int k;
    double w, t, y, v, theta;
    static double a[8] = {
        -1.4810349e-13, 3.363594618e-11, 
        -5.65140051697e-9, 6.7816840144764e-7, 
        -5.425347222188379e-5, 0.00260416666666662438, 
        -0.06249999999999999799, 0.49999999999999999998
    };
    static double b[65] = {
        2.43721316e-12, -9.400554763e-11, 
        3.0605338998e-9, -8.287270492518e-8, 
        1.83020515991344e-6, -3.219783841164382e-5, 
        4.3795830161515318e-4, -0.00442952351530868999, 
        0.03157908273375945955, -0.14682160488052520107, 
        0.39309619054093640008, -0.4795280821510107028, 
        0.1414899934402712514, 
        1.82119257e-12, -6.862117678e-11, 
        2.1732790836e-9, -5.69359291782e-8, 
        1.20771046483277e-6, -2.020151799736374e-5, 
        2.5745933218048448e-4, -0.00238514907946126334, 
        0.01499220060892984289, -0.05707238494868888345, 
        0.10375225210588234727, -0.02721551202427354117, 
        -0.06420643306727498985, 
        1.352611196e-12, -4.9706947875e-11, 
        1.527944986332e-9, -3.8602878823401e-8, 
        7.82618036237845e-7, -1.23499947484511e-5, 
        1.45508295194426686e-4, -0.001203649737425854162, 
        0.006299092495799005109, -0.016449840761170764763, 
        0.002106328565019748701, 0.05852741000686073465, 
        -0.031896615709705053191, 
        9.97982124e-13, -3.5702556073e-11, 
        1.062332772617e-9, -2.5779624221725e-8, 
        4.96382962683556e-7, -7.310776625173004e-6, 
        7.8028107569541842e-5, -5.50624088538081113e-4, 
        0.002081442840335570371, -7.71292652260286633e-4, 
        -0.019541271866742634199, 0.033361194224480445382, 
        0.017516628654559387164, 
        7.31050661e-13, -2.5404499912e-11, 
        7.29360079088e-10, -1.6915375004937e-8, 
        3.06748319652546e-7, -4.151324014331739e-6, 
        3.8793392054271497e-5, -2.11180556924525773e-4, 
        2.74577195102593786e-4, 0.003378676555289966782, 
        -0.013842821799754920148, -0.002041834048574905921, 
        0.032167266073736023299
    };
    static double c[70] = {
        -1.185964494e-11, 3.9110295657e-10, 
        1.80385519493e-9, -5.575391345723e-8, 
        -1.8635897017174e-7, 5.42738239401869e-6, 
        1.181490114244279e-5, -3.300031939852107e-4, 
        -3.7717832892725053e-4, 0.01070685852970608288, 
        0.00356629346707622489, -0.13524776185998074716, 
        0.00980725611657523952, 0.27312196367405374425, 
        -3.029591097e-11, 9.259293559e-11, 
        4.96321971223e-9, -1.518137078639e-8, 
        -5.7045127595547e-7, 1.71237271302072e-6, 
        4.271400348035384e-5, -1.2152454198713258e-4, 
        -0.00184155714921474963, 0.00462994691003219055, 
        0.03671737063840232452, -0.06863857568599167175, 
        -0.21090395092505707655, 0.16126443075752985095, 
        -2.19760208e-11, -2.7659100729e-10, 
        3.74295124827e-9, 3.684765777023e-8, 
        -4.5072801091574e-7, -3.27941630669276e-6, 
        3.5713715545163e-5, 1.7664005411843533e-4, 
        -0.00165119297594774104, -0.00485925381792986774, 
        0.03593306985381680131, 0.04997877588191962563, 
        -0.22913866929783936544, -0.07885001422733148814, 
        5.16292316e-12, -3.9445956763e-10, 
        -6.6220021263e-10, 5.511286218639e-8, 
        5.01257940078e-8, -5.22111059203425e-6, 
        -1.34311394455105e-6, 3.0612891890766805e-4, 
        -7.103391195326182e-5, -0.00949316714311443491, 
        0.00455036998246516948, 0.11540391585989614784, 
        -0.04779493761902840455, -0.2283786206653234746, 
        2.697817493e-11, -1.6633326949e-10, 
        -4.3313486035e-9, 2.508404686362e-8, 
        4.8528284780984e-7, -2.58267851112118e-6, 
        -3.521049080466759e-5, 1.6566324273339952e-4, 
        0.00146474737522491617, -0.00565140892697147306, 
        -0.028338820556793004, 0.07580744376982855057, 
        0.16012275906960187978, -0.16548380461475971845
    };
    static double d[52] = {
        -1.272346002224188092e-14, 3.370464692346669075e-13, 
        -1.144940314335484869e-11, 6.863141561083429745e-10, 
        -9.491933932960924159e-8, 5.301676561445687562e-5, 
        0.162867503967639974, -3.652982212914147794e-13, 
        1.151126750560028914e-11, -5.165585095674343486e-10, 
        4.657991250060549892e-8, -1.186794704692706504e-5, 
        0.01562499999999994026, 
        -8.713069680903981555e-15, 3.140780373478474935e-13, 
        -1.139089186076256597e-11, 6.862299023338785566e-10, 
        -9.491926788274594674e-8, 5.301676558106268323e-5, 
        0.162867503967646622, -2.792555727162752006e-13, 
        1.108650207651756807e-11, -5.156745588549830981e-10, 
        4.657894859077370979e-8, -1.186794650130550256e-5, 
        0.01562499999987299901, 
        -6.304859171204770696e-15, 2.857249044208791652e-13, 
        -1.124956921556753188e-11, 6.858482894906716661e-10, 
        -9.49186795351689846e-8, 5.301676509057781574e-5, 
        0.1628675039678191167, -2.185193490132496053e-13, 
        1.048820673697426074e-11, -5.132819367467680132e-10, 
        4.65740943737299422e-8, -1.186794150862988921e-5, 
        0.01562499999779270706, 
        -4.74041720979200985e-15, 2.578715253644144182e-13, 
        -1.104148898414138857e-11, 6.850134201626289183e-10, 
        -9.49167823417491964e-8, 5.301676277588728159e-5, 
        0.1628675039690033136, -1.75512205749384229e-13, 
        9.848723331445182397e-12, -5.094535425482245697e-10, 
        4.656255982268609304e-8, -1.186792402114394891e-5, 
        0.01562499998712198636
    };

    w = fabs(x);
    if (w < 1) {
        t = w * w;
        y = (((((((a[0] * t + a[1]) * t + 
            a[2]) * t + a[3]) * t + a[4]) * t + 
            a[5]) * t + a[6]) * t + a[7]) * w;
    } else if (w < 8.5) {
        t = w * w * 0.0625;
        k = (int) t;
        t -= k + 0.5;
        k *= 13;
        y = ((((((((((((b[k] * t + b[k + 1]) * t + 
            b[k + 2]) * t + b[k + 3]) * t + b[k + 4]) * t + 
            b[k + 5]) * t + b[k + 6]) * t + b[k + 7]) * t + 
            b[k + 8]) * t + b[k + 9]) * t + b[k + 10]) * t + 
            b[k + 11]) * t + b[k + 12]) * w;
    } else if (w < 12.5) {
        k = (int) w;
        t = w - (k + 0.5);
        k = 14 * (k - 8);
        y = ((((((((((((c[k] * t + c[k + 1]) * t + 
            c[k + 2]) * t + c[k + 3]) * t + c[k + 4]) * t + 
            c[k + 5]) * t + c[k + 6]) * t + c[k + 7]) * t + 
            c[k + 8]) * t + c[k + 9]) * t + c[k + 10]) * t + 
            c[k + 11]) * t + c[k + 12]) * t + c[k + 13];
    } else {
        v = 24 / w;
        t = v * v;
        k = 13 * ((int) t);
        y = ((((((d[k] * t + d[k + 1]) * t + 
            d[k + 2]) * t + d[k + 3]) * t + d[k + 4]) * t + 
            d[k + 5]) * t + d[k + 6]) * sqrt(v);
        theta = (((((d[k + 7] * t + d[k + 8]) * t + 
            d[k + 9]) * t + d[k + 10]) * t + d[k + 11]) * t + 
            d[k + 12]) * v - 0.78539816339744830962;
        y *= sin(w + theta);
    }
    return x < 0 ? -y : y;
}

int
main()
{
	double  d;

	//FILE * fp = fopen("Release/bestest.dat","w");
	for (d = 0.0; d < 1<<12; d += 0.001) {
		double m0 = j0(d);
		double m1 = j1(d);
		double d0 = dbesj0(d);
		double d1 = dbesj1(d);
		if (fabs(m0-d0) > 1e-15 || fabs(m1-d1) > 1e-15)
			printf("%9.4f\t%0.16e\t%0.16e\t%0.16e\t%0.16e\n",d,m0,m1,d0,d1);
	}
	//fclose(fp);
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

#ifdef NOBUILD
int
main()
{
	fixed<6> a;
	fixed<6> b(0);
	fixed<6> c(12);
	fixed<6> d(12.34);
	fixed<6> e(b);
	fixed<6> f(-12);

	a = -c;
	b = b + c;
	e += d;
	e -= 1;
	e = (fixed<6>)(0.06) + e;
	++e;
	printf(a == f ? "true\n" : "false\n");
	printf(d > c ? "true\n" : "false\n");
	printf((double)(d) > 11.0 ? "true\n" : "false\n");
	printf("%g\n",(double)(b));
	printf("%g\n",(double)(e));
}
#endif

#ifdef NOBUILD
int
main()
{
	int i;

	for (i = 0; i < 10000000; i++) {
		double u = stdnormal_rnd();
		printf("%0.16e\n",u);
	}
}
#endif

#ifdef NOBUILD
int
main()
{
	double p;

	for (p = 2.5; p <= 6.0; p += 0.001) {
		double u = -pow(2,p);
		double p0 = stdnormal_pdf(u);
		double p1 = stdnormal_cdf(u);
		double p2 = quad8_stdnormal_pdf(-100.0,u,0.0);
		double u0 = stdnormal_inv(p1);
		printf("%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\t%0.16e\n",u,p0,p1,p2,u0,fabs(p1-p2));
	}
}
#endif

#ifdef NOBUILD
#define n 10
#define m 9
int
main()
{
	int i, j, k, r = 0;
	double s, dot = 0.0;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * b = new double[n];
	double * x = new double[n];
	if (A == 0 || B == 0 || b == 0 || x == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}

again:
	printf(".");
	memset(A,0,sizeof(double)*n*n);
	for (i = 0; i < n; i++) {
		for (j = i; j <= i+m && j < n; j++)
			A[i*n+j] = RAND(-1.0,1.0);
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			B[i*n+j] = 0.0;
			for (k = 0; k < n; k++)
				B[i*n+j] += A[i*n+k]*A[j*n+k];
		}
	}
	for (i = 0; i < n; i++)
		b[i] = RAND(0.0,1.0);
	printf("(%d",r++);
	fflush(NULL);
	
	memcpy(A,B,sizeof(double)*n*n);
	equ_gauss(n,A,b,x);
	//equ_lu(n,A,b,x);
	//equ_chol(n,A,b,x);
	//equ_ldl(n,A,b,x);
	//equ_svd(n,A,b,x);
	//equ_eig(n,A,b,x);

	//for (i = 0; i < n; i++) {
	//	for (j = i; j <= i+m && j < n; j++)
	//		A[B_IDX(n,m,i,j)] = B[i*n+j];
	//}
	//equ_chol(n,m,A,b,x,n*n*10e-12);
	
	printf(")");
	double c1, y1, t1, c2, y2, t2;
	for (i = 0, dot = 0.0, c1 = 0.0; i < n; i++) {
		for (j = 0, s = -b[i], c2 = 0.0; j < n; j++) {
			y2 = B[i*n+j]*x[j] - c2; t2 = s + y2;
			c2 = t2 - s - y2; s = t2;
		}
		y1 = s*s - c1; t1 = dot + y1;
		c1 = t1 - dot - y1; dot = t1;
	}
	if (sqrt(dot) > n*n*1e-12) {
		printf("\n%g\nA = [ ",sqrt(dot));
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++)
				printf("%.60e%s ...\n",B[i*n+j],(j == n-1 ? ";" : ","));
		}
		printf("];\n");
		printf(" b = [");
		for (i = 0; i < n; i++) {
			printf("%.60e; ...\n",b[i]);
		}
		printf("];\n x = [");
		for (i = 0; i < n; i++) {
			printf("%.60e; ...\n",x[i]);
		}
		printf("];\n");
		exit(1);
	}
	goto again;

abort:
	delete [] x;
	delete [] b;
	delete [] B;
	delete [] A;
}
#endif

#ifdef NOBUILD
int
main()
{
	int n, m, i, j, k;
	int * A;

	for (n = 0; n < 1000; n++) {
		for (m = 0; m < n; m++) {
			A = new int[B_SIZE(n,m)];
			memset(A,0,sizeof(int)*B_SIZE(n,m));
			for (i = 0; i < n; i++) {
				for (j = i; j >= i-m && j >= 0; j--) {
					k = B_IDX(n,m,j,i);
					if (k < 0 || k >= B_SIZE(n,m)) {
						printf("%d\t%d\t%d\t%d\t%d\t%d\n",n,m,i,j,B_SIZE(n,m),k);
						exit(1);
					}
					printf("%d\t%d\t%d\t%d\t%d\n",n,m,i,j,k);
					if (A[k] == 1) {
						printf("OVERLAP\n");
						exit(1);
					}
					A[k] = 1;
				}
			}
			for (i = 0; i < B_SIZE(n,m); i++) {
				if (A[i] != 1) {
					printf("MISSED\n");
					exit(1);
				}
			}
			delete [] A;
		}
	}
}
#endif

#ifdef NOBUILD
#define n 10
int
main()
{
	bool rv = false;
	int i, j, k, iter = 0;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * I = new double[n*n];
	if (A == 0 || B == 0 || I == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}

again:
	printf(".");
	for (i = 0; i < n*n; i++)
		B[i] = RAND(0.0,1.0), I[i] = A[i] = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				A[i*n+j] += B[i*n+k]*B[j*n+k];
		}
	}
	for (i = 0; i < n*n; i++)
		B[i] = A[i];
	printf("(");
	//inv_lu(n,A);
	//inv_chol(n,A);
	//inv_svd(n,A);
	inv_eig(n,A);
	printf("%d)",++iter);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				I[i*n+j] += A[i*n+k]*B[k*n+j];
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			//printf("%f\t",I[i*n+j]-(i==j?1.0:0.0));
			if (fabs(I[i*n+j]-(i==j?1.0:0.0)) > 1e-6)
				rv = true;
		}
		//printf("\n");
	}
	if (rv)
		exit(1);
	goto again;

abort:
	if (A != 0)
		delete [] A;
	if (B != 0)
		delete [] B;
	if (I != 0)
		delete [] I;

}
#endif

#ifdef NOBUILD
#define n 5
int
main()
{
	int i, j, k;
	double * A = new double[n*n];
	double * B = new double[n*n];
	double * I = new double[n*n];
	double * d = new double[n];
	double * e = new double[n];
	if (A == 0 || d == 0 || e == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		goto abort;
	}

again:
	//printf(".");
	for (i = 0; i < n*n; i++)
		B[i] = RAND(0.0,1.0), I[i] = A[i] = 0.0;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++)
				A[i*n+j] += B[i*n+k]*B[j*n+k];
		}
	}
	for (i = 0; i < n*n; i++)
		B[i] = A[i];
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%8.4f%s",A[i*n+j],(j==n-1?"\n":"\t"));
	}
	printf("\n");
	//printf("(");
	tridiag_hh(n,A,d,e);
	eig_tri_ql(n,d,e,A);
	//printf(")");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%8.4f%s",A[i*n+j],(j==n-1?"\n":"\t"));
	}
	printf("\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == j)
				printf("%8.4f%s",d[i],(j==n-1?"\n":"\t"));
			else if (i == j+1)
				printf("%8.4f%s",e[i],(j==n-1?"\n":"\t"));
			else if (i == j-1)
				printf("%8.4f%s",e[i],(j==n-1?"\n":"\t"));
			else
				printf("%8.4f%s",0.0,(j==n-1?"\n":"\t"));
		}
	}
	//goto again;
abort:
	if (e != 0)
		delete [] e;
	if (d != 0)
		delete [] d;
	if (B != 0)
		delete [] B;
	if (A != 0)
		delete [] A;
	if (I != 0)
		delete [] I;
}
#endif

#ifdef BUILD
#define n 5
int
main()
{
	int i, j;

	double * A = new double[n*n];
	if (A == 0) {
		event_msg(EVENT_ERROR,"Out of memory!");
		exit(1);
	}
	for (i = 0; i < n*n; i++)
		A[i] = RAND(0.0,1.0);

	matrix_dense *s = new matrix_dense(n,n,A);
	matrix a(s);
	matrix b;
	
	b = -a;
	matrix c(!a);

	printf("%g\t%g\t%g\n",a(1,1),b(1,1),c(1,1));

	delete [] A;
}
#endif
