/**************************************************************************

	TRAFFIC.H - Interface for traffic related classes.

	$OpenPave$

	Purpose:
		This file implements various traffic related classes.  Most
		importantly it can convert SANRAL style .RSA WIM data files into
		our .WIM binary file format, and provides facilities to save and
		load these files.

	Future:
		In the future this file will contain classes to handle LEF
		calculations, and determination of traffic growth and cumulative
		ESALs with time.

	History:
		Based on AXW-FIN.PAS from Coen Coetzee.
		2002/10/05 - Created by Jeremy Lea <jlea@csir.co.za>

**************************************************************************/

#ifndef __TRAFFIC_H
#define __TRAFFIC_H

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <time.h>

// Known axle types.
#define STEER		0
#define DUAL		1
#define TANDEM		2
#define TRIDEM		3
#define BEDA		4	// Bogey equivalent dual axles (for all > 3).
#define AXLE_TYPES	5

#define MAX_AXLES	15	// Max number of axles per vehicle.
#define MAX_LANES	6	// Maximum number of lanes in survey.
#define LOADS		21
const double load[LOADS-1] = {
	 10,  20,  30,  40,  50,  60,  70,  80,  90, 100,
	110, 120, 130, 140, 150, 160, 170, 180, 190, 200
};

/*
 * Daily count structure.
 */
struct WIMday {
	time_t day;								// Day of the count
	int SV;									// Short vehicles (ignored)
	int LVPL[MAX_LANES];					// Light Vehicles per lane
	int HVPL[MAX_LANES];					// Heavy Vehicles per lane
	int WG[MAX_LANES][AXLE_TYPES][LOADS];	// The actual WIM counts
	double AAS[AXLE_TYPES];					// Average axle spacing

	// Your basic constructor...
	inline WIMday(time_t d = -1) {
		Clear(d);
	};
	// Clear the strucuture
	inline void Clear(time_t d = -1) {
		memset(this,0,sizeof(*this));
		day = d;
	};
	// Write binary file...
	inline bool Write(FILE * bp) {
		return (fwrite(this,sizeof(*this),1,bp) == 1);
	};
	// Read binary file...
	inline bool Read(FILE * bp) {
		return (fread(this,sizeof(*this),1,bp) == 1);
	};
	static time_t GetRSADate(const char * fname);
	bool AddRSAVehicle(const char * buf);
	WIMday & operator += (const WIMday & daily);
	
	// The total vehicle count.
	inline int TV() const {
		return TLV()+THV();
	};
	// The total light vehicle count.
	inline int TLV() const {
		int l, rv;
		for (l = 0, rv = 0; l < MAX_LANES; l++)
			rv += LVPL[l];
		return rv;
	};
	// The total heavy vehicle count.
	inline int THV() const {
		int l, rv;
		for (l = 0, rv = 0; l < MAX_LANES; l++)
			rv += HVPL[l];
		return rv;
	};
	// The total vehicle count in lane L.
	inline int VPL(const int l) const {
		return LVPL[l]+HVPL[l];
	};
	// The total axle count in lane L.
	inline int TAC(const int l) const {
		int t, rv;
		for (t = 0, rv = 0; t < AXLE_TYPES; t++)
			rv += AC(l,t);
		return rv;
	};
	// The total axle count for axle type T.
	inline int AxC(const int t) const {
		int l, rv;
		for (l = 0, rv = 0; l < MAX_LANES; l++)
			rv += AC(l,t);
		return rv;
	};
	// The axle count in lane L for axle type T.
	inline int AC(const int l, const int t) const {
		int i, rv;
		for (i = 0, rv = 0; i < LOADS; i++)
			rv += WG[l][t][i];
		return rv;
	};
};

/*
 * Cummulative WIM survey results
 */
class WIMsurvey : public WIMday {
public:
	// Your basic constructor...
	inline WIMsurvey(time_t s = -1, time_t e = -1,const char * fname = 0) {
		memset(this,0,sizeof(*this));
		if (s != -1)
			start = s;
		if (e != -1)
			end = e;
		if (fname != 0)
			Read(fname);
	};
	bool ProcessRSAFile(const char * fname, FILE * bp, WIMday * day);
	bool ProcessRSADir(const char * dir, const char * bname);
	bool Read(const char * bname);

private:
	time_t start;							// First day of count
	time_t end;								// Last day of count
};

int str2int(const char * c, int w);
time_t str2time(const char * s, struct tm * date);
int fexists(const char * fname);
void fgetln(char * b, int l, FILE * fp);

#endif // TRAFFIC_H

