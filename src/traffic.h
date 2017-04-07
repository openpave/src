/**************************************************************************

	TRAFFIC.H - Interface for traffic related classes.

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
		2002/10/05 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __TRAFFIC_H
#define __TRAFFIC_H

#include <stdio.h>
#include <string.h>
#include <time.h>

namespace OP {

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
	WIMday(time_t d = -1) {
		Clear(d);
	}
	// Clear the strucuture
	void Clear(time_t d = -1) {
		memset(this,0,sizeof(*this));
		day = d;
	}
	// Write binary file...
	bool Write(FILE * bp) {
		return (fwrite(this,sizeof(*this),1,bp) == 1);
	}
	// Read binary file...
	bool Read(FILE * bp) {
		return (fread(this,sizeof(*this),1,bp) == 1);
	}
	static time_t GetRSADate(const char * fname);
	bool AddRSAVehicle(const char * buf);
	WIMday & operator += (const WIMday & daily);

	// The total vehicle count.
	int TV() const {
		return TLV()+THV();
	}
	// The total light vehicle count.
	int TLV() const {
		int l, rv;
		for (l = 0, rv = 0; l < MAX_LANES; l++)
			rv += LVPL[l];
		return rv;
	}
	// The total heavy vehicle count.
	int THV() const {
		int l, rv;
		for (l = 0, rv = 0; l < MAX_LANES; l++)
			rv += HVPL[l];
		return rv;
	}
	// The total vehicle count in lane L.
	int VPL(int l) const {
		return LVPL[l]+HVPL[l];
	}
	// The total axle count in lane L.
	int TAC(int l) const {
		int t, rv;
		for (t = 0, rv = 0; t < AXLE_TYPES; t++)
			rv += AC(l,t);
		return rv;
	}
	// The total axle count for axle type T.
	int AxC(int t) const {
		int l, rv;
		for (l = 0, rv = 0; l < MAX_LANES; l++)
			rv += AC(l,t);
		return rv;
	}
	// The axle count in lane L for axle type T.
	int AC(int l, int t) const {
		int i, rv;
		for (i = 0, rv = 0; i < LOADS; i++)
			rv += WG[l][t][i];
		return rv;
	}
};

/*
 * Cummulative WIM survey results
 */
class WIMsurvey : public WIMday {
public:
	// Your basic constructor...
	WIMsurvey(time_t s = -1, time_t e = -1, const char * fname = 0) {
		memset(this,0,sizeof(*this));
		if (s != -1)
			start = s;
		if (e != -1)
			end = e;
		if (fname != 0)
			Read(fname);
	}
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

} // namespace OP

#endif // TRAFFIC_H
