/*************************************************************************

	TRAFFIC.CPP - Implementations for traffic related classes.

	$OpenPave$

	See TRAFFIC.H.

	History:
		2002/11/05 - Created by Jeremy Lea <jlea@csir.co.za>

*************************************************************************/

#include "traffic.h"
#include "../include/set.h"
#include <ctype.h>
#if defined(_MSC_VER)
#include <direct.h>
#include <io.h>
#define strdup		_strdup
#define chdir 		_chdir
#define getcwd		_getcwd
#define stat		_stat
#define S_IFREG		_S_IFREG
#else
#include <unistd.h>
#include <glob.h>
#endif
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

/*
 * Make an integer from a fixed width string.
 */
int str2int(const char * c, int w)
{
	int n = 0, r = 0;

	if (c == NULL || w <= 0)
		return 0;
	while (w > 0 && isspace(*c))
		w--, c++;
	while (w-- > 0) {
		if (r == 0 && *c == '-')
			n = 1;
		else if (!isdigit(*c))
			break;
		r = r*10 + *c++ - '0';
	}
	return (n ? -r : r);
}

/*
 * Convert a date in YYMMDD format to a tm strucuture, and
 * return the Unix style time.
 */
time_t str2time(const char * s, struct tm * date)
{
	static struct tm fake;
	int d;
	
	d = str2int(s,6);
	if (s == NULL || d == 0)
		return -1;
	if (d < 690000)
		d += 1000000;
	if (date == NULL)
		date = &fake;
	date->tm_isdst = 0;
	date->tm_hour = date->tm_min = date->tm_sec = 0;
	date->tm_year = d/10000, d -= (d/10000)*10000;
	date->tm_mon = d/100 - 1, d -= (d/100)*100;
	date->tm_mday = d;
	return mktime(date);
}

/*
 * Quick check to see if a file exists.
 */
int fexists(const char * fname)
{
	struct stat sb;

	return (stat(fname, &sb) == 0 && sb.st_mode & S_IFREG);
}

/*
 * Really read a line from a file, cleaning as we go...
 * b is the buffer, which should be l long.  Only l-1 chars are placed
 * into b, and the last character is always '\0'.  However, the procedure
 * keeps reading until it finds the end of the line.  Only normal chars
 * are placed into the string.  Anything funny is replaced by a space.
 */
void fgetln(char * b, int l, FILE * fp)
{
	int c;

	if (b == NULL || l <= 0)
		return;
	while (!feof(fp)) {
		c = fgetc(fp);
		if (c == '\n') {
			if (*b == '\r')
				*b = '\0';
			break;
		}
		if (l > 1)
			*b++ = (isprint(c) ? c : ' '), l--;
	}
	while (l > 0)
		*b++ = '\0', l--;
}

/*
 * Read a .RSA file and retrieve the first date from the file.
 */
time_t WIMday::GetRSADate(const char * fname)
{
	FILE * fp;
	char buf[256];
	struct tm d;

	if (fname == NULL)
		return (time_t)-1;
	if ((fp = fopen(fname, "r")) == NULL)
		return (time_t)-1;
	while (!feof(fp)) {
		fgetln(buf,256,fp);
		if (!strncmp(buf,"H9",2))
			break;
	}
	while (!feof(fp)) {
		fgetln(buf,256,fp);
		if (!strncmp(buf,"13",2))
			break;
	}
	if (feof(fp)) {
		fclose(fp);
		return (time_t)-1;
	}
	fclose(fp);
	return str2time(&buf[5],&d);
}

/*
 * Process a '13' line from a .RSA file and add the vehicle's axles
 * to the various bins.
 */
bool WIMday::AddRSAVehicle(const char * buf)
{
	int i, l, ia, ib;
	int il = 0, na = 0, at = STEER;
	double L = 0.0, GVM = 0.0;
	double AS[MAX_AXLES],AW[MAX_AXLES];

	if (str2time(&buf[5],NULL) != day)
		return false;
	il = str2int(&buf[24],2)-1;			// Lane of measurement
	if (il < 0 || il >= MAX_LANES)		// which must be in range
		return false;
	L = 0.01*str2int(&buf[33],4);		// Vehicle length in metres
	if (L < 2.5) {						// which needs to be realistic
		SV++;
		return true;
	}
	na = str2int(&buf[51],2);			// No. of axles
	if (na < 0 || na > MAX_AXLES)		// which must be in range
		return false;
	for (ia = 0; ia < na; ia++) {
		AW[ia] = 0.1*str2int(&buf[54+ia*9],3);	// Axle weight in Tonne
		GVM += AW[ia];
	}
	for (ia = 1, AS[0] = 0.0; ia < na; ia++)
		AS[ia] = 0.01*str2int(&buf[49+ia*9],4);	// Axle spacing in metres
	// XXX: Correction????
	if (na == 2 && AS[1] <= 2.9 && GVM >= 3.0)
		AW[0] *= 2.5/GVM, AW[1] *= 2.5/GVM, GVM = 2.5;
	// Sort Light and Heavy vehicles
	if ((na == 0) 
	 || (na == 1 && GVM < 2.0 && L <  8.5)
	 || (na == 2 && GVM < 3.5 && L < 10.8)
	 || (na == 3 && GVM < 5.5)
	 || (na == 4 && GVM < 7.0)) {
		LVPL[il]++;
		return true;
	} else {
		HVPL[il]++;
	}
	// Determine axle combinations
	ia = 0;								// First axle in bogey		
	while (ia < na) {
		ib = 1;							// No. of axles in bogey
		while (ia+ib < na && AS[ia+ib] <= 1.7)
			ib++;
		if (ib > 3)
			at = BEDA;
		else if (ib == 3)
			at = TRIDEM;
		else if (ib == 2)
			at = TANDEM;
		else if (ia == 0 && na > 1)
			at = STEER;
		else
			at = DUAL;
		if (at == BEDA) {
			// For more than three axles, add each individually.
			for (i = 0; i < ib; i++) {
				// Calculate correct load group.
				for (l = 0; l < LOADS-1; l++)
					if (AW[ia+i] < load[l]/10)
						break;
				WG[il][at][l]++;
			}
		} else {
			// Calculate correct load group.
			double W = 0.0;
			for (i = 0; i < ib; i++)
				W += AW[ia+i];
			for (l = 0, W /= ib; l < LOADS-1; l++)
				if (W < load[l]/10.0)
					break;
			WG[il][at][l]++;
			// Add in the axle spacings, which will be averaged later.
			for (i = 1; i < ib; i++)
				AAS[at] += AS[ia+i]/(ib-1);
		}
		ia += ib;
	}
	return true;
}

/*
 * Add the totals for two days.
 */
WIMday & WIMday::operator += (const WIMday & daily)
{
	int il, ia, i;

	day = daily.day;
	SV += daily.SV;
	for (ia = 0; ia < AXLE_TYPES; ia++) {
		int ac1 = AxC(ia), ac2 = daily.AxC(ia);
		AAS[ia] = (ac1*AAS[ia]+ac2*daily.AAS[ia])/(ac1+ac2);
	}
	for (il = 0; il < MAX_LANES; il++) {
		LVPL[il] += daily.LVPL[il];
		HVPL[il] += daily.HVPL[il];
		for (ia = 0; ia < AXLE_TYPES; ia++) {
			for (i = 0; i < LOADS; i++)
				WG[il][ia][i] += daily.WG[il][ia][i];
		}
	}
	return *this;
}

/*
 * Process a .RSA file, reading the vehicles into day
 */
bool WIMsurvey::ProcessRSAFile(const char * fname, FILE * bp, WIMday * day)
{
	char buf[256];
	time_t current;
	FILE * fp;
	
	if ((fp = fopen(fname, "r")) == NULL)
		return false;
	while (!feof(fp)) {
		fgetln(buf,256,fp);
		// H9 always starts the data
		if (!strncmp(buf,"H9",2))
			break;
	}
	if (feof(fp))
		return false;
	while (!feof(fp)) {
		fgetln(buf,256,fp);
		// H0 always ends the data
		if (!strncmp(buf,"H0",2))
			break;
		// 13 is a vehicle
		if (!strncmp(buf,"13",2)) {
			current = str2time(&buf[5],NULL);
			// New day.  Save away the old data and clear.
			if (current != day->day) {
				day->Write(bp);
				if (start <= day->day && end >= day->day)
					*this += *day;
				day->Clear(current);
			}
			day->AddRSAVehicle(buf);
		}
	}
    fclose(fp);
	return true;
}

/*
 * Read in a binary datafile and sum the data between start and end dates.
 */
bool WIMsurvey::Read(const char * bname)
{
	FILE * bp;
	WIMday day;
	
	if ((bp = fopen(bname, "rb")) == NULL)
		return false;
	while (!feof(bp)) {
		day.Read(bp);
		if (start <= day.day && end >= day.day)
			*this += day;
	}
    fclose(bp);
	return true;
}

bool WIMsurvey::ProcessRSADir(const char * dir, const char * bname)
{
	int i;
	char pwd[FILENAME_MAX], * fname;
#if defined(_MSC_VER)
    _finddata_t found_file;
    long fn;
#else
	glob_t g;
	const char * fn;
#endif
	FILE * bp;
	aoset<char *,time_t> RSAfiles;
	WIMday day;

	if (getcwd(pwd,FILENAME_MAX) == NULL) {
		event_msg(EVENT_ERROR,"Unable to get current working directory!");
		return false;
	}
	if (chdir(dir)) {
		event_msg(EVENT_ERROR,"Unable to change to directory '%s'!",dir);
		return false;
	}
#if defined(_MSC_VER)
	if ((fn = _findfirst("*.RSA",&found_file)) == -1L) {
		event_msg(EVENT_ERROR,"No .RSA files found in direcory '%s'!\n",dir);
		return false;
	} else {
		do {
			if ((fname = strdup(found_file.name)) == NULL
			 || !RSAfiles.add(fname,WIMday::GetRSADate(fname))) {
				event_msg(EVENT_FATAL,"Out of memeory in WIMsurvey::ProcessRSADir()!");
				return false;
			}
		} while (_findnext(fn,&found_file) == 0);
		_findclose(fn);
	}
#else
	memset(&g,0,sizeof(glob_t));
	if (glob("*.RSA",0,NULL,&g) || g.gl_pathc == 0) {
		event_msg(EVENT_ERROR,"No .RSA files found in direcory '%s'!\n",dir);
		return false;
	} else {
		fn = *(g.gl_pathv);
		do {
			if ((fname = strdup(fn)) == NULL
			 || !RSAfiles.add(fname,WIMday::GetRSADate(fname))) {
				event_msg(EVENT_FATAL,"Out of memeory in WIMsurvey::ProcessRSADir()!");
				return false;
			}
		} while (++fn != 0);
	}
	globfree(&g);
#endif
	RSAfiles.sort();
	bp = fopen(bname,"wb");
	day.day = RSAfiles.getvalue(1);
	for (i = 0; i < RSAfiles.length(); i++) {
		fname = RSAfiles.getkey(i+1);
		WIMsurvey::ProcessRSAFile(fname,bp,&day);
		free(fname);
	}
	day.Write(bp);
	if (start <= day.day && end >= day.day)
		*this += day;
	fclose(bp);
	chdir(pwd);
	return true;
}

