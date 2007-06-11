/*************************************************************************

	LIBOP.CPP - DLL wrapper for OpenPave classes.

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

	See LIBOP.H.

	History:
		2007/06/07 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#if defined(XP_PC)
#include <windows.h>

BOOL WINAPI
DllMain(HANDLE hModule, DWORD fdwReason, LPVOID lpReserved)
{
	switch (fdwReason) {
	case DLL_PROCESS_ATTACH:
		break;
	case DLL_THREAD_ATTACH:
		break;
	case DLL_THREAD_DETACH:
		break;
	case DLL_PROCESS_DETACH:
		break;
	}
	return TRUE;
}
#endif

#define _EVENT_IMP
#define _PROGRESS_IMP
#include "event.h"
#include "pavement.h"
#include "libop.h"

int OP_EXPORT
OP_LE_Calc(const int flags,
           const int nl, const double * h, const double * E,
             const double * v, const double * f,
           const int na, const double * ax, const double * ay,
             const double * al, const double * ap, const double * ar,
		   const int np, const double * px, const double * py,
		     const double * pz, double (* res)[27])
{
	LEsystem pave;
	int i;
	bool rv = false;

	for (i = 0; i < nl; i++) {
		pave.addlayer(h[i],E[i],v[i],f[i]);
	}
	for (i = 0; i < na; i++) {
		pave.addload(point2d(ax[i],ay[i]),al[i],ap[i],ar[i]);
	}
	for (i = 0; i < np; i++) {
		pave.addpoint(point3d(px[i],py[i],pz[i]));
	}
	switch (flags & 0xFF) {
	case 0x00:
		rv = !pave.calculate(flags > 0xFF ? LEsystem::disp : LEsystem::all);
		break;
	case 0x01:
		rv = !pave.calculate(flags > 0xFF ? LEsystem::fastdisp : LEsystem::fast);
		break;
	case 0x02:
		rv = !pave.calculate(LEsystem::dirty);
		break;
	case 0x03:
		rv = !pave.calc_odemark();
		break;
	case 0x04:
		rv = !pave.calc_fastnum();
		break;
	case 0xFF:
		rv = !pave.calc_accurate();
		break;
	default:
		rv = !pave.calculate(flags > 0xFF ? LEsystem::disp : LEsystem::all);
		break;
	}
	for (i = 0; i < np; i++) {
		point3d p(px[i],py[i],pz[i]);
		const pavedata & d = pave.result(p);
		res[i][ 0] = d.result(pavedata::stress,pavedata::xx);
		res[i][ 1] = d.result(pavedata::stress,pavedata::yy);
		res[i][ 2] = d.result(pavedata::stress,pavedata::zz);
		res[i][ 3] = d.result(pavedata::stress,pavedata::xy);
		res[i][ 4] = d.result(pavedata::stress,pavedata::xz);
		res[i][ 5] = d.result(pavedata::stress,pavedata::yz);
		res[i][ 6] = d.result(pavedata::stress,pavedata::p1);
		res[i][ 7] = d.result(pavedata::stress,pavedata::p2);
		res[i][ 8] = d.result(pavedata::stress,pavedata::p3);
		res[i][ 9] = d.result(pavedata::stress,pavedata::s1);
		res[i][10] = d.result(pavedata::stress,pavedata::s2);
		res[i][11] = d.result(pavedata::stress,pavedata::s3);
		res[i][12] = d.result(pavedata::deflct,pavedata::xx);
		res[i][13] = d.result(pavedata::deflct,pavedata::yy);
		res[i][14] = d.result(pavedata::deflct,pavedata::zz);
		res[i][15] = d.result(pavedata::strain,pavedata::xx);
		res[i][16] = d.result(pavedata::strain,pavedata::yy);
		res[i][17] = d.result(pavedata::strain,pavedata::zz);
		res[i][18] = d.result(pavedata::strain,pavedata::xy);
		res[i][19] = d.result(pavedata::strain,pavedata::xz);
		res[i][20] = d.result(pavedata::strain,pavedata::yz);
		res[i][21] = d.result(pavedata::strain,pavedata::p1);
		res[i][22] = d.result(pavedata::strain,pavedata::p2);
		res[i][23] = d.result(pavedata::strain,pavedata::p3);
		res[i][24] = d.result(pavedata::strain,pavedata::s1);
		res[i][25] = d.result(pavedata::strain,pavedata::s2);
		res[i][26] = d.result(pavedata::strain,pavedata::s3);
	}
#if defined(_MSC_VER) || defined(__MINGW32__)
	_clearfp();
#endif
	return rv;
}
