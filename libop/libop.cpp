/*************************************************************************

	LIBOP.CPP - DLL wrapper for OpenPave classes.

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

	Portions Copyright (C) 2006-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	See LIBOP.H.

	History:
		2007/06/07 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#if defined(XP_PC)
#include <windows.h>

BOOL WINAPI
DllMain(HANDLE, DWORD fdwReason, LPVOID)
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
#else
#define _EVENT_IMP
#endif

#define _PROGRESS_IMP
#include "event.h"
#include "tree.h"
#include "pavement.h"
#include "thermal.h"
#include "libop.h"

using namespace OP;

#if defined(XP_PC)
/*
 * Windows error event handler.
 */
void OP::event_msg(const int level, const char * fmt, ...)
{
	va_list args;
	char * buf;
	int len;

	va_start(args,fmt);
	if (level < EVENT_DEBUG && fmt != NULL) {
		len = _vscprintf_p(fmt,args)+1;
		buf = (char *)malloc(len*sizeof(char));
		if (buf == NULL) {
			::MessageBoxA(0,"Out of Memory!","LIBOP.DLL Error",MB_OK|MB_ICONERROR);
		} else {
			_vsprintf_p(buf,len,fmt,args);
			::MessageBoxA(0,buf,"LIBOP.DLL Error",MB_OK|MB_ICONERROR);
			free(buf);
		}
	}
	va_end(args);
}
#endif

int OP_EXPORT
OP_LE_Calc(const unsigned flags,
           const unsigned nl, const double * h, const double * E,
             const double * v, const double * f,
           const unsigned na, const double * ax, const double * ay,
             const double * al, const double * ap, const double * ar,
           const unsigned np, const double * px, const double * py,
             const double * pz, const unsigned * pl, double (* res)[27])
{
	LEsystem pave;
	unsigned i;
	bool rv = false;

	for (i = 0; i < nl; i++) {
		pave.addlayer(h[i],E[i],v[i],f[i]);
	}
	for (i = 0; i < na; i++) {
		pave.addload(point2d(ax[i],ay[i]),al[i],ap[i],ar[i]);
	}
	for (i = 0; i < np; i++) {
		if (pl[i] == 0)
			pave.addpoint(point3d(px[i],py[i],pz[i]));
		else
			pave.addpoint(point3d(px[i],py[i],pz[i]),pl[i]-1);
	}
	switch (flags & 0xFF) {
	case 0x00:
	default:
		rv = !pave.calculate(flags > 0xFF
				? LEsystem::disp : LEsystem::all);
		break;
	case 0x01:
		rv = !pave.calculate(flags > 0xFF
				? LEsystem::fastdisp : LEsystem::fast);
		break;
	case 0x02:
		rv = !pave.calculate(flags > 0xFF
				? LEsystem::dirtydisp : LEsystem::dirty);
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
	}
	for (i = 0; i < np; i++) {
		point3d p(px[i],py[i],pz[i]);
		unsigned l = (pl[i] == 0 ? UINT_MAX : pl[i]-1);
		const pavedata & d = pave.result(p,l);
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

int OP_EXPORT
OP_LE_Calc_CalME(const unsigned flags,
           const unsigned nl, const double * h, const double * E,
             const double * v, const double * f,
           const unsigned na, const double * ax, const double * ay,
             const double * ap, const double * ar,
           const unsigned np, const double * px, const double * py,
             const double * pz, const unsigned * pl, double * res)
{
	LEsystem pave;
	unsigned i;
	bool rv = false;

	for (i = 0; i < nl; i++)
		pave.addlayer(h[i],E[i]*1000,v[i],f[i]);
	for (i = 0; i < na; i++)
		pave.addload(point2d(ax[i],ay[i]),0.0,ap[i]*1000,ar[i]);
	for (i = 0; i < np; i++)
		pave.addpoint(point3d(px[i],py[i],pz[i]),pl[i]-1);
	switch (flags & 0xFF) {
	case 0x00:
	default:
		rv = !pave.calculate(flags > 0xFF
				? LEsystem::disp : LEsystem::all);
		break;
	case 0x01:
		rv = !pave.calculate(flags > 0xFF
				? LEsystem::fastdisp : LEsystem::fast);
		break;
	case 0x02:
		rv = !pave.calculate(flags > 0xFF
				? LEsystem::dirtydisp : LEsystem::dirty);
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
	}
	for (i = 0; i < np; i++) {
		const pavedata & d = pave.result(point3d(px[i],py[i],pz[i]),pl[i]-1);
		res[ 1*150+i] = -d.result(pavedata::stress,pavedata::yy)/1000;
		res[ 2*150+i] = -d.result(pavedata::stress,pavedata::zz)/1000;
		res[ 4*150+i] = -d.result(pavedata::stress,pavedata::xz)/1000;
		res[15*150+i] = -d.result(pavedata::strain,pavedata::xx);
		res[16*150+i] = -d.result(pavedata::strain,pavedata::yy);
		res[17*150+i] = -d.result(pavedata::strain,pavedata::zz);
		res[ 9*150+i] = -d.result(pavedata::strain,pavedata::p1);
		res[10*150+i] = -d.result(pavedata::strain,pavedata::p2);
		res[11*150+i] = -d.result(pavedata::strain,pavedata::p3);
	}
#if defined(_MSC_VER) || defined(__MINGW32__)
	_clearfp();
#endif
	return rv;
}

ktree_avl<long, FEMthermal *> tokens;

long OP_EXPORT
OP_HT_Init(const unsigned nl, const double * h, const double * D,
           const unsigned nn, const double * nd, const double * nt,
           const unsigned nw, const double dt)
{
	unsigned i;
	long token = 0;

	if (!(nw == 1 || nw == 2 || nw == 4)) {
		event_msg(EVENT_WARN, "Invalid bandwidth!");
		return 0;
	}
	if (nn < nw+1 || (nn-1)%nw != 0) {
		event_msg(EVENT_WARN, "Bad temperature mesh! Nodes must match bandwidth!");
		return 0;
	}
	// Set up layer top and bottom.
	for (i = 0; i < nl-1; i++) {
		if (h[i] <= 0.0) {
			event_msg(EVENT_WARN,"Invalid Layer thickness!");
			return 0;
		}
	}
	FEMthermal * rv = new FEMthermal(nl,h,D,nn,nd,nt,dt,nw);
	while (tokens.haskey(token))
		token++;
	tokens.add(token, rv);

#if defined(_MSC_VER) || defined(__MINGW32__)
	_clearfp();
#endif
	return token;
}

void OP_EXPORT
OP_HT_Step(const long token, const unsigned nt,
           const double * tt, const double tb)
{
	if (!tokens.haskey(token)) {
		event_msg(EVENT_ERROR, "Invalid token!");
		return;
	}
	FEMthermal * system = tokens.get(token);

	for (unsigned i = 0; i < nt; i++)
		system->step(tt[i],tb);
#if defined(_MSC_VER) || defined(__MINGW32__)
	_clearfp();
#endif
}

void OP_EXPORT
OP_HT_Interpolate(const long token, const unsigned np, const double * pd,
                  double * pt)
{
	if (!tokens.haskey(token)) {
		event_msg(EVENT_ERROR, "Invalid token!");
		return;
	}
	FEMthermal * system = tokens.get(token);

	system->interpolate(np,pd,pt);
#if defined(_MSC_VER) || defined(__MINGW32__)
	_clearfp();
#endif
}

void OP_EXPORT
OP_HT_Reset(const long token)
{
	if (!tokens.haskey(token)) {
		event_msg(EVENT_ERROR, "Invalid token!");
		return;
	}
	FEMthermal * system = tokens.get(token);
	tokens.remove(token);
	delete system;
}
