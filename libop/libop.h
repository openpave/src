/*************************************************************************

	LIBOP.H - DLL wrapper for OpenPave classes.

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

	Purpose:
		LIBOP consists of various thin C wrappers around the OpenPave
		classes, which allow them to be used in from C, C++ or Visual
		Basic, or any other language for that matter.

	History:
		2007/06/07 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#ifndef __LIBOP_H
#define __LIBOP_H


BEGIN_C_DECLS

/*
 * OP_LE_Calc - Layered elastic calculation with circular loads.
 *
 * Res will hold 27 values on return
 *   1- 6: sxx, syy, szz, sxy, sxz, syz
 *   7-12: sp1, sp2, sp3, ss1, ss2, ss3
 *  13-15: dxx, dyy, dzz
 *  16-21: exx, eyy, ezz, exy, exz, eyz
 *  22-27: ep1, ep2, ep3, es1, es2, es3
 * 
 * For Visual Basic use:
Declare Function OP_LE_Calc Lib "libop.dll" Alias "_OP_LE_Calc@68" ( _
    ByVal flags as Long,
    ByVal nl As Long, h As Double, E As Double, v As Double, f As Double_
    ByVal na As Long, Ax As Double, Ay As Double, Al As Double, Ap As Double, Ar as Double_
    ByVal np As Long, Px As Double, Py As Double, Pz As Double, Pl As Long, Res As Double) _
  As Long
 *
 * Note: VB uses Fortan array storage (cols first) so the results are Res(1,i).
 */
int OP_EXPORT
OP_LE_Calc(
	const unsigned flags,           // Flags to choose method
	const unsigned nl,              // Number of layers
	const double * h,               // Layer thickness (0 for semi-inf)
	const double * E,               // Elastic modulus
	const double * v,               // Poisson's ratio
	const double * f,               // Friction (0.0 to 1.0)
	const unsigned na,              // Number of loads
	const double * ax,              // Center X location
	const double * ay,              // Center Y location
	const double * al,              // Load (0 for auto) 
	const double * ap,              // Pressure (0 for auto)
	const double * ar,              // Radius (0 for auto)
	const unsigned np,              // Number of evaluation points
	const double * px,              // Point X
	const double * py,              // Point Y
	const double * pz,              // Point Z
	const unsigned * pl,            // Point layer (0 for auto)
	double (* res)[27]);            // Results

int OP_EXPORT
OP_LE_Calc_CalME(
	const unsigned flags,           // Flags to choose method
	const unsigned nl,              // Number of layers
	const double * h,               // Layer thickness (0 for semi-inf)
	const double * E,               // Elastic modulus
	const double * v,               // Poisson's ratio
	const double * f,               // Friction (0.0 to 1.0)
	const unsigned na,              // Number of loads
	const double * ax,              // Center X location
	const double * ay,              // Center Y location
	const double * ap,              // Pressure
	const double * ar,              // Radius
	const unsigned np,              // Number of evaluation points
	const double * px,              // Point X
	const double * py,              // Point Y
	const double * pz,              // Point Z
	const unsigned * pl,            // Point layer (0 for auto)
	double * res);                  // Results

long OP_EXPORT
OP_HT_Init(const unsigned nl, const double * h, const double * D,
           const unsigned nn, const double * nd, const double * nt,
           const unsigned nw, const double dt);

void OP_EXPORT
OP_HT_Step(const long token, const unsigned nt,
           const double * tt, const double tb);

void OP_EXPORT
OP_HT_Interpolate(const long token, const unsigned np, const double * pd,
                  double * pt);

void OP_EXPORT
OP_HT_Reset(const long token);

END_C_DECLS

#endif // LIBOP_H
