/*************************************************************************

	LIBOP.H - DLL wrapper for OpenPave classes.

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
	double (* res)[27]              // Results
) noexcept;

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
	double * res                    // Results
) noexcept;

long OP_EXPORT
OP_LE_Create() noexcept;

int OP_EXPORT
OP_LE_AddLayer(const long token, const double h, const double E,
	const double v, const double f) noexcept;

int OP_EXPORT
OP_LE_InsertLayer(const long token, const unsigned p, const double h,
	const double E, const double v, const double f) noexcept;

int OP_EXPORT
OP_LE_RemoveLayer(const long token, const unsigned p) noexcept;

int OP_EXPORT
OP_LE_RemoveLayers(const long token) noexcept;

int OP_EXPORT
OP_LE_AddLoad(const long token, const double ax, const double ay,
	const double al, const double ap, const double ar) noexcept;

int OP_EXPORT
OP_LE_RemoveLoad(const long token, const unsigned l) noexcept;

int OP_EXPORT
OP_LE_RemoveLoads(const long token) noexcept;

int OP_EXPORT
OP_LE_AddGroupLoad(const long token, const unsigned g, const double ax,
	const double ay, const double al, const double ap, const double ar) noexcept;

int OP_EXPORT
OP_LE_RemoveGroupLoad(const long token, const unsigned g, const unsigned l) noexcept;

int OP_EXPORT
OP_LE_RemoveGroup(const long token, const unsigned g) noexcept;

int OP_EXPORT
OP_LE_RemoveAllLoadGroups(const long token) noexcept;

int OP_EXPORT
OP_LE_AddPoint(const long token, const double px, const double py,
	const double pz, const unsigned pl) noexcept;

int OP_EXPORT
OP_LE_RemovePoint(const long token, const double px, const double py,
	const double pz, const unsigned pl) noexcept;

int OP_EXPORT
OP_LE_RemovePoints(const long token) noexcept;

int OP_EXPORT
OP_LE_Calculate(const long token, const unsigned flags) noexcept;

int OP_EXPORT
OP_LE_Result(const long token, const unsigned i, double res[27]) noexcept;

int OP_EXPORT
OP_LE_Results(const long token, double (*res)[27]) noexcept;

int OP_EXPORT
OP_LE_GroupResult(const long token, const unsigned g, const unsigned i,
	double res[27]) noexcept;

int OP_EXPORT
OP_LE_GroupResults(const long token, const unsigned g, double (* res)[27]) noexcept;

int OP_EXPORT
OP_LE_Reset(const long token) noexcept;

long OP_EXPORT
OP_HT_Init(const unsigned nl, const double * h, const double * D,
	       const unsigned nn, const double * nd, const double * nt,
	       const unsigned nw, const double dt) noexcept;

int OP_EXPORT
OP_HT_Step(const long token, const unsigned nt,
	       const double * tt, const double tb) noexcept;

int OP_EXPORT
OP_HT_Interpolate(const long token, const unsigned np, const double * pd,
	              double * pt) noexcept;

int OP_EXPORT
OP_HT_Reset(const long token) noexcept;

END_C_DECLS

#endif // LIBOP_H
