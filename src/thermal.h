/**************************************************************************

	THERMAL.H - Interface for the pavement class.

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

	Portions Copyright (C) 2008 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This class implements a 1-D FEM heat transfer solution combined with
		a forward finite difference time step to solve the temperatures
		within a pavement.

		The code is run by calling the constructor with all of the initial
		parameters and then calling step to calculate at each new timestep.
		Step makes use of the first and last temperatures in the passed in
		array (the surface and at-depth temperatures).

		Interpolate can be used to get temperatures at non-nodal depths.

	Future:
		This class needs to be designed more carefully to accommodate
		thermal and moisture properties.

	History:
		2008/11/21 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __THERMAL_H
#define __THERMAL_H

namespace OP {

class FEMthermal
{
public:
	FEMthermal(unsigned nl, const double * lh, const double * ld,
			unsigned _n, const double * _nd, const double * _nt,
			double _dt = 1, unsigned _w = 1);
	~FEMthermal() {
		delete [] nd;
		delete [] nt;
		delete [] ng;
		delete [] KK;
		delete [] FF;
		delete [] Kt;
		delete [] Kb;
	}
	// Perform the FEM analysis for a single timestep, based in the top and
	// bottom temperature at the desired time.  This function can be called
	// repeatedly to obtain temperatures at sucessive time steps.
	void step(double tt, double tb);
	// Interpolate temperatures at abitrary depths. pd must be sorted.
	void interpolate(unsigned np, const double * pd, double * pt);

private:
	const unsigned n;  			// Number of nodes
	const unsigned w;			// Number of bands (element size)
	const double dt;			// Time delta
	double * nd;				// Nodal depths
	double * nt;				// Nodal temperatures
	double * ng;				// Nodal temp. gradient
	double * KK;				// Nodal diffusivity matrix
	double * FF;				// Finite diff matrix
	double * Kt;				// Saved first col.
	double * Kb;				// Saved last col.
};

} // namespace OP

#endif // THERMAL_H
