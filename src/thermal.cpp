/*************************************************************************

	THERMAL.CPP - Implementation for the thermal class.

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

	Portions Copyright (C) 2008-2022 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	See THERMAL.H.

	History:
		2008/11/21 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#include <memory>
#include <cstring>
#include <stdexcept>
#include "autodelete.h"
#include "linalg.h"
#include "thermal.h"

namespace OP {

/*
 * Thermal element matrices.
 */
const double Ke1[2][2] = {
	{ 1.0/2.0,-1.0/2.0 },
	{-1.0/2.0, 1.0/2.0 }
};
const double Fe1[2][2] = {
	{ 2.0/3.0, 1.0/3.0 },
	{ 1.0/3.0, 2.0/3.0 }
};

const double Ke2[3][3] = {
	{ 7.0/6.0,-4.0/3.0, 1.0/6.0 },
	{-4.0/3.0, 8.0/3.0,-4.0/3.0 },
	{ 1.0/6.0,-4.0/3.0, 7.0/6.0 }
};
const double Fe2[3][3] = {
	{ 4.0/15.0, 2.0/15.0,-1.0/15.0 },
	{ 2.0/15.0,16.0/15.0, 2.0/15.0 },
	{-1.0/15.0, 2.0/15.0, 4.0/15.0 }
};

const double Ke4[5][5] = {
	{  985.0/ 378.0,-3424.0/945.0,  508.0/315.0, -736.0/945.0,  347.0/1890.0 },
	{-3424.0/ 945.0, 1664.0/189.0,-2368.0/315.0, 2944.0/945.0, -736.0/ 945.0 },
	{  508.0/ 315.0,-2368.0/315.0,  248.0/ 21.0,-2368.0/315.0,  508.0/ 315.0 },
	{ -736.0/ 945.0, 2944.0/945.0,-2368.0/315.0, 1664.0/189.0,-3424.0/ 945.0 },
	{  347.0/1890.0, -736.0/945.0,  508.0/315.0,-3424.0/945.0,  985.0/ 378.0 }
};
const double Fe4[5][5] = {
	{ 292.0/2835.0, 296.0/2835.0, -58.0/945.0,   8.0/ 405.0, -29.0/2835.0 },
	{ 296.0/2835.0, 256.0/ 405.0,-128.0/945.0, 256.0/2835.0,   8.0/ 405.0 },
	{ -58.0/ 945.0,-128.0/ 945.0, 208.0/315.0,-128.0/ 945.0, -58.0/ 945.0 },
	{   8.0/ 405.0, 256.0/2835.0,-128.0/945.0, 256.0/ 405.0, 296.0/2835.0 },
	{ -29.0/2835.0,   8.0/ 405.0, -58.0/945.0, 296.0/2835.0, 292.0/2835.0 }
};

FEMthermal::FEMthermal(unsigned nl, const double * lh, const double * ld,
		unsigned _n, const double * _nd, const double * _nt,
		double _dt, unsigned _w)
  : n(_n), w(_w), dt(_dt), nd(nullptr), nt(nullptr), ng(nullptr),
		KK(nullptr), FF(nullptr), Kt(nullptr), Kb(nullptr)
{
	unsigned il, i, j, k;

	if (n < w+1)
		throw std::logic_error("Cannot have less elements than bands!");
	if ((n-1)%w > 0)
		throw std::logic_error("Number of elements must be a multiple of the number of bands!");
	autodelete<double> lb(new double[nl]);	// Layer bottom
	nd = new double[n];
	nt = new double[n];
	ng = new double[n];
	KK = new double[B_SIZE(n,w)];
	FF = new double[B_SIZE(n,w)];
	Kt = new double[n];
	Kb = new double[n];
	memcpy(nd,_nd,n*sizeof(double));
	memcpy(nt,_nt,n*sizeof(double));
	for (i = 0; i < nl; i++) {
		if (lh[i] <= 0.0)
			throw std::range_error("Layer thickness must be positive!");
		lb[i] = (i == 0 ? nd[0] : lb[i-1]) + lh[i];
	}

	// Set up the element 'stiffness' matrices.  Since we are working in
	// 1-D we can integrate the shape functions directly, so there is no
	// gaussian integration involved.
	memset(KK,0,B_SIZE(n,w)*sizeof(double));
	memset(FF,0,B_SIZE(n,w)*sizeof(double));
	memset(Kt,0,n*sizeof(double));
	memset(Kb,0,n*sizeof(double));
	for (j = 0, il = 0; j < n - w; j += w) {
		while (lb[il] <= nd[j])
			il++;
		if (lb[il] - lh[il] > nd[j] || lb[il] < nd[j+w])
			throw std::logic_error("Node coordinates and layers do not align!");
		const double d = ld[il]*2/(nd[j+w]-nd[j]);
		const double f = (nd[j+w]-nd[j])/2;
		switch (w) {
		case 1:
			for (i = 0; i <= 1; i++) {
				for (k = i; k <= 1; k++) {
					KK[B_IDX(n,w,j+i,j+k)] += d*Ke1[i][k];
					FF[B_IDX(n,w,j+i,j+k)] += f*Fe1[i][k];
				}
			}
			break;
		case 2:
			for (i = 0; i <= 2; i++) {
				for (k = i; k <= 2; k++) {
					KK[B_IDX(n,w,j+i,j+k)] += d*Ke2[i][k];
					FF[B_IDX(n,w,j+i,j+k)] += f*Fe2[i][k];
				}
			}
			break;
		case 4:
			for (i = 0; i <= 4; i++) {
				for (k = i; k <= 4; k++) {
					KK[B_IDX(n,w,j+i,j+k)] += d*Ke4[i][k];
					FF[B_IDX(n,w,j+i,j+k)] += f*Fe4[i][k];
				}
			}
			break;
		default:
			throw std::invalid_argument("Invalid choice of bandwidth in temperature!");
		}
	}

	// Combine the two 'stiffness' matrices to get the final LHS of
	// the system.
	for (i = 0; i < B_SIZE(n,w); i++)
		KK[i] += 2*FF[i]/dt;
	// Strip off the vectors for the top and bottom nodes.
	for (i = 0; i < w+1; i++) {
		Kt[i] = KK[B_IDX(n,w,0,i)];
		KK[B_IDX(n,w,0,i)] = (i == 0 ? 1.0 : 0.0);
	}
	for (i = n-w-1; i < n; i++) {
		Kb[i] = KK[B_IDX(n,w,i,n-1)];
		KK[B_IDX(n,w,i,n-1)] = (i == n-1 ? 1.0 : 0.0);
	}
	// Perform cholesky decompostion
	decmp_chol(n,w,KK);
	// Setup the initial gradient vector
	for (i = 0; i < n; i++)
		ng[i] = 2*nt[i]/dt;
}

void
FEMthermal::step(double tt, double tb) noexcept
{
	unsigned i, j;

	// All of the information is already in the gradient.
	memset(nt,0,n*sizeof(double));
	nt[0] = tt;
	// Make corrections for first nodes
	for (i = 1; i < w+1; i++)
		nt[i] -= tt*Kt[i];
	for (i = 1; i < n-1; i++) {
		for (j = (i > w ? i-w : 0); j < i; j++)
			nt[i] += FF[B_IDX(n,w,j,i)]*ng[j];
		for (j = i; j <= i+w && j < n; j++)
			nt[i] += FF[B_IDX(n,w,i,j)]*ng[j];
	}
	// Make corrections for last nodes
	for (i = n-w-1; i < n-1; i++)
		nt[i] -= tb*Kb[i];
	nt[n-1] = tb;
	// Solve.
	bksub_chol(n,w,KK,nt);
	// Update the gradient
	for (i = 0; i < n; i++)
		ng[i] = -ng[i] + 4*nt[i]/dt;
}

void
FEMthermal::interpolate(unsigned np, const double * pd, double * pt)
{
	unsigned i, j;

	for (i = 0, j = 0; i < np; i++) {
		if (pd[i] <= nd[0]) {
			pt[i] = nt[0];
			continue;
		}
		if (pd[i] >= nd[n-1]) {
			pt[i] = nt[n-1];
			continue;
		}
		// Find element end node depth greater than or equal to
		// current interolation depth
		while (pd[i] > nd[j])
			j += w;
		if (pd[i] == nd[j]) {
			pt[i] = nt[j];
			continue;
		}
		// Interpolate shape functions over [-1 1]
		const double z = 2*(pd[i]-nd[j-w])/(nd[j]-nd[j-w])-1;
		switch (w) {
		case 1:
			pt[i] = ((z+1)*nt[j]-(z-1)*nt[j-1])/2;
			break;
		case 2:
			pt[i] = (z-1)*z*nt[j-2]/2
				  + (1-z*z)*nt[j-1]
				  + (z+1)*z*nt[j  ]/2;
			break;
		case 4:
			pt[i] = z*(z-1)*(4*z*z-1)*nt[j-4]/6 -
				  z*(2*z-1)*(4*z*z-4)*nt[j-3]/3 +
				    (z*z-1)*(4*z*z-1)*nt[j-2]     -
				  z*(2*z+1)*(4*z*z-4)*nt[j-1]/3 +
				    z*(z+1)*(4*z*z-1)*nt[j  ]/6;
			break;
		default:
			throw std::invalid_argument("Invalid choice of bandwidth in temperature!");
		}
	}
}

} // namespace OP
