/**************************************************************************

	RNG.H - A random number generator

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
		This header file implements a random number generator.

	Design:
		These are based on DSFMT, which has the license below.  The SIMD
		and Altivec code have been stripped, and only one exponent remains.

	History:
		2008/12/12 - Initial check-in.

**************************************************************************/

#ifndef __RNG_H
#define __RNG_H

#include <mathplus.h>

/*

Copyright (c) 2007, 2008 Mutsuo Saito, Makoto Matsumoto and Hiroshima
University.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the name of the Hiroshima University nor the names of
      its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#if defined(_MSC_VER)
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
#else
#include <inttypes.h>
#endif

#define DSFMT_MEXP        19937
#define DSFMT_N           ((DSFMT_MEXP - 128) / 104 + 1)

#define DSFMT_LOW_MASK    0x000FFFFFFFFFFFFFULL
#define DSFMT_HIGH_CONST  0x3FF0000000000000ULL
#define DSFMT_SR          12
#define DSFMT_POS1        117
#define DSFMT_SL1         19
#define DSFMT_MSK1        0x000ffafffffffb3fULL
#define DSFMT_MSK2        0x000ffdfffc90fffdULL
#define DSFMT_FIX1        0x90014964b32f4329ULL
#define DSFMT_FIX2        0x3b8d12ac548a7c7aULL
#define DSFMT_PCV1        0x3d84e1ac0dc82880ULL
#define DSFMT_PCV2        0x0000000000000001ULL

/*
 * class rng - self contained random number generator.
 */
class rng {
public:
	// Construct and seed the RNG.
	inline rng(uint32_t seed = 12345)
	  : ni(1) {
		uint64_t inner;

		uint32_t * psfmt32 = &(status[0].u32[0]);
		psfmt32[0] = seed;
		for (unsigned i = 1; i < (DSFMT_N+1)*4; i++) {
			psfmt32[i] = 1812433253UL 
			    * (psfmt32[i-1] ^ (psfmt32[i-1] >> 30)) + i;
		}
		uint64_t * psfmt64 = &(status[0].u[0]);
		for (unsigned i = 0; i < DSFMT_N*2; i++)
			psfmt64[i] = (psfmt64[i] & DSFMT_LOW_MASK) | DSFMT_HIGH_CONST;
		inner  = (status[DSFMT_N].u[0] ^ DSFMT_FIX1) & DSFMT_PCV1;
		inner ^= (status[DSFMT_N].u[1] ^ DSFMT_FIX2) & DSFMT_PCV2;
		for (unsigned i = 32; i > 0; i >>= 1) 
			inner ^= inner >> i;
		inner &= 1;
		if (inner != 1)
			status[DSFMT_N].u[1] ^= 1;
		idx = DSFMT_N*2;
	}
	inline double c1o2() {
		double *psfmt = &(status[0].d[0]);
		if (idx >= DSFMT_N*2) {
			regen();
			idx = 0;
		}
		return psfmt[idx++];
	}
	inline double c0o1() {
		return c1o2()-1.0;
	}
	inline double o0c1() {
		return 2.0-c1o2();
	}
	inline double o0o1() {
		union {
			double   d;
			uint64_t u;
		} r;
		r.d = c1o2();
		r.u |= 1;
		return r.d-1.0;
	}
	inline double stdnormal() {
		double r[2];

		if (ni == 1) {
			r[0] = sqrt(-2*log(o0c1()));
			r[1] = M_2PI*o0c1();
			n[0] = r[0]*sin(r[1]);
			n[1] = r[0]*cos(r[1]);
			ni = 0;
		} else
			ni = 1;
		return n[ni];
	}

private:
	int idx;
	int ni;
	union {
		uint64_t u[2];
		uint32_t u32[4];
		double   d[2];
	} status[DSFMT_N+1];
	double n[2];

	inline void regen() {
		for (unsigned i = 0; i < DSFMT_N; i++) {
			uint64_t t0 = status[i].u[0];
			uint64_t t1 = status[i].u[1];
			uint64_t L0 = status[DSFMT_N].u[0];
			uint64_t L1 = status[DSFMT_N].u[1];
			status[DSFMT_N].u[0] = (t0 << DSFMT_SL1) ^ (L1 >> 32)
				^ (L1 << 32) ^ status[(i+DSFMT_POS1)%DSFMT_N].u[0];
			status[DSFMT_N].u[1] = (t1 << DSFMT_SL1) ^ (L0 >> 32)
				^ (L0 << 32) ^ status[(i+DSFMT_POS1)%DSFMT_N].u[1];
			status[i].u[0] = (status[DSFMT_N].u[0] >> DSFMT_SR)
				^ (status[DSFMT_N].u[0] & DSFMT_MSK1) ^ t0;
			status[i].u[1] = (status[DSFMT_N].u[1] >> DSFMT_SR)
				^ (status[DSFMT_N].u[1] & DSFMT_MSK2) ^ t1;
		}
	}
};

#endif // RNG_H
