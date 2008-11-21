/**************************************************************************

	TEST_FFT.CPP - A test harness for fft.h.

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

	History:
		2008/09/12 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include "event.h"
#include <stdio.h>
#include <math.h>
#include "fft.h"

double
myrand()
{
  static int pos = 0;
  static unsigned in[4] = {0,0,0,0};
  static unsigned t[16];
  unsigned x, y, sum;
  int rounds;
  
  if (!pos) {
    x = 0;
    sum = 0;
    rounds = 3;
    t[1] =t[2] =t[3] =0; t[0] =in[0];
    t[5] =t[6] =t[7] =0; t[4] =in[1];
    t[9] =t[10]=t[11]=0; t[8] =in[2];
    t[13]=t[14]=t[15]=0; t[12]=in[3];
    do {
      sum += 0x9e3779b9; y=x; y+=sum;
      x = (((x) << (5)) | ((x) >> (32 - (5)))); x^=y; y=t[0]; x+=y; y=sum; t[0] =x; y+=x;
      x = (((x) << (7)) | ((x) >> (32 - (7)))); x^=y; y=t[1]; x+=y; y=sum; t[1] =x; y+=x;
      x = (((x) << (9)) | ((x) >> (32 - (9)))); x^=y; y=t[2]; x+=y; y=sum; t[2] =x; y+=x;
      x = (((x) << (13)) | ((x) >> (32 - (13)))); x^=y; y=t[3]; x+=y; y=sum; t[3] =x; y+=x;
      x = (((x) << (5)) | ((x) >> (32 - (5)))); x^=y; y=t[4]; x+=y; y=sum; t[4] =x; y+=x;
      x = (((x) << (7)) | ((x) >> (32 - (7)))); x^=y; y=t[5]; x+=y; y=sum; t[5] =x; y+=x;
      x = (((x) << (9)) | ((x) >> (32 - (9)))); x^=y; y=t[6]; x+=y; y=sum; t[6] =x; y+=x;
      x = (((x) << (13)) | ((x) >> (32 - (13)))); x^=y; y=t[7]; x+=y; y=sum; t[7] =x; y+=x;
      x = (((x) << (5)) | ((x) >> (32 - (5)))); x^=y; y=t[8]; x+=y; y=sum; t[8] =x; y+=x;
      x = (((x) << (7)) | ((x) >> (32 - (7)))); x^=y; y=t[9]; x+=y; y=sum; t[9] =x; y+=x;
      x = (((x) << (9)) | ((x) >> (32 - (9)))); x^=y; y=t[10]; x+=y; y=sum; t[10]=x; y+=x;
      x = (((x) << (13)) | ((x) >> (32 - (13)))); x^=y; y=t[11]; x+=y; y=sum; t[11]=x; y+=x;
      x = (((x) << (5)) | ((x) >> (32 - (5)))); x^=y; y=t[12]; x+=y; y=sum; t[12]=x; y+=x;
      x = (((x) << (7)) | ((x) >> (32 - (7)))); x^=y; y=t[13]; x+=y; y=sum; t[13]=x; y+=x;
      x = (((x) << (9)) | ((x) >> (32 - (9)))); x^=y; y=t[14]; x+=y; y=sum; t[14]=x; y+=x;
      x = (((x) << (13)) | ((x) >> (32 - (13)))); x^=y; y=t[15]; x+=y;
      t[15] = x;
    } while (--rounds);
    if (!++in[0]) if (!++in[1]) if (!++in[2]) ++in[3];
    pos = 16;
  }
  return 0.0000000004656612873077392578125 * long(t[--pos] & 0x7fffffff);
}

template<unsigned N>
void
doitr8()
{
  timeme();

  double x8[N];
  double y8[N];
  double z8[N];

  for (unsigned i = 0; i < N; i++) x8[i] = myrand();
  for (unsigned i = 0; i < N; i++) y8[i] = myrand();
  for (unsigned i = 0; i < N; i++) z8[i] = 0;
  for (unsigned i = 0; i < N; i += 2) {
    for (unsigned j = 0; j < N; j += 2) {
      if ((i+j) < N) {
        z8[i+j  ] += x8[i  ] * y8[j  ];
        z8[i+j+1] += x8[i  ] * y8[j+1];
        z8[i+j+1] += x8[i+1] * y8[j  ];
        z8[i+j  ] += x8[i+1] * y8[j+1];
      } else {
        z8[i+j-N+1] += x8[i  ] * y8[j  ];
        z8[i+j-N  ] += x8[i  ] * y8[j+1];
        z8[i+j-N  ] += x8[i+1] * y8[j  ];
        z8[i+j-N+1] += x8[i+1] * y8[j+1];
      }
    }
  }
  fftr<N>(y8);
  fftr_scale<N>(y8);
  fftr<N>(x8);
  fftr_mul<N>(x8,y8);
  fftr_un<N>(x8);

  double diff, error = 0;
  for (unsigned i = 0; i < N; i++) {
    diff = x8[i] - z8[i]; error += diff * diff;
  }
  printf("%6d r8 %.30f ",N,sqrt(error/(N/2))/(N/2));
  timeme("\n");
}

template<unsigned N>
void
doitc8()
{
  timeme();
  
  complex x8[N];
  complex y8[N];
  complex z8[N];

  for (unsigned i = 0; i < N; i++) x8[i].re = myrand();
  for (unsigned i = 0; i < N; i++) x8[i].im = myrand();
  for (unsigned i = 0; i < N; i++) y8[i].re = myrand();
  for (unsigned i = 0; i < N; i++) y8[i].im = myrand();
  for (unsigned i = 0; i < N; i++) z8[i].re = 0;
  for (unsigned i = 0; i < N; i++) z8[i].im = 0;
  for (unsigned i = 0; i < N; i++) {
    for (unsigned j = 0; j < N; j++) {
      z8[(i + j) & (N - 1)].re += x8[i].re * y8[j].re;
      z8[(i + j) & (N - 1)].im += x8[i].re * y8[j].im;
      z8[(i + j) & (N - 1)].im += x8[i].im * y8[j].re;
      z8[(i + j) & (N - 1)].re -= x8[i].im * y8[j].im;
    }
  }
  fftc<N>(y8);
  fftc_scale<N>(y8);
  fftc<N>(x8);
  fftc_mul<N>(x8,y8);
  fftc_un<N>(x8);
  
  double diff, error = 0;
  for (unsigned i = 0; i < N; i++) {
    diff = x8[i].re - z8[i].re; error += diff * diff;
    diff = x8[i].im - z8[i].im; error += diff * diff;
  }
  printf("%6d c8 %.30f ",N,sqrt(error/N)/N);
  timeme("\n");
}

#ifdef NOBUILD
int
main()
{
  doitr8<2>();
  doitc8<2>();
  doitr8<4>();
  doitc8<4>();
  doitr8<8>();
  doitc8<8>();
  doitr8<16>();
  doitc8<16>();
  doitr8<32>();
  doitc8<32>();
  doitr8<64>();
  doitc8<64>();
  doitr8<128>();
  doitc8<128>();
  doitr8<256>();
  doitc8<256>();
  doitr8<512>();
  doitc8<512>();
  doitr8<1024>();
  doitc8<1024>();
  doitr8<2048>();
  doitc8<2048>();
  doitr8<4096>();
  doitc8<4096>();
  doitr8<8192>();
  doitc8<8192>();
  return 0;
}
#endif
