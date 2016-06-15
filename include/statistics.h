/**************************************************************************

	STATISTICS.H - Prototypes for basic statistical functions.

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
		This file contains prototypes for the basic statistical functions.

	Design:
		These are close to pure C functions.
	
	Status:
		Only the standard normal distribution is implemented.
  
	History:
		2002/01/10 - Created by Jeremy Lea <reg@openpave.org>
		2016/02/02 - Seperated from RELIABILITY.H

**************************************************************************/

#ifndef __STATISTICS_H
#define __STATISTICS_H

namespace OP {

double stdnormal_rnd();
double stdnormal_pdf(double);
double quad8_stdnormal_pdf(double, double, double);
double stdnormal_cdf(double);
double stdnormal_inv(double);

} // namespace OP

#endif // STATISTICS_H
