/**************************************************************************

	RELIABILIY.H - Interface for the reliability class.

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
		The reliability class serves two roles.  The first is to act
		as a factory class for random variables, and the second is to
		perform reliability analysis with these rv's.

	Design:
		The general procedure for performing reliability analysis is
		to instantiate a number of random variables, within a
		reliability problem, then set their distribution parameters. 
		Once these have been set, then any cross correllations
		between the variables must be set.

		Once the random variables are defined, then failure criteria
		are added, which use the random variables to determine any
		failure criteria.  Several gfunction's can be added, and they can
		return multiple results (to run several problems in
		parallel).  Each gfunction can also return two variables (one
		load, one capacity) for statistical analysis on a per
		function basis.

		Once this has been done, the reliability calculation is
		performed, and the results can be accessed.  The type of
		analysis depends on the type of reliability class created.
	
	Status:
		The current design is limited to doing normal and log normal
		variables, and to doing adaptive importance sampling.  The
		analysis cannot be tuned.
  
	History:
		2002/01/10 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#ifndef __RELIABILITY_H
#define __RELIABILITY_H

double stdnormal_rnd();
double stdnormal_pdf(double);
double quad8_stdnormal_pdf(double, double, double);
double stdnormal_cdf(double);
double stdnormal_inv(double);

class reliability;

struct randomvar  
{
#define RV_DIST_PARAM		4			// # of distribution parameters
public:
	enum distribution {					// Possible types of distribution
		NORMAL,
		LOGNORMAL
	};

	virtual distribution type() = 0;	// Distribution type.
	virtual double x() = 0;				// Current value.
	virtual double mean() = 0;			// Mean value.
	virtual double stddev() = 0;		// Standard Deviation.
	double z() {				// Std correllated value.
		return stdnormal;
	}
	double u() {				// Std uncorrellated value.
		return ustdnormal;
	}
	double param(int i) {	 	// Get distribution parmater i. 
		return (i < 0 || i > (RV_DIST_PARAM-1) ? 0.0 : d[i]);
	}

protected:
	// The reliability class does most of the work for us...
	friend class reliability;

	double d[RV_DIST_PARAM];			// Four distribution parmeters.
	double stdnormal;					// The current std normal value.
	double ustdnormal;					// Uncorrelated std normal value.

	double param(int i, double dp) { // Set distribution parmater i. 
		return (i < 0 || i >  (RV_DIST_PARAM-1) ? 0.0 : d[i] = dp);
	}
	randomvar(reliability * /* owner */, randomvar * /* prev */);
	virtual ~randomvar();

private:
	// List management stuff.
	reliability * owner;
	randomvar * prev;
	randomvar * next;
};

struct gfunction  
{
	virtual double g() = 0;				// Current value.

protected:
	// The reliability class does most of the work for us...
	friend class reliability;

	gfunction(reliability * owner, gfunction * prev);
	virtual ~gfunction();

private:
	// List management stuff.
	reliability * owner;
	gfunction * prev;
	gfunction * next;
};


/*
 * This is a class to perform a reliability calculation.
 */
class reliability  
{
private:
	randomvar * rv_head;		// The list of rv's.
	gfunction * gf_head;		// The list of rv's.

public:
	randomvar * NewRV(randomvar::distribution type,
		double mean = 0.0, double stddev = 1.0);
	void AddGFunc(gfunction * gfunction);

	reliability();
	~reliability();
};

#endif // RELIABILITY_H
