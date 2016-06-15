/**************************************************************************

	RELIABILITY.CPP - Implementation for the reliability class.

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

	See RELIABILITY.H

	History:
		2002/01/10 - Created by Jeremy Lea <reg@openpave.org>

**************************************************************************/

#include <stdlib.h>
#include "mathplus.h"
#include "statistics.h"
#include "reliability.h"

namespace OP {

/*
 *	Simple constructor.  We always have an owner, and optionally a previous list member.
 */
randomvar::randomvar(reliability * o, randomvar * p) {
	owner = o;
	prev = p, next = 0;
	if (prev != 0) {
		if (prev->next != 0) {
			next = prev->next;
			next->prev = this;
		}
		prev->next = this;
	}

	stdnormal = 0.0;
	ustdnormal = 0.0;
	for (int i = 0; i < RV_DIST_PARAM; i++)
		d[i] = 0.0;
}

/*
 *	A simple destructor which does a bit of list management...
 */
randomvar::~randomvar() {
	if (prev != 0) {
		prev->next = next;
	}
	if (next != 0) {
		next->prev = prev;
	}
}

/*
 * An implementation of a normally distributed random variable.
 */
struct rv_normal : public randomvar
{
public:
	virtual distribution type() {
		return randomvar::NORMAL;
	}
	virtual double x() {
		return d[0] + d[1]*z();
	}
	virtual double mean() {
		return d[0];
	}
	virtual double stddev() {
		return d[1];
	}

protected:
	friend class reliability;

	rv_normal(reliability * o, randomvar * p,
			const double m, const double s)
		: randomvar(o,p) {
		param(0,m);
		param(1,s);
	}
	virtual ~rv_normal() {
	}
};

/*
 * An implementation of a natural log normally distributed random variable.
 */
struct rv_lognormal : public randomvar
{
public:
	virtual distribution type() {
		return randomvar::LOGNORMAL;
	}
	virtual double x() {
		return exp(d[0] + d[1]*z());
	}
	virtual double mean() {
		return exp(d[0] + d[1]*d[1]/2);
	}
	virtual double stddev() {
		return mean()*sqrt(exp(d[1]*d[1])-1);
	}

protected:
	friend class reliability;

	rv_lognormal(reliability * o, randomvar * p,
			const double m, const double s)
		: randomvar(o,p) {
		param(1,sqrt(log(1+s*s/m/m)));
		param(0,log(m) - d[1]*d[1]/2);
	}
	virtual ~rv_lognormal() {
	}
};

/*
 *	Simple constructor.  We always have an owner, and optionally a previous list member.
 */
gfunction::gfunction(reliability * o, gfunction * p) {
	owner = o;
	prev = p, next = 0;
	if (prev != 0) {
		if (prev->next != 0) {
			next = prev->next;
			next->prev = this;
		}
		prev->next = this;
	}
}

/*
 *	A simple destructor which does a bit of list management...
 */
gfunction::~gfunction() {
	if (prev != 0) {
		prev->next = next;
	}
	if (next != 0) {
		next->prev = prev;
	}
}

/*
 *	A simple construcutor for a reliability problem.
 */
reliability::reliability() {
	rv_head = 0;
	gf_head = 0;
}

/*
 *	Simple destuctor.  delete's all of our lists.
 */
reliability::~reliability() {
	if (rv_head != 0) {
		while (rv_head->next != 0) {
			delete rv_head->next;
		}
		delete rv_head;
	}
	rv_head = 0;
	if (gf_head != 0) {
		while (gf_head->next != 0) {
			delete gf_head->next;
		}
		delete gf_head;
	}
	gf_head = 0;
}

/*
 * Creates a new random variable of the type specified.
 */
randomvar *
reliability::NewRV(randomvar::distribution t, double m, double s) {
	randomvar * rv;

	rv = rv_head;
	while (rv != 0)
		rv = rv->next;
	switch (t) {
	case (randomvar::NORMAL):
		rv = new rv_normal(this, rv, m, s);
		break;
	case (randomvar::LOGNORMAL):
		rv = new rv_lognormal(this, rv, m, s);
		break;
	}
	if (!rv) {
		// XXX: Error.
		return 0;
	}
	if (rv_head == 0)
		rv_head = rv;

	return rv;
}

} // namespace OP
