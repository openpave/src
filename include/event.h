/**************************************************************************

	EVENT.H - A set of callbacks for handling events.

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

	The Original Code is OpenPave.org Core Libraries.

	The Initial Developer of the Original Code is OpenPave.org.

	Portions Copyright (C) 2006 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This simple C header defines a global external error handler, and
		provides two simple implementations.

	Design:
		Code sometimes needs to complain.  These functions provide a way
		for that to happen.

		There are a number of defines which determine what it's
		complaining about, and also a basic DOS implementation.

	History:
		2002/10/03 - Created by Jeremy Lea <reg@openpave.org>
		2002/11/05 - Rewritten into C from C++...

**************************************************************************/

#ifndef __EVENT_H
#define __EVENT_H

/* Error levels */
#define	EVENT_FATAL	0		/* Fatal error condition */
#define	EVENT_ERROR	1		/* Non-fatal error */
#define	EVENT_WARN	2		/* Non-ignorable Warning */
#define	EVENT_NAG	3		/* Ingorable Warning */
#define	EVENT_NOTE	4		/* Important message */
#define EVENT_MSG	5		/* General message */
#define EVENT_DEBUG	6		/* Debugging message */

#define EVENT_PROGRESS_START	0	/* Start a progress bar */
#define EVENT_PROGRESS_UPDATE	1	/* Update progress */
#define EVENT_PROGRESS_STOP		2	/* Finish a progress bar */

extern void event_msg(const int level, const char *fmt,...) OP_PRINTF(2,3);
extern void event_progress_bar(const int level, const double p,
							   const char *fmt,...) OP_PRINTF(3,4);
extern void event_progress(const int type, const int marker,
						   const char *fmt,...) OP_PRINTF(3,4);

#ifdef _EVENT_IMP

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/*
 * Default error event handler.
 */
void event_msg(const int level, const char * fmt, ...)
{
	va_list args;

	va_start(args,fmt);
	if (level < EVENT_DEBUG && fmt != NULL)
		vfprintf(stderr,fmt,args);
	va_end(args);
};

#undef _EVENT_IMP
#endif

#ifdef _PROGRESS_IMP

#include <stdarg.h>
#include <stdio.h>
#include <string.h>

/*
 * Default progress bar handler.
 */
void event_progress_bar(const int level, const double p,
						const char * fmt, ...)
{
	va_list args;
	static char buf[7] = "";

	va_start(args,fmt);
	if (strlen(buf) > 0) {
		memset(buf,'\b',strlen(buf));
		event_msg(EVENT_NOTE,"%s",buf);
		buf[0] = '\0';
	}
	if (level == 0)
		event_msg(EVENT_NOTE,fmt,args);
	sprintf(buf,"(%3.0f%%)",(p<0.0?0.0:(p>1.0?1.0:p))*100);
	event_msg(EVENT_NOTE,"%s",buf);
	va_end(args);
};

/*
 * Default error event handler.
 */
void event_progress(const int type, const int marker,
					const char * fmt, ...)
{
	va_list args;
	static int level = -1;
	static int max[3], mark[3];

	va_start(args,fmt);
	switch (type) {
	case EVENT_PROGRESS_START:
		level++;
		if (level < 3)
			max[level] = marker, mark[level] = 0;
		event_progress_bar(level,0.0,fmt,args);
		break;
	case EVENT_PROGRESS_UPDATE:
		if (level < 3) {
			double p = 0.0, m = 1.0;
			mark[level] = marker;
			for (int i = 0; i <= level; i++) {
				p += m*double(mark[i])/double(max[i]);
				m *= 1.0/double(max[i]);
			}
			event_progress_bar(level,p,fmt,args);
		}
		break;
	case EVENT_PROGRESS_STOP:
		event_progress_bar(level,1.0,fmt,args);
		if (level < 3)
			max[level] = 0, mark[level] = 0;
		if (level == 0)
			event_msg(EVENT_NOTE,"\n");
		break;
		level--;
	default:
		event_msg(EVENT_ERROR,"Invalid progress indication!");
	};
	va_end(args);
};
#undef _PROGRESS_IMP
#endif

#endif // EVENT_H
