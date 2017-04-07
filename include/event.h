/**************************************************************************

	EVENT.H - A set of callbacks for handling events.

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

#if defined(_EVENT_IMP) || defined(_PROGRESS_IMP)
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#if defined(DARWIN)
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif
#include <exception>
#include <stdexcept>
#include <new>

namespace OP {

/* Error levels */
#define EVENT_FATAL             0       /* Fatal error condition */
#define EVENT_ERROR             1       /* Non-fatal error */
#define EVENT_WARN              2       /* Non-ignorable Warning */
#define EVENT_NAG               3       /* Ingorable Warning */
#define EVENT_NOTE              4       /* Important message */
#define EVENT_MSG               5       /* General message */
#define EVENT_DEBUG             6       /* Debugging message */

#define EVENT_PROGRESS_START    0   /* Start a progress bar */
#define EVENT_PROGRESS_UPDATE   1   /* Update progress */
#define EVENT_PROGRESS_STOP     2   /* Finish a progress bar */

extern void event_msg(int level, const char *fmt, ...) OP_PRINTF(2,3);
extern void event_progress_bar(int level, double p,
                               const char *fmt, ...) OP_PRINTF(3,4);
extern void event_progress(int type, int marker,
                           const char *fmt, ...) OP_PRINTF(3,4);
extern void timeme(const char * msg = nullptr);

#ifdef _EVENT_IMP

/*
 * Default error event handler.
 */
void
event_msg(int level, const char * fmt, ...)
{
	va_list args;

	va_start(args,fmt);
	if (level < EVENT_DEBUG && fmt != nullptr) {
		vfprintf(stderr,fmt,args);
		fprintf(stderr,"\n");
	}
	va_end(args);
}

void
timeme(const char * msg)
{
#if !defined(_MSC_VER) && !defined(DARWIN)
#if !defined(__FreeBSD__)
#define CLOCK_PROF CLOCK_PROCESS_CPUTIME_ID
#endif
	static struct timespec start;
	struct timespec stop;
	double run_time;
	if (msg == nullptr) {
		clock_gettime(CLOCK_PROF,&start);
		return;
	}
	clock_gettime(CLOCK_PROF,&stop);
	run_time = double(stop.tv_sec - start.tv_sec) + double(stop.tv_nsec - start.tv_nsec) / 1000000000.0;
#elif defined(DARWIN)
	static struct timeval start;
	struct timeval stop;
	double run_time;
	if (msg == nullptr) {
		gettimeofday(&start,NULL);
		return;
	}
	gettimeofday(&stop,NULL);
	run_time = double(stop.tv_sec - start.tv_sec) + double(stop.tv_usec - start.tv_usec) / 1000000.0;
#elif defined(_MSC_VER)
	static clock_t start;
	clock_t stop;
	double run_time;
	if (msg == nullptr) {
		start = clock();
		return;
	}
	stop = clock();
	run_time = double(stop - start) / CLOCKS_PER_SEC;
#endif
	printf(" %f%s",run_time,msg);
}

#undef _EVENT_IMP
#endif

#ifdef _PROGRESS_IMP

/*
 * Default progress bar handler.
 */
void
event_progress_bar(int level, double p, const char * fmt, ...)
{
	va_list args;
	static char buf[7] = "";

	va_start(args,fmt);
	if (strlen(buf) > 0) {
		memset(buf,'\b',strlen(buf));
		fprintf(stderr,"%s",buf);
		buf[0] = '\0';
	}
	if (level == 0 && fmt != nullptr)
		vfprintf(stderr,fmt,args);
	sprintf(buf,"(%3.0f%%)",(p<0.0?0.0:(p>1.0?1.0:p))*100);
	fprintf(stderr,"%s",buf);
	va_end(args);
}

/*
 * Default error event handler.
 */
void
event_progress(int type, int marker, const char * fmt, ...)
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
		if (level < 0) {
			throw std::logic_error("Event update called on unstarted progress indicator!");
		} else if (level < 3) {
			double p = 0.0, m = 1.0;
			mark[level] = marker;
			for (int i = 0; level >= 0 && i < level; i++) {
				p += m*double(mark[i])/double(max[i]);
				m *= 1.0/double(max[i]);
			}
			p += m*double(mark[level])/double(max[level]);
			event_progress_bar(level,p,fmt,args);
		}
		break;
	case EVENT_PROGRESS_STOP:
		event_progress_bar(level,1.0,fmt,args);
		if (level < 3)
			max[level] = 0, mark[level] = 0;
		if (level == 0)
			fprintf(stderr,"\n");
		level--;
		break;
	default:
		throw std::invalid_argument("Invalid progress indication!");
	}
	va_end(args);
}

#undef _PROGRESS_IMP
#endif

} // namespace OP

#endif // EVENT_H
