/**************************************************************************

	AUTODELETE.H - Simple RAII smart pointers

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

	Portions Copyright (C) 2014 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements templated C++ classes for RAII (Resource
		Allocation Is Initialization), without the bulk of the STL's
		implementation.  These are not reference counted or container safe,
		they are only intended to free allocated memory when throwing.

	History:
		2014/02/06 - Created a basic implementation.

**************************************************************************/

#ifndef __AUTODELETE_H
#define __AUTODELETE_H

template<typename T>
class autodelete
{
public:
	autodelete(T * p)
		: ptr(p) {
	}
	~autodelete() {
		delete [] ptr;
	}
	T & operator[] (const int i) {
		return ptr[i];
	}
	T & operator[] (const unsigned i) {
		return ptr[i];
	}
	T & operator* () {
		return *ptr;
	}
	T * operator-> () {
		return ptr;
	}
	operator T* () {
		return ptr;
	}
	operator void* () {
		return ptr;
	}

private:
	T * ptr;
};

#endif // AUTODELETE_H
