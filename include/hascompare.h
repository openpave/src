/**************************************************************************

	HASCOMPARE.H - Simple SFINAE class to check if we have compare()

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

	Portions Copyright (C) 2016 OpenPave.org.

	Contributor(s): Jeremy Lea <reg@openpave.org>.

	Purpose:
		This header implements templated C++ classes for SFINAE
		(Substitution failure is not an error) checking if a class has
		a compare() method that returns -1/0/1.  It is used to specialise
		various classes that need to to order comparisons.

	History:
		2016/03/08 - Created a basic implementation.

**************************************************************************/

#ifndef __HASCOMPARE_H
#define __HASCOMPARE_H

#include <type_traits>

template<typename T>
struct has_compare
{
private:
	template<typename C>
	static typename std::is_same<
		decltype(std::declval<C>().compare(std::declval<C>())),int>::type
	check(C *) { return std::true_type(); };
	template<typename>
	static std::false_type check(...) { return std::false_type(); };

	typedef decltype(check<T>(0)) type;
public:
	static const bool value = type::value;
};

#endif // HASCOMPARE_H
