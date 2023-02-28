/*************************************************************************

	WASM.CPP - WASM wrapper for OpenPave classes.

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
		2022/07/20 - Created by Jeremy Lea <reg@openpave.org>

*************************************************************************/

#include <cstdarg>
#include <memory>
#define _EVENT_IMP
#include "event.h"
#include "tree.h"
#include "pavement.h"
#include "thermal.h"
#include "emscripten/emscripten.h"
#include "emscripten/bind.h"

using namespace OP;
using namespace emscripten;

EMSCRIPTEN_BINDINGS(openpave) {
	class_<point2d>("point2d")
	.constructor<double,double>()
	.constructor<const point2d &>()
	.property("x",&point2d::x)
	.property("y",&point2d::y)
	;
	class_<paveload,emscripten::base<point2d>>("paveload")
	.constructor<const point2d &,double,double,double>()
	.constructor<const paveload &>()
	.property("force",&paveload::force)
	.property("pressure",&paveload::pressure)
	.property("radius",&paveload::radius)
	.property("area",&paveload::area)
	;
	class_<point3d,emscripten::base<point2d>>("point3d")
	.constructor<double,double,double>()
	.constructor<const point3d &>()
	.property("z",&point3d::z)
	;
	class_<pavepoint,emscripten::base<point3d>>("pavepoint")
	.constructor<double,double,double>()
	.constructor<double,double,double,unsigned>()
	.property("il",&pavepoint::il)
	;
	enum_<pavedata::type>("pavedata.type")
	.value("deflct",pavedata::deflct)
	.value("stress",pavedata::stress)
	.value("strain",pavedata::strain)
	;
	enum_<pavedata::direction>("pavedata.direction")
	.value("xx",pavedata::xx)
	.value("yy",pavedata::yy)
	.value("zz",pavedata::zz)
	.value("xy",pavedata::xy)
	.value("xz",pavedata::xz)
	.value("yz",pavedata::yz)
	.value("p1",pavedata::p1)
	.value("p2",pavedata::p2)
	.value("p3",pavedata::p3)
	.value("s1",pavedata::s1)
	.value("s2",pavedata::s2)
	.value("s3",pavedata::s3)
	;
	class_<pavedata,emscripten::base<pavepoint>>("pavedata")
	.function("result",&pavedata::result)
	;
	class_<LElayer>("LElayer")
	.property("bottom",&LElayer::bottom)
	.property("top",&LElayer::top)
	.property("thickness",
		select_overload<double() const>(&LElayer::thickness),
		select_overload<double(double)>(&LElayer::thickness))
	.property("emod",
		select_overload<double() const>(&LElayer::emod),
		select_overload<double(double)>(&LElayer::emod))
	.property("poissons",
		select_overload<double() const>(&LElayer::poissons),
		select_overload<double(double)>(&LElayer::poissons))
	.property("slip",
		select_overload<double() const>(&LElayer::slip),
		select_overload<double(double)>(&LElayer::slip))
	;
	enum_<LEsystem::resulttype>("LEsystem.resulttype")
	.value("all"      ,LEsystem::all)
	.value("fast"     ,LEsystem::fast)
	.value("dirty"    ,LEsystem::dirty)
	.value("odemark"  ,LEsystem::odemark)
	.value("fastnum"  ,LEsystem::fastnum)
	.value("accurate" ,LEsystem::accurate)
	.value("mask"     ,LEsystem::mask)
	.value("disp"     ,LEsystem::disp)
	.value("grad"     ,LEsystem::grad)
	.value("dispgrad" ,LEsystem::dispgrad)
	.value("fastdisp" ,LEsystem::fastdisp)
	.value("dirtydisp",LEsystem::dirtydisp)
	.value("fastgrad" ,LEsystem::fastgrad)
	.value("failure"  ,LEsystem::failure)
	;
	class_<LEsystem>("LEsystem")
	.constructor<>()
	.function("addlayer",optional_override(
			[](LEsystem& this_, double h, double e, double v, double s,
					unsigned p) -> LElayer * {
				return &(this_.addlayer(h,e,v,s,p));
			}
		),allow_raw_pointers())
	.function("addlayer",optional_override(
			[](LEsystem& this_, double h, double e, double v, double s)
					-> LElayer * {
				return &(this_.addlayer(h,e,v,s));
			}
		),allow_raw_pointers())
	.function("addlayer",optional_override(
			[](LEsystem& this_, double h, double e, double v) -> LElayer * {
				return &(this_.addlayer(h,e,v));
			}
		),allow_raw_pointers())
	.function("removelayer",&LEsystem::removelayer)
	.function("removelayers",&LEsystem::removelayers)
	.property("layers",&LEsystem::layers)
	.function("layer",optional_override(
			[](LEsystem& this_, unsigned i)
					-> LElayer * {
				return &(this_.layer(i));
			}
		),allow_raw_pointers())
	.property("defaultgroup",
		select_overload<unsigned() const>(&LEsystem::defaultgroup),
		select_overload<unsigned(unsigned)>(&LEsystem::defaultgroup))
	.property("groups",&LEsystem::groups)
	.function("addload",
		select_overload<void(const point2d &,double,double,double)>(
			&LEsystem::addload))
	.function("addload",
		select_overload<void(const paveload &)>(&LEsystem::addload))
	.function("addload",
		select_overload<void(unsigned,const point2d &,double,double,double)>(
			&LEsystem::addload))
	.function("addload",
		select_overload<void(unsigned,const paveload &)>(&LEsystem::addload))
	.function("removeload",
		select_overload<void(const unsigned)>(&LEsystem::removeload))
	.function("removeload",
		select_overload<void(const unsigned,const unsigned)>(
			&LEsystem::removeload))
	.function("removeloads",
		select_overload<void()>(&LEsystem::removeloads))
	.function("removeloads",
		select_overload<void(const unsigned)>(
			&LEsystem::removeloads))
	.function("loads",
		select_overload<unsigned() const>(&LEsystem::loads))
	.function("loads",
		select_overload<unsigned(const unsigned) const>(
			&LEsystem::loads))
	.function("addpoint",optional_override(
			[](LEsystem& this_, const point3d & z) -> void {
				this_.addpoint(z);
			}
		))
	.function("addpoint",&LEsystem::addpoint)
	.function("removepoint",optional_override(
			[](LEsystem& this_, const point3d & z) -> void {
				this_.removepoint(z);
			}
		))
	.function("removepoint",&LEsystem::removepoint)
	.function("removepoints",&LEsystem::removepoints)
	.function("results",&LEsystem::results)
	.function("result",optional_override(
			[](const LEsystem& this_, const point3d & p)
					-> const pavedata * {
				return &(this_.result(p));
			}
		),allow_raw_pointers())
	.function("result",optional_override(
			[](const LEsystem& this_, const point3d & p, const unsigned l)
					-> const pavedata * {
				return &(this_.result(p,l));
			}
		),allow_raw_pointers())
	.function("result",optional_override(
			[](const LEsystem& this_, const unsigned g, const point3d & p,
				const unsigned l) -> const pavedata * {
				return &(this_.result(g,p,l));
			}
		),allow_raw_pointers())
	.function("calc_accurate",&LEsystem::calc_accurate)
	.function("calculate",optional_override(
			[](LEsystem& this_, LEsystem::resulttype r) -> bool {
				return this_.calculate(r);
			}
		))
	.function("calculate",optional_override(
			[](LEsystem& this_) -> bool {
				return this_.calculate();
			}
		))
	.function("calc_odemark",&LEsystem::calc_odemark)
	.function("calc_fastnum",&LEsystem::calc_fastnum)
	;
	class_<defldata,emscripten::base<point3d>>("defldata")
	.constructor<const unsigned,const point3d &,const double>()
	.property("lg",&defldata::lg)
	.property("measured",&defldata::measured)
	.property("calculated",&defldata::calculated)
	;
	class_<LEbackcalc,emscripten::base<LEsystem>>("LEbackcalc")
	.constructor<>()
	.function("adddefl",
		select_overload<void(const point3d &,double)>(
			&LEbackcalc::adddefl))
	.function("adddefl",
		select_overload<void(const unsigned,const point3d &,const double)>(
			&LEbackcalc::adddefl))
	.function("adddefl",
		select_overload<void(const defldata &)>(
			&LEbackcalc::adddefl))
	.property("deflections",&LEbackcalc::deflections)
	.function("deflection",optional_override(
			[](LEbackcalc& this_, unsigned i)
					-> const defldata * {
				return &(this_.getdefl(i));
			}
		),allow_raw_pointers())
	.function("removedeflections",&LEbackcalc::removedeflections)
	.function("setup",&LEbackcalc::setup)
	.function("backcalc",&LEbackcalc::backcalc)
	;
	EM_ASM(
		Module['pavedata']['type'] = Module['pavedata.type'];
		delete Module['pavedata.type'];
		Module['pavedata']['direction'] = Module['pavedata.direction'];
		delete Module['pavedata.direction'];
		Module['LEsystem']['resulttype'] = Module['LEsystem.resulttype'];
		delete Module['LEsystem.resulttype'];
	);
}
