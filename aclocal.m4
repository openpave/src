dnl
dnl $OpenPave$
dnl 
dnl READ_OPCONFIG() - Read in 'opconfig' file
AC_DEFUN([READ_OPCONFIG],
[AC_REQUIRE([AC_INIT_BINSH])dnl
# Read in 'opconfig' script to set the initial options.
# See the mozconfig2configure script for more details.
_AUTOCONF_TOOLS_DIR=`dirname [$]0`/[$1]/build
. $_AUTOCONF_TOOLS_DIR/opconfig2configure])

dnl This gets inserted at the top of the configure script
READ_OPCONFIG(.)
