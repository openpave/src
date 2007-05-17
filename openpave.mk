#
# 	OPENPAVE.MK - Global makefile for the Openpave source
#
#	$OpenPave$
#
#	The OpenPave build system is based on the Mozilla build system and is
#	thus sibject to their license, as given below.
#
#	Purpose:
#		This file controls the entire build process.
#

# ***** BEGIN LICENSE BLOCK *****
# Version: MPL 1.1/GPL 2.0/LGPL 2.1
#
# The contents of this file are subject to the Mozilla Public License Version
# 1.1 (the "License"); you may not use this file except in compliance with
# the License. You may obtain a copy of the License at
# http://www.mozilla.org/MPL/
#
# Software distributed under the License is distributed on an "AS IS" basis,
# WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
# for the specific language governing rights and limitations under the
# License.
#
# The Original Code is mozilla.org code.
#
# The Initial Developer of the Original Code is
# Netscape Communications Corporation.
# Portions created by the Initial Developer are Copyright (C) 1998
# the Initial Developer. All Rights Reserved.
#
# Contributor(s):
#   Stephen Lamm
#   Benjamin Smedberg <bsmedberg@covad.net>
#   Chase Phillips <chase@mozilla.org>
#   Mark Mentovai <mark@moxienet.com>
#
# Alternatively, the contents of this file may be used under the terms of
# either the GNU General Public License Version 2 or later (the "GPL"), or
# the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
# in which case the provisions of the GPL or the LGPL are applicable instead
# of those above. If you wish to allow use of your version of this file only
# under the terms of either the GPL or the LGPL, and not to allow others to
# use your version of this file under the terms of the MPL, indicate your
# decision by deleting the provisions above and replace them with the notice
# and other provisions required by the GPL or the LGPL. If you do not delete
# the provisions above, a recipient may use your version of this file under
# the terms of any one of the MPL, the GPL or the LGPL.
#
# ***** END LICENSE BLOCK *****

# How to build a openpave.org application:
#
# To checkout and build a tree,
#    1. cvs co openpave/openpave.mk
#    2. cd openpave
#    3. gmake -f openpave.mk 
#

AVAILABLE_PROJECTS = \
  core \
  test \
  $(NULL)

BOOTSTRAP_core :=                                \
  openpave/.cvsignore                            \
  openpave/COPYING-ADDL-1.0                      \
  openpave/configure.in                          \
  openpave/Makefile.in                           \
  $(NULL)

ifdef RUN_AUTOCONF
BOOTSTRAP_core += openpave/configure
endif

MODULES_core :=                                  \
  openpave/build                                 \
  openpave/include                               \
  openpave/src                                   \
  $(NULL)

REQUIRES_test :=                                 \
  core                                           \
  $(NULL)

MODULES_test :=                                  \
  openpave/test                                  \
  $(NULL)

#######################################################################
# Defines
#
CVS ?= cvs
comma := ,
empty :=
space := $(empty) $(empty)

CWD := $(shell pwd)
ifneq (1,$(words $(CWD)))
$(error The openpave directory cannot be located in a path with spaces.)
endif

ifeq "$(CWD)" "/"
CWD   := /.
endif

ifneq (, $(wildcard openpave.mk))
# Ran from openpave directory
ROOTDIR   := $(shell dirname $(CWD))
TOPSRCDIR := $(CWD)
else
# Ran from openpave/.. directory (?)
ROOTDIR   := $(CWD)
TOPSRCDIR := $(CWD)/openpave
endif

# if ROOTDIR equals only drive letter (i.e. "C:"), set to "/"
DIRNAME := $(shell echo "$(ROOTDIR)" | sed -e 's/^.://')
ifeq ($(DIRNAME),)
ROOTDIR := /.
endif

AUTOCONF ?= autoconf
MKDIR ?= mkdir
SH ?= /bin/sh
MAKE ?= gmake

CONFIG_GUESS_SCRIPT := $(wildcard $(TOPSRCDIR)/build/config.guess)
ifdef CONFIG_GUESS_SCRIPT
  CONFIG_GUESS = $(shell $(CONFIG_GUESS_SCRIPT))
else
  _IS_FIRST_CHECKOUT := 1
endif

####################################
# CVS

# Add the CVS root to CVS_FLAGS if needed
CVS_ROOT_IN_TREE := $(shell cat $(TOPSRCDIR)/CVS/Root 2>/dev/null)
ifneq ($(CVS_ROOT_IN_TREE),)
ifneq ($(CVS_ROOT_IN_TREE),$(CVSROOT))
  CVS_FLAGS := -d $(CVS_ROOT_IN_TREE)
endif
endif

CVS_CO_DATE_FLAGS = $(if $(OP_CO_DATE),-D "$(OP_CO_DATE)")
CVSCO = $(CVS) $(CVS_FLAGS) co $(OP_CO_FLAGS) $(if $(OP_CO_TAG),-r $(OP_CO_TAG),-A) $(CVS_CO_DATE_FLAGS)

CVSCO_LOGFILE := $(ROOTDIR)/op-cvs.log
CVSCO_LOGFILE := $(shell echo $(CVSCO_LOGFILE) | sed s%//%/%)

####################################
# Load opconfig Options

OPCONFIG_LOADER := openpave/build/opconfig2make
OPCONFIG_FINDER := openpave/build/opconfig-find 
run_for_side_effects := \
  $(shell cd $(ROOTDIR); \
     if test "$(_IS_FIRST_CHECKOUT)"; then \
        $(CVSCO) $(OPCONFIG_FINDER) $(OPCONFIG_LOADER); \
     else true; \
     fi; \
     $(OPCONFIG_LOADER) $(TOPSRCDIR) openpave/.opconfig.mk > openpave/.opconfig.out)
include $(TOPSRCDIR)/.opconfig.mk

####################################
# Options that may come from opconfig

OP_PROJECT_LIST := $(subst $(comma), ,$(OP_CO_PROJECT))
OP_PROJECT_LIST += $(foreach project,$(OP_PROJECT_LIST),$(REQUIRES_$(project)))

ifneq (,$(filter-out $(AVAILABLE_PROJECTS),$(OP_PROJECT_LIST)))
$(error OP_CO_PROJECT contains an unrecognized project.)
endif

ifeq (all,$(filter all,$(OP_PROJECT_LIST)))
  OP_PROJECT_LIST := $(AVAILABLE_PROJECTS)
endif

OP_MODULE_LIST := $(subst $(comma), ,$(OP_CO_MODULE)) $(foreach project,$(OP_PROJECT_LIST),$(MODULES_$(project)))

OPCONFIG_CONFIG += $(foreach project,$(OP_PROJECT_LIST),$(BOOTSTRAP_$(project)))

# Using $(sort) here because it also removes duplicate entries.
OP_MODULE_LIST := $(sort $(OP_MODULE_LIST))
OPCONFIG_CONFIG := $(sort $(OPCONFIG_CONFIG))

# OP_CVS_FLAGS - Basic CVS flags
ifeq "$(origin OP_CVS_FLAGS)" "undefined"
  CVS_FLAGS := $(CVS_FLAGS) -q -z 3 
else
  CVS_FLAGS := $(OP_CVS_FLAGS)
endif

ifdef OP_BUILD_PROJECTS

ifndef OP_OBJDIR
  $(error When OP_BUILD_PROJECTS is set, you must set OP_OBJDIR)
endif
ifdef OP_CURRENT_PROJECT
  OBJDIR = $(OP_OBJDIR)/$(OP_CURRENT_PROJECT)
  OP_MAKE = $(MAKE) $(OP_MAKE_FLAGS) -C $(OBJDIR)
  BUILD_PROJECT_ARG = OP_BUILD_APP=$(OP_CURRENT_PROJECT)
else
  OBJDIR = $(error Cannot find the OBJDIR when OP_CURRENT_PROJECT is not set.)
  OP_MAKE = $(error Cannot build in the OBJDIR when OP_CURRENT_PROJECT is not set.)
endif

else # OP_BUILD_PROJECTS

ifdef OP_OBJDIR
  OBJDIR = $(OP_OBJDIR)
  OP_MAKE = $(MAKE) $(OP_MAKE_FLAGS) -C $(OBJDIR)
else
  OBJDIR := $(TOPSRCDIR)
  OP_MAKE := $(MAKE) $(OP_MAKE_FLAGS)
endif

endif # OP_BUILD_PROJECTS

###################################
# Checkout main modules
#

# sort is used to remove duplicates.
OP_MODULE_LIST := $(sort $(OP_MODULE_LIST))

MODULES_CO_FLAGS := -P
ifdef OP_CO_FLAGS
  MODULES_CO_FLAGS := $(OP_CO_FLAGS)
endif
MODULES_CO_FLAGS := $(MODULES_CO_FLAGS) $(if $(OP_CO_TAG),-r $(OP_CO_TAG),-A)

ifeq (,$(strip $(OP_MODULE_LIST)))
CHECKOUT_MODULES   = $(error No modules or projects were specified. Use OP_CO_PROJECT to specify a project for checkout.)
else
CHECKOUT_MODULES   := cvs_co $(CVS) $(CVS_FLAGS) co $(MODULES_CO_FLAGS) $(CVS_CO_DATE_FLAGS) $(OP_MODULE_LIST);
endif

#######################################################################
# Rules
# 

# Print out any options loaded from opconfig.
all build checkout clean depend distclean export libs install realclean::
	@if test -f .opconfig.out; then \
	  cat .opconfig.out; \
	  rm -f .opconfig.out; \
	else true; \
	fi

ifdef _IS_FIRST_CHECKOUT
all:: checkout build
else
all:: checkout alldep
endif

# Do everything from scratch
everything: checkout clean build

####################################
# CVS checkout
#
checkout::
#	@: Backup the last checkout log.
	@if test -f $(CVSCO_LOGFILE) ; then \
	  mv $(CVSCO_LOGFILE) $(CVSCO_LOGFILE).old; \
	else true; \
	fi
	@echo "checkout start: "`date` | tee $(CVSCO_LOGFILE)
	@echo '$(CVSCO) openpave/openpave.mk $(OPCONFIG_CONFIG)'; \
        cd $(ROOTDIR) && \
	$(CVSCO) openpave/openpave.mk $(OPCONFIG_CONFIG)
	@cd $(ROOTDIR) && $(MAKE) -f openpave/openpave.mk real_checkout

#	Start the checkout. Split the output to the tty and a log file.

real_checkout:
	@set -e; \
	cvs_co() { set -e; echo "$$@" ; \
	  "$$@" 2>&1 | tee -a $(CVSCO_LOGFILE); }; \
	$(CHECKOUT_MODULES)
	@echo "checkout finish: "`date` | tee -a $(CVSCO_LOGFILE)
#	@: Check the log for conflicts. ;
	@conflicts=`egrep "^C " $(CVSCO_LOGFILE)` ;\
	if test "$$conflicts" ; then \
	  echo "$(MAKE): *** Conflicts during checkout." ;\
	  echo "$$conflicts" ;\
	  echo "$(MAKE): Refer to $(CVSCO_LOGFILE) for full log." ;\
	  false; \
	else true; \
	fi

#####################################################
# First Checkout

ifdef _IS_FIRST_CHECKOUT
# First time, do build target in a new process to pick up new files.
build::
	@cd $(ROOTDIR) && $(MAKE) -f openpave/openpave.mk build
else

#####################################################
# After First Checkout

# If we're building multiple projects, but haven't specified which project,
# loop through them.

ifeq (,$(OP_CURRENT_PROJECT)$(if $(OP_BUILD_PROJECTS),,1))
configure depend build install export libs clean realclean distclean alldep::
	set -e; \
	for app in $(OP_BUILD_PROJECTS); do \
	  @cd $(ROOTDIR) && $(MAKE) -f openpave/openpave.mk $@ OP_CURRENT_PROJECT=$$app; \
	done

else

# OP_CURRENT_PROJECT: either doing a single-project build, or building an
# individual project in a multi-project build.

####################################
# Configure

CONFIG_STATUS = $(wildcard $(OBJDIR)/config.status)
CONFIG_CACHE  = $(wildcard $(OBJDIR)/config.cache)

EXTRA_CONFIG_DEPS := \
	$(NULL)

ifdef RUN_AUTOCONF
$(TOPSRCDIR)/configure: $(TOPSRCDIR)/configure.in $(EXTRA_CONFIG_DEPS)
	@echo Generating $@ using autoconf
	cd $(TOPSRCDIR); $(AUTOCONF)
endif

CONFIG_STATUS_DEPS := \
	$(TOPSRCDIR)/openpave.mk \
	$(TOPSRCDIR)/configure.in \
	$(TOPSRCDIR)/configure \
	$(TOPSRCDIR)/.opconfig.mk \
	$(NULL)

CONFIG_STATUS_OUTS := \
	$(OBJDIR)/config.status \
	$(OBJDIR)/Makefile \
	$(OBJDIR)/build/autoconf.mk \
	$(NULL)

# configure uses the program name to determine @srcdir@. Calling it without
#   $(TOPSRCDIR) will set @srcdir@ to "."; otherwise, it is set to the full
#   path of $(TOPSRCDIR).
ifeq ($(TOPSRCDIR),$(OBJDIR))
  CONFIGURE = ./configure
else
  CONFIGURE = $(TOPSRCDIR)/configure
endif

configure:: $(TOPSRCDIR)/configure.in $(CONFIG_STATUS_DEPS)
ifdef OP_BUILD_PROJECTS
	@if test ! -d $(OP_OBJDIR); then $(MKDIR) $(OP_OBJDIR); else true; fi
endif
	@if test ! -d $(OBJDIR); then $(MKDIR) $(OBJDIR); else true; fi
	@echo cd $(OBJDIR);
	@echo $(CONFIGURE) $(CONFIGURE_ARGS)
	@cd $(OBJDIR) && $(BUILD_PROJECT_ARG) $(CONFIGURE_ENV_ARGS) $(CONFIGURE) $(CONFIGURE_ARGS) \
	  || ( echo "*** Fix above errors and then restart with\
               \"$(MAKE) -f openpave/openpave.mk build\"" && exit 1 )
	@touch $(OBJDIR)/Makefile

$(CONFIG_STATUS_OUTS): $(CONFIG_STATUS_DEPS)
	@cd $(ROOTDIR) && $(MAKE) -f openpave/openpave.mk configure

####################################
# Depend

depend:: $(CONFIG_STATUS_OUTS)
	$(OP_MAKE) export && $(OP_MAKE) depend

####################################
# Build it

build::  $(CONFIG_STATUS_OUTS)
	$(OP_MAKE)

####################################
# Other targets

# Pass these target onto the real build system
install export libs clean realclean distclean alldep:: $(OBJDIR)/Makefile
	$(OP_MAKE) $@

endif # OP_CURRENT_PROJECT
# (! IS_FIRST_CHECKOUT)
endif

cleansrcdir:
	@cd $(TOPSRCDIR); \
	if [ -f Makefile ]; then \
	  $(MAKE) distclean ; \
	else \
	  echo "Removing object files from srcdir..."; \
	  rm -fr `find . -type d \( -name .deps -print -o -name CVS \
	          -o -exec test ! -d {}/CVS \; \) -prune \
	          -o \( -name '*.[ao]' -o -name '*.so' \) -type f -print`; \
	fi; \
	rm -f config-defs.h \
	      config.cache \
	      config.log \
	      config.status \
	      .opconfig.out \
	      .opconfig.mk \
	;


.PHONY: checkout real_checkout depend build export libs alldep install clean realclean distclean cleansrcdir everything configure
