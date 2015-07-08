#
#    OPENPAVE.MK - Global makefile for the Openpave source
#
#    The OpenPave build system is based on the Mozilla build system and is
#    thus sibject to their license, as given below.
#
#    Purpose:
#        This file controls the entire build process.
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
  core  \
  test  \
  libop \
  fem3d \
  fcgi  \
  $(NULL)

BOOTSTRAP_core :=                       \
  .cvsignore                            \
  COPYING-ADDL-1.0                      \
  aclocal.m4                            \
  configure.in                          \
  Makefile.in                           \
  build                                 \
  $(NULL)

ifndef OP_AUTOCONF
BOOTSTRAP_core += configure
endif

MODULES_core :=                         \
  include                               \
  src                                   \
  $(NULL)

define default_mods =
ifndef MODULES_$(1)
	MODULES_$(1) := $(1)
endif
ifndef REQUIRES_$(1)
	REQUIRES_$(1) := core
endif
endef
$(foreach p,$(filter-out core,$(AVAILABLE_PROJECTS)),$(eval $(call default_mods,$(p))))

#######################################################################
# Defines
#
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
CVS ?= cvs
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

# OP_CVS_ARGS - Basic CVS flags
ifdef OP_CVS_ARGS
  CVS_ARGS := $(OP_CVS_ARGS)
else
# Add the CVS root to CVS_ARGS if needed
CVS_ROOT_IN_TREE := $(shell cat $(TOPSRCDIR)/CVS/Root 2>/dev/null)
ifneq ($(CVS_ROOT_IN_TREE),)
  CVS_ARGS := -d $(CVS_ROOT_IN_TREE)
  ifneq ($(subst :ext:,,$(CVS_ROOT_IN_TREE)),$(CVS_ROOT_IN_TREE))
    OP_CVS_USER := $(subst @cvs.openpave.org:/home/cvs,,$(subst :ext:,,$(CVS_ROOT_IN_TREE)))
  endif
endif
  CVS_ARGS := $(CVS_ARGS) -q -z 3 
endif

# OP_CO_ARGS - Checkout flags
ifdef OP_CO_ARGS
  CVSCO_ARGS := $(OP_CO_ARGS)
else
  CVSCO_ARGS := -P
endif
CVSCO_ARGS := $(CVSCO_ARGS) $(if $(OP_CO_DATE),-D "$(OP_CO_DATE)")
CVSCO_ARGS := $(CVSCO_ARGS) $(if $(OP_CO_TAG),-r $(OP_CO_TAG),-A)
CVSCO = $(CVS) $(CVS_ARGS) co $(CVSCO_ARGS)

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
OPCONFIG := $(shell cd $(ROOTDIR) && $(OPCONFIG_FINDER) $(TOPSRCDIR))

####################################
# Options that may come from opconfig

OP_PROJECT_LIST := $(subst $(comma), ,$(OP_PROJECTS))
ifeq ($(OP_PROJECT_LIST),)
	OP_PROJECT_LIST := all
endif
ifeq (all,$(filter all,$(OP_PROJECT_LIST)))
	OP_PROJECT_LIST := $(AVAILABLE_PROJECTS)
endif
ifneq (,$(filter-out $(AVAILABLE_PROJECTS),$(OP_PROJECT_LIST)))
$(error OP_PROJECTS contains an unrecognized project.  Options are $(strip $(AVAILABLE_PROJECTS)))
endif
OP_PROJECT_TMP := $(empty)
define requires =
	OP_PROJECT_TMP += $(foreach r,$(REQUIRES_$(1)),$(eval $(call requires,$(r)))) $(filter-out $(OP_PROJECT_TMP),$(1))
endef
$(foreach project,$(OP_PROJECT_LIST),$(eval $(call requires,$(project))))
OP_PROJECT_LIST := $(strip $(OP_PROJECT_TMP))

OP_MODULE_LIST := $(strip $(foreach project,$(OP_PROJECT_LIST),$(MODULES_$(project))))
OP_BOOTSTRAP_LIST := $(strip $(foreach project,$(OP_PROJECT_LIST),$(BOOTSTRAP_$(project))))

OP_PROJECT_TMP := $(empty)
define requires_mod =
	OP_MODULE_REQUIRES_$(1) := $(strip $(foreach project,$(OP_PROJECT_REQUIRES_$(2)),$(filter-out $(1),$(MODULES_$(project)))))
endef
define requires_proj =
	$(eval $(call requires,$(1)))
	$(eval OP_PROJECT_REQUIRES_$(1) := $(strip $(foreach project,$(OP_PROJECT_TMP),$(project))))
	OP_PROJECT_TMP := $(empty)
	$(foreach module,$(MODULES_$(1)),$(eval $(call requires_mod,$(module),$(1))))
endef
$(foreach project,$(OP_PROJECT_LIST),$(eval $(call requires_proj,$(project))))
OP_MODULE_REQUIRES := $(empty)
define requires_list =
	OP_MODULE_REQUIRES += OP_MODULE_REQUIRES_$(1)=$(subst $(space),:,$(OP_MODULE_REQUIRES_$(1)))
endef
$(foreach module,$(OP_MODULE_LIST),$(eval $(call requires_list,$(module))))

ifdef OP_OBJDIR
  OBJDIR = $(OP_OBJDIR)
  OP_MAKE = $(MAKE) $(OP_MAKE_ARGS) -C $(OBJDIR)
else
  OBJDIR := $(TOPSRCDIR)
  OP_MAKE := $(MAKE) $(OP_MAKE_ARGS)
endif

#######################################################################
# Rules
# 

# Print out any options loaded from opconfig.
all build scan-build checkout clean depend distclean install realclean::
	@if test -f .opconfig.out; then \
	  cat .opconfig.out; \
	  rm -f .opconfig.out; \
	else true; \
	fi

.DEFAULT_GOAL := all

.PHONY: all
all:: checkout build

# Do everything from scratch
.PHONY: everything
everything: checkout clean build

####################################
# CVS checkout
#
.PHONY: checkout
checkout::
	@set -e; \
	echo "*** CVS checkout started: "`date` | tee $(CVSCO_LOGFILE); \
	echo '*** Checking out bootstrap files...'; \
	cd $(ROOTDIR) && \
	    $(CVSCO) openpave/openpave.mk $(addprefix openpave/,$(sort $(OP_BOOTSTRAP_LIST))) 2>&1 | tee -a $(CVSCO_LOGFILE); \
	cd $(ROOTDIR) && \
	    $(MAKE) $(OP_MAKE_ARGS) -f openpave/openpave.mk real_checkout

#    Start the checkout. Split the output to the tty and a log file.

###################################
# Checkout main modules
#

.PHONY: real_checkout
ifeq (,$(strip $(OP_MODULE_LIST)))
real_checkout:
	$(error No modules or projects were specified. Use OP_PROJECTS to specify a project for checkout.)
else
real_checkout:
	@echo '*** Checking out project files...'; \
	    $(CVSCO) $(addprefix openpave/,$(sort $(OP_MODULE_LIST))) 2>&1 | tee -a $(CVSCO_LOGFILE); \
	echo "*** CVS checkout finished: "`date` | tee -a $(CVSCO_LOGFILE); \
	conflicts=`egrep "^C " $(CVSCO_LOGFILE)`; \
	if test "$$conflicts"; then \
	  echo "*** Conflicts during checkout."; \
	  echo "$$conflicts"; \
	  echo "*** Refer to $(CVSCO_LOGFILE) for full log."; \
	  false; \
	else true; \
	fi
endif

#####################################################
# First Checkout

ifdef _IS_FIRST_CHECKOUT
# First time, do build target in a new process to pick up new files.
.PHONY: build
build::
	@cd $(TOPSRCDIR) && \
	    $(MAKE) $(OP_MAKE_ARGS) -f openpave.mk build

.PHONY: scan-build
scan-build::
	@cd $(TOPSRCDIR) && \
	    $(MAKE) $(OP_MAKE_ARGS) -f openpave.mk scan-build
else

#####################################################
# After First Checkout

####################################
# Configure

EXTRA_CONFIG_DEPS := \
	$(TOPSRCDIR)/aclocal.m4 \
	$(NULL)

ifdef OP_AUTOCONF
$(TOPSRCDIR)/configure: $(TOPSRCDIR)/configure.in $(EXTRA_CONFIG_DEPS)
	@echo Generating $@ using $(AUTOCONF)
	cd $(TOPSRCDIR); $(AUTOCONF)
endif

CONFIG_STATUS_DEPS := \
	$(TOPSRCDIR)/openpave.mk \
	$(TOPSRCDIR)/configure.in \
	$(TOPSRCDIR)/configure \
	$(TOPSRCDIR)/.opconfig.mk \
	$(TOPSRCDIR)/build/autoconf.mk.in \
	$(OPCONFIG)

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

# OP_CONFIGURE_ENV - Configure environment
ifdef OP_CONFIGURE_ENV
  CONFIGURE_ENV := $(OP_CONFIGURE_ENV)
else
  CONFIGURE_ENV := $(NULL)
endif
CONFIGURE_ENV   += OP_PROJECTS="$(OP_PROJECT_LIST)"
CONFIGURE_ENV   += OP_MODULES="$(OP_MODULE_LIST)"
CONFIGURE_ENV   += OP_MODULE_REQUIRES="$(OP_MODULE_REQUIRES)"
CONFIGURE_ENV   += OBJDIR="$(OBJDIR)"
ifdef OP_CVS_USER
  CONFIGURE_ENV += OP_CVS_USER="$(OP_CVS_USER)"
endif

# OP_CONFIGURE_ARGS - Configure arguments
ifdef OP_CONFIGURE_ARGS
  CONFIGURE_ARGS := $(OP_CONFIGURE_ARGS)
else
  CONFIGURE_ARGS := $(NULL)
endif

.PHONY: configure
configure:: $(CONFIG_STATUS_DEPS)
	@if test ! -d $(OBJDIR); then $(MKDIR) $(OBJDIR); else true; fi; \
	echo '*** Running configure...' && \
	cd $(OBJDIR) && $(CONFIGURE_ENV) $(CONFIGURE) $(CONFIGURE_ARGS) \
	  || ( echo "*** Fix above errors and then restart with\
	           \"$(value MAKE) -f openpave.mk build\"" && exit 1 )

$(CONFIG_STATUS_OUTS): $(CONFIG_STATUS_DEPS)
	@cd $(TOPSRCDIR) && \
	    $(MAKE) $(OP_MAKE_ARGS) -f openpave.mk configure

####################################
# Build it

.PHONY: build
build:: $(CONFIG_STATUS_OUTS)
	@echo '*** Building...' && \
	    $(OP_MAKE)

.PHONY: scan-build
scan-build:: $(CONFIG_STATUS_OUTS)
	@echo '*** Building...' && \
	    scan-build34 -v $(OP_MAKE)

# Pass these target onto the real build system
.PHONY: depend install install clean realclean distclean
depend install clean realclean distclean:: $(CONFIG_STATUS_OUTS)
	$(OP_MAKE) $@

endif # (! IS_FIRST_CHECKOUT)

distclean::
	@cd $(TOPSRCDIR) && \
	rm -f config-defs.h \
	      config.cache \
	      config.log \
	      config.status \
	      .opconfig.out \
	      .opconfig.mk \
	;

.SUFFIXES:
