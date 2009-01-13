################################################################################
#
# $OpenPave$
#
# all - Create libraries.
#       Create programs. 
#
# install - Install headers, libraries, and programs on the system.
#
# Parameters to this makefile (set these before including):
#
# a)
#	TARGETS	-- the target to create 
#			(defaults to $LIBRARY $PROGRAM)
# b)
#	DIRS	-- subdirectories for make to recurse on
#			(the 'all' rule builds $TARGETS $DIRS)
# c)
#	CSRCS   -- .c files to compile
#			(used to define $OBJS)
# d)
#	PROGRAM	-- the target program name to create from $OBJS
# e)
#	LIBRARY	-- the target library name to create from $OBJS
#
################################################################################

PWD := $(shell pwd)

.DEFAULT_GOAL := all

ifeq ($(NSDISTMODE),copy)
# copy files, but preserve source mtime
INSTALL		= $(INSTALL) -t
else
# install using relative symbolic links
INSTALL		= $(INSTALL) -R
endif

GARBAGE		+= .depend core $(wildcard core.[0-9]*)
DIST_GARBAGE	+= Makefile

#
# This makefile contains rules for building the following kinds of
# libraries:
# - LIBRARY: a static (archival) library
# - SHARED_LIBRARY: a shared (dynamic link) library
# - IMPORT_LIBRARY: an import library, used only on Windows and OS/2
#
# The names of these libraries can be generated by simply specifying
# LIBRARY_NAME.
#

ifdef LIBRARY_NAME
LIBRARY		= lib$(LIBRARY_NAME).$(LIB_SUFFIX)
ifndef ONLY_STATIC_LIB
OS_CPPFLAGS += -DBUILD_DLL
ifdef WIN32
LIBRARY		= lib$(LIBRARY_NAME)_s.$(LIB_SUFFIX)
SHARED_LIBRARY	= lib$(LIBRARY_NAME).$(DLL_SUFFIX)
IMPORT_LIBRARY	= lib$(LIBRARY_NAME).$(LIB_SUFFIX)
else
SHARED_LIBRARY	= lib$(LIBRARY_NAME).$(DLL_SUFFIX)
endif
endif

echo-lib-deps:
ifdef ONLY_STATIC_LIB
	@echo $(LIB_PATH)/$(LIBRARY)
else
ifdef WIN32
	@echo $(LIB_PATH)/$(IMPORT_LIBARY)
else
	@echo $(LIB_PATH)/$(SHARED_LIBARY)
endif
endif

echo-libs:
	@echo -n "-l$(LIBRARY_NAME) ";
endif

ifdef LIBS
_LIBS = $(shell \
	for l in $(LIBS); do \
		echo -n "-L$(MOD_DEPTH)/$$l ";\
	done; \
	for l in $(LIBS); do \
		$(MAKE) -s --no-print-directory -C $(MOD_DEPTH)/$$l echo-libs; \
	done;)

_LIB_DEPS = $(shell \
	for l in $(LIBS); do \
		$(MAKE) -s --no-print-directory -C $(MOD_DEPTH)/$$l \
			LIB_PATH=$(MOD_DEPTH)/$$l echo-lib-deps; \
	done;)
endif

ifdef PROGRAM_NAME
PROGRAM		= $(PROGRAM_NAME)$(EXE_SUFFIX)
endif

ifndef TARGETS
ifdef LIBRARY
TARGETS		+= $(LIBRARY) $(SHARED_LIBRARY) $(IMPORT_LIBRARY)
endif
ifdef PROGRAM
TARGETS		+= $(PROGRAM)
endif
endif

#
# OBJS is the list of object files.  It can be constructed by
# specifying CSRCS (list of C source files) and ASFILES (list
# of assembly language source files).
#

ifndef OBJS
OBJS	= $(strip \
	$(CSRCS:.c=.$(OBJ_SUFFIX)) \
	$(CXXSRCS:.cpp=.$(OBJ_SUFFIX)) \
	$(ASFILES:.$(ASM_SUFFIX)=.$(OBJ_SUFFIX)))
endif

ALL_TRASH	= $(TARGETS) $(OBJS) $(RES) $(GARBAGE)

ifdef DIRS
LOOP_OVER_DIRS = \
	@set -e; \
	for d in $(DIRS); do \
		$(MAKE) -C $$d $@; \
	done;
endif

CONFIG_DEPS = \
	$(topsrcdir)/build/rules.mk		\
	$(MOD_DEPTH)/config.status		\
	$(MOD_DEPTH)/build/autoconf.mk		\
	Makefile

.PHONY: all
all:: $(CONFIG_DEPS)
	+$(LOOP_OVER_DIRS)

all:: $(TARGETS)

ifneq (,$(TARGETS))
$(TARGETS): $(CONFIG_DEPS)
endif

ifneq (,$(OBJS))
$(OBJS): $(CONFIG_DEPS)
endif

.PHONY: clean
clean:: $(CONFIG_DEPS)
	rm -rf $(OBJS) $(RES) $(GARBAGE)
	+$(LOOP_OVER_DIRS)

.PHONY: realclean
realclean:: $(CONFIG_DEPS)
	rm -rf $(ALL_TRASH)
	+$(LOOP_OVER_DIRS)

.PHONY: distclean
distclean:: $(CONFIG_DEPS)
	rm -rf $(ALL_TRASH) $(DIST_GARBAGE)
	+$(LOOP_OVER_DIRS)

.PHONY: install
install:: $(RELEASE_BINS) $(RELEASE_HEADERS) $(RELEASE_LIBS)
ifdef RELEASE_BINS
	$(INSTALL) -t -m 0755 $(RELEASE_BINS) $(DESTDIR)$(bindir)
endif
ifdef RELEASE_HEADERS
	$(INSTALL) -t -m 0644 $(RELEASE_HEADERS) $(DESTDIR)$(includedir)/$(include_subdir)
endif
ifdef RELEASE_LIBS
	$(INSTALL) -t -m 0755 $(RELEASE_LIBS) $(DESTDIR)$(libdir)/$(lib_subdir)
endif
	+$(LOOP_OVER_DIRS)

.PHONY: release
release:: all
ifdef RELEASE_BINS
	@echo "Copying executable programs and scripts to release directory"
	@if test ! -d $(RELEASE_BIN_DIR); then \
		rm -rf $(RELEASE_BIN_DIR); \
		$(INSTALL) -D $(RELEASE_BIN_DIR);\
	else \
		true; \
	fi
	cp $(RELEASE_BINS) $(RELEASE_BIN_DIR)
endif
ifdef RELEASE_LIBS
	@echo "Copying libraries to release directory"
	@if test ! -d $(RELEASE_LIB_DIR); then \
		rm -rf $(RELEASE_LIB_DIR); \
		$(INSTALL) -D $(RELEASE_LIB_DIR);\
	else \
		true; \
	fi
	cp $(RELEASE_LIBS) $(RELEASE_LIB_DIR)
endif
ifdef RELEASE_HEADERS
	@echo "Copying header files to release directory"
	@if test ! -d $(RELEASE_HEADERS_DEST); then \
		rm -rf $(RELEASE_HEADERS_DEST); \
		$(INSTALL) -D $(RELEASE_HEADERS_DEST);\
	else \
		true; \
	fi
	cp $(RELEASE_HEADERS) $(RELEASE_HEADERS_DEST)
endif
	+$(LOOP_OVER_DIRS)

Makefile: $(srcdir)/Makefile.in
	$(MAKE) -C $(topsrcdir) -f openpave.mk configure

$(PROGRAM): $(OBJS) $(_LIB_DEPS)
ifdef WIN32
	@sh $(topsrcdir)/build/cygwin-wrapper \
		$(CC) $(OBJS) -Fe$@ -link $(OS_LDFLAGS) \
		$(patsubst -l%,lib%.$(LIB_SUFFIX),$(subst -L,/LIBPATH:,$(_LIBS))) \
		$(OS_LIBS)
else
ifdef CXXSRCS
	$(CXX) $(OS_LDFLAGS) $(OBJS) $(_LIBS) $(OS_LIBS) -o $@
else
	$(CC) $(OS_LDFLAGS) $(OBJS) $(_LIBS) $(OS_LIBS) -o $@
endif
ifndef BUILD_DBG
ifndef BUILD_PRF
ifdef BUILD_OPT
	$(STRIP) $@
endif
endif
endif
endif

$(LIBRARY): $(OBJS)
	@rm -f $@
ifdef WIN32
	@sh $(topsrcdir)/build/cygwin-wrapper \
		$(AR) $(ARFLAGS) -OUT:"$@" $(OBJS) $(AR_EXTRA_ARGS)
else
	$(AR) $(ARFLAGS) $(OBJS) $(AR_EXTRA_ARGS)
ifdef RANLIB
	$(RANLIB) $@
endif
endif

$(SHARED_LIBRARY): $(_LIB_DEPS) $(OBJS) $(RES)
	@rm -f $@
ifdef WIN32
	@sh $(topsrcdir)/build/cygwin-wrapper \
		$(LD) $(DSO_LDFLAGS) $(patsubst -l%,lib%.$(LIB_SUFFIX),$(subst -L,-LIBPATH:,$(_LIBS))) \
	      -OUT:"$@" $(OBJS) $(RES)
	@sh $(topsrcdir)/build/cygwin-wrapper \
		$(RANLIB) -MANIFEST $@.manifest -OUTPUTRESOURCE:"$@;2"
else
ifdef CXXSRCS
	$(CXX) $(DSO_LDFLAGS) $(_LIBS) $(OBJS) $(RES) -o $@
else
	$(CC) $(DSO_LDFLAGS) $(_LIBS) $(OBJS) $(RES) -o $@
endif
ifndef BUILD_DBG
ifndef BUILD_PRF
ifdef BUILD_OPT
	$(STRIP) $@
endif
endif
endif
endif

ifdef RC
$(RES): $(RESNAME)
# The resource compiler does not understand the -U option.
ifdef WIN32
	@sh $(topsrcdir)/build/cygwin-wrapper \
		$(RC) $(RCFLAGS) $(filter-out -U%,$(DEFINES)) $(INCLUDES) -Fo$@ $<
else
	$(RC) $(RCFLAGS) $(filter-out -U%,$(DEFINES)) $(INCLUDES:-I%=--include-dir %) -o $@ $<
endif
endif

#
# Translate source filenames to absolute paths. This is required for
# debuggers under Windows to find source files automatically.
#

ifdef WIN32
abspath = $(if $(findstring :,$(1)),$(1),$(if $(filter /%,$(1)),$(1),$(PWD)/$(1)))
endif

%.$(OBJ_SUFFIX): %.cpp
ifdef WIN32
	@echo $@ ": \\" > $@.d 
	@sh $(topsrcdir)/build/cygwin-wrapper -quiet \
		$(CXX) -EP -showIncludes $(OS_CXXFLAGS) $(DEFINES) $(INCLUDES) $(call abspath,$<) 2>&1 \
			| grep "including" | sed -e 's/.*file: //g' -e 's/^[ ]*//g' -e 's^\\^/^g' -e 's/^/\t/g' -e 's/$$/ \\/g' | grep -v " [^\]" | grep -v ";" >> $@.d 
	@sh $(topsrcdir)/build/cygwin-wrapper \
		$(CXX) -Fo$@ -c $(OS_CXXFLAGS) $(DEFINES) $(INCLUDES) $(call abspath,$<)
else
	$(CXX) -c -MMD -MP $(OS_CXXFLAGS) $(DEFINES) $(INCLUDES) $< -o $@ 
endif

%.$(OBJ_SUFFIX): %.c
ifdef WIN32
	@echo $< ": \\" > $@.d 
	@sh $(topsrcdir)/build/cygwin-wrapper -quiet \
		$(CC) -EP -showIncludes $(OS_CFLAGS) $(DEFINES) $(INCLUDES) $(call abspath,$<) 2>&1 \
			| grep "including" | sed -e 's/.*ile: //g' -e 's/^[ ]*//g' -e 's^\\^/^g' -e 's/^/\t/g' -e 's/$$/ \\/g' | grep -v " " >> $@.d 
	| grep "including" | sed -e 's/.*ile: //g' -e 's/ //g' > $($(call abspath,$@):.$(OBJ_SUFFIX)=.d)
else
	$(CC) -c -MMD -MP $(OS_CFLAGS) $(DEFINES) $(INCLUDES) $< -o $@ 
endif

%.$(OBJ_SUFFIX): %.$(ASM_SUFFIX)
	$(AS) -c $(ASFLAGS) $<-o $@ 

%.i: %.c
	$(CPP) -C $(OS_CPPFLAGS) $(DEFINES) $(INCLUDES) $< > $@

ifneq (,$(OBJS))
ifdef WIN32
-include $(OBJS:.$(OBJ_SUFFIX)=.$(OBJ_SUFFIX).d)
else
-include $(OBJS:.$(OBJ_SUFFIX)=.d)
endif
endif

.SUFFIXES:
.SUFFIXES: .$(LIB_SUFFIX) .$(DLL_SUFFIX) .$(OBJ_SUFFIX) .c .cpp .$(ASM_SUFFIX) .h .i
