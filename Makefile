SHARP_TARGET?=auto
ifndef SHARP_TARGET
  SHARP_TARGET:=$(error SHARP_TARGET undefined. Please see README.compilation for help)UNDEFINED
endif

default: compile_all
SRCROOT:=$(shell pwd)
include $(SRCROOT)/config/config.$(SHARP_TARGET)
include $(SRCROOT)/config/rules.common

all_hdr:=
all_lib:=
all_cbin:=

FULL_INCLUDE:=

include c_utils/planck.make
include libfftpack/planck.make
include libsharp/planck.make
include docsrc/planck.make

$(all_lib): %: | $(LIBDIR)_mkdir
	@echo "#  creating library $*"
	$(ARCREATE) $@ $^

$(all_cbin): %: | $(BINDIR)_mkdir
	@echo "#  linking C binary $*"
	$(CL) -o $@ $^ $(CLFLAGS)
#	$(CXX) -o $@ $^ $(CLFLAGS)

compile_all: $(all_cbin) hdrcopy

autotune: sharp_bench
	$(BINDIR)/sharp_bench
	mv sharp_oracle.inc $(SRCROOT)/libsharp
	$(MAKE)

hdrclean:
	@if [ -d $(INCDIR) ]; then rm -rf $(INCDIR)/* ; fi

hdrcopy: | $(INCDIR)_mkdir
	@if [ "$(all_hdr)" ]; then cp -p $(all_hdr) $(INCDIR); fi

$(notdir $(all_cbin)) : % : $(BINDIR)/%

test: compile_all
	$(BINDIR)/sharp_acctest && \
	$(BINDIR)/sharp_test healpix 2048 1024 1 0 1 && \
	$(BINDIR)/sharp_test ecp 2047 4096 0 2 1 && \
	$(BINDIR)/sharp_test gauss 2047 4096 0 0 2
