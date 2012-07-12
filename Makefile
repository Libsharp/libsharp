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

compile_all: $(all_cbin) hdrcopy

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

perftest: compile_all
	$(BINDIR)/sharp_test healpix 2048 1024 0 0 1 && \
	$(BINDIR)/sharp_test gauss 63 128 0 0 1 && \
	$(BINDIR)/sharp_test gauss 127 256 0 0 1 && \
	$(BINDIR)/sharp_test gauss 255 512 0 0 1 && \
	$(BINDIR)/sharp_test gauss 511 1024 0 0 1 && \
	$(BINDIR)/sharp_test gauss 1023 2048 0 0 1 && \
	$(BINDIR)/sharp_test gauss 2047 4096 0 0 1 && \
	$(BINDIR)/sharp_test gauss 4095 8192 0 0 1 && \
	$(BINDIR)/sharp_test gauss 8191 16384 0 0 1
