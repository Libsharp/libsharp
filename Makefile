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

CYTHON_MODULES=python/libsharp/libsharp.so $(if $(MPI_CFLAGS), python/libsharp/libsharp_mpi.so)

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
	$(BINDIR)/sharp_testsuite acctest && \
	$(BINDIR)/sharp_testsuite test healpix 2048 -1 1024 -1 0 1 && \
	$(BINDIR)/sharp_testsuite test fejer1 2047 -1 -1 4096 2 1 && \
	$(BINDIR)/sharp_testsuite test gauss 2047 -1 -1 4096 0 2

perftest: compile_all
	$(BINDIR)/sharp_testsuite test healpix 2048 -1 1024 -1 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 63 -1 -1 128 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 127 -1 -1 256 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 255 -1 -1 512 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 511 -1 -1 1024 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 1023 -1 -1 2048 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 2047 -1 -1 4096 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 4095 -1 -1 8192 0 1 && \
	$(BINDIR)/sharp_testsuite test gauss 8191 -1 -1 16384 0 1

%.c: %.c.in
# Only do this if the md5sum changed, in order to avoid Python and Jinja
# dependency when not modifying the c.in file
	grep `md5sum $< | cut -d ' ' -f 1` $@ || ./runjinja.py < $< > $@

genclean:
	rm libsharp/sharp_legendre.c || exit 0

$(CYTHON_MODULES): %.so: %.pyx
ifndef PIC_CFLAGS
	$(error Python extension must be built using the --enable-pic configure option.)
endif
	cython $<
	$(CC) $(DEBUG_CFLAGS) $(OPENMP_CFLAGS) $(PIC_CFLAGS) `python-config --cflags` -I$(INCDIR) -o $(<:.pyx=.o) -c $(<:.pyx=.c)
	$(CL) -shared $(<:.pyx=.o) $(OPENMP_CFLAGS) $(CYTHON_OBJ) -L$(LIBDIR) -lsharp -lfftpack -lc_utils -L`python-config --prefix`/lib `python-config --ldflags` -o $@

python: $(all_lib) hdrcopy $(CYTHON_MODULES)

pytest: python
	cd python && nosetests --nocapture libsharp/tests/test_sht.py
