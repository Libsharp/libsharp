PKG:=libsharp

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libsharp.a
BIN:=sharp_testsuite
LIBOBJ:=sharp_ylmgen_c.o sharp.o sharp_announce.o sharp_geomhelpers.o sharp_almhelpers.o sharp_core.o sharp_legendre.o sharp_legendre_roots.o
ALLOBJ:=$(LIBOBJ) sharp_testsuite.o
LIBOBJ:=$(LIBOBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_libfftpack) $(HDR_c_utils)
$(OD)/sharp_core.o: $(SD)/sharp_core_inchelper.c $(SD)/sharp_core_inc.c $(SD)/sharp_core_inc2.c
$(OD)/sharp.o: $(SD)/sharp_mpi.c
BDEP:=$(LIB_$(PKG)) $(LIB_libfftpack) $(LIB_c_utils)

$(LIB_$(PKG)): $(LIBOBJ)

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
BIN:=$(BIN:%=$(BINDIR)/%)
$(BIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
all_cbin+=$(BIN)
