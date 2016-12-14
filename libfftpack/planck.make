PKG:=libfftpack

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libfftpack.a
OBJ_$(PKG):=fftpack.o bluestein.o ls_fft.o
OBJ_$(PKG):=$(OBJ_$(PKG):%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_c_utils)

$(OD)/fftpack.o: $(SD)/fftpack_inc.c

$(OBJ_$(PKG)): $(ODEP) | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ_$(PKG))

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
