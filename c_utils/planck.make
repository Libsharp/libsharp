PKG:=c_utils

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libc_utils.a
OBJ_$(PKG):=trig_utils.o c_utils.o walltime_c.o memusage.o
OBJ_$(PKG):=$(OBJ_$(PKG):%=$(OD)/%)

$(OBJ_$(PKG)): $(HDR_$(PKG)) | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ_$(PKG))

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
