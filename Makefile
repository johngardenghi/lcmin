# Set the executables
AR := ar
FC := gfortran

# Set the compilation flags
FFLAGS := -O3 -fcheck=all

# Set LCMIN to the absolute path of the LCMIN installation files.
# The default value is the current directory of this Makefile.
LCMIN := $(CURDIR)

BIN := $(LCMIN)/bin
INC := $(LCMIN)/include
LIB := $(LCMIN)/lib
OBJ := $(LCMIN)/obj
SRC := $(LCMIN)/sources

HSL   := $(SRC)/hsl
INTER := $(SRC)/interfaces
SOLV  := $(SRC)/solver

# Fortran 90 problems codification path
F90PPATH := $(INTER)/f90

###
# Please stop your modifications here
###

export

all:
	$(MAKE) -C $(OBJ) lib

%:
	$(MAKE) -C $(OBJ) $* install-$*

clean:
	$(MAKE) -C $(OBJ) clean

distclean:
	$(MAKE) -C $(OBJ) distclean
