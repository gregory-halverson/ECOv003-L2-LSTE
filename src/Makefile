# !gmake
# DEBUG = 1
#
# !Makefile Name: Makefile
#
# !Description: Makefile for L2_PGE Land Surface Temperature and Emissivity.
#               This makefile used to generate l2_tes_pge
#
# !Env Variables:
#
#      ENV VARIABLES DESCRIPTION
#      ~~~~~~~~~~~~~ ~~~~~~~~~~~
#      CC            The c compiler
#      CFLAGS        Compiler flags
#      HDFINC        Include directory of HDF
#      HDF5INC       Include directory of HDF5
#      HDFLIB        Library for HDF
#      HDF5LIB       Library for HDF5
#
# !Revision History:
#  Date:        PR#     Author          Description
#  -------      -----   ------          -----------
#  05/03/16             R. Freepartner  Initial Version
###################################################################
# If using conda environment, add its lib/pkgconfig to PKG_CONFIG_PATH
ifdef CONDA_PREFIX
PKG_CONFIG_PATH := $(CONDA_PREFIX)/lib/pkgconfig:$(PKG_CONFIG_PATH)
export PKG_CONFIG_PATH
endif

# Auto-discover dependencies using pkg-config
# This makes the build portable across any system with installed packages

PKG_CHECK = $(shell pkg-config --cflags $(1) 2>/dev/null)
PKG_LIBS = $(shell pkg-config --libs $(1) 2>/dev/null)

# HDF4 (required) - provides hdf.h, mfhdf.h, hdfi.h
# Try pkg-config first, then check conda/homebrew paths (HDF4 often doesn't have .pc file)
HDF4_PREFIX := $(shell if [ -n "$$CONDA_PREFIX" ]; then echo $$CONDA_PREFIX; fi)
ifeq ($(HDF4_PREFIX),)
  HDF4_PREFIX := $(shell command -v brew >/dev/null 2>&1 && brew --prefix hdf4 2>/dev/null || echo "")
endif

ifeq ($(HDF4_PREFIX),)
  HDFINC := $(call PKG_CHECK,hdf)
  HDFLIB := $(call PKG_LIBS,hdf)
else
  HDFINC := -I$(HDF4_PREFIX)/include
  HDFLIB := -L$(HDF4_PREFIX)/lib -lmfhdf -ldf
endif

ifeq ($(HDFINC),)
$(error HDF4 not found. Use: conda create -n ecostress -c conda-forge hdf4 && conda activate ecostress)
endif

# HDF5 (required) - provides hdf5.h
HDF5INC := $(call PKG_CHECK,hdf5)
HDF5LIB := $(call PKG_LIBS,hdf5)
ifeq ($(HDF5INC),)
$(error HDF5 not found. Install with: brew install hdf5 (macOS) or apt-get install libhdf5-dev (Ubuntu))
endif

# LibXML2 (required) - provides libxml/xmlwriter.h
# LibXML2 often doesn't provide .pc file in conda, so check conda/homebrew paths
XMLINC := $(call PKG_CHECK,libxml-2.0)
ifeq ($(XMLINC),)
  # Try conda/homebrew paths
  ifdef CONDA_PREFIX
    XMLINC := -I$(CONDA_PREFIX)/include/libxml2
    XMLLIB := -L$(CONDA_PREFIX)/lib -lxml2
  else
    XMLINC := $(shell if command -v brew >/dev/null 2>&1; then echo -I$$(brew --prefix libxml2)/include/libxml2; fi)
    XMLLIB := $(shell if command -v brew >/dev/null 2>&1; then echo -L$$(brew --prefix libxml2)/lib -lxml2; fi)
  endif
endif
ifeq ($(XMLINC),)
$(error libxml2 not found. Use: conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config && conda activate ecostress)
endif
XMLLIB_VAL := $(XMLLIB)
ifeq ($(XMLLIB_VAL),)
  XMLLIB_VAL := $(call PKG_LIBS,libxml2-2.0)
endif
ifeq ($(XMLLIB_VAL),)
  XMLLIB_VAL := -lxml2
endif
XMLLIB := $(XMLLIB_VAL)

# Expat (required) - provides expat.h
# Expat often doesn't provide .pc file in conda, so check conda/homebrew paths
EXPATLIB := $(call PKG_LIBS,expat)
ifeq ($(EXPATLIB),)
  # Try conda/homebrew paths
  ifdef CONDA_PREFIX
    EXPATLIB := -L$(CONDA_PREFIX)/lib -lexpat
  else
    EXPATLIB := $(shell if command -v brew >/dev/null 2>&1; then echo -L$$(brew --prefix expat)/lib -lexpat; fi)
  endif
endif
ifeq ($(EXPATLIB),)
$(error expat not found. Use: conda create -n ecostress -c conda-forge hdf4 hdf5 libxml2 eccodes pkg-config && conda activate ecostress)
endif

# GRIB API / ecCodes (required) - provides grib_api.h
# Try ecCodes first (modern), fall back to grib_api (legacy)
GRIBINC := $(call PKG_CHECK,eccodes)
GRIBLIB := $(call PKG_LIBS,eccodes)
ifeq ($(GRIBLIB),)
  GRIBINC := $(call PKG_CHECK,grib_api)
  GRIBLIB := $(call PKG_LIBS,grib_api)
endif
ifeq ($(GRIBLIB),)
$(error ecCodes/GRIB API not found. Install with: brew install eccodes (macOS) or apt-get install libeccodes-dev (Ubuntu))
endif

##################
# EXPAT settings
# Adds support for checking lexical aspects of namespaces
XP_NS = -DXML_NS
# Adds support for processing DTDs
XP_DTD = -DXML_DTD
# Byte order macro
#This is for little endian machines
XP_BO = -DXML_BYTE_ORDER=12

XP_DEFINES = $(XP_NS) $(XP_DTD) $(XP_BO) $(XP_UNI) $(XP_UNI_WC) $(XP_MM)
##################

ifdef DEBUG
ADD_CFLAGS = -g -Wall
else
ADD_CFLAGS = -O3 -Wall
endif

# Object files
TES_OBJ = \
    tes_main.o \
    asterged.o \
    cloud.o \
    config.o \
    describe.o \
    error.o \
    fileio.o \
    interps.o \
    lste_lib.o \
    matrix.o \
    maptogrid.o \
    metadata.o \
    rttov_util.o \
    smooth2d.o \
    tes_util.o \
    tg_wvs.o 

ALLOBJ = $(TES_OBJ)

# Combine the includes using pkg-config results
INC = -I. $(HDFINC) $(HDF5INC) $(XMLINC) $(GRIBINC) -I/usr/include

# Header files
HDR = asterged.h \
      cloud.h \
      config.h \
      describe.h \
      fileio.h \
      interps.h \
      lste_lib.h \
      maptogrid.h \
      matrix.h \
      metadata.h \
      rttov_util.h \
      smooth2d.h \
      tes_util.h \
      tg_wvs.h
      
# Combine the libraries using pkg-config results
LIB = $(HDFLIB) $(HDF5LIB) $(GRIBLIB) $(XMLLIB) $(EXPATLIB) -ldl -ldf -ljpeg -lz -lm 

# Define executable name
TARGET = L2_PGE	

all:	$(TARGET)
$(TARGET):	$(TES_OBJ) $(HDR)
	 $(CC) $(ADD_CFLAGS) $(XP_DEFINES) $(TES_OBJ) $(LIB) -o $(TARGET)

.c.o:
	$(CC) $(CFLAGS) $(ADD_CFLAGS) $(INC) -c $< -o $@
	
# Delete object files
clean: clean_exe clean_o

clean_o:
	rm -f $(ALLOBJ)

# Delete the executable file $TARGET
clean_exe:
	rm -f $(TARGET)

###################################################################
