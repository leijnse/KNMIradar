#! /bin/bash
#
# Makefile
#

#Setting names of program

SOURCE="reduce_ODIM_h5_volume.c"
BINARY="reduce_ODIM_h5_volume.x"

#Setting compiler and options

COMPILE="gcc -Wall -O"

#Setting include and library directories

IPATH=""
LPATH=""

#Setting libraries

LIBRARIES="-lhdf5_hl -lhdf5 -lz -lm"

#Compiling of source

$COMPILE $SOURCE -o $BINARY $IPATH $LPATH $LIBRARIES
