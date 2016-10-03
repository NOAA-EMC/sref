#!/bin/sh

FOPTS="-g -C -traceback"

cc -c -D_UNDERSCORE fbioc.c
ftn -c $FOPTS dio.f90
ftn -c $FOPTS wrfheader.f90
