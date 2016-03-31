#!/bin/sh

FOPTS="-g -C -traceback"

icc -c -D_UNDERSCORE fbioc.c
ifort -c $FOPTS dio.f90
ifort -c $FOPTS wrfheader.f90
