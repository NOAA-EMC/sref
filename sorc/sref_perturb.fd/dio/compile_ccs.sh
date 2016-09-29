#!/bin/sh

FOPTS="-g -qflttrap=zerodivide:invalid:enable -qsigtrap -C -qinitauto=FF911299"

xlc -c fbioc.c 
xlf90 -c $FOPTS dio.f90
xlf90 -c $FOPTS wrfheader.f90
