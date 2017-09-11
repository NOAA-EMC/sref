#! /usr/bin/env bash
module load nco-gnu-sandybridge/4.4.4

set -eux

stretch_factor=1.5
target_lon=-97.5
target_lat=35.5

cube_res=768
parent_tile=6
refine_ratio=3
istart_nest=193
jstart_nest=329
iend_nest=1344
jend_nest=1288

for file in "$@"
do

echo "${file}"
ncatted -O -a stretch_factor,global,c,d,${stretch_factor} -h "${file}"
ncatted -O -a target_lat,global,c,d,${target_lat} -h "${file}"
ncatted -O -a target_lon,global,c,d,${target_lon} -h "${file}"

ncatted -O -a cube_res,global,c,i,${cube_res} -h "${file}"
ncatted -O -a parent_tile,global,c,i,${parent_tile} -h "${file}"
ncatted -O -a refine_ratio,global,c,i,${refine_ratio} -h "${file}"
ncatted -O -a istart_nest,global,c,i,${istart_nest} -h "${file}"
ncatted -O -a jstart_nest,global,c,i,${jstart_nest} -h "${file}"
ncatted -O -a iend_nest,global,c,i,${iend_nest} -h "${file}"
ncatted -O -a jend_nest,global,c,i,${jend_nest} -h "${file}"

done
