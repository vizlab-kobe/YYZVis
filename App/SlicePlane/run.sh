#!/bin/sh
DIM_RAD=201
DIM_LAT=204
DIM_LON=608
DIM_ZHONG=222
DATA_DIR=~/Work/Data/MHD

YIN_DATA=${DATA_DIR}/Jun28b.000.wyin.vx.n000550000.t00067
YANG_DATA=${DATA_DIR}/Jun28b.000.wyng.vx.n000550000.t00067
ZHONG_DATA=${DATA_DIR}/Jun28b.000.icore_3d.vx.n000550000.t00067

./SlicePlane -dim_rad ${DIM_RAD} -dim_lat ${DIM_LAT} -dim_lon ${DIM_LON} -dim_zhong ${DIM_ZHONG} -yin ${YIN_DATA} -yang ${YANG_DATA} -zhong ${ZHONG_DATA}
