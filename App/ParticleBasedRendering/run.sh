#!/bin/sh
DIM_RAD=201
DIM_LAT=204
DIM_LON=608
DIM_ZHONG=222
REPEAT=1
DATA_DIR=~/Work/Test/tsubota/bx_vx

YIN_DATA=${DATA_DIR}/oct09b.011.wyin.bx.n000250000.t00302
YANG_DATA=${DATA_DIR}/oct09b.011.wyng.bx.n000250000.t00302
ZHONG_DATA=${DATA_DIR}/oct09b.011.icore_3d.bx.n000250000.t00302

./ParticleBasedRendering -dim_rad ${DIM_RAD} -dim_lat ${DIM_LAT} -dim_lon ${DIM_LON} -dim_zhong ${DIM_ZHONG} -repeat ${REPEAT} -yin ${YIN_DATA} -yang ${YANG_DATA} -zhong ${ZHONG_DATA}
