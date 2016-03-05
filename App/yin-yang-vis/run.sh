#!/bin/sh
REPEAT_LEVEL=1
YIN_DATA_FILE=../../../bx_vx/oct09b.011.wyin.vx.n000250000.t00302
YANG_DATA_FILE=../../../bx_vx/oct09b.011.wyng.vx.n000250000.t00302
ZONG_DATA_FILE=../../../bx_vx/oct09b.011.icore_3d.vx.n000250000.t00302
./yin-yang-vis $REPEAT_LEVEL $YIN_DATA_FILE $YANG_DATA_FILE $ZONG_DATA_FILE
