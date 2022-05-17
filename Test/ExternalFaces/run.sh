#!/bin/sh
PROGRAM=${PWD##*/}

INPUT_FILE=./Jun28b.000.n000550000.t00067.vx.json
#TEMP_FILES=output.kvsml output_value.dat

./$PROGRAM ${INPUT_FILE}
