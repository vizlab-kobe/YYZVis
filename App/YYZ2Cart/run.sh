#!/bin/sh
PROGRAM=${PWD##*/}

DIM=200
INPUT_FILE=./Jun28b.000.n000550000.t00067.vx.json
OUTPUT_FILE=output.kvsml
TEMP_FILES=output*

./$PROGRAM -dim ${DIM} -output ${OUTPUT_FILE} ${INPUT_FILE}
