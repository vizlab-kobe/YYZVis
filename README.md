# YinYangVis
Yin-Yang Grid Data Visualization.
This program requires that KVS has been installed.


###compiling Library
```
cd YinYangVis/Lib
./kvsmake.py
```
###compiling Application
```
cd YinYangVis/App/yin-yang-vis
kvsmake -G
kvsmake 
```
###run Application
```
sh run.sh
```
+ Put data in the YinYangVis directory with the same height.
+ In run.sh YIN_DATA_FILE,YANG_DATA_FILE,ZHONG_DATA_FILE,you write the data name.
+ You can decide the repeat count by changing run.sh NUMBER_OF_REPETITIONS.
