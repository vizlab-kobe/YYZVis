# YYZVis
YYZVis is a visualization toolkit based on KVS for Yin-Yang-Zhong volume dataset.

## Installation

### Prerequisite
* [KVS](https://github.com/naohisas/KVS)

### Build
Clone YYZVis repository from GitHub as follows:
```bash
$ git clone https://github.com/vizlab-kobe/YYZVis.git
```

#### Lib
Build YYZVis library required to compile each application in App and Test.
```bash
$ cd YYZVis
$ cd Lib
$ ./kvsmake.py
```

If necessary, rebuild the library as follows:
```bash
$ ./kvsmake.py rebuild
```

#### Test
Some of the test programs are in the Test directory. All of these programs can be built using the kvsmake command in each test program directory (XXX).
```bash
$ cd Test
$ cd XXX
$ kvsmake -G
$ kvsmake
$ ./run.sh
```

#### App
The App directory contains several application programs. To build and run these application programs, refer to the ReadMe files in each directories.

## Class list
The library includes the following classes used in the test programs and application programs. All of classes included in the library are defined in the namespace of `YYZVis`.

* `YYZVis::DensityMap`

* `YYZVis::Edge::CreateLineMeshObject`

* `YYZVis::Edge::CreateLineEdgeObject`

* `YYZVis::YinYangGrid`

* `YYZVis::YinYangGridSampling`

* `YYZVis::YinYangVolumeObject`

* `YYZVis::ZhongGrid`

* `YYZVis::ZhongGridSampling`

* `YYZVis::ZhongVolumeObject`

## Publications
1. Kaoru Uemori, Naohisa Sakamoto, Nobuaki Ohno, Akira Kageyama, An Interactive Visualization Toolkit for Yin-Yang-Zhong Grid Dataset, The 38th JSST Annual International Conference on Simulation Technology (JSST2019), pp.27-30, 2019.11. [[repo](https://github.com/vizlab-kobe-paper/2019_JSST__KaoruUemori)]
1. Kaoru Uemori, Naohisa Sakamoto, Nobuaki Ohno, Akira Kageyama, [YYZVis: An Efficient Visualization Toolkit for Yin-Yang-Zhong Grid Dataset](https://www.jstage.jst.go.jp/article/jasse/7/1/7_15/_article/-char/en), Journal of Advanced Simulation in Science and Engineering (JASSE), Vol.7, Issue 1, pp.15-33, 2020.3 [[repo](https://github.com/vizlab-kobe-paper/2020_JASSE__KaoruUemori))]
