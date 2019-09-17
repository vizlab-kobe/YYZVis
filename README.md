# YinYangVis
YinYangVis is a visualization toolkit based on KVS for Yin-Yang-Zhong volume dataset.

## Requirements
* [KVS](https://github.com/naohisas/KVS)

## Compiling
Clone the repository.
```bash
$ git clone https://github.com/vizlab-kobe/YinYangVis.git
```

Build YinYangVis library required to compile each application in App.
```bash
$ cd Lib
$ ./kvsmake.py
```

Build applications with kvsmake in each application directory. e.g.) App/Isosurface
```bash
$ cd App/Isosurface
$ kvsmake -G
$ kvsmake
$ ./run.sh
```

## Publications

1. Kaoru Uemori, Naohisa Sakamoto, Nobuaki Ohno, Akira Kageyama, An Interactive Visualization Toolkit for Yin-Yang-Zhong Grid Dataset, The 38th JSST Annual International Conference on Simulation Technology (JSST2019), (Accepted), 2019. [[repo](https://github.com/vizlab-kobe-paper/2019_JSST__KaoruUemori)]
