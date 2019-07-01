# YinYangVis
YinYangVis ia a set of visualization modues based on KVS for Yin-Yang-Zhong volume dataset.

## Requirements
* [KVS](https://github.com/naohisas/KVS)

## Building from sources
Clone the repository.
```bash
$ git clone https://github.com/vizlab-kobe/YinYangVis.git
```

Build YinYangVis library required to compile each application in App.
```bash
$ cd Lib
$ ./kvsmake.py
```

Build applications with kvsmake in each application directory.
```bash
e.g.) App/Isosurface
$ cd App/Isosurface
$ kvsmake -G
$ kvsmake
$ ./run.sh
```
