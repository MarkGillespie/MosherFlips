# Flipping between triangulations

This repository implements Mosher's algorithm for computing a sequence of edge flips to transform an input triangulation into any other triangulation specified by normal coordinates.

The algorithm is described beginning at the bottom of page 37 of [the following paper](https://doi.org/10.2307/2000830):

Lee Mosher (Mar. 1988). “Tiling the projective foliation space of a punctured surface”.
_Transactions of the American Mathematical Society_ 306.1, pp. 1–70.
doi: https://doi.org/10.2307/2000830.




## Getting started
On mac/linux, you can set up this project with the following commands.
```bash
git clone --recursive https://github.com/MarkGillespie/MosherFlips.git
cd MosherFlips
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j7
```

Then run the code with
```
bin/run /path/to/a/mesh
```

Run the tests with
```
bin/test
```
