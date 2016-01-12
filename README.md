parADMM Engine
==============

This beta release is associated with *Testing fine-grained parallelism for the ADMM on a factor-graph* (Hao, Oghbaee, Rostami, Derbinsky, Bento).

The code includes a demo of how to use parADMM for a GPU or multiple cpu cores that is based on the circle packing example found in the numerical section of the paper above.

To compile, run: `scons`

Available flags:
* `--verbose` shows debug output (including timing data)
* `--openmp` compiles with OpenMP support
* `--cuda=/path/to/cuda` compiles with CUDA support

To run:
* `out/packing` (cpu and OpenMP)
* `out/packing-gpu` (gpu)
