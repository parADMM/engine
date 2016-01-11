parADMM Engine
==============

This beta release is associated with *Testing fine-grained parallelism for the ADMM on a factor-graph* (Hao, Oghbaee, Rostami, Derbinsky, Bento).

As of right now it only includes circle packing demo programs comparing single CPU core vs. multi-CPU vs. GPU.

To compile, run: `scons`

Available flags:
* `--verbose` shows debug output (including timing data)
* `--openmp` compiles with OpenMP support
* `--cuda=/path/to/cuda` compiles with CUDA support

To run:
* `out/packing` (cpu and OpenMP)
* `out/packing-gpu` (gpu)
