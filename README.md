HPC Benchmarks
==============

This repository contains a basic suite of benchmarks that I use for for testing 
HPC performance across different system architectures.  This is an ongoing 
project that is still far from complete or self-contained; notably, 
documentation is severely lacking.  

At present, this repository is comprised of a few utility functions for manually
working on a cluster of nodes without a resource manager, a few simple 
synthetic benchmarks, and a few assorted files to make life easier.

Framework for this Benchmark Suite
----------------------------------
* Makefile - build the microbenchmarks and their nodefiles
* NODES - the nodelist containing all unique nodes to be used
* gen_nodefiles.sh - generate all permutations of node pairs from the NODES file.  Needed by Makefile
* run_bench.sh - run all microbenchmarks

Utility Scripts
---------------
* doall - Execute a command on all of the nodes in the NODES file
* pushout - Distribute a file to all of the nodes in the NODES file

High-Performance Linpack (HPL)
------------------------------
* Make.gordon - Build HPL using Intel + MKL (tested on Gordon, a Sandy Bridge-based cluster)
* HPL.dat - Input for HPL to get best performance on Gordon (Sandy Bridge).  Significantly different from default HPL.dat that ships with HPL.
* run_hpl.qsub - Torque submit script (or standalone bash script) to run HPL

Microbenchmarks, Point-to-Point
-------------------------------
* osu_bw.c - OSU bandwidth microbenchmark
* osu_latency.c - OSU latency microbenchmark
* osu_bibw.c - OSU bidirectional bandwidth benchmark

Microbenchmarks, Collective
---------------------------
* osu_alltoallv.c - OSU alltoallv microbenchmark
* osu_coll.c - needed by osu_alltoallv
* osu_coll.h - needed by osu_alltoallv
* ring_bw.c - Natural- and random-ring bandwidth and latency benchmark.  I don't really like this benchmark and am in the process of updating it to be a little more robust and meaningful.
