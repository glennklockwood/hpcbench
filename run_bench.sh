#!/bin/bash
################################################################################
#  Controller script for launching pairwise and collective benchmark 
#  applications in this benchmark suite
#
#  Glenn K. Lockwood, San Diego Supercomputer Center            December 2013
################################################################################

MPIRUN=mpirun
export OMPI_MCA_btl=self,sm,tcp

### Run all pairwise benchmarks
for benchmark in osu_bw.x osu_bibw.x osu_latency.x ring_bw.x
do
  for iteration in 2.1 2.2 2.3 2.4 2.5 2.6 2.1 2.2 2.3 2.4 2.5 2.6
  do
    $MPIRUN -machinefile ./hostfile.$iteration ./$benchmark | tee $(/bin/basename $benchmark .x).out.$iteration
  done
done

### Run collective benchmarks
for benchmark in ring_bw.x osu_alltoallv.x
do
  for iteration in $(seq 1 10)
  do
    $MPIRUN -machinefile ./hostfile.master ./$benchmark | tee $(/bin/basename $benchmark .x).out.$iteration
  done
done
