#!/bin/bash
#
#  Generate nodefiles of all node pairs contained in the master NODES file
#

uniq NODES > hostfile.master
nodect=$(cat hostfile.master | wc -l)

seqnum=1
for i in $(seq 1 $((nodect-1)))
do
    n1=$(sed -e "${i}q;d" hostfile.master)
    for j in $(seq $((i+1)) $nodect)
    do
       n2=$(sed -e "${j}q;d" hostfile.master)
       echo $n1 > hostfile.2.$seqnum
       echo $n2 >> hostfile.2.$seqnum
       let "seqnum++"
    done
done
