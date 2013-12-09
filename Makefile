MPICC=mpicc

all: osu_bw.x osu_bibw.x osu_latency.x osu_alltoallv.x ring_bw.x hostfile.master

hostfile.master: $(NODEFILE)
	./gen_hostfiles.sh

%.x: %.c
	$(MPICC) -o $@ $*.c -lm
