################################################################################
Cray XE/6
################################################################################
o With GNU PrgEnv

module swap PrgEnv-pgi PrgEnv-gnu

./configure CC=cc --prefix=$PLATHOME/local/gasnet/1.22.4-gnu --enable-gemini \
--disable-mpi --enable-par --enable-segment-fast --disable-aligned-segments \
--disable-pshm --with-segment-mmap-max=2GB \
CPPFLAGS="-I/opt/cray/gni-headers/default/include \
-I/opt/cray/pmi/default/include" LDFLAGS="-L/opt/cray/pmi/default/lib64"

o NOTE: no --enable-cross-compile
o export CROSS_PAGESIZE=4096 required for cross compile

Legion Build:
export GASNET_ROOT=$PLATHOME/local/gasnet/1.22.4-gnu

export MPIRUN_CMD='aprun -n %N %P %A'

gasnetrun_gemini -n 1 -N 1 ./lgn-hpcg -nx 128 -ny 128 -nz 128 -s 16 -ll:cpu 16 -no-precond

Not really needed -- just use aprun

Wish For GASNet Builds of Legion. Problem: No cpu_set_t.
