CC=nvcc
NVCC=nvcc
EXEC=hydro
HEADERS=hydro.h hydro_struct.h
OBJS=main.o dev_funcs.o hydro.o outfile.o
LIBS=-L"/projects/opt/mpi/openmpi/1.6.4-xrc_centos_gcc-4.4.7/lib/" -lmpi -lm
INCS=-I"/projects/opt/mpi/openmpi/1.6.4-xrc_centos_gcc-4.4.7/include/"
CUFLAGS=-arch=sm_13


all: ${EXEC}

${EXEC}:${OBJS} ${HEADERS}
	${CC} ${CFLAGS} ${LIBS} ${OBJS} -o ${EXEC}

debug:CFLAGS+=-g
debug: all

optim:CFLAGS+=-O3
optim: all

%.o: %.cu
	${NVCC} ${INCS} ${CFLAGS} ${CUFLAGS} -c $<

.PHONY: clean
clean:
	rm *.o ${EXEC}
