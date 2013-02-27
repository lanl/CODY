CC=nvcc
NVCC=nvcc
EXEC=hydro
HEADERS=hydro.h hydro_struct.h
OBJS=main.o dev_funcs.o hydro.o outfile.o
LIBS=-L/opt-/mpi/openmpi/1.5.1/pgi/11.1/lib/ -lmpi -lm
INCS=-I"/opt-/mpi/openmpi/1.5.1/pgi/11.1/include/"
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
