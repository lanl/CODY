CC=nvcc
NVCC=nvcc
EXEC=hydro
HEADERS=hydro.h hydro_struct.h
OBJS=main.o dev_funcs.o hydro.o outfile.o
LIBS=-lm
CUFLAGS=-arch=sm_20


all: ${EXEC}

${EXEC}:${OBJS} ${HEADERS}
	${CC} ${CFLAGS} ${LIBS} ${OBJS} -o ${EXEC}

debug:CFLAGS+=-G -g
debug: all

optim:CFLAGS+=-O3
optim: all

%.o: %.cu
	${NVCC} ${CFLAGS} ${CUFLAGS} -c $<

.PHONY: clean
clean:
	rm *.o ${EXEC}
