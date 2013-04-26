CC=pgcc -acc
EXEC=hydro
HEADERS=hydro.h hydro_struct.h
OBJS=main.o hydro.o outfile.o
CFLAGS+=
LIBS=-lm

all: ${EXEC}

${EXEC}:${OBJS} ${HEADERS}
	${CC} ${CFLAGS} ${LIBS} ${OBJS} -o ${EXEC}

debug:CFLAGS+=-g
debug: all

optim:CFLAGS+=-O3
optim: all

.PHONY: clean
clean:
	rm *.o ${EXEC}
