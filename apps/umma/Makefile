#OBJ=stack/micro-app-aos-serial.o \
    stack/micro-app-soa-serial.o \
    stack/micro-app-aos-openmp.o \
    stack/micro-app-soa-openmp.o \
    stack/cuda/micro-app-aos-cuda.o \
    stack/cuda/micro-app-soa-cuda.o \
    stack/cuda/micro-app-cuda.o \
    stack/ispc/micro-app-soa-ispc.o \
    stack/ispc/micro-app-soa.o \
    stack/ispc/micro-app-aos-ispc.o \
    stack/ispc/micro-app-aos.o

#--- local machine
CFLAGS=-I/path/to/lua -I/path/to/micro-app/stack \
       -L/usr/lib/x86_64-linux-gnu -L/usr/local/cuda/lib64 \
       -fopenmp

ISPC=/path/to/ispc

ISPC_FLAGS=--wno-perf

LIBS=-llua -lm -ldl -lcudart

NVCFLAGS=-I/path/to/micro-app/stack -L/home/cuda/cuda4.2/lib64 \
	 -arch=sm_20

.SUFFIXES: .c .cu .ispc

.c.o:
	gcc $(CFLAGS) -c $< -o $@

.cu.o:
	nvcc $(NVCFLAGS) -c $< -o $@

.ispc.o:
	$(ISPC) $(ISPC_FLAGS) $< -o $@


all: micro-app-aos-serial micro-app-aos-openmp \
    micro-app-soa-serial micro-app-soa-openmp \
    micro-app-aos-cuda micro-app-soa-cuda \
    micro-app-soa-ispc micro-app-aos-ispc

micro-app-aos-serial: stack/micro-app-aos-serial.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-aos-openmp: stack/micro-app-aos-openmp.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-aos-cuda: stack/cuda/micro-app-aos-cuda.o stack/cuda/micro-app-cuda.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-aos-ispc: stack/ispc/micro-app-aos-ispc.o \
    stack/ispc/micro-app-aos.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-serial: stack/micro-app-soa-serial.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-openmp: stack/micro-app-soa-openmp.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-cuda: stack/cuda/micro-app-soa-cuda.o stack/cuda/micro-app-cuda.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-ispc: stack/ispc/micro-app-soa-ispc.o \
    stack/ispc/micro-app-soa.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

headers:
	$(ISPC) $(ISPC_FLAGS) stack/ispc/micro-app-aos.ispc \
	    -h stack/ispc/micro-app-aos.h
	$(ISPC) $(ISPC_FLAGS) stack/ispc/micro-app-soa.ispc \
	    -h stack/ispc/micro-app-soa.h

.PHONY:  clean

clean:
	rm -f micro-app-aos-serial micro-app-aos-openmp \
	    micro-app-soa-serial micro-app-soa-openmp \
	    micro-app-aos-cuda micro-app-soa-cuda \
	    micro-app-soa-ispc micro-app-aos-ispc \
	    stack/*.o stack/cuda/*.o stack/ispc/*.o
	
