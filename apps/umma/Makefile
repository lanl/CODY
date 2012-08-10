OBJ=stack/micro-app-aos-serial.o \
    stack/micro-app-soa-serial.o \
    stack/micro-app-aos-openmp.o \
    stack/micro-app-soa-openmp.o \
    stack/micro-app-aos-cuda.o \
    stack/micro-app-soa-cuda.o \
    stack/micro-app-cuda.o

# local machine
#CFLAGS=-I/path/to/lua -I/path/to/micro-app/stack \
       -L/usr/lib/x86_64-linux-gnu -L/usr/local/cuda/lib64 \
       -fopenmp

NVCFLAGS=-I/path/to/micro-app/stack -L/home/cuda/cuda4.2/lib64 \
	 -arch=sm_20

LIBS=-llua -lm -ldl -lcudart

.SUFFIXES: .c .cu

.c.o:
	gcc $(CFLAGS) -c $< -o $@

.cu.o:
	nvcc $(NVCFLAGS) -c $< -o $@


all: micro-app-aos-serial micro-app-aos-openmp \
    micro-app-soa-serial micro-app-soa-openmp \
    micro-app-aos-cuda micro-app-soa-cuda

micro-app-aos-serial: stack/micro-app-aos-serial.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-aos-openmp: stack/micro-app-aos-openmp.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-aos-cuda: stack/micro-app-aos-cuda.o stack/micro-app-cuda.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-serial: stack/micro-app-soa-serial.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-openmp: stack/micro-app-soa-openmp.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

micro-app-soa-cuda: stack/micro-app-soa-cuda.o stack/micro-app-cuda.o
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY:  clean

clean:
	rm -f micro-app-aos-serial micro-app-aos-openmp \
	    micro-app-soa-serial micro-app-soa-openmp \
	    micro-app-aos-cuda micro-app-soa-cuda stack/*.o
