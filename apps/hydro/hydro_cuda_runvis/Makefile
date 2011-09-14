PGM=hydro
UTILS= hydro_utils
EXECS=$(PGM)_god $(PGM)_lw $(UTILS)

SRC_DIR=$(PWD)
IN_DIR=$(PWD)/input
OUT_DIR=$(PWD)/output

#Source files

H_SRCC=

H_SRCCU= \
    main.cu \
    hydro_args.cu \
    hydro_utils.cu \
    hydro_dmp.cu

HG_SRCC= $(H_SRCC)

HG_SRCCU= $(H_SRCCU) \
    engine_godunov.cu \
    device_funcs_god.cu

HU_SRCC=

HU_SRCCU= \
  utils.cu \
  hydro_dmp.cu

#Implied files
HEADERS=$(wildcard *.h)
HG_OBJS=$(patsubst %.c,%.o,$(HG_SRCC)) $(patsubst %.cu,%.o,$(HG_SRCCU))
HG_DEPS=$(patsubst %.c,%.d,$(HG_SRCC)) $(patsubst %.cu,%.d,$(HG_SRCCU))
HU_OBJS=$(patsubst %.c,%.o,$(HU_SRCC)) $(patsubst %.cu,%.o,$(HU_SRCCU))
HU_DEPS=$(patsubst %.c,%.d,$(HU_SRCC)) $(patsubst %.cu,%.d,$(HU_SRCCU))

#Compiler configurations
LIBS = -lglut -lGLU -lGL


CC=gcc
NVCC=nvcc
CUFLAGS=-arch=sm_13 --maxrregcount=32 $(LIBS)
#CUFLAGS=-Xptxas="-v" -arch=sm_20

#Executable targets
.PHONY: all
all:	$(EXECS)

$(PGM)_god: $(HG_DEPS) $(HG_OBJS)
	$(NVCC) $(CFLAGS) ${CUFLAGS} -o $(PGM)_god $(OPT) $(HG_OBJS) -lm

$(UTILS): $(HU_DEPS) $(HU_OBJS)
	$(NVCC) $(CFLAGS) ${CUFLAGS} -o $(UTILS) $(OPT) $(HU_OBJS) -lm

#Dependencies
-include $(DEP)

#Implicit targets
.SUFFIXES= .c .o .d .cu
%.o:	%.c
	$(CC) ${CFLAGS} -c $<

%.o:	%.cu
	$(NVCC) ${CFLAGS} ${CUFLAGS} -c $<

%.d:	%.c
	@$(CC) ${CPPFLAGS} $(CFLAGS) $(HFLAGS) -M $< -o $@
	@sed -i.bak "s|\($*\)\.o[ :]*|\1.o $@ : |g" $@
	@rm $@.bak

%.d:	%.cu
	@$(NVCC) ${CPPFLAGS} $(CFLAGS) ${CUFLAGS} $(HFLAGS) -M $< -o $@
	@sed -i.bak "s|\($*\)\.o[ :]*|\1.o $@ : |g" $@
	@rm $@.bak

#Clean targets
.PHONY: clean cleanall

clean:
	-rm *.o

cleanall:
	-rm *.o *.d *.bak $(EXECS)

#Run targets
.PHONY: run crun vrun tests debug

TEST_OPTS=-i $(IN_DIR)/test.nml

.PHONY: run_god
run: run_god
RUN_OPTS=-O3

run_god: CFLAGS+=$(RUN_OPTS)
run_god: $(PGM)_god
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p \
		../$(PGM)_god $(TEST_OPTS) &> test_god.log

.PHONY: vrun_god
vrun: vrun_god

vrun_god: $(PGM)_god
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p valgrind --tool=memcheck \
		--leak-check=full -v \
		../$(PGM)_god $(TEST_OPTS) &> test_vrun_god.log

.PHONY: crun_god
crun: crun_god
CRUN_FLAGS=-g -G

crun_god: CFLAGS+=$(CRUN_FLAGS)
crun_god: clean $(PGM)_god
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p cuda-memcheck \
		../$(PGM)_god $(TEST_OPTS) &> test_crun_god.log

.PHONY: debug_god
debug:	CFLAGS+=-g -G
debug:	clean debug_god

debug_god: CFLAGS+=-g -G
debug_god: clean $(PGM)_god

tests:	crun

