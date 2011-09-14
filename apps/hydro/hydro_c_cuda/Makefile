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

HLW_SRCC= $(H_SRCC)

HLW_SRCCU= $(H_SRCCU) \
    engine_laxWend.cu \
    device_funcs_lw.cu

HU_SRCC=

HU_SRCCU= \
  utils.cu \
  hydro_dmp.cu

#Implied files
HEADERS=$(wildcard *.h)
HG_OBJS=$(patsubst %.c,%.o,$(HG_SRCC)) $(patsubst %.cu,%.o,$(HG_SRCCU))
HG_DEPS=$(patsubst %.c,%.d,$(HG_SRCC)) $(patsubst %.cu,%.d,$(HG_SRCCU))
HLW_OBJS=$(patsubst %.c,%.o,$(HLW_SRCC)) $(patsubst %.cu,%.o,$(HLW_SRCCU))
HLW_DEPS=$(patsubst %.c,%.d,$(HLW_SRCC)) $(patsubst %.cu,%.d,$(HLW_SRCCU))
HU_OBJS=$(patsubst %.c,%.o,$(HU_SRCC)) $(patsubst %.cu,%.o,$(HU_SRCCU))
HU_DEPS=$(patsubst %.c,%.d,$(HU_SRCC)) $(patsubst %.cu,%.d,$(HU_SRCCU))

#Compiler configurations

CC=gcc
NVCC=nvcc
CUFLAGS=-arch=sm_13 --maxrregcount=32
#CUFLAGS=-Xptxas="-v" -arch=sm_20

#Executable targets
.PHONY: all
all:	$(EXECS)

$(PGM)_god: $(HG_DEPS) $(HG_OBJS)
	echo $(HG_DEPS)
	$(NVCC) $(CFLAGS) ${CUFLAGS} -o $(PGM)_god $(OPT) $(HG_OBJS) -lm

$(PGM)_lw: $(HLW_DEPS) $(HLW_OBJS)
	$(NVCC) $(CFLAGS) ${CUFLAGS} -o $(PGM)_lw $(OPT) $(HLW_OBJS) -lm

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

.PHONY: run_god run_lw
run: run_god run_lw
RUN_OPTS=-O3
run_lw: CFLAGS+=$(RUN_OPTS)
run_lw: $(PGM)_lw
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p \
		../$(PGM)_lw $(TEST_OPTS) &> test_lw.log

run_god: CFLAGS+=$(RUN_OPTS)
run_god: $(PGM)_god
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p \
		../$(PGM)_god $(TEST_OPTS) &> test_god.log

.PHONY: vrun_god vrun_lw
vrun: vrun_god vrun_lw
vrun_lw: $(PGM)_lw
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p valgrind --tool=memcheck \
		--leak-check=full -v \
		../$(PGM)_lw $(TEST_OPTS) &> test_vrun_lw.log

vrun_god: $(PGM)_god
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p valgrind --tool=memcheck \
		--leak-check=full -v \
		../$(PGM)_god $(TEST_OPTS) &> test_vrun_god.log

.PHONY: crun_god crun_lw
crun: crun_god crun_lw
CRUN_FLAGS=-g -G
crun_lw: CFLAGS+=$(CRUN_FLAGS)
crun_lw: clean $(PGM)_lw
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p cuda-memcheck \
		../$(PGM)_lw $(TEST_OPTS) &> test_crun_lw.log

crun_god: CFLAGS+=$(CRUN_FLAGS)
crun_god: clean $(PGM)_god
	mkdir -p $(OUT_DIR)
	cd $(OUT_DIR) && /usr/bin/time -p cuda-memcheck \
		../$(PGM)_god $(TEST_OPTS) &> test_crun_god.log

.PHONY: debug_god debug_lw
debug:	CFLAGS+=-g -G
debug:	clean debug_god debug_lw
debug_lw: CFLAGS+=-g -G
debug_lw: clean $(PGM)_lw

debug_god: CFLAGS+=-g -G
debug_god: clean $(PGM)_god

tests:	crun

