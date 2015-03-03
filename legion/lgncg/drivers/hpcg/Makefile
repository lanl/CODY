#
# Copyright 2013 Stanford University
# Copyright (c) 2014      Los Alamos National Security, LLC
#                         All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# LA-CC 10-123
#

ifndef LG_RT_DIR
$(error LG_RT_DIR variable is not defined, aborting build)
endif

# Flag indicating whether or not run-time offset density checks will be
# performed. See: offsetsAreDense in src/utils.h. Calls to offsetsAreDense are
# not free, but are useful when writing code and moving to new architectures.
LGNCG_DO_DENSITY_CHECKS = 0

ifeq ($(LGNCG_DO_DENSITY_CHECKS), 1)
LGNCG_CC_FLAGS = -DLGNCG_CHECK_OFF_R_DENSE
endif

# Flag indicating whether or not to enable Legion tracing within CG.
LGNCG_DO_CG_TRACE = 0
ifeq ($(LGNCG_DO_CG_TRACE), 1)
LGNCG_CC_FLAGS = -DLGNCG_DO_CG_TRACE
endif

# for now just point to top-level src dir of lgncg
LGCG_DIR := ../../src

#Flags for directing the runtime makefile what to include
DEBUG=0                   # Do NOT Include debugging symbols
#DEBUG=1                   # Include debugging symbols
#OUTPUT_LEVEL=LEVEL_DEBUG  # Compile time print level
OUTPUT_LEVEL=LEVEL_NONE # Compile time print level
SHARED_LOWLEVEL=1	  # Use the shared low level
#ALT_MAPPERS=1		  # Compile the alternative mappers

# Put the binary file name here
OUTFILE	:= lgn-hpcg
# List all the application source files here
GEN_SRC	:= \
$(wildcard *.cc) \
$(wildcard $(LGCG_DIR)/*.cc)
GEN_GPU_SRC	:=
# make sure that we depend on the headers of lgncg and our driver
includes = \
$(wildcard *.h) \
$(wildcard $(LGCG_DIR)/*.h)

# You can modify these variables, some will be appended to by the runtime makefile
INC_FLAGS    := -I$(LGCG_DIR)
CC_FLAGS := $(LGNCG_CC_FLAGS)
# for gcc on OS X
#CC_FLAGS     := -Wall -O2 -Wa,-q
# optimized
#CC_FLAGS     += -Wall -Wextra -Ofast -ftree-vectorize -ffast-math -Wa,-q
# optimized/w profiling make sure OUTPUT_LEVEL=LEVEL_DEBUG
CC_FLAGS     += -Wall -Ofast -ffast-math -DLEGION_PROF -DMAX_FIELDS=64 -DFIELD_LOG2=6
# debug
#CC_FLAGS += -Wall -Wextra -O0 -DPRIVILEGE_CHECKS -DBOUNDS_CHECKS -DINORDER_EXECUTION -DFULL_SIZE_INSTANCES -Wa,-q
#CC_FLAGS += -Wall -Wextra -O0 -DPRIVILEGE_CHECKS -DBOUNDS_CHECKS -Wa,-q
#CC_FLAGS     += -Wall -Ofast -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=0 -Wa,-q
# debug
#CC_FLAGS     += -Wall -O0 -g
NVCC_FLAGS   :=
GASNET_FLAGS :=

###########################################################################
#
#   Don't change anything below here
#
###########################################################################

# All these variables will be filled in by the runtime makefile
LOW_RUNTIME_SRC	:=
HIGH_RUNTIME_SRC:=
GPU_RUNTIME_SRC	:=
MAPPER_SRC	:=

include $(LG_RT_DIR)/runtime.mk

# General shell commands
SHELL	:= /bin/sh
SH	:= sh
RM	:= rm -f
LS	:= ls
MKDIR	:= mkdir
MV	:= mv
CP	:= cp
SED	:= sed
ECHO	:= echo
TOUCH	:= touch
MAKE	:= make
ifndef GCC
GCC	:= g++
# for gcc on OS X
#GCC	:= g++-4.8 -mavx -march=corei7-avx
#GCC	:= g++-4.8
endif
ifndef NVCC
NVCC	:= $(CUDA)/bin/nvcc
endif
SSH	:= ssh
SCP	:= scp

common_all : all

.PHONY	: common_all

GEN_OBJS	:= $(GEN_SRC:.cc=.o)
LOW_RUNTIME_OBJS:= $(LOW_RUNTIME_SRC:.cc=.o)
HIGH_RUNTIME_OBJS:=$(HIGH_RUNTIME_SRC:.cc=.o)
MAPPER_OBJS	:= $(MAPPER_SRC:.cc=.o)
# Only compile the gpu objects if we need to
ifndef SHARED_LOWLEVEL
GEN_GPU_OBJS	:= $(GEN_GPU_SRC:.cu=.o)
GPU_RUNTIME_OBJS:= $(GPU_RUNTIME_SRC:.cu=.o)
else
GEN_GPU_OBJS	:=
GPU_RUNTIME_OBJS:=
endif

ALL_OBJS	:= $(GEN_OBJS) $(GEN_GPU_OBJS) $(LOW_RUNTIME_OBJS) $(HIGH_RUNTIME_OBJS) $(GPU_RUNTIME_OBJS) $(MAPPER_OBJS)

all:
	$(MAKE) $(OUTFILE)

# If we're using the general low-level runtime we have to link with nvcc
$(OUTFILE) : $(ALL_OBJS)
	@echo "---> Linking objects into one binary: $(OUTFILE)"
ifdef SHARED_LOWLEVEL
	$(GCC) -o $(OUTFILE) $(ALL_OBJS) $(LD_FLAGS) $(GASNET_FLAGS)
else
	$(NVCC) -o $(OUTFILE) $(ALL_OBJS) $(LD_FLAGS) $(GASNET_FLAGS)
endif

$(GEN_OBJS) : %.o : %.cc $(includes)
	$(GCC) -o $@ -c $< $(INC_FLAGS) $(CC_FLAGS)

$(LOW_RUNTIME_OBJS) : %.o : %.cc
	$(GCC) -o $@ -c $< $(INC_FLAGS) $(CC_FLAGS)

$(HIGH_RUNTIME_OBJS) : %.o : %.cc
	$(GCC) -o $@ -c $< $(INC_FLAGS) $(CC_FLAGS)

$(MAPPER_OBJS) : %.o : %.cc
	$(GCC) -o $@ -c $< $(INC_FLAGS) $(CC_FLAGS)

$(GEN_GPU_OBJS) : %.o : %.cu
	$(NVCC) -o $@ -c $< $(INC_FLAGS) $(NVCC_FLAGS)

$(GPU_RUNTIME_OBJS): %.o : %.cu
	$(NVCC) -o $@ -c $< $(INC_FLAGS) $(NVCC_FLAGS)

clean:
	$(RM) -rf $(GEN_OBJS) $(OUTFILE)

# very clean
vclean:
	$(RM) -rf $(ALL_OBJS) $(OUTFILE)
