################################################################################
#
# Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
#
################################################################################

# Flags to detect 32-bit or 64-bit OS platform
OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

#OS_SIZE = 64
#OS_ARCH = x86_64

# Location of the CUDA Toolkit binaries and libraries
CUDA_PATH       ?= /usr/local/cuda-5.5
CUDA_INC_PATH   ?= $(CUDA_PATH)/include
CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib64

# Common binaries
NVCC            ?= $(CUDA_BIN_PATH)/nvcc
GCC             ?= g++

# Extra user flags
EXTRA_NVCCFLAGS ?=
EXTRA_LDFLAGS   ?=

# CUDA code generation flags
GENCODE_SM10    := -gencode arch=compute_10,code=sm_10
GENCODE_SM20    := -gencode arch=compute_20,code=sm_20
GENCODE_SM30    := -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35
GENCODE_FLAGS   := $(GENCODE_SM10) $(GENCODE_SM20) $(GENCODE_SM30)

# OS-specific build flags
LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart
CCFLAGS   := -m64

# OS-architecture specific flags
NVCCFLAGS := -m64

# Debug build flags
ifeq ($(dbg),1)
      CCFLAGS   += -g
      NVCCFLAGS += -g -G
endif


# Common includes and paths for CUDA
INCLUDES      := -I$(CUDA_INC_PATH) -I.

# Additional parameters
MAXRREGCOUNT  :=  -po maxrregcount=16

# Target rules
all: build

build: blackscholes

blackscholes.o: blackscholes.cu
	$(NVCC) $(NVCCFLAGS) $(EXTRA_NVCCFLAGS) $(GENCODE_FLAGS) $(MAXRREGCOUNT) $(INCLUDES) -o $@ -c $<

blackscholes_gold.o: blackscholes_gold.cpp
	$(GCC) $(CCFLAGS) $(INCLUDES) -o $@ -c $<

blackscholes: blackscholes.o blackscholes_gold.o
	$(GCC) $(CCFLAGS) -o $@ $+ $(LDFLAGS) $(EXTRA_LDFLAGS)

run: build
	./blackscholes

clean:
	rm -f blackscholes blackscholes.o blackscholes_gold.o

