# -*- makefile-gmake -*-

# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

CFLAGS = -O3 -g -I.
PREFIX = ${HOME}/infer/gkylsoft

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_STLIB = ${PREFIX}/gkylzero/lib/libgkylzero.a
CUDA_LIB = 

ifeq ($(CC), nvcc)
	CUDA_LIB = -lcublas -lcusparse -lcusolver
endif

G0_LIBS = ${G0_STLIB} ${CUDA_LIB} -lm -lpthread

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = ${PREFIX}/OpenBLAS/include
LAPACK_LIB = ${PREFIX}/OpenBLAS/lib/libopenblas.a

# SuperLU includes and librararies
SUPERLU_INC = ${PREFIX}/superlu/include
SUPERLU_LIB_DIR = ${PREFIX}/superlu/lib64
SUPERLU_LIB = ${PREFIX}/superlu/lib64/libsuperlu.a

# determine OS we are running on
UNAME = $(shell uname)

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_INC = . # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
endif

INCLUDES = -I${G0_INC_DIR} -I${LAPACK_INC} -I${SUPERLU_INC}

all: jets

jets: jets.c
	 ${CC} ${CFLAGS} ${INCLUDES} jets.c -o jets -L${G0_LIB_DIR} ${G0_LIBS} ${LAPACK_LIB} ${SUPERLU_LIB}

clean:
	rm -rf jets jets.dSYM

