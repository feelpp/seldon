# C++ compiler and options
CC = g++
FLAGS = -g

# file to compile and executable
#MAIN_CC = direct_test.cpp
MAIN_CC = iterative_test.cpp
BIN = main.x

# directory where Seldon is 
INCLUDE := -I../../

# basic libraries : BLAS, LAPACK
LIB := -llapack -L$(HOME)/CBLAS -lcblas -lblas -lg2c

# Mumps
SOLVE = $(HOME)/SOLVE/lib_PC
DIR_MUMPS = $(SOLVE)/MUMPS_4.7.3
INCLUDE := -I$(DIR_MUMPS)/include -I$(DIR_MUMPS)/libseq \
           -I$(DIR_METIS)/metis-4.0/Lib $(INCLUDE)
LIB := -L$(DIR_MUMPS)/lib -ldmumps -lzmumps -lpord \
     -L$(DIR_MUMPS)/libseq -lmpiseq \
    -L$(HOME)/common/lib -lgfortran -lpthread \
    -L$(SOLVE)/metis-4.0 -lmetis $(LIB)

# UmfPack
DIR_UMF = $(SOLVE)/UMFPACKv5.02
INCLUDE := -I$(DIR_UMF)/AMD/Include -I$(DIR_UMF)/UFconfig \
           -I$(DIR_UMF)/UMFPACK/Include $(INCLUDE)
LIB := -L$(DIR_UMF)/UMFPACK/Lib -lumfpack -L$(DIR_UMF)/AMD/Lib -lamd $(LIB)

# SuperLU
DIR_SUPERLU = $(SOLVE)/SuperLU
INCLUDE := -I$(DIR_SUPERLU)/SRC $(INCLUDE)
LIB := -L$(DIR_SUPERLU) -lsuperlu $(LIB)

# compilation directive
all :
	$(CC) $(FLAGS) $(MAIN_CC) -o test.x $(INCLUDE) $(LIB)
	mv test.x $(BIN)
