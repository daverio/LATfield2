# programming environment
COMPILER     := CC
INCLUDE      := -I../
LIB          := -lfftw3 -lm -lhdf5 -lgsl -lgslcblas

# target and source
EXEC         := mulit
SOURCE       := multigrid.cpp
HEADERS      := $(wildcard *.hpp ../*.hpp ../*.h)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DHDF5 -DMULTIGRID

# optional compiler settings (LATfield2)
DLATFIELD2   += -DH5_HAVE_PARALLEL -DH5_HAVE_PIXIE
#DLATFIELD2   += -DEXTERNAL_IO

# compiler settings
DNBODY := -DDEBUG

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DNBODY) $(INCLUDE) $(LIB)
