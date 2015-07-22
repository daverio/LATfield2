#---------------- parameters -----------------------------------------------#
EXEC                      := poisson wave_test gettingStarted #fft_test
TEST_WITH_EVERY_BUILD     := false
COMPILER                  := CC
INC                       := -I./../
#LIB                       := -lfftw3 -lm
LIB_OPENACC  			  := -lcufftw -lcufft
LIB_CPU                   := -lfftw3 -lm
DEF_LATFIELD_OPENACC      := -DFFT3D_ACC
DEF_LATFIELD_CPU         := -DFFT3D
ifdef DEBUG
OPT_OPENACC               :=
OPT_CPU                   :=
CFLAG                     := -g -std=c++11 -Mbounds -traceback -w
else
OPT_OPENACC               := -O3 -acc -ta=tesla:cc35,pin -Mcuda #-Minfo=accel
OPT_CPU                   := -O3
#CFLAG                     := -std=c++11 -w
CFLAG                     := -std=c++11 -w -DNO_OUTPUT -DLOOP_COUNT=4 -VH15.7

endif

#---------------- setup -----------------------------------------------------#
HEADER                    := $(wildcard ../*.hpp) $(wildcard LATfield2d/*.h)
EXEC_CPU                  := $(addsuffix _cpu,$(EXEC))
EXEC_OPENACC              := $(addsuffix _openacc,$(EXEC))
DEF                       := $(DEF_LATFIEILD)
ifeq ($(TEST_WITH_EVERY_BUILD),true)
ADDITIONAL_BUILD_TARGETS=tests
else
ADDITIONAL_BUILD_TARGETS=
endif
.PHONY: all tests clean openacc cpu

#---------------- target build rules ----------------------------------------#
#mpic++ main.cpp -I./LATfield2d -I../local/include/gsl/  -DFFT3D -DPHINONLINEAR -DCHECK_B -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -std\
# =c++0x -O3  -o test

all: openacc cpu ${ADDITIONAL_BUILD_TARGETS}

openacc: ${EXEC_OPENACC}

cpu: ${EXEC_CPU}

tests: ${EXEC_OPENACC} ${EXEC_CPU}
	./RunTest.sh

%_cpu: %.cpp $(HEADER) makefile
	$(COMPILER) $< $(INC) $(DEF_LATFIELD_CPU) $(LIB_CPU) $(OPT_CPU) $(CFLAG) -o $@

%_openacc: %.cpp $(HEADER) makefile
	$(COMPILER) $< $(INC) $(DEF_LATFIELD_OPENACC) $(LIB_OPENACC) $(OPT_OPENACC) $(CFLAG) -o $@

clean:
	rm -f $(EXEC_CPU) $(EXEC_OPENACC) *.o *~
