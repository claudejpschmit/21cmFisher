#################################################
# This file should just make a bunch of objects,#
# to simulate what class does.                  #
#################################################

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
LIBRARIES = lib/
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

vpath %.cpp src:main:$(LIBRARIES)ALGLIB_source:$(LIBRARIES)GLOBAL21CM_source
#:$(LIBRARIES)ODEsolver_source
#vpath %.c $(LIBRARIES)CLASS_source:$(LIBRARIES)CLASS_tools:$(LIBRARIES)CLASS_main
vpath %.o build
vpath .base build

# Compiler required for c++ code.
# including -ffast-math may not be as bad as anticipated.
CXX = g++ -Wall -std=c++11 -ffast-math -s -Wno-deprecated -fopenmp 
# Compiler required for c code.
CC = gcc -Wall -s
# Compiler required for fortran RECFAST
FF = gfortran

OPTFLAG = -O4
ARMAFLAGS = -larmadillo
GSLFLAGS = -lgsl -lgslcblas

# This decreases the precision of the bessel functions significantly if 
# added to the compilation of the files containing boost->sph_bess(l,x).
OPTFLAG_CLASS = -ffast-math
OMPFLAG = -fopenmp
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

# leave blank to compile without HyRec, or put path to HyRec directory
# (with no slash at the end: e.g. CLASS_hyrec or ../CLASS_hyrec)
HYREC = $(LIBRARIES)CLASS_hyrec 
# automatically add external programs if needed. First, initialize to blank.
EXTERNAL =

CCFLAG += -D__CLASSDIR__='"$(MDIR)"'
INCLUDES = -I../include
#INCLUDES += -I../$(LIBRARIES)CLASS_include
INCLUDES += -I../$(LIBRARIES)ALGLIB_include
INCLUDES += -I../$(LIBRARIES)GLOBAL21CM_include
#INCLUDES += -I../$(LIBRARIES)ODEsolver_include

# These lines seem to be unnecessary, but I leave them in anyways.
INCLUDES += -I/usr/include/boost
LINKER = -L/usr/include/boost #-lboost_filesystem

# Linking GSL
INCLUDES += -I/usr/include
LINKER += -L/usr/lib -lgsl -lgslcblas


# eventually update flags for including HyRec
ifneq ($(HYREC),)
vpath %.c $(HYREC)
CCFLAG += -DHYREC
#LDFLAGS += -DHYREC
#INCLUDES += -I../$(LIBRARIES)CLASS_hyrec
EXTERNAL += hyrectools.o helium.o hydrogen.o history.o
endif

%.o: %.cpp .base
	cd $(WRKDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o $(ARMAFLAGS) $(GSLFLAGS)

# This line creates the CLASS objects.
%.o: %.c .base
	cd $(WRKDIR);$(CC) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o

TOOLS = growTable.o dei_rkck.o sparse.o evolver_rkck.o  evolver_ndf15.o arrays.o parser.o quadrature.o hyperspherical.o common.o
SOURCE = input.o background.o thermodynamics.o perturbations.o primordial.o nonlinear.o transfer.o spectra.o lensing.o
CLASS = class.o
OUTPUT = output.o
ALGLIB = alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o fasttransforms.o integration.o interpolation.o linalg.o optimization.o solvers.o specialfunctions.o statistics.o
GLOBAL21CM = dnumrecipes.o dcomplex.o dcosmology.o astrophysics.o twentyonecm.o spline.o spline2D.o
ODE = ODE_Solver.o ODEs.o
SRC = Integrator.o CosmoBasis.o CosmologyCalculatorClass.o CosmologyWriterClass.o FisherClass.o Engine.o ClassEngine.o CAMB_interface.o Global21cmInterface.o SanityChecker.o LevinIntegrator.o ARES_interface.o
MAIN = Main.o
ANALYSE = Analyser.o Analyse.o
TOMOGRAPHY = calc_tomography.o Tomography.o 

RESTRUCT = CosmoBasis.o Models.o AnalysisInterface.o Cosmology3D.o Tomography2D.o Fisher.o Fisher1.o Fisher_Santos.o Integrator.o CAMB_interface.o ARES_interface.o Global21cmInterface.o


all: calc analyse 
#class_test analyse calc_tomography 

analyse: $(ALGLIB) $(ANALYSE)
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o analyse $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS)

#calc_tomography: $(TOMOGRAPHY) $(SRC) $(SOURCE) $(TOOLS) $(OUTPUT) $(EXTERNAL) $(ALGLIB) $(GLOBAL21CM)
#	cd $(MDIR);$(CXX) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o calc_tomography $(addprefix# build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS)


#calc: $(SRC) $(SOURCE) $(TOOLS) $(OUTPUT) $(EXTERNAL) $(ALGLIB) $(GLOBAL21CM) $(ODE) $(MAIN) 
#	cd $(MDIR);$(CXX) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o calc $(addprefix build/, $(#notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS)
calc: $(RESTRUCT) $(MAIN) $(ALGLIB) $(GLOBAL21CM)
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o calc $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS)

#class_test: $(SOURCE) $(TOOLS) $(OUTPUT) $(EXTERNAL) $(CLASS) 
#	cd $(MDIR);$(CC) $(OPTFLAG) $(OPTFLAG_CLASS) $(OMPFLAG) $(LDFLAG) $(LINKER) -o class_test $(addprefix build#/, $(notdir $^)) -lm

#install:
#	$(FF) -o $(LIBRARIES)GLOBAL21CM_dependencies/RECFAST_CODE/recfast $(LIBRARIES)GLOBAL21CM_dependencies/RECFA#ST_CODE/recfast.for
#	cd $(MDIR)
#	make all

#clean_integration: .base
#	rm $(WRKDIR)/Integrator.o;
#	rm $(WRKDIR)/Main.o;
#	rm $(WRKDIR)/CosmologyCalculatorClass.o;
#	rm calc

clean: .base
	rm -rf $(WRKDIR);
#	rm -f libclass.a;
	rm calc;
#	rm class_test;
#	rm analyse;
#	rm calc_tomography

