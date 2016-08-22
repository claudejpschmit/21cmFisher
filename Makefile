####################################################
# This makefile creates all the objects in a build #
# directory and links them into executables.       #
####################################################

MDIR := $(shell pwd)
WRKDIR = $(MDIR)/build
LIBRARIES = lib/
.base:
	if ! [ -e $(WRKDIR) ]; then mkdir $(WRKDIR) ; mkdir $(WRKDIR)/lib; fi;
	touch build/.base

# determines the location of the .cpp files
vpath %.cpp src:\
	main:\
	$(LIBRARIES)ALGLIB_source:\
	$(LIBRARIES)GLOBAL21CM_source:\
	$(LIBRARIES)ODE_source

vpath %.o build
vpath .base build

# Compiler required for c++ code.
# including -ffast-math may not be as bad as anticipated.
CXX = g++ -Wall -std=c++11 -ffast-math -s -Wno-deprecated -fopenmp 

OPTFLAG = -O4
ARMAFLAGS = -larmadillo
GSLFLAGS = -lgsl -lgslcblas
WIGNERFLAGS = -lwignerSymbols
CUBAFLAGS = -lcuba

# This decreases the precision of the bessel functions significantly if 
# added to the compilation of the files containing boost->sph_bess(l,x).
OPTFLAG_CLASS = -ffast-math
OMPFLAG = -fopenmp
CCFLAG = -g -fPIC
LDFLAG = -g -fPIC

# Header files for local libraries
INCLUDES = -I../include
INCLUDES += -I../$(LIBRARIES)ALGLIB_include
INCLUDES += -I../$(LIBRARIES)GLOBAL21CM_include
INCLUDES += -I../$(LIBRARIES)ODE_include

INCLUDES += -I/usr/include/boost
BOOSTFLAGS = -lboost_filesystem -lboost_system

# Including GSL
INCLUDES += -I/usr/include
LINKER += -L/usr/lib

%.o: %.cpp .base
	cd $(WRKDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(INCLUDES) -c ../$< -o $*.o $(ARMAFLAGS) $(GSLFLAGS) $(CUBAFLAGS) -lm

ALGLIB = alglibinternal.o alglibmisc.o ap.o dataanalysis.o diffequations.o\
		 fasttransforms.o integration.o interpolation.o linalg.o optimization.o\
		 solvers.o specialfunctions.o statistics.o
GLOBAL21CM = dnumrecipes.o dcomplex.o dcosmology.o astrophysics.o twentyonecm.o spline.o spline2D.o
MAIN = Main.o
ANALYSE = Analyser.o Analyse.o
SORTING = SortFiles.o
INIREADER = iniReader.o

SOURCE = CosmoBasis.o Models.o AnalysisInterface.o Cosmology3D.o Tomography2D.o\
		 IntensityMapping.o HighZAnalysis.o Fisher.o Fisher1.o Fisher_Santos.o Integrator.o\
		 CAMB_interface.o ARES_interface.o Global21cmInterface.o Zygelman.o
ODE = ODEs.o ODE_Solver.o
BISPECTRUM = Bispectrum_main.o Bispectrum.o LISW.o Bispectrum_Fisher.o

# ---------------------------------------------------------------------------------------------------------------#

all: calc analyse sortFiles bispectrum 

analyse: $(ALGLIB) $(INIREADER) $(ANALYSE) 
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(LINKER)\
		-o analyse $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)

calc: $(SOURCE) $(INIREADER) $(MAIN) $(ALGLIB) $(GLOBAL21CM)
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(LINKER)\
		-o calc $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS) $(BOOSTFLAGS)

sortFiles: $(SORTING)
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(LINKER)\
		-o sortFiles $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS)

bispectrum: $(BISPECTRUM) $(INIREADER) $(SOURCE) $(ALGLIB) $(GLOBAL21CM) $(ODE)
	cd $(MDIR);$(CXX) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(LINKER)\
		-o bispectrum $(addprefix build/, $(notdir $^)) -lm $(ARMAFLAGS) $(GSLFLAGS) $(WIGNERFLAGS) $(CUBAFLAGS) $(BOOSTFLAGS)


clean: .base
	rm -rf $(WRKDIR);
	rm calc;
	rm analyse;
	rm sortFiles;
	rm bispectrum;

