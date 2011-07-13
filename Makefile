# Petter Strandmark 2011

# You might have to change these to g++-4.5 or something similar
COMPILER = g++ -c
LINKER = g++
# Temporary directory to store object files; needs to exist
OBJDIR = ~/obj/surfaces/
# Location of Clp if not installed where gcc can find it
CLPLIBDIR = ~/Programming/coin-Clp/lib/
CLPINCLUDEDIR = ~/Programming/coin-Clp/include/
# Locations of Armadillo
ARMADILLOINCLUDEDIR = ~/Programming/libraries/armadillo-1.0.0/include

##################################################################################

OBJFILES := $(OBJDIR)Petter-Color.o $(OBJDIR)curv3d.o $(OBJDIR)curv3d_lp.o $(OBJDIR)qpbo_3dcurv.o $(OBJDIR)mesh3d.o 
QPBOOBJ := $(OBJDIR)QPBO.o $(OBJDIR)QPBO_extra.o $(OBJDIR)QPBO_maxflow.o $(OBJDIR)QPBO_postprocessing.o 
LIBDIR = $(OBJDIR)
INCLUDE = -I software/cv/curv3d  -I software/thirdparty/HOCR -I software/thirdparty/QPBO -I software/cv/mesh -I software/cv/points
OPTIONS = -std=c++0x
TARGET = bin/surfaces_gcc

$(TARGET): $(OBJFILES) $(QPBOOBJ)
	$(LINKER) -L $(CLPLIBDIR) $(OBJFILES) $(QPBOOBJ) -o $(TARGET) -lClp -lCoinUtils -llapack -lblas
	
clean:
	rm -f $(TARGET)*
	rm -f $(OBJDIR)*.o
	rm -f $(OBJDIR)*.la
	
$(OBJDIR)main_program.o: source/main_program.cpp
	$(COMPILER) $(OPTIONS) $(INCLUDE) source/main_program.cpp -o $(OBJDIR)main_program.o

$(LIBDIR)libpetter.a : $(SUBMOBJ) $(OBJDIR)Petter-Color.o
	ar rcs $(LIBDIR)libpetter.a $(OBJDIR)Petter-Color.o 

$(OBJDIR)Petter-Color.o: software/cv/curv3d/Petter-Color.cc software/cv/curv3d/Petter-Color.h
	$(COMPILER) $(OPTIONS) $(INCLUDE) software/cv/curv3d/Petter-Color.cc -o $(OBJDIR)Petter-Color.o	
	
$(OBJDIR)%.o: software/cv/curv3d/%.cc software/cv/curv3d/3dcurv.hh software/cv/mesh/mesh3d.hh software/cv/mesh/lp_constraints.hh
	$(COMPILER) $(OPTIONS) $(INCLUDE) -I $(ARMADILLOINCLUDEDIR) -I $(CLPINCLUDEDIR) $< -o $@

$(OBJDIR)%.o: software/cv/points/%.cc software/cv/curv3d/points_lp.hh software/cv/mesh/mesh3d.hh
	$(COMPILER) $(OPTIONS) $(INCLUDE) -I $(ARMADILLOINCLUDEDIR) -I $(CLPINCLUDEDIR) $< -o $@
	
$(OBJDIR)%.o: software/cv/mesh/%.cc software/cv/mesh/mesh3d.hh 
	$(COMPILER) $(OPTIONS) $(INCLUDE) -I $(ARMADILLOINCLUDEDIR) -I $(CLPINCLUDEDIR) $< -o $@
	
$(OBJDIR)%.o: software/thirdparty/QPBO/%.cpp 
	$(COMPILER) $(OPTIONS) $(INCLUDE) $< -o $@

%.hh :
	#nothing
	
