CC = $(CC_ENV)
#CC = g++
#CC = /usr/bin/g++-4.4
CCFLAGS = -Wall -g
SOURCES =
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA
DELPHES = -I.

# root
ROOTCFLAGS = $(shell root-config --cflags)
ROOTGLIBS = $(shell root-config --glibs)

# boost
BOOSTFLAGS = $(BOOSTFLAGS_ENV)
BOOSTLIBS = $(BOOSTLIBS_ENV)

# File names
SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
EXEC = $(SOURCES:.cc=.exe)

# Main target
# $(EXEC): $(OBJECTS) ExRootTreeReader.o
# $(CC) $(OBJECTS) $(ROOTGLIBS) -L`pwd` -lDelphes ExRootTreeReader.o -o $(EXEC)

# To compile
%.exe: %.o ExRootTreeReader.o
# 	$(CC) $(CCFLAGS) $(ROOTGLIBS) -L`pwd` -lDelphes ExRootTreeReader.o $< -o $@
	$(CC) $(CCFLAGS) $(ROOTGLIBS) $(BOOSTLIBS) -L`pwd` -lDelphes ExRootTreeReader.o $< -o $@

# To obtain object files
%.o: %.cc
# 	$(CC) -c $(CCFLAGS) $(ROOTCFLAGS) $(DELPHES) $< -o $@
	$(CC) -c $(CCFLAGS) $(ROOTCFLAGS) $(BOOSTFLAGS) $(DELPHES) $< -o $@

# needed everywhere: ExRootTreeReader
ExRootTreeReader.o: ExRootAnalysis/ExRootTreeReader.cc ExRootAnalysis/ExRootTreeReader.h
	$(CC) $(CCFLAGS) $(ROOTCFLAGS) $(ROOTGLIBS) $(DELPHES) -c ExRootAnalysis/ExRootTreeReader.cc

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

