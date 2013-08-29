CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA
DELPHES = -I.

# root
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)


# File names
SOURCES = $(wildcard *.cc)
OBJECTS = $(SOURCES:.cc=.o)
EXEC = $(SOURCES:.cc=.exe)

# Main target
# $(EXEC): $(OBJECTS) ExRootTreeReader.o
# 	$(CC) $(OBJECTS) $(ROOTGLIBS) -L`pwd` -lDelphes ExRootTreeReader.o  -o $(EXEC)

# To compile
%.exe: %.o ExRootTreeReader.o
	$(CC) $(ROOTGLIBS) -L`pwd` -lDelphes ExRootTreeReader.o $< -o $@

# To obtain object files
%.o: %.cc
	$(CC) -c $(CC_FLAGS) $(ROOTCFLAGS) $(ROOTGLIBS) $(DELPHES) -L`pwd` -lDelphes $< -o $@

# needed everywhere: ExRootTreeReader
ExRootTreeReader.o: ExRootAnalysis/ExRootTreeReader.cc ExRootAnalysis/ExRootTreeReader.h
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(ROOTGLIBS) $(DELPHES) -c ExRootAnalysis/ExRootTreeReader.cc

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

