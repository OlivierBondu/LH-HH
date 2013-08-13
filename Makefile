CC        = g++
CCFLAGS   = -Wall -g
SOURCES   =
ROOTFLAGS = `root-config --cflags`
ROOTLIBS  = `root-config --libs --ldflags`
ROOFITLIBS = -lRooFit -lRooFitCore -lMinuit -lFoam
ROOSTATSLIBS = -lRooStats
TMVA = -L${ROOTSYS}lib -lTMVA
DELPHES = -I.

all: selection.exe
#all:ExRootTreeReader.o selection.o

ExRootTreeReader.o: ExRootAnalysis/ExRootTreeReader.cc ExRootAnalysis/ExRootTreeReader.h
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(DELPHES) -c ExRootAnalysis/ExRootTreeReader.cc

selection.o: selection.cc
	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(DELPHES) -c selection.cc	

selection.exe: selection.o ExRootTreeReader.o
	$(CC) $(ROOTLIBS) -L`pwd` -lDelphes selection.o ExRootTreeReader.o -o $@

#selection.exe:
#	$(CC) $(CCFLAGS) $(ROOTFLAGS) $(DELPHES) $(ROOTLIBS) -o $@ selection.cc ExRootAnalysis/ExRootTreeReader.cc

clean:
	rm *.exe; rm *.o
