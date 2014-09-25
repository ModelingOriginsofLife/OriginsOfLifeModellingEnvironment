CC = g++
CFLAGS = -O6 -m64 -std=c++11 -w -fpermissive
LIBS = -lm -lboost_system -lboost_filesystem -larmadillo

ANALYSES = analyses/analysisDistribution.o analyses/analysisHeredity.o analyses/analysisReactionNet.o
SIMULATIONS = simulations/simulationTimeDependent.o simulations/simulationReactionChain.o

OBJ = simulation.o random.o output.o weighted_tree.o file_parser.o \
reaction_parser.o reaction_engine.o analysis.o $(ANALYSES) $(SIMULATIONS) 

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all: oolme

oolme: $(OBJ)
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o analyses/*.o simulations/*.o *~ oolme
