CC = g++
CFLAGS = -O6 -m64 -std=c++0x -w -fpermissive -g
LIBS = -lm -lboost_system -lboost_filesystem -larmadillo

ANALYSES = analyses/analysisDistribution.o analyses/analysisHeredity.o analyses/analysisReactionNet.o analyses/analysisOutputTimeseries.o analyses/analysisSaveKnockouts.o
SIMULATIONS = simulations/simulationTimeDependent.o simulations/simulationReactionChain.o

OBJ = simulation.o random.o output.o weighted_tree.o file_parser.o \
reaction_parser.o reaction_engine.o analysis.o linear.o $(ANALYSES) $(SIMULATIONS) 

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

all: oolme

oolme: $(OBJ)
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS)

.PHONY: clean

clean:
	rm -f *.o analyses/*.o simulations/*.o *~ oolme
