#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <algorithm>
#include <sstream>

#include <boost/filesystem.hpp>

#include <armadillo>

using namespace std;
using namespace arma;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <memory.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "weighted_tree.h"
#include "reaction_parser.h"
#include "simulation.h"
#include "analysis.h"
#include "file_parser.h"
#include "random.h"
#include "output.h"
#include "linear.h"

#include "simulations/simulationTimeDependent.h"
#include "simulations/simulationReactionChain.h"

#include "analyses/analysisDistribution.h"
#include "analyses/analysisOutputTimeseries.h"
#include "analyses/analysisHeredity.h"
#include "analyses/analysisReactionNet.h"
#include "analyses/analysisSaveKnockouts.h"
