#include "../includes.h"

void AnalysisDistribution::onIteration(SimulationRequest *SR)
{
	SimulationTimeDependent *S = (SimulationTimeDependent *)SR;
	string datadir;
	string outputSubDirectory = strParams["OUTPUT_DIR"];
	char Str[512];
		
	if (!directoryExists(outputSubDirectory))
		makeDirectory(outputSubDirectory);
	
	datadir = outputSubDirectory + "iter" + to_string(S->iter);
	
	if (!directoryExists(datadir))
		makeDirectory(datadir);
		
	S->System.writePopulationText(datadir.c_str());
}

AnalysisDistribution::AnalysisDistribution()
{	
	simType = "TimeDependent";
	analysisType = "Distribution";
	iterateCallback = 1;

	numParams["PERIOD"] = 100;
	strParams["OUTPUT_DIR"] = "";
}

AnalysisRequest *AnalysisDistribution::clone()
{
	AnalysisDistribution *T = new AnalysisDistribution;
	
	*T = *this;
	
	return T;
}
