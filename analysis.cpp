#include "includes.h"

const string AnalysisRequest::simType = "";
const string AnalysisRequest::analysisType = "";

const string AnalysisDistribution::simType = "TimeDependent";
const string AnalysisDistribution::analysisType = "Distribution";

void AnalysisDistribution::onIteration(SimulationTimeDependent *S)
{
	string datadir;
	
	char Str[512];
	
	if (!directoryExists(outputSubDirectory))
		makeDirectory(outputSubDirectory);
	
	datadir = outputSubDirectory + "iter" + to_string(S->iter);
	
	if (!directoryExists(datadir))
		makeDirectory(datadir);
		
	S->System.writePopulationText(datadir.c_str());
}
