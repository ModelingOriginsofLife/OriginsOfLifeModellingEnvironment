#include "../includes.h"

void AnalysisHeredity::onIteration(SimulationRequest *SR)
{
}

AnalysisHeredity::AnalysisHeredity()
{	
}

AnalysisRequest *AnalysisHeredity::clone()
{
	AnalysisHeredity *T = new AnalysisHeredity;
	
	*T = *this;
	
	return T;
}
