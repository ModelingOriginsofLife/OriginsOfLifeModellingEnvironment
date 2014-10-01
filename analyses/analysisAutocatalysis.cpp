#include "../includes.h"

void AnalysisAutocatalysis::onSimulationEnd(SimulationRequest *SR)
{
}

void AnalysisAutocatalysis::onReaction(vector<string> &reac, vector<string> &prod, Region *Reg)
{
}

AnalysisAutocatalysis::AnalysisAutocatalysis()
{	
}

AnalysisRequest *AnalysisAutocatalysis::clone()
{
	AnalysisAutocatalysis *T = new AnalysisAutocatalysis;
	
	*T = *this;
	
	return T;
}
