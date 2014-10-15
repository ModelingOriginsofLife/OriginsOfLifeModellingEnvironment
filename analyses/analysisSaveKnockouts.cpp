#include "../includes.h"

void AnalysisSaveKnockouts::onSimulationBegin(SimulationRequest *SR)
{
	SimulationTimeDependent *S = (SimulationTimeDependent *)SR;
	string datadir;
	string outputSubDirectory = S->getSubdirectory(strParams["OUTPUT_DIR"]);
	Simulation *Sim = &S->System;
	char Str[512];
			
	if (!directoryExists(outputSubDirectory))
		makeDirectory(outputSubDirectory);
	
	unordered_set<string>::iterator it;
	
	string filename = outputSubDirectory + "/full_system.txt";
	FILE *f=fopen(filename.c_str(),"wb");
	
	for (it=Sim->knockouts.begin();it!=Sim->knockouts.end();++it)
	{
		fprintf(f,"%s\n",it->c_str());
	}
	fclose(f);
	
	unordered_map<string, int>::iterator it2;
	
	for (it2 = Sim->regionMap.begin();it2 != Sim->regionMap.end();++it2)
	{
		filename = outputSubDirectory + "/region_"+it2->first;
		
		f = fopen(filename.c_str(),"wb");
		Region *R = &Sim->regions[it2->second];
		for (it = R->knockouts.begin();it!=R->knockouts.end();++it)
			fprintf(f,"%s\n",it->c_str());
		fclose(f);
	}
}

AnalysisSaveKnockouts::AnalysisSaveKnockouts()
{	
	simType = "TimeDependent";
	analysisType = "SaveKnockouts";
	beginCallback = 1;

	numParams["PERIOD"] = 100;
	strParams["OUTPUT_DIR"] = "";
}

AnalysisRequest *AnalysisSaveKnockouts::clone()
{
	AnalysisSaveKnockouts *T = new AnalysisSaveKnockouts;
	
	*T = *this;
	
	return T;
}
