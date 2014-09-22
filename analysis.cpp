#include "includes.h"

vector <AnalysisRequest*> registeredAnalyses;

void registerAnalyses()
{
	registeredAnalyses.push_back(new AnalysisDistribution);
}

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

AnalysisRequest::AnalysisRequest()
{
	numParams["PERIOD"] = 1;
}

void AnalysisRequest::parse(ifstream &file)
{
	int istop = 0;
	string str, substr;
	int place = 0;
	
	do
	{
		getCollapsedLine(file,str);				
		str = getFirstToken(str, substr);
		if (substr == "ENDPARAMS") istop = 1;
		else
		{
			if (((str[0] >= '0')&&(str[0] <= '9')) || (str[0] == '-') || (str[0] == '.')) // Numerical parameter
			{
				numParams[substr] = atof(str.c_str());
			}
			else
			{
				if (str[0] == '\"') // remove a leading quote mark
				{
					str = str.substr(1,str.length()-1);
					if (str[str.length()-1] == '\"') // closing the quotes is natural, so remove trailing if present
						str = str.substr(0,str.length()-1);					
				}
				
				strParams[substr] = str;
			}
		}
	} while (!istop);
}
