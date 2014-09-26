#include "../includes.h"

void sumOverCompounds(Node *N, unordered_map<char, double> &symbols);

void sumOverCompounds(Node *N, unordered_map<char, double> &symbols)
{
	if (N==NULL) return;
	
	if (N->isleaf)
	{
		for (int i=0;i<N->value.length();i++)
		{
			symbols[N->value[i]]+=N->weight;
		}
	}
	else
	{
		sumOverCompounds(N->child[0], symbols);
		sumOverCompounds(N->child[1], symbols);
	}
}

void AnalysisOutputTimeseries::onIteration(SimulationRequest *SR)
{
	SimulationTimeDependent *S = (SimulationTimeDependent *)SR;
	string datadir;
	string outputSubDirectory = strParams["OUTPUT_DIR"];
	char Str[512];
	FILE *f;
	
	if (!directoryExists(outputSubDirectory))
		makeDirectory(outputSubDirectory);
	
	datadir = outputSubDirectory;
	
	unordered_map<string, int>::iterator it;
	int i,j;
	
	int numValues = 1;
	vector<double> values;
	vector<string> cNames;
	
	if (strParams["OUTPUT_NUMBER"] != "false")
	{
		numValues += 1;
	}
	
	if (strParams["OUTPUT_LENGTH"] != "false")
	{
		numValues += 1;
	}

	if (strParams["OUTPUT_DERIVED"] != "false")
	{
	}
	
	if (strParams["OUTPUT_SYMBOLS"] != "false")
	{
		numValues += S->System.parentChem->L.atomTypes.size();
	}
		
	if (strParams["OUTPUT_COMPOUNDS"].length())
	{
		f=fopen(strParams["OUTPUT_COMPOUNDS"].c_str(),"rb");
		if (f!=NULL)
		{
			while (fscanf(f,"%s\n",Str)!=EOF)
			{
				string sStr = Str;
				
				cNames.push_back(sStr);
				numValues++;
			}
			
			fclose(f);
		}
	}
		
	values.resize(numValues,0);
		
	for (it=S->System.regionMap.begin();it!=S->System.regionMap.end();++it)
	{
		i = it->second;
		Region *R = &S->System.regions[i];
		j=0;
		
		values[j++] = S->iter;
				
		if (strParams["OUTPUT_NUMBER"] != "false")
		{
			if (R->population.Root != NULL)
				values[j++] += R->population.Root->weight;
		}
	
		if (strParams["OUTPUT_LENGTH"] != "false")
		{
			if (R->population.Root != NULL)
				values[j++] += R->population.Root->lweight;
		}

		if (strParams["OUTPUT_DERIVED"] != "false")
		{
		}
	
		if (strParams["OUTPUT_SYMBOLS"] != "false")
		{
			unordered_map<char, double> symbols;
			sumOverCompounds(R->population.Root, symbols);
			
			for (int k=0;k<S->System.parentChem->L.atomTypes.size();++k)
			{
				values[j++] += symbols[S->System.parentChem->L.atomTypes[k]];
			}
		}
			
		for (int k=0;k<cNames.size();k++)
		{
			unordered_map<string, Node*>::iterator it = R->population.accessHash.find(cNames[k]);
			
			if (it != R->population.accessHash.end())
			{
				values[j+k] += it->second->weight;
			}
		}
		j+=cNames.size();
		
		if (strParams["PER_REGION"] != "false")
		{
			string fname = datadir + it->first + ".csv";
			
			f=fopen(fname.c_str(),"a");
			
			for (int k=0;k<values.size();k++)
				fprintf(f,"%.9g%s", values[k], k<values.size()-1 ? ", " : "\n");
			
			fclose(f);
			
			std::fill(values.begin(), values.end(), 0);
		}
	}	
	
	if (strParams["PER_REGION"] == "false")
	{
		string fname = datadir + "total.csv";
			
		f=fopen(fname.c_str(),"a");
			
		for (int k=0;k<values.size();k++)
			fprintf(f,"%.9g%s", values[k], k<values.size()-1 ? ", " : "\n");
			
		fclose(f);
	}
		
	values.clear();
}

AnalysisOutputTimeseries::AnalysisOutputTimeseries()
{	
	simType = "TimeDependent";
	analysisType = "OutputTimeseries";
	iterateCallback = 1;

	numParams["PERIOD"] = 1;
	strParams["PER_REGION"] = "true";
	strParams["OUTPUT_DIR"] = "";

	strParams["OUTPUT_NUMBER"] = "false";
	strParams["OUTPUT_LENGTH"] = "false";
	strParams["OUTPUT_DERIVED"] = "false";	
	strParams["OUTPUT_SYMBOLS"] = "false";
	strParams["OUTPUT_COMPOUNDS"] = "";
}

AnalysisRequest *AnalysisOutputTimeseries::clone()
{
	AnalysisOutputTimeseries *T = new AnalysisOutputTimeseries;
	
	*T = *this;
	
	return T;
}
