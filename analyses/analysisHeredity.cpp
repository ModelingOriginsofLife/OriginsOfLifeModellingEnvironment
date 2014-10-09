#include "../includes.h"

int HeredityRunData::getFeatureCount(int vCols, double threshold, unordered_map<string, int> &cIndex)
{
	int maxval = frames.size();
	
	mat data(maxval, vCols, fill::ones);
	
	data = -1 * data;
	
	mat cov;
	
	for (int j=0;j<frames.size();j++)
	{
		for (unordered_set<string>::iterator it = frames[j].compounds.begin(); it != frames[j].compounds.end(); ++it)
		{
			int idx = cIndex[*it]; 
	
			if (idx >= 0)
			{
				data(j,idx) = 1;
			}
		}
	}
	
	cov = data.t() * data / maxval;
	
	for (int i=0;i<cov.n_rows;i++)
	{
		bool stop=false;
		for (int j=0;(j<i)&&(!stop);j++)
		{
			if (fabs(cov(i,j)) > threshold)
			{
				cov.shed_row(i);
				cov.shed_col(i);
				i--;
				stop=true;
			}
		}
	}
	
	return cov.n_rows;
}

double HeredityRunData::getEntropy()
{
	unordered_map<string, int> counts;
	
	for (int i=0;i<frames.size();++i)
	{
		for (unordered_set<string>::iterator it = frames[i].compounds.begin();it != frames[i].compounds.end(); ++it)
		{
			counts[*it]++;
		}
	}
	
	double H = 0, total = frames.size();
	
	unordered_map<string, int>::iterator it;
	
	for (it=counts.begin();it!=counts.end();++it)
	{
		double P = it->second/total;
		
		H -= P*log(P);
	}
	
	return H;
}

void AnalysisHeredity::onIteration(SimulationRequest *SR)
{
	SimulationTimeDependent *S = (SimulationTimeDependent*)SR;
	Simulation *Sim = &S->System;
	double threshold = numParams["DETECTION_THRESHOLD"];
	
	// Regions are all lumped together - lets do this differently in the future, or give multiple options
	HereditySnapshot HS;
	ReactionRule RRule;
	
	if (strParams["RANDOM_KNOCKOUTS"] != "false")
	{
		RRule.rule = strParams["KNOCKOUT_RULE"];
		RRule.parseRule(Sim->parentChem->L);
	}
	
	for (int i=0;i<Sim->regions.size();++i)
	{
		Region *R = &Sim->regions[i];		
		unordered_map<string, Node*>::iterator it;
				
		for (it=R->population.accessHash.begin();it!=R->population.accessHash.end();++it)
		{
			bool accept = true;
			// Exclude compounds that match the knockout rule
			if (strParams["RANDOM_KNOCKOUTS"] != "false")
			{
				if (RRule.matchCompound(it->first, Sim->parentChem->L))
				{
					accept = false;
				}
			}

			if (accept)
				if (it->second->weight >= threshold)
				{
					if (!cIndex.count(it->first))
					{
						cIndex[it->first] = cidx++;
					}
								
					HS.compounds.insert(it->first);
				}			
		}
	}
	
	runs[simidx].frames.push_back(HS);
}

void AnalysisHeredity::onSimulationBegin(SimulationRequest *SR)
{
	SimulationTimeDependent *S = (SimulationTimeDependent*)SR;
	ChemistryComputation *C = S->System.parentChem;
	
	iterCounter = 0;
	
	if (strParams["RANDOM_KNOCKOUTS"]!="false")
	{
		int nKnock = numParams["KNOCKOUT_COUNT"];
		double wildlength = numParams["KNOCKOUT_MEAN_WILDCARD_LENGTH"];
		string rule = strParams["KNOCKOUT_RULE"];
		
		for (int i=0;i<nKnock;i++)
		{
			string compound = generateRandomCompound(rule, *C, wildlength);
			
			if (compound.length())
			{
				S->System.knockouts.insert(compound);
			}
		}
	}
	
	HeredityRunData HR;	
	runs.push_back(HR);
}	

void AnalysisHeredity::PCAHeredity()
{
	int vCols = cidx;
	int maxval = endStates.frames.size();
	vector<int> cCount;
	vector<int> cMap;
/*	
	cCount.resize(cidx, 0);	
	
	for (int j=0;j<maxval;j++)
	{
		for (unordered_set<string>::iterator it = endStates.frames[j].compounds.begin(); it != endStates.frames[j].compounds.end(); ++it)
		{
			int idx = cIndex[*it];
			cCount[idx]++;
		}
	}
	
	int k=0;
	
	for (int j=0;j<cidx;j++)
	{
		if ((cCount[j] == 0)||(cCount[j] == maxval))
		{
			cMap.push_back(-1);
			vCols--;
		}
		else
		{
			cMap.push_back(k);
			k++;
		}
	}
	*/
	mat data(maxval, vCols, fill::zeros);
	
	for (int j=0;j<maxval;j++)
	{
		for (unordered_set<string>::iterator it = endStates.frames[j].compounds.begin(); it != endStates.frames[j].compounds.end(); ++it)
		{
			int idx = cIndex[*it]; // Use cMap to reduce the matrix
	
			if (idx >= 0)
			{
				data(j,idx) = 1;
			}
		}
	}
	
	vec eigvals;
	mat princomps;
	mat scores;
	
	int nVals = vCols;
	if (maxval < nVals) nVals = maxval;
	
	printf("Performing PCA of %d x %d matrix\n", data.n_rows, data.n_cols);
	
	doPCA(data, eigvals, princomps, scores, nVals);

	string outputSubDirectory = strParams["OUTPUT_DIR_PCA"];
		
	if (!directoryExists(outputSubDirectory))
		makeDirectory(outputSubDirectory);
		
	string filename = outputSubDirectory + "/eigenvalues.txt";
	
	FILE *f=fopen(filename.c_str(),"wb");
	
	for (int i=0;i<eigvals.n_rows;i++)
		fprintf(f,"%.9g\n", eigvals(i));
	fclose(f);

	double threshold = eigvals(0), weighted = 0, cutoff[4];
	double mu = 0;
	int count = 0;
	
	memset(cutoff,0,sizeof(double)*4);
	
	for (int i=0;i<eigvals.n_rows;i++)
	{
		if (eigvals(i) > threshold*0.5) cutoff[0]++;
		if (eigvals(i) > threshold*0.25) cutoff[1]++;
		if (eigvals(i) > threshold*0.1) cutoff[2]++;
		if (eigvals(i) > threshold*0.05) cutoff[3]++;
	}

	filename = outputSubDirectory + "/heredity_measure.csv";	
	
	f=fopen(filename.c_str(),"rb");
	
	if (f == NULL)
	{
		f=fopen(filename.c_str(), "wb");
		fprintf(f,"Number of Simulations, Cliff-0.5, Cliff-0.25, Cliff-0.1, Cliff-0.05\n");
		fclose(f);
	}
	else fclose(f);
	
	f=fopen(filename.c_str(),"a");
	fprintf(f,"%d, %.9g, %.9g, %.9g, %.9g\n", simidx, cutoff[0]+1, cutoff[1]+1, cutoff[2]+1, cutoff[3]+1);
	fclose(f);

	filename = outputSubDirectory + "/pca_scores.txt";	
	
	f=fopen(filename.c_str(),"wb");
	for (int i=0;i<scores.n_rows;i++)
	{
		for (int j=0;j<scores.n_cols;j++)
		{
			fprintf(f,"%.9g ", scores(i,j));
		}
		fprintf(f,"\n");
	}
	fclose(f);

	filename = outputSubDirectory + "/rawdata.txt";
	f = fopen(filename.c_str(),"wb");
	for (int i=0;i<data.n_rows;i++)
	{
		for (int j=0;j<data.n_cols;j++)
		{
			fprintf(f,"%.9g ", data(i,j));
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

void AnalysisHeredity::FeatureEliminationHeredity()
{
	int endFeatures;
	double runFeatures = 0, runFeatures_max = 0;	
	double threshold = numParams["FEATURE_ELIMINATION_THRESHOLD"];
	
	for (int i=0;i<runs.size();i++)
	{
		double tmp = runs[i].getFeatureCount(cidx, threshold, cIndex);
		runFeatures += tmp;
		if (tmp > runFeatures_max) runFeatures_max = tmp;
	}
	
	runFeatures /= (double)runs.size();
	
	endFeatures = endStates.getFeatureCount(cidx, threshold, cIndex);
	
	string filename = strParams["FEATURE_ELIMINATION_OUTPUT"];	
	
	FILE *f = fopen(filename.c_str(), "rb");
	
	if (f == NULL)
	{
		f=fopen(filename.c_str(), "wb");
		fprintf(f,"Runs, runFeatures_avg, runFeatures_max, finalFeatures, difference\n");
		fclose(f);
	} else fclose(f);
	
	f = fopen(filename.c_str(), "a");
	fprintf(f,"%d, %.6g, %.6g, %.6g\n", simidx, runFeatures, runFeatures_max, (double)endFeatures, (double)(endFeatures - runFeatures));
	fclose(f);
}

void AnalysisHeredity::EntropyHeredity()
{
	double Hbar = 0, Hcross = 0;
	
	for (int i=0;i<runs.size();i++)
		Hbar += runs[i].getEntropy();
		
	Hbar /= (double)runs.size();
	
	Hcross = endStates.getEntropy();
	
	string filename = strParams["ENTROPY_OUTPUT"];
	
	FILE *f = fopen(filename.c_str(), "a");
	fprintf(f,"%d, %.6g, %.6g, %.6g\n", simidx, Hbar, Hcross, Hcross-Hbar);
	fclose(f);
}

void AnalysisHeredity::onSimulationEnd(SimulationRequest *SR)
{
	SimulationTimeDependent *S = (SimulationTimeDependent*)SR;
	Simulation *Sim = &S->System;
	double threshold = numParams["DETECTION_THRESHOLD"];
	
	// Regions are all lumped together - lets do this differently in the future, or give multiple options
	HereditySnapshot HS;
	
	for (int i=0;i<Sim->regions.size();++i)
	{
		Region *R = &Sim->regions[i];		
		unordered_map<string, Node*>::iterator it;
				
		for (it=R->population.accessHash.begin();it!=R->population.accessHash.end();++it)
		{
			if (it->second->weight >= threshold)
			{
				if (!cIndex.count(it->first))
				{
					cIndex[it->first] = cidx++;
				}
								
				HS.compounds.insert(it->first);
			}			
		}
	}
	
	endStates.frames.push_back(HS);
	
	int oPeriod = numParams["OUTPUT_PERIOD"];
	
	if (endStates.frames.size() % oPeriod != 0) return;
	// Lets do output up to this point now

	if (strParams["PCA_ANALYSIS"] != "false")
	{
		PCAHeredity();
	}
	
	if (strParams["FEATURE_ELIMINATION_ANALYSIS"] != "false")
	{
		FeatureEliminationHeredity();
	}
	
	if (strParams["ENTROPY_ANALYSIS"] != "false")
	{
		EntropyHeredity();
	}
	
	simidx++;
}

AnalysisHeredity::AnalysisHeredity()
{	
	simType = "TimeDependent";
	analysisType = "Heredity";
	iterateCallback = 1;
	endCallback = 1;
	beginCallback = 1;
	
	numParams["OUTPUT_PERIOD"] = 1;
	strParams["RANDOM_KNOCKOUTS"] = "false";
	numParams["KNOCKOUT_COUNT"] = 10;
	numParams["KNOCKOUT_MEAN_WILDCARD_LENGTH"] = 1;
	numParams["DETECTION_THRESHOLD"] = 0.5;
	strParams["KNOCKOUT_RULE"] = "";
	strParams["PCA_ANALYSIS"] = "false";
	strParams["FEATURE_ELIMINATION_ANALYSIS"] = "false";
	strParams["FEATURE_ELIMINATION_OUTPUT"] = "heredity_elimination.txt";
	numParams["FEATURE_ELIMINATION_THRESHOLD"] = 0.9;
	strParams["OUTPUT_DIR_PCA"] = "heredity_PCA";
	strParams["ENTROPY_ANALYSIS"] = "false";
	strParams["ENTROPY_OUTPUT"] = "heredity_entropy.txt";
	
	cidx = 0;
	simidx = 0;
}

AnalysisRequest *AnalysisHeredity::clone()
{
	AnalysisHeredity *T = new AnalysisHeredity;
	
	*T = *this;
	
	return T;
}
