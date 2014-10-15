#include "includes.h"

bool pairComparator(pair<string,int> A, pair<string,int> B) { return A.second > B.second; }

void Region::writePopulationText(char *fname)
{
	FILE *f=fopen(fname,"wb");
	unordered_map<string, Node*>::iterator it;
	vector< pair<string,int> > cList;
	
	for (it=population.accessHash.begin();it!=population.accessHash.end();it++)
	{
		pair<string,int> P;
		
		P.first = it->first; 
		P.second = floor(it->second->weight);
		cList.push_back(P);
	}
	
	sort(cList.begin(),cList.end(), pairComparator);
	
	for (int i=0;i<cList.size();i++)
		fprintf(f,"%s %d\n",cList[i].first.c_str(),cList[i].second);
		
	fclose(f);
	cList.clear();
}

void Simulation::writePopulationText(char *directory)
{
	char Str[4096];
	unordered_map<string,int>::iterator it;
	
	for (it=regionMap.begin();it!=regionMap.end();++it)
	{
		snprintf(Str,4096,"%s/%s.txt",directory,it->first.c_str());
		
		regions[it->second].writePopulationText(Str);
	}
}

bool directoryExists(string dir)
{
	using namespace boost::filesystem;
	return exists(dir);
}

void makeDirectory(string dir)
{
	using namespace boost::filesystem;
	create_directory(dir);
}

string SimulationRequest::getSubdirectory(string dir)
{
	string curheader = parentChem->getSubdirectory(headerDirectory);
	
	if (!directoryExists(curheader))
		makeDirectory(curheader);		
	
	return curheader + dir;
}

string ChemistryComputation::getSubdirectory(string dir)
{
	string totaldir = headerDirectory + dir;
	
	if (!directoryExists(headerDirectory))
		makeDirectory(headerDirectory);		
	
	return totaldir;
}
