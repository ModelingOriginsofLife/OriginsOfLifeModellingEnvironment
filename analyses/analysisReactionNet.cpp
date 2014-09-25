#include "../includes.h"

void parseReactionSubstring(string str, vector<string> &array)
{
	int pos;
	
	do
	{
		pos = str.find_first_of("+");
		
		if (pos != string::npos)
		{
			string sstr = str.substr(0,pos);
			if (sstr.length())
				array.push_back(sstr);
			str = str.substr(pos+1);
		}
		else
		{
			if (str.length())
				array.push_back(str);
		}
	} while (pos != string::npos);
}

void parseReactionString(string str, vector<string> &reac, vector<string> &prod)
{
	string tok, rstr, pstr;			

	reac.clear(); prod.clear();

	str.erase (std::remove (str.begin(), str.end(), ' '), str.end());
			
	int pos = str.find_first_of("=");

	rstr = str.substr(0,pos);
	pstr = str.substr(pos+1);
	
	parseReactionSubstring(rstr, reac);
	parseReactionSubstring(pstr, prod);					
}

void AnalysisReactionNet::onSimulationEnd(SimulationRequest *SR)
{
	unordered_map<string,long int>::iterator it;
	FILE *f = fopen(strParams["OUTPUT_FILENAME"].c_str(), "wb");
	FILE *f2;
	int j;
	int doGraphviz = 0;
	
	if (strParams["GRAPHVIZ"]!="false")
	{
		doGraphviz = 1;
		f2 = fopen(strParams["OUTPUT_GRAPHVIZ_FILENAME"].c_str(), "wb");
		fprintf(f2,"digraph G {\n");
	}
	
	int ridx = 0;
	unordered_map<string, int> compounds;
	
	for (it=reactionNetwork.begin();it!=reactionNetwork.end();++it)
	{
		fprintf(f,"%s, %ld\n",it->first.c_str(), it->second);
		
		if (doGraphviz)
		{
			vector<string> reac, prod;
			string str = it->first;
			
			parseReactionString(str,reac,prod);
			
			fprintf(f2,"\trNode%d [shape=point];\n", ridx);
			
			for (j=0;j<reac.size();j++)
			{				
				if (!compounds.count(reac[j]))
				{
					fprintf(f2,"\t\"%s\";\n", reac[j].c_str());
					compounds[reac[j]]=1;
				}
				fprintf(f2,"\t\"%s\" -> rNode%d;\n", reac[j].c_str(), ridx);
			}
			
			for (j=0;j<prod.size();j++)
			{
				if (!compounds.count(prod[j]))
				{
					fprintf(f2,"\t\"%s\";\n", prod[j].c_str());
					compounds[prod[j]]=1;
				}
				fprintf(f2,"\trNode%d -> \"%s\";\n", ridx, prod[j].c_str());
			}
			
			ridx++;
		}
	}
	
	if (doGraphviz) 
	{
		fprintf(f2,"}\n");
		fclose(f2);
	}
	fclose(f); 
}

void AnalysisReactionNet::onReaction(vector<string> &reac, vector<string> &prod)
{
	string rxnString = "";
	
	for (int i=0;i<reac.size();i++)
	{
		rxnString = rxnString + reac[i]; 
		if (i<reac.size()-1)
			rxnString = rxnString + " + ";
	}

	rxnString = rxnString + " = ";
	
	for (int i=0;i<prod.size();i++)
	{
		rxnString = rxnString + prod[i]; 
		if (i<prod.size()-1)
			rxnString = rxnString + " + ";
	}
	
	if (reactionNetwork.count(rxnString))
		reactionNetwork[rxnString]++;
	else reactionNetwork[rxnString] = 1;
}

AnalysisReactionNet::AnalysisReactionNet()
{	
	simType = "TimeDependent";
	analysisType = "ReactionNet";
	iterateCallback = 0;
	reactionCallback = 1;
	endCallback = 1;

    strParams["GRAPHVIZ"] = "false";
	strParams["OUTPUT_GRAPHVIZ_FILENAME"] = "reaction_network.dot";
	strParams["OUTPUT_FILENAME"] = "reaction_network.txt";
}

AnalysisRequest *AnalysisReactionNet::clone()
{
	AnalysisReactionNet *T = new AnalysisReactionNet;
	
	*T = *this;
	
	return T;
}
