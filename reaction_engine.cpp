#include "includes.h"

string ChemistryComputation::getReactionString(vector<string> &reactants)
{
	string S = "";
	
	for (int i=0;i<reactants.size();i++)
		S = S + reactants[i]+" ";
		
	return S;
}

Outcome ReactionList::randomOutcome()
{
	Outcome O;
	
	if (rList.size())
	{
		int idx = irand(rList.size());	
		O.reacted = true;
		O.products = rList[idx].products; // Currently this ignores reaction rates in sampling the list, so it forces an accept/reject step at the terminal point
		O.rate = rList[idx].weight;
	}
	else O.reacted = false;
	
	return O;
}

Outcome ChemistryComputation::getReactionProducts(vector<string> &reactants)
{
	string reactionString = getReactionString(reactants);
	ReactionList R;
	
	if (reactionHash.count(reactionString))
	{
		return reactionHash[reactionString].randomOutcome();
	}

	for (int i=0;i<reactions.size();i++)
	{
		if (reactants.size() == reactions[i].Nreac)
		{
			// This is a possible result...
			
			vector<vector<string> > S = reactions[i].getAllProducts(reactants, L);
			
			for (int j=0;j<S.size();j++)
			{
				Reaction reac;
			
				reac.weight = reactions[i].k;
				reac.typeidx = i;
				reac.products = S[j];
				
				R.rList.push_back(reac);
			}
		}
	}
	
/*	for (int i=0;i<reactants.size();i++)
	{
		printf("%s%s", i>0 ? " + " : "", reactants[i].c_str());
	}
	
	printf(" = ");
	
	for (int i=0;i<R.rList.size();i++)
	{
		printf("%s", i>0 ? ", " : "");
		for (int j=0;j<R.rList[i].products.size();j++)
		{
			printf("%s%s", j>0 ? " + " : "", R.rList[i].products[j].c_str());
		}
	}
	
	printf("\n");*/
	reactionHash[reactionString] = R;
	
	return R.randomOutcome();
}

int main(int argc, char **argv)
{
	registerSimulations();
	registerAnalyses();
	
	ifstream file(argv[1]);	
	ChemistryComputation C = parseConfigFile(file);
	file.close();
	
	int jobid;
	
	for (jobid = 0; jobid < C.jobs.size(); jobid++)
	{
		C.jobs[jobid]->executeSimulation(C);
	}
	
	/*	printf("%d: %.3g %.3g\n",iter, C.simTemplate.regions[0].getTotalPopulation(), C.simTemplate.regions[0].getTotalLength());
		
		for (subiter=0;subiter<10000;subiter++)
			C.simTemplate.Iterate(&C);
	*/
}
