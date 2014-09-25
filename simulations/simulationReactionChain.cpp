#include "../includes.h"

/* SimulationReactionChain
 * 
 *	This type of simulation takes a single chemical compound and progressively reacts it with things in the bath.
 *  The result is effectively an exploration of 'first order chemistry', which represents the case where there is
 *  a large separation in concentration between compounds in the bath and compounds in some newly formed autocatalytic
 *  set.
 * 
 *  Reaction chain simulations don't generally play well with multiple regions because they do not collect separate per-region statistics
 *  However, multiple regions can be used to probe the properties of multiple separate compounds in the same simulation
 */

SimulationReactionChain::SimulationReactionChain()
{
	simType = "ReactionChain";

	numParams["DEPTH"] = 20; // How many transformations deep should the simulation explore starting from the initial state?
	numParams["MAXITER"] = 10000; // How many different runs should be executed?

	strParams["RANDOM_STARTING_COMPOUND"] = "false"; // Use a random starting compound to initialize each run
	
	// If using a random starting compound, generate it using these characters and wildcards
	// Empty rules are ignored
	numParams["MEAN_WILDCARD_LENGTH"] = 1; // This is the mean of a Poisson distribution of lengths if a variable-length wildcard is used in these rules
	strParams["STARTING_COMPOUND_RULE"] = ""; 
}

string generateRandomCompound(string rulestr, ChemistryComputation &C, double meanLength)
{
	ReactionRule R;
	unordered_map<char, string> bound;
	double continueChance = meanLength/(1.0+meanLength);
	string compound = "";
	
	R.rule = rulestr;
	R.parseRule(C.L);
	
	for (int i=0;i<R.reacRules[0].size();i++)
	{
		Symbol S = R.reacRules[0][i];
		
		if (S.isWild)
		{
			if (S.isBound)
			{
				if (bound.count(S.label))
				{
					if (S.isReversed)
						compound = compound + reverseString(bound[S.label]);
					else
						compound = compound + bound[S.label];
				}
				else
				{
					string substr = "";
					
					while (prand(continueChance))
					{
						substr = substr + C.L.getMemberOfSet(S.label, S.doConjugate);
					} 
					
					if (S.isReversed) substr = reverseString(substr);
					
					compound = compound + substr;
					bound[S.label] = substr;
				}
			}
			else
			{
				while (prand(continueChance))
				{
					compound = compound + C.L.getMemberOfSet(S.label, S.doConjugate);
				}
			}
		}
		else
		{
			if (S.isBound)
			{
				if (bound.count(S.label))
				{
					compound = compound + bound[S.label];
				}
				else
				{
					string substr = "";
					
					substr = substr + C.L.getMemberOfSet(S.label, S.doConjugate);
					compound = compound + substr;
					
					bound[S.label] = substr;
				}
			}
			else
			{
				compound = compound + C.L.getMemberOfSet(S.label, S.doConjugate);
			}
		}
	}
		
	return compound;
}

bool SimulationReactionChain::Iterate(ChemistryComputation &C)
{
	int depth = numParams["DEPTH"];
	int maxiter = numParams["MAXITER"];
	IterationParams I;	
	int i, j;
	
	I.bathReactions = false;
	
	System = C.simTemplate; // Reset the system each iteration for this kind of simulation
	
	for (i=0;i<System.regions.size();i++)
	{
		if (strParams["RANDOM_STARTING_COMPOUND"] != "false")
		{
			string compound;
			
			// Generate a random compound according to the rules
			compound = generateRandomCompound(strParams["STARTING_COMPOUND_RULE"], C, numParams["MEAN_WILDCARD_LENGTH"]);
			
			System.regions[i].addCompound(compound, 1);
		}
	}
	
	for (int j=0;j<depth;j++)
	{
		System.Iterate(&C, I);
		for (i=0;i<System.regions.size();i++)
		{			
			// Restrict the region back down to a single compound
			if ((System.regions[i].population.Root != NULL)&&(System.regions[i].population.Root->weight > 0))
			{
				string compound = System.regions[i].pickRandomCompound(WEIGHT_HEAVY, 0);
				
				Tree newTree;
				System.regions[i].population = newTree;
				System.regions[i].addCompound(compound, 1);
			}
		}
	}
	
	iter++;
	
	if (iter>=maxiter) return true;
	return false;
}

void SimulationReactionChain::setupSimulation(ChemistryComputation &C)
{
	System = C.simTemplate;
}

SimulationRequest *SimulationReactionChain::clone()
{
	SimulationReactionChain *T = new SimulationReactionChain;
	
	*T = *this;
	
	return T;
}
