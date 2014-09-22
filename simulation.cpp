#include "includes.h"

vector <SimulationRequest*> registeredSimulations;

void registerSimulations()
{
	registeredSimulations.push_back(new SimulationTimeDependent);
	registeredSimulations.push_back(new SimulationExplore);
}

/* SimulationTimeDependent */
SimulationTimeDependent::SimulationTimeDependent()
{
	simType = "TimeDependent";
	numParams["SUBITER"] = 10000;
	numParams["MAXITER"] = 1000;
}

bool SimulationTimeDependent::Iterate(ChemistryComputation &C)
{
	int subIters = numParams["SUBITER"];
	int maxiter = numParams["MAXITER"];
	
	for (int i=0;i<subIters;i++)
		System.Iterate(&C);
		
	iter++;
	
	if (iter>=maxiter) return true;
	return false;
}

void SimulationTimeDependent::setupSimulation(ChemistryComputation &C)
{
	System = C.simTemplate;
}

SimulationRequest *SimulationTimeDependent::clone()
{
	SimulationTimeDependent *T = new SimulationTimeDependent;
	
	*T = *this;
	
	return T;
}

/* SimulationExplore */
SimulationExplore::SimulationExplore()
{
	simType = "Explore";
	numParams["MAXITER"] = 1000;
}

bool SimulationExplore::Iterate(ChemistryComputation &C)
{
	return true;
}

void SimulationExplore::setupSimulation(ChemistryComputation &C)
{
}

SimulationRequest *SimulationExplore::clone()
{
	SimulationExplore *T = new SimulationExplore;
	
	*T = *this;
	
	return T;
}

/* SimulationRequest */
/*SimulationRequest *SimulationRequest::clone()
{
	SimulationRequest *T = new SimulationRequest;
	
	*T = *this;
	
	printf("I shouldn't be running this.\n");
	return T;
}

bool SimulationRequest::Iterate(ChemistryComputation &C)
{
	printf("This code should never run.\n");
	return true;
}

void SimulationRequest::setupSimulation(ChemistryComputation &C)
{
}
*/
void SimulationRequest::doSimulation(ChemistryComputation &C)
{
	C.curSimType = simType;
	
	setupSimulation(C);
	
	while (!Iterate(C))
	{
		// do any analyses requested
		
		for (int i=0;i<C.analyses.size();i++)
		{
			if ((C.analyses[i]->iterateCallback)&&(C.analyses[i]->simType == simType))
			{
				C.analyses[i]->iterCounter++;
				
				if (C.analyses[i]->iterCounter >= C.analyses[i]->numParams["PERIOD"])
				{
					C.analyses[i]->onIteration(this);
					C.analyses[i]->iterCounter = 0;
				}
			}
		}
	}
	
	for (int i=0;i<C.analyses.size();i++)
	{
		if ((C.analyses[i]->endCallback)&&(C.analyses[i]->simType == simType))
		{
			C.analyses[i]->onSimulationEnd(this);
		}
	}
}

void SimulationRequest::parse(ifstream &file)
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

/* Region */

Region::Region()
{
}

double Region::getTotalPopulation()
{
	double total = 0;
	
	if (population.Root != NULL) total += population.Root->weight;
	if (bath.Root != NULL) total += bath.Root->weight;
	
	return total;
}

double Region::getTotalLength()
{
	double total = 0;
	
	if (population.Root != NULL) 
		total += population.Root->lweight;
	if (bath.Root != NULL) 
		total += bath.Root->lweight;
	
	return total;
}

void Region::addCompound(string str, int count)
{
	Node *leaf;
	
	if (bath.accessHash.count(str)) { return; } // Can't change the concentration of the bath
	
	if (population.accessHash.count(str))
	{
		leaf = population.accessHash[str];
		
		leaf->modifyWeight(leaf->weight  + count, 
						   leaf->lweight + count * str.length());
	}
	else
	{
		population.addNewLeaf(str, count);
	}
}

void Region::removeCompound(string str, int count)
{
	Node *leaf;
	
	if (bath.accessHash.count(str)) return; // Can't change the concentration of the bath

	if (population.accessHash.count(str))
	{
		leaf = population.accessHash[str];
		double curPop = leaf->weight;
		
		if (curPop-count<=1e-3)
		{
			population.deleteLeaf(leaf);
		}
		else 
		{
			leaf->modifyWeight(curPop-count, leaf->lweight - count * str.length());
		}
	}
}

string Region::pickRandomCompound(WeightingType wType)
{
	if (population.Root == NULL) 
	{ 
		if (bath.Root == NULL) return ""; 
		else return bath.Root->findRandomLeaf(wType)->value; 
	}
	else if (bath.Root == NULL) 
		return population.Root->findRandomLeaf(wType)->value;
	
	double total;
	
	if (wType == WEIGHT_HEAVY) total = population.Root->weight + bath.Root->weight;
	if (wType == WEIGHT_HEAVYLENGTH) total = population.Root->lweight + bath.Root->lweight;
	
	if (prand(population.Root->weight / total)) return population.Root->findRandomLeaf(wType)->value;
	else return bath.Root->findRandomLeaf(wType)->value;
}

void Region::doRandomSinglet(ChemistryComputation *C)
{	
	vector<string> rList; 
	
	rList.push_back(pickRandomCompound(C->singletRate));
		
	Outcome P = C->getReactionProducts(rList);
	
	if (P.reacted)
	{
		if (prand(P.rate))
		{
			for (int i=0;i<P.products.size();i++)
			{
				addCompound(P.products[i],1);
			}
			removeCompound(rList[0],1);
		}
	}
}

double Region::getConcentration(string compound)
{
	if (bath.accessHash.count(compound)) return bath.accessHash[compound]->weight;
	if (population.accessHash.count(compound)) return population.accessHash[compound]->weight;
	
	return 0;
}

void Region::doRandomDoublet(ChemistryComputation *C)
{	
	vector<string> rList; 

	rList.push_back(pickRandomCompound(WEIGHT_HEAVY));
	rList.push_back(pickRandomCompound(WEIGHT_HEAVY));
	
	if (rList[0] == rList[1])
	{
		if (getConcentration(rList[0])<2-1e-3) return; // Not enough reactant to make this go!
	}
	
	Outcome P = C->getReactionProducts(rList);
	
	if (P.reacted)
	{
		if (prand(P.rate))
		{
			/* If there are lots of analyses, this will be slow and we 
			 * should pre-generate the list of applicable ones. However,
			 * it is unclear whether that is actually the usual case so 
			 * lets do some profiling before taking that step */
			 
			for (int i=0;i<C->analyses.size();i++)
			{
				if ((C->analyses[i]->reactionCallback)&&(C->analyses[i]->simType == C->curSimType)) // This is a bit of a hack, but otherwise we have to thread in info about the simulationrequest to here
				{
					C->analyses[i]->onReaction(rList, P.products);
				}
			}
			
			for (int i=0;i<P.products.size();i++)
			{
				addCompound(P.products[i],1);
			}
			
			for (int i=0;i<rList.size();i++)
			{
				removeCompound(rList[i],1);
			}
		}
	}	
}

void Simulation::Iterate(ChemistryComputation *C)
{
	int i;
	
	for (i=0;i<regions.size();i++)
	{
		double Pdouble, Psingle, P=0;
		
		if ((regions[i].population.Root != NULL)||(regions[i].bath.Root != NULL))
		{
			if (C->singletRate == WEIGHT_HEAVY)
				Psingle = regions[i].getTotalPopulation();
			else if (C->singletRate == WEIGHT_HEAVYLENGTH)
				Psingle = regions[i].getTotalLength();
				
			if (C->doubletRate == WEIGHT_HEAVY)
				Pdouble = regions[i].getTotalPopulation();
			else if (C->doubletRate == WEIGHT_HEAVYLENGTH)
				Pdouble = regions[i].getTotalLength();
				
			P=Psingle+Pdouble+1e-8;
			Psingle/=P; Pdouble/=P;
		
			if (prand(Psingle))
				regions[i].doRandomSinglet(C);
			else
				regions[i].doRandomDoublet(C);
		}
	}
}

ChemistryComputation::ChemistryComputation()
{
	singletRate = WEIGHT_HEAVY;
	doubletRate = WEIGHT_HEAVY;
}
