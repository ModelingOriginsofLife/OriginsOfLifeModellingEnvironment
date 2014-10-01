#include "includes.h"

vector <SimulationRequest*> registeredSimulations;

void registerSimulations()
{
	registeredSimulations.push_back(new SimulationTimeDependent);
	registeredSimulations.push_back(new SimulationReactionChain);
}

void SimulationRequest::doSimulation(ChemistryComputation &C)
{
	C.curSimType = simType;
	
	setupSimulation(C);
	
	for (int i=0;i<C.analyses.size();i++)
	{
		if ((C.analyses[i]->beginCallback)&&(C.analyses[i]->simType == simType))
		{
			C.analyses[i]->onSimulationBegin(this);					
		}
	}
	
	int globalIterCounter = 0;
	while (!Iterate(C))
	{
		// do any analyses requested
		for (int i=0;i<C.analyses.size();i++)
		{
			if ((C.analyses[i]->iterateCallback)&&(C.analyses[i]->simType == simType))
			{
				C.analyses[i]->iterCounter++;
				
				if ((C.analyses[i]->iterCounter >= C.analyses[i]->numParams["PERIOD"]) && (C.analyses[i]->iterCounter >= C.analyses[i]->numParams["BEGIN_ITER"]))
				{
					int end = C.analyses[i]->numParams["END_ITER"];
					
					if ((end<0)||(globalIterCounter < end))
					{
						C.analyses[i]->onIteration(this);
						C.analyses[i]->iterCounter = 0;
					}
				}
			}
		}
		
		globalIterCounter++;
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

/* IterationParams */
IterationParams::IterationParams()
{
	adjustConcentrations = true;
	noRejections = false;
	bathReactions = true;
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

bool Region::compoundExists(string str)
{
	if (bath.accessHash.count(str)) { return true; }
	if (population.accessHash.count(str)) { return true; }
	
	return false;
}

void Region::addCompound(string str, int count)
{
	Node *leaf;
	
	if (str.length() == 0) return; // ignore zero-length strings
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

string Region::pickRandomCompound(WeightingType wType, bool includeBath)
{
	if (population.Root == NULL) 
	{ 
		if (!includeBath || (bath.Root == NULL)) return ""; 
		else return bath.Root->findRandomLeaf(wType)->value; 
	}
	else if (!includeBath || (bath.Root == NULL))
		return population.Root->findRandomLeaf(wType)->value;
	
	double total;
	
	if (wType == WEIGHT_HEAVY) total = population.Root->weight + bath.Root->weight;
	if (wType == WEIGHT_HEAVYLENGTH) total = population.Root->lweight + bath.Root->lweight;
	
	if (prand(population.Root->weight / total)) return population.Root->findRandomLeaf(wType)->value;
	else return bath.Root->findRandomLeaf(wType)->value;
}

void Region::doRandomSinglet(IterationParams &I)
{	
	vector<string> rList; 
	ChemistryComputation *C = parentSimulation->parentChem;
	
	rList.push_back(pickRandomCompound(C->singletRate, I.bathReactions));
	
	if (!rList[0].length()) return;	
/*	if (!I.bathReactions)
	{
		if (bath.accessHash.count(rList[0])) return; // No bath reactions permitted
	}*/
	
	Outcome P = C->getReactionProducts(rList);
	
	if (P.reacted)
	{
		if (prand(P.rate))
		{
			for (int i=0;i<P.products.size();i++)
			{
				if (I.adjustConcentrations || !compoundExists(P.products[i])) // Never increase concentrations during network exploration
					addCompound(P.products[i],1);
			}
			
			if (I.adjustConcentrations) // Never remove compounds during network exploration
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

bool Region::checkKnockout(string &compound)
{
	if (knockouts.count(compound)) return true;
	if (parentSimulation->knockouts.count(compound)) return true;
	
	return false;
}

void Region::doRandomDoublet(IterationParams &I)
{	
	vector<string> rList; 
	ChemistryComputation *C = parentSimulation->parentChem;
	
	rList.push_back(pickRandomCompound(WEIGHT_HEAVY, true));
	rList.push_back(pickRandomCompound(WEIGHT_HEAVY, I.bathReactions));
	
	if (!rList[0].length() || !rList[1].length()) return;
/*	if (!I.bathReactions)
	{
		if (bath.accessHash.count(rList[0]) && bath.accessHash.count(rList[1])) return; // No bath reactions permitted
	}*/
	
	if (I.adjustConcentrations && (rList[0] == rList[1]))
	{
		if (getConcentration(rList[0])<2-1e-3) return; // Not enough reactant to make this go!
	}
	
	Outcome P = C->getReactionProducts(rList);
	
	for (int i=0;i<rList.size();i++)
		if (checkKnockout(rList[i])) return; // Product compound is specifically excluded
	
	if (P.reacted)
	{
		if (I.noRejections || prand(P.rate))
		{
			/* If there are lots of analyses, this will be slow and we 
			 * should pre-generate the list of applicable ones. However,
			 * it is unclear whether that is actually the usual case so 
			 * lets do some profiling before taking that step */
			 
			for (int i=0;i<C->analyses.size();i++)
			{
				if ((C->analyses[i]->reactionCallback)&&(C->analyses[i]->simType == C->curSimType)) // This is a bit of a hack, but otherwise we have to thread in info about the simulationrequest to here
				{
					C->analyses[i]->onReaction(rList, P.products, this);
				}
			}
			
			for (int i=0;i<P.products.size();i++)
			{
				if (I.adjustConcentrations || !compoundExists(P.products[i])) // Never increase concentrations during network exploration
					addCompound(P.products[i],1);
			}
			
			if (I.adjustConcentrations) // Never remove compounds during network exploration
			{
				for (int i=0;i<rList.size();i++)
				{
					removeCompound(rList[i],1);
				}
			}
		}
	}	
}

void Simulation::Iterate(IterationParams &I)
{
	int i;
	ChemistryComputation *C = parentChem;
	 
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
				regions[i].doRandomSinglet(I);
			else
				regions[i].doRandomDoublet(I);
		}
	}
}

void Simulation::connectToChemistry(ChemistryComputation *C)
{
	parentChem = C;
	
	for (int i=0;i<regions.size();i++)
		regions[i].parentSimulation = this;
}

ChemistryComputation::ChemistryComputation()
{
	singletRate = WEIGHT_HEAVY;
	doubletRate = WEIGHT_HEAVY;
}
