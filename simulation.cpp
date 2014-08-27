#include "includes.h"

void Region::addCompound(string str, int count)
{
	Node *leaf;
	
	if (bath.accessHash.count(str)) return; // Can't change the concentration of the bath
	
	if (population.accessHash.count(str))
	{
		leaf = population.accessHash[str];
		leaf->modifyWeight(leaf->weight + count);
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
		
		if (curPop<=count)
		{
			population.deleteLeaf(leaf);
		}
		else leaf->modifyWeight(curPop-count);
	}
}

string Region::pickRandomCompound()
{
	if (population.Root == NULL) { if (bath.Root == NULL) return ""; else return bath.Root->findRandomLeaf(false)->value; }
	else if (bath.Root == NULL) return population.Root->findRandomLeaf(false)->value;
	
	double total = population.Root->weight + bath.Root->weight;
	
	if (prand(population.Root->weight / total)) return population.Root->findRandomLeaf(false)->value;
	else return bath.Root->findRandomLeaf(false)->value;
}

void Region::doRandomSinglet(ChemistryComputation *C)
{	
	vector<string> rList; rList.push_back(pickRandomCompound());
	
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

	rList.push_back(pickRandomCompound());
	rList.push_back(pickRandomCompound());
	
	if (rList[0] == rList[1])
	{
		if (getConcentration(rList[0])<2) return; // Not enough reactant to make this go!
	}
	
	Outcome P = C->getReactionProducts(rList);
	
	if (P.reacted)
	{
		if (prand(P.rate))
		{
			for (int i=0;i<P.products.size();i++)
			{
				addCompound(P.products[i],1);
			}
			
			for (int i=0;i<rList.size();i++)
				removeCompound(rList[i],1);
		}
	}
}

void Simulation::Iterate(ChemistryComputation *C)
{
	int i;
	
	for (i=0;i<regions.size();i++)
	{
		int N=0;
		
		if ((regions[i].population.Root != NULL)||(regions[i].bath.Root != NULL))
		{
			regions[i].doRandomSinglet(C);
			regions[i].doRandomDoublet(C);
		}
	}
}
