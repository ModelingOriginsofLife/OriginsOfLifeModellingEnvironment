#include "../includes.h"

/* SimulationTimeDependent */
SimulationTimeDependent::SimulationTimeDependent()
{
	simType = "TimeDependent";
	numParams["SUBITER"] = 10000;
	numParams["MAXITER"] = 1000;
	numParams["REPEAT"] = 1;
	strParams["FIXED_CONCENTRATION"] = "false";
	strParams["NO_REJECTIONS"] = "false";
}

bool SimulationTimeDependent::Iterate(ChemistryComputation &C)
{
	int subIters = numParams["SUBITER"];
	int maxiter = numParams["MAXITER"];
	IterationParams I;
		
	if (strParams["FIXED_CONCENTRATION"] == "false")
		I.adjustConcentrations = true;
	else I.adjustConcentrations = false;
	
	if (strParams["NO_REJECTIONS"] == "false")
		I.noRejections = false;
	else I.noRejections = true;

	for (int i=0;i<subIters;i++)
		System.Iterate(I);
		
	iter++;
	
	if (iter>=maxiter) return true;
	return false;
}

void SimulationTimeDependent::setupSimulation(ChemistryComputation &C)
{
	System = C.simTemplate;
	System.connectToChemistry(&C);
}

SimulationRequest *SimulationTimeDependent::clone()
{
	SimulationTimeDependent *T = new SimulationTimeDependent;
	
	*T = *this;
	
	return T;
}

void SimulationTimeDependent::executeSimulation(ChemistryComputation &C)
{
	printf("Overload succeeded!\n");
	for (int i=0;i<numParams["REPEAT"];i++)
		doSimulation(C);
}
