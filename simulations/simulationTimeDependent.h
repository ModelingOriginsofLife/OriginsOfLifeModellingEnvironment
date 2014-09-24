/* Do a full time-dependent simulation */
class SimulationTimeDependent : public SimulationRequest
{
	public:		
		void setupSimulation(ChemistryComputation &C);
		bool Iterate(ChemistryComputation &C);
		SimulationRequest *clone();
		
		SimulationTimeDependent();

		int iter;
		Simulation System;
};
