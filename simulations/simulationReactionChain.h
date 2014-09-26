/* Do a simulation which follows the transformation of a single 'backbone' compound */
class SimulationReactionChain : public SimulationRequest
{
	public:		
		void setupSimulation(ChemistryComputation &C);
		bool Iterate(ChemistryComputation &C);
		SimulationRequest *clone();
		
		SimulationReactionChain();

		int iter;
		Simulation System;
		
		unordered_map<string, long int> normalizationMatrix;
		unordered_map<string, long int> transitionMatrix;
};
