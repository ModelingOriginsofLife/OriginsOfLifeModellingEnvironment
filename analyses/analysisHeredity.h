class HereditySnapshot
{
	public:
		unordered_set<string> compounds;
};

class HeredityRunData
{
	public:
		vector<HereditySnapshot> frames;
		
		double getEntropy();
};

class AnalysisHeredity: public AnalysisRequest
{
	public:
		void onIteration(SimulationRequest *SR);
		void onSimulationBegin(SimulationRequest *SR);
		void onSimulationEnd(SimulationRequest *SR);
		void onReaction(vector<string> &reac, vector<string> &prod, Region *Reg) {}; // Pointers make everything faster...
		void PCAHeredity();
		void EntropyHeredity();
	
		AnalysisRequest *clone();
		AnalysisHeredity();		
		
	private:
		int simidx; // This keeps track of the simulation index in a sheaf of multiple simulations
		int cidx; // This keeps track of novel compounds
		
		unordered_map<string, int> cIndex;
		vector<HeredityRunData> runs;
		
		HeredityRunData endStates;
};
