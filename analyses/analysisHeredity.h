class HereditySnapshot
{
	public:
		unordered_set<string> compounds;
};

class HeredityRunData
{
	public:
		vector<HereditySnapshot> frames;
		int getFeatureCount(int vCols, double threshold, unordered_map<string, int> &cIndex, string outfile, vector<char> &features);
		
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
		void FeatureEliminationHeredity();
		void EntropyHeredity();
		bool checkAgainstKnockouts(string compound, Library &L);
		string AnalysisHeredity::getSubdirectory(string dir);
	
		AnalysisRequest *clone();
		AnalysisHeredity();		
		
	private:
		int simidx; // This keeps track of the simulation index in a sheaf of multiple simulations
		int cidx; // This keeps track of novel compounds
		string subdirectory;
		
		ReactionRule RRule;
		unordered_set<string> knockoutList;
		
		unordered_map<string, int> cIndex;
		unordered_map<int, string> rIndex;
		vector<HeredityRunData> runs;
		
		HeredityRunData endStates;
};
