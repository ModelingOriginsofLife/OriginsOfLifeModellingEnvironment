class AnalysisHeredity: public AnalysisRequest
{
	public:
		void onIteration(SimulationRequest *SR);
		void onSimulationEnd(SimulationRequest *SR) {};
		void onReaction(vector<string> &reac, vector<string> &prod) {}; // Pointers make everything faster...
	
		AnalysisRequest *clone();
		AnalysisHeredity();		
};
