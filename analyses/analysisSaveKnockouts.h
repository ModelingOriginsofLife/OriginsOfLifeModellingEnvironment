class AnalysisSaveKnockouts : public AnalysisRequest
{
	public:		
		void onIteration(SimulationRequest *SR) {};
		void onSimulationBegin(SimulationRequest *SR);
		void onSimulationEnd(SimulationRequest *SR) {};
		void onReaction(vector<string> &reac, vector<string> &prod, Region *Reg) {}; // Pointers make everything faster...
	
		AnalysisRequest *clone();
		AnalysisSaveKnockouts();
};

