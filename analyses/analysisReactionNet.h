class AnalysisReactionNet : public AnalysisRequest
{
	public:		
		void onIteration(SimulationRequest *SR) {};
		void onSimulationBegin(SimulationRequest *SR) {};
		void onSimulationEnd(SimulationRequest *SR);
		void onReaction(vector<string> &reac, vector<string> &prod, Region *Reg); 
	
		AnalysisRequest *clone();
		AnalysisReactionNet();
	
	private:
		unordered_map<string, long int> reactionNetwork;
};

