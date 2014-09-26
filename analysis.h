
/* AnalysisRequest
 * 
 * This class references a particular sort of analysis method which the user has requested be applied to the simulation(s)
 * 
 * The main way these classes function is by callbacks. Each class defines onIteration, onSimulationEnd, onReaction, etc 
 * functions which contribute to the analysis. References are passed to the simulation state as well as any detailed
 * information of relevance to the particular kind of callback.
 * 
 * For efficiency, static values are defined to determine if a given callback is needed, so that we don't call each thing on each type
 * of analysis.
 */

class AnalysisRequest
{
	public:
		string simType;
		string analysisType;
	
		unordered_map<string, double> numParams;
		unordered_map<string, string> strParams;
		
		int iterCounter; // Used to apply onIteration based on PERIOD numParam
	
		virtual void onIteration(SimulationRequest *SR)=0;
		virtual void onSimulationBegin(SimulationRequest *SR)=0;
		virtual void onSimulationEnd(SimulationRequest *SR)=0;
		virtual void onReaction(vector<string> &reac, vector<string> &prod)=0; // Pointers make everything faster...
	
		int iterateCallback;
		int beginCallback;
		int endCallback;
		int reactionCallback;
		
		AnalysisRequest();
		virtual AnalysisRequest *clone()=0;
		void parse(ifstream &file);
	private:
		
};

extern vector<AnalysisRequest*> registeredAnalyses;
extern void registerAnalyses();
