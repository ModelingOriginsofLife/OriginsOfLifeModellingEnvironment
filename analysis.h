
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
		static const string simType;
		static const string analysisType;
	
		int iterFrequency;
		int iterCounter; // Used to apply onIteration based on iterFrequency
	
		void onIteration(SimulationRequest *S);
		void onSimulationEnd(SimulationRequest *S);
		void onReaction(vector<string> &reac, vector<string> &prod); // Pointers make everything faster...
	
		static const int iterateCallback = 0;
		static const int endCallback = 0;
		static const int reactionCallback = 0;
		
	private:
		
};

class AnalysisDistribution : public AnalysisRequest
{
	public:
		static const string simType;
		static const string analysisType;
		
		string outputSubDirectory;
		
		void onIteration(SimulationTimeDependent *S);
	
		static const int iterateCallback = 1;
};

