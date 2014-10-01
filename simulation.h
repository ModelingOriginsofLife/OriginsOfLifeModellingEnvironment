class ChemistryComputation;

/* IterationParams allows various hard-coded options to be passed to the chemistry engine */

class IterationParams
{
	public:
		bool adjustConcentrations; // If this is false, only keep track of whether a compound is present or absent
		bool noRejections; // Never reject a reaction on the basis of rate
		bool bathReactions; // Allow reactions between pairs of compounds in the bath
		IterationParams(); 
};

class Simulation;

class Region
{
	public:
		Tree population;
		Tree bath; // Compounds that are assumed to always be available, and their relative concentrations
		Simulation *parentSimulation;
		
		unordered_set<string> knockouts; // This adds to any global knockouts defined for the simulation
		bool checkKnockout(string& compound);
		
		bool compoundExists(string str);
		void addCompound(string str, int count);
		void removeCompound(string str, int count);
		string pickRandomCompound(WeightingType wType, bool includeBath);
		void writePopulationText(char *FName);
		double getConcentration(string compound);
		
		void doRandomSinglet(IterationParams &I);
		void doRandomDoublet(IterationParams &I);
		
		double getTotalPopulation();
		double getTotalLength();
		
		Region();
};

/* Simulation
 * 
 * This class contains all the variables needed to step forward a given simulation run in time. Keeping it separate like this allows for cloning simulations,
 * doing branched runs, etc.
 * 
 * Optimizations like the reaction hash tables are stored globally however, to avoid duplicating large data structures and to share speedups across many runs.
 * Because of this, the iterator associated with simulations needs access to the parent ChemistryComputation
 * 
 */
 
class Simulation
{
	public:
		unordered_map<string, int> regionMap;		
		vector<Region> regions;
		ChemistryComputation *parentChem;
		unordered_set<string> knockouts;
		
		void connectToChemistry(ChemistryComputation *C);
		void Iterate(IterationParams &I);
		void writePopulationText(char *directory);
};

class Outcome
{
	public:
		bool reacted;
		vector<string> products;
		double rate;
};

/* ChemicalData
 */
 
class ChemicalData
{
	public:
		float energy;
};

/* Reaction
 * 
 * A given chemical reaction
 */
 
class Reaction
{
	public:
		vector<string> products;
		float weight; // selection weight for this reaction
		int typeidx; // Index for the type of this reaction
};

/* ReactionList
 * 
 * Class containing a list of reactions which a given pair of compounds can undergo
 * 
 */
 
class ReactionList
{
	public:
		vector<Reaction> rList;
		
		Outcome randomOutcome();
};

/* SimulationRequest
 * 
 * This class contains a specific kind of simulation, as well as the relevant data from the simulation to be used in processing
 * 
 * Each AnalysisRequest can specify a certain kind of SimulationRequest. Shared requests are grouped when possible
 * 
 * Question: is it possible to set up things like 'do a simulation for X time, then use that as the initial condition for Y other simulation'?
 * 
 */

class SimulationRequest
{
	public:
		string simType;
		virtual void setupSimulation(ChemistryComputation &C)=0;
		virtual bool Iterate(ChemistryComputation &C)=0; // Run a single pass of this simulation type. Returns true if the simulation has ended. Some analyses run every iteration (or every N iterations), whereas others run on completion
		void doSimulation(ChemistryComputation &C);
		virtual void executeSimulation(ChemistryComputation &C) { doSimulation(C); }
		virtual SimulationRequest *clone()=0;
		void parse(ifstream &file);
		
		unordered_map<string, double> numParams;
		unordered_map<string, string> strParams;
};

/* ChemistryComputation
 * 
 * This class contains all the various hashes and data structures used to speed up computing the properties of compounds and their reactions
 * A given simulation should probably only have one instance of this, as there aren't many situations in which separating two different
 * sets of chemistries would be all that useful.
 * 
 * Because this structure is central to storing the model of the chemistry, it contains the functions which evaluate compounds and their
 * reactions, so a simulation must start by preparing and retaining an instance of this class in a global manner (for that simulation).
 */

class AnalysisRequest;

class ChemistryComputation 
{
	public:
		Library L;
		vector<ReactionRule> reactions;
		WeightingType singletRate; // Whether singlet reactions are weighted by molecule count or total length
		WeightingType doubletRate; // Whether doublet reactions are weighted by molecule count or total length
		
		Simulation simTemplate; 
		string curSimType;
		
		unordered_map<string, ChemicalData> compoundHash; // For a given compound string, looks up the energy, chemical vector, etc
		unordered_map<string, ReactionList> reactionHash; // The accessor to this hash should be the two compounds, separated by ',', listed in lexographic order

		vector<SimulationRequest*> jobs;
		vector<AnalysisRequest*> analyses;
		
		void parseSingleConjugateSet(ifstream& file, int setidx);
		void parseConjugates(ifstream& file);
		void parseLibrary(ifstream& file);
		void parseRules(ifstream& file);
		void parseSingleRule(ifstream& file);
		void parseInitial(ifstream& file);
		void parseControl(ifstream& file);
		void parseConfig(ifstream& file);
		int lookupRegion(string str);
		ChemistryComputation();
		
		Outcome getReactionProducts(vector<string> &reactants); // Using a reference here means we don't have to allocate/free memory as much
		string getReactionString(vector<string> &reactants); // Using a reference here means we don't have to allocate/free memory as much
};

extern vector<SimulationRequest*> registeredSimulations;
extern void registerSimulations();
