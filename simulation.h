class ChemistryComputation;

class Region
{
	public:
		Tree population;
		Tree bath; // Compounds that are assumed to always be available, and their relative concentrations
		
		void addCompound(string str, int count);
		void removeCompound(string str, int count);
		string pickRandomCompound();
		void writePopulationText(char *FName);
		double getConcentration(string compound);
		
		void doRandomSinglet(ChemistryComputation *C);
		void doRandomDoublet(ChemistryComputation *C);
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
		
		void Iterate(ChemistryComputation *C);
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

/* ChemistryComputation
 * 
 * This class contains all the various hashes and data structures used to speed up computing the properties of compounds and their reactions
 * A given simulation should probably only have one instance of this, as there aren't many situations in which separating two different
 * sets of chemistries would be all that useful.
 * 
 * Because this structure is central to storing the model of the chemistry, it contains the functions which evaluate compounds and their
 * reactions, so a simulation must start by preparing and retaining an instance of this class in a global manner (for that simulation).
 */
  
class ChemistryComputation 
{
	public:
		Library L;
		vector<ReactionRule> reactions;

		Simulation simTemplate; 
		unordered_map<string, ChemicalData> compoundHash; // For a given compound string, looks up the energy, chemical vector, etc
		unordered_map<string, ReactionList> reactionHash; // The accessor to this hash should be the two compounds, separated by ',', listed in lexographic order

		void parseLibrary(ifstream& file);
		void parseRules(ifstream& file);
		void parseSingleRule(ifstream& file);
		void parseInitial(ifstream& file);
		void parseControl(ifstream& file);
		int lookupRegion(string str);
		
		Outcome getReactionProducts(vector<string> reactants);
		string getReactionString(vector<string> reactants);
};
