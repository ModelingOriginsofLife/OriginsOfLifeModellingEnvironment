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
		unordered_map<string, ChemicalData> compoundHash; // For a given compound string, looks up the energy, chemical vector, etc
		unordered_map<string, ReactionList> reactionHash; // The accessor to this hash should be the two compounds, separated by ',', listed in lexographic order
};
