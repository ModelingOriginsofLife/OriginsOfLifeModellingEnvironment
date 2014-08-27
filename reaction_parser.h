class Library
{
	public:
		vector<char> atomTypes;
		
		unordered_map<char, vector<char>> atomClasses;

		string conjugateString(string S, int idx);
		char applyConjugation(char C, int idx);
		bool isMemberOfSet(char C, char label, int idx);		
};

class Symbol
{
	public:
		bool isWild;
		bool isSingleton;
		bool isReversed;
		bool isBound;
		int doConjugate;
		
		char label;
		
		void clear();
};

class SearchSubtree;

class SearchSubtree
{
	public:
		int stringpos, rulepos, isNull, isEnd;
		unordered_map<char, string> bound;
		vector<SearchSubtree> subTrees;
				
		bool isMatch(Symbol S, char C, Library L);		
		void parseOnLeaves(vector<Symbol> rule, string matchstr, Library L);
		vector<vector<string>> generateProducts(vector<vector<Symbol>> rule, Library L);
		void applyRule(vector<Symbol> rule, string matchstr, Library L);
		void branchTree(int offset, int strpos, vector<Symbol> rule, string matchstr, Library L);
		void writeValidLeaves();
		string getProduct(vector<Symbol> rule, Library L);
};

class ReactionRule
{
	public:
		string rule;
		vector<vector<Symbol>> reacRules, prodRules;
		int Nreac, Nprod;
		double k; // Underlying rate constant
		
		vector<vector<string>> (*overrideProducts)(vector<string> reactants, Library L); // Support for a custom function here; still takes a library though
		
		vector<vector<string>> getAllProducts(vector<string> reactants, Library L);
		void parseRule(Library L);
};
