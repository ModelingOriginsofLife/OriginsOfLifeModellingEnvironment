class Library
{
	public:
		vector<char> atomTypes;
		vector< vector<char> > conjugates; // Conjugation 0 is always the identity, so is a special case. A vector is used here for efficient access compared to a hash.
		unordered_map<char, vector<char>> atomClasses;

		void expandConjugates(int idx);
		string conjugateString(string S, int idx);
		char applyConjugation(char C, int idx);
		bool isMemberOfSet(char C, char label, int idx);		
		char getMemberOfSet(char label, int idx);		
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
				
		bool isMatch(Symbol S, char C, Library &L);		
		void parseOnLeaves(vector<Symbol> &rule, string matchstr, Library &L);
		vector<vector<string>> generateProducts(vector<vector<Symbol>> &rule, Library &L);
		void applyRule(vector<Symbol> &rule, string matchstr, Library &L);
		void branchTree(int offset, int strpos, vector<Symbol> &rule, string matchstr, Library &L);
		void writeValidLeaves();
		string getProduct(vector<Symbol> &rule, Library &L);
};

class ReactionRule
{
	public:
		string rule;
		vector<vector<Symbol>> reacRules, prodRules;
		int Nreac, Nprod;
		double k; // Underlying rate constant
		
		vector<vector<string>> (*overrideProducts)(vector<string> &reactants, Library &L); // Support for a custom function here; still takes a library though
		
		vector<vector<string>> getAllProducts(vector<string> &reactants, Library &L);
		void parseRule(Library &L);
};

extern string reverseString(string S);
