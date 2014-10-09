enum WeightingType { WEIGHT_HEAVY, WEIGHT_LIGHT, WEIGHT_HEAVYLENGTH };

class Node
{
	public:
		string value;
		double weight, lweight;
		bool isleaf;
		Node *parent;
		Node *child[2];
		
		Node *findRandomLeaf(WeightingType wType);
		void modifyWeight(double nw, double nwl);
		Node *getOtherChild(Node *T);
		void writeContents();
		Node();
};

class Tree
{
	public:
		Node *Root;
		unordered_map<string, Node*> accessHash;
		void deleteLeaf(Node *N);
		Node *addNewLeaf(string &val, double w);
		void addSubtree(Node *N);
		
		Tree();
		Tree(const Tree &other);
		~Tree();
		
		Tree &operator=(const Tree &other);
};
