class Node
{
	public:
		string value;
		double weight;
		bool isleaf;
		Node *parent;
		Node *child[2];
		
		Node *findRandomLeaf(bool invert);
		void modifyWeight(double nw);
		Node *getOtherChild(Node *T);
};

class Tree
{
	public:
		Node *Root;
		unordered_map<string, Node*> accessHash;
		void deleteLeaf(Node *N);
		Node *addNewLeaf(string val, double w);
		
		Tree();
		~Tree();
};
