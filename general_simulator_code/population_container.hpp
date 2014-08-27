class cObject
{
	public:
		string label;
		Object *prototype;
		int count;
};

class Population_Container: public Container
{
	public:
		unordered_map<string, cObject> contents;
		
		cObject selectRandom();
};
