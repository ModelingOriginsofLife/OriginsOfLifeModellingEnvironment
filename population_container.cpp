class cObject
{
	public:
		Object *prototype;
		int count;
};

class Population_Container: public Container
{
	public:
		unordered_map<string, cObject> contents;
		
		Accessor selectRandom();
};
