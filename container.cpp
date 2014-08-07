class Accessor
{
	public:
		Object *prototype; // If this is a regular container, writing to this is okay; otherwise, don't!
		string label; // Only if population_container!
		
		void (*decrement); // Function pointer to reduce the count of this object
		void (*increment); // Function pointer to increase the count of this object
};

class Container
{
	public:
		vector<Object *> contents;
		
		Accessor selectRandom();
		
};
