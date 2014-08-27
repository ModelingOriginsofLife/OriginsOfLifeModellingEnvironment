class DataType
{
	public:
		string typelabel();
		DataType();
		
		void *accessContents();
		DataType *newMember();
};

extern unordered_map<string, (DataType *)(*)()> registeredDataTypes; // New data types should add themselves to this list at runtime

extern DataType *createNewDataType(string label);
