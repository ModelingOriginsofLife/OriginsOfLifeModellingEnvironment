class ObjectSpec
{
	public:
		string label; // Used to identify the objectspec for dependency between modules
		vector<string> provides; // List of properties that the spec provides
		
		DataType *DataContainer; // Generic data container associated with the ObjectSpec (use templates instead?)
};
