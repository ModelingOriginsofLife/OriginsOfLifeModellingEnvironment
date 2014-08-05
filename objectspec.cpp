#include "includes.h"

/*
 * ObjectSpec
 * 
 * - Stores details of what constitutes a simulation object. This is constructed from the choice of Reactor and Arranger
 * - Each ObjectSpec contains certain specific information
 * - ObjectSpecs can interact along shared 'models'. For example, a string chemistry can associate something with every place on the string
 * - These interactions are encoded ???
 */
 
class ObjectSpec
{
	public:
		string label; // Used to identify the objectspec for dependency between modules
		vector<string> provides; // List of properties that the spec provides
		
		DataType *DataContainer; // Generic data container associated with the ObjectSpec (use templates instead?)
};
