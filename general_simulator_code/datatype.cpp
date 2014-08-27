#include "includes.h"

/* 
 * DataType
 * 
 * An abstract container type for different kinds of data.
 * Sub-classes are specific types of container
 * 
 */

DataType::DataType()
{
}

string DataType::typelabel()
{
	return "empty";
}

void *DataType::accessContents()
{
	return NULL;
}

DataType *DataType::newMember()
{
	return new DataType;
}

/*
 * List of registered data types
 */
 
unordered_map<string, DataType *(*)()> registeredDataTypes; // New data types should add themselves to this list at runtime

DataType *createNewDataType(string label)
{
	return registeredDataTypes[label]();
}
