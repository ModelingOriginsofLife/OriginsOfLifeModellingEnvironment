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

