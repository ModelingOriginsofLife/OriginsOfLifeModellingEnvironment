#include "includes.h"

cObject Population_Container::selectRandom()
{
	unordered_map<string, cObject>::iterator it;
	int total=0;
	
	for (it=contents.begin();it!=contents.end();++it)
	{
		total+=it->second.count;
	}
	
	int ridx = irandom(total);
	
	it=contents.begin();
	
	while (ridx>=0)
	{
		ridx -= it->second.count;
		if (ridx>=0) ++it;
	}
	
	return it->second;
}
