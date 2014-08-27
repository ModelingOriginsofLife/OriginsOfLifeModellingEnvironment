#include "includes.h"

int Container::selectRandom()
{
	return contents[irandom(contents.size())];
}
