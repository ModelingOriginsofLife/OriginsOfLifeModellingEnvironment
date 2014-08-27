#include "includes.h"

int irandom(int scale)
{
	return rand()%scale;
}

float frandom()
{
	return (rand()%1000001)/1000000.0;
}
