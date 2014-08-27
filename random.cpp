#include "includes.h"

int mrandom()
{
	return rand();
}

int irand(int x)
{
	return mrandom()%x;
}

float frand(float x)
{
	return x*irand(1000001)/1000000.0;
}

bool prand(double P)
{
	if (frand(1) < P) return true;
	return false;
}
