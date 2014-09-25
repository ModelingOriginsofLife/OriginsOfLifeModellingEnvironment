#include "includes.h"

/* This file contains various routines for linear algebra analysis of systems - PCA, LDA, finding eigenvalues of the transition matrix, etc */

void normalizeData(mat &data)
{
	vec mu = mean(data, 0);
	vec std = stddev(data, 1, 0);
	
	data.each_col() -= mu;
	data.each_col() /= std;
}

void doPCA(mat &data, vec &eigvals, mat &princomps, mat &scores, int maxVals)
{
	svd_econ(scores, eigvals, princomps, data);
}
