#include "includes.h"

/* This file contains various routines for linear algebra analysis of systems - PCA, LDA, finding eigenvalues of the transition matrix, etc */

void normalizeData(mat &data)
{
	rowvec mu = mean(data, 0);
	rowvec std = stddev(data, 1, 0);
	
	for (int i=0;i<std.n_cols;i++)
		if (std(i)==0) std(i)=1e-6;
		
	data.each_row() -= mu;
	data.each_row() /= std;
}

void doPCA(mat &data, vec &eigvals, mat &princomps, mat &scores, int maxVals)
{
	svd_econ(scores, eigvals, princomps, data);
	
	printf("%d\n",eigvals.n_rows);
	for (int i=0;i<eigvals.n_rows;i++)
		eigvals(i) = pow(eigvals(i),2)/(double)data.n_rows;
}
