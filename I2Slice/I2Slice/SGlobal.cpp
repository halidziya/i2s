#include "FastMat.h"

double kappa0 = 0;
double kappa1 = 0;
int MAXCOMP = 20;

 Normal priormean;
 IWishart priorcov;

Vector stickBreaker(double ustar, double betastar = 1.0, double alpha = 1)
{
	Dirichlet  ds(v({ 1, alpha })); // Beta
	Vector lengths(MAXCOMP);
	for (int i = 0; i < MAXCOMP; i++)
		lengths[i] = ds.rnd()[0];
	Vector beta = zeros(MAXCOMP + 1);
	int i = 0;
	double totallength = betastar;
	double betasum = 0;
	for (i = 0; i<MAXCOMP; i++)
	{
		beta[i] = betastar*lengths[i];
		betasum += beta[i];
		betastar = betastar*(1 - lengths[i]);
		if (betastar < ustar)
			break;
	}
	beta.resize(i + 2);
	beta[i + 1] = totallength - betasum;
	return beta;
}