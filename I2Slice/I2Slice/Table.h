#pragma once
#include "FastMat.h"
class Restaurant;
class Table
{
public:
	Normal dist;
	int n;
	Vector sum;
	Matrix scatter;
	double u=0;
	
	Restaurant* cls=NULL;
	int id;
	Table();
	~Table();
	Table(Vector& mu, Matrix& sigma);
	Table(Restaurant* cls,int n,Vector& sum,Matrix& scatter);
	void sampleMean();
	void reset();
};

