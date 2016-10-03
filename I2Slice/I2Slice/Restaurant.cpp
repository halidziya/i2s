#include "Restaurant.h"



Restaurant::Restaurant()
{
}

void Restaurant::reset()
{
	tables.resize(0);
}

void  Restaurant::addTable(Table& table)
{
	tables.push_back(table);
}


void Restaurant::sampleTables(double ustar)
{
	// Create Betas
	Vector valpha = zeros(tables.size()+1);
	int i = 0;
	for (auto ti = tables.begin(); i < tables.size(); i++, ti++)
		valpha[i] = ti->n;
	valpha[i] = gamma;
	Dirichlet dr(valpha);
	beta = dr.rnd();
	//New Sticks
	Vector	newsticks = stickBreaker(ustar, beta[beta.n - 1], gamma);
	beta.resize(beta.n - 1);
	beta = beta.append(newsticks);
	Normal priormean(this->dist.mu, this->sigma / kappa1);
	for (i = 0; i < newsticks.n; i++)
	{
		tables.push_back(Table(priormean.rnd(), this->sigma));
	}
}

void Restaurant::setDist(Vector & mu, Matrix & sigma)
{
	this->sigma = sigma;
	this->dist = Normal(mu, sigma);
}


Restaurant::~Restaurant()
{
}
