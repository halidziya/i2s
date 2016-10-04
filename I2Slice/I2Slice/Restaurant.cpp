#include "Restaurant.h"



Restaurant::Restaurant()
{
	reset();
}

void Restaurant::reset()
{
	tables.resize(0);
	resetStats();
}

void Restaurant::resetStats()
{
	scatter = zeros(d, d);
	sum = zeros(d);
	n = 0;
}

void Restaurant::sampleParams()
{
	Matrix otherscatter = zeros(d, d);
	sigma = IWishart(Psi + scatter + otherscatter, n + m).rnd();
	dist = Normal(Normal((mu0*kappa0 + sum) / (kappa0 + tables.size()), sigma / kappa1).rnd(), sigma);
}


void  Restaurant::addTable(Table& table)
{
	tables.push_back(table);
	addStats(table);
}

void Restaurant::addStats(Table& table)
{
	scatter += table.scatter;
	sum += table.dist.mu;
	n += table.n;
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
		Table t = Table(priormean.rnd(), this->sigma);
		t.cls = this;
		tables.push_back(t);
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
