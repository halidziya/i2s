#include "Restaurant.h"
#include <iostream>


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
	nt = 0;
}

void Restaurant::sampleParams()
{
	Matrix otherscatter = zeros(d, d);
	if (tables.size() > 0)
	{
		Vector& diff = (mu0 - (sum / tables.size()));
		otherscatter = (diff >> diff)*(kappa0*tables.size() / (kappa0 + tables.size()));
		for (auto& atable : tables)
		{
			otherscatter += (atable->dist.mu - dist.mu) >> (atable->dist.mu - dist.mu);
		}
	}
	sigma = IWishart(Psi + scatter + otherscatter , n + m + tables.size()).rnd();
	dist = Normal(Normal((mu0*kappa0 + sum) / (kappa0 + tables.size()), sigma / (kappa0 + tables.size())).rnd(), sigma/kappa1);
	collapseddist = Stut((mu0*kappa0 + sum) / (kappa0 + tables.size()), (Psi + scatter + otherscatter) * ((kappa0 +1 )/(kappa0*(n + m + tables.size() -d +1 ))),n+m+tables.size()-d+1);
}


void  Restaurant::addTable(Table* table)
{
	tables.push_back(table);
	addStats(table);
}

void Restaurant::addStats(Table* table)
{
	scatter += table->scatter;
	sum += table->dist.mu;
	n += table->n;
	nt += 1;
	Matrix otherscatter = zeros(d, d);
	if (tables.size() > 0)
	{
		Vector& diff = (mu0 - (sum / tables.size()));
		otherscatter = (diff >> diff)*(kappa0*tables.size() / (kappa0 + tables.size()));
		for (auto& atable : tables)
		{
			otherscatter += (atable->dist.mu - dist.mu) >> (atable->dist.mu - dist.mu);
		}
	}
	collapseddist = Stut((mu0*kappa0 + sum) / (kappa0 + tables.size()), (Psi + scatter + otherscatter) * ((kappa0 + 1) / (kappa0*(n + m + tables.size() - d + 1))), n + m + tables.size() - d + 1);
}

void Restaurant::remStats(Table* table)
{
	scatter -= table->scatter;
	sum -= table->dist.mu;
	n -= table->n;
	nt -= 1;
	Matrix otherscatter = zeros(d, d);
	if (tables.size() > 0)
	{
		Vector& diff = (mu0 - (sum / tables.size()));
		otherscatter = (diff >> diff)*(kappa0*tables.size() / (kappa0 + tables.size()));
		for (auto& atable : tables)
		{
			otherscatter += (atable->dist.mu - dist.mu) >> (atable->dist.mu - dist.mu);
		}
	}
	collapseddist = Stut((mu0*kappa0 + sum) / (kappa0 + tables.size()), (Psi + scatter + otherscatter) * ((kappa0 + 1) / (kappa0*(n + m + tables.size() - d + 1))), n + m + tables.size() - d + 1);
}



void Restaurant::sampleTables(list<Table>& mainlist,double ustar)
{
	// Create Betas
	Vector valpha = zeros(tables.size()+1);
	int i = 0;
	for (auto ti = tables.begin(); i < tables.size(); i++, ti++)
		valpha[i] = (*ti)->n;
	valpha[i] = alpha;
	Dirichlet dr(valpha);
	beta = dr.rnd();
	//New Sticks
	Vector	newsticks = stickBreaker(ustar, beta[beta.n - 1], gamma);
	beta.resize(beta.n - 1);
	beta = beta.append(newsticks);
	Normal priormean(this->dist.mu, this->sigma / kappa1);
	for (i = 0; i < newsticks.n; i++)
	{
		mainlist.push_back(Table(priormean.rnd(), this->sigma/kappa2));
		mainlist.back().cls = this;
		tables.push_back(&mainlist.back());
	}
	i = 0;
	for (auto& atable : tables)
	{
		atable->beta = beta[i++];
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
