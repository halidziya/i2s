#include "Restaurant.h"
#include "Table.h"

#include "Restaurant.h"
#include <iostream>


Restaurant::Restaurant()
{
	reset();
}

void Restaurant::reset()
{
	tables.resize(0);
	collapsedtables.resize(0);
	resetStats();
}

void Restaurant::resetStats()
{
	scatter = zeros(d, d);
	sum = zeros(d);
	n = 0;
	nt = 0;
	dist.mu = mu0;
}

double Restaurant::sampleParams()
{
	Matrix otherscatter = zeros(d, d);
	if (tables.size() > 0)
	{
		Vector& diff = (mu0 - (sum / tables.size()));
		otherscatter = (diff >> diff)*(kappa*tables.size() / (kappa + tables.size()));
		for (auto& atable : tables)
		{
			otherscatter += (atable->dist.mu - dist.mu) >> (atable->dist.mu - dist.mu);
		}
	}
	sigma = IWishart(Psi + scatter + otherscatter, n + m + tables.size()).rnd();
	Normal ndist((mu0*kappa + sum) / (kappa + tables.size()), sigma / (kappa + tables.size()));
	dist = Normal(ndist.rnd(), sigma / kappa1);
	return ndist.likelihood(dist.mu);
}


void Restaurant::calculateDist()
{
	double kap = harmean(nt + kappa, kappa1);
	tdist.eta = m + 2 + this->n - d - nt; //dpar0
	Vector& diff = (mu0 - (sum / tables.size()));
	Matrix otherscatter = (diff >> diff)*(kappa*tables.size() / (kappa + tables.size()));
	tdist.mu = (sum + mu0*kappa) / (nt + kappa); // dpar1
	tdist.cholsigma = ((Psi + scatter + otherscatter)*((kap + 1) / (kap*tdist.eta))).chol(); // dpar2// mu_ti is deleted because of it is same with dpar0
	tdist.calculateNormalizer();
}

void  Restaurant::addTable(Table* table)
{
	tables.push_back(table);
	addStats(table);
}

void  Restaurant::addCollapsedTable(Table* table)
{
	collapsedtables.push_back(table);
}

void Restaurant::addStats(Table* table)
{
	scatter += table->scatter;
	sum += table->dist.mu;
	n += table->n;
	nt += 1;
}

void Restaurant::remStats(Table* table)
{
	scatter -= table->scatter;
	sum -= table->dist.mu;
	n -= table->n;
	nt -= 1;

}



void Restaurant::sampleTables(list<Table>& mainlist)
{
	// Create Betas
	Vector valpha = zeros(tables.size() + 1);
	int i = 0;
	for (auto ti = tables.begin(); i < tables.size(); i++, ti++)
		valpha[i] = (*ti)->n;
	valpha[i] = alpha;
	Dirichlet dr(valpha);
	beta = dr.rnd();
	//New Sticks
	//Vector	newsticks = stickBreaker(ustar, beta[beta.n - 1], gam);
	//beta.resize(beta.n - 1);
	//beta = beta.append(newsticks);
	//Normal priormean(this->dist.mu, this->sigma / kappa1);
	//for (i = 0; i < newsticks.n; i++)
	//{
	//	mainlist.push_back(Table(priormean.rnd(), this->sigma));
	//	mainlist.back().cls = this;
	//	tables.push_back(&mainlist.back());
	//}
	i = 0;
	for (auto& atable : tables)
	{
		atable->beta = beta[i++];
	}
	ustar = 1; //For next iteration
}

void Restaurant::setDist(Vector & mu, Matrix & sigma)
{
	this->sigma = sigma;
	this->dist = Normal(mu, sigma);
}


Restaurant::~Restaurant()
{
}
