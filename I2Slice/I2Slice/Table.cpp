#include "Table.h"
#include "Restaurant.h"


Table::Table()
{
	dist = Normal(mu0, Psi);
	reset();
}

Table::Table(Vector& mu,Matrix& sigma)
{
	dist = Normal(mu,sigma);
	reset();

}

void Table::reset()
{
	this->n = 0;
	this->sum = zeros(d);
	this->scatter = zeros(d,d);
}

Table::Table(Restaurant* cls, int n, Vector& sum, Matrix& scatter)
{
	this->n = n;
	this->sum = sum;
	this->scatter = scatter;
	this->cls = cls;
	sampleMean();
}


void Table::sampleMean()
{
	Vector& mu = Normal((cls->dist.mu*kappa1 + sum) / (n + kappa1), cls->sigma / (n + kappa1)).rnd();
	dist = Normal(mu,cls->sigma);
}


Table::~Table()
{
}
