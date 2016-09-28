#include "Table.h"
#include "Restaurant.h"


Table::Table()
{
	// Sample from Prior
}

Table::Table(Vector& mu,Matrix& sigma)
{
	dist = Normal(mu,sigma);
	n = 0;

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
