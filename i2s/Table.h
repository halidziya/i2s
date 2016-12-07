#pragma once
#include <list>
#include "FastMat.h"
#include "GMMBase.h"
class Restaurant;

using namespace std;

class Table{
public:
	Normal dist;
	int n;
	Vector sum;
	Matrix scatter;
	double beta;
	double loglik0;
	list<int> plist;
	Restaurant* cls = NULL;
	int id;
	Table();
	~Table();
	Table(Vector& mu, Matrix& sigma);
	Table(Restaurant* cls, int n, Vector& sum, Matrix& scatter);
	double sampleMean();
	void reset();

	//Stut dist;
	//int npoints;
	//double logprob; // Used in sampling
	//double loglikelihood;
	//int tableid;

	//Matrix sampleScatter;
	//Vector sampleMean;
	//list<Table>::iterator copy; // To facilate copy


	//Table(void);
	//Table(int dim);
	//Table(const Table& t);
	//Table(list<Dish>::iterator d);
	//~Table(void);

	//list<Dish>::iterator dishp;
	//void addPoint(Vector& v);
	//void removePoint(Vector& v);
	//void addInitPoint(Vector& v);

	//void operator=(const Table& t);
	//void calculateCov(); // Fill covariance sampleCov , again scatter actually
	//void calculateDist();

	//friend ostream& operator<<(ostream& os, const Table& t);
};

