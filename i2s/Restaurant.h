#pragma once
#include <list>
#include "FastMat.h"
#include "Table.h"
#include "SGlobal.h"
using namespace std;


class Restaurant
{
public:
	Vector beta;
	Normal dist;
	Stut tdist;
	Matrix sigma;
	int id = 0;
	Vector sum;
	int n;
	int nt;
	double ustar;
	Matrix scatter;
	vector<Table*> tables;
	void sampleTables(list<Table>& mainlist);
	void addStats(Table* table);
	void remStats(Table* table);
	void addTable(Table* table);
	void setDist(Vector& mu, Matrix& sigma);
	void reset();
	void resetStats();
	void calculateDist();
	void sampleParams();
	Restaurant();
	~Restaurant();
};


//class Restaurant : public virtual Task // Collection of clusters and data points
//{
//public:
//
//	int Restaurantid;
//	double likelihood;
//
//
//	list<Table> tables;
//	vector<Customer> customers;
//
//
//	void addTable(Table& t);
//	Restaurant(void);
//	
//	~Restaurant(void);
//	void run(int id);
//
//	void operator=(Restaurant& r);
//	Restaurant(Restaurant& r); // Copy constructor
//	friend ostream& operator<<(ostream& os, Restaurant& v);
//};

