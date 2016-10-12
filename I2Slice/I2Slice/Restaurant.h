#pragma once
#include "FastMat.h"
#include "Table.h"
#include "SGlobal.h"
class Restaurant 
{
public:
	Vector beta;
	Normal dist;
	Matrix sigma;
	int id = 0;
	Vector sum;
	int n;
	int nt;
	Matrix scatter;
	vector<Table*> tables;
	void sampleTables(list<Table>& mainlist,double minu);
	void addStats(Table* table);
	void addTable(Table* table);
	void setDist(Vector& mu, Matrix& sigma);
	void reset();
	void resetStats();
	void sampleParams();
	Restaurant();
	~Restaurant();
};

