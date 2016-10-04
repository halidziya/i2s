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
	Matrix scatter;
	list<list<Table>::iterator> tables;
	void sampleTables(double minu);
	void addStats(list<Table>::iterator table);
	void addTable(list<Table>::iterator table);
	void setDist(Vector& mu, Matrix& sigma);
	void reset();
	void resetStats();
	void sampleParams();
	list<Table>& tablestorage;
	Restaurant(list<Table>& tablestorage);
	~Restaurant();
};

