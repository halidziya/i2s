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
	list<Table> tables;
	void sampleTables(double minu);
	void addTable(Table& table);
	void setDist(Vector& mu, Matrix& sigma);
	void reset();
	Restaurant();
	~Restaurant();
};

