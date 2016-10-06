#pragma once
#include "FastMat.h"
#include "Restaurant.h"
class CompTask : public Task
{
public:
	atomic<int> taskid;
	int nchunks;
	vector<Restaurant*>& c; //Cluster Labels
	vector<Table*>& z; //Component Labels
	Vector& u; //Component u
	Matrix& x; //Data Matrix
	CompTask(Matrix& x, vector<Table*>& z, Vector& u, vector<Restaurant*>& c);
	~CompTask();
	void run(int id);
	void reset(int nchunks);

};



