#pragma once
#include "FastMat.h"
#include "Restaurant.h"
class CompTask : public Task
{
public:
	atomic<int> taskid;
	int nchunks;
	Vector& c; //Cluster Labels
	Vector& z; //Component Labels
	Vector& u; //Component u
	Matrix& x; //Data Matrix
	vector<Restaurant>& clusters;
	CompTask(Matrix& x, Vector& z, Vector& u, Vector& c,vector<Restaurant>& clusters);
	~CompTask();
	void run(int id);
	void reset(int nchunks);

};



