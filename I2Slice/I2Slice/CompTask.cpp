#include "CompTask.h"
#include <iostream>


CompTask::CompTask(Matrix& x, Vector& z, Vector& u,Vector& c, vector<Restaurant>& clusters) : x(x),c(c),u(u),z(z),clusters(clusters)
{
}


CompTask::~CompTask()
{
}

void CompTask::run(int id) {
{
		SETUP_ID()
		int taskid = this->taskid++; // Oth thread is the main process
		auto range = trange(n, nchunks, taskid); // 2xNumber of Threads chunks		
		for (auto i = range[0]; i< range[1]; i++) // Operates on its own chunk
		{
			Restaurant& cl = clusters[c[i]];
			int NTABLE = cl.tables.size();
			Vector likelihoods(NTABLE);
			int j = 0;
			for (auto& t : cl.tables)
			{
				if (cl.beta[j] >= u[i])
				{
					likelihoods[j] = t.dist.likelihood(x(i)); //**
				}
				else
				{
					likelihoods[j] = -INFINITY;
				}
				j++;
			}
			z[i] = sampleFromLog(likelihoods); //**
		}
	}

}

void CompTask::reset(int nchunks)
{
	this->nchunks = nchunks;
	taskid = 0;
}
