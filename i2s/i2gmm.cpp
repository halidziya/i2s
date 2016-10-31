#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "FastMat.h"
#include "Restaurant.h"
#include <thread>  
#include <iostream>
#include <algorithm>


using namespace std;


int MAX_SWEEP=1500;
int BURNIN=1400;
int NSAMPLE = 10;
int STEP=(MAX_SWEEP-BURNIN)/NSAMPLE; // Default value is 10 sample + 1 post burnin
char* result_dir = "./";
int NINITIAL = 2;
// Variables
double kep,eta;
Vector u;
vector<Restaurant*> c;
vector<Table*> z;
Matrix x;
list<Table> tables;
list<Restaurant> clusters;
Vector loglik0;

class CompTask : public Task
{
public:
	atomic<int> taskid;
	int nchunks;
	vector<Table* > tables;
	void run(int id) {
		{
			SETUP_ID()
			int taskid = this->taskid++; // Oth thread is the main process
			auto range = trange(n, nchunks, taskid); // 2xNumber of Threads chunks		

			for (auto i = range[0]; i< range[1]; i++) // Operates on its own chunk
			{
				Restaurant& cl = *c[i];
				int NTABLE = cl.tables.size();
				Vector likelihoods(NTABLE);
				int j = 0;
				for (auto& t : cl.tables)
				{
					if ( (t->beta >= u[i]))
					{
						likelihoods[j] = t->dist.likelihood(x(i)); //**
					}
					else
					{
						likelihoods[j] = -INFINITY;
					}
					j++;
				}
				z[i] = cl.tables[sampleFromLog(likelihoods)]; //**
			}
		}

	}

	void reset(int nchunks)
	{
		this->nchunks = nchunks;
		taskid = 0;
		tables.resize(0);
		for (auto& t : ::tables) // Global tables
		{
			tables.push_back(&t);
		}
	}
};

class Collector : public Task
{
public:
	atomic<int> taskid;
	int ntable;
	Collector(int ntable) : ntable(ntable) {
	}
	void reset() {
		for (auto& atable : tables)
			atable.reset();
		taskid = 0;
	}
	void subscatter() {
		for (auto& atable : tables)
			atable.scatter -= ((atable.sum >> atable.sum) / atable.n); // Actual Scatter
		taskid = 0;
	}
	void run(int id) {
		// Use thread specific buffer
		SETUP_ID()
			int taskid = this->taskid++;
		for (int i = 0; i < n; i++)
		{
			if (z[i]->id == taskid) // Distributed according to label id
			{
				z[i]->n += 1;
				z[i]->sum += x(i);
				z[i]->scatter += x(i) >> x(i); // Second moment actually
				z[i]->loglik0 += loglik0[i];
			}
		}

	}
};


class Cluster : public Task
{
public:
	atomic<int> taskid;
	Vector& logprobs;
	Table& atable;
	vector<Restaurant*>& cp;
	Cluster(Vector& logprobs,Table& atable, vector<Restaurant*>& cp) : logprobs(logprobs),atable(atable),cp(cp) {
		taskid = 0;
	}
	void run(int id) {
		// Use thread specific buffer
		SETUP_ID()
		int taskid = this->taskid++;
		double logprob;
			auto cc = cp[taskid];
			if (cc->nt == 0)
				logprobs[taskid] = -INFINITY;
			else
			{
				logprob = 0;
				for (auto points = atable.plist.begin(); points != atable.plist.end(); points++)
				{
					logprob += cc->tdist.likelihood(x[*points]);
				}
				logprob = logprob + log(cc->n);//+log(atable.n); //Prior
				logprobs[taskid] = logprob;
			}
	}
};

void reid(list<Table>& tables)
{
	int i = 0;
	for (auto& atable : tables)
		atable.id = i++;
}

void reid(list<Restaurant>& clusters)
{
	int i = 0;
	for (auto& acul : clusters)
		acul.id = i++;
}

Matrix SliceSampler(Matrix& data, ThreadPool& workers, Matrix& superlabels)
{

	Vector priormean = mu0;
	Matrix priorvariance(d, d);
	priorvariance = Psi*((kep + 1) / ((kep)*eta));
	Stut stt(priormean, priorvariance, eta);
	loglik0 = stt.likelihood(data);
	// Point level variables
	x = data;
	int NTABLE = NINITIAL;
	int k = 0, i = 0, j = 0;
	u = ones(n);
	z = vector<Table*>(n); // Component labels
	c = vector<Restaurant*>(n);

	clusters = list<Restaurant>(NINITIAL);
	for (auto& cluster : clusters)
	{
		cluster.sampleParams();
		cluster.ustar = 0.05;
		cluster.sampleTables(tables);
	}

	vector<Restaurant*> cp = vector<Restaurant*>(clusters.size());
	i = 0;
	for (auto& cluster : clusters)
	{
		cp[i] = &cluster;
		i++;
	}

	for (int i = 0; i < n; i++)
	{
		c[i] = cp[rand()%NINITIAL];
		j = rand() % c[i]->tables.size();
		u[i] = c[i]->beta[j]*urand();
		z[i] = c[i]->tables[j];
	}

	Matrix zi((MAX_SWEEP - BURNIN) / STEP + 1, n);
	CompTask cmsampler;
	Collector  collector(NTABLE);

	for (int iter = 0; iter <= MAX_SWEEP; iter++) {
		reid(clusters);
		cmsampler.reset(2 * nthd);
		for (auto i = 0; i < cmsampler.nchunks; i++) {
			workers.submit(cmsampler);
		}
		workers.waitAll();

		reid(tables);
		NTABLE = tables.size();
		collector.reset();
		for (i = 0; i < NTABLE; i++) {
			workers.submit(collector);
		}
		workers.waitAll();
		collector.subscatter();


		// Remove if tables are not used 
		for (auto& cluster : clusters)
			for (auto tt = cluster.tables.begin(); tt != cluster.tables.end();)
			{
				if ((*tt)->n == 0)
					tt = cluster.tables.erase(tt);
				else
					tt++;
			}

		for (auto tt = tables.begin(); tt != tables.end();)
		{
			if (tt->n == 0)
				tt = tables.erase(tt);
			else
				tt++;
		}

		// Update table means
		for (auto& table : tables)
			table.sampleMean();

		reid(clusters);
		// Sample tables

		cp = vector<Restaurant*>(clusters.size());
		vector<int> cidx = vector<int>(tables.size());
		i = 0;
		for (auto& cluster : clusters)
		{
			cp[i] = &cluster;
			i++;
		}


		for (auto& cluster : clusters)
			cluster.reset();

		for (auto& atable : tables)
		{
			atable.cls->addTable(&atable);
		}


		if ((iter % 20) == 0)
		{
			cout << iter << endl;
			for (auto cc = clusters.begin(); cc != clusters.end(); cc++)
				cout << " " << cc->n;
			cout << endl << "Tables :" << tables.size() << endl;
			cout << endl;
		}

		// 3rd Loop
		int kal = 1;
		double newdishprob = 0,maxdishprob,logprob;
		
		for (auto& atable : tables)
			atable.plist.resize(0);
		for (int i = 0; i < n; i++)
		{
			z[i]->plist.push_back(i);
		}
		
		
		cp = vector<Restaurant*>(clusters.size());
		i = 0;
		for (auto& cluster : clusters)
		{
			cp[i] = &cluster;
			i++;
		}

		for (auto& cluster : clusters)
			cluster.reset();
		for (auto& table : tables)
			table.cls->addTable(&table);
		for (auto& cluster : clusters)
			cluster.sampleParams();
		for (auto& cluster : clusters)
			cluster.calculateDist();
			
			for (auto& atable : tables)
			{
				Vector logprobs(clusters.size() + 1);
				atable.cls->remStats(&atable);
				atable.cls->calculateDist();

				k = 0;
				newdishprob =  log(gamma) + atable.loglik0;
				Cluster cls(logprobs, atable,cp);
				for (auto i = 0; i < clusters.size(); i++) {
					workers.submit(cls);
				}
				workers.waitAll();


				//for (auto cc = clusters.begin(); cc != clusters.end();cc++)
				//{
				//	if (cc->nt == 0)
				//		logprobs[k] = -INFINITY;
				//	else
				//	{
				//		logprob = 0;
				//		// Parallelize here !!!! 
				//		for (auto points = atable.plist.begin(); points != atable.plist.end(); points++)
				//		{
				//			logprob += cc->tdist.likelihood(x[*points]);
				//		}
				//		logprob = logprob + log(cc->nt) * atable.n; //Prior
				//		logprobs[k] = logprob;
				//		
				//	}

				//	if (false)
				//	{
				//		cc->tdist.mu.print();
				//		cout << cc->id << ":";
				//		cout << logprobs[k] << " ";
				//		cout << (cc->tdist.mu - atable.dist.mu).norm() << " ";
				//		cout << ((cc->sum/n) - atable.dist.mu).norm() << " ";
				//		cout << endl;
				//	}
				//	k++;
				//}
				logprobs[clusters.size()] = newdishprob;
				//if (iter > BURNIN)
				//{
				//	logprobs.print();
				//	atable.dist.mu.print();
				//	cout << newdishprob;

				//}

				int idx = sampleFromLog(logprobs);
				//if (iter > BURNIN)
				//	cout << endl << "Selected " << idx << endl;
	
				if (idx == clusters.size())
				{
					clusters.push_back(Restaurant());
					clusters.back().addTable(&atable);
					clusters.back().id = clusters.size() - 1;
					clusters.back().calculateDist();
					clusters.back().sampleParams();
					atable.cls = &clusters.back();
					cp.push_back(&clusters.back());
					
				}
				else
				{
					atable.cls = cp[idx];
					
					atable.cls->addStats(&atable);
					atable.cls->calculateDist();
					atable.cls->sampleParams();
				}
			}

			for (auto cc = clusters.begin(); cc != clusters.end();)
			{

				if (cc->n == 0)
				{
					cc = clusters.erase(cc);
				}
				else
					cc++;
			}



			for (auto& cluster : clusters)
				cluster.reset();
			for (auto& table : tables)
				table.cls->addTable(&table);

		for (auto& cluster : clusters)
			cluster.sampleParams();
		reid(clusters);
		reid(tables);
		for (int i = 0; i < n; i++)
			c[i] = z[i]->cls;

		if (iter > BURNIN && ((iter - BURNIN) % STEP == 0))
			for (int i = 0; i < n; i++)
			{
				zi((iter - BURNIN) / STEP)[i] = z[i]->id;
				superlabels((iter - BURNIN) / STEP)[i] = c[i]->id;
			}

		// Sample Tables
		for (auto& cluster : clusters)
			if (cluster.n > 0)
				cluster.sampleTables(tables);

		for (i = 0; i < n; i++)
		{
			u[i] = z[i]->beta * urand();
			if (z[i]->cls->ustar < u[i])
				z[i]->cls->ustar = u[i];
		}
	}
	return zi;
}




PILL_DEBUG
int main(int argc,char** argv)
{
	Matrix x;
	Matrix hyperparams;
	Matrix initialLabels;
	generator.seed(time(NULL));
	srand(time(NULL));
	system("dir");

	CBLAS = 0; // Do not use vector library in small dimensions

	if (argc > 1)
	{
		x.readBin(argv[1]);
		cout << argv[1] << endl;
	}
	else
	{
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (d,m,kappa0,kappa1,gamma)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
		return -1;
	}
	cout << "NPOINTS :" << x.r << " NDIMS:" << x.m << endl;
	nthd = thread::hardware_concurrency();
	n = x.r; // Number of Points
	d = x.m;

	init_buffer(nthd, x.m);
	cout << " Available number of threads : " << nthd << endl;
	precomputeGammaLn(2 * n + 100 * d);


	// Hyper-parameters with default values
	if (x.data == NULL)
	{
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (d,m,kappa,gamma)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
		return -1;
	}

	if (argc > 2)
	{
		Matrix mu;
		mu.readBin(argv[2]);
		mu0 = mu;
	}
	else
		mu0 = x.mean().copy();


	if (argc > 4)
	{
		hyperparams.readBin(argv[4]);
		hyperparams.print();
		m = hyperparams.data[1];
		kappa0 = hyperparams.data[2];
		kappa1 = hyperparams.data[3];
		alpha = hyperparams.data[4];
		gamma = hyperparams.data[5];
		cout << m << " " << kappa0 << " " << kappa1 << " " << alpha << " " << gamma << endl;
	}
	else
	{
		m = x.m + 3;
		kappa0 = 1;
		kappa1 = 1;
		gamma = 1;
		alpha = 1;
	}

	if (argc > 3)
		Psi.readBin(argv[3]);
	else
		Psi = (eye(d)).copy();




	if (argc > 5)
		MAX_SWEEP = atoi(argv[5]);

	if (argc > 6)
		BURNIN = atoi(argv[6]);

	if (argc > 7)
	{
		NSAMPLE = atoi(argv[7]);
		STEP = (MAX_SWEEP - BURNIN) / NSAMPLE;
	}

	if (argc > 8)
		CBLAS = atoi(argv[8]);


	printf("Reading...\n");
	kep = kappa0*kappa1 / (kappa0 + kappa1);
	eta = m - d + 2;
	precomputeGammaLn(2 * (n + d) + 1);  // With 0.5 increments
	init_buffer(thread::hardware_concurrency(), d);

	ThreadPool tpool(nthd);
	Matrix superlabels((MAX_SWEEP - BURNIN) / STEP + 1, n);
	auto labels = SliceSampler(x, tpool, superlabels); // data,m,kappa,gamma,mean,cov 
	string filename = argv[1];
	labels.writeBin(filename.append(".labels").c_str());
	filename = argv[1];
	superlabels.writeBin(filename.append(".superlabels").c_str());

}
