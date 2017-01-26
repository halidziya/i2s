#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "FastMat.h"
#include "Restaurant.h"
#include <thread>  
#include <iostream>
#include <algorithm>
#include "Algorithms.h"
#include "GMMBase.h"


using namespace std;


int MAX_SWEEP = 1500;
int BURNIN = 1000;
int SAMPLE = 20;
int STEP = (MAX_SWEEP - BURNIN) / SAMPLE;

char* result_dir = "./";
int NINITIAL = 3; // Should be smaller than d in current matrix implemetation
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
	atomic<int> nnew;
	int nchunks;
	vector<Table* > tables;

	void run(int id) {
		{
			SETUP_ID()
			int taskid = this->taskid++; // Oth thread is the main process
			auto range = trange(n, nchunks, taskid); // 2xNumber of Threads chunks		
			int NTABLE = 100;
			Vector likelihoods(NTABLE + 1);
			int j = 0;
			for (auto i = range[0]; i< range[1]; i++) // Operates on its own chunk
			{
				Restaurant& cl = *c[i];
				NTABLE  = cl.tables.size();
				likelihoods.resize(NTABLE + 1);
				j = 0;
				for (auto& t : cl.tables)
				{
					likelihoods[j] = t->dist.likelihood(x(i)) +  log(t->beta); //**
					j++;
				}
				
				likelihoods[NTABLE] = loglik0[i] + log(cl.beta[NTABLE]);
				int idx = sampleFromLog(likelihoods);
				if (idx < NTABLE)
				{
					z[i] = cl.tables[idx]; //**
				}
				else
				{
					z[i] = NULL;
					nnew++;
				}
			}
		}

	}

	void reset(int nchunks)
	{
		this->nchunks = nchunks;
		taskid = 0;
		nnew = 0;
		tables.resize(0);
		for (auto& t : ::tables) // Global tables
		{
			tables.push_back(&t);
		}
	}
};

class ListCollect : public Task
{
public:
	atomic<int> taskid;
	atomic<double> loglikelihood;
	int nchunks;
	void reset()
	{
		taskid = 0;
	}
	void run(int id) {
		SETUP_ID()
		int taskid = this->taskid++;
		for (int i = 0; i < n; i++)
		{
			if (z[i]->id == taskid) // Distributed according to label id
			{
				z[i]->plist.push_back(i);
			}
		}
	}

};

class CurrentLikelihood : public Task
{
public:
	atomic<int> taskid;
	atomic<double> loglikelihood;
	int nchunks;
	void run(int id) {
			double llikelihood=0;
			SETUP_ID()
			int taskid = this->taskid++; // Oth thread is the main process
			auto range = trange(n, nchunks, taskid); // 2xNumber of Threads chunks

			for (auto i = range[0]; i < range[1]; i++) // Operates on its own chunk
			{
				llikelihood = llikelihood + z[i]->dist.likelihood(x(i));
			}
			loglikelihood = loglikelihood + llikelihood;
		}
		void reset(int nc){
			loglikelihood = 0;
			nchunks = nc;
			taskid = 0;
		}

};

double getFullLikelihood(ThreadPool& tp)
{
	CurrentLikelihood cl;
	cl.reset(nthd);
	for (int i = 0; i < nthd; i++) {
		tp.submit(cl);
	}
	tp.waitAll();
	return cl.loglikelihood;
}

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
				#pragma parallel
				for (auto apoint : atable.plist)
				{
					logprob += cc->tdist.likelihood(x[apoint]);
				}
				//logprob += cc->tdist.likelihood(atable.dist.mu);
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

bool compare_clusters(Restaurant& c1, Restaurant& c2)
{
	return c1.n > c2.n;
}


void reid(list<Restaurant>& clusters)
{
	clusters.sort(compare_clusters);
	int i = 0;
	for (auto& acul : clusters)
		acul.id = i++;
}


Matrix SliceSampler(Matrix& data, ThreadPool& workers, Matrix& superlabels)
{

	if (NINITIAL > d)
		NINITIAL = d;

	Vector priormean = mu0;
	Matrix priorvariance(d, d);
	priorvariance = Psi*((kep + 1) / ((kep)*eta));
	Stut stt(priormean, priorvariance, eta);
	loglik0 = stt.likelihood(data);
	// Point level variables
	x = data;
	int NTABLE = NINITIAL;
	int k = 0, i = 0, j = 0;
	//u = ones(n);
	z = vector<Table*>(n); // Component labels
	c = vector<Restaurant*>(n);

	Vector initiallabels = kmeans(data, NINITIAL);
	clusters = list<Restaurant>(NINITIAL);
	for (auto& cluster : clusters)
		cluster.resetStats();

	vector<Restaurant*> cp = vector<Restaurant*>(clusters.size());
	i = 0;
	for (auto& cluster : clusters)
	{
		cp[i] = &cluster;
		i++;
	}

	for (auto& cluster : clusters)
	{
		cluster.sampleParams();
		Normal priormean(cluster.dist.mu, cluster.sigma / kappa1);
		tables.push_back(Table(priormean.rnd(), cluster.sigma));
		tables.back().cls = &cluster;
		cluster.tables.push_back(&tables.back());
		cluster.beta = v({ 0.99,0.01 });
		tables.back().beta = 0.99;
		cluster.dist.mu = zeros(d);
		cluster.dist.cholsigma = eye(d);
		//cluster.sampleTables(tables);
	}

	for (int i = 0; i < n; i++)
	{
		c[i] = cp[initiallabels[i]];
		//u[i] = c[i]->beta[0]*urand();
		z[i] = c[i]->tables[0];
	}

	Matrix zi(SAMPLE, n);
	CompTask cmsampler;
	Collector  collector(NTABLE);

	reid(tables);
	NTABLE = tables.size();
	collector.reset();
	for (i = 0; i < NTABLE; i++) {
		workers.submit(collector);
	}
	workers.waitAll();
	collector.subscatter();
	// Update table means
	for (auto& table : tables)
		table.sampleMean();

	reid(clusters);
	// Sample tables

	for (auto& cluster : clusters)
		cluster.reset();
	for (auto& table : tables)
		table.cls->addTable(&table);
	for (auto& cluster : clusters)
		cluster.sampleParams();
	for (auto& cluster : clusters)
		cluster.calculateDist();
	for (auto& table : tables)
		table.sampleMean();

	//Vector kappas = v({ 0.01,0.02,0.05,0.1 });

	for (int iter = 0; iter < MAX_SWEEP; iter++) {

		reid(clusters);
		cmsampler.reset(10 * nthd);
		for (auto i = 0; i < cmsampler.nchunks; i++) {
			workers.submit(cmsampler);
		}
		workers.waitAll();
		
		//kappa = kappas[rand() % kappas.n];
		//kappa1 = 10 * kappa;
		Vector likes;
		likes.resize(100);
		for (int i = 0; i < n; i++)
		{
			if (z[i] == NULL)
			{
				likes.resize(c[i]->collapsedtables.size() + 1);
				j = 0;
				for (auto& atable : c[i]->collapsedtables)
				{
					Normal dist = Normal((atable->sum + c[i]->dist.mu*kappa1) / (kappa1 + atable->n) , (Psi + c[i]->scatter + atable->scatter)/(atable->n + c[i]->n + m )) ;
					likes[j] = dist.likelihood(x[i])  + log(atable->n);
					j++;
				}
				//mu0.print();
				likes[j] = loglik0[i] + alpha;
				int idx = sampleFromLog(likes);
				

				if (idx == c[i]->collapsedtables.size())
				{
					tables.push_back(Table());
					c[i]->addCollapsedTable(&tables.back());
					tables.back().cls = c[i];
				}
				z[i] = c[i]->collapsedtables[idx];
				z[i]->n++;
				z[i]->sum += x[i];
				z[i]->scatter += ((x[i]-c[i]->dist.mu) >> (x[i] - c[i]->dist.mu));
			}
			//cout << i << endl;
		}
		//cout << cmsampler.nnew << endl;

		for (auto& cl : clusters)
		{
			cl.tables.insert(cl.tables.begin(), cl.collapsedtables.begin(), cl.collapsedtables.end());
			cl.collapsedtables.resize(0);
		}


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
			fflush(stdout);
		}

		// 3rd Loop
		int kal = 1;
		double newdishprob = 0,maxdishprob,logprob;
		
		for (auto& atable : tables)
			atable.plist.resize(0);

		//ListCollect lc;
		//NTABLE = tables.size();
		//lc.reset();
		//for (i = 0; i < NTABLE; i++) {
			//workers.submit(lc);
		//}
		//workers.waitAll();

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
			cluster.calculateDist();
			
			for (auto& atable : tables)
			{
				Vector logprobs(clusters.size() + 1);
				atable.cls->remStats(&atable);
				atable.cls->calculateDist();

				k = 0;
				
				newdishprob = log(gam) + atable.loglik0; //stt.likelihood(atable.dist.mu);
				Cluster cls(logprobs, atable,cp);
				//newdishprob = log(gam) + stt.likelihood(atable.dist.mu); //stt.likelihood(atable.dist.mu);
				for (auto i = 0; i < clusters.size(); i++) {
					//if (cp[i]->nt == 0)
					//	logprobs[i] = -INFINITY;
					//else
					//	logprobs[i]= cp[i]->tdist.likelihood(atable.dist.mu) + log(cp[i]->n);
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

		//for (auto& cluster : clusters)
		//	cluster.sampleParams();






			if (false){//iter % 5000 == 1) {
				Vector loglikelihood(20);
				i = 0;
				double kapparatio = kappa1 / kappa;
				for (i=0;i<20;i++)
				{
					loglikelihood[i] = 0;
					kappa = 0.005 + 0.01*i;
					kappa1 = kappa * kapparatio;
					for (auto& cluster : clusters)
						loglikelihood[i] += cluster.sampleParams();
					//for (auto& table : tables)
					//	loglikelihood[i] += table.sampleMean()/tables.size();
					//loglikelihood[i] += getFullLikelihood(workers)/n;
					i++;
				}
				//loglikelihood.print();

				kappa = (sampleFromLog(loglikelihood)*0.01 + 0.005);
				kappa1 = kappa * kapparatio;
				for (auto& cluster : clusters)
					cluster.sampleParams();
				for (auto& table : tables)
					table.sampleMean();

				//i = 0;
				//for (kapparatio = 6; kapparatio < 16; kapparatio += 1)
				//{
				//	kappa1 = kapparatio*kappa;
				//	loglikelihood[i] = 0;
				//	for (auto& cluster : clusters)
				//		loglikelihood[i] += cluster.sampleParams();
				//	//for (auto& table : tables)
				//	//	loglikelihood[i] += table.sampleMean()/tables.size();
				//	//loglikelihood[i] += getFullLikelihood(workers)/n;
				//	i++;
				//}
				//kapparatio = (sampleFromLog(loglikelihood)*1 + 6);
				//kappa1 = kapparatio*kappa;
				//for (auto& cluster : clusters)
				//	cluster.sampleParams();
				//for (auto& table : tables)
				//	table.sampleMean();

				//i = 0;
				//Matrix Psioz = (Psi / (m - d - 1)).copy();
				//for (m = d + 2; m < 101 * d + 2; m += 10 * d)
				//{
				//	Psi = Psioz*(m - d - 1);
				//	loglikelihood[i] = 0;
				//	for (auto& cluster : clusters)
				//		loglikelihood[i] += cluster.sampleParams();
				//	for (auto& table : tables)
				//		loglikelihood[i] += table.sampleMean();
				//	loglikelihood[i] += getFullLikelihood(workers) / n;
				//	i++;
				//}

				//m = (sampleFromLog(loglikelihood) * 10 * d + d + 2);
				//Psi = Psioz*(m - d - 1);
				//for (auto& cluster : clusters)
				//	cluster.sampleParams();
				//for (auto& table : tables)
				//	table.sampleMean();
				//cout << m << " " << kappa1 << " " << kappa;
				//Psioz = (Psi / (m - d - 1));
				//Psi = eye(d)*(m - d - 1) * 1;
				//for (i = 0; i < d; i++)
				//{
				//	k = 0;
				//	for (auto s = 1.;s < 11; s += 1)
				//	{
				//		Psi.data[i*d + i] = s*(m - d - 1);
				//		loglikelihood[k] = 0;
				//		for (auto& cluster : clusters)
				//			loglikelihood[k] += cluster.sampleParams();
				//		for (auto& table : tables)
				//			loglikelihood[k] += table.sampleMean();
				//		loglikelihood[k] += getFullLikelihood(workers) / n;
				//		k++;
				//	}
				//	Psi.data[i*d + i] = ((sampleFromLog(loglikelihood) + 1) * 1)*(m - d - 1);
				//	for (auto& cluster : clusters)
				//		cluster.sampleParams();
				//	for (auto& table : tables)
				//		table.sampleMean();
				//}

			}





		reid(clusters);
		reid(tables);
		for (int i = 0; i < n; i++)
			c[i] = z[i]->cls;

		if (((MAX_SWEEP - iter - 1) % STEP) == 0 && iter >= BURNIN)
		{
			int sampleno = (MAX_SWEEP - iter - 1) / STEP;
			if (sampleno<SAMPLE)
			for (int i = 0; i < n; i++)
			{
				zi(sampleno)[i] = z[i]->id;
				superlabels(sampleno)[i] = c[i]->id;
			}
		}
		// Sample Tables
		for (auto& cluster : clusters)
			if (cluster.n > 0)
				cluster.sampleTables(tables);

		//for (i = 0; i < n; i++)
		//{
		//	u[i] = z[i]->beta * urand();
		//	if (z[i]->cls->ustar < u[i])
		//		z[i]->cls->ustar = u[i];
		//}

		//Psi = eye(d);
		//for (auto& cluster : clusters)
		//	Psi += cluster.scatter;
		//Psi /= (n/(2*(m-d-1)));

	}
	return zi;
}




PILL_DEBUG
int main(int argc,char** argv)
{

	Matrix initialLabels;
	generator.seed(time(NULL));
	srand(time(NULL));

	CBLAS = 0; // Do not use vector library in small dimensions
	if (argc < 1)
	{
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (m,kappa,kappai,alpha,gam)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
		return -1;
	}
	nthd = thread::hardware_concurrency();
	DataSet ds(argc,argv);
	Matrix& x = ds.data;
	cout << "NPOINTS :" << x.r << " NDIMS:" << x.m << endl;
	n = x.r; // Number of Points
	d = x.m;
	//init_buffer(nthd, x.m);
	cout << " Available number of threads : " << nthd << endl;
	precomputegamLn(2 * n + 100 * d);
	
	// Hyper-parameters with default values
	if (x.data == NULL)
	{
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (m,kappa,kappai,alpha,gam)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
		return -1;
	}
	cout << m << " " << kappa << " " << kappa1 << " " << alpha << " " << gam << endl;

	if (argc>4)
		MAX_SWEEP = atoi(argv[4]);
	if (argc>5)
		BURNIN = atoi(argv[5]);
	if (argc > 6)
		result_dir = argv[6];
	else
	{
		string str(argv[1]);
		result_dir = (char*)str.substr(0, str.find_last_of("/\\")).c_str(); // Datafile folder
	}
	if (argc > 7)
	{
		SAMPLE = atoi(argv[7]);
	}
	STEP = (MAX_SWEEP - BURNIN) / SAMPLE;
	if (BURNIN >= MAX_SWEEP | STEP == 0) // Housekeeping
	{
		BURNIN = MAX_SWEEP - 2;
		SAMPLE = 1; STEP = 1;
	}
	printf("Reading...\n");
	kep = kappa*kappa1 / (kappa + kappa1);
	eta = m - d + 2;

	ThreadPool tpool(nthd);
	Matrix superlabels(SAMPLE, n);
	cout << "Starting sampling ...";
	auto labels = SliceSampler(x, tpool, superlabels); // data,m,kappa,gam,mean,cov 
	string filename = result_dir;
	labels.writeBin(filename.append("Sublabels.matrix").c_str());
	filename = result_dir;
	superlabels.writeBin(filename.append("Labels.matrix").c_str());

}
