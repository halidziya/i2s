#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "FastMat.h"
#include "Restaurant.h"
#include <thread>  
#include <iostream>
#include <algorithm>
#include "Algorithms.h"


using namespace std;


int MAX_SWEEP=1500;
int BURNIN=1400;
int NSAMPLE = 50;
int STEP=(MAX_SWEEP-BURNIN)/NSAMPLE; // Default value is 10 sample + 1 post burnin
char* result_dir = "./";
int NINITIAL = 5; // Should be smaller than d in current matrix implemetation
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
	u = ones(n);
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
		u[i] = c[i]->beta[0]*urand();
		z[i] = c[i]->tables[0];
	}

	Matrix zi((MAX_SWEEP - BURNIN) / STEP + 1, n);
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

	for (int iter = 0; iter <= MAX_SWEEP; iter++) {

		reid(clusters);
		cmsampler.reset(2 * nthd);
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
				newdishprob =  log(gam) + atable.loglik0;
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

		//for (auto& cluster : clusters)
		//	cluster.sampleParams();






			if (false){//iter % 5000 == 1) {
				Vector loglikelihood(50);
				i = 0;
				for (kappa = 0.005; kappa < 0.505; kappa += 0.01)
				{
					loglikelihood[i] = 0;
					kappa1 = kappa * 10;
					for (auto& cluster : clusters)
						loglikelihood[i] += cluster.sampleParams();
					//for (auto& table : tables)
					//	loglikelihood[i] += table.sampleMean()/tables.size();
					//loglikelihood[i] += getFullLikelihood(workers)/n;
					i++;
				}
				//loglikelihood.print();

				kappa = (sampleFromLog(loglikelihood)*0.01 + 0.005);
				kappa1 = kappa * 10;
				for (auto& cluster : clusters)
					cluster.sampleParams();
				for (auto& table : tables)
					table.sampleMean();


				i = 0;
				Matrix Psioz = (Psi / (m - d - 1)).copy();
				for (m = d + 2; m < 101 * d + 2; m += 2 * d)
				{
					Psi = Psioz*(m - d - 1);
					loglikelihood[i] = 0;
					for (auto& cluster : clusters)
						loglikelihood[i] += cluster.sampleParams();
					for (auto& table : tables)
						loglikelihood[i] += table.sampleMean();
					loglikelihood[i] += getFullLikelihood(workers) / n;
					i++;
				}


				m = (sampleFromLog(loglikelihood) * 2 * d + d + 2);
				Psi = Psioz*(m - d - 1);
				for (auto& cluster : clusters)
					cluster.sampleParams();
				for (auto& table : tables)
					table.sampleMean();

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
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (d,m,kappa,kappa1,gam)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
		return -1;
	}
	cout << "NPOINTS :" << x.r << " NDIMS:" << x.m << endl;
	nthd = thread::hardware_concurrency();
	n = x.r; // Number of Points
	d = x.m;

	init_buffer(nthd, x.m);
	cout << " Available number of threads : " << nthd << endl;
	precomputegamLn(2 * n + 100 * d);


	// Hyper-parameters with default values
	if (x.data == NULL)
	{
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (d,m,kappa,gam)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
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
		kappa = hyperparams.data[2];
		kappa1 = hyperparams.data[3];
		alpha = hyperparams.data[4];
		gam = hyperparams.data[5];
		cout << m << " " << kappa << " " << kappa1 << " " << alpha << " " << gam << endl;
	}
	else
	{
		m = x.m + 3;
		kappa = 0.05;
		kappa1 = 0.5;
		gam = 1;
		alpha = 1;
	}

	if (argc > 3)
		Psi.readBin(argv[3]);
	else
		Psi = (eye(d)*2*(m-d-1)).copy();

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
	kep = kappa*kappa1 / (kappa + kappa1);
	eta = m - d + 2;
	precomputegamLn(2 * (n + d) + 1);  // With 0.5 increments
	init_buffer(thread::hardware_concurrency(), d);

	ThreadPool tpool(nthd);
	Matrix superlabels((MAX_SWEEP - BURNIN) / STEP + 1, n);
	cout << "Starting sampling ...";
	auto labels = SliceSampler(x, tpool, superlabels); // data,m,kappa,gam,mean,cov 
	string filename = argv[1];
	labels.writeBin(filename.append(".labels").c_str());
	filename = argv[1];
	superlabels.writeBin(filename.append(".superlabels").c_str());

}
