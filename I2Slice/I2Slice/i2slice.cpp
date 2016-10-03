#include <iostream>
#include "FastMat.h"
#include "SGlobal.h"
#include "CompTask.h"
#include <string>
#include <map>


int MAX_SWEEP = 500;
int NINITIAL = 1;
int BURNIN = 300;
int STEP = (MAX_SWEEP - BURNIN) / 10; // Default value is 10 sample + 1 post burnin



class Collector : public Task
{
public:
	Matrix& x;
	Vector& labels;
	Vector count;
	vector<Matrix> scatter;
	vector<Vector> sum;
	atomic<int> taskid;

	Collector(Matrix& x, Vector& labels) : x(x), labels(labels) {
	}

	void reset() {
		int vsize = labels.maximum() + 1;
		count = zeros(vsize);
		scatter = vector<Matrix>(vsize, zeros(d, d));
		sum = vector<Vector>(vsize, zeros(d));
		taskid = 0;
	}

	void run(int id) {
		// Use thread specific buffer
		SETUP_ID()
		int taskid = this->taskid++;
		int maxlab = labels.maximum();
		int label;
		for (int i = 0; i < n; i++)
		{
			label = labels(i);
			if (label == taskid) // Distributed according to label id
			{
				count[label] += 1;
				sum[label] += x(i);
				scatter[label] += x(i) >> x(i); // Second moment actually
			}
		}

		scatter[taskid] -= ((sum[taskid] >> sum[taskid]) / count[taskid]); // Actual Scatter
	}
};


vector<int> relabel(Vector& labels,Vector& superlabels)
{
	for (int i = 0; i < labels.n; i++)
	{
		labels[i] = labels[i] + 1. / superlabels.n; // Dirty way, maybe I can create hierarhical unique labels later
	}
	Vector ulabels = labels.unique();
	vector<int> parents(ulabels.n);
	map<double, int> dict;
	int j = 0;
	for (int i = 0; i < ulabels.n; i++)
	{
		dict[ulabels(i)] = j++;
	}

	for (int i = 0; i < labels.n; i++)
	{
		labels[i] = dict[labels(i)];
		parents[labels[i]] = superlabels[i];
	}
	return parents;
}



Vector SliceSampler(Matrix& x, ThreadPool& workers)
{
	// Point level variables
	int NTABLE = NINITIAL;
	Vector u = ones(n);
	Vector c = zeros(n); // Cluster labels
	Vector beta = ones(NINITIAL);
	Vector z = zeros(n); // Component labels
	vector<Restaurant> clusters = vector<Restaurant>(1);
	vector<int> parents; // Parents of components
	clusters[0].setDist(mu0, eye(d)/100);
	clusters[0].sampleTables(0.1);
	u = urand(n);
	u *= clusters[0].beta.minimum();

	CompTask cmsampler(x,z,u,c,clusters);
	Collector  collector(x,z);
	
	for (int iter = 0; iter < 100; iter++) {
		cmsampler.reset(2 * nthd);
		for (auto i = 0; i < cmsampler.nchunks; i++) {
			workers.submit(cmsampler);
		}
		workers.waitAll();

		// relabel, remove empty ones
		parents = relabel(z,c);
		NTABLE = parents.size();
		collector.reset();
		for (auto i = 0; i < NTABLE; i++) {
			workers.submit(collector);
		}
		workers.waitAll();

		for (int i = 0; i < clusters.size(); i++)
			clusters[i].reset();
		for (int i = 0; i < NTABLE; i++)
		{
			Restaurant* cls = &clusters[0];
			cls->addTable(Table(cls, collector.count[i], collector.sum[i], collector.scatter[i]));
		}
		// Sample Tables
		for (int i = 0; i < clusters.size(); i++)
			clusters[i].sampleTables(u.minimum());

		// Auxilary variables
		for (int i = 0; i < n; i++)
			u[i] = clusters[c[i]].beta[z[i]] * urand();


		// Cluster level beta variable
		Vector valpha = zeros(clusters.size() + 1);
		int i = 0;
		for (auto ti = clusters.begin(); i < clusters.size(); i++, ti++)
			valpha[i] = ti->tables.size();
		valpha[i] = gamma;
		Dirichlet dr(valpha);
		beta = dr.rnd();


	}
	/*
	// Create Betas
	Vector alpha = collector.count.append(gamma);
	Dirichlet dr(alpha);
	beta = dr.rnd();

	// Sample U
	u = rand(n);
	u *= beta[z];
	//New Sticks
	double ustar = u.minimum();
	Vector newsticks = stickBreaker(ustar, beta[beta.n - 1], gamma);

	// Distribute B by assigning to its table
	// Add new tables, zero scatter
	// 
	*/
	return z;
}


int main(int argc,char** argv)
{

	Matrix x;
	Matrix hyperparams;
	Matrix initialLabels;
	generator.seed(time(NULL));
	srand(time(NULL));
	system("dir");
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
	nthd =  thread::hardware_concurrency();
	n = x.r; // Number of Points
	d = x.m;
	init_buffer(nthd, x.m);
	cout << " Available number of threads : " << nthd << endl;
	precomputeGammaLn(2 * n + 100 * d);


	// Hyper-parameters with default values
	if (x.data == NULL)
	{
		cout << "Usage: " << "dpsl.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (d,m,kappa,gamma)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
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
		kappa1 = hyperparams.data[2];
		gamma = hyperparams.data[3];
		cout << m << " " << kappa0 << " "  << kappa1 << " " << gamma << endl;
	}
	else
	{
		m = x.m + 3;
		kappa0 = 1;
		kappa1 = 1;
		gamma = 1;
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
		STEP = (MAX_SWEEP - BURNIN) / atoi(argv[7]);


	ThreadPool tpool(nthd);
	auto labels = SliceSampler(x,tpool); // data,m,kappa,gamma,mean,cov 
	string filename = argv[1];
	labels.writeBin(filename.append(".labels").c_str());
}