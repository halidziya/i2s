#include <iostream>
#include "FastMat.h"
#include "Restaurant.h"
#include "SGlobal.h"
#include <string>
#include <map>


int MAX_SWEEP = 1800;
int NINITIAL = 1;
int BURNIN = 1600;
int NSAMPLE = 200;
int STEP = (MAX_SWEEP - BURNIN) / NSAMPLE; // Default value is 10 sample + 1 post burnin


Vector u0;
Vector u1;
vector<Restaurant*> c;
vector<Table*> z;
Matrix x;
list<Table> tables;
Vector beta;
list<Restaurant> clusters;

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
				int NTABLE = tables.size();
				Vector likelihoods(NTABLE);
				int j = 0;
				for (auto& t : tables)
				{
					if ((u0[i] <= beta[t->cls->id]) & (t->beta >= u1[i]))
					{
						likelihoods[j] = t->dist.likelihood(x(i)); //**
					}
					else
					{
						likelihoods[j] = -INFINITY;
					}
					j++;
				}
				z[i] = tables[sampleFromLog(likelihoods)]; //**

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
			}
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
	// Point level variables
	x = data;
	int NTABLE = NINITIAL;
	int k = 0, i = 0, j = 0;
	u0 = ones(n);
	u1 = ones(n);
	beta = ones(NINITIAL);
	beta /= NINITIAL;
	z = vector<Table*>(n); // Component labels
	c = vector<Restaurant*>(n);
	clusters = list<Restaurant>(NINITIAL);
	vector<Restaurant*> cp = vector<Restaurant*>(clusters.size());
	Stut priorcollapsed = Stut(mu0, Psi * ((kappa0 + 1) / (kappa0*(m - d + 1))), m-d+1);
	for (auto& cluster : clusters)
	{
		cluster.sampleParams();
		cluster.sampleTables(tables, 0.005);
		cp[i++] = &cluster;
	}
	for (int i = 0; i < n; i++)
	{
		j = rand() % NINITIAL;
		c[i] = cp[j];
		u0[i] = beta[j] * urand();
		j = rand() % c[i]->tables.size();
		z[i] = c[i]->tables[j];
		u1[i] = c[i]->beta[j] * urand();
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
		cout << iter << endl;
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

		vector<Restaurant*> cp = vector<Restaurant*>(clusters.size());
		vector<int> cidx = vector<int>(tables.size());
		i = 0;
		for (auto& cluster : clusters)
		{
			cp[i] = &cluster;
			i++;
		}

		i = 0;
		for (auto& cluster : clusters)
			cluster.resetStats();
		for (auto& atable : tables)
			atable.cls->addStats(&atable);
		for (auto& cluster : clusters)
			cluster.sampleParams();
		for (auto& atable : tables)
		{
			j = 0;
			Vector likelihood(clusters.size() );
			atable.cls->remStats(&atable);
			for (auto& cluster : clusters)
			{

				likelihood[j] = cluster.dist.likelihood(atable.sum / atable.n) + log(cluster.n);
				//Vector dd = absbuffer.get();
				//dd <<= atable.dist.mu - cluster.dist.mu;
				//cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, atable.dist.cholsigma.r, atable.dist.cholsigma.data, atable.dist.cholsigma.r, dd.data, 1);
				//double dist = cblas_dnrm2(dd.n, dd.data, 1);
				//dist = dist*dist;
				//likelihood[j] = -0.5 * (((atable.cls->sigma).inverse()*cluster.sigma).diag().sum() + dist - d + atable.dist.cholsigma.sumlogdiag() - cluster.dist.cholsigma.sumlogdiag());
                //likelihood[j] += log(cluster.n) ;
				//likelihood[j] += Wishart((cluster.scatter+Psi)/(cluster.n + m), cluster.n + m).likelihood(atable.scatter + Psi);
				j++;
			}
			//Vector dd = absbuffer.get();
			//dd <<= atable.dist.mu - mu0;
			//cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, atable.dist.cholsigma.r, atable.dist.cholsigma.data, atable.dist.cholsigma.r, dd.data, 1);
			//double dist = cblas_dnrm2(dd.n, dd.data, 1);
			//dist = dist*dist;
			//Matrix asigmaj = priorcov.rnd();


			//likelihood[j] = -0.5 * (((atable.cls->sigma).inverse()*asigmaj).diag().sum() + dist - d + atable.dist.cholsigma.sumlogdiag() - asigmaj.chol().sumlogdiag());
			//likelihood[j] +=  log(gam);// +Wishart(Psi, m).likelihood(atable.scatter + Psi);

			cidx[i] = sampleFromLog(likelihood);
			//if (cidx[i] >= cp.size()) {
			//	clusters.push_back(Restaurant());
			//	cp.push_back(&clusters.back());
			//}

			atable.cls = cp[cidx[i]];
			atable.cls->addStats(&atable);
			atable.cls->sampleParams();
			i++;
		}

		for (auto& cluster : clusters)
			cluster.reset();

		for (auto& atable : tables)
		{
			atable.cls->addTable(&atable);
		}
		// Remove empty clusters
		for (auto cc = clusters.begin(); cc != clusters.end();)
		{

			if (cc->n == 0)
				cc = clusters.erase(cc);
			else
				cc++;
		}

		for (auto cc = clusters.begin(); cc != clusters.end(); cc++)
			cout << " " << cc->n;
		cout << endl;


		// Cluster level beta variable
		Vector valpha = zeros(clusters.size() + 1);
		int i = 0;
		for (auto ti = clusters.begin(); i < clusters.size(); i++, ti++)
			valpha[i] = ti->n; // ti->tables.size();
		valpha[i] = gam;
		Dirichlet dr(valpha);
		beta = dr.rnd();

		Vector	newsticks = stickBreaker(u0.minimum(), beta[beta.n - 1], alpha);
		beta.resize(beta.n - 1);
		beta = beta.append(newsticks);

		for (int i = clusters.size(); i < beta.n; i++)
			clusters.push_back(Restaurant());

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
				cluster.sampleTables(tables, u1.minimum());

		for (i = 0; i < n; i++)
			u1[i] = z[i]->beta * urand();

		for (i = 0; i < n; i++)
			u0[i] = beta[c[i]->id] * urand();

	}
	return zi;
}



int main(int argc, char** argv)
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
		cout << "Usage: " << "i2slice.exe datafile.matrix [hypermean.matrix] [hyperscatter.matrix] [params.matrix (d,m,kappa0,kappa1,gam)]  [#ITERATION] [#BURNIN] [#SAMPLE]  [initiallabels.matrix]: In fixed order";
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
		kappa0 = hyperparams.data[2];
		kappa1 = hyperparams.data[2];
		gam = hyperparams.data[3];
		cout << m << " " << kappa0 << " " << kappa1 << " " << gam << " " << alpha << endl;
	}
	else
	{
		m = x.m + 3;
		kappa0 = 1;
		kappa1 = 1;
		kappa2 = 1;
		gam = 1;
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


	ThreadPool tpool(nthd);
	Matrix superlabels((MAX_SWEEP - BURNIN) / STEP + 1, n);
	priorcov = IWishart(Psi, m);
	priormean = Normal(mu0, priorcov.rnd() / kappa0);
	auto labels = SliceSampler(x, tpool, superlabels); // data,m,kappa,gam,mean,cov 
	string filename = argv[1];
	labels.writeBin(filename.append(".labels").c_str());
	filename = argv[1];
	superlabels.writeBin(filename.append(".superlabels").c_str());
}
