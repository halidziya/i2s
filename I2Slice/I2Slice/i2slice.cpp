#include <iostream>
#include "FastMat.h"
#include "SGlobal.h"
#include "CompTask.h"
#include <string>
#include <map>


int MAX_SWEEP = 1700;
int NINITIAL = 4;
int BURNIN = 1600;
int NSAMPLE = 100;
int STEP = (MAX_SWEEP - BURNIN)/ NSAMPLE; // Default value is 10 sample + 1 post burnin



class Collector : public Task
{
public:
	Matrix& x;
	vector<Table*>& labels;
	atomic<int> taskid;
	int ntable;

	Collector(Matrix& x, vector<Table*>& labels,int ntable) : x(x), labels(labels) , ntable(ntable) {
	}

	void reset(list<Table>& tables) {
		for (auto& atable : tables)
			atable.reset();
		taskid = 0;
	}

	void subscatter(list<Table>& tables) {
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
			if (labels[i]->id == taskid) // Distributed according to label id
			{
				labels[i]->n += 1;
				labels[i]->sum += x(i);
				labels[i]->scatter += x(i) >> x(i); // Second moment actually
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


Matrix SliceSampler(Matrix& x, ThreadPool& workers,Matrix& superlabels)
{

	
	// Point level variables
	int NTABLE = NINITIAL;
	int k = 0, i = 0, j = 0;
	Vector u = ones(n);
	Vector beta = ones(NINITIAL);
	beta /= NINITIAL;
	vector<Table*> z = vector<Table*>(n); // Component labels
	vector<Restaurant*> c = vector<Restaurant*>(n);
	list<Table> tables; // Component labels
	list<Restaurant> clusters = list<Restaurant>(NINITIAL);
	vector<Restaurant*> cp = vector<Restaurant*>(clusters.size());
	vector<int> parents; // Parents of components
	

	for (auto& cluster : clusters)
	{
		cluster.sampleParams();
		cluster.sampleTables(tables, 0.0005);
		cp[i++] = &cluster;
	}
	for (int i = 0; i < n; i++)
	{
		c[i] = cp[rand()%NINITIAL];
		j = rand() % c[i]->tables.size();
		z[i] = c[i]->tables[j];
		u[i] = c[i]->beta[j] * urand();
	}

	Matrix zi((MAX_SWEEP - BURNIN) / STEP + 1,n);
	CompTask cmsampler(x,z,u,c);
	Collector  collector(x,z,tables.size());

	for (int iter = 0; iter <= MAX_SWEEP; iter++) {
		cmsampler.reset(2 * nthd);
		for (auto i = 0; i < cmsampler.nchunks; i++) {
			workers.submit(cmsampler);
		}
		workers.waitAll();
		cout << iter << endl;
		reid(tables);
		// Auxilary variables

		// relabel, remove empty ones
		//parents = relabel(z,c);
		//NTABLE = parents.size();
		NTABLE = tables.size();
		collector.reset(tables);
		for (i = 0; i < NTABLE; i++) {
			workers.submit(collector);
		}
		workers.waitAll();
		collector.subscatter(tables);

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


		for (auto& table : tables)
			table.sampleMean();


		i = 0;
		for (auto& cluster : clusters)
			cluster.id = i++;
		// Upper Layer
		// Auxilary variables


		double ustar = 1;
		for (auto& cluster : clusters)
			for (auto& atable : cluster.tables) {
				atable->u = beta[atable->cls->id] * urand();
				k = k + 1;
				if ((atable->n) >0 &&  (atable->u < ustar))
					ustar = atable->u;
			}
		// Sample tables

		Vector likelihood(clusters.size());
		vector<Restaurant*> cp = vector<Restaurant*>(clusters.size());
		vector<int> cidx = vector<int>(tables.size());
		i = 0;
		for (auto& cluster : clusters)
		{
			cp[i] = &cluster;
			i++;
		}

		i = 0;
//		int maxi;
//		double maxli = -INFINITY;
		for (auto& atable : tables)
		{
			j = 0;
			//maxli = -INFINITY;
			for (auto& cluster : clusters)
			{
				if (atable.u <= beta[j])
				{
					likelihood[j] = cluster.dist.likelihood(atable.dist.mu);
					likelihood[j] += Wishart(cluster.sigma, atable.n+m).likelihood(atable.scatter+ Psi);
					//for (auto& btable : cluster.tables)
					//{
					//	likelihood[j] += btable->dist.likelihood(atable.dist.mu);
					//	//likelihood[j] += Wishart((btable->scatter+ Psi)/(btable->n + m), atable.n + m).likelihood(atable.scatter + Psi);

					//}
					//likelihood[j] /= (cluster.tables.size()+1);
				}
				else
				{
					likelihood[j] = -INFINITY;
				}
				//if (likelihood[j] > maxli)
				//{
				//	maxli = likelihood[j];
				//	maxi = j;
				//}
				j++;
			}
			//if (iter > 500)
			//	cidx[i] = maxi;//sampleFromLog(likelihood);
			//else
				cidx[i] = sampleFromLog(likelihood);
			atable.cls = cp[cidx[i]];
			
			i++;
		}


		//for (auto& cluster : clusters)
		//	cluster.reset();

		//for (auto& atable : tables)
		//{
		//	atable.cls->addStats(&atable);
		//}

		////Metropolis Hasting Step
		//i = 0;
		//for (auto& atable : tables)
		//{
		//	int idx = rand() % cp.size();
		//	Restaurant* randclust = cp[idx];
		//	if (randclust != atable.cls)
		//	{
		//		priormean = Normal(mu0, priorcov.rnd() / kappa0);
		//		if (randclust->n == 0)
		//			randclust->sigma = IWishart(atable.scatter + Psi, atable.n + m).rnd();
		//		Normal dist1 = Normal(Normal((mu0*kappa0 + randclust->sum + atable.dist.mu) / (kappa0 + randclust->nt + 1), randclust->sigma / (kappa0 + randclust->nt + 1)).rnd(), randclust->sigma / kappa1);
		//		double mhratio = dist1.likelihood(atable.dist.mu) + log(beta[idx]) + priormean.likelihood(dist1.mu) + Wishart(randclust->sigma, atable.n + m).likelihood(atable.scatter + Psi)
		//			- atable.cls->dist.likelihood(atable.dist.mu) - log(beta[cidx[i]]) - priormean.likelihood(randclust->dist.mu) -  Wishart(atable.cls->sigma, atable.n + m).likelihood(atable.scatter + Psi);;
		//		if (rand() < exp(mhratio)) // Accept
		//		{
		//			cidx[i] = idx;
		//			atable.cls =randclust;
		//		}
		//	}
		//	i++;
		//}
		//cout << endl;


		// Move tables in higher level
		for (auto& cluster : clusters)
			cluster.reset();

		for (auto& atable : tables)
		{
			atable.cls->addTable(&atable);
		}




		for (auto cc = clusters.begin(); cc != clusters.end();)
		{
			
			if (cc->n == 0)
				cc = clusters.erase(cc);
			else
				cc++;
		}

		for (auto cc = clusters.begin(); cc != clusters.end();cc++)
		cout << " " << cc->n;
		cout << endl;


		// Cluster level beta variable
		Vector valpha = zeros(clusters.size() + 1);
		int i = 0;
		for (auto ti = clusters.begin(); i < clusters.size(); i++, ti++)
			valpha[i] = ti->n; // ti->tables.size();
		valpha[i] = gamma;
		Dirichlet dr(valpha);
		beta = dr.rnd();

		Vector	newsticks = stickBreaker(ustar, beta[beta.n - 1], alpha);
		beta.resize(beta.n - 1);
		beta = beta.append(newsticks);

		//Remove unused clusters
		for (int i = clusters.size(); i < beta.n; i++)
			clusters.push_back(Restaurant());

		i = 0;
		for (auto& cluster : clusters)
		{
			cluster.sampleParams();
			cluster.id = i++;
		}

		for (int i = 0; i < n; i++)
		{
			c[i] = z[i]->cls;
		}

		if (iter > BURNIN && ((iter - BURNIN) % STEP == 0))
			for (int i = 0; i < n; i++)
			{
				zi((iter - BURNIN) / STEP)[i] = z[i]->id;
				superlabels((iter - BURNIN) / STEP)[i] = c[i]->id;
			}

		// Sample Tables
		for (auto& cluster : clusters)
			if (cluster.n > 0)
				cluster.sampleTables(tables,u.minimum());
		for (auto& cluster : clusters) {
			i = 0;
			for (auto& table : cluster.tables)
				table->id = i++;
		}

		for (auto& table : tables)
		{
			table.dist = Normal(table.dist.mu, table.cls->sigma);
			cout << table.n << " ";
		}
		cout << endl;

		for (i = 0; i < n; i++)
			u[i] = c[i]->beta[z[i]->id] * urand();


	}	
	return zi;
}


int imain(int argc,char** argv)
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
		cout << m << " " << kappa0 << " "  << kappa1 << " " << gamma << " " << alpha << endl;
	}
	else
	{
		m = x.m + 2;
		kappa0 = .05;
		kappa1 = .5;
		gamma = 1;
		alpha = 1;
	}

	if (argc > 3)
		Psi.readBin(argv[3]);
	else
		Psi = (x.cov()/100).copy();




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
	priorcov  = IWishart(Psi,m);
	priormean = Normal(mu0, priorcov.rnd()/kappa0);
	auto labels = SliceSampler(x,tpool,superlabels); // data,m,kappa,gamma,mean,cov 
	string filename = argv[1];
	labels.writeBin(filename.append(".labels").c_str());
	filename = argv[1];
	superlabels.writeBin(filename.append(".superlabels").c_str());
}
