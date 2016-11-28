#include "Restaurant.h"
#include "Table.h"

#include "Restaurant.h"
#include <iostream>


Restaurant::Restaurant()
{
	reset();
}

void Restaurant::reset()
{
	tables.resize(0);
	collapsedtables.resize(0);
	resetStats();
}

void Restaurant::resetStats()
{
	scatter = zeros(d, d);
	sum = zeros(d);
	n = 0;
	nt = 0;
	dist.mu = mu0;
}

void Restaurant::sampleParams()
{
	Matrix otherscatter = zeros(d, d);
	if (tables.size() > 0)
	{
		Vector& diff = (mu0 - (sum / tables.size()));
		otherscatter = (diff >> diff)*(kappa*tables.size() / (kappa + tables.size()));
		for (auto& atable : tables)
		{
			otherscatter += (atable->dist.mu - dist.mu) >> (atable->dist.mu - dist.mu);
		}
	}
	sigma = IWishart(Psi + scatter + otherscatter, n + m + tables.size()).rnd();
	dist = Normal(Normal((mu0*kappa + sum) / (kappa + tables.size()), sigma / (kappa + tables.size())).rnd(), sigma / kappa1);
}


void Restaurant::calculateDist()
{
	double kap = harmean(nt + kappa, kappa1);
	tdist.eta = m + 2 + this->n - d - nt; //dpar0
	Vector& diff = (mu0 - (sum / tables.size()));
	Matrix otherscatter = (diff >> diff)*(kappa*tables.size() / (kappa + tables.size()));
	tdist.mu = (sum + mu0*kappa) / (nt + kappa); // dpar1
	tdist.cholsigma = ((Psi + scatter + otherscatter)*((kap + 1) / (kap*tdist.eta))).chol(); // dpar2// mu_ti is deleted because of it is same with dpar0
	tdist.calculateNormalizer();
}

void  Restaurant::addTable(Table* table)
{
	tables.push_back(table);
	addStats(table);
}

void  Restaurant::addCollapsedTable(Table* table)
{
	collapsedtables.push_back(table);
}

void Restaurant::addStats(Table* table)
{
	scatter += table->scatter;
	sum += table->dist.mu;
	n += table->n;
	nt += 1;
}

void Restaurant::remStats(Table* table)
{
	scatter -= table->scatter;
	sum -= table->dist.mu;
	n -= table->n;
	nt -= 1;

}



void Restaurant::sampleTables(list<Table>& mainlist)
{
	// Create Betas
	Vector valpha = zeros(tables.size() + 1);
	int i = 0;
	for (auto ti = tables.begin(); i < tables.size(); i++, ti++)
		valpha[i] = (*ti)->n;
	valpha[i] = alpha;
	Dirichlet dr(valpha);
	beta = dr.rnd();
	//New Sticks
	//Vector	newsticks = stickBreaker(ustar, beta[beta.n - 1], gamma);
	//beta.resize(beta.n - 1);
	//beta = beta.append(newsticks);
	//Normal priormean(this->dist.mu, this->sigma / kappa1);
	//for (i = 0; i < newsticks.n; i++)
	//{
	//	mainlist.push_back(Table(priormean.rnd(), this->sigma));
	//	mainlist.back().cls = this;
	//	tables.push_back(&mainlist.back());
	//}
	i = 0;
	for (auto& atable : tables)
	{
		atable->beta = beta[i++];
	}
	ustar = 1; //For next iteration
}

void Restaurant::setDist(Vector & mu, Matrix & sigma)
{
	this->sigma = sigma;
	this->dist = Normal(mu, sigma);
}


Restaurant::~Restaurant()
{
}



//
//Restaurant::Restaurant(void)
//{
//}
//
//
//void Restaurant::operator=(Restaurant& r)
//{
//	tables = r.tables;
//	customers = r.customers;
//	Restaurantid = r.Restaurantid;
//	list<Table>::iterator it,tit;
//	vector<Customer>::iterator cit,rit;
//	int i;
//	for(it=tables.begin(),tit=r.tables.begin();tit!=r.tables.end();tit++,it++)
//	{
//		tit->copy = it;
//		it->dishp = tit->dishp->copy;
//	}
//	for(cit=customers.begin(),rit=r.customers.begin();cit!=customers.end();cit++,rit++)
//	{
//		cit->table = rit->table->copy; // Point to copied object
//	}
//
//}
//
//
//Restaurant::Restaurant(Restaurant& r)
//{
//	operator=(r);
//}
//
//Restaurant::~Restaurant(void)
//{
//}
//
//
//void Restaurant::addTable(Table& t)
//{
//	tables.push_back(t);
//}
//
//void Restaurant::run(int id)
//{
//	// Use thread specific buffer
//	SETUP_ID();
//
//	int i,n,ti,npts;
//	double sum,val,max,subrest,newclust; // Likelihood sum for normalization
//	n = customers.size();
//	list<Table>::iterator tit,oldtable; // Table iterator
//	list<Dish>::iterator ki;
//	Table copy;
//	likelihood=0;
//	for (i=0;i< n;i++)
//	{
//		Vector& x = customers[i].data;
//		oldtable = customers[i].table;
//		copy = *oldtable;
//		ki = customers[i].table->dishp;
//
//		customers[i].table->removePoint(x);
//		// printf("%d\n",customers[i].table->npoints);
//		sum = 0;
//		
//		
//		npts = customers[i].table->npoints;
//
//		if (customers[i].table->npoints ==0 )
//			tables.erase(customers[i].table);
//
//		/*for(tit=tables.begin();tit!=tables.end();tit++)
//			if (ki == tit->dishp)
//				subrest+= tit->npoints;*/
//		
//		newclust = customers[i].loglik0 + log(alpha); ///(subrest+alpha)
//		max = newclust;
//
//		for(tit=tables.begin();tit!=tables.end();tit++)
//		{
//
//				tit->loglikelihood = tit->dist.likelihood(x);
//				tit->logprob = tit->loglikelihood  + log(tit->npoints); ///(subrest+alpha)
//
//			if (tit->logprob > max)
//				max = tit->logprob;
//		}
//
//
//		for(tit=tables.begin();tit!=tables.end();tit++)
//		{
//				tit->logprob = exp(tit->logprob-max); // No longer in logarithm actually
//				sum += tit->logprob;
//		}
//
//		sum += exp(newclust - max); // New class
//
//
//		val = urand()*sum;
//		for(tit=tables.begin();tit!=tables.end();tit++)
//		{
//			if ( (tit->logprob) >= val )
//			{
//				break; // Find it
//			}
//			val -= (tit->logprob);
//		}
//
//		if (tit==tables.end()) // Not in current tables add new one
//		{
//			tables.emplace_front(ki); // New empty table
//			tit = tables.begin();
//			tit->loglikelihood = customers[i].loglik0;
//		}
//
//		customers[i].table = tit;
//
//		likelihood += tit->loglikelihood;
//		if (tit==oldtable)
//			*tit = copy;
//		else
//			tit->addPoint(x); // Add point to selected one
//
//	}
//}
//
//
//ostream& operator<<(ostream& os, Restaurant& r)
//{
//	os.write((char*) &r.Restaurantid,sizeof(int));
//	os.write((char*) &r.likelihood,sizeof(double));
//	int ncustomers = r.customers.size();
//	int i;
//
//	int ntables = r.tables.size();
//	os.write((char*) &ntables,sizeof(int));
//	list<Table>::iterator tit;
//	for(i=0,tit=r.tables.begin();tit!=r.tables.end();tit++,i++)
//	{
//		tit->tableid = i+1;
//		os << *tit;
//	}
//
//	os.write((char*) &ncustomers,sizeof(int));
//	for (i=0;i<ncustomers;i++)
//	{
//	/*	if (r.customers[i].table->tableid>10)
//		{
//			printf("Tables\n");
//			for(tit=r.tables.begin();tit!=r.tables.end();tit++)
//			{
//				printf("%p %d %d\n",tit._Ptr,tit->tableid,tit->npoints);
//			}
//			printf("Customer\n");
//			printf("%p\n",r.customers[i].table._Ptr);
//			printf("%d\n",r.customers[i].table->tableid);
//			printf("%d\n",r.customers[i].table->npoints);
//		}*/
//		os << r.customers[i];
//	}
//
//
//	return os;
//}