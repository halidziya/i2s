#include "FastMat.h"


class Tester : public Task
{
public:
	Matrix& a;
	Matrix& b;
	Tester(Matrix& a, Matrix &b) : a(a), b(b)
	{
	}
	void run(int id)
	{
		a*b;
	}

};
int main()
{
	
	PILL_DEBUG;
	debugMode(true);
	
	Matrix m(1000, 1000);
	Matrix l(1000, 1000);
	step();
	Tester t(m,l);

	ThreadPool tp(12);
	init_buffer(12, 1000);
	for (int i = 0; i < 24;i++)
	{
		tp.submit(t);
	}
	tp.waitAll();
	step();
	pause();
}