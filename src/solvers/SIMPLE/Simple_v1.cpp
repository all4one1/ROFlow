#include "Simple_v1.h"

SIMPLE::SIMPLE(const Configuration& config)
{
	if (dim == 3)
	{
		Nvx = (nx + 1) * ny * nz;
		Nvy = nx * (ny + 1) * nz;
		Nvz = nx * ny * (nz + 1);
	}
	if (dim == 2)
	{
		Nvx = (nx + 1) * ny;
		Nvy = nx * (ny + 1);
		Nvz = 0;
	}
	if (dim == 1)
	{
		Nvx = (nx + 1);
		Nvy = 0;
		Nvz = 0;
	}

	NV = Nvx + Nvy + Nvz;
	stride = Nvx;
	stride2 = Nvx + Nvy;

	auto alloc = [](double** f, unsigned int N)
	{
		*f = new double[N];
		memset(*f, 0, sizeof(double) * N);
		//return sizeof(double) * N;
	};

	mk = MatrixMaker(*this);

	T = Variable(*this, true);
	T0 = Variable(*this, true);
	
}


