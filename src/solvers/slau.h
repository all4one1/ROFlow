#pragma once
#include "../FromOuterSparse/SparseMatrix.h"
#include <iostream>


inline void solveGS(double* f, double* f0, double* bb, int NN, SparseMatrix& M)
{
	unsigned int k = 0;
	double eps_iter = 1e-6;

	for (k = 1; k < 100000; k++)
	{
		double s = 0;
		for (int j = 0; j < NN; j++)
		{
			s = M.line(j, f);
			f[j] = f[j] + (bb[j] - s) / M[j][j];
		}


		double max = 0;
		double dif;
		for (int i = 0; i < NN; i++)
		{
			dif = abs(f0[i] - f[i]);
			if (dif > max)
				max = dif;
		}
		for (int j = 0; j < NN; j++)
			f0[j] = f[j];

		if (max < eps_iter)	break;

		if (k % 10000 == 0) std::cout << "stream_cpu k = " << k << ", eps = " << max << std::endl;
	}

}
