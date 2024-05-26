#pragma once
#include <vector>

struct Checker
{
	std::vector<double> Ek;
	size_t N, i, itotal;
	double eps = 1e-10;
	Checker(size_t n = 3)
	{
		N = n;
		itotal = 0;
		i = 0;
		Ek.resize(N);
		for (size_t k = 0; k < N; k++)
		{
			Ek[k] = 0;
		}
	}

	bool stop(double f, bool print = false)
	{
		itotal++;
		i = itotal % N;

		Ek[i] = f;
		if (itotal < N) return false;

		double av = 0.0;
		for (size_t k = 0; k < N; k++)
		{
			av += Ek[k];
		} av = av / N;

		if (print)
		{
			std::cout << "Ek = " << Ek[i] << ", dEk = " << Ek[i] - Ek[(itotal - 1) % N];
			std::cout << ", converged: " << 100 * tanh(eps / abs(av - Ek[i])) << " %";
			std::cout << std::endl;
		}

		if (abs(av - Ek[i]) < eps) return true;
		else return false;
	}
};