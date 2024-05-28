#pragma once
#include <vector>
#include <map>
#include <sstream>
#include <string>
struct Checker
{
	std::vector<double> Ek;
	size_t N, i, itotal;
	double eps = 1e-15;
	double dEk = 0;
	Checker(size_t n = 10)
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

		dEk = Ek[i] - Ek[(itotal - 1) % N];

		//double res = tanh(eps / abs(av - Ek[i]));

		double d = (abs(dEk) / Ek[i]);
		double res = tanh(eps / d);
		if (print)
		{
			std::cout << "Ek = " << Ek[i] << ", dEk = " << dEk;
			std::cout << ", converged: " << 100 * res << " %";
			std::cout << std::endl;
		}

		if (res > 0.99999) return true;
		else return false;
	}
};


struct StateOut
{
	std::map<std::string, double> m;
	std::map<std::string, double*> mp;
	std::ofstream file;
	std::string defname = "stats.dat";
	size_t iter = 0;
	StateOut()
	{
		file.open(defname, std::ios::app);
	}

	void write_head(std::string h)
	{
		if (iter == 0)
			file << h << std::endl;
		iter++;
	}

	void write(std::vector<double> v)
	{
		for (auto& it : v) 
			file << it << " "; 
		file << std::endl; 
		iter++;
	}

	void write()
	{
		iter++;
		if (iter == 1)
		{
			for (auto& it : mp)	file << it.first << " ";
			file << std::endl;
		}

		for (auto& it : mp)	file << *(it.second) << " ";
		file << std::endl;
	}
	void add(std::string name, double &f)
	{
		mp[name] = &f;
	}

	//std::string int_to_str(int i)
	//{
	//	std::stringstream ss;
	//	std::string name;
	//	ss.str("");
	//	ss.clear();
	//	ss << i;
	//	name = ss.str();
	//	return name;
	//}
};