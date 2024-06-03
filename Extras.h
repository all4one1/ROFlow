#pragma once
#include <vector>
#include <map>
#include <sstream>
#include <string>
#ifdef _WIN32
#include "Windows.h"
#endif // 

struct Checker
{
	std::vector<double> Ek;
	size_t N, i, itotal;
	double eps = 1e-7;
	double dEk = 0;
	double dif = 0, res = 0, long_dif = 0;
	Checker(size_t n = 20)
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
		if (itotal < 2) return false;

		dEk = Ek[i] - Ek[(itotal - 1) % N];
		dif = (abs(dEk) / Ek[i]);
		res = tanh(eps / dif);
		if (print)
		{
			std::cout << "Ek = " << Ek[i] << ", dEk = " << dEk;
			std::cout << ", converged: " << 100 * res << " %";
			std::cout << std::endl;
		}
		if (res > 0.99999) return true;



		if (itotal > N*2)
		{
			double av = 0.0;
			for (size_t k = 0; k < N; k++)
			{
				av += Ek[k];
			} av = av / N;
			//long_dif = tanh(eps / abs(av - Ek[i]));
			long_dif = abs(av - Ek[i]);

			if (long_dif < 1e-7 && av > 1e-5) return true;
		}
		
		return false;
	}
};


struct StateOut
{
	std::map<std::string, double> m;
	std::map<std::string, double*> mp;
	std::ofstream file;
	size_t iter = 0;
	bool header = false;
	StateOut()
	{
		create_folder("results");
	}

	void open(bool append = true, std::string path = "results\\stats.dat")
	{
		if (append)
		{
			file.open(path, std::ios_base::app);
			if (file_empty(path))
			{
				header = false;
			}
			else
			{
				header = true;
				file << std::endl;
			}
		}
		else
		{
			file.open(path);
			header = false;
		}	
	}

	void create_folder(std::string name)
	{
		#ifdef __linux__
		std::string str = "mkdir -p " + name + "/";
		system(str.c_str());
		#endif 

		#ifdef _WIN32
		CreateDirectoryA(name.c_str(), NULL);
		#endif 
	}

	void write_header(std::string h, bool forced = false)
	{
		if (header == false || forced == true)
		{
			file << h << std::endl;
			header = true;
		}
	}

	void write(std::vector<double> v)
	{
		for (auto& it : v) 
			file << it << " "; 
		file << std::endl; 
		iter++;
	}

	bool file_empty(std::string s)
	{
		std::ifstream pFile(s);
		return pFile.peek() == std::ifstream::traits_type::eof();
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


struct FuncTimer
{
private:
	std::map <std::string, double> timer, t1, t2;
	std::string active;
	double dt, total;
public:
	FuncTimer()
	{
		dt = total = 0;
	}
	void start(std::string s)
	{
		active = s;
		t1[s] = clock();
	}
	void end(std::string s)
	{
		if (t1.find(s) == t1.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return;
		}
		else
		{
			t2[s] = clock();
			dt = (t2[s] - t1[s]) / CLOCKS_PER_SEC;
			//if (s == active)
			{
				timer[s] += dt;
			}
		}
	}

	double get(std::string s)
	{
		if (timer.find(s) == timer.end())
		{
			std::cout << s + " trigger not started" << std::endl;
			return 0.0;
		}
		else
		{
			return timer[s];
		}
	}

	std::string get_info()
	{
		int n = timer.size();
		std::ostringstream oss;
		oss << "Calculation time in seconds. Number of cases: " << n << ".\n";

		for (auto& it : timer)
		{
			oss << it.first << ": " << it.second <<  std::endl;
		}
		return oss.str();
	}

	void show_info()
	{
		std::cout << get_info() << std::endl;
	}
	
	void write_info(std::string path = "results\\report.dat")
	{
		std::ofstream file(path);
		file << get_info() << std::endl;
	}

};