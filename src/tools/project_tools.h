#pragma once
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <windows.h>

#include "types.h"
#include "common_tools.h"

using std::string;
using std::ifstream;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::cout;
using std::endl;

//#define __FILENAME__ (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)
//#define print(message) {std::cout << __FILENAME__  << " (Line #" << __LINE__ << "): " << message << std::endl;}
//#define pause {std::cout << "Pause at line: " << __LINE__ << " in file: " << __FILENAME__  << std::endl; system("pause"); }

std::map<string, string> folders;
inline void init_parameters(Configuration& c)
{
	ReadingFile par("parameters.txt");
	//addFolder(folders, { "fields", "final", "lines", "matrices", "pictures" });

	std::string str;
	par.reading_string(str, "domain", "closed_box");
	if (str == "closed_box") c.domain = closed_box;
	if (str == "open_tube") c.domain = open_tube;

	par.reading_string(str, "method", "FV");
	if (str == "FV") c.disc = FV;
	if (str == "FD") c.disc = FD;

	if (c.disc == FV) c.q = 0;
	if (c.disc == FD) c.q = 1;


	int bc;
	par.reading<int>(bc, "xbc", 0); c.xbc = static_cast<bc_type>(bc);
	par.reading<int>(bc, "ybc", 0); c.ybc = static_cast<bc_type>(bc);

	par.reading<double>(c.Pr, "Pr", 10);
	par.reading<double>(c.Ra, "Ra", 1000);

	par.reading<double>(c.Ly, "Ly", 1.0);
	par.reading<double>(c.Lx, "Lx", 1.0);
	par.reading<unsigned int>(c.nx, "nx", 20);
	par.reading<unsigned int>(c.ny, "ny", 20);
	c.nx = (unsigned int)(c.Lx * c.nx);
	c.ny = (unsigned int)(c.Ly * c.ny);
	par.reading<double>(c.grav_y, "grav_y", 1.0);
	par.reading<double>(c.incr_parameter, "increment", -20);


	par.reading<double>(c.Re, "Re", 1);
	par.reading<double>(c.tau, "tau", 3e-7);



	par.reading<double>(c.Cn, "Cn", 1e-4);
	par.reading<double>(c.Pe, "Pe", 10000);


	par.reading<double>(c.p_in, "p_in", 0.0);
	par.reading<double>(c.p_out, "p_out", 0.0);


	if (c.nx > 0) c.dim = 1;
	if (c.ny > 0) c.dim = 2;
	if (c.nz > 0) c.dim = 3;

	c.hx = c.Lx / c.nx;
	c.hy = c.Ly / c.ny;
	c.hz = c.Lz / c.nz;



	if (c.dim == 1)
	{
		c.ny = c.nz = 1 - c.q;
		c.hy = c.hz = 0;
		c.Sx = 1;
		c.Sy = c.Sz = 0;
		c.dV = c.hx;
		c.N = (c.nx + c.q);
	}

	if (c.dim == 2)
	{
		c.nz = 1 - c.q;
		c.hz = 0;
		c.Sx = c.hy;
		c.Sy = c.hx;
		c.Sz = 0;
		c.dV = c.hx * c.hy;
		c.N = (c.nx + c.q) * (c.ny + c.q);
		c.offset = c.nx + c.q;
		c.offset2 = 0;
	}

	if (c.dim == 3)
	{
		c.Sx = c.hy * c.hz;
		c.Sy = c.hx * c.hz;
		c.Sz = c.hx * c.hy;
		c.dV = c.hx * c.hy * c.hz;
		c.N = (c.nx + c.q) * (c.ny + c.q) * (c.nz + c.q);
		c.offset = c.nx + c.q;
		c.offset2 = (c.nx + c.q) * (c.ny + c.q);
	}



	c.ox = 0.5 * c.hx * (1 - c.q);
	c.oy = 0.5 * c.hy * (1 - c.q);
	c.oz = 0.5 * c.hz * (1 - c.q);

	c.Nbytes = c.N * sizeof(double);

	double pi = 3.1415926535897932384626433832795;
	auto make_vector = [&pi](double angle, double (*func)(double))
	{	return std::floor(func(angle * pi / 180.0) * 1e+10) / 1e+10;	};
	auto set_angles = [&make_vector, &c](double a, double b)
	{
		c.grav_y = make_vector(a, cos);
		c.grav_x = make_vector(a, sin);

		c.vibr_y = make_vector(b, cos);
		c.vibr_x = make_vector(b, sin);

		c.density_y = make_vector(a, cos);
		c.density_x = make_vector(a, sin);
	};
	set_angles(c.alpha, c.beta);
}



struct BackUp
{
private:
	struct OneField
	{
		int n = 0;
		double* d = nullptr;
		std::string name;
		OneField(int n_ = 0, double* d_ = nullptr, std::string name_ = "") : n(n_), d(d_), name(name_) {}
	};
	std::vector<OneField> fields;

	std::ofstream write[2];
	std::ifstream read;
	int copy = 1;

public:
	//BackUp(){}
	BackUp(bool withCopy = false)
	{
		if (withCopy) copy = 2;
	}

	void add(int n, double* d, std::string str)
	{
		fields.emplace_back(n, d, str);
	}
	void save_fields()
	{
		open_write();
		for (auto& it : fields)
		{
			for (int q = 0; q < copy; q++)
			{
				write[q] << it.name << " " << it.n << " ";
				for (int i = 0; i < it.n; i++)
					write[q] << it.d[i] << " ";
				write[q] << std::endl;
			}
		}
		close_write();
	}
	void recover(std::string str, double* f)
	{
		open_read();
		std::stringstream ss;
		std::string line;


		while (getline(read, line))
		{
			std::string symbol;
			ss.str("");
			ss.clear();
			ss << line;
			ss >> symbol;
			if (str == symbol)
			{
				int n = 0;
				ss >> n;
				for (int i = 0; i < n; i++)
					ss >> f[i];
			}
		}
	}



private:
	void open_write()
	{
		write[0].open("recovery.dat");
		if (copy == 2)	write[1].open("recovery2.dat");
	}
	void close_write()
	{
		write[0].close();
		if (copy == 2)	write[1].close();
	}
	void open_read()
	{
		read.open("recovery.dat");
	}
	void close_read()
	{
		read.close();
	}
};
