#pragma once
#include <random>
#include "toolkit.h"
#include <iostream>

struct Configuration
{
	double Lx, Ly, Lz, hx, hy, hz, Sx, Sy, Sz, dV, tau;
	double rho, nu, chi, betaT, gravity, typicalU, typicalL;
	int nx, ny, nz, off, off2, N, dim;

	double* vx, * vy, * vz, * p, * T;

	Configuration()
	{
		rho = 1000;
		nu = 8.9e-4;
		chi = 0.143e-6;
		betaT = 3.02e-4;
		gravity = 9.81;
		typicalU = 1;
		typicalL = 1;
		tau = 1.0;
		Lx = Ly = Lz = 1;
		hx = hy = hz = 1;
		dV = 1;
		nx = ny = nz = N = 1;
		off = off2 = 0;
		Sx = Sy = Sz = 0;
		dim = 3;
		if (nz > 1)
		{
			dim = 3;
			Sx = hy * hz;
			Sy = hx * hz;
			Sz = hx * hy;
			dV = hx * hy * hz;
		}

		if (ny > 1 && nz == 1)
		{
			dim = 2;
			Sx = hy;
			Sy = hx;
			Sz = 0.0;
			dV = hx * hy;
		}

		if (ny == 1 && nz == 1)
		{
			dim = 1;
			Sx = 1.0;
			Sy = 0.0;
			Sz = 0.0;
			dV = hx;
		}

		alloc(&vx, N);
		alloc(&vy, N);
		alloc(&vz, N);
		alloc(&p, N);
		alloc(&T, N);
	}



	void set_cell_number(int nx_, int ny_ = 1, int nz_ = 1)
	{
		nx = nx_;
		ny = ny_;
		nz = nz_;
		off = nx;
		off2 = nx * ny;
		N = nx * ny * nz;

		hx = Lx / nx;
		hy = Ly / ny;
		hz = Lz / nz;		
		set_surfaces();
	}
	void set_domain_LxLyLz(double Lx_, double Ly_ = 1, double Lz_ = 1)
	{
		Lx = Lx_;
		Ly = Ly_;
		Lz = Lz_;

		hx = Lx / nx;
		hy = Ly / ny;
		hz = Lz / nz;
		set_surfaces();
	}
	void set_surfaces()
	{
		if (nz > 1)
		{
			dim = 3;
			Sx = hy * hz;
			Sy = hx * hz;
			Sz = hx * hy;
			dV = hx * hy * hz;
		}

		if (ny > 1 && nz == 1)
		{
			dim = 2;
			Sx = hy;
			Sy = hx;
			Sz = 0.0;
			dV = hx * hy;
		}

		if (ny == 1 && nz == 1)
		{
			dim = 1;
			Sx = 1.0;
			Sy = 0.0;
			Sz = 0.0;
			dV = hx;
		}
	}
	void initFields()
	{
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::uniform_int_distribution<> distrib(0, 1);
		auto p_lin = [this](int i)
		{
			double x = hx * 0.5 + i * hx;
			return 5.0 * (1 - x / Lx);
		};
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + off * j;
					double x = i * hx + 0.5 * hx;
					double y = j * hy + 0.5 * hy;
					double z = k * hz + 0.5 * hz;

					vx[l] = 0;
					vy[l] = 0;
					vz[l] = 0;
					p[l] = 0;
					T[l] = 1 - y;

				}
			}
		}
	}

	void show_parameters()
	{
		using std::cout;
		using std::endl;

		cout << "dim = " << dim << endl;
		cout << "Lx,Ly,Lz = " << Lx << ", " << Ly << ", " << Lz << endl;
		cout << "hx,hy,hz = " << hx << ", " << hy << ", " << hz << endl;
		cout << "Sx,Sy,Sz = " << Sx << ", " << Sy << ", " << Sz << endl;
		cout << "Qx,Qy,Qz = " << Sx / hx << ", " << Sy / hy << ", " << Sz / hz << endl;
		cout << "nx,ny,nz = " << nx << ", " << ny << ", " << nz << endl;
		cout << "N = " << N << ", off = " << off << ", off2 = " << off2 << endl;
	}
};
