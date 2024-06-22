#pragma once
#include <iostream>
#include <fstream>
#include "toolkit.h"
#include <map>

enum class Dimension
{
	D1, D2, D3
};

enum class Component
{
	x, y, z
};

enum class Side
{
	center,
	west, east,
	south, north,
	front, back
};

enum class MathBoundary
{
	Dirichlet,
	Neumann,
	Periodic,
	Inlet,
	Outlet
};


//think of it in the framework of the finite volume approximation
//index i,j,k is a cell index
//get access to a velocity value on a side

#define GETVALUE(i, j, k) v[(i) + off * (j) + off2 * (k)]

#define GETFX(i2,i1) (GETVALUE((i2), j, k) + GETVALUE((i1), j, k)) * 0.5
#define GETFY(j2,j1_) (GETVALUE(i, (j2), k) + GETVALUE(i, (j1_), k)) * 0.5
#define GETFZ(k2,k1) (GETVALUE(i, j, (k2)) + GETVALUE(i, j, (k1))) * 0.5;

#define GETDX(i2,i1) (GETVALUE((i2), j, k) - GETVALUE((i1), j, k)) / hx
#define GETDY(j2,j1_) (GETVALUE(i, (j2), k) - GETVALUE(i, (j1_), k)) / hy
#define GETDZ(k2,k1) (GETVALUE(i, j, (k2)) - GETVALUE(i, j, (k1))) / hz

#define DIFFX(i2,i1) (GETVALUE((i2), j, k) - GETVALUE((i1), j, k))
#define DIFFY(j2,j1_) (GETVALUE(i, (j2), k) - GETVALUE(i, (j1_), k))
#define DIFFZ(k2,k1) (GETVALUE(i, j, (k2)) - GETVALUE(i, j, (k1)))

struct Boundary
{
	struct One
	{
		double value;
		double sign;
		double h;
		bool shifted = false;
		MathBoundary type;	
		One(double v = 0, double s = 0, double H = 0, MathBoundary t = MathBoundary::Neumann)
			: value(v), sign(s), h(H), type(t) {}

	};
	std::map<Side, One> b;
	int dim = 1;
	double* v = nullptr;
	int nx, ny, nz, off, off2;

	Boundary(MathBoundary def_type = MathBoundary::Neumann, double def_value = 0.0,
		double hx = 0, double hy = 0, double hz = 0)
	{
		b[Side::west]  = One(def_value, -1, hx, def_type);
		b[Side::east]  = One(def_value, +1, hx, def_type);
		b[Side::south] = One(def_value, -1, hy, def_type);
		b[Side::north] = One(def_value, +1, hy, def_type);
		b[Side::front] = One(def_value, -1, hz, def_type);
		b[Side::back]  = One(def_value, +1, hz, def_type);
		nx = ny = nz = off = off2 = 0;
	}
	void shifted_on(Side s)
	{
		b[s].shifted = true;
	}
	
	void bind_with_field(double* p, int nx_, int ny_, int nz_, int off_, int off2_)
	{
		v = p;
		nx = nx_;
		ny = ny_;
		nz = nz_;
		off = off_;
		off2 = off2_;
	}
	auto iterator(Side side)
	{
		return &b.find(side)->second;
	}
	MathBoundary type(Side side)
	{
		return b[side].type;
	}

	void set_boundary(Side side, MathBoundary bc, double v = 0.0, double h = 0.0)
	{
		auto it = iterator(side);
		it->type = bc;
		if (v != 0.0) it->value = v;
		if (h != 0.0) it->h = h;
 	}
	void set_period_pair(Side s1, Side s2)
	{
		bool found = false;

		auto check = [this, &found, &s1, &s2](bool ok)
		{
			if (ok)
			{
				double h1 = b[s1].h, h2 = b[s2].h;
				double h = h1 + h2;
				set_boundary(s1, MathBoundary::Periodic, 0.0, h);
				set_boundary(s2, MathBoundary::Periodic, 0.0, h);
				found = true;
			}
		};

		check(s1 == Side::west && s2 == Side::east);
		check(s2 == Side::west && s1 == Side::east);

		check(s1 == Side::south && s2 == Side::north);
		check(s2 == Side::south && s1 == Side::north);

		check(s1 == Side::front && s2 == Side::back);
		check(s2 == Side::front && s1 == Side::back);

		if (!found) std::cout << "not adequate pair" << std::endl;
	}

	void set_steps(double hx_, double hy_ = 0, double hz_ = 0)
	{
		b[Side::west].h = b[Side::east].h = hx_;
		b[Side::south].h = b[Side::north].h = hy_;
		b[Side::front].h = b[Side::back].h = hz_;
	}
	double retrieve_from_deriv(Side side, double f)
	{
		auto it = iterator(side);
		return f + it->sign * it->h * it->value;
	}
	double normal_deriv(Side side, int i = 0, int j = 0, int k = 0)
	{
		auto it = iterator(side);
		if (it->type == MathBoundary::Neumann)
			return it->value;
		if (it->type == MathBoundary::Periodic)
		{
			if (it->shifted == false)
			{
				switch (side)
				{
				case Side::west:		
					return (GETVALUE(0, j, k) - GETVALUE(nx - 1, j, k)) / (it->h);
					break;
				case Side::east:
					return (GETVALUE(0, j, k) - GETVALUE(nx - 1, j, k)) / (it->h);
					break;
				case Side::south:
					return (GETVALUE(i, 0, k) - GETVALUE(i, ny - 1, k)) / (it->h);
					break;
				case Side::north:
					return (GETVALUE(i, 0, k) - GETVALUE(i, ny - 1, k)) / (it->h);
					break;
				case Side::front:
					return (GETVALUE(i, j, 0) - GETVALUE(i, j, nz - 1)) / (it->h);
					break;
				case Side::back:
					return (GETVALUE(i, j, 0) - GETVALUE(i, j, nz - 1)) / (it->h);
					break;
				default:
					break;
				}
			} 
			else
			{
				switch (side)
				{
				case Side::west:
					return (GETVALUE(0, j, k) - GETVALUE(nx - 1, j, k)) / (it->h);
					break;
				case Side::east:
					return (GETVALUE(1, j, k) - GETVALUE(nx, j, k)) / (it->h);
					break;
				case Side::south:
					return (GETVALUE(i, 0, k) - GETVALUE(i, ny - 1, k)) / (it->h);
					break;
				case Side::north:
					return (GETVALUE(i, 1, k) - GETVALUE(i, ny, k)) / (it->h);
					break;
				case Side::front:
					return (GETVALUE(i, j, 0) - GETVALUE(i, j, nz - 1)) / (it->h);
					break;
				case Side::back:
					return (GETVALUE(i, j, 1) - GETVALUE(i, j, nz)) / (it->h);
					break;
				default:
					break;
				}
			}
		}

		return it->sign * (it->value - GETVALUE(i, j, k)) / (it->h);
	}

	double normal_difference(Side side, int i = 0, int j = 0, int k = 0)
	{
		return normal_deriv(side, i, j, k) * b[side].h;
	}

	double normal_deriv_oriented(Side side, int i = 0, int j = 0, int k = 0)
	{
		auto it = iterator(side);
		return normal_deriv(side, i, j, k) * it->sign;
	}
	void show_boundary(Side side)
	{
		std::cout << std::endl << "Boundary: " << std::endl;
		auto it = iterator(side);
		std::cout << "side: " << static_cast<int>(side) << ", ";
		std::cout << "value: " << it->value << ", ";
		std::cout << "type: " << static_cast<int>(it->type) << ", ";
		std::cout << "sign: " << static_cast<int>(it->sign) << std::endl;
		std::cout << "step,h: " << it->h << ", dim: " << dim << std::endl;
	}
	double operator()(Side side, int i = 0, int j = 0, int k = 0)
	{
		auto it = iterator(side);
		if (it->type == MathBoundary::Neumann)
			return retrieve_from_deriv(side, GETVALUE(i, j, k));
		else if (it->type == MathBoundary::Dirichlet)
			return b[side].value;
		else if (it->type == MathBoundary::Periodic)
		{
			if (it->shifted == false)
			{
				switch (side)
				{
				case Side::west:					
					return 0.5 * (GETVALUE(0, j, k) + GETVALUE(nx - 1, j, k));
					break;
				case Side::east:
					return 0.5 * (GETVALUE(0, j, k) + GETVALUE(nx - 1, j, k));
					break;
				case Side::south:
					return 0.5 * (GETVALUE(i, 0, k) + GETVALUE(i, ny - 1, k));
					break;
				case Side::north:
					return 0.5 * (GETVALUE(i, 0, k) + GETVALUE(i, ny - 1, k));
					break;
				case Side::front:
					return 0.5 * (GETVALUE(i, j, 0) + GETVALUE(i, j, nz - 1));
					break;
				case Side::back:
					return 0.5 * (GETVALUE(i, j, 0) + GETVALUE(i, j, nz - 1));
					break;
				default:
					break;
				}
			}
			else
			{
				switch (side)
				{
				case Side::west:					
					return 0.5 * (GETVALUE(0, j, k) + GETVALUE(nx - 1, j, k));
					break;
				case Side::east:
					return 0.5 * (GETVALUE(1, j, k) + GETVALUE(nx, j, k));
					break;
				case Side::south:
					return 0.5 * (GETVALUE(i, 0, k) + GETVALUE(i, ny - 1, k));
					break;
				case Side::north:
					return 0.5 * (GETVALUE(i, 1, k) + GETVALUE(i, ny, k));
					break;
				case Side::front:
					return 0.5 * (GETVALUE(i, j, 0) + GETVALUE(i, j, nz - 1));
					break;
				case Side::back:
					return 0.5 * (GETVALUE(i, j, 1) + GETVALUE(i, j, nz));
					break;
				default:
					break;
				}
			}
		}
		MYERROR("type operator()");
		return NAN;
	}
	double get_value_from_deriv(Side side, int i = 0, int j = 0, int k = 0)
	{
		if (b[side].type == MathBoundary::Neumann)
			return retrieve_from_deriv(side, GETVALUE(i, j, k));
		else
			MYERROR("incorrect bc type");
	}
	double get_fixed_value(Side side)
	{
		if (b[side].type == MathBoundary::Dirichlet)
			return b[side].value;
		//else if (b[side].type == MathBoundary::Periodic)
		else
			MYERROR("incorrect bc type");
	}

	One& operator[](Side side)
	{
		return b[side];
	}
};




struct ScalarVariable
{
	double* v = nullptr; 
	double hx = 1, hy = 1, hz = 1;
	double x0, y0, z0;
	double Sx = 0, Sy = 0, Sz = 0;
	int nx = 1, ny = 1, nz = 1;
	int off = 0, off2 = 0, N = 0, Nb = 0;
	int xExtra = 0, yExtra = 0, zExtra = 0;
	char dim = 1;
	Boundary boundary;
	ScalarVariable()
	{
		default_(false, nullptr);
	}

	ScalarVariable(int nx_, int ny_, int nz_, double hx_, double hy_, double hz_, bool newPointer = true, double* ptr = nullptr)
		: nx(nx_), ny(ny_), nz(nz_), hx(hx_), hy(hy_), hz(hz_)
	{
		dim = 3;
		N = nx * ny * nz;
		off = nx;
		off2 = nx * ny;

		default_(newPointer, ptr);
	}

	void default_(bool newPointer, double *ptr)
	{
		if (nz == 1)
		{
			off2 = 0;
			dim = 2;
		}
		if (ny == 1)
		{
			off = 0;
			dim = 1;
		}

		Sx = hy * hz;
		Sy = hx * hz;
		if (dim < 2) Sy = 0.0;
		Sz = hx * hy;
		if (dim < 3) Sz = 0.0;

		x0 = 0.5 * hx;
		y0 = 0.5 * hy;
		z0 = 0.5 * hz;
		boundary.set_boundary(Side::west,  MathBoundary::Neumann, 0.0);
		boundary.set_boundary(Side::east,  MathBoundary::Neumann, 0.0);
		boundary.set_boundary(Side::south, MathBoundary::Neumann, 0.0);
		boundary.set_boundary(Side::north, MathBoundary::Neumann, 0.0);
		boundary.set_boundary(Side::front, MathBoundary::Neumann, 0.0);
		boundary.set_boundary(Side::back,  MathBoundary::Neumann, 0.0);
		boundary.set_steps(0.5 * hx, 0.5 * hy, 0.5 * hz);


		if (newPointer)
		{
			v = new double[N];
			for (int l = 0; l < N; l++)
				v[l] = 0.0;
		}
		else if (!newPointer)
		{
			if (ptr != nullptr)
				v = ptr;
		}

		boundary.bind_with_field(v, nx, ny, nz, off, off2);
	}

	double get_dx(Side side, int i, int j = 0, int k = 0)
	{
		if (side == Side::west)
		{
			if (i == 0)
				return boundary.normal_deriv(side, i, j, k);
			else
				return GETDX(i, i - 1);
		}
		
		if (side == Side::east)
		{
			if (i == nx - 1)
				return boundary.normal_deriv(side, i, j, k);
			else
				return GETDX(i + 1, i);
		}
		if (side == Side::south)
		{
			if (j == 0)
			{
				if (i == 0 && boundary.type(Side::west) == MathBoundary::Periodic)
					return (boundary(side, 1, j, k) - (boundary(side, nx - 1, j, k))) / (2 * hx);
				if (i == nx - 1 && boundary.type(Side::east) == MathBoundary::Periodic)
					return (boundary(side, 0, j, k) - (boundary(side, nx - 2, j, k))) / (2 * hx);
				if (i == 0)				
					return (boundary(side, i + 1, j, k) - (boundary(side, i, j, k))) / (hx);
				else if (i == nx - 1)	
					return (boundary(side, i, j, k) - (boundary(side, i - 1, j, k))) / (hx);
				else                    
					return (boundary(side, i + 1, j, k) - (boundary(side, i - 1, j, k))) / (2 * hx);
			}
			else
			{
				double f1 = 0.5 * (get_at_side(Side::west, i, j, k) + get_at_side(Side::west, i, j - 1, k));
				double f2 = 0.5 * (get_at_side(Side::east, i, j, k) + get_at_side(Side::east, i, j - 1, k));
				return (f2 - f1) / hx;
			}
		}

		if (side == Side::north)
		{
			if (j == ny - 1)
			{
				if (i == 0 && boundary.type(Side::west) == MathBoundary::Periodic)
					return (boundary(side, 1, j, k) - (boundary(side, nx - 1, j, k))) / (2 * hx);
				if (i == nx - 1 && boundary.type(Side::east) == MathBoundary::Periodic)
					return (boundary(side, 0, j, k) - (boundary(side, nx - 2, j, k))) / (2 * hx);				
				if (i == 0)				return (boundary(side, i + 1, j, k) - (boundary(side, i, j, k))) / (hx);
				else if (i == nx - 1)	return (boundary(side, i, j, k) -     (boundary(side, i - 1, j, k))) / (hx);
				else                    return (boundary(side, i + 1, j, k) - (boundary(side, i - 1, j, k))) / (2 * hx);
			}
			else
			{
				double f1 = 0.5 * (get_at_side(Side::west, i, j, k) + get_at_side(Side::west, i, j + 1, k));
				double f2 = 0.5 * (get_at_side(Side::east, i, j, k) + get_at_side(Side::east, i, j + 1, k));
				return (f2 - f1) / hx;
			}
		}


		if (side == Side::center)
		{
			return (get_at_side(Side::east, i, j, k) - get_at_side(Side::west, i, j, k)) / (hx);
		}

		throw std::invalid_argument("wrong side");
		return NAN;
	}
	double get_dy(Side side, int i, int j, int k = 0)
	{
		if (side == Side::south)
		{
			if (j == 0)
				return boundary.normal_deriv(side, i, j, k);
			else
				return GETDY(j, j - 1);
		}

		if (side == Side::north)
		{
			if (j == ny - 1)
				return boundary.normal_deriv(side, i, j, k);
			else
				return GETDY(j + 1, j);
		}

		if (side == Side::west)
		{
			if (i == 0)
			{
				if (j == 0 && boundary.type(Side::south) == MathBoundary::Periodic)
					return (boundary(side, i, 1, k) - (boundary(side, i, ny - 1, k))) / (2 * hy);
				if (j == ny - 1 && boundary.type(Side::north) == MathBoundary::Periodic)
					return (boundary(side, i, 0, k) - (boundary(side, i, ny - 2, k))) / (2 * hy);

				if (j == 0)				return (boundary(side, i, j + 1, k) - (boundary(side, i, j, k))) / (hy);
				else if (j == ny - 1)	return (boundary(side, i, j, k)		- (boundary(side, i , j - 1, k))) / (hy);
				else                    return (boundary(side, i, j + 1, k) - (boundary(side, i, j - 1, k))) / (2 * hy);
			}
			else
			{
				double f1 = 0.5 * (get_at_side(Side::south, i, j, k) + get_at_side(Side::south, i - 1, j , k));
				double f2 = 0.5 * (get_at_side(Side::north, i, j, k) + get_at_side(Side::north, i - 1, j , k));
				return (f2 - f1) / hy;
			}
		}

		if (side == Side::east)
		{
			if (i == nx - 1)
			{
				if (j == 0 && boundary.type(Side::south) == MathBoundary::Periodic)
					return (boundary(side, i, 1, k) - (boundary(side, i, ny - 1, k))) / (2 * hy);
				if (j == ny - 1 && boundary.type(Side::north) == MathBoundary::Periodic)
					return (boundary(side, i, 0, k) - (boundary(side, i, ny - 2, k))) / (2 * hy);

				if (j == 0)				return (boundary(side, i, j + 1, k) - (boundary(side, i, j, k))) / (hy);
				else if (j == ny - 1)	return (boundary(side, i, j, k) - (boundary(side, i, j - 1, k))) / (hy);
				else                    return (boundary(side, i, j + 1, k) - (boundary(side, i, j - 1, k))) / (2 * hy);
			}
			else
			{
				double f1 = 0.5 * (get_at_side(Side::south, i, j, k) + get_at_side(Side::south, i + 1, j, k));
				double f2 = 0.5 * (get_at_side(Side::north, i, j, k) + get_at_side(Side::north, i + 1, j, k));
				return (f2 - f1) / hy;
			}
		}


		if (side == Side::center)
		{
			return (get_at_side(Side::north, i, j, k) - get_at_side(Side::south, i, j, k)) / (hy);
		}

		throw std::invalid_argument("wrong side");
		return NAN;
	}
	double get_dz(Side side, int i, int j, int k)
	{
		if (side == Side::front)
		{
			if (k == 0)
				return boundary.normal_deriv(side, i, j, k);
			else
				return GETDZ(k, k - 1);
		}
		if (side == Side::back)
		{
			if (k == nz - 1)
				return boundary.normal_deriv(side, i, j, k);
			else
				return GETDZ(k + 1, k);
		}
		if (side == Side::center)
		{
			return (get_at_side(Side::back, i, j, k) - get_at_side(Side::front, i, j, k)) / hz;
		}

		throw std::invalid_argument("wrong side");
		return NAN;
	}

	double get_diff_x(Side side, int i, int j = 0, int k = 0)
	{
		if (side == Side::west)
		{
			if (i <= 0)
				return boundary.normal_difference(side, 0, j, k);
			else
				return DIFFX(i, i - 1);
		}

		if (side == Side::east)
		{
			if (i >= nx - 1)
				return boundary.normal_difference(side, nx - 1, j, k);
			else
				return DIFFX(i + 1, i);
		}
		if (side == Side::center)
		{
			return (get_at_side(Side::east, i, j, k) - get_at_side(Side::west, i, j, k));
		}

		throw std::invalid_argument("wrong side");
		return NAN;
	}
	double get_diff_y(Side side, int i, int j, int k = 0)
	{
		if (side == Side::south)
		{
			if (j <= 0)
				return boundary.normal_difference(side, i, 0, k);
			else
				return DIFFY(j, j - 1);
		}

		if (side == Side::north)
		{
			if (j >= ny - 1)
				return boundary.normal_difference(side, i, ny - 1, k);
			else
				return DIFFY(j + 1, j);
		}
		if (side == Side::center)
		{
			return (get_at_side(Side::north, i, j, k) - get_at_side(Side::south, i, j, k));
		}

		throw std::invalid_argument("wrong side");
		return NAN;
	}
	double get_diff_z(Side side, int i, int j, int k)
	{
		if (side == Side::front)
		{
			if (k <= 0)
				return boundary.normal_difference(side, i, j, 0);
			else
				return DIFFZ(k, k - 1);
		}
		if (side == Side::back)
		{
			if (k >= nz - 1)
				return boundary.normal_difference(side, i, j, nz - 1);
			else
				return DIFFZ(k + 1, k);
		}
		if (side == Side::center)
		{
			return (get_at_side(Side::back, i, j, k) - get_at_side(Side::front, i, j, k));
		}

		throw std::invalid_argument("wrong side");
		return NAN;
	}




	double get_at_side(Side side, int i, int j = 0, int k = 0)
	{
		switch (side)
		{
		case Side::west:
			if (i == 0) 
				return boundary(side, i, j, k);
			return GETFX(i, i - 1);
			break;
		case Side::east:
			if (i == nx - 1)
				return boundary(side, i, j, k);
			return GETFX(i, i + 1);
			break;
		case Side::south:
			if (j == 0)
				return boundary(side, i, j, k);
			return GETFY(j, j - 1);
			break;
		case Side::north:
			if (j == ny - 1)
				return boundary(side, i, j, k);
			return GETFY(j, j + 1);
			break;
		case Side::front:
			if (k == 0)
				return boundary(side, i, j, k);
			return GETFZ(k, k - 1);
			break;
		case Side::back:
			if (k == nz - 1)
				return boundary(side, i, j, k);
			return GETFZ(k, k + 1);
			break;
		default:
			MYERROR("get at side");
			break;
		}
	}
	double get_shifted(Component comp, int i, int j = 0, int k = 0)
	{
		switch (comp)
		{
		case Component::x:
			if (i == 0)
				return boundary(Side::west, 0, j, k);
			else if (i == nx)
				return boundary(Side::east, nx - 1, j, k);
			else
				return GETFX(i, i - 1);
			break;
		case Component::y:
			if (j == 0)
				return boundary(Side::south, i, 0, k);
			else if (j == ny)
				return boundary(Side::north, i, ny - 1, k);
			else
				return GETFY(j, j - 1);
			break;
		case Component::z:
			if (k == 0)
				return boundary(Side::front, i, j, 0);
			else if (k == nz)
				return boundary(Side::back, i, j, nz - 1);
			else
				return GETFZ(k, k - 1);
			break;
		default:
			return 0;
			break;
		}
	}
	double get_shifted_diff(Component comp, int i, int j = 0, int k = 0)
	{
		switch (comp)
		{
		case Component::x:
			if (i == 0)
				return boundary.normal_difference(Side::west, 0, j, k);
			else if (i == nx)
				return boundary.normal_difference(Side::east, nx - 1, j, k);
			else
				return DIFFX(i, i - 1);
			break;
		case Component::y:
			if (j == 0)
				return boundary.normal_difference(Side::south, i, 0, k);
			else if (j == ny)
				return boundary.normal_difference(Side::north, i, ny - 1, k);
			else
				return DIFFY(j, j - 1);
			break;
		case Component::z:
			if (k == 0)
				return boundary.normal_difference(Side::front, i, j, 0);
			else if (k == nz)
				return boundary.normal_difference(Side::back, i, j, nz - 1);
			else
				return DIFFZ(k, k - 1);
			break;
		default:			
			return 0;
			break;
		}
	}
	double get_shifted_deriv(Component comp, int i, int j = 0, int k = 0)
	{
		switch (comp)
		{
		case Component::x:
			if (i == 0)
				return boundary.normal_deriv(Side::west, 0, j, k);
			else if (i == nx)
				return boundary.normal_deriv(Side::east, nx - 1, j, k);
			else
				return GETDX(i, i - 1);
			break;
		case Component::y:
			if (j == 0)
				return boundary.normal_deriv(Side::south, i, 0, k);
			else if (j == ny)
				return boundary.normal_deriv(Side::north, i, ny - 1, k);
			else
				return GETDY(j, j - 1);
			break;
		case Component::z:
			if (k == 0)
				return boundary.normal_deriv(Side::front, i, j, 0);
			else if (k == nz)
				return boundary.normal_deriv(Side::back, i, j, nz - 1);
			else
				return GETDZ(k, k - 1);
			break;
		default:
			return 0;
			break;
		}
	}

	double get_neighbor(Side side, int i, int j, int k)
	{
		switch (side)
		{
		case Side::west:			
			if (i == 0) return boundary(side, i, j, k);
			else return (*this)(i - 1, j, k);
			break;
		case Side::east:
			if (i == nx - 1) return boundary(side, i, j, k);
			else return (*this)(i + 1, j, k);
			break;
		case Side::south:			
			if (j == 0) return boundary(side, i, j, k);
			else return (*this)(i, j - 1, k);
			break;
		case Side::north:
			if (j == ny - 1) return boundary(side, i, j, k);
			else return (*this)(i, j + 1, k);
			break;
		case Side::front:			
			if (k == 0) return boundary(side, i, j, k);
			else return (*this)(i, j, k - 1);
			break;
		case Side::back:
			if (k == nz - 1) return boundary(side, i, j, k);
			else return (*this)(i, j, k + 1);
			break;
		default:
			break;
		}
	}

	double* get_ptr()
	{
		return v;
	}

	int get_l(int i, int j = 0, int k = 0)
	{
		return i + off * j + off2 * k;
	}

	double laplace_finite_volume(int i, int j, int k)
	{
		double lapl = 0.0;
		lapl += Sx * (get_dx(Side::east, i, j, k) -  get_dx(Side::west, i, j, k));
		lapl += Sy * (get_dy(Side::north, i, j, k) - get_dy(Side::south, i, j, k));
		lapl += Sz * (get_dz(Side::back, i, j, k) -  get_dz(Side::front, i, j, k));
		return lapl;
	}

	double& operator()(int i, int j = 0, int k = 0)
	{
		return v[i + off * j + off2 * k];
	}
	double& operator[](int l)
	{
		return v[l];
	}

	void reset_to_null()
	{
		for (int l = 0; l < N; l++)
			v[l] = 0.0;
	}

	static void swap(ScalarVariable& V1, ScalarVariable& V2)
	{
		double* temp;
		temp = V1.v;
		V1.v = V2.v;
		V2.v = temp;
	}

	void transfer_data_to(ScalarVariable& target)
	{
		int N_ = target.N;
		if (N_ != N)
		{
			std::cout << "not equal size" << std::endl;
			return;
		}
		
		for (int l = 0; l < N; l++)
		{
			target[l] = this->v[l];
		}
	}

	double x_(int i)
	{
		return i * hx + x0;
	}
	double y_(int j)
	{
		return j * hy + y0;
	}
	double z_(int k)
	{
		return k * hz + z0;
	}

	void write3D()
	{
		std::ofstream w("one_field.dat");
		w << "i, j, k, f" << std::endl;
		w << std::fixed;
		for (int k = 0; k < nz + zExtra; k++) {
			for (int j = 0; j < ny + yExtra; j++) {
				for (int i = 0; i < nx + xExtra; i++)
				{
					w << x_(i) << " " << y_(j) << " " << z_(k) << " " << (*this)(i, j, k) << std::endl;


				}
			}
		}
	}
	void write_line_x(int nx, int j, int k = 0)
	{
		std::ofstream w("line_x.dat");
		w << "i, f" << std::endl;

		for (int i = 0; i < nx; i++)
		{
			w << x_(i) << " " << (*this)(i, j, k) << std::endl;
		}

	}
	void write_line_y(int i, int k = 0)
	{
		std::ofstream w("line_y.dat");
		w << "j, f" << std::endl;

		for (int j = 0; j < ny; j++)
		{
			w << j << " " << (*this)(i, j, k) << std::endl;
		}

	}

	void set_all_boundaries(MathBoundary type, double val = 0.0)
	{
		boundary.set_boundary(Side::west,  type, val);
		boundary.set_boundary(Side::east,  type, val);
		boundary.set_boundary(Side::south, type, val);
		boundary.set_boundary(Side::north, type, val);
		boundary.set_boundary(Side::front, type, val);
		boundary.set_boundary(Side::back, type, val);
	}
	void set_linear()
	{
		double FW = 0, FE = 0;
		double FS = 0, FN = 0;
		double FF = 0, FB = 0;

		if (boundary.type(Side::west) == MathBoundary::Dirichlet)
			FW = boundary.get_fixed_value(Side::west);
		if (boundary.type(Side::east) == MathBoundary::Dirichlet)
			FE = boundary.get_fixed_value(Side::east);

		if (boundary.type(Side::south) == MathBoundary::Dirichlet)
			FS = boundary.get_fixed_value(Side::south);
		if (boundary.type(Side::north) == MathBoundary::Dirichlet)
			FN = boundary.get_fixed_value(Side::north);

		if (boundary.type(Side::front) == MathBoundary::Dirichlet)
			FF = boundary.get_fixed_value(Side::front);
		if (boundary.type(Side::back) == MathBoundary::Dirichlet)
			FB = boundary.get_fixed_value(Side::back);

		double Lx = hx * nx;
		double Ly = hy * ny;
		double Lz = hz * nz;

		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {

					(*this)(i, j, k) = 
						(FE - FW) / Lx * x_(i) + FW +
						(FN - FS) / Ly * y_(j) + FS +
						(FB - FF) / Lz * z_(k) + FF;

				}
			}
		}
	}

	void show_max_min()
	{
		double max = -1e+10;
		double min = +1e+10;

		for (int l = 0; l < N; l++)
		{
			if (v[l] > max) max = v[l];
			if (v[l] < min) min = v[l];
		}
		std::cout << "max = " << max << ", min = " << min << std::endl;
	}
	void show_info()
	{
		using std::cout;
		using std::endl;
		cout << "dim = " << (int)dim << ", (nx=" << nx << ",ny=" << ny << ",nz=" << nz << "), N = " << N << endl;
		cout << "offset = " << off << ", offset2 = " << off2 << endl;
		cout << "hx, hy, hz = " << hx << " " << hy << " " << hz << endl;
		cout << "Sx, Sy, Sz = " << Sx << " " << Sy << " " << Sz << endl;

		cout << endl << "Boundaries: " << endl;
		cout << "West: "; boundary.show_boundary(Side::west);
		cout << "East: "; boundary.show_boundary(Side::east);

		if (dim > 1)
		{
			cout << "South: "; boundary.show_boundary(Side::south);
			cout << "North: "; boundary.show_boundary(Side::north);
		}
		if (dim > 2)
		{
			cout << "Front: "; boundary.show_boundary(Side::front);
			cout << "Back: "; boundary.show_boundary(Side::back);
		}
	}

	
};


struct Velocity : public ScalarVariable
{
	Component type = Component::x;
	Velocity()
	{

	}

	Velocity(Component type_, int nx_, int ny_, int nz_,
		double hx_, double hy_, double hz_,
		bool newPointer = true, double* ptr = nullptr)
		: ScalarVariable(nx_, ny_, nz_, hx_, hy_, hz_, newPointer, ptr), type(type_)
	{
		dim = 3;
		if (nz == 1) dim = dim - 1;
		if (ny == 1) dim = dim - 1;

		if (type == Component::x)
		{
			if (dim >= 2)	off = nx + 1;
			if (dim >= 3)	off2 = (nx + 1) * ny;
			N = (nx + 1) * ny * nz;
			xExtra = 1;
		}

		if (type == Component::y)
		{
			if (dim >= 2)	off = nx;
			if (dim >= 3)	off2 = (ny + 1) * nz;
			N = nx * (ny + 1) * nz;
			yExtra = 1;
		}

		if (type == Component::z)
		{
			if (dim >= 2)	off = nx;
			if (dim >= 3)	off2 = nx * ny;
			N = nx * ny * (nz + 1);
			zExtra = 1;
		}

		default_(newPointer, ptr);

		boundary.set_boundary(Side::west,  MathBoundary::Dirichlet, 0.0);
		boundary.set_boundary(Side::east,  MathBoundary::Dirichlet, 0.0);
		boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 0.0);
		boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0.0);
		boundary.set_boundary(Side::front, MathBoundary::Dirichlet, 0.0);
		boundary.set_boundary(Side::back,  MathBoundary::Dirichlet, 0.0);
		boundary.set_steps(hx, hy, hz);

		if (type == Component::x) boundary.set_steps(0.5 * hx, 0.5 * hy, 0.5 * hz);
		if (type == Component::y) boundary.set_steps(0.5 * hx, 0.5 * hy, 0.5 * hz);
		if (type == Component::z) boundary.set_steps(0.5 * hx, 0.5 * hy, 0.5 * hz);

		if (type == Component::x) x0 = 0;
		if (type == Component::y) y0 = 0;
		if (type == Component::z) z0 = 0;

		if (type == Component::x)
		{
			boundary.shifted_on(Side::west); 
			boundary.shifted_on(Side::east);
		}
		if (type == Component::y)
		{
			boundary.shifted_on(Side::south);
			boundary.shifted_on(Side::north);
		}		
		if (type == Component::z)
		{
			boundary.shifted_on(Side::front);
			boundary.shifted_on(Side::back);
		}

		//if (type == Component::x) boundary.bind_with_field(v, nx + 1, ny, nz, off, off2);
		//if (type == Component::y) boundary.bind_with_field(v, nx, ny + 1, nz, off, off2);
		//if (type == Component::z) boundary.bind_with_field(v, nx, ny, nz + 1, off, off2);
		boundary.bind_with_field(v, nx, ny, nz, off, off2);
	}

	void get_at_side(Side side, int i, int j = 0, int k = 0) = delete;
	void laplace_finite_volume(Side side, int i, int j = 0, int k = 0) = delete;

	double get_for_centered_cell(Side side, int i, int j = 0, int k = 0)
	{
		if (type == Component::x)
		{
			if (side == Side::west)
			{
				return GETVALUE(i, j, k);
			}
			if (side == Side::east)
			{
				return GETVALUE(i + 1, j, k);
			}
			if (side == Side::south)
			{
				if (j == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i + 1, j, k) +
						GETVALUE(i, j - 1, k) +
						GETVALUE(i + 1, j - 1, k));
			}
			if (side == Side::north)
			{
				if (j == ny - 1)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i + 1, j, k) +
						GETVALUE(i, j + 1, k) +
						GETVALUE(i + 1, j + 1, k));
			}
			if (side == Side::front)
			{
				if (k == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i + 1, j, k) +
						GETVALUE(i, j, k - 1) +
						GETVALUE(i + 1, j, k - 1));
			}
			if (side == Side::back)
			{
				if (k == nz - 1)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i + 1, j, k) +
						GETVALUE(i, j, k + 1) +
						GETVALUE(i + 1, j, k + 1));
			}
			if (side == Side::center)
			{
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i + 1, j, k));
			}
		}
		if (type == Component::y)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j + 1, k) +
						GETVALUE(i - 1, j + 1, k) +
						GETVALUE(i - 1, j, k));
			}
			if (side == Side::east)
			{
				if (i == nx - 1)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j + 1, k) +
						GETVALUE(i + 1, j + 1, k) +
						GETVALUE(i + 1, j, k));
			}
			if (side == Side::south)
			{
				return GETVALUE(i, j, k);
			}
			if (side == Side::north)
			{
				return GETVALUE(i, j + 1, k);
			}
			if (side == Side::front)
			{
				if (k == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j + 1, k) +
						GETVALUE(i, j, k - 1) +
						GETVALUE(i, j + 1, k - 1));
			}
			if (side == Side::back)
			{
				if (k == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j + 1, k) +
						GETVALUE(i, j, k + 1) +
						GETVALUE(i, j + 1, k + 1));
			}
			if (side == Side::center)
			{
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j + 1, k));
			}
		}
		if (type == Component::z)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j, k + 1) +
						GETVALUE(i - 1, j, k) +
						GETVALUE(i - 1, j, k + 1));
			}
			if (side == Side::east)
			{
				if (i == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j, k + 1) +
						GETVALUE(i + 1, j, k) +
						GETVALUE(i + 1, j, k + 1));
			}
			if (side == Side::south)
			{
				if (j == 0)
					return boundary(side, i, j, k);
				else
					return 0.25 * (
						GETVALUE(i, j, k) +
						GETVALUE(i, j - 1, k) +
						GETVALUE(i, j, k + 1) +
						GETVALUE(i, j - 1, k + 1));
			}
			if (side == Side::front)
			{
				return GETVALUE(i, j, k);
			}
			if (side == Side::back)
			{
				return GETVALUE(i, j, k + 1);
			}
			if (side == Side::center)
			{
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k + 1));
			}
		}

		return 0.0;
	}
	double get_for_vx_cell(Side side, int i, int j = 0, int k = 0)
	{
		if (type == Component::x)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i - 1, j, k));
			}
			if (side == Side::east)
			{
				if (i == nx)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i + 1, j, k));
			}
			if (side == Side::south)
			{
				if (j == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j - 1, k));
			}
			if (side == Side::north)
			{
				if (j == ny - 1)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j + 1, k));
			}
			if (side == Side::front)
			{
				if (k == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k - 1));
			}
			if (side == Side::back)
			{
				if (k == nz - 1)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k + 1));
			}
			if (side == Side::center)
			{
				return GETVALUE(i, j, k);
			}
		}
		if (type == Component::y)
		{
			//if (side == Side::west)
			//{
			//	if (i == 0) 
			//		return 0.5 * (boundary(side, i, j, k) + boundary(side, i, j + 1, k));
			//	return 0.5 * (GETVALUE(i - 1, j, k) + GETVALUE(i - 1, j + 1, k));
			//}
			//if (side == Side::east)
			//{
			//	if (i == nx)
			//		return 0.5 * (boundary(side, i, j, k) + boundary(side, i, j + 1, k));
			//	return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j + 1, k));
			//}
			if (side == Side::south)
			{
				//wtf
				if (i == 0) 
					return boundary(Side::west, i, j, k);
				else if (i == nx)
					return boundary(Side::east, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i - 1, j, k));
			}
			if (side == Side::north)
			{
				if (i == 0)
					return boundary(Side::west, i, j + 1, k);
				else if (i == nx)
					return boundary(Side::east, i, j + 1, k);
				return 0.5 * (GETVALUE(i, j + 1, k) + GETVALUE(i - 1, j + 1, k));
			}

			//if (side == Side::front)
			//{
			//	if (k == 0)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//				GETVALUE(i, j, k) +
			//				GETVALUE(i, j + 1, k) +
			//				GETVALUE(i - 1, j + 1, k) +
			//				GETVALUE(i - 1, j, k) +
			//				GETVALUE(i, j, k - 1) +
			//				GETVALUE(i, j + 1, k - 1) +
			//				GETVALUE(i - 1, j + 1, k - 1) +
			//				GETVALUE(i - 1, j, k - 1));
			//}
			//if (side == Side::back)
			//{
			//	if (k == nz - 1)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j + 1, k) +
			//			GETVALUE(i - 1, j + 1, k) +
			//			GETVALUE(i - 1, j, k) +
			//				GETVALUE(i, j, k + 1) +
			//				GETVALUE(i, j + 1, k + 1) +
			//				GETVALUE(i - 1, j + 1, k + 1) +
			//				GETVALUE(i - 1, j, k + 1));
			//}
			if (side == Side::center)
			{
				return 0.25 * (
					GETVALUE(i, j, k) +
					GETVALUE(i, j + 1, k) +
					GETVALUE(i - 1, j + 1, k) +
					GETVALUE(i - 1, j, k));
			}
		}
		if (type == Component::z)
		{
			//if (side == Side::west)
			//{
			//	return 0.5 * (GETVALUE(i - 1, j, k) + GETVALUE(i - 1, j, k + 1));
			//}
			//if (side == Side::east)
			//{
			//	return 0.5 * (GETVALUE(i + 1, j, k) + GETVALUE(i + 1, j, k + 1));
			//}
			//if (side == Side::south)
			//{
			//	if (j == 0)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j, k + 1) +
			//			GETVALUE(i - 1, j, k + 1) +
			//			GETVALUE(i - 1, j, k) +
			//			GETVALUE(i, j - 1, k) +
			//			GETVALUE(i, j - 1, k + 1) +
			//			GETVALUE(i - 1, j - 1, k + 1) +
			//			GETVALUE(i - 1, j - 1, k)
			//			);
			//}
			//if (side == Side::north)
			//{
			//	if (j == ny - 1)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j, k + 1) +
			//			GETVALUE(i - 1, j, k + 1) +
			//			GETVALUE(i - 1, j, k) +
			//			GETVALUE(i, j + 1, k) +
			//			GETVALUE(i, j + 1, k + 1) +
			//			GETVALUE(i - 1, j + 1, k + 1) +
			//			GETVALUE(i - 1, j + 1, k)
			//			);
			//}
			if (side == Side::front)
			{
				if (i == 0)
					return boundary(Side::west, i, j, k);
				else if (i == nx)
					return boundary(Side::east, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i - 1, j, k));
			}
			if (side == Side::back)
			{
				if (i == 0)
					return boundary(Side::west, i, j, k + 1);
				else if (i == nx)
					return boundary(Side::east, i, j, k + 1);
				return 0.5 * (GETVALUE(i, j, k + 1) + GETVALUE(i - 1, j, k + 1));
			}
			if (side == Side::center)
			{
				return 0.25 * (
					GETVALUE(i, j, k) +
					GETVALUE(i, j, k + 1) +
					GETVALUE(i - 1, j, k + 1) +
					GETVALUE(i - 1, j, k));
			}
		}
		
		std::cout << "should be here" << std::endl;
		return 0.0;
	}
	double get_for_vy_cell(Side side, int i, int j, int k = 0)
	{
		if (type == Component::x)
		{
			if (side == Side::west)
			{
				if (j == 0)
					return boundary(Side::south, i, j, k);
				else if (j == ny)
					return boundary(Side::north, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j - 1, k));
			}
			if (side == Side::east)
			{
				if (j == 0)
					return boundary(Side::south, i + 1, j, k);
				else if (j == ny)
					return boundary(Side::north, i + 1, j, k);
				return 0.5 * (GETVALUE(i + 1, j, k) + GETVALUE(i + 1, j - 1, k));
			}
			//if (side == Side::south)
			//{
			//	//if (j == 0)
			//	//	return boundary(side, i, j, k);
			//	//return 0.5 * (GETVALUE(i, j - 1, k) + GETVALUE(i + 1, j - 1, k));
			//}
			//if (side == Side::north)
			//{
			//	//if (j == ny)
			//	//	return boundary(side, i, j, k);
			//	//return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i + 1, j, k));
			//}
			//if (side == Side::front)
			//{
			//	//if (k == 0)
			//	//	return boundary(side, i, j, k);
			//	//else
			//	//	return 0.125 * (
			//	//		GETVALUE(i, j, k) +
			//	//		GETVALUE(i, j - 1, k) +
			//	//		GETVALUE(i + 1, j, k) +
			//	//		GETVALUE(i + 1, j - 1, k) +
			//	//		GETVALUE(i, j, k - 1) +
			//	//		GETVALUE(i, j - 1, k - 1) +
			//	//		GETVALUE(i + 1, j, k - 1) +
			//	//		GETVALUE(i + 1, j - 1, k - 1));
			//}
			//if (side == Side::back)
			//{
			//	//if (k == nz - 1)
			//	//	return boundary(side, i, j, k);
			//	//else
			//	//	return 0.125 * (
			//	//		GETVALUE(i, j, k) +
			//	//		GETVALUE(i, j - 1, k) +
			//	//		GETVALUE(i + 1, j, k) +
			//	//		GETVALUE(i + 1, j - 1, k) +
			//	//		GETVALUE(i, j, k + 1) +
			//	//		GETVALUE(i, j - 1, k + 1) +
			//	//		GETVALUE(i + 1, j, k + 1) +
			//	//		GETVALUE(i + 1, j - 1, k + 1));
			//}
			if (side == Side::center)
			{
				return 0.25 * (
					GETVALUE(i, j, k) +
					GETVALUE(i, j - 1, k) +
					GETVALUE(i + 1, j, k) +
					GETVALUE(i + 1, j - 1, k));
			}
		}
		if (type == Component::y)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i - 1, j, k));
			}
			if (side == Side::east)
			{
				if (i == nx - 1)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i + 1, j, k));
			}
			if (side == Side::south)
			{
				if (j == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j - 1, k));
			}
			if (side == Side::north)
			{
				if (j == ny)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j + 1, k));
			}
			if (side == Side::front)
			{
				if (k == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k - 1));
			}
			if (side == Side::back)
			{
				if (k == nz - 1)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k + 1));
			}
			if (side == Side::center)
			{
				return GETVALUE(i, j, k);
			}
		}
		if (type == Component::z)
		{
			//if (side == Side::west)
			//{
			//	if (i == 0)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j - 1, k + 1) +
			//			GETVALUE(i, j, k + 1) +
			//			GETVALUE(i, j - 1, k) +
			//			GETVALUE(i - 1, j, k) +
			//			GETVALUE(i - 1, j - 1, k + 1) +
			//			GETVALUE(i - 1, j, k + 1) +
			//			GETVALUE(i - 1, j - 1, k));
			//}
			//if (side == Side::east)
			//{
			//	if (i == nx - 1)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j - 1, k + 1) +
			//			GETVALUE(i, j, k + 1) +
			//			GETVALUE(i, j - 1, k) +
			//			GETVALUE(i + 1, j, k) +
			//			GETVALUE(i + 1, j - 1, k + 1) +
			//			GETVALUE(i + 1, j, k + 1) +
			//			GETVALUE(i + 1, j - 1, k));
			//}
			//if (side == Side::south)
			//{
			//	return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k + 1));
			//}
			//if (side == Side::north)
			//{
			//	return 0.5 * (GETVALUE(i, j + 1, k) + GETVALUE(i, j + 1, k + 1));
			//}
			if (side == Side::front)
			{
				if (j == 0)
					return boundary(Side::south, i, j, k);
				else if (j == ny)
					return boundary(Side::north, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j - 1, k));
			}
			if (side == Side::back)
			{
				if (j == 0)
					return boundary(Side::south, i, j, k + 1);
				else if (j == ny)
					return boundary(Side::north, i, j, k + 1);
				return 0.5 * (GETVALUE(i, j, k + 1) + GETVALUE(i, j - 1, k + 1));
			}
			if (side == Side::center)
			{
				return 0.25 * (
					GETVALUE(i, j, k) +
					GETVALUE(i, j - 1, k + 1) +
					GETVALUE(i, j, k + 1) +
					GETVALUE(i, j - 1, k));
			}
		}
		

		std::cout << "should be here" << std::endl;
		return 0.0;
	}
	double get_for_vz_cell(Side side, int i, int j, int k)
	{
		if (type == Component::x)
		{
			if (side == Side::west)
			{
				if (k == 0)
					return boundary(Side::front, i, j, k);
				else if (k == nz)
					return boundary(Side::back, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k - 1));
			}
			if (side == Side::east)
			{
				if (k == 0)
					return boundary(Side::front, i + 1, j, k);
				else if (k == nz)
					return boundary(Side::back, i + 1, j, k);
				return 0.5 * (GETVALUE(i + 1, j, k) + GETVALUE(i + 1, j, k - 1));
			}
		//	if (side == Side::south)
		//	{
		//		if (j == 0)
		//			return boundary(side, i, j, k);
		//		else
		//			return 0.125 * (
		//				GETVALUE(i, j, k) +
		//				GETVALUE(i, j, k - 1) +
		//				GETVALUE(i + 1, j, k - 1) +
		//				GETVALUE(i + 1, j, k) +
		//				GETVALUE(i, j - 1, k) +
		//				GETVALUE(i, j - 1, k - 1) +
		//				GETVALUE(i + 1, j - 1, k - 1) +
		//				GETVALUE(i + 1, j - 1, k));
		//	}
		//	if (side == Side::north)
		//	{
		//		if (j == ny - 1)
		//			return boundary(side, i, j, k);
		//		else
		//			return 0.125 * (
		//				GETVALUE(i, j, k) +
		//				GETVALUE(i, j, k - 1) +
		//				GETVALUE(i + 1, j, k - 1) +
		//				GETVALUE(i + 1, j, k) +
		//				GETVALUE(i, j + 1, k) +
		//				GETVALUE(i, j + 1, k - 1) +
		//				GETVALUE(i + 1, j + 1, k - 1) +
		//				GETVALUE(i + 1, j + 1, k));
		//	}
		//	if (side == Side::front)
		//	{
		//		if (k == 0)
		//			return boundary(side, i, j, k);
		//		else
		//			return 0.5 * (GETVALUE(i, j, k - 1) + GETVALUE(i + 1, j, k - 1));
		//	}
		//	if (side == Side::back)
		//	{
		//		if (k == nz - 1)
		//			return boundary(side, i, j, k);
		//		else
		//			return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i + 1, j, k));
		//	}
		//	if (side == Side::center)
		//	{
		//		return 0.25 * (
		//			GETVALUE(i, j, k) +
		//			GETVALUE(i, j, k - 1) +
		//			GETVALUE(i + 1, j, k - 1) +
		//			GETVALUE(i + 1, j, k));
		//	}
		}
		if (type == Component::y)
		{
			//if (side == Side::west)
			//{
			//	if (i == 0)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k - 1) +
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j + 1, k - 1) +
			//			GETVALUE(i, j + 1, k) +
			//			GETVALUE(i - 1, j, k - 1) +
			//			GETVALUE(i - 1, j, k) +
			//			GETVALUE(i - 1, j + 1, k - 1) +
			//			GETVALUE(i - 1, j + 1, k));
			//}
			//if (side == Side::east)
			//{
			//	if (i == nx - 1)
			//		return boundary(side, i, j, k);
			//	else
			//		return 0.125 * (
			//			GETVALUE(i, j, k - 1) +
			//			GETVALUE(i, j, k) +
			//			GETVALUE(i, j + 1, k - 1) +
			//			GETVALUE(i, j + 1, k) +
			//			GETVALUE(i + 1, j, k - 1) +
			//			GETVALUE(i + 1, j, k) +
			//			GETVALUE(i + 1, j + 1, k - 1) +
			//			GETVALUE(i + 1, j + 1, k));
			//}
			if (side == Side::south)
			{
				if (k == 0)
					return boundary(Side::front, i, j, k);
				else if (k == nz)
					return boundary(Side::back, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k - 1));
			}
			if (side == Side::north)
			{
				if (k == 0)
					return boundary(Side::front, i, j + 1, k);
				else if (k == nz)
					return boundary(Side::back, i, j + 1, k);
				return 0.5 * (GETVALUE(i, j + 1, k) + GETVALUE(i, j + 1, k - 1));
			}
			//if (side == Side::front)
			//{
			//	return 0.5 * (GETVALUE(i, j, k - 1) + GETVALUE(i, j + 1, k - 1));
			//}
			//if (side == Side::back)
			//{
			//	return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j + 1, k));
			//}
			if (side == Side::center)
			{
				return 0.25 * (
					GETVALUE(i, j, k - 1) +
					GETVALUE(i, j, k) +
					GETVALUE(i, j + 1, k - 1) +
					GETVALUE(i, j + 1, k));
			}
		}
		if (type == Component::z)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i - 1, j, k));
			}
			if (side == Side::east)
			{
				if (i == nx - 1)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i + 1, j, k));
			}
			if (side == Side::south)
			{
				if (j == 0)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j - 1, k));
			}
			if (side == Side::north)
			{
				if (j == ny - 1)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j + 1, k));
			}
			if (side == Side::front)
			{
				if (k == 0) 
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k - 1));
			}
			if (side == Side::back)
			{
				if (k == nz)
					return boundary(side, i, j, k);
				return 0.5 * (GETVALUE(i, j, k) + GETVALUE(i, j, k + 1));
			}
			if (side == Side::center)
			{
				return GETVALUE(i, j, k);
			}
		}
		
		std::cout << "should be here" << std::endl;
		return 0.0;
	}


	double get_dx_for_vx_cell(Side side, int i, int j = 0, int k = 0)
	{
		if (type == Component::x)
		{
			if (side == Side::west)
			{
				if (i == 0)
				{
					if (boundary[side].type == MathBoundary::Periodic)
						return GETDX(0, nx - 1);
				}
				return GETDX(i, i - 1);
			}
			if (side == Side::east)
			{
				if (i == nx)
				{
					if (boundary[side].type == MathBoundary::Periodic)
						return GETDX(1, nx);
				}
				return GETDX(i + 1, i);
			}
		}
		MYERROR("Velocity");
		return 0.0;
	}
	double get_dy_for_vx_cell(Side side, int i, int j, int k = 0)
	{
		if (type == Component::x)
		{
			if (side == Side::south)
			{
				if (j == 0)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDY(j, j - 1);
			}

			if (side == Side::north)
			{
				if (j == ny - 1)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDY(j + 1, j);
			}
		}
		MYERROR("Velocity");
		return 0.0;
	}
	double get_dz_for_vx_cell(Side side, int i, int j, int k)
	{
		if (type == Component::x)
		{
			if (side == Side::front)
			{
				if (k == 0)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDZ(k, k - 1);
			}
			if (side == Side::back)
			{
				if (k == nz - 1)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDZ(k + 1, k);
			}
		}
		MYERROR("Velocity");
		return 0.0;
	}

	double get_dx_for_vy_cell(Side side, int i, int j, int k = 0)
	{
		if (type == Component::y)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDX(i, i - 1);
			}
			if (side == Side::east)
			{
				if (i == nx - 1)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDX(i + 1, i);
			}
			MYERROR("Velocity");
			return 0.0;
		}
	}
	double get_dy_for_vy_cell(Side side, int i, int j, int k = 0)
	{
		if (type == Component::y)
		{
			if (side == Side::south)
				return GETDY(j, j - 1);

			if (side == Side::north)
				return GETDY(j + 1, j);
		}
		MYERROR("Velocity");
		return 0.0;
	}
	double get_dz_for_vy_cell(Side side, int i, int j, int k)
	{
		if (type == Component::y)
		{
			if (side == Side::front)
			{
				if (k == 0)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDZ(k, k - 1);
			}
			if (side == Side::back)
			{
				if (k == nz - 1)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDZ(k + 1, k);
			}
		}
		MYERROR("Velocity");
		return 0.0;
	}

	double get_dx_for_vz_cell(Side side, int i, int j, int k)
	{
		if (type == Component::z)
		{
			if (side == Side::west)
			{
				if (i == 0)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDX(i, i - 1);
			}
			if (side == Side::east)
			{
				if (i == nx - 1)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDX(i + 1, i);
			}
			MYERROR("Velocity");
			return 0.0;
		}
	}
	double get_dy_for_vz_cell(Side side, int i, int j, int k)
	{
		if (type == Component::z)
		{
			if (side == Side::south)
			{
				if (j == 0)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDY(j, j - 1);
			}

			if (side == Side::north)
			{
				if (j == ny - 1)
					return boundary.normal_deriv(side, i, j, k);
				else
					return GETDY(j + 1, j);
			}
		}
		MYERROR("Velocity");
		return 0.0;
	}
	double get_dz_for_vz_cell(Side side, int i, int j, int k)
	{
		if (type == Component::z)
		{
			if (side == Side::front)
				return GETDZ(k, k - 1);
			if (side == Side::back)
				return GETDZ(k + 1, k);
		}
	}

	//double laplace_for_vx_cell(Side side, int i, int j, int k)
	//{
	//	double lapl = 0.0;
	//	lapl += (get_dx_for_vx_cell(Side::east, i, j, k)
	//		- get_dx_for_vx_cell(Side::west, i, j, k)) * Sx;
	//	lapl += (get_dy_for_vx_cell(Side::north, i, j, k)
	//		- get_dy_for_vx_cell(Side::south, i, j, k)) * Sy;
	//	lapl += (get_dz_for_vx_cell(Side::back, i, j, k)
	//		- get_dz_for_vx_cell(Side::front, i, j, k)) * Sz;
	//	return lapl;
	//}
	//double laplace_for_vy_cell(Side side, int i, int j, int k)
	//{
	//	double lapl = 0.0;
	//	lapl += (get_dx_for_vy_cell(Side::east, i, j, k)
	//			- get_dx_for_vy_cell(Side::west, i, j, k)) * Sx;
	//	lapl += (get_dy_for_vy_cell(Side::north, i, j, k)
	//			- get_dy_for_vy_cell(Side::south, i, j, k)) * Sy;
	//	lapl += (get_dz_for_vy_cell(Side::back, i, j, k)
	//			- get_dz_for_vy_cell(Side::front, i, j, k)) * Sz;
	//	return lapl;
	//}
	//double laplace_for_vz_cell(Side side, int i, int j, int k)
	//{
	//	double lapl = 0.0;
	//	lapl += (get_dx_for_vz_cell(Side::east, i, j, k)
	//		   - get_dx_for_vz_cell(Side::west, i, j, k)) * Sx;
	//	lapl += (get_dy_for_vz_cell(Side::north, i, j, k)
	//		   - get_dy_for_vz_cell(Side::south, i, j, k)) * Sy;
	//	lapl += (get_dz_for_vz_cell(Side::back, i, j, k)
	//		   - get_dz_for_vz_cell(Side::front, i, j, k)) * Sz;
	//	return lapl;
	//}

	//static double divergence_finite_volume(Velocity &FX, Velocity&FY, Velocity&FZ,	int i, int j, int k)
	//{
	//	double div = 0.0;
	//	div += FX.Sx * (FX.get_for_centered_cell(Side::east, i, j, k) - FX.get_for_centered_cell(Side::west, i, j, k));
	//	div += FY.Sy * (FY.get_for_centered_cell(Side::north, i, j, k) - FY.get_for_centered_cell(Side::south, i, j, k));
	//	div += FZ.Sz * (FZ.get_for_centered_cell(Side::back, i, j, k) - FZ.get_for_centered_cell(Side::front, i, j, k));
	//	return div;
	//}
	//static double divergence_finite_volume(Velocity& FX, Velocity& FY, 	int i, int j)
	//{
	//	double div = 0.0;
	//	div += FX.Sx * (FX.get_for_centered_cell(Side::east, i, j) - FX.get_for_centered_cell(Side::west, i, j));
	//	div += FY.Sy * (FY.get_for_centered_cell(Side::north, i, j) - FY.get_for_centered_cell(Side::south, i, j));
	//	return div;
	//}
	//static double divergence_finite_volume(Velocity& FX, int i)
	//{
	//	double div = 0.0;
	//	div += FX.Sx * (FX.get_for_centered_cell(Side::east, i) - FX.get_for_centered_cell(Side::west, i));		
	//	return div;
	//}

	int range_out(int i, int j = 0, int k = 0)
	{
		if (i < 0 || j < 0 || k < 0) return 1;
		if (type == Component::x && i >= nx + 1)
			return 2;
		else if (type != Component::x && i >= nx)
			return 3;

		if (type == Component::y && j >= ny + 1)
			return 4;
		else if (type != Component::y && j >= ny)
			return 5;

		if (type == Component::z && k >= nz + 1)
			return 6;
		else if (type != Component::z && k >= nz)
			return 7;

		return 0;

	}

	int index_shift(Component c, int i)
	{
	}

};


struct StaticVector
{
	double a = 0;
	double b = 0;
	double x, y, z;
	StaticVector(double a_ = 0, double b_ = 0) : a(a_), b(b_)
	{
		set_by_degree(a, b);
	}
	void set_by_degree(double a_ = 0, double b_ = 0)
	{
		double rad = acos(-1.0) / 180;
		a = a_;
		b = b_;

		x = cos(a * rad) * cos(b * rad);
		y = sin(a * rad) * cos(b * rad);
		z = sin(b * rad);
	}
	void set_directly_xyz(double x_, double y_ = 0, double z_ = 0)
	{
		x = x_;
		y = y_;
		z = z_;
		double l = sqrt(x * x + y * y + z * z);
		if (l == 0)
		{
			x = y = z = 0;
		}
		{
			x = x / l;
			y = y / l;
			z = z / l;
		}
	}

	double operator() (Component c) const
	{
		switch (c)
		{
		case Component::x:
			return x;
			break;
		case Component::y:
			return y;
			break;
		case Component::z:
			return z;
			break;
		default:
			return 0;
			break;
		}
	}
};