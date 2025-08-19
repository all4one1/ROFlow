#pragma once
#include <iostream>
#include <fstream>
#include <map>

#include "tools/types.h"
#pragma warning(disable: 4018)

#define GETVALUE(i, j, k) v[(i) + offset * (j) + offset2 * (k)]

#define GETFX(i2,i1) (GETVALUE((i2), j, k) + GETVALUE((i1), j, k)) * 0.5
#define GETFY(j2,j1_) (GETVALUE(i, (j2), k) + GETVALUE(i, (j1_), k)) * 0.5
#define GETFZ(k2,k1) (GETVALUE(i, j, (k2)) + GETVALUE(i, j, (k1))) * 0.5;

#define GETDX(i2,i1) (GETVALUE((i2), j, k) - GETVALUE((i1), j, k)) / hx
#define GETDY(j2,j1_) (GETVALUE(i, (j2), k) - GETVALUE(i, (j1_), k)) / hy
#define GETDZ(k2,k1) (GETVALUE(i, j, (k2)) - GETVALUE(i, j, (k1))) / hz

#define DIFFX(i2,i1) (GETVALUE((i2), j, k) - GETVALUE((i1), j, k))
#define DIFFY(j2,j1_) (GETVALUE(i, (j2), k) - GETVALUE(i, (j1_), k))
#define DIFFZ(k2,k1) (GETVALUE(i, j, (k2)) - GETVALUE(i, j, (k1)))

namespace newvar
{
	struct Boundary
	{
		double value = 0.0;
		MathBoundary type = MathBoundary::Neumann;
		double* v = nullptr;
		Boundary() {};
		Boundary(MathBoundary t, double v = 0.0) : type(t), value(v) {};
		Boundary(MathBoundary t, bool new_array = false, unsigned int n = 0, double fixed_val = NAN) : type(t)
		{
			if (new_array) v = new double[n];
			if (fixed_val == fixed_val)
			{
				for (unsigned int l = 0; l < n; l++)
					v[l] = fixed_val;
			}
		};
		//Boundary(MathBoundary t, bool new_array, unsigned int n, double (*func)(double)) : type(t)
		double& operator()(unsigned int q) { return v[q]; }
	};
	using walls = std::vector<Boundary>;
	struct Variable
	{
		double* f = nullptr;
		walls wall;

		Variable()
		{
			wall.resize(4);
		};
		Variable(double* ptr) : f(ptr)
		{
			wall.resize(4);
		}
		double& operator[](int i) { return f[i]; }
		operator double* () const { return f; }

		void set_walls(Boundary west_, Boundary east_, Boundary south_, Boundary north_)
		{
			wall[west] = west_;
			wall[east] = east_;
			wall[south] = south_;
			wall[north] = north_;
		}
	};
}

//namespace newvar2
//{
	struct Boundary
	{
		struct OneWall
		{
			double value;
			double sign;
			double h;
			bool shifted = false;
			MathBoundary type;
			OneWall(double v = 0, double s = 0, double H = 0, MathBoundary t = MathBoundary::Neumann)
				: value(v), sign(s), h(H), type(t) {
			}
	
		};
		std::map<Side, OneWall> wall; //via vector?
		int dim = 1;
		double* v = nullptr;
		int nx, ny, nz, offset, offset2;
	
		Boundary(MathBoundary def_type = MathBoundary::Neumann, double def_value = 0.0,
			double hx = 0, double hy = 0, double hz = 0)
		{
			wall[Side::west] = OneWall(def_value, -1, hx, def_type);
			wall[Side::east] = OneWall(def_value, +1, hx, def_type);
			wall[Side::south] = OneWall(def_value, -1, hy, def_type);
			wall[Side::north] = OneWall(def_value, +1, hy, def_type);
			wall[Side::front] = OneWall(def_value, -1, hz, def_type);
			wall[Side::back] = OneWall(def_value, +1, hz, def_type);
			nx = ny = nz = offset = offset2 = 0;
		}
		void shifted_on(Side s)
		{
			wall[s].shifted = true;
		}
	
		void bind_with_field(double* p, int nx_, int ny_, int nz_, int off_, int off2_)
		{
			v = p;
			nx = nx_;
			ny = ny_;
			nz = nz_;
			offset = off_;
			offset2 = off2_;
		}
		auto iterator(Side side)
		{
			return &wall.find(side)->second;
		}
		MathBoundary type(Side side)
		{
			return wall[side].type;
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
					double h1 = wall[s1].h, h2 = wall[s2].h;
					double h = h1 + h2;
					set_boundary(s1, MathBoundary::periodic, 0.0, h);
					set_boundary(s2, MathBoundary::periodic, 0.0, h);
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
			wall[Side::west].h = wall[Side::east].h = hx_;
			wall[Side::south].h = wall[Side::north].h = hy_;
			wall[Side::front].h = wall[Side::back].h = hz_;
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
			if (it->type == MathBoundary::periodic)
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
			return normal_deriv(side, i, j, k) * wall[side].h;
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
				return wall[side].value;
			else if (it->type == MathBoundary::periodic)
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
			return NAN;
		}
		double get_value_from_deriv(Side side, int i = 0, int j = 0, int k = 0)
		{
			if (wall[side].type == MathBoundary::Neumann)
				return retrieve_from_deriv(side, GETVALUE(i, j, k));
			else
				cout << "incorrect bc type" << endl;
		}
		double get_fixed_value(Side side)
		{
			if (wall[side].type == MathBoundary::Dirichlet)
				return wall[side].value;
			//else if (wall[side].type == MathBoundary::periodic)
			else
			{
				cout << "incorrect bc type" << endl;
			}
		}
	
		OneWall& operator[](Side side)
		{
			return wall[side];
		}
	};

	struct Variable : Configuration
	{
		double* v = nullptr;
		int xExtra = 0, yExtra = 0, zExtra = 0;
		Boundary boundary;

		Variable() : Configuration()
		{

		}

		Variable(Configuration &config, bool newPointer = false, double* ptr = nullptr) : Configuration(config)
		{
			default_(newPointer, ptr);
		}

		void default_(bool newPointer, double* ptr)
		{
			boundary.set_boundary(Side::west, MathBoundary::Neumann, 0.0);
			boundary.set_boundary(Side::east, MathBoundary::Neumann, 0.0);
			boundary.set_boundary(Side::south, MathBoundary::Neumann, 0.0);
			boundary.set_boundary(Side::north, MathBoundary::Neumann, 0.0);
			boundary.set_boundary(Side::front, MathBoundary::Neumann, 0.0);
			boundary.set_boundary(Side::back, MathBoundary::Neumann, 0.0);

			boundary.set_steps(0.5 * hx, 0.5 * hy, 0.5 * hz);


			if (newPointer)
			{
				v = new double[N];
				for (unsigned int l = 0; l < N; l++)
					v[l] = 0.0;
			}
			else if (!newPointer)
			{
				if (ptr != nullptr)
					v = ptr;
			}

			boundary.bind_with_field(v, nx, ny, nz, offset, offset2);
		}


		double& operator()(int i, int j = 0, int k = 0)	{return v[i + offset * j + offset2 * k];}
		double& operator[](int l) {return v[l];}
		operator double* () const { return v; }
		double* get_ptr() {return v;}




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
					if (i == 0 && boundary.type(Side::west) == MathBoundary::periodic)
						return (boundary(side, 1, j, k) - (boundary(side, nx - 1, j, k))) / (2 * hx);
					if (i == nx - 1 && boundary.type(Side::east) == MathBoundary::periodic)
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
					double f1 = 0.5 * (get_value(Side::west, i, j, k) + get_value(Side::west, i, j - 1, k));
					double f2 = 0.5 * (get_value(Side::east, i, j, k) + get_value(Side::east, i, j - 1, k));
					return (f2 - f1) / hx;
				}
			}

			if (side == Side::north)
			{
				if (j == ny - 1)
				{
					if (i == 0 && boundary.type(Side::west) == MathBoundary::periodic)
						return (boundary(side, 1, j, k) - (boundary(side, nx - 1, j, k))) / (2 * hx);
					if (i == nx - 1 && boundary.type(Side::east) == MathBoundary::periodic)
						return (boundary(side, 0, j, k) - (boundary(side, nx - 2, j, k))) / (2 * hx);
					if (i == 0)				return (boundary(side, i + 1, j, k) - (boundary(side, i, j, k))) / (hx);
					else if (i == nx - 1)	return (boundary(side, i, j, k) - (boundary(side, i - 1, j, k))) / (hx);
					else                    return (boundary(side, i + 1, j, k) - (boundary(side, i - 1, j, k))) / (2 * hx);
				}
				else
				{
					double f1 = 0.5 * (get_value(Side::west, i, j, k) + get_value(Side::west, i, j + 1, k));
					double f2 = 0.5 * (get_value(Side::east, i, j, k) + get_value(Side::east, i, j + 1, k));
					return (f2 - f1) / hx;
				}
			}


			if (side == Side::center)
			{
				return (get_value(Side::east, i, j, k) - get_value(Side::west, i, j, k)) / (hx);
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
					if (j == 0 && boundary.type(Side::south) == MathBoundary::periodic)
						return (boundary(side, i, 1, k) - (boundary(side, i, ny - 1, k))) / (2 * hy);
					if (j == ny - 1 && boundary.type(Side::north) == MathBoundary::periodic)
						return (boundary(side, i, 0, k) - (boundary(side, i, ny - 2, k))) / (2 * hy);

					if (j == 0)				return (boundary(side, i, j + 1, k) - (boundary(side, i, j, k))) / (hy);
					else if (j == ny - 1)	return (boundary(side, i, j, k) - (boundary(side, i, j - 1, k))) / (hy);
					else                    return (boundary(side, i, j + 1, k) - (boundary(side, i, j - 1, k))) / (2 * hy);
				}
				else
				{
					double f1 = 0.5 * (get_value(Side::south, i, j, k) + get_value(Side::south, i - 1, j, k));
					double f2 = 0.5 * (get_value(Side::north, i, j, k) + get_value(Side::north, i - 1, j, k));
					return (f2 - f1) / hy;
				}
			}

			if (side == Side::east)
			{
				if (i == nx - 1)
				{
					if (j == 0 && boundary.type(Side::south) == MathBoundary::periodic)
						return (boundary(side, i, 1, k) - (boundary(side, i, ny - 1, k))) / (2 * hy);
					if (j == ny - 1 && boundary.type(Side::north) == MathBoundary::periodic)
						return (boundary(side, i, 0, k) - (boundary(side, i, ny - 2, k))) / (2 * hy);

					if (j == 0)				return (boundary(side, i, j + 1, k) - (boundary(side, i, j, k))) / (hy);
					else if (j == ny - 1)	return (boundary(side, i, j, k) - (boundary(side, i, j - 1, k))) / (hy);
					else                    return (boundary(side, i, j + 1, k) - (boundary(side, i, j - 1, k))) / (2 * hy);
				}
				else
				{
					double f1 = 0.5 * (get_value(Side::south, i, j, k) + get_value(Side::south, i + 1, j, k));
					double f2 = 0.5 * (get_value(Side::north, i, j, k) + get_value(Side::north, i + 1, j, k));
					return (f2 - f1) / hy;
				}
			}


			if (side == Side::center)
			{
				return (get_value(Side::north, i, j, k) - get_value(Side::south, i, j, k)) / (hy);
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
				return (get_value(Side::back, i, j, k) - get_value(Side::front, i, j, k)) / hz;
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
				return (get_value(Side::east, i, j, k) - get_value(Side::west, i, j, k));
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
				return (get_value(Side::north, i, j, k) - get_value(Side::south, i, j, k));
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
				return (get_value(Side::back, i, j, k) - get_value(Side::front, i, j, k));
			}

			throw std::invalid_argument("wrong side");
			return NAN;
		}

		double get_value(Side side, int i, int j = 0, int k = 0)
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
				break;
			}
		}

		//used within vec variables
		double get_shifted_value(Component comp, int i, int j = 0, int k = 0)
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

		double laplace_finite_volume(int i, int j, int k)
		{
			double lapl = 0.0;
			lapl += Sx * (get_dx(Side::east, i, j, k) - get_dx(Side::west, i, j, k));
			lapl += Sy * (get_dy(Side::north, i, j, k) - get_dy(Side::south, i, j, k));
			lapl += Sz * (get_dz(Side::back, i, j, k) - get_dz(Side::front, i, j, k));
			return lapl;
		}


		void reset_to_null()
		{
			for (int l = 0; l < N; l++)
				v[l] = 0.0;
		}

		static void swap(Variable& V1, Variable& V2)
		{
			double* temp;
			temp = V1.v;
			V1.v = V2.v;
			V2.v = temp;
		}

		void transfer_data_to(Variable& target)
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
			return i * hx + ox;
		}
		double y_(int j)
		{
			return j * hy + oy;
		}
		double z_(int k)
		{
			return k * hz + oz;
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
			boundary.set_boundary(Side::west, type, val);
			boundary.set_boundary(Side::east, type, val);
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
			cout << "offset = " << offset << ", offset2 = " << offset2 << endl;
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

//}

