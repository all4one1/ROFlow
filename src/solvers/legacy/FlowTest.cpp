#pragma once
#include "FlowSolver.h"

#define DIRICHLET u.boundary.type(side) == MathBoundary::Dirichlet
#define PERIODIC u.boundary.type(side) == MathBoundary::Periodic


void FlowSolver::form_matrix_test(SparseMatrix& M, double* b)
{
	timer.start("matrix_formation");

	M.reset();
	struct Contribution
	{
		std::map<int, double> m;
		int l, i, j, k, nx, ny, nz, off, off2;
		int x1, y1, z1, stride;
		Contribution(ScalarVariable& F, int stride_, int l_, int i_ = 0, int j_ = 0, int k_ = 0)
			: stride(stride_), l(l_), i(i_), j(j_), k(k_)
		{
			off = F.off;
			off2 = F.off2;
			nx = F.nx;
			ny = F.ny;
			nz = F.nz;
			x1 = F.xExtra;
			y1 = F.yExtra;
			z1 = F.zExtra;
		}
		int index(Side side)
		{
			switch (side)
			{
			case Side::center:
				return l;
				break;
			case Side::west:
				return l - 1;
				break;
			case Side::east:
				return l + 1;
				break;
			case Side::south:
				return l - off;
				break;
			case Side::north:
				return l + off;
				break;
			case Side::front:
				return l - off2;
				break;
			case Side::back:
				return l + off2;
				break;
			default:
				return -1;
				break;
			}
			return -1;
		};

		int index_period(Side side)
		{
			switch (side)
			{
			case Side::center:
				return l;
				break;
			case Side::west:
				return (nx - 1) + off * j + off2 * k + stride;
				break;
			case Side::east:
				return x1 + off * j + off2 * k + stride;
				break;
			case Side::south:
				return i + (ny - 1) * off + off2 * k + stride;
				break;
			case Side::north:
				return i + off * y1 + off2 * k + stride;
				break;
			case Side::front:
				return i + j * off + (nz - 1) * off2 + stride;
				break;
			case Side::back:
				return i + j * off + z1 * off2 + stride;
				break;
			default:
				return -1;
				break;
			}
			return -1;
		};

		double& operator()(Side side, bool period = false)
		{
			if (period == false) return (*this)[index(side)];
			else return (*this)[index_period(side)];
		}
		double& operator[](int l)
		{
			return m[l];
		}
		map<int, double>& get_map()
		{
			return m;
		}
	};


	auto general = [](Contribution &a, Velocity& v, Side side, double S, double h, double coef, bool border = false)
	{
		if (border) //half step from a boundary
		{
			if (v.boundary.type(side) == MathBoundary::Dirichlet)
			{
				a(side) = 0.0;
				a(Side::center) += coef * S / (0.5 * h);
				//b[l] += v.boundary.get_value(side) * coef * S / (0.5 * h);
			}
			else if (v.boundary.type(side) == MathBoundary::Neumann)
			{ /*b[l] += v.boundary.normal_deriv_oriented(side) * S * coef;*/	}
			else if (v.boundary.type(side) == MathBoundary::Periodic)
			{
				a(Side::center) += coef * S / h;
				a(side, true) += -coef * S / h;
			}
		}
		else
		{
			a(Side::center) += coef * S / h;
			a(side) += -coef * S / h;
		}
	};
	auto fixed_node = [b](Velocity& v, Side side, int l, Contribution& a)
	{
		a[l] = 1;
		b[l] = v.boundary.get_fixed_value(side);
	};
	auto deriv_at_boundary = [b, &M](Velocity& v, Side side, int shift, double S, double h, double coef, int l, Contribution& a)
	{
		//if (v.boundary.type(side) == MathBoundary::Neumann)
		{
			a[l] = coef * S / h;
			a[l + shift] = -coef * S / h;;
			//b[l] += v.boundary.normal_deriv_oriented(side);
		}
	};
	auto period_boundary = [&M](Velocity& v, Side period_side, Side normal_side, double S, double h, double coef, int l, Contribution& a)
	{
		a(Side::center) += coef * S / h;
		a(period_side, true) = -coef * S / h;

		a(Side::center) += coef * S / h;
		a(normal_side) = -coef * S / h;
	};
	auto pre_boundary = [b, &M](Velocity& v, Side side, double S, double h, double PhysCoef, int l, Contribution& a)
	{
		if (v.boundary.type(side) == MathBoundary::Dirichlet)
		{
			a(Side::center) += PhysCoef * S / h;
		//	b[l] += v.boundary.get_fixed_value(side) * S / h * PhysCoef;
		}
	};
	auto temporal = [this](double dV, Contribution& a)
	{
		a(Side::center) += dV / tau;
		a(Side::center) /= alpha_relax;
	};

	if (1)
	{
		Velocity& u = ux;
		double PhysCoef = 1.0 / Re;
		if (u.type == Component::x)
		{
			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i <= nx; i++) {
						int l = i + u.off * j + u.off2 * k;
						Contribution a(u, 0, l, i, j, k);
						double reduced = 1.0;
						double S_reduced = 1.0;
						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;


						if ((i == 0 || i == nx)
							&& !(u.boundary.type(Side::west) == MathBoundary::Periodic))
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SY = SY * 0.5;
						}

						//Side::west
						if (i < nx)
						{
							Side side = Side::west;
							if (i == 0)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									fixed_node(u, side, l, a);
									M.add_line_with_map(a.get_map(), l);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, +1, SX, hx, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									Side opposide = Side::east;
									period_boundary(u, side, opposide, SX, hx, PhysCoef, l, a);
								}
							}
							else if (i == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SX, hx, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann
									|| u.boundary.type(side) == MathBoundary::Periodic)
								{
									general(a, u, side, SX, hx, PhysCoef, false);
								}
							}
							else
							{
								general(a, u, side, SX, hx, PhysCoef, false);
							}
						}

						//Side::east
						if (i > 0)
						{
							Side side = Side::east;
							if (i == nx)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									fixed_node(u, side, l, a);
									M.add_line_with_map(a.get_map(), l);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, -1, SX, hx, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									Side opposide = Side::west;
									period_boundary(u, side, opposide, SX, hx, PhysCoef, l, a);
								}
							}
							else if (i == nx - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SX, hx, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann
									|| u.boundary.type(side) == MathBoundary::Periodic)
								{
									general(a, u, side, SX, hx, PhysCoef, false);
								}
							}
							else
							{
								general(a, u, side, SX, hx, PhysCoef, false);
							}
						}

						if (dim > 1)
						{
							general(a, u, Side::south, SY, hy, PhysCoef, j == 0);
							general(a, u, Side::north, SY, hy, PhysCoef, j == ny - 1);
						}
						if (dim > 2)
						{
							general(a, u, Side::front, SZ, hz, PhysCoef, k == 0);
							general(a, u, Side::back, SZ, hz, PhysCoef, k == nz - 1);
						}

						temporal(DV, a);


						M.add_line_with_map(a.get_map(), l);
					}
				}
			}
		}
	}

	if (1 && dim > 1)
	{
		Velocity& u = uy;
		double PhysCoef = 1.0 / Re;
		if (u.type == Component::y)
		{
			int qq = 0;
			for (int k = 0; k < nz; k++) {
				for (int j = 0; j <= ny; j++) {
					for (int i = 0; i < nx; i++) {
						int l = i + u.off * j + u.off2 * k + stride;
						Contribution a(u, stride, l, i, j, k);
						double reduced = 1.0;
						double S_reduced = 1.0;
						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;

						if ((j == 0 || j == ny)
							&& !(u.boundary.type(Side::south) == MathBoundary::Periodic))
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SX = SX * 0.5;
						}

						//Side::south
						if (j < ny)
						{
							Side side = Side::south;
							if (j == 0)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									fixed_node(u, side, l, a);
									M.add_line_with_map(a.get_map(), l);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, +off, SY, hy, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									Side opposide = Side::north;
									period_boundary(u, side, opposide, SY, hy, PhysCoef, l, a);
								}
							}
							else if (j == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SY, hy, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann
									|| u.boundary.type(side) == MathBoundary::Periodic)
								{
									general(a, u, side, SY, hy, PhysCoef, false);
								}
							}
							else
							{
								general(a, u, side, SY, hy, PhysCoef, false);
							}
						}

						//Side::north
						if (j > 0)
						{
							Side side = Side::north;
							if (j == ny)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									fixed_node(u, side, l, a);
									M.add_line_with_map(a.get_map(), l);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, -off, SY, hy, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									Side opposide = Side::south;
									period_boundary(u, side, opposide, SY, hy, PhysCoef, l, a);
								}
							}
							else if (j == ny - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SY, hy, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann
									|| u.boundary.type(side) == MathBoundary::Periodic)
								{
									general(a, u, side, SY, hy, PhysCoef, false);
								}
							}
							else
							{
								general(a, u, side, SY, hy, PhysCoef, false);
							}
						}


						//if (dim > 1)
						{
							general(a, u, Side::west, SX, hx, PhysCoef, i == 0);
							general(a, u, Side::east, SX, hx, PhysCoef, i == nx - 1);
						}
						if (dim > 2)
						{
							general(a, u, Side::front, SZ, hz, PhysCoef, k == 0);
							general(a, u, Side::back, SZ, hz, PhysCoef, k == nz - 1);
						}

						temporal(DV, a);


						M.add_line_with_map(a.get_map(), l);
					}
				}
			}
		}
	}

	if (1 && dim > 2)
	{
		Velocity& u = uz;
		double PhysCoef = 1.0 / Re;
		if (u.type == Component::z)
		{
			for (int k = 0; k <= nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						int l = i + u.off * j + u.off2 * k + stride2;
						Contribution a(u, stride2, l, i, j, k);
						double reduced = 1.0;
						double S_reduced = 1.0;
						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;

						if ((k == 0 || k == nz)
							&& !(u.boundary.type(Side::front) == MathBoundary::Periodic))
						{
							DV = DV * 0.5;
							SY = SY * 0.5;
							SX = SX * 0.5;
						}


						//Side::front
						if (k < nz)
						{
							Side side = Side::front;
							if (k == 0)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									fixed_node(u, side, l, a);
									M.add_line_with_map(a.get_map(), l);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, +off2, SZ, hz, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									Side opposide = Side::back;
									period_boundary(u, side, opposide, SZ, hz, PhysCoef, l, a);
								}
							}
							else if (k == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SZ, hz, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann
									|| u.boundary.type(side) == MathBoundary::Periodic)
								{
									general(a, u, side, SZ, hz, PhysCoef, false);
								}
							}
							else
							{
								general(a, u, side, SZ, hz, PhysCoef, false);
							}
						}

						//Side::north
						if (k > 0)
						{
							Side side = Side::back;
							if (k == nz)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									fixed_node(u, side, l, a);
									M.add_line_with_map(a.get_map(), l);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, -off2, SZ, hz, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									Side opposide = Side::front;
									period_boundary(u, side, opposide, SZ, hz, PhysCoef, l, a);
								}
							}
							else if (k == nz - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SZ, hz, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann
									|| u.boundary.type(side) == MathBoundary::Periodic)
								{
									general(a, u, side, SZ, hz, PhysCoef, false);
								}
							}
							else
							{
								general(a, u, side, SZ, hz, PhysCoef, false);
							}
						}


						//if (dim > 1)
						{
							general(a, u, Side::west, SX, hx, PhysCoef, i == 0);
							general(a, u, Side::east, SX, hx, PhysCoef, i == nx - 1);
						}
						//if (dim > 2)
						{
							general(a, u, Side::south, SY, hy, PhysCoef, j == 0);
							general(a, u, Side::north, SY, hy, PhysCoef, j == ny - 1);
						}

						temporal(DV, a);


						M.add_line_with_map(a.get_map(), l);
					}
				}
			}
		}
	}


	timer.end("matrix_formation");
}
void FlowSolver::form_rhs_test(double* b, bool reset)
{
	double PhysCoef = 1.0 / Re;

	auto half_border = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
	{
		double res = 0.0;
		if (V.boundary.type(side) == MathBoundary::Neumann)
			res = V.boundary.normal_deriv_oriented(side) * S * coef;
		else if (V.boundary.type(side) == MathBoundary::Dirichlet)
			res = V.boundary.get_fixed_value(side) * coef * S / (0.5 * h);
		else if (V.boundary.type(side) == MathBoundary::Periodic)
		{	/*there should be nothing to do*/  }
		else
		{
			print("half_step error");
		}
		return res;
	};
	auto pre_fixed_node = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
	{
		double res = 0.0;
		if (V.boundary.type(side) == MathBoundary::Dirichlet)
			res = V.boundary.get_fixed_value(side) * coef * S / (h);
		//else if (V.boundary.type(side) == MathBoundary::Neumann)
		return res;
	};
	auto FluxVX = [this](int i, int j, int k, double Sx, double Sy, double Sz)
	{
		double FX = Sx * (vx.get_for_vx_cell(Side::east, i, j, k) * vx.get_for_vx_cell(Side::east, i, j, k)
						- vx.get_for_vx_cell(Side::west, i, j, k) * vx.get_for_vx_cell(Side::west, i, j, k));
		double FY = Sy * (vx.get_for_vx_cell(Side::north, i, j, k) * vy.get_for_vx_cell(Side::north, i, j, k)
						- vx.get_for_vx_cell(Side::south, i, j, k) * vy.get_for_vx_cell(Side::south, i, j, k));
		double FZ = Sz * (vx.get_for_vx_cell(Side::back, i, j, k) * vz.get_for_vx_cell(Side::back, i, j, k)
						- vx.get_for_vx_cell(Side::front, i, j, k) * vz.get_for_vx_cell(Side::front, i, j, k));
		return (FX + FY + FZ);
	};
	auto FluxVY = [this](int i, int j, int k, double Sx, double Sy, double Sz)
	{
		double FX = Sx * (vx.get_for_vy_cell(Side::east, i, j, k) * vy.get_for_vy_cell(Side::east, i, j, k)
						- vx.get_for_vy_cell(Side::west, i, j, k) * vy.get_for_vy_cell(Side::west, i, j, k));
		double FY = Sy * (vy.get_for_vy_cell(Side::north, i, j, k) * vy.get_for_vy_cell(Side::north, i, j, k)
						- vy.get_for_vy_cell(Side::south, i, j, k) * vy.get_for_vy_cell(Side::south, i, j, k));
		double FZ = Sz * (vy.get_for_vy_cell(Side::back, i, j, k) * vz.get_for_vy_cell(Side::back, i, j, k)
						- vy.get_for_vy_cell(Side::front, i, j, k) * vz.get_for_vy_cell(Side::front, i, j, k));
		return (FX + FY + FZ);
	};
	auto FluxVZ = [this](int i, int j, int k, double Sx, double Sy, double Sz)
	{
		double FX = Sx * (vx.get_for_vz_cell(Side::east, i, j, k) * vz.get_for_vz_cell(Side::east, i, j, k)
			- vx.get_for_vz_cell(Side::west, i, j, k) * vz.get_for_vz_cell(Side::west, i, j, k));
		double FY = Sy * (vy.get_for_vz_cell(Side::north, i, j, k) * vz.get_for_vz_cell(Side::north, i, j, k)
			- vy.get_for_vz_cell(Side::south, i, j, k) * vz.get_for_vz_cell(Side::south, i, j, k));
		double FZ = Sz * (vz.get_for_vz_cell(Side::back, i, j, k) * vz.get_for_vz_cell(Side::back, i, j, k)
			- vz.get_for_vz_cell(Side::front, i, j, k) * vz.get_for_vz_cell(Side::front, i, j, k));
		return (FX + FY + FZ);
	};
	auto general = [this, b](int i, int j, int k, Velocity& u, Velocity& v0, Component comp, double tau, double dV,  double S, double h)
		{
			double res = 0.0;
			double T_KC = T.get_shifted(comp, i, j, k) + K * C.get_shifted(comp, i, j, k);
			double dT_KdC = T.get_shifted_deriv(comp, i, j, k) + K * C.get_shifted_deriv(comp, i, j, k);
			res+= v0(i, j, k) * dV / tau
				- P.get_shifted_diff(comp, i, j, k) * S
				+ Ra / Pr * grav(comp) * T_KC * dV
				- Rav / Pr * bufferVibr.get_shifted(comp, i, j, k) * dT_KdC * dV;
			res += (1.0 - alpha_relax) / alpha_relax * Aii(u, i, j, k) * v0(i, j, k);
			
			return res;
		};


	if (1)
	{
		Velocity& u = ux;
		Velocity& v = vx;
		if (u.type == Component::x)
		{
			for (int k = 0; k < nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i <= nx; i++) {
						int l = i + u.off * j + u.off2 * k;
						if (reset) b[l] = 0;

						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;
						double rx = 1, ry = 1, rz = 1;


						if ((i == 0 || i == nx) 
							&& !(u.boundary.type(Side::west) == MathBoundary::Periodic))
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SY = SY * 0.5;
						}

						// west
						if (i < nx)
						{
							Side side = Side::west;
							if (i == 0)
							{
								//fixed node
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									b[l] = u.boundary.get_fixed_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SX;
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									// nothing to do
								}
							}
							else if (i == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet) //dont like it
									b[l] += half_border(u, l, side, SX, 2.0 * hx, PhysCoef);
							}
						}

						// east
						if (i > 0)
						{
							Side side = Side::east;
							if (i == nx)
							{
								//fixed node
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									b[l] = u.boundary.get_fixed_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SX;
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									// nothing to do
								}
							}
							else if (i == nx - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
									b[l] += half_border(u, l, side, SX, 2.0 * hx, PhysCoef);
							}
						}


						if (dim > 1)
						{
							//south
							if (j == 0)
							{
								Side side = Side::south;
								b[l] += half_border(u, l, side, SY, hy, PhysCoef);
							}

							//north
							if (j == ny - 1)
							{
								Side side = Side::north;
								b[l] += half_border(u, l, side, SY, hy, PhysCoef);
							}
						}


						if (dim > 2)
						{
							//front
							if (k == 0)
							{
								Side side = Side::front;
								b[l] += half_border(u, l, side, SZ, hz, PhysCoef);
							}
							//back
							if (k == nz - 1)
							{
								Side side = Side::back;
								b[l] += half_border(u, l, side, SZ, hz, PhysCoef);
							}
						}
						//rx = ry = rz = 1;
						b[l] += general(i, j, k, u, v, Component::x, tau, DV, SX, hx);
						b[l] += -FluxVX(i, j, k, SX * rx, SY * ry, SZ * rz);
					}
				}
			}
		}
	}

	if (1 && dim > 1)
	{
		Velocity& u = uy;
		Velocity& v = vy;
		if (u.type == Component::y)
		{
			for (int k = 0; k < nz; k++) {
				for (int j = 0; j <= ny; j++) {
					for (int i = 0; i < nx; i++) {
						int l = i + u.off * j + u.off2 * k + stride;
						if (reset) b[l] = 0;

						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;
						double rx = 1, ry = 1, rz = 1;


						if ((j == 0 || j == ny)
							&& !(u.boundary.type(Side::south) == MathBoundary::Periodic))
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SX = SX * 0.5;
						}


						// south
						if (j < ny)
						{
							Side side = Side::south;
							if (j == 0)
							{
								//fixed node
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									b[l] = u.boundary.get_fixed_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SY;
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									// nothing to do
								}
							}
							else if (j == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
									b[l] += half_border(u, l, side, SY, 2.0 * hy, PhysCoef);
							}
						}

						// north
						if (j > 0)
						{
							Side side = Side::north;
							if (j == ny)
							{
								//fixed node
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									b[l] = u.boundary.get_fixed_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SY;
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									// nothing to do
								}
							}
							else if (j == ny - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
									b[l] += half_border(u, l, side, SY, 2.0 * hy, PhysCoef);
							}
						}


						//if (dim > 1)
						{
							//west
							if (i == 0)
							{
								Side side = Side::west;
								b[l] += half_border(u, l, side, SX, hx, PhysCoef);
							}

							//east
							if (i == nx - 1)
							{
								Side side = Side::east;
								b[l] += half_border(u, l, side, SX, hx, PhysCoef);
							}
						}


						if (dim > 2)
						{
							//front
							if (k == 0)
							{
								Side side = Side::front;
								b[l] += half_border(u, l, side, SZ, hz, PhysCoef);
							}
							//back
							if (k == nz - 1)
							{
								Side side = Side::back;
								b[l] += half_border(u, l, side, SZ, hz, PhysCoef);
							}
						}
						//rx = ry = rz = 1;
						b[l] += general(i, j, k, u, v, Component::y, tau, DV, SY, hy);
						b[l] += -FluxVY(i, j, k, SX * rx, SY * ry, SZ * rz);
					}
				}
			}
		}
	}

	if (1 && dim > 2)
	{
		Velocity& u = uz;
		Velocity& v = vz;
		if (u.type == Component::z)
		{
			for (int k = 0; k <= nz; k++) {
				for (int j = 0; j < ny; j++) {
					for (int i = 0; i < nx; i++) {
						int l = i + u.off * j + u.off2 * k + stride2;
						if (reset) b[l] = 0;

						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;
						double rx = 1, ry = 1, rz = 1;

						if ((k == 0 || k == nz)
							&& !(u.boundary.type(Side::front) == MathBoundary::Periodic))
						{
							DV = DV * 0.5;
							SY = SY * 0.5;
							SX = SX * 0.5;
						}


						// front
						if (k < nz)
						{
							Side side = Side::front;
							if (k == 0)
							{
								//fixed node
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									b[l] = u.boundary.get_fixed_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SZ;
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									// nothing to do
								}
							}
							else if (k == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
									b[l] += half_border(u, l, side, SZ, 2.0 * hz, PhysCoef);
							}
						}

						// back
						if (k > 0)
						{
							Side side = Side::back;
							if (k == nz)
							{
								//fixed node
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									b[l] = u.boundary.get_fixed_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SZ;
									//continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Periodic)
								{
									// nothing to do
								}
							}
							else if (k == nz - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
									b[l] += half_border(u, l, side, SZ, 2.0 * hz, PhysCoef);
							}
						}


						//if (dim > 1)
						{
							//west
							if (i == 0)
							{
								Side side = Side::west;
								b[l] += half_border(u, l, side, SX, hx, PhysCoef);
							}

							//east
							if (i == nx - 1)
							{
								Side side = Side::east;
								b[l] += half_border(u, l, side, SX, hx, PhysCoef);
							}
						}


						//if (dim > 2)
						{
							//front
							if (j == 0)
							{
								Side side = Side::south;
								b[l] += half_border(u, l, side, SY, hy, PhysCoef);
							}
							//back
							if (j == ny - 1)
							{
								Side side = Side::north;
								b[l] += half_border(u, l, side, SY, hy, PhysCoef);
							}
						}
						//rx = ry = rz = 1;
						b[l] += general(i, j, k, u, v, Component::z, tau, DV, SZ, hz);
						b[l] += -FluxVZ(i, j, k, SX * rx, SY * ry, SZ * rz);
					}
				}
			}
		}
	}
}

void FlowSolver::correction_for_simple_test()
{
	timer.start("correction");
	double alpha = alpha_relax;
	double max = 0;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				P(i, j, k) = P(i, j, k) + alpha * p_prime(i, j, k);
			}
		}
	}

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i <= nx; i++) {
				int l = i + ux.off * j + ux.off2 * k;
				double v_prime = -(p_prime.get_shifted_diff(Component::x, i, j, k)) * Sx / SM(l, l);
				if (abs(v_prime) > max) max = v_prime;
				ux(i, j, k) = alpha * ux(i, j, k) + (1.0 - alpha) * vx(i, j, k) + v_prime;
			}
		}
	}

	if (1 && dim > 1)
	{
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j <= ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + uy.off * j + uy.off2 * k + stride;
					double v_prime = -(p_prime.get_shifted_diff(Component::y, i, j, k)) * Sy / SM(l, l);
					if (abs(v_prime) > max) max = v_prime;
					uy(i, j, k) = alpha * uy(i, j, k) + (1.0 - alpha) * vy(i, j, k) + v_prime;
				}
			}
		}
	}

	if (1 && dim > 2)
	{
		for (int k = 0; k <= nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + uz.off * j + uz.off2 * k + stride2;
					double v_prime = -(p_prime.get_shifted_diff(Component::z, i, j, k)) * Sz / SM(l, l);
					if (abs(v_prime) > max) max = v_prime;
					uz(i, j, k) = alpha * uz(i, j, k) + (1.0 - alpha) * vz(i, j, k) + v_prime;
				}
			}
		}
	}
	test_var = max;

	timer.end("correction");
}
