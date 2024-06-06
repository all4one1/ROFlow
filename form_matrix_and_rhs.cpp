#include "FlowSolver.h"

void FlowSolver::form_big_matrix(SparseMatrix& M, double* b)
{
	timer.start("matrix_formation");

	M.reset();
	struct Contribution
	{
		std::map<int, double> m;
		int l, off, off2;
		Contribution(int l_, int off_ = 0, int off2_ = 0) : l(l_), off(off_), off2(off2_)
		{

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

		double& operator()(Side side)
		{
			return (*this)[index(side)];
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

	auto fixed_node = [b](Velocity& v, Side side, int l, Contribution& a)
		{
			a[l] = 1;
			b[l] = v.boundary.get_fixed_value(side);
		};

	auto deriv_at_boundary = [b, &M](Velocity& v, Side side, int shift, double S, double h, double coef, int l, Contribution& a)
		{
			if (v.boundary.type(side) == MathBoundary::Neumann)
			{
				a[l] = coef * S / h;
				a[l + shift] = -coef * S / h;;
				b[l] += v.boundary.normal_deriv_oriented(side);
			}
		};

	auto pre_boundary = [b, &M](Velocity& v, Side side, double S, double h, double PhysCoef, int l, Contribution& a)
		{
			if (v.boundary.type(side) == MathBoundary::Dirichlet)
			{
				a(Side::center) += PhysCoef * S / h;
				b[l] += v.boundary.get_fixed_value(side) * S / h * PhysCoef;
			}
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
						Contribution a(l, u.off, u.off2);
						double reduced = 1.0;
						double S_reduced = 1.0;
						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;



						auto general = [b, &a, i, j, k, l](Velocity& v, Side side, double S, double h, double coef, bool border = false)
							{
								if (border) //half step from a boundary
								{
									if (v.boundary.type(side) == MathBoundary::Dirichlet)
									{
										a(side) = 0.0;
										a(Side::center) += coef * S / (0.5 * h);
										//b[l] += v.boundary.get_value(side) * coef * S / (0.5 * h);
									}
									if (v.boundary.type(side) == MathBoundary::Neumann)
										b[l] += v.boundary.normal_deriv_oriented(side) * S * coef;
								}
								else
								{
									a(Side::center) += coef * S / h;
									a(side) += -coef * S / h;
								}
							};


						if (i == 0 || i == nx)
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SY = SY * 0.5;
						}

						if (REDUCED)
						{


							if (dim > 1)
							{
								if ((j == 0 && u.boundary.type(Side::south) == MathBoundary::Dirichlet)
									|| (j == (ny - 1) && u.boundary.type(Side::north) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SX = SX * 0.75;
									SZ = SZ * 0.75;
								}
							}
							if (dim > 2)
							{
								if ((k == 0 && u.boundary.type(Side::front) == MathBoundary::Dirichlet)
									|| (k == (nz - 1) && u.boundary.type(Side::back) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SX = SX * 0.75;
									SY = SY * 0.75;
								}
							}
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
								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, +1, SX, hx, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
							}
							else if (i == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SX, hx, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									general(u, side, SX, hx, PhysCoef, false);
								}
							}
							else
							{
								general(u, side, SX, hx, PhysCoef, false);
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

								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, -1, SX, hx, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
							}
							else if (i == nx - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SX, hx, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									general(u, side, SX, hx, PhysCoef, false);
								}
							}
							else
							{
								general(u, side, SX, hx, PhysCoef, false);
							}
						}

						if (dim > 1)
						{
							general(u, Side::south, SY, hy, PhysCoef, j == 0);
							general(u, Side::north, SY, hy, PhysCoef, j == ny - 1);
						}
						if (dim > 2)
						{
							general(u, Side::front, SZ, hz, PhysCoef, k == 0);
							general(u, Side::back, SZ, hz, PhysCoef, k == nz - 1);
						}

						a(Side::center) += DV / tau;


						M.add_line_with_map(a.get_map(), l);
					}
				}
			}
		}
	}


	if (1)
	{
		Velocity& u = uy;
		double PhysCoef = 1.0 / Re;
		if (u.type == Component::y)
		{
			for (int k = 0; k < nz; k++) {
				for (int j = 0; j <= ny; j++) {
					for (int i = 0; i < nx; i++) {
						int l = i + u.off * j + u.off2 * k + stride;
						Contribution a(l, u.off, u.off2);
						double reduced = 1.0;
						double S_reduced = 1.0;
						double SX = this->Sx;
						double SY = this->Sy;
						double SZ = this->Sz;
						double DV = this->dV;


						auto general = [b, &a, i, j, k, l](Velocity& v, Side side, double S, double h, double coef, bool border = false)
							{
								if (border) //half step from a boundary
								{
									if (v.boundary.type(side) == MathBoundary::Dirichlet)
									{
										a(side) = 0.0;
										a(Side::center) += coef * S / (0.5 * h);
										//b[l] += v.boundary.get_value(side) * coef * S / (0.5 * h);
									}
									if (v.boundary.type(side) == MathBoundary::Neumann)
										b[l] += v.boundary.normal_deriv_oriented(side) * S * coef;
								}
								else
								{
									a(Side::center) += coef * S / h;
									a(side) += -coef * S / h;
								}
							};

						if (j == 0 || j == ny)
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SX = SX * 0.5;
						}

						if (REDUCED)
						{


							//if (dim > 1)
							{
								if ((i == 0 && u.boundary.type(Side::west) == MathBoundary::Dirichlet)
									|| (i == (nx - 1) && u.boundary.type(Side::east) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SY = SY * 0.75;
									SZ = SZ * 0.75;
								}
							}
							if (dim > 2)
							{
								if ((k == 0 && u.boundary.type(Side::front) == MathBoundary::Dirichlet)
									|| (k == (nz - 1) && u.boundary.type(Side::back) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SX = SX * 0.75;
									SY = SY * 0.75;
								}
							}
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
								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, +off, SY, hy, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
							}
							else if (j == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SY, hy, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									general(u, side, SY, hy, PhysCoef, false);
								}
							}
							else
							{
								general(u, side, SY, hy, PhysCoef, false);
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

								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									deriv_at_boundary(u, side, -off, SY, hy, PhysCoef, l, a);
									//M.add_line_with_map(a.get_map(), l);
									//continue;
								}
							}
							else if (j == ny - 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
								{
									pre_boundary(u, side, SY, hy, PhysCoef, l, a);
								}
								if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									general(u, side, SY, hy, PhysCoef, false);
								}
							}
							else
							{
								general(u, side, SY, hy, PhysCoef, false);
							}
						}


						//if (dim > 1)
						{
							general(u, Side::west, SX, hx, PhysCoef, i == 0);
							general(u, Side::east, SX, hx, PhysCoef, i == nx - 1);
						}
						if (dim > 2)
						{
							general(u, Side::front, SZ, hz, PhysCoef, k == 0);
							general(u, Side::back, SZ, hz, PhysCoef, k == nz - 1);
						}

						a(Side::center) += DV / tau;


						M.add_line_with_map(a.get_map(), l);
					}
				}
			}
		}
	}

	timer.end("matrix_formation");
}

void FlowSolver::form_big_rhs(double* b, bool reset)
{
	auto half_border = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
		{
			double res = 0.0;
			if (V.boundary.type(side) == MathBoundary::Neumann)
				res = V.boundary.normal_deriv_oriented(side) * S * coef;
			else if (V.boundary.type(side) == MathBoundary::Dirichlet)
			{
				res = V.boundary.get_fixed_value(side) * coef * S / (0.5 * h);
			}
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
			return FX + FY + FZ;
		};
	auto FluxVY = [this](int i, int j, int k, double Sx, double Sy, double Sz)
		{
			double FX = Sx * (vx.get_for_vy_cell(Side::east, i, j, k) * vy.get_for_vy_cell(Side::east, i, j, k)
				- vx.get_for_vy_cell(Side::west, i, j, k) * vy.get_for_vy_cell(Side::west, i, j, k));
			double FY = Sy * (vy.get_for_vy_cell(Side::north, i, j, k) * vy.get_for_vy_cell(Side::north, i, j, k)
				- vy.get_for_vy_cell(Side::south, i, j, k) * vy.get_for_vy_cell(Side::south, i, j, k));
			double FZ = Sz * (vy.get_for_vy_cell(Side::back, i, j, k) * vz.get_for_vy_cell(Side::back, i, j, k)
				- vy.get_for_vy_cell(Side::front, i, j, k) * vz.get_for_vy_cell(Side::front, i, j, k));
			return FX + FY + FZ;
		};
	auto FluxVZ = [this](int i, int j, int k, double Sx, double Sy, double Sz)
		{
			double FX = Sx * (vx.get_for_vz_cell(Side::east, i, j, k) * vz.get_for_vz_cell(Side::east, i, j, k)
				- vx.get_for_vz_cell(Side::west, i, j, k) * vz.get_for_vz_cell(Side::west, i, j, k));
			double FY = Sy * (vy.get_for_vz_cell(Side::north, i, j, k) * vz.get_for_vz_cell(Side::north, i, j, k)
				- vy.get_for_vz_cell(Side::south, i, j, k) * vz.get_for_vz_cell(Side::south, i, j, k));
			double FZ = Sz * (vz.get_for_vz_cell(Side::back, i, j, k) * vz.get_for_vz_cell(Side::back, i, j, k)
				- vz.get_for_vz_cell(Side::front, i, j, k) * vz.get_for_vz_cell(Side::front, i, j, k));
			return FX + FY + FZ;
		};


	double PhysCoef = 1.0 / Re;
	double PhysCoef2 = Ra / Pr;
	double PhysCoef3 = Rav / Pr;
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

						auto general = [this, b, i, j, k, l, PhysCoef2, PhysCoef3](Velocity& v, Velocity& v0, Side side, double dp, double dV, double tau, double S, double h, double coef, double reduced = 1)
							{
								Component comp = Component::x;
								double T_KC = T.get_shifted(comp, i, j, k) + K * C.get_shifted(comp, i, j, k);
								double dT_KdC = T.get_shifted_deriv(comp, i, j, k) + K * C.get_shifted_deriv(comp, i, j, k);
								return v0(i, j, k) * dV / tau
									- P.get_shifted_diff(comp, i, j, k) * S
									+ PhysCoef2 * grav.x * T_KC * dV
									- PhysCoef3 * bufferVibr.get_shifted(comp, i, j, k) * dT_KdC * dV;
							};

						if (i == 0 || i == nx)
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SY = SY * 0.5;
						}

						if (REDUCED)
						{
							if (dim > 1)
							{
								if ((j == 0 && u.boundary.type(Side::south) == MathBoundary::Dirichlet)
									|| (j == (ny - 1) && u.boundary.type(Side::north) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SX = SX * 0.75;
									SZ = SZ * 0.75;
									ry = 0.75;
								}
							}
							if (dim > 2)
							{
								if ((k == 0 && u.boundary.type(Side::front) == MathBoundary::Dirichlet)
									|| (k == (nz - 1) && u.boundary.type(Side::back) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SX = SX * 0.75;
									SY = SY * 0.75;
									rz = 0.75;
								}
							}
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
									// nothin to do
								}
							}
							else if (i == 1)
							{
								if (u.boundary.type(side) == MathBoundary::Dirichlet)
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
						b[l] += general(u, v, Side::center, 1, DV, tau, SX, hx, PhysCoef, 1.0);
						b[l] += -FluxVX(i, j, k, SX * rx, SY * ry, SZ * rz);
					}
				}
			}
		}
	}

	if (1)
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

						auto general = [this, b, i, j, k, l, PhysCoef2, PhysCoef3](Velocity& v, Velocity& v0, Side side, double dp, double dV, double tau, double S, double h, double coef, double reduced = 1)
							{
								Component comp = Component::y;
								double T_KC = T.get_shifted(comp, i, j, k) + K * C.get_shifted(comp, i, j, k);
								double dT_KdC = T.get_shifted_deriv(comp, i, j, k) + K * C.get_shifted_deriv(comp, i, j, k);
								return v0(i, j, k) * dV / tau
									- P.get_shifted_diff(comp, i, j, k) * S
									+ PhysCoef2 * grav.y * T_KC * dV
									- PhysCoef3 * bufferVibr.get_shifted(comp, i, j, k) * dT_KdC * dV;
							};

						if (j == 0 || j == ny)
						{
							DV = DV * 0.5;
							SZ = SZ * 0.5;
							SX = SX * 0.5;
						}

						if (REDUCED)
						{

							//if (dim > 1)
							{
								if ((i == 0 && u.boundary.type(Side::west) == MathBoundary::Dirichlet)
									|| (i == (nx - 1) && u.boundary.type(Side::east) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SY = SY * 0.75;
									SZ = SZ * 0.75;
									rx = 0.75;
								}
							}
							if (dim > 2)
							{
								if ((k == 0 && u.boundary.type(Side::front) == MathBoundary::Dirichlet)
									|| (k == (nz - 1) && u.boundary.type(Side::back) == MathBoundary::Dirichlet))
								{
									DV = DV * 0.75;
									SX = SX * 0.75;
									SY = SY * 0.75;
									rz = 0.75;
								}
							}
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
						b[l] += general(u, v, Side::center, 1, DV, tau, SY, hy, PhysCoef, 1);
						b[l] += -FluxVY(i, j, k, SX * rx, SY * ry, SZ * rz);
					}
				}
			}
		}
	}
}