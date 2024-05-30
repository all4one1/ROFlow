#include "FlowSolver.h"


void FlowSolver::form_matrix_test(SparseMatrix& M, double* b)
{

	struct Contribution2
	{
		std::map<int, double> m;
		int l, off, off2;
		Contribution2(int l_, int off_ = 0, int off2_ = 0) : l(l_), off(off_), off2(off2_)
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

	if (ux.type == Component::x)
	{
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i <= nx; i++) {
					int l = i + ux.off * j + ux.off2 * k;
					Contribution2 a(l, ux.off, ux.off2);
					double reduced = 1.0;
					double S_reduced = 1.0;
					double SX = this->Sx;
					double SY = this->Sy;
					double SZ = this->Sz;
					double DV = this->dV;

					auto reduce = [](Velocity& V, Side side, double& reduced)
					{
						if (V.boundary.type(side) == MathBoundary::Dirichlet)
						{
							reduced *= 0.75;
						}
					};
					double PhysCoef = 1.0 / Re;

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
					auto fixed_node = [&a, l, b](Velocity& v, Side side)
					{
						a[l] = 1;
						b[l] = v.boundary.get_value(side);
					};
					auto deriv_at_boundary = [&a, l, b, &M](Velocity& v, Side side, int shift, double S, double h, double coef)
					{
						if (v.boundary.type(side) == MathBoundary::Neumann)
						{
							a[l] = coef * S / h;
							a[l + shift] = -coef * S / h;;
							b[l] += v.boundary.normal_deriv_oriented(side);
						}
					};
					auto pre_boundary = [&a, l, b, &M, general](Velocity& v, Side side, double S, double h, double PhysCoef)
					{
						if (v.boundary.type(side) == MathBoundary::Dirichlet)
						{
							a(Side::center) += PhysCoef * S / h;
							b[l] += v.boundary.get_value(side) * S / h * PhysCoef;
						}
						else if (v.boundary.type(side) == MathBoundary::Neumann)
						{
							general(v, side, S, h, PhysCoef, false);
						}
						else
						{
							print("wrong pre_boundary");
						}

					};


					//Side::west
					if (i < nx)
					{
						Side side = Side::west;
						if (i == 0)
						{
							if (ux.boundary.type(side) == MathBoundary::Dirichlet)
							{
								fixed_node(ux, side);
								M.add_line_with_map(a.get_map(), l);
								continue;
							}
							if (ux.boundary.type(side) == MathBoundary::Neumann)
							{
								deriv_at_boundary(ux, side, +1, Sx, hx, PhysCoef);
								//M.add_line_with_map(a.get_map(), l);
								//continue;
							}
						}
						else if (i == 1)
						{
							pre_boundary(ux, side, Sx, hx, PhysCoef);
						}
						else
						{
							general(ux, side, Sx, hx, PhysCoef, false);
						}
					}

					//Side::east
					if (i > 0)
					{
						Side side = Side::east;
						if (i == nx)
						{
							if (ux.boundary.type(side) == MathBoundary::Dirichlet)
							{
								fixed_node(ux, side);							
								M.add_line_with_map(a.get_map(), l);
								continue;
							}

							if (ux.boundary.type(side) == MathBoundary::Neumann)
							{
								deriv_at_boundary(ux, side, -1, Sx, hx, PhysCoef);
								//M.add_line_with_map(a.get_map(), l);
								//continue;
							}
						}
						else if (i == nx - 1)
						{
							pre_boundary(ux, side, Sx, hx, PhysCoef);
						}
						else
						{
							general(ux, side, Sx, hx, PhysCoef, false);
						}
					}

					if (dim > 1)
					{
						if (i == 0 || i == nx) S_reduced = 0.5;
						general(ux, Side::south, Sy * S_reduced, hy, PhysCoef, j == 0);
						general(ux, Side::north, Sy * S_reduced, hy, PhysCoef, j == ny - 1);

						if (j == 0)		 reduce(ux, Side::south, reduced);
						if (j == ny - 1) reduce(ux, Side::north, reduced);
					}
					if (dim > 2)
					{
						if (i == 0 || i == nx) S_reduced = 0.5;
						general(ux, Side::front, Sz * S_reduced, hz, PhysCoef, k == 0);
						general(ux, Side::back, Sz * S_reduced, hz, PhysCoef, k == nz - 1);

						if (k == 0)		 reduce(ux, Side::front, reduced);
						if (k == nz - 1) reduce(ux, Side::back, reduced);
					}
					if (i == 0 || i == nx) reduced *= 0.5;

					a(Side::center) += dV / tau * reduced;


					M.add_line_with_map(a.get_map(), l);
				}
			}
		}
	}

}

void FlowSolver::form_rhs_test(double* b, bool reset)
{

	if (ux.type == Component::x)
	{
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i <= nx; i++) {
					int l = i + ux.off * j + ux.off2 * k;
					if (reset) b[l] = 0;

					double SX = this->Sx;
					double SY = this->Sy;
					double SZ = this->Sz;
					double DV = this->dV;
					double dp = 0;
					double reduced = 1;
					double PhysCoef = 1.0 / Re;

					auto fixed_node = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
					{
						double res = 0.0;
						if (V.boundary.type(side) == MathBoundary::Dirichlet)
							res = V.boundary.get_value(side);
						else if (V.boundary.type(side) == MathBoundary::Neumann)
							res = V.boundary.normal_deriv_oriented(side);

						return res;
					};
					auto half_border = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
					{
						double res = 0.0;
						if (V.boundary.type(side) == MathBoundary::Neumann)
							res = V.boundary.normal_deriv_oriented(side) * S * coef;
						else if (V.boundary.type(side) == MathBoundary::Dirichlet)
						{
							res = V.boundary.get_value(side) * coef * S / (0.5 * h);
						}
						else
						{
							print("half_step error");
						}
						return res;
					};
					auto reduce = [](Velocity& V, Side side, double& reduced)
					{
						if (V.boundary.type(side) == MathBoundary::Dirichlet)
						{
							reduced *= 0.75;
						}
					};
					auto pre_fixed_node = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
					{
						double res = 0.0;
						if (V.boundary.type(side) == MathBoundary::Dirichlet)
							res = V.boundary.get_value(side) * coef * S / (h);
						//else if (V.boundary.type(side) == MathBoundary::Neumann)
						return res;
					};

					auto general = [this, b, i, j, k, l](Velocity& v, Velocity& v0, Side side, double dp, double dV, double tau, double S, double h, double coef, double reduced = 1)
					{
						return v0(i, j, k) * dV * reduced / tau
							- P.get_shifted_diff(Component::x, i, j, k) * S * reduced
							+ Ra * grav.x * T.get_shifted(Component::x, i, j, k) * dV * reduced;
					};


					// west
					if (i < nx)
					{
						Side side = Side::west;
						if (i == 0)
						{
							//fixed node
							if (ux.boundary.type(side) == MathBoundary::Dirichlet)
							{
								b[l] = ux.boundary.get_value(side);
								continue;
							}
							else if (ux.boundary.type(side) == MathBoundary::Neumann)
							{
								b[l] = PhysCoef * ux.boundary.normal_deriv_oriented(side) * Sx;
								//continue;
							}
						}
						else if (i == 1)
						{
							if (ux.boundary.type(side) == MathBoundary::Dirichlet)
								b[l] += half_border(ux, l, side, Sx, 2.0 * hx, PhysCoef);
						}
					}

					// east
					if (i > 0)
					{
						Side side = Side::east;
						if (i == nx)
						{
							//fixed node
							if (ux.boundary.type(side) == MathBoundary::Dirichlet)
							{
								b[l] = ux.boundary.get_value(side);
								continue;
							}
							else if (ux.boundary.type(side) == MathBoundary::Neumann)
							{
								b[l] = PhysCoef * ux.boundary.normal_deriv_oriented(side) * Sx;
								//continue;
							}
						}
						else if (i == nx - 1)
						{
							if (ux.boundary.type(side) == MathBoundary::Dirichlet)
								b[l] += half_border(ux, l, side, Sx, 2.0 * hx, PhysCoef);
						}
					}


					if (dim > 1)
					{
						//south
						if (j == 0)
						{
							Side side = Side::south;
							if (i == 0 || i == nx - 1)
								SY = SY * 0.5;
							b[l] += half_border(ux, l, side, SY, hy, PhysCoef);
							reduce(ux, side, reduced);
						}

						//north
						if (j == ny - 1)
						{
							Side side = Side::north;
							if (i == 0 || i == nx - 1)
								SY = SY * 0.5;
							b[l] += half_border(ux, l, Side::north, SY, hy, PhysCoef);
							reduce(ux, side, reduced);
						}
					}


					if (dim > 2)
					{
						//front
						if (k == 0)
						{
							Side side = Side::front;
							if (i == 0 || i == nx - 1)
								SZ = SZ * 0.5;
							b[l] += half_border(ux, l, side, SZ, hz, PhysCoef);
							reduce(ux, side, reduced);
						}
						//back
						if (k == nz - 1)
						{
							Side side = Side::back;
							if (i == 0 || i == nx - 1)
								SZ = SZ * 0.5;
							b[l] += half_border(ux, l, side, Sz, hz, PhysCoef);
							reduce(ux, side, reduced);
						}
					}

					//general
					if (i == 0)
					{
						dp = P(0, j, k) - P.boundary(Side::west, P(0, j, k));
						//dp = P.get_diff_x(Side::west, i, j, k);
						DV = DV * 0.5;
					}

					else if (i == nx)
					{
						dp = P.boundary(Side::east, P(nx - 1, j, k)) - P(nx - 1, j, k);
						//dp = P.get_diff_x(Side::east, i, j, k);
						DV = DV * 0.5;
					}
					else
						dp = P(i, j, k) - P(i - 1, j, k);


					//dp = P.get_shifted_diff(Component::x, i, j, k);

					b[l] += general(ux, vx, Side::center, dp, DV, tau, Sx, hx, PhysCoef, reduced);

				}
			}
		}
	}
}


#define DIRICHLET(side) boundary.type(side) == MathBoundary::Dirichlet
#define NEUMANN(side) boundary.type(side) == MathBoundary::Neumann

void FlowSolver::form_big_matrix(SparseMatrix& M, double* b)
{
	timer.start("matrix_formation");
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

	auto fixed_node = [b](Velocity& v, Side side, int l, Contribution &a)
	{
		a[l] = 1;
		b[l] = v.boundary.get_value(side);
	};

	auto deriv_at_boundary = [b, &M](Velocity& v, Side side, int shift, double S, double h, double coef, int l, Contribution &a)
	{
		if (v.boundary.type(side) == MathBoundary::Neumann)
		{
			a[l] = coef * S / h;
			a[l + shift] = -coef * S / h;;
			b[l] += v.boundary.normal_deriv_oriented(side);
		}
	};

	auto pre_boundary = [b, &M](Velocity& v, Side side, double S, double h, double PhysCoef, int l, Contribution &a)
	{
		if (v.boundary.type(side) == MathBoundary::Dirichlet)
		{
			a(Side::center) += PhysCoef * S / h;
			b[l] += v.boundary.get_value(side) * S / h * PhysCoef;
		}
	};

	if(1)
	{
		Velocity& u = ux;
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

						double PhysCoef = 1.0 / Re;

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

						double PhysCoef = 1.0 / Re;

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
	auto fixed_node = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
	{
		double res = 0.0;
		if (V.boundary.type(side) == MathBoundary::Dirichlet)
			res = V.boundary.get_value(side);
		else if (V.boundary.type(side) == MathBoundary::Neumann)
			res = V.boundary.normal_deriv_oriented(side);

		return res;
	};
	auto half_border = [](Velocity& V, int l, Side side, double S, double h, double coef, double f = 0)
	{
		double res = 0.0;
		if (V.boundary.type(side) == MathBoundary::Neumann)
			res = V.boundary.normal_deriv_oriented(side) * S * coef;
		else if (V.boundary.type(side) == MathBoundary::Dirichlet)
		{
			res = V.boundary.get_value(side) * coef * S / (0.5 * h);
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
			res = V.boundary.get_value(side) * coef * S / (h);
		//else if (V.boundary.type(side) == MathBoundary::Neumann)
		return res;
	};

	auto FluxVX = [this](int i, int j, int k, double Sx, double Sy, double Sz)
	{
		double FX = Sx * (vx.get_for_vx_cell(Side::east, i, j, k) *  vx.get_for_vx_cell(Side::east, i, j, k)
						- vx.get_for_vx_cell(Side::west, i, j, k) *  vx.get_for_vx_cell(Side::west, i, j, k));
		double FY = Sy * (vx.get_for_vx_cell(Side::north, i, j, k) * vy.get_for_vx_cell(Side::north, i, j, k)
						- vx.get_for_vx_cell(Side::south, i, j, k) * vy.get_for_vx_cell(Side::south, i, j, k));
		double FZ = Sz * (vx.get_for_vx_cell(Side::back, i, j, k) *  vz.get_for_vx_cell(Side::back, i, j, k)
						- vx.get_for_vx_cell(Side::front, i, j, k) * vz.get_for_vx_cell(Side::front, i, j, k));
		return FX + FY + FZ;
	};
	auto FluxVY = [this](int i, int j, int k, double Sx, double Sy, double Sz)
	{
		double FX = Sx * (vx.get_for_vy_cell(Side::east, i, j, k) *  vy.get_for_vy_cell(Side::east, i, j, k)
						- vx.get_for_vy_cell(Side::west, i, j, k) *  vy.get_for_vy_cell(Side::west, i, j, k));
		double FY = Sy * (vy.get_for_vy_cell(Side::north, i, j, k) * vy.get_for_vy_cell(Side::north, i, j, k)
						- vy.get_for_vy_cell(Side::south, i, j, k) * vy.get_for_vy_cell(Side::south, i, j, k));
		double FZ = Sz * (vy.get_for_vy_cell(Side::back, i, j, k) *  vz.get_for_vy_cell(Side::back, i, j, k)
						- vy.get_for_vy_cell(Side::front, i, j, k) * vz.get_for_vy_cell(Side::front, i, j, k));
		return FX + FY + FZ;
	};
	auto FluxVZ = [this](int i, int j, int k, double Sx, double Sy, double Sz)
	{
		double FX = Sx * (vx.get_for_vz_cell(Side::east, i, j, k) *  vz.get_for_vz_cell(Side::east, i, j, k)
						- vx.get_for_vz_cell(Side::west, i, j, k) *  vz.get_for_vz_cell(Side::west, i, j, k));
		double FY = Sy * (vy.get_for_vz_cell(Side::north, i, j, k) * vz.get_for_vz_cell(Side::north, i, j, k)
						- vy.get_for_vz_cell(Side::south, i, j, k) * vz.get_for_vz_cell(Side::south, i, j, k));
		double FZ = Sz * (vz.get_for_vz_cell(Side::back, i, j, k) *  vz.get_for_vz_cell(Side::back, i, j, k)
						- vz.get_for_vz_cell(Side::front, i, j, k) * vz.get_for_vz_cell(Side::front, i, j, k));
		return FX + FY + FZ;
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


						double PhysCoef = 1.0 / Re;

						auto general = [this, b, i, j, k, l](Velocity& v, Velocity& v0, Side side, double dp, double dV, double tau, double S, double h, double coef, double reduced = 1)
						{
							return v0(i, j, k) * dV * reduced / tau
								- P.get_shifted_diff(Component::x, i, j, k) * S * reduced
								+ Ra * grav.x * T.get_shifted(Component::x, i, j, k) * dV * reduced;
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
									b[l] = u.boundary.get_value(side);
									continue;
								}
								else if (u.boundary.type(side) == MathBoundary::Neumann)
								{
									b[l] = PhysCoef * u.boundary.normal_deriv_oriented(side) * SX;
									//continue;
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
									b[l] = u.boundary.get_value(side);
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

						b[l] += general(u, v, Side::center, 1, DV, tau, SX, hx, PhysCoef, 1.0);
						b[l] += -FluxVX(i, j, k, SX*rx, SY*ry, SZ*rz);
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
						double PhysCoef = 1.0 / Re;

						auto general = [this, b, i, j, k, l](Velocity& v, Velocity& v0, Side side, double dp, double dV, double tau, double S, double h, double coef, double reduced = 1)
						{
							//print(i << " " << j << " " << T.get_shifted(Component::y, i, j, k))
							return v0(i, j, k) * dV * reduced / tau
								- P.get_shifted_diff(Component::y, i, j, k) * S * reduced
								+ Ra * grav.y * T.get_shifted(Component::y, i, j, k) * dV * reduced;
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
									b[l] = u.boundary.get_value(side);
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
									b[l] = u.boundary.get_value(side);
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

						b[l] += general(u, v, Side::center, 1, DV, tau, SY, hy, PhysCoef, 1);
						b[l] += -FluxVY(i, j, k, SX*rx, SY*ry, SZ*rz);
					}
				}
			}
		}
	}
}
