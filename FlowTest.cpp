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



void FlowSolver::poisson_equation_pulsation_stream_function()
{
	timer.start("SF poisson");
	double tau_p = 0.25 * pow(fmin(fmin(hx, hy), hz), 2);
	int iter_p = 0;
	double res, res0 = 0;
	double abs_eps = 1, rel_eps = 1, eps0 = 1e-6;

	
	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
		{
			double lapl = 0.0;

			lapl += Sx * (f.get_diff_x(Side::east, i, j, k)  - f.get_diff_x(Side::west, i, j, k));
			lapl += Sy * (f.get_diff_y(Side::north, i, j, k)  - f.get_diff_y(Side::south, i, j, k));

			return lapl;
		};

	//double aw, ae, as, an, af, ab;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				double div = 0.0;
				div +=  (T.get_at_side(Side::east, i, j, k) -  T.get_at_side(Side::west, i, j, k)) * vibr.y * Sx;
				div += -(T.get_at_side(Side::north, i, j, k) - T.get_at_side(Side::south, i, j, k)) * vibr.x * Sy;

				div +=  K * (C.get_at_side(Side::east, i, j, k) -  C.get_at_side(Side::west, i, j, k)) * vibr.y * Sx;
				div += -K * (C.get_at_side(Side::north, i, j, k) - C.get_at_side(Side::south, i, j, k)) * vibr.x * Sy;

				bufferU(i, j, k) = div;
			}
		}
	}

	//SF0.reset_to_null();
	//buffer.boundary = F.boundary;



	//while (eps > eps0 * res0)
	while (true)
	{
		iter_p = iter_p + 1;
		//if (iter_p > iter_limit) break;
		res = 0;
		for (int k = 0; k < nz; k++) {
			#pragma omp parallel for num_threads(4) reduction(+:res)
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + off * j + off2 * k;

					SF[l] = SF0[l] + tau_p / dV * (
						laplace(SF0, i, j, k) - bufferU(i, j, k)
						);

					res = res + abs(SF[l]);
				}
			}
		}

		abs_eps = abs(res - res0);
		rel_eps = abs_eps / res0;
		res0 = res;

		if (abs_eps < eps0 * 10 || rel_eps < eps0) break;


		SF.transfer_data_to(SF0);
	}

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				double VX = +SF.get_dy(Side::center, i, j, k);
				double VY = -SF.get_dx(Side::center, i, j, k);
				bufferVibr[l] = VX * vibr.x + VY * vibr.y;
			}
		}
	}








	timer.end("SF poisson");
	//cout << "SF poisson: " << iter_p << endl;
}