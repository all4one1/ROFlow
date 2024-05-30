#include "FlowSolver.h"


void FlowSolver::form_u_matrix_for_simple(SparseMatrix &M, double *b)
{
	struct Contribution
	{
		std::map<Side, double> m;
		Contribution()
		{
			m[Side::center] = 0.0;
			m[Side::west] = 0.0;
			m[Side::east] = 0.0;
			m[Side::south] = 0.0;
			m[Side::north] = 0.0;
			m[Side::front] = 0.0;
			m[Side::back] = 0.0;
		}

		void send_to_sparse(map<int, double>& line, int l, int dim, int off, int off2)
		{
			line[l] = (*this)[Side::center];
			line[l - 1] = (*this)[Side::west];
			line[l + 1] = (*this)[Side::east];

			if (dim > 1) {
				line[l - off] = (*this)[Side::south];
				line[l + off] = (*this)[Side::north];
			}
			if (dim > 2) {
				line[l - off2] = (*this)[Side::front];
				line[l + off2] = (*this)[Side::back];
			}
		}
		double& operator[](Side side)
		{
			return m[side];
		}
	};
	struct Contribution2
	{
		std::map<int, double> m;
		int l, off, off2;
		Contribution2(int l_, int off_ = 0, int off2_ = 0) : l(l_), off(off_), off2(off2_)
		{

		}
		int index (Side side)
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
		map<int, double> &get_map()
		{
			return m;
		}
	};
	auto rhs = [](bool fixed, Velocity &V, int l, Side side, double S, double h, double coef, double f = 0)
	{
		if (fixed)
		{
			if (V.boundary.type(side) == MathBoundary::Neumann)
				return V.boundary.normal_deriv_oriented(side); //?
			if (V.boundary.type(side) == MathBoundary::Dirichlet)
				return V.boundary.get_value(side);
		}
		else 
		{
			if (V.boundary.type(side) == MathBoundary::Neumann)
				return V.boundary.normal_deriv_oriented(side) * S * coef; //?
			if (V.boundary.type(side) == MathBoundary::Dirichlet)
				return V.boundary.get_value(side) * coef * S / (0.5 * h);
		}

		MYERROR("wrong bc");
		return 0.0;
	};
	auto index = [](Side side, int l, int off, int off2)
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
			break;
		}
	};

	if (ux.type == Component::x)
	{
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i <= nx; i++) {
					int l = i + ux.off * j + ux.off2 * k;
					//map<int, double> line;
					Contribution2 a(l, ux.off, ux.off2);
					double reduced = 1.0;
					auto reduce = [](Velocity& V, Side side, double& reduced)
					{
						if (V.boundary.type(side) == MathBoundary::Dirichlet)
						{
							reduced *= 0.75;
						}
					};
					double PhysCoef = 1.0 / Re;

					auto general = [&rhs, b, &a, i, j, k, l](Velocity &v, Side side, double S, double h, double coef, bool border = false)
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
							}
							if (ux.boundary.type(side) == MathBoundary::Neumann)
							{
								deriv_at_boundary(ux, side, +1, Sx, hx, PhysCoef);

								if (dim > 1)
								{
									general(ux, Side::south, 0.5 * Sy, hy, PhysCoef, j == 0);
									general(ux, Side::north, 0.5 * Sy, hy, PhysCoef, j == ny - 1);
								}
								if (dim > 2)
								{
									general(ux, Side::front, 0.5 * Sz, hz, PhysCoef, k == 0);
									general(ux, Side::back, 0.5 * Sz, hz, PhysCoef, k == nz - 1);
								}
							}
							M.add_line_with_map(a.get_map(), l);
							continue;
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
							}

							if (ux.boundary.type(side) == MathBoundary::Neumann)
							{
								deriv_at_boundary(ux, side, -1, Sx, hx, PhysCoef);

								if (dim > 1)
								{
									general(ux, Side::south, 0.5 * Sy, hy, PhysCoef, j == 0);
									general(ux, Side::north, 0.5 * Sy, hy, PhysCoef, j == ny - 1);
								}
								if (dim > 2)
								{
									general(ux, Side::front, 0.5 * Sz, hz, PhysCoef, k == 0);
									general(ux, Side::back, 0.5 * Sz, hz, PhysCoef, k == nz - 1);
								}
							}
							M.add_line_with_map(a.get_map(), l);
							continue;
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
						general(ux, Side::south, Sy, hy, PhysCoef, j == 0);
						general(ux, Side::north, Sy, hy, PhysCoef, j == ny - 1);

						if (j == 0)		 reduce(ux, Side::south, reduced);
						if (j == ny - 1) reduce(ux, Side::north, reduced);
					}
					if (dim > 2)
					{
						general(ux, Side::front, Sz, hz, PhysCoef, k == 0);
						general(ux, Side::back, Sz, hz, PhysCoef, k == nz - 1);

						if (k == 0)		 reduce(ux, Side::front, reduced);
						if (k == nz - 1) reduce(ux, Side::back, reduced);
					}

					a(Side::center) += dV / tau * reduced;


					M.add_line_with_map(a.get_map(), l);
				}
			}
		}
	}
}

void FlowSolver::form_u_rhs_for_simple(double* b, bool reset)
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

					auto general = [b, i, j, k, l](Velocity& v, Velocity& v0, Side side, double dp, double dV, double tau, double S, double h, double coef, double reduced = 1)
					{
						return v0(i, j, k) * dV * reduced / tau - dp * S * reduced;
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
							b[l] += half_border(ux, l, side, SZ, hz, PhysCoef);
							reduce(ux, side, reduced);
						}
					}

					//general
					if (i == 0)
					{
						
						dp = P(0, j, k) - P.boundary(Side::west, P(0, j, k)); 
						DV = DV * 0.5;
					}

					else if (i == nx)
					{
						dp = P.boundary(Side::east, P(nx - 1, j, k)) - P(nx - 1, j, k);
						DV = DV * 0.5;
					}
					else
						dp = P(i, j, k) - P(i - 1, j, k);


					b[l] += general(ux, vx, Side::center, dp, DV, tau, Sx, hx, PhysCoef, reduced);

				}
			}
		}
	}

}


void FlowSolver::guessed_velocity_for_simple()
{
	timer.start("guessed_u");
	
	form_big_rhs(B, true);
	itsol.solveGS(U, U0, B, NV, SM);

	timer.end("guessed_u");
}

void FlowSolver::poisson_equation_for_p_prime()
{
	timer.start("p_prime");
	double tau_p = 30.25 * pow(fmin(fmin(hx, hy), hz), 2);
	int iter_p = 0;
	double res, res0 = 0;
	double eps = 1, eps0 = 1e-5;

	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
	{
		double lapl = 0.0;
		double aw = SM.get_diag(ux.get_l(i, j, k));
		double ae = SM.get_diag(ux.get_l(i + 1, j, k));
		double as = SM.get_diag(uy.get_l(i, j, k) + stride);
		double an = SM.get_diag(uy.get_l(i, j + 1, k) + stride);

		lapl += Sx * (f.get_diff_x(Side::east, i, j, k) / ae - f.get_diff_x(Side::west, i, j, k) / aw);
		lapl += Sy * (f.get_diff_y(Side::north, i, j, k) / an - f.get_diff_y(Side::south, i, j, k) / as);

		//lapl += Sx * (f.get_dx(Side::east, i, j, k) / ae - f.get_dx(Side::west, i, j, k) / aw);
		//lapl += Sy * (f.get_dy(Side::north, i, j, k) - f.get_dy(Side::south, i, j, k));
		//lapl += Sz * (f.get_dz(Side::back, i, j, k) - f.get_dz(Side::front, i, j, k));
		return lapl;
	};
	auto divU = [this](int i, int j = 0, int k = 0)
	{
		double div = 0.0;
		div += Sx * (ux.get_for_centered_cell(Side::east, i, j, k) - ux.get_for_centered_cell(Side::west, i, j, k));
		div += Sy * (uy.get_for_centered_cell(Side::north, i, j, k) - uy.get_for_centered_cell(Side::south, i, j, k));
		div += Sz * (uz.get_for_centered_cell(Side::back, i, j, k) - uz.get_for_centered_cell(Side::front, i, j, k));
		return div;
	};

	//double aw, ae, as, an, af, ab;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k; 
				double div = 0.0;
				div += Sx * (ux(i + 1, j, k) - ux(i, j, k));
				div += Sy * (uy(i, j + 1, k) - uy(i, j, k));

				bufferU(i, j, k) = div;
			}
		}
	}



	buffer.reset_to_null();
	buffer.boundary = p_prime.boundary;



	//while (eps > eps0 * res0)
	while(true)
	{
		iter_p = iter_p + 1;
		//if (iter_p > iter_limit) break;
		res = 0;
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + off * j + off2 * k;

					p_prime[l] = buffer[l] + tau_p / dV * (
						laplace(buffer, i, j, k) 
						- bufferU(i,j,k)
						);

					res = res + abs(p_prime[l]);
				}
			}
		}
		eps = abs(res - res0);
		res0 = res;


		double max = 0;
		double dif;
		for (int l = 0; l < N; l++)
		{
			dif = abs(p_prime[l] - buffer[l]);
			if (dif > max)
				max = dif;
		}
		if (max < eps0)	break;

		if (eps < 1e-12) break;

		p_prime.transfer_data_to(buffer);

		//if (iter_p % 10000 == 0)	cout << "iter_p: " << iter_p << " res: " << res << " eps: " << eps << endl;
		//if (iter_p % 10000 == 0)	cout << "eps: " << abs(res - res0) << " " << eps << endl;
		
		//if (k > 1000) break;

		
	}
	//iter_limit = iter_p;
	
	timer.end("p_prime");

	if (iter == 1) ofstream f("p_prime.dat");
	ofstream f("p_prime.dat", ios_base::app);
	f << iter << " " << iter_p << " " << eps << " " << res << endl;

}

//void FlowSolver::poisson_equation_for_p_prime2()
//{
//	double tau_p = 10.25 * pow(fmin(fmin(hx, hy), hz), 2);
//	int iter_p = 0;
//	double res, res0 = 0;
//	double eps = 1, eps0 = 1e-6;
//
//	buffer.reset_to_null();
//	buffer.boundary = p_prime.boundary;
//
//
//	auto div = [this](int i, int j, int k)
//	{
//		return Sx * (ux(i + 1, j, k) - ux(i, j, k));
//	};
//
//	while (eps > eps0 * res0 && eps > 1e-12)
//	{
//		if (iter_p > iter_limit) break;
//		res = 0;
//		for (int k = 0; k < nz; k++) {
//			for (int j = 0; j < ny; j++) {
//				for (int i = 0; i < nx; i++) {
//					int l = i + off * j + off2 * k;
//
//
//					//if (i == 0 || i == nx - 1)	c = 0.5;
//					//double PE = ddx(i + 1, j, k) * (buffer(i + 1, j, k) - buffer(i, j, k));
//					//double PW = ddx(i, j, k) * (buffer(i, j, k) - buffer(i - 1, j, k));
//					double PE = ddx(i + 1, j, k) * (buffer.get_diff_x(Side::east, i, j, k));
//					double PW = ddx(i, j, k) * (buffer.get_diff_x(Side::west, i, j, k));
//
//					
//					p_prime(i, j, k) = buffer(i, j, k) + tau_p / dV * ((PE - PW)  - div(i, j, k));
//					
//					
//					res = res + abs(p_prime[l]);
//				}
//			}
//		}
//		eps = abs(res - res0);
//		res0 = res;
//
//		p_prime.transfer_data_to(buffer);
//		iter_p = iter_p + 1;
//		//if (iter_p % 10000 == 0)	cout << "iter_p: " << iter_p << " res: " << res << " eps: " << eps << endl;
//		if (iter_p % 100000 == 0)	cout << "eps: " << abs(res - res0) << " " << eps << endl;
//
//		//if (k > 1000) break;
//	}
//	iter_limit = iter_p;
//
//}
//void FlowSolver::correction2()
//{
//	double alpha = 0.1;
//
//	for (int k = 0; k < nz; k++) {
//		for (int j = 0; j < ny; j++) {
//			for (int i = 0; i < nx; i++) {
//				P(i, j, k) = P(i, j, k) + alpha * p_prime(i, j, k);
//			}
//		}
//	}
//
//	for (int k = 0; k < nz; k++) {
//		for (int j = 0; j < ny; j++) {
//			for (int i = 0; i <= nx; i++) {
//				auto FF = [this](double y)
//				{
//					return -Re / Lx / 2.0 * (y * y - y * Ly);
//				};
//				double vx_prime = 0;
//				double dp = 0.0;
//				int l = i + ux.off * j + ux.off2 * k;
//				if (i > 0 && i < nx)
//				{
//					dp = p_prime(i, j, k) - p_prime(i - 1, j, k);
//				}
//				else if (i == 0) dp = (p_prime(i, j, k) - 0);
//				else if (i == nx) dp = (0 - p_prime(i-1, j, k));
//
//				ux(i, j, k) = ux(i, j, k) - dp * ddx(i, j, k);
//			}
//		}
//	}
//}
//


void FlowSolver::correction_for_simple()
{
	timer.start("correction");
	double alpha = 0.8;

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				P(i, j, k) = P(i, j, k) + alpha * p_prime(i, j, k);
			}
		}
	}
	//auto FF = [this](double y)
	//{
	//	return -Re / Lx / 2.0 * (y * y - y * Ly);
	//};

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i <= nx; i++) {
				int l = i + ux.off * j + ux.off2 * k;

				double v_prime = -(p_prime.get_shifted_diff(Component::x, i, j, k)) * Sx / SM(l, l);
				ux(i, j, k) = ux(i, j, k) + v_prime;
			}
		}
	}

	if (1)
	{
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j <= ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + uy.off * j + uy.off2 * k + stride;
					double v_prime = -(p_prime.get_shifted_diff(Component::y, i, j, k)) * Sy / SM(l, l);
					uy(i, j, k) = uy(i, j, k) + v_prime;
				}
			}
		}
	}

	timer.end("correction");
}


void FlowSolver::solve_simple(int steps_at_ones)
{
	//form_matrix_test(SM, B);
	form_big_matrix(SM, B);

	for (int q = 0; q < steps_at_ones; q++)
	{
		iter++;
		total_time += tau;

		//form_rhs_test(B, true);
		form_big_rhs(B, true);

		int sq = 0;
		while (true)
		{
			sq++;
			form_big_rhs(B, true);
			itsol.solveGS(U, U0, B, Nvx + Nvy, SM);
			poisson_equation_for_p_prime();
			correction_for_simple();

			double div = check_div();
			if (sq % 100 == 0) print(div)
				if (div < 1e-4) break;
		}

		ux.transfer_data_to(vx);
		uy.transfer_data_to(vy);
		if (q % 100 == 0) print(ux(nx / 2, ny / 2, 0) << " " << check_div());
	}


	auto FF = [this](int j)
		{
			double dp = P.boundary.get_value(Side::east) - P.boundary.get_value(Side::west);
			double y = vx.y_(j);
			return Re * 0.5 * (dp) / Lx * (y * y - y * Ly);
		};

	for (int j = 0; j < ny; j++)
	{
		cout << fixed << vx.y_(j) << " " << vx(0, j, 0) << " " << vx(nx / 2, j, 0) << " " << FF(j) << endl;
	}

}

void FlowSolver::solve_system(int steps_at_ones)
{

	if (iter == 0)
	{
		form_big_matrix(SM, B);
		form_matrix_for_heat_equation(sm, b);
	}
	
	for (int q = 0; q < steps_at_ones; q++)
	{
		iter++;
		total_time += tau;
		
		int sq = 0;
		while (true)
		{
			sq++;
			guessed_velocity_for_simple();
			poisson_equation_for_p_prime();
			correction_for_simple();

			double div = check_div();
			if (sq % 100 == 0) print("div = " << div)
				if (div < 1e-7) break;
		}

		ux.transfer_data_to(vx);
		uy.transfer_data_to(vy);

		form_rhs_for_heat_equation(b, true);
		itsol.solveGS(T.get_ptr(), T0.get_ptr(), b, N, sm);
		//solve_heat_equation_explicitly();


		//if (q % 100 == 0) print(uy(nx / 4, ny / 2, 0) << " " << check_div());
		
		if ((q + 1) % 100 == 0)
		{
			double Ek, Vmax;
			statistics(Ek, Vmax);

			print("t = " << total_time << ", Ek = " << Ek << ", Vmax = " << Vmax);


			//for (auto &it : m_timer)	cout << it.first << ": " << it.second << endl;
			if (kinetic_check.stop(Ek, true))
			{
				stats.write({total_time, Ra, Ek, Vmax });
				break;
			}
		}
	}
}