#include "FlowSolver.h"


void FlowSolver::solve_heat_equation_explicitly(int steps_at_ones)
{
	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
	{
		double lapl = 0.0;
		lapl += Sx * (f.get_dx(Side::east, i, j, k) - f.get_dx(Side::west, i, j, k));
		lapl += Sy * (f.get_dy(Side::north, i, j, k) - f.get_dy(Side::south, i, j, k));
		lapl += Sz * (f.get_dz(Side::back, i, j, k) - f.get_dz(Side::front, i, j, k));
		return lapl;
	};

	auto vF = [this](ScalarVariable& f, Velocity& VX, Velocity& VY, Velocity& VZ, int i, int j, int k)
	{
		double vxF = Sx * (f.get_at_side(Side::east, i, j, k) * VX.get_for_centered_cell(Side::east, i, j, k)
			- f.get_at_side(Side::west, i, j, k) * VX.get_for_centered_cell(Side::west, i, j, k));
		double vyF = Sy * (f.get_at_side(Side::north, i, j, k) * VY.get_for_centered_cell(Side::north, i, j, k)
			- f.get_at_side(Side::south, i, j, k) * VY.get_for_centered_cell(Side::south, i, j, k));
		double vzF = Sz * (f.get_at_side(Side::back, i, j, k) * VZ.get_for_centered_cell(Side::back, i, j, k)
			- f.get_at_side(Side::front, i, j, k) * VZ.get_for_centered_cell(Side::front, i, j, k));

		return vxF + vyF + vzF;
	};

	for (int iter = 0; iter < steps_at_ones; iter++)
	{
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++)
				{
					T(i, j, k) = T0(i, j, k) + tau / dV *
						(laplace(T0, i, j, k) / Pr
							- vF(T0, vx, vy, vz, i, j, k));
				}
			}
		}
		ScalarVariable::swap(T, T0);
	}
}

void FlowSolver::form_matrix_for_heat_equation(SparseMatrix &M, double *b)
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
		double& operator[](Side side)
		{
			return m[side];
		}
	};

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				map<int, double> line;
				Contribution a;

				auto flux = [this, b, &a, i, j, k, l](Side side, double S, double h, double coef, bool border)
				{
					if (border)
					{
						if (T.boundary.type(side) == MathBoundary::Neumann)
						{
							b[l] += T.boundary.normal_deriv_oriented(side) * S * coef;
						}
						if (T.boundary.type(side) == MathBoundary::Dirichlet)
						{
							a[Side::center] += coef * S / (0.5 * h);
							b[l] += T.boundary(side) * coef * S / (0.5 * h);
						}
					}
					else
					{
						a[Side::center] +=  coef * S / h;
						a[side] += -coef * S / h;
					}
				};

				double PhysCoef = 1.0 / Pr;
				flux(Side::west, Sx, hx, PhysCoef, i == 0);
				flux(Side::east, Sx, hx, PhysCoef, i == nx - 1);
				

				if (dim > 1)
				{
					flux(Side::south, Sy, hy, PhysCoef, j == 0);
					flux(Side::north, Sy, hy, PhysCoef, j == ny - 1);
				}
				if (dim > 2)
				{
					flux(Side::front, Sz, hz, PhysCoef, k == 0);
					flux(Side::back, Sz, hz, PhysCoef, k == nz - 1);
				}

				a[Side::center] += dV / tau;

				line[l] = a[Side::center];
				line[l - 1] = a[Side::west];
				line[l + 1] = a[Side::east];

				if (dim > 1) {
					line[l - off] = a[Side::south];
					line[l + off] = a[Side::north];
				}
				if (dim > 2) {
					line[l - off2] = a[Side::front];
					line[l + off2] = a[Side::back];
				}
				M.add_line_with_map(line, l);
			}
		}
	}
	//nval = (int)val.size();

}

void FlowSolver::form_rhs_for_heat_equation(double* b, bool reset)
{
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				map<int, double> line;
				if (reset) b[l] = 0.0;
				b[l] = T0(i, j, k) * dV / tau;


				auto bc = [this, b, i, j, k, l](Side side, double S, double h, double coef, bool border)
				{
					if (border)
					{
						if (T.boundary.type(side) == MathBoundary::Neumann)
						{
							b[l] += T.boundary.normal_deriv_oriented(side) * S * coef;
						}
						if (T.boundary.type(side) == MathBoundary::Dirichlet)
						{
							b[l] += T.boundary(side) * coef * S / (0.5 * h);
						}
					}
				};

				double PhysCoef = 1.0 / Pr;
				bc(Side::west, Sx, hx, PhysCoef, i == 0);
				bc(Side::east, Sx, hx, PhysCoef, i == nx - 1);


				if (dim > 1)
				{
					bc(Side::south, Sy, hy, PhysCoef, j == 0);
					bc(Side::north, Sy, hy, PhysCoef, j == ny - 1);
				}
				if (dim > 2)
				{
					bc(Side::front, Sz, hz, PhysCoef, k == 0);
					bc(Side::back, Sz, hz, PhysCoef, k == nz - 1);
				}

			}
		}
	}
}


void FlowSolver::solve_heat_equation(int steps_at_ones)
{
	if (sm.nval == 0)
	{
		form_matrix_for_heat_equation(sm, b);
	}

	form_rhs_for_heat_equation(b);
	itsol.solveGS(T.get_ptr(), T0.get_ptr(), b, N, sm);

}