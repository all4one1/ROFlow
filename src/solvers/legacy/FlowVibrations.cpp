#pragma once
#include "FlowSolver.h"
#define OMP 4

void FlowSolver::poisson_equation_pulsation_stream_function()
{
	timer.start("SF poisson");
	double tau_p = pow(fmin(fmin(hx, hy), hz), 2) / (dim * 2.0) * 0.99;
	int iter_p = 0;
	double res, res0 = 0;
	double abs_eps = 1, rel_eps = 1, eps0 = 1e-5;


	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
		{
			double lapl = 0.0;

			lapl += Sx * (f.get_dx(Side::east, i, j, k) - f.get_dx(Side::west, i, j, k));
			lapl += Sy * (f.get_dy(Side::north, i, j, k) - f.get_dy(Side::south, i, j, k));
			//if (dim > 1)
			//if (dim > 2) lapl += Sz * (f.get_dz(Side::back, i, j, k) - f.get_dz(Side::front, i, j, k));
			return lapl;
		};

	//double aw, ae, as, an, af, ab;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				double div = 0.0;
				div += (T.get_at_side(Side::east, i, j, k) - T.get_at_side(Side::west, i, j, k)) * vibr.y * Sx;
				div += -(T.get_at_side(Side::north, i, j, k) - T.get_at_side(Side::south, i, j, k)) * vibr.x * Sy;

				div += K * (C.get_at_side(Side::east, i, j, k) - C.get_at_side(Side::west, i, j, k)) * vibr.y * Sx;
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
		#pragma omp parallel for num_threads(OMP) reduction(+:res)
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + off * j + off2 * k;

					SF[l] = SF0[l] + tau_p / dV * (
						laplace(SF0, i, j, k) + bufferU(i, j, k)
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

void FlowSolver::poisson_equation_pulsation_velocity_x()
{
	timer.start("vibr poisson");
	double tau_p = 0.25 * pow(fmin(fmin(hx, hy), hz), 2);
	int iter_p = 0;
	double res, res0 = 0;
	double abs_eps = 1, rel_eps = 1, eps0 = 1e-6;


	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
		{
			double lapl = 0.0;
			//lapl += Sx * (f.get_diff_x(Side::east, i, j, k) - f.get_diff_x(Side::west, i, j, k));
			//lapl += Sy * (f.get_diff_y(Side::north, i, j, k) - f.get_diff_y(Side::south, i, j, k));
			lapl += Sx * (f.get_dx(Side::east, i, j, k) - f.get_dx(Side::west, i, j, k));
			lapl += Sy * (f.get_dy(Side::north, i, j, k) - f.get_dy(Side::south, i, j, k));
			return lapl;
		};

	//double aw, ae, as, an, af, ab;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				double div = 0.0;
				div += laplace(T, i, j, k) * vibr.x;

				div += -vibr.x * Sx * (T.get_dx(Side::east, i, j, k) - T.get_dx(Side::west, i, j, k));
				div += -vibr.y * Sy * (T.get_dx(Side::north, i, j, k) - T.get_dx(Side::south, i, j, k));

				bufferU(i, j, k) = div;
			}
		}
	}

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

					VibrX[l] = VibrX0[l] + tau_p / dV * (
						laplace(VibrX0, i, j, k) - bufferU(i, j, k)
						);

					if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
						res = res + abs(VibrX[l]);
				}
			}
		}

		abs_eps = abs(res - res0);
		rel_eps = abs_eps / res0;
		res0 = res;

		if (abs_eps < eps0 * 10 || rel_eps < eps0) break;


		VibrX.transfer_data_to(VibrX0);
	}


	timer.end("vibr poisson");
	//cout << "SF poisson: " << iter_p << endl;
}
void FlowSolver::poisson_equation_pulsation_velocity_y()
{
	timer.start("vibr poisson");
	double tau_p = 0.25 * pow(fmin(fmin(hx, hy), hz), 2);
	int iter_p = 0;
	double res, res0 = 0;
	double abs_eps = 1, rel_eps = 1, eps0 = 1e-6;


	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
		{
			double lapl = 0.0;
			//lapl += Sx * (f.get_diff_x(Side::east, i, j, k) - f.get_diff_x(Side::west, i, j, k));
			//lapl += Sy * (f.get_diff_y(Side::north, i, j, k) - f.get_diff_y(Side::south, i, j, k));
			lapl += Sx * (f.get_dx(Side::east, i, j, k) - f.get_dx(Side::west, i, j, k));
			lapl += Sy * (f.get_dy(Side::north, i, j, k) - f.get_dy(Side::south, i, j, k));
			return lapl;
		};

	//double aw, ae, as, an, af, ab;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;
				double div = 0.0;
				div += laplace(T, i, j, k) * vibr.y;

				div += -vibr.x * Sx * (T.get_dy(Side::east, i, j, k) - T.get_dy(Side::west, i, j, k));
				div += -vibr.y * Sy * (T.get_dy(Side::north, i, j, k) - T.get_dy(Side::south, i, j, k));

				bufferU(i, j, k) = div;
			}
		}
	}

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

					VibrY[l] = VibrY0[l] + tau_p / dV * (
						laplace(VibrY0, i, j, k) - bufferU(i, j, k)
						);

					if (i > 0 && i < nx - 1 && j > 0 && j < ny - 1)
						res = res + abs(VibrY[l]);
				}
			}
		}

		abs_eps = abs(res - res0);
		rel_eps = abs_eps / res0;
		res0 = res;

		if (abs_eps < eps0 * 10 || rel_eps < eps0) break;


		VibrY.transfer_data_to(VibrY0);
	}

	
	timer.end("vibr poisson");
	//cout << "SF poisson: " << iter_p << endl;
}
void FlowSolver::make_vibr_buffer()
{

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k;

				bufferVibr[l] = VibrX(i,j,k) * vibr.x + VibrY(i, j, k) * vibr.y;
			}
		}
	}

}