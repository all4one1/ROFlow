#include "FlowSolver.h"
#include <omp.h>


void FlowSolver::guessed_velocity_for_simple()
{
	timer.start("guessed_u");
	form_rhs_test(B, true);
	//form_big_rhs(B, true);
	//itsol.solveGS(U, U0, B, NV, SM);
	itsol.solveGS(U, U0, B, NV, SM);

	timer.end("guessed_u");
}
void FlowSolver::poisson_equation_for_p_prime()
{
	timer.start("p_prime");
	double tau_p = 20.25 * pow(fmin(fmin(hx, hy), hz), 2);
	int iter_p = 0;
	double res, res0 = 0;
	double abs_eps = 1, rel_eps = 1, eps0 = 1e-4;

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
		//if (iter_p > 4000) break;
		//if (iter_p > iter_limit) break;
		res = 0;
		for (int k = 0; k < nz; k++) {
			#pragma omp parallel for num_threads(4) reduction(+:res)
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

		abs_eps = abs(res - res0);
		rel_eps = abs_eps / res0;
		res0 = res;
		
		
		if (abs_eps < eps0*10 || rel_eps < eps0) break;


		//double max = 0;
		//for (int l = 0; l < N; l++)
		//{
		//	double dif = abs(p_prime[l] - buffer[l]);
		//	if (dif > max)	max = dif;
		//}
		//if (max < eps0)	break;
		//print(iter << " " << max)
		//print(iter_p)

		p_prime.transfer_data_to(buffer);

		//if (iter_p % 10000 == 0)	cout << "iter_p: " << iter_p << " res: " << res << " eps: " << eps << endl;
		//if (iter_p % 10000 == 0)	cout << "eps: " << abs(res - res0) << " " << eps << endl;
	}
	timer.end("p_prime");

	if (iter == 1) ofstream f("p_prime.dat");
	if (iter % 1 == 0)
	{
		ofstream f("p_prime.dat", ios_base::app);
		f << iter << " " << iter_p << " " << abs_eps << " " << rel_eps << " " << res << endl;
	}
}
void FlowSolver::correction_for_simple()
{
	timer.start("correction");
	double alpha = 1;

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


void FlowSolver::solve_simple(size_t steps_at_ones)
{
	//form_matrix_test(SM, B);
	form_big_matrix(SM, B);

	for (size_t q = 0; q < steps_at_ones; q++)
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
			double dp = P.boundary.get_fixed_value(Side::east) - P.boundary.get_fixed_value(Side::west);
			double y = vx.y_(j);
			return Re * 0.5 * (dp) / Lx * (y * y - y * Ly);
		};

	for (int j = 0; j < ny; j++)
	{
		cout << fixed << vx.y_(j) << " " << vx(0, j, 0) << " " << vx(nx / 2, j, 0) << " " << FF(j) << endl;
	}

}

void FlowSolver::solve_system(size_t steps_at_ones)
{
	if (iter == 0)
	{
		//form_big_matrix(SM, B);
		form_matrix_test(SM, B);
		form_matrix_for_heat_equation(SMt, b);
	}
	
	for (size_t q = 0; q < steps_at_ones; q++)
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

			double div = check_div2();
			if (sq % 100 == 0) print("div = " << div << " " << sq)
				if (div < 1e-5) break;
			if (sq > 10000) { print("bad div"); break; }
		}

		ux.transfer_data_to(vx);
		uy.transfer_data_to(vy);

		timer.start("T,C equations");
		form_rhs_for_heat_equation(b, true);
		itsol.solveGS(T.get_ptr(), T0.get_ptr(), b, N, SMt);


		//form_rhs_for_concentration_equation(b, true);
		//itsol.solveGS(C.get_ptr(), C0.get_ptr(), b, N, SMc);
		timer.end("T,C equations");


		if (Rav != 0)
			poisson_equation_pulsation_stream_function();


		//if (q % 100 == 0) print(uy(nx / 4, ny / 2, 0) << " " << check_div());
		
		if ((q + 1) % 100 == 0)
		{
			double Ek, Vmax;
			statistics(Ek, Vmax);

			print(iter << ", t = " << total_time << ", Ek = " << Ek << ", Vmax = " << Vmax);
			
			temporal.write_header("Ra, Rav, total_time, Ek, Vmax, check_Ek.dif, check_Ek.long_dif");
			temporal.write({ Ra, Rav, total_time, Ek, Vmax, check_Ek.dif, check_Ek.long_dif });


			if (check_Ek.stop(Ek, false))
			{
				break;
			}
		}
		if ((q + 1) % 1000 == 0)
		{
			write_fields("F.dat");
		}

	}
}