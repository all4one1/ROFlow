#include "FlowSolver.h"
#include <omp.h>


void FlowSolver::guessed_velocity_for_simple()
{
	timer.start("guessed_u");
	form_rhs_test(B, true);
	//form_rhs_Uxyz(B, true);
	itsol.solveGS(U, U0, B, NV, SM);

	timer.end("guessed_u");
}
void FlowSolver::poisson_equation_for_p_prime()
{
	timer.start("p_prime");
	double hmin = hx;
	if (dim > 1) hmin = fmin(hmin, hy);
	if (dim > 2) hmin = fmin(hmin, hz);

	double tau_p = 20.25 * pow(hmin, 2);
	iter_p = 0;
	double res = 0, res0 = 0;
	double abs_eps = 1, rel_eps = 1, eps0 = 1e-5;

	

	double a[6] = { 0,0,0,0,0,0 };

	double max_ = 0;
	double min_ = 10000;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				a[0] = SM.get_diag(ux.get_l(i, j, k));
				a[1] = SM.get_diag(ux.get_l(i + 1, j, k));
				if (dim > 1) {
					a[2] = SM.get_diag(uy.get_l(i, j, k));
					a[3] = SM.get_diag(uy.get_l(i, j + 1, k));
				}
				if (dim > 2) {
					a[4] = SM.get_diag(uz.get_l(i, j, k));
					a[5] = SM.get_diag(uz.get_l(i, j, k + 1));
				}
				for (int i = 0; i < 6; i++)
				{
					if (a[i] > max_) max_ = a[i];
					if (a[i] < min_ && a[i] > 0) min_ = a[i];
				}
			}
		}
	}
	//print(min_ * hmin / (2 * dim));
	tau_p = min_ * hmin / (2 * dim) * 0.99;

	
	auto laplace = [&a, this](ScalarVariable& f, int i, int j = 0, int k = 0)
	{
		double lapl = 0.0;
		a[0] = SM.get_diag(ux.get_l(i, j, k));
		a[1] = SM.get_diag(ux.get_l(i + 1, j, k));
		lapl += Sx * (f.get_diff_x(Side::east, i, j, k) / a[1] - f.get_diff_x(Side::west, i, j, k) / a[0]);

		if (dim > 1) {
			a[2] = SM.get_diag(uy.get_l(i, j, k));
			a[3] = SM.get_diag(uy.get_l(i, j + 1, k));
			lapl += Sy * (f.get_diff_y(Side::north, i, j, k) / a[3] - f.get_diff_y(Side::south, i, j, k) / a[2]);
		}
		if (dim > 2) {
			a[4] = SM.get_diag(uz.get_l(i, j, k));
			a[5] = SM.get_diag(uz.get_l(i, j, k + 1));
			lapl += Sz * (f.get_diff_z(Side::back, i, j, k) / a[5] - f.get_diff_z(Side::front, i, j, k) / a[4]);
		}
		return lapl;
	};


	//double aw, ae, as, an, af, ab;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				int l = i + off * j + off2 * k; 
				double div = 0.0;
				div += Sx * (ux(i + 1, j, k) - ux(i, j, k));
				if (dim > 1) div += Sy * (uy(i, j + 1, k) - uy(i, j, k));
				if (dim > 2) div += Sz * (uz(i, j, k + 1) - uz(i, j, k));

				bufferU(i, j, k) = div;
			}
		}
	}

	buffer.reset_to_null();
	buffer.boundary = p_prime.boundary;

	

	//while (eps > eps0 * res0)
	while(true)
	{
		iter_p++;
		//if (iter_p > 4000) break;
		//if (iter_p > iter_limit) break;
		res = 0;
		#pragma omp parallel for num_threads(4) reduction(+:res)
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

		abs_eps = abs(res - res0);
		rel_eps = abs_eps / res0;
		res0 = res;
		
		
		if (abs_eps < eps0*10 || rel_eps < eps0) break;

		p_prime.transfer_data_to(buffer);

		if (iter_p % 10000 == 0)	cout << "iter_p: " << iter_p << " res: " << res << " eps: " << abs_eps << endl;
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
	double alpha = alpha_relax;

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

	if (1 && dim > 1)
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

	if (1 && dim > 2)
	{
		for (int k = 0; k <= nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int l = i + uz.off * j + uz.off2 * k + stride2;
					double v_prime = -(p_prime.get_shifted_diff(Component::z, i, j, k)) * Sz / SM(l, l);
					uz(i, j, k) = uz(i, j, k) + v_prime;
				}
			}
		}
	}


	timer.end("correction");
}


void FlowSolver::solve_system(size_t steps_at_ones, bool inTimeUnits)
{
	if (iter == 0)
	{
		form_matrix_test(SM, B);
		form_matrix_for_heat_equation(SMt, b);
	}
	if (inTimeUnits)
		steps_at_ones = size_t((1.0 / tau) * steps_at_ones);

	auto tt = [this](double each)
		{
			size_t iterPerTimeUnit = tau > 1.0 ? 1 : size_t(1.0 / tau);
			int d = int(each * iterPerTimeUnit);
			if (d == 0)	return 1;
			else return d;
		};
	#define EVERY(each) (iter) % tt(each)== 0




	while (true)
	{
		iter++;
		total_time += tau;
		if (iter >= steps_at_ones) stop_signal = 1;

		timer.start("T,C equations");
		form_rhs_for_heat_equation(b, true);
		itsol.solveGS(T.get_ptr(), T0.get_ptr(), b, N, SMt);
		//form_rhs_for_concentration_equation(b, true);
		//itsol.solveGS(C.get_ptr(), C0.get_ptr(), b, N, SMc);
		timer.end("T,C equations");


		//alpha_relax = 0.75;
		iter_div = 0;
		while (1)
		{
			iter_div++;
			guessed_velocity_for_simple();
			poisson_equation_for_p_prime();
			correction_for_simple_test();

			double div = check_div2();
			if (iter_div % 100 == 0) print("div = " << div << " " << iter_div)
				if (div < 1e-5) break;
			if (iter_div > 20000) { print("bad div"); break; }
		}
		StateOut::simple_write(iter, "iterations.dat", { double(iter), double(iter_div), double(iter_p)});

		ux.transfer_data_to(vx);
		uy.transfer_data_to(vy);
		uz.transfer_data_to(vz);


		if (Rav != 0)
		{
			poisson_equation_pulsation_stream_function();
			//poisson_equation_pulsation_velocity_x();
			//poisson_equation_pulsation_velocity_y();
			//make_vibr_buffer();
		}



		if (EVERY(1))
		{
			double t_ = timer.get("total");
			double Ek, Vmax;
			statistics(Ek, Vmax);
			cout << iter << ", t = " << total_time << ", t_c = " << t_ << " (" << t_/60 << "), Ek = " << Ek << ", Vmax = " << Vmax << endl;
		}

		
		if (EVERY(1))
		{
			double t_ = timer.get("total");
			double Ek, Vmax;
			statistics(Ek, Vmax);
			
			temporal.write_header("Ra, Rav, total_time, compute_time, Ek, Vmax, check_Ek.dif, check_Ek.long_dif, iter_div, iter_p");
			temporal.write({ Ra, Rav, total_time, t_, Ek, Vmax, check_Ek.dif, check_Ek.long_dif, double(iter_div), double(iter_p)});



			if (check_Ek.stop(Ek, false))
			{
				stop_signal = 1;
			}
		}

		if (EVERY(10))
		{
			write_fields("F.dat");
			write_section_xz(ny / 2);
			back.save_fields();
		}


		if (stop_signal)
		{
			finalize();
			break;
		}
	}
}