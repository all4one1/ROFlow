//#include "FlowSolver.h"
//
//
//void FlowSolver::solve_projection(int steps_at_ones)
//{
//	for (int q = 0; q < steps_at_ones; q++)
//	{
//		iter++;
//		//projection_quasi_velocity();
//		//projection_poisson_for_pressure();
//		//projection_true_velocity();
//		auto FFj = [this](int j)
//		{
//			double y = ux.y_(j);
//			return -Re * dP / Lx / 2.0 * y * (y - Ly);
//		};
//		
//
//		for (int k = 0; k < nz; k++) {
//			for (int j = 0; j < ny; j++) {
//				for (int i = 1; i < nx; i++)
//				{
//					int l = i + off * j + off2 * k;
//
//					double V = dV;
//
//					if (j == 0 || j == ny - 1) V = dV * 0.75;
//
//					ux(i, j, k) = vx(i, j, k) + tau / (V) *
//						(
//							+vx.get_dy_for_vx_cell(Side::north, i, j, k) * Sy / Re
//							- vx.get_dy_for_vx_cell(Side::south, i, j, k) * Sy / Re
//
//							- (P(i, j, k) - P(i - 1, j, k)) / hx * V
//							);
//							
//						
//				}
//
//
//
//				ux(0, j, k) = ux(1, j, k) - hx * (P(0, j, k) - 1);
//				ux(nx, j, k) = ux(nx - 1, j, k) - hx * (0 - P(nx - 1, j, k));
//			}
//		}
//
//
//		for (int k = 0; k < nz; k++) {
//			for (int j = 0; j < ny; j++) {
//				for (int i = 0; i <= nx; i++)
//				{
//					vx(i, j, k) = ux(i, j, k);
//				}
//			}
//		}
//
//	}
//}
//void FlowSolver::projection_quasi_velocity()
//{
//	auto FluxVX = [this](Velocity& VX, Velocity& VY, Velocity& VZ, int i, int j, int k)
//	{
//		double FX = Sx * (VX.get_for_vx_cell(Side::east, i, j, k) * VX.get_for_vx_cell(Side::east, i, j, k)
//			- VX.get_for_vx_cell(Side::west, i, j, k) * VX.get_for_vx_cell(Side::west, i, j, k));
//		double FY = Sy * (VX.get_for_vx_cell(Side::north, i, j, k) * VY.get_for_vx_cell(Side::north, i, j, k)
//			- VX.get_for_vx_cell(Side::south, i, j, k) * VY.get_for_vx_cell(Side::south, i, j, k));
//		double FZ = Sz * (VX.get_for_vx_cell(Side::back, i, j, k) * VZ.get_for_vx_cell(Side::back, i, j, k)
//			- VX.get_for_vx_cell(Side::front, i, j, k) * VZ.get_for_vx_cell(Side::front, i, j, k));
//		return FX + FY + FZ;
//	};
//
//
//	//vx
//	for (int k = 0; k < nz; k++) {
//		for (int j = 0; j < ny; j++) {
//			for (int i = 1; i < nx; i++)
//			{
//				int l = i + off * j + off2 * k;
//
//				ux(i, j, k) = vx(i, j, k) + tau / (dV) *
//					(
//						//-FluxVX(vx, vy, vz, i, j, k)
//						//+ vx.get_dx_for_vx_cell(Side::east, i, j, k) * Sx / Re
//						//- vx.get_dx_for_vx_cell(Side::west, i, j, k) * Sx / Re
//
//						+ vx.get_dy_for_vx_cell(Side::north, i, j, k) * Sy / Re
//						- vx.get_dy_for_vx_cell(Side::south, i, j, k) * Sy / Re
//
//						+ vx.get_dz_for_vx_cell(Side::back, i, j, k) * Sz / Re
//						- vx.get_dz_for_vx_cell(Side::front, i, j, k) * Sz / Re
//
//						);
//
//
//			}
//			ux(0, j, k) = ux(1, j, k);
//			ux(nx, j, k) = ux(nx - 1, j, k);
//		}
//	}
//}
//void FlowSolver::projection_poisson_for_pressure()
//{
//	double tau_p = 0.1 * pow(fmin(fmin(hx, hy), hz), 2);
//	int iter_p = 0;
//	double res, res0 = 0;
//	double eps = 1, eps0 = 1e-6;
//
//	auto laplace = [this](ScalarVariable& f, int i, int j = 0, int k = 0)
//	{
//		double lapl = 0.0;
//		lapl += Sx * (f.get_dx(Side::east, i, j, k) - f.get_dx(Side::west, i, j, k));
//		lapl += Sy * (f.get_dy(Side::north, i, j, k) - f.get_dy(Side::south, i, j, k));
//		lapl += Sz * (f.get_dz(Side::back, i, j, k) - f.get_dz(Side::front, i, j, k));
//		return lapl;
//	};
//	auto divU = [this](int i, int j = 0, int k = 0)
//	{
//		double div = 0.0;
//		div += Sx * (ux.get_for_centered_cell(Side::east, i, j, k) - ux.get_for_centered_cell(Side::west, i, j, k));
//		div += Sy * (uy.get_for_centered_cell(Side::north, i, j, k) - uy.get_for_centered_cell(Side::south, i, j, k));
//		div += Sz * (uz.get_for_centered_cell(Side::back, i, j, k) - uz.get_for_centered_cell(Side::front, i, j, k));
//		return div;
//	};
//
//	while (eps > eps0 * res0)
//	{
//		if (iter_p > iter_limit) break;
//		res = 0;
//		for (int k = 0; k < nz; k++) {
//			for (int j = 0; j < ny; j++) {
//				for (int i = 0; i < nx; i++) {
//					int l = i + off * j + off2 * k;
//
//					//P[l] = P0[l] + tau_p / dV * (fullFlux(P0, i, j, k));
//
//					P[l] = P0[l] + tau_p / dV * (laplace(P0, i, j, k) - divU(i, j, k) / tau);
//
//					res = res + abs(P[l]);
//				}
//			}
//		}
//		eps = abs(res - res0);
//		res0 = res;
//
//		ScalarVariable::swap(P, P0);
//		//for (int l = 0; l < N; l++)			P0[l] = P[l];
//		iter_p = iter_p + 1;
//		if (iter_p % 10000 == 0) cout << "iter_p: " << iter_p << " res: " << res << " eps: " << eps << endl;
//		//if (k > 1000) break;
//	}
//	iter_limit = iter_p;
//	//cout << k << endl;
//}
//void FlowSolver::projection_true_velocity()
//{
//	for (int k = 0; k < nz; k++) {
//		for (int j = 0; j < ny; j++) {
//			for (int i = 1; i < nx; i++) {
//
//				int l = i + off * j + off * 2;
//				vx(i, j, k) = ux(i, j, k) - tau * Sx / dV * (P(i, j, k) - P(i - 1, j, k));
//			}
//			vx(0, j, k) = vx(1, j, k);
//			vx(nx, j, k) = vx(nx - 1, j, k);
//		}
//	}
//
//}
