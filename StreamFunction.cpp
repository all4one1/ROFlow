//#include "FlowSolver.h"
//
//
//void FlowSolver::solve_stream_function(size_t steps_at_ones)
//{
//	if (iter == 0)
//	{
//		ksi = ScalarVariable(nx, ny, nz, hx, hy, hz);
//		ksi0 = ScalarVariable(nx, ny, nz, hx, hy, hz);
//		curl = ScalarVariable(nx, ny, nz, hx, hy, hz);
//		curl0 = ScalarVariable(nx, ny, nz, hx, hy, hz);
//	}
//
//	for (size_t q = 0; q < steps_at_ones; q++)
//	{
//	
//			for (int j = 1; j < ny - 1; j++) {
//				for (int i = 1; i < nx - 1; i++) {
//					curl(i, j) = curl0(i, j); //+
//						//tau / hx / hy / 4.0 * ((ksi(i + 1, j) - ksi(i - 1, j)) * (curl(i, j + 1) - curl(i, j - 1)) - (ksi(i, j + 1) - ksi(i, j - 1)) * (curl(i + 1, j) - curl(i - 1, j))) + // nonlinear term
//						//tau / hx / hy * (curl(i + 1, j) + curl(i - 1, j) + curl(i, j + 1) + curl(i, j - 1) - 4. * curl(i, j)) + // laplacian term
//						//-dt * Pr * Ra * (-(-cosG * a + (tem(i + 1, j) - tem(i - 1, j)) / hx / 2.0 + kappa * (-cosG + (con1(i + 1, j) - con1(i - 1, j)) / hx / 2.0)) * sinA + (-sinG * a + (tem(i, j + 1) - tem(i, j - 1)) / hy / 2.0 + kappa * (-sinG + (con1(i, j + 1) - con1(i, j - 1)) / hy / 2.0)) * cosA) + // gravity term
//						//Rav * Pr * dt * ((-(-ksib(i - 1, j + 1) + ksib(i + 1, j + 1) + ksib(i - 1, j - 1) - ksib(i + 1, j - 1)) / hx / hy * cosb / 0.4 + (ksib(i + 1, j) + ksib(i - 1, j) - 2 * ksib(i, j)) / hx * hx * sinb) * (-sinG * a + (tem(i, j + 1) - tem(i, j - 1)) / hy / 2.0 + Kappa * (-sinG + (con1(i, j + 1) - con1(i, j - 1)) / hy / 2.0)) - (-(ksib(i, j + 1) + ksib(i, j - 1) - 2 * ksib(i, j)) / hy / hy * cosb + (-ksib(i - 1, j + 1) + ksib(i + 1, j + 1) + ksib(i - 1, j - 1) - ksib(i + 1, j - 1)) / hx / hy * sinb / 0.4) * (-cosG * a + (tem(i + 1, j) - tem(i - 1, j)) / hx / 2.0 + Kappa * (-cosG + (con1(i + 1, j) - con1(i - 1, j)) / hx / 2.0))); // vibration term
//
//
//
//				}
//			}
//		
//	}
//
//
//
//}