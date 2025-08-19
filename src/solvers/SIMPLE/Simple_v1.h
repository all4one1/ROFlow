#pragma once
#include "../FromOuterSparse/SparseMatrix.h"
#include "tools/types.h"
#include "tools/common_tools.h"
#include "newvar.h"

#include "solvers/MatrixMaker.h"
#include "solvers/slau.h"

#include <vector>
#include <map>


struct SIMPLE : Configuration
{
	unsigned int Nvx = 0, Nvy = 0, Nvz = 0, NV = 0, stride = 0, stride2 = 0;
	Variable T, T0;
	
	FuncTimer ftimer;
	MatrixMaker mk;

	SIMPLE() : Configuration() {}
	SIMPLE(const Configuration& config);



	//void form_matrix_for_laplace_centered(SparseMatrix &M, Variable &F, double ax, double ay, double phys_coef, double time_coef)
	//{
	//	ftimer.start("matrix_formation");

	//	M.reset();
	//	struct Contribution
	//	{
	//		std::map<int, double> m;
	//		int l, i, j, k, nx, ny, nz, off, off2;
	//		Contribution(Variable& F, int l_, int i_ = 0, int j_ = 0, int k_ = 0)
	//			: l(l_), i(i_), j(j_), k(k_)
	//		{
	//			off = F.off;
	//			off2 = F.off2;
	//			nx = F.nx;
	//			ny = F.ny;
	//			nz = F.nz;
	//		}
	//		int index(Side side)
	//		{
	//			switch (side)
	//			{
	//			case Side::center:
	//				return l;
	//				break;
	//			case Side::west:
	//				return l - 1;
	//				break;
	//			case Side::east:
	//				return l + 1;
	//				break;
	//			case Side::south:
	//				return l - off;
	//				break;
	//			case Side::north:
	//				return l + off;
	//				break;
	//			case Side::front:
	//				return l - off2;
	//				break;
	//			case Side::back:
	//				return l + off2;
	//				break;
	//			default:
	//				return -1;
	//				break;
	//			}
	//			return -1;
	//		};
	//		int index_period(Side side)
	//		{
	//			switch (side)
	//			{
	//			case Side::center:
	//				return l;
	//				break;
	//			case Side::west:
	//				return (nx - 1) + off * j + off2 * k;
	//				break;
	//			case Side::east:
	//				return (0) + off * j + off2 * k;
	//				break;
	//			case Side::south:
	//				return i + (ny - 1) * off + off2 * k;
	//				break;
	//			case Side::north:
	//				return i + off * 0 + off2 * k;
	//				break;
	//			case Side::front:
	//				return i + j * off + (nz - 1) * off2;
	//				break;
	//			case Side::back:
	//				return i + j * off + 0 * off2;
	//				break;
	//			default:
	//				return -1;
	//				break;
	//			}
	//			return -1;
	//		};

	//		double& operator()(Side side, bool period = false)
	//		{
	//			if (period == false) return (*this)[index(side)];
	//			else  return (*this)[index_period(side)];
	//		}
	//		double& operator[](int l)
	//		{
	//			return m[l];
	//		}
	//		std::map<int, double>& get_map()
	//		{
	//			return m;
	//		}
	//	};


	//	for (int k = 0; k < nz; k++) {
	//		for (int j = 0; j < ny; j++) {
	//			for (int i = 0; i < nx; i++) {
	//				int l = i + offset * j + offset2 * k;
	//				Contribution a(F, l, i, j, k);

	//				auto flux = [&F, &a, i, j, k, l](Side side, double S, double h, double coef, bool border)
	//				{
	//					if (border)
	//					{
	//						if (F.boundary.type(side) == MathBoundary::Neumann)
	//						{

	//						}
	//						if (F.boundary.type(side) == MathBoundary::Dirichlet)
	//						{
	//							a(Side::center) += coef * S / (0.5 * h);
	//						}
	//						if (F.boundary.type(side) == MathBoundary::periodic)
	//						{
	//							a(Side::center) += coef * S / h;
	//							a(side, true) += -coef * S / h;
	//						}
	//					}
	//					else
	//					{
	//						a(Side::center) += coef * S / h;
	//						a(side) += -coef * S / h;
	//					}
	//				};


	//				flux(Side::west, Sx, hx, phys_coef, i == 0);
	//				flux(Side::east, Sx, hx, phys_coef, i == nx - 1);

	//				if (dim > 1)
	//				{
	//					flux(Side::south, Sy, hy, phys_coef, j == 0);
	//					flux(Side::north, Sy, hy, phys_coef, j == ny - 1);
	//				}
	//				if (dim > 2)
	//				{
	//					flux(Side::front, Sz, hz, phys_coef, k == 0);
	//					flux(Side::back, Sz, hz, phys_coef, k == nz - 1);
	//				}

	//				a(Side::center) += time_coef;
	//				M.add_line_with_map(a.get_map(), l);
	//			}
	//		}
	//	}


	//	ftimer.end("matrix_formation");
	//}
};

