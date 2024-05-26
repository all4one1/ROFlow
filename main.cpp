#define DEBUG

#include "toolkit.h"
#include "Configuration.h"
#include "Variable.h"
#include "FlowSolver.h"
#include "IterativeSolver.h"
#include <iostream>
using std::cout;
using std::endl;

#define callNtimes(N, func) for (int i = 0; i < N; i++) {func; cout << "iter: " << i << endl;}
int main(int argc, char** argv)
{
	Configuration config;

	int nx = 20;
	int ny = 20;
	int nz = 1;
	config.set_domain_LxLyLz(1, 1, 0.4);
	config.set_cell_number(nx, ny, nz);

	config.tau = 1e+100;
	config.tau = 0.0005;
	config.dim = 3;

	FlowSolver solver(config);
	config.show_parameters();


	solver.Ra = 5000;
	solver.Pr = 1;
	solver.grav.x = 0;
	solver.grav.y = 1;
	//solver.T.boundary.set_boundary(Side::west, MathBoundary::Dirichlet, 1);
	//solver.T.boundary.set_boundary(Side::east, MathBoundary::Dirichlet, 0);
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);
	//solver.solve_heat_equation();
	//solver.T.set_linear();
	solver.T0.boundary = solver.T.boundary;
	

	solver.solve_system(10000);
	solver.uy.show_max_min();
	solver.uy.write3D(0,1,0);

	solver.SM.recover_full_with_rhs(5, solver.B);

	ofstream w("check.dat");


	


	cout << "end" << endl;
	return 0;
}



//
//auto FF = [&solver, dp](double y)
//{
//	return -solver.Re * dp / solver.Lx / 2.0 * (y * y - y * solver.Ly);
//};
//for (int i = 0; i < solver.ny; i++)
//{
//	double y = solver.ux.y_(i);
//	cout << fixed << y << " " << solver.ux(0, i, 0) << " " << solver.ux(nx / 2, i, 0) << " " << solver.ux(nx, i, 0) << " " << FF(y) << endl;
//}
//
//for (int i = 0; i < solver.nx + 1; i++)
//{
//	double x = solver.ux.x_(i);
//	//cout << fixed << x << " " << solver.ux(i, 0, 0) << endl;
//	//w << fixed << x << " " << solver.ux(i, 0, 0) << endl;
//}
