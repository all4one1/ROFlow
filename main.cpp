#define DEBUG

#include "toolkit.h"
#include "Configuration.h"
#include "Variable.h"
#include "FlowSolver.h"
#include "IterativeSolver.h"
#include <iostream>
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
	Configuration config;

	int nx = 20;
	int ny = 20;
	int nz = 1;
	config.set_domain_LxLyLz(1, 1, 0.4);
	config.set_cell_number(nx, ny, nz);

	config.tau = 1e+100;
	config.tau = 0.001;
	config.dim = 3;

	FlowSolver solver(config);
	config.show_parameters();


	solver.Ra = 2550;
	solver.Pr = 1;
	solver.grav.x = 0;
	solver.grav.y = 1;
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);
	solver.T0.boundary = solver.T.boundary;


	

	solver.solve_system(1000000);


	solver.uy.show_max_min();
	solver.uy.write3D(0, 1, 0);

	solver.SM.recover_full_with_rhs(5, solver.B);

	cout << "end" << endl;
	return 0;
}
