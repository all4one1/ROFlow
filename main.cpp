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

	int nx = 40;
	int ny = 40;
	int nz = 1;
	config.set_domain_LxLyLz(1,1, 0.4);
	config.set_cell_number(nx, ny, nz);

	config.tau = 1e+100;
	config.tau = 0.001;
	config.dim = 3;

	FlowSolver solver(config);
	config.show_parameters();


	solver.Ra = 5000;
	solver.Pr = 1;
	solver.grav.x = 0;
	solver.grav.y = 1;
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);
	solver.T0.boundary = solver.T.boundary;

	//solver.set_boundary(Side::west, PhysBoundary::Inlet, 1);
	//solver.set_boundary(Side::east, PhysBoundary::Outlet, 0);


	//solver.solve_simple(1000);

	solver.solve_system(10000);
	solver.write_fields();
	solver.finalize();


	solver.uy.show_max_min();
	//solver.SM.recover_full_with_rhs(5, solver.B);

	cout << "end" << endl;

	return 0;
}
