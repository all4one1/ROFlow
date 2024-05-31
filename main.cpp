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
	//config.set_domain_LxLyLz(1,1, 0.4);
	//config.set_cell_number(nx, ny, nz);

	config.set_uniform(20, "2D", 10.0, 1.0);


	config.tau = 1e+100;
	config.tau = 0.001;
	

	FlowSolver solver(config);
	config.show_parameters();


	solver.Ra = 2100;
	solver.Pr = 1;
	solver.grav.x = 0;
	solver.grav.y = 1;
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);
	solver.T0.boundary = solver.T.boundary;

	//solver.set_boundary(Side::west, PhysBoundary::Inlet, 1);
	//solver.set_boundary(Side::east, PhysBoundary::Outlet, 0);


	//solver.solve_simple(1000);

	solver.timer.start("total");

	double R = 1750;
	//for (double R = 2000; R > 1500; R = R - 50)
	{
		solver.reset();
		solver.Ra = R;
		solver.solve_system(3000);
		solver.finalize();
	}

	solver.timer.end("total");

	solver.finalize();



	solver.uy.show_max_min();
	//solver.SM.save_full_matrix_with_rhs(5, solver.B);

	cout << "end" << endl;

	return 0;
}
