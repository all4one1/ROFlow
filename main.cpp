#define CODEVERSION 010624

#include "toolkit.h"
#include "Configuration.h"
#include "FlowSolver.h"

#include <iostream>
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
	Configuration config;

	// n, Dim, Lx, Ly, Lz
	config.set_uniform(20, "2D", 15.0, 1.0); 
	//config.tau = 1e+100;
	config.tau = 0.0005;
	

	FlowSolver solver(config);

	solver.K = 0;
	solver.Pr = 10;
	solver.Rav = 320;

	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);

	solver.grav.set_directly_xyz(0, 1, 0);
	solver.vibr.set_directly_xyz(1, 0, 0);


	solver.timer.start("total");

	double R = 2000;
	for (double R = 2000; R > 1000; R = R - 50)
	{
		solver.reset();
		solver.Ra = R;
		solver.solve_system();  //size_t(1.0 / config.tau) * 2000
		solver.finalize();
	}

	solver.timer.end("total");
	solver.timer.write_info();

	solver.uy.show_max_min();
	//solver.SM.save_full_matrix_with_rhs(5, solver.B);

	cout << "end" << endl;
	return 0;
}
