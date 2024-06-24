#define CODEVERSION 240624

#include "toolkit.h"
#include "Configuration.h"
#include "FlowSolver.h"

#include <iostream>
using std::cout;
using std::endl;

int main(int argc, char** argv)
{
	Configuration config;

	
	config.set_uniform(20, "3D", 2.0, 1.0, 2.0); // n, Dim, Lx, Ly, Lz
	config.tau = 0.1;
	

	FlowSolver solver(config);
	config.show_parameters();


	solver.K = 0;
	solver.Pr = 7;
	solver.Rav = 0;
	solver.grav.set_directly_xyz(0, 1, 0);
	solver.vibr.set_directly_xyz(0, 1, 0);

	solver.set_period_pair(Side::west, Side::east);			    
	solver.set_period_pair(Side::front, Side::back);
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);


	solver.timer.start("total");

	double R = 4000;
	//for (double R = 3000; R > 2000; R = R - 100)
	{
		solver.Ra = R;
		solver.solve_system(200, true);  
		solver.reset();
	}


	solver.timer.end("total");
	solver.timer.write_info();

	if (solver.N < 50)
		solver.SM.save_full_matrix_with_rhs(5, solver.B);

	cout << "end" << endl;
	return 0;
}
