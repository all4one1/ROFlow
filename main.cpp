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
	config.set_uniform(20, "2D", 1.0, 1.0); 
	//config.set_cell_number(3, 3); 	config.set_domain_LxLyLz(1, 2, 1);
	config.tau = 1e+100;
	config.tau = 0.01;
	

	FlowSolver solver(config);

	solver.K = 0;
	solver.Pr = 10;
	solver.Rav = 0;
	solver.vibr.set_directly_xyz(1, 0, 0);


	solver.set_period_pair(Side::west, Side::east);			    
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);
	solver.grav.set_directly_xyz(0, 1, 0);


	//solver.set_period_pair(Side::south, Side::north);
	//solver.T.boundary.set_boundary(Side::west, MathBoundary::Dirichlet, 1);
	//solver.T.boundary.set_boundary(Side::east, MathBoundary::Dirichlet, 0);
	//solver.grav.set_directly_xyz(1, 0, 0);
	


	config.show_parameters();

	solver.timer.start("total");

	double R = 5000;
	//for (double R = 5000; R > 1000; R = R - 50)
	{
		solver.reset();
		solver.Ra = R;
		solver.solve_system(3000);  //size_t(1.0 / config.tau) * 2000
		solver.finalize();
	}



	solver.timer.end("total");
	solver.timer.write_info();

	solver.uy.write3D();
	solver.ux.show_max_min();
	solver.uy.show_max_min();	
	cout << solver.uy(0, config.ny / 2, 0) << endl;
	//solver.SM.save_full_matrix_with_rhs(5, solver.B);

	cout << "end" << endl;
	return 0;
}
