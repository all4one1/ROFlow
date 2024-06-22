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
	config.set_uniform(20, "3D", 1.0, 1.0, 1.0); 
	//config.set_cell_number(3, 3); 	config.set_domain_LxLyLz(1, 2, 1);
	config.tau = 1e+100;
	config.tau = 0.01;
	

	FlowSolver solver(config);

	solver.K = 0;
	solver.Pr = 7;
	solver.Rav = 4000;

	solver.grav.set_directly_xyz(0, 1, 0);
	solver.vibr.set_directly_xyz(1, 0, 0);

	//solver.set_period_pair(Side::west, Side::east);			    
	solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);
	
	solver.solve_heat_equation(1000);
	solver.T.write3D();




	//solver.T.set_linear();
	//for (int i = 0; i < 10; i++)
	//{
	//	solver.poisson_equation_pulsation_stream_function();
	//	//solver.poisson_equation_pulsation_velocity_x();
	//}
	//for (int k = 0; k < solver.nz; k++) {
	//	for (int j = 0; j < solver.ny; j++) {
	//		for (int i = 0; i < solver.nx; i++) {
	//			int l = i + solver.off * j + solver.off2 * k;
	//			double VX = +solver.SF.get_dy(Side::center, i, j, k);
	//			double VY = -solver.SF.get_dx(Side::center, i, j, k);
	//			solver.VibrX(i, j, k) = VX;
	//			solver.VibrY(i, j, k) = VY;
	//			solver.buffer[l] = VX /** solver.vibr.x + VY * solver.vibr.y*/;
	//		}
	//	}
	//}



	config.show_parameters();

	solver.timer.start("total");

	double R = 0;
	solver.Ra = 0;
	
	for (double R = 3000; R > 1000; R = R - 50)
	{
		solver.reset();
		solver.Rav = R;
		#define TT(i) size_t(1.0 / config.tau) * i
		solver.solve_system(TT(0));  //size_t(1.0 / config.tau) * 2000
		solver.finalize();
	}


	solver.timer.end("total");
	solver.timer.write_info();

	
	solver.ux.show_max_min();
	solver.uy.show_max_min();	
	cout << solver.uy(0, config.ny / 2, 0) << endl;
	//solver.SM.save_full_matrix_with_rhs(5, solver.B);

	cout << "end" << endl;
	return 0;
}
