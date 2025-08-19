#define VERSION 0_0_1

#include "tools/common_tools.h"
#include "tools/project_tools.h"
#include "newvar.h"
#include "solvers/MatrixMaker.h"


//#include "toolkit.h"
//#include "Configuration.h"
//#include "solvers/FlowSolver.h"


#include <iostream>
using std::cout;
using std::endl;



int main(int argc, char** argv)
{
	Configuration config{};

	init_parameters(config);

	Variable T(config);

	SparseMatrix A(config.N), B(config.N), SM(config.N*2);


	MatrixMaker mk(config);



	//FlowSolver solver(config);
	//config.show_parameters();


	//solver.K = 0;
	//solver.Pr = 7;
	//solver.Rav = 0;
	//solver.grav.set_directly_xyz(0, 1, 0);


	////solver.set_period_pair(Side::west, Side::east);			    
	////solver.set_period_pair(Side::front, Side::back);	
	//solver.T.boundary.set_boundary(Side::east, MathBoundary::Neumann, 0);
	//solver.T.boundary.set_boundary(Side::west, MathBoundary::Neumann, 0);
	//solver.T.boundary.set_boundary(Side::south, MathBoundary::Dirichlet, 1);
	//solver.T.boundary.set_boundary(Side::north, MathBoundary::Dirichlet, 0);


	//solver.timer.start("total");

	//double R = 6000;
	////for (double R = 6000; R > 2000; R = R - 100)
	//{
	//	solver.Ra = R;
	//	solver.solve_system(200, true);  
	//	solver.reset();
	//}


	//solver.timer.end("total");
	//solver.timer.write_info();

	//if (solver.N < 50)
	//	solver.SM.save_full_matrix_with_rhs(5, solver.B);

	cout << "end" << endl;
	return 0;
}
