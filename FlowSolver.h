#pragma once
#include "FromOuterSparse/SparseMatrix.h"
#include "IterativeSolver.h"
#include "Configuration.h"
#include "Variable.h"
#include "Extras.h"
#include <map>

#define REDUCED 0

using namespace std;
enum class PhysBoundary
{
	Closed,
	Inlet,
	Outlet,
	Periodic,
	HeatFixed,
	HeatFlux
};


struct FlowSolver
{
	
public:
	//physical fields:
	Velocity vx, vy, vz;
	Velocity ux, uy, uz;
	ScalarVariable P, P0, T, T0, C, C0, SF, SF0, p_prime;
	ScalarVariable buffer, bufferU, bufferVibr;
	ScalarVariable ksi, ksi0, curl, curl0;
	ScalarVariable VibrX, VibrX0, VibrY, VibrY0;
	double* vx_, * vy_, * vz_;

	double* V, * V0, *U, *U0, * B, *Ap,  *b;
	SparseMatrix SM, SMt, SMc, SMtest;
	int nx, ny, nz, off, off2, N;
	int dim, stride, stride2, Nvx, Nvy, Nvz, NV;
	double hx, hy, hz, Sx, Sy, Sz, Lx, Ly, Lz, dV;
	double Re = 1, Ra = 0, Rav = 0, Pr = 1, Gr = 0, Le = 1, K = 0;
	double P_in = 0, P_out = 0, dP = 0.0;
	double alpha_relax = 1.0;
	StaticVector grav, vibr, dens;
	double tau, total_time = 0.0;
	size_t iter = 0, iter_div = 0, iter_p = 0, iter_limit = 1000, bytes_allocated = 0;
	int stop_signal = 0;
	std::map <Side, PhysBoundary> phys_bc;
	FuncTimer timer;
	IterativeSolver itsol;
	Checker check_Ek;
	StateOut final, temporal;
	BackUp back;

	FlowSolver(Configuration config);

	void set_preset_geometry();
	void set_boundary_for_pressure(Side s, MathBoundary t, double v);
	void set_boundary_for_temperature(Side s, MathBoundary t, double v);
	void set_boundary_for_velocity(Component c, Side s, MathBoundary t, double v = 0.0);
	void set_linear_pressure(double px0, double pxN, double py0 = 0, double pyN = 0, double pz0 = 0, double pzN = 0);
	void set_boundary(Side s, PhysBoundary t, double v = 0.0);
	void set_period_pair(Side s1, Side s2);

	void solve_projection(int steps_at_ones = 1);
	void projection_quasi_velocity();
	void projection_poisson_for_pressure();
	void projection_true_velocity();

	void form_matrix_for_heat_equation(SparseMatrix& M, double* b);
	void form_rhs_for_heat_equation(double* b, bool reset = false);
	void form_matrix_for_concentration_equation(SparseMatrix& M, double* b);
	void form_rhs_for_concentration_equation(double* b, bool reset = false);
	void solve_heat_equation(size_t steps_at_ones = SIZE_MAX);
	void solve_heat_equation_explicitly(int steps_at_ones = 1);

	void guessed_velocity_for_simple();
	void poisson_equation_for_p_prime();
	void correction_for_simple();
	void solve_simple(size_t steps_at_ones = SIZE_MAX);
	void solve_system(size_t steps_at_ones = SIZE_MAX, bool inTimeUnits = false);
	
	double Aii(Velocity &v, int i, int j, int k);
	double check_div();
	double check_div2();
	void statistics(double& Ek, double& Vmax);
	void write_fields(std::string path = "results\\field.dat");
	void write_section_xz(int j, std::string path = "results\\section_xy.dat");
	void recover();

	void finalize();
	void reset();


	void form_matrix_test2(SparseMatrix& M, double* b);

	void form_rhs_test2(double* b, bool reset);

	void form_matrix_test(SparseMatrix& M, double* b);
	void form_rhs_test(double* b, bool reset);
	void correction_for_simple_test();
	
	void form_matrix_Uxyz(SparseMatrix& M, double* b);
	void form_rhs_Uxyz(double* b, bool reset);

	void poisson_equation_pulsation_stream_function();
	void poisson_equation_pulsation_velocity_x();
	void poisson_equation_pulsation_velocity_y();
	void make_vibr_buffer();

	void solve_stream_function(size_t steps_at_ones = SIZE_MAX);
};


