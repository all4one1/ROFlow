#pragma once
#include "FromOuterSparse/SparseMatrix.h"
#include "IterativeSolver.h"
#include "Configuration.h"
#include "Variable.h"
#include "Extras.h"
#include <map>


#define REDUCED 0


#define timer(name, func) {double start_ = clock();  func;  double end_ = clock(); m_timer[name] += (end_ - start_) / CLOCKS_PER_SEC;}

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
	Velocity vx,  vy, vz;
	Velocity vx_prime, vy_prime, vz_prime;
	Velocity ux, uy, uz;
	Velocity ddx, ddy, ddz;
	ScalarVariable P, P0, T, T0, p_prime, buffer, bufferU;
	double* vx_, * vy_, * vz_, * p;

	double* p0, * dp, * dp0;
	double* V, * V0, *U, *U0, * B, *Ap,  *b;
	SparseMatrix SM, sm;
	int nx, ny, nz, off, off2, N;
	int dim, stride, stride2, Nvx, Nvy, Nvz, NV;
	double hx, hy, hz, Sx, Sy, Sz, Lx, Ly, Lz, dV;
	double Re = 1, Ra = 1, Pr = 1, Gr = 1;
	double P_in = 0, P_out = 0, dP = 0.0;
	StaticVector grav, vibr, dens;
	double tau, total_time = 0.0, compute_time = 0.0;
	size_t iter = 0;
	int iter_limit = 1000;
	std::map <Side, PhysBoundary> phys_bc;
	std::map <std::string, double> m_timer;
	IterativeSolver itsol;
	Checker kinetic_check;
	StateOut stats;

	FlowSolver(Configuration config);

	void set_preset_geometry();
	void set_boundary_for_pressure(Side s, MathBoundary t, double v);
	void set_boundary_for_temperature(Side s, MathBoundary t, double v);
	void set_boundary_for_velocity(Component c, Side s, MathBoundary t, double v = 0.0);
	void set_linear_pressure(double px0, double pxN, double py0 = 0, double pyN = 0, double pz0 = 0, double pzN = 0);
	void set_boundary(Side s, PhysBoundary t, double v = 0.0);

	void solve_projection(int steps_at_ones = 1);
	void projection_quasi_velocity();
	void projection_poisson_for_pressure();
	void projection_true_velocity();

	void solve_heat_equation_explicitly(int steps_at_ones = 1);
	void form_matrix_for_heat_equation(SparseMatrix& M, double* b);
	void form_rhs_for_heat_equation(double* b, bool reset = false);

	void solve_heat_equation(int steps_at_ones = 1);

	void form_u_matrix_for_simple(SparseMatrix& M, double* b);
	void form_u_rhs_for_simple(double* b, bool reset = false);
	void poisson_equation_for_p_prime();
	void poisson_equation_for_p_prime2();
	void correction2();
	void correction_for_simple();
	void solve_simple(int steps_at_ones = 1);
	void solve_system(int steps_at_ones = 1);
	double check_div();

	void statistics(double& Ek, double& Vmax);


	void form_matrix_test(SparseMatrix& M, double* b);
	void form_rhs_test(double* b, bool reset);

	void form_big_matrix(SparseMatrix& M, double* b);

	void form_big_rhs(double* b, bool reset);



};


