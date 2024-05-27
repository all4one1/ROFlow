#include "FlowSolver.h"
#include "toolkit.h"

FlowSolver::FlowSolver(Configuration config)
{
	nx = config.nx;
	ny = config.ny;
	nz = config.nz;

	N = nx * ny * nz;
	
	off = ny == 1 ? 0 : nx;
	off2 = nz == 1 ? 0 : nx * ny;
	
	Lx = config.Lx;
	Ly = config.Ly;
	Lz = config.Lz;
	hx = config.hx;
	hy = config.hy;
	hz = config.hz;
	tau = config.tau;

	Sx = 0.0;
	Sy = 0.0;
	Sz = 0.0;
	dV = hx * hy * hz;

	dim = 3;
	if (nz > 1)
	{
		dim = 3;
		Sx = hy * hz;
		Sy = hx * hz;
		Sz = hx * hy;
		dV = hx * hy * hz;

		Nvx = (nx + 1) * ny * nz;
		Nvy = nx * (ny + 1) * nz;
		Nvz = nx * ny * (nz + 1);
	}

	if (ny > 1 && nz == 1)
	{
		dim = 2;
		Sx = hy;
		Sy = hx;
		Sz = 0.0;
		dV = hx * hy;

		Nvx = (nx + 1) * ny;
		Nvy = nx * (ny + 1);
		Nvz = 0;
	}

	if (ny == 1 && nz == 1)
	{
		dim = 1;
		Sx = 1.0;
		Sy = 0.0;
		Sz = 0.0;
		dV = hx;

		Nvx = (nx + 1);
		Nvy = 0;
		Nvz = 0;
	}

	NV  = Nvx + Nvy + Nvz;



	stride  = Nvx;
	stride2 = Nvx + Nvy;

	SM.resize(NV);
	alloc(&V, NV);
	alloc(&V0, NV);
	alloc(&U, NV);
	alloc(&U0, NV);
	alloc(&B, NV);
	alloc(&Ap, NV);

	sm.resize(N);
	alloc(&b, N);

	vx_ = &V[0];
	vy_ = &V[stride];
	vz_ = &V[stride2];

	vx = Velocity(Component::x, nx, ny, nz, hx, hy, hz, false, vx_);
	vy = Velocity(Component::y, nx, ny, nz, hx, hy, hz, false, vy_);
	vz = Velocity(Component::z, nx, ny, nz, hx, hy, hz, false, vz_);

	ux = Velocity(Component::x, nx, ny, nz, hx, hy, hz, false, &U[0]);
	uy = Velocity(Component::y, nx, ny, nz, hx, hy, hz, false, &U[stride]);
	uz = Velocity(Component::z, nx, ny, nz, hx, hy, hz, false, &U[stride2]);


	vx_prime = Velocity(Component::x, nx, ny, nz, hx, hy, hz);
	vy_prime = Velocity(Component::y, nx, ny, nz, hx, hy, hz);
	vz_prime = Velocity(Component::z, nx, ny, nz, hx, hy, hz);

	ddx = Velocity(Component::x, nx, ny, nz, hx, hy, hz);
	ddy = Velocity(Component::y, nx, ny, nz, hx, hy, hz);
	ddz = Velocity(Component::z, nx, ny, nz, hx, hy, hz);

	buffer = ScalarVariable(nx, ny, nz, hx, hy, hz); 
	bufferU = ScalarVariable(nx, ny, nz, hx, hy, hz);
	P = ScalarVariable(nx, ny, nz, hx, hy, hz);
	P0 = ScalarVariable(nx, ny, nz, hx, hy, hz);
	p_prime = ScalarVariable(nx, ny, nz, hx, hy, hz);
	
	T = ScalarVariable(nx, ny, nz, hx, hy, hz);
	T0 = ScalarVariable(nx, ny, nz, hx, hy, hz);

	grav.set_directly(0, 1, 0);

}


void FlowSolver::set_boundary_for_pressure(Side s, MathBoundary t, double v)
{
	P.boundary.set_boundary(s, t, v);	
	P0.boundary[s] = P.boundary[s];
}
void FlowSolver::set_boundary_for_temperature(Side s, MathBoundary t, double v)
{
	T.boundary.set_boundary(s, t, v);
	T0.boundary[s] = T.boundary[s];
}
void FlowSolver::set_boundary_for_velocity(Component c, Side s, MathBoundary t, double v)
{
	if (c == Component::x)
	{
		vx.boundary.set_boundary(s, t, v);
		ux.boundary.set_boundary(s, t, v);
	}	
	if (c == Component::y)
	{
		vy.boundary.set_boundary(s, t, v);
		uy.boundary.set_boundary(s, t, v);
	}	
	if (c == Component::z)
	{
		vz.boundary.set_boundary(s, t, v);
		uz.boundary.set_boundary(s, t, v);
	}
}
void FlowSolver::set_boundary(Side s, PhysBoundary t, double v)
{
	phys_bc[s] = t;
	if (t == PhysBoundary::Closed)
	{
		P.boundary.set_boundary(s, MathBoundary::Neumann, 0.0);
		P0.boundary[s] = P.boundary[s];

		p_prime.boundary.set_boundary(s, MathBoundary::Neumann, 0.0);

		vx.boundary.set_boundary(s, MathBoundary::Dirichlet, 0.0);
		vy.boundary[s] = vz.boundary[s] = vx.boundary[s];

		ux.boundary.set_boundary(s, MathBoundary::Dirichlet, 0.0);
		uy.boundary[s] = uz.boundary[s] = ux.boundary[s];
	}
	if (t == PhysBoundary::Inlet || t == PhysBoundary::Outlet)
	{
		P.boundary.set_boundary(s, MathBoundary::Dirichlet, v);
		P0.boundary[s] = P.boundary[s];

		p_prime.boundary.set_boundary(s, MathBoundary::Dirichlet, 0.0);

		vx.boundary.set_boundary(s, MathBoundary::Neumann, 0.0);
		vy.boundary[s] = vz.boundary[s] = vx.boundary[s];

		ux.boundary.set_boundary(s, MathBoundary::Neumann, 0.0);
		uy.boundary[s] = uz.boundary[s] = ux.boundary[s];
	}
}
void FlowSolver::set_linear_pressure(double px0, double pxN, double py0, double pyN, double pz0, double pzN)
{
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				double x = hx * i + 0.5 * hx;
				double y = hy * j + 0.5 * hy;
				double z = hz * k + 0.5 * hz;
				P(i, j, k) = P0(i, j, k) =
					(pxN - px0) / Lx * x + px0 +
					(pyN - py0) / Ly * y + py0 +
					(pzN - pz0) / Lz * z + pz0;
				P(i, j, k) *= 1.0;
			}
		}
	}

}

double FlowSolver::check_div()
{
	double div = 0.0;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++)
			{
				div += (ux(i + 1, j, k) - ux(i, j, k)) / hx; 
				div += (uy(i, j + 1, k) - uy(i, j, k)) / hy;
				div += (uz(i, j, k + 1) - uz(i, j, k)) / hz;
			}
		}
	}
	return abs(div);
}

void FlowSolver::statistics(double &Ek, double &Vmax)
{
	Ek = 0;
	Vmax = 0;
	double ux_ = 0, uy_ = 0, uz_ = 0;
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++)
			{
				ux_ = 0.5 * (ux(i + 1, j, k) + ux(i, j, k));
				if (dim > 1) uy_ = 0.5 * (uy(i, j + 1, k) + uy(i, j, k));
				if (dim > 2) uz_ = 0.5 * (uz(i, j, k + 1) + uz(i, j, k));

				double V = sqrt(pow(ux_, 2) + pow(uy_, 2) + pow(uz_, 2));
				if (V > Vmax) Vmax = V;
				Ek += 0.5 * (pow(V, 2)) * dV;
			}
		}
	}

}