#pragma once

#pragma once
#include <type_traits>
enum discretization { FV = 0, FD = 1 };
enum bc_type { closed = 0, periodic = 1, Dirichlet = 2, Neumann = 3, not_boundary = 4 };
enum domain_config { closed_box, open_tube, periodic_cell };
enum field_name {
	field_T, field_vx, field_vy, field_p, field_p_prime, field_vx_prime, field_vy_prime,
	field_C, field_mu, field_C2, field_mu2, field_C3, field_mu3, field_stream, field_omega, field_count
};
enum enumside { west, east, south, north, front, back, inner };

enum class MathBoundary { Dirichlet, Neumann, periodic, not_boundary, array_value, complex };
enum class Side { center, west, east, south, north, front, back };
enum class Component { x, y, z };

#define ENUM_DEFINED

struct Configuration
{
	double hx, hy, hz, Sx, Sy, Sz, Lx, Ly, Lz, dV;
	double ox, oy, oz;
	double tau, tau_p;
	double alpha, beta, gamma;
	double sinA, cosA, sinB, cosB, sinG, cosG;
	double vibr_x, vibr_y, vibr_z;
	double grav_x, grav_y, grav_z;
	double p_in, p_out;
	double density_x, density_y, density_z;
	double Sc11, Sc12, Sc21, Sc22, psi1, psi2, psiS, Q, Ra, Rad, Rav, K, Pr, Le;
	double Re, eps, Pe, Mach, Mach2, Mach3, M0, A, Cn, Gr, Gr2, Gr3, Sigma1, Sigma2, Sigma3;
	double incr_parameter;
	unsigned int q, dim, nx, ny, nz, N, offset, offset2, Nbytes;
	unsigned int heatflux;
	bc_type xbc, ybc, zbc;
	domain_config domain;
	discretization disc;
};

struct PhysicalValues
{
	double Ek, Vmax, AM, AMx, AMy, AMz;
};

struct StatValues
{
	double ksi_max = 0, ksi_sum = 0, omega_max = 0, omega_sum = 0, C_sum = 0, C_sum_signed = 0;
	double NuTop = 0, NuDown = 0, NuLeft = 0, NuRight = 0, ShrTop = 0, ShrDown = 0;
	double Vmax = 0, Vx = 0, Vy = 0, Ek = 0;
	double Cu = 0, Pe = 0;
};


struct Arrays
{
	double *buffer, * buffer2, * buffer3, * rhs;
	double* p, * p0, * p_prime, * p_prime0, * ux, * uy, * uz, * vx, * vy, * vz;
	double* ux_buf, * uy_buf, * Hvisc;
	double* T, * T0, * C, * C0, *C2, * C20, * C3, * C30, * mu, * mu0, * mu2, * mu20, * mu3, * mu30;
	double* omega, * omega0, * ksi, * ksi0;
	double* src_x, * src_y, * src_z;
	double* U, * U0, * Urhs;
	double* F, * F0, * Frhs;
};
