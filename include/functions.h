#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

#define _USE_MATH_DEFINES

#include "types.h"
#include "border.h"
#include "globals.h"
#include <set>

static const int FIRST_ORDER 		= 0;
static const int SECOND_ORDER 		= 1;
static const int FIRST_ORDER_NS		= 2;
static const int SECOND_ORDER_NS	= 3;

static const int IDEAL_GAS			= 0;
static const int ABEL_DUPRE			= 1;
static const int POWDER_EQ			= 2;
static const int PISTON				= 3;

double euler_Usn(double P_sn, double S, double F, double dt, double U_prev);

double euler_Xsn(double x_prev, double U_new);

void euler_proj_broder(double array[5], int j, double Xsn, double dx,
		double dr);

double euler_bar_Vx(cell2d& cell, int n, int i, int j, double dt, double dx,
		double dr, int var);

double euler_bar_Vr(cell2d& cell, int n, int i, int j, double dt, double dx,
		double dr, int var);

void rotateVectors(double& Vx, double& Vr, LineAngle2D angle);

double * smoothSpeed(double * Vx, double * Vr, LineAngle2D angle);

double euler_bar_e(cell2d& cell, int n, int i, int j, double dt, double dx,
		double dr, int var);

double lagrange_rho(gasCell * cell, gasCell * prevCell, int i, int j, double dt,
		double dx, double dr);

double lagrange_m(gasCell * cell);

void lagrange_mass(double array[21], cell2d& cell, int i, int j, int n,
		double dx, double dr, double dt);

double lagrange_e(gasCell * prevCell, gasCell * cell);

double final_calc_Vx(cell2d& cell, int i, int j, int n, double dx, double dr,
		double dt);

double final_calc_Vr(cell2d& cell, int i, int j, int n, double dx, double dr,
		double dt);

double final_calc_e(cell2d& cell, int i, int j, int n, double dx, double dr,
		double dt);

double final_calc_z(cell2d * previousCell, gasCell * cell, int n, int i, int j);

double final_calc_psi(cell2d * previousCell, gasCell * cell, int n, int i,
		int j);

double final_calc_p(gasCell * prevCell, gasCell * cell, int var);

double euler_z(cell2d * previousCell, gasCell * cell, int n, int i, int j);

double euler_psi(cell2d * previousCell, gasCell * cell, int n, int i, int j);

double new_final_z(cell2d& cell, int i, int j, int n, double dx, double dr,
		double dt);

double new_final_psi(cell2d& cell, int i, int j, int n, double dx, double dr,
		double dt);

double smooth_Vr(cell2d * cell);

#endif	/* FUNCTIONS_H */
