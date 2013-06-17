#ifndef DEBUG_H
#define	DEBUG_H

#include "types.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>


void debug_projectile_par(int i_sn, int j, double borderP, 
	double newRho, double newVx, double newE, double newP, 
	double final_psi, double x_sn);
	
void debug_weights(int i, int j, double P,
	std::vector < WeightPart > weightVector);
	
void debug_equality_Vx_e(int i_sn, int max_j, int n, cell2d & cell);

void debug_p_output(int n, int nArray, int *i, int *j, cell2d & cell);

void debug_dM_rho_output(int n, int nArray, int *i, int *j, cell2d & cell);

void debug_Vx_Vr_P_A_barVx_output(int n, int nArray, int *i, int *j, cell2d & cell);

void debug_final_output(int n, int nArray, int *i, int *j, cell2d & cell);

#endif	/* DEBUG_H */
