/*
 * projectile.cpp
 *
 *  Created on: May 22, 2013
 *      Author: frost
 */

#include "projectile.h"

void projPCalc(cell2dStatic & cell, double & P_sn,
		int & top_j, int & bottom_j, int borderI) {
	double count = 0;

	for (int j = 0; j < max_j; j++) {
		if (cell.at(borderI-1).at(j).type != 18) {
			count++;
			top_j = j;
			P_sn += cell.at(borderI-1).at(j).P[0];
		} else {
			if (bottom_j == 0) bottom_j++;
		}
	}
	bottom_j++;

	P_sn /= count;

//
//	gasCell * curCell = &cell.at(i_sn-1).at(10);
//	// Local speed of sound
//	double ai = sqrt(k * curCell->P[0] / curCell->rho);
//	// Pressure at the border
//	P_sn = curCell->P[0] + ai*curCell->rho *
//		(U_sn.back() - curCell->Vx[0]);
}

void projSpeedPositionCalc(cell2dStatic & cell, double P_sn,
		int top_j, int bottom_j, bool PROJECTILE) {
	if (PROJECTILE) {
		U_sn.push_back(euler_Usn(P_sn, M_PI*pow(top_j*dr,2) - M_PI*pow(bottom_j*dr,2),
			0, dt, U_sn.back()));
		x_sn.push_back(euler_Xsn(x_sn.back(), U_sn.back()));
		i_sn = floor(x_sn.back() / dx);
	} else {
		U_pist.push_back(euler_Usn(P_sn, M_PI*pow(top_j*dr,2) - M_PI*pow(bottom_j*dr,2),
			0, dt, U_pist.back()));
		x_pist.push_back(euler_Xsn(x_pist.back(), U_pist.back()));
		i_pist = floor(x_pist.back() / dx);
	}
}

void projCheckIfChanged(cell2dStatic & cell, cell2dStatic & nextTCell,
		int borderI_prev, int borderI) {
	if (borderI_prev != borderI) {
		for (int j = 0; j < max_j; j++) {
			/* Return cell to its original shape */
			// For n
			double tempArray[5];
			pre_cell_geometry(tempArray, cell.at(borderI_prev-1).at(j), borderI_prev-1, j);
			gasCell * oldCell = &cell.at(borderI_prev-1).at(j);
			gasCell * oldTCell = &nextTCell.at(borderI_prev-1).at(j);
			gasCell * newCell = &cell.at(borderI-1).at(j);

			oldCell->A[0] = tempArray[0];
			oldCell->A[1] = tempArray[1];
			oldCell->A[2] = tempArray[2];
			oldCell->A[3] = tempArray[3];
			oldCell->A[4] = tempArray[4];
			oldTCell->A[0] = tempArray[0];
			oldTCell->A[1] = tempArray[1];
			oldTCell->A[2] = tempArray[2];
			oldTCell->A[3] = tempArray[3];
			oldTCell->A[4] = tempArray[4];

			newCell->P[0] = oldCell->P[0];
			newCell->P[1] = oldCell->P[1];
			newCell->P[2] = oldCell->P[2];
			newCell->P[3] = oldCell->P[3];
			newCell->P[4] = oldCell->P[4];
			newCell->rho = oldCell->rho;
			newCell->e = oldCell->e;
			newCell->bar_Vx[0] = oldCell->bar_Vx[0];
			newCell->bar_Vx[1] = oldCell->bar_Vx[1];
			newCell->bar_Vx[2] = oldCell->bar_Vx[2];
			newCell->bar_Vx[3] = oldCell->bar_Vx[3];
			newCell->bar_Vx[4] = oldCell->bar_Vx[4];
			newCell->bar_Vr[0] = oldCell->bar_Vr[0];
			newCell->bar_Vr[1] = oldCell->bar_Vr[1];
			newCell->bar_Vr[2] = oldCell->bar_Vr[2];
			newCell->bar_Vr[3] = oldCell->bar_Vr[3];
			newCell->bar_Vr[4] = oldCell->bar_Vr[4];
			newCell->Vx[0] = oldCell->Vx[0];
			newCell->Vx[1] = oldCell->Vx[1];
			newCell->Vx[2] = oldCell->Vx[2];
			newCell->Vx[3] = oldCell->Vx[3];
			newCell->Vx[4] = oldCell->Vx[4];
			newCell->Vr[0] = oldCell->Vr[0];
			newCell->Vr[1] = oldCell->Vr[1];
			newCell->Vr[2] = oldCell->Vr[2];
			newCell->Vr[3] = oldCell->Vr[3];
			newCell->Vr[4] = oldCell->Vr[4];
			newCell->final_z = oldCell->final_z;
			newCell->final_psi = oldCell->final_psi;
		}
	}
}

void projBorderMove(cell2dStatic & cell, int borderI) {
	double arrayT[5] = {0};
	for (int j = 0; j < max_j; j++) {
		if (cell.at(borderI-1).at(j).type != 18) {
			euler_proj_broder(arrayT, j, x_sn.back(), dx, dr);
			// For n
			for (int iter = 0; iter < 5; iter++) {
				cell.at(borderI-1).at(j).A[iter] = arrayT[iter];
			}
			if (cell.at(borderI-1).at(j+1).type == 18) cell.at(borderI-1).at(j).A[4] = 0;
			if (cell.at(borderI-1).at(j-1).type == 18) cell.at(borderI-1).at(j).A[3] = 0;

			for (int iter = 0; iter < 5; iter++) {
				if (cell.at(borderI-1).at(j).A[iter] >= 2) {
					cell.at(borderI-1).at(j).A[iter] -= 1;
				}
			}
		}
	}
}

void projParCalc(cell2d & cell, int borderI_prev, int borderI, int var, bool debug) {
	for (int j = 0; j < max_j; j++) {
		if (cell.at(n).at(borderI-1).at(j).type != 18) {
			gasCell * curCell = &cell.at(n).at(borderI-1).at(j);
			double barQi;
			double Qi;
			double newP = 0;

			if (borderI_prev == borderI) {
				if (debug) {
					if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
					if (j == 2) printf("Prev A[0] = %16.16f\n",
							cell.at(n-1).at(borderI-1).at(j).A[0]);
				}
				barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(n-1).at(borderI-1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			} else {
				if (debug) {
					if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
					if (j == 2) printf("Prev A[0] = %16.16f\n",
							cell.at(n-1).at(borderI-2).at(j).A[0]);
				}
				barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(n-1).at(borderI-2).at(j).A[0]-1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			}

//			Qi = M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
//			barQi = (1 + curCell->A[0] - cell.at(n-1).at(i_sn-1).at(j).A[0]) * Qi;
//			if (i_sn_prev != i_sn)
//				barQi = (2 + curCell->A[0] - cell.at(n-1).at(i_sn-2).at(j).A[0]) * Qi;

			// Local speed of sound
			double ai = sqrt(k * curCell->P[0] / curCell->rho);
			// Pressure at the border
			double borderP = curCell->P[0] + ai*curCell->rho *
					(U_sn.back() - curCell->Vx[0]);
			// Density at the center of the cell
			double newRho = curCell->rho * Qi / barQi;
			// Gas velocity at the center of the cell
			double newVx = curCell->rho / newRho * Qi / barQi * curCell->Vx[0] +
					(borderP - curCell->P[0]) / newRho / barQi *
					dt * M_PI*(2*(j-axis_j)+1)*pow(dr,2);
			// Gas full energy at the center of the cell
			double newE = curCell->e + borderP * (Qi - barQi) / curCell->rho / Qi;
			// Gas pressure at the center of the cell
			switch (var) {
			case IDEAL_GAS:
				newP = (k-1) * (newE - (pow(newVx,2) + pow(curCell->Vr[0],2))/2) * newRho;
				break;

			case POWDER_EQ:
				newP = (k-1) * (newE - (pow(newVx,2) + pow(curCell->Vr[0],2))/2 -
						f/(k-1)*(1-curCell->final_psi)) /
					( 1/newRho - (1 - curCell->final_psi)/delta -
							alpha_k * curCell->final_psi);
				break;

			case PISTON:
				newP = PISTON_B * newRho/PISTON_RHO * (newRho/PISTON_RHO - 1) /
					pow(PISTON_C - newRho/PISTON_RHO, 2);
				break;

			default:
				break;
			}

			/** DEBUG **/
			if (debug) {
				if (j == 2) printf("Old rho = %8.8f\n", curCell->rho);
				if (j == 2) printf("Qi = %10.10f\n",Qi);
				if (j == 2) printf("barQi = %10.10f\n",barQi);
				if (j == 2) printf("delta rho = %10.10f\n",curCell->rho - newRho);
				if (j == 2) printf("delta P = %10.10f\n",curCell->P[0] - newP);
				if (j == 2) debug_projectile_par(borderI, j, borderP, newRho,
						newVx, newE, newP, curCell->final_psi, x_sn.back());
			}

			curCell->P[0] = newP;
			curCell->e = newE;
			curCell->Vx[0] = newVx;
			curCell->rho = newRho;
		}
	}
}

void projCalc(cell2d & cell, int var, int borderI,
		bool PROJECTILE, bool debug) {
	double P_sn = 0; int top_j = 0; int bottom_j = 0;
	int borderI_prev = borderI;
	projPCalc(cell.at(n), P_sn, top_j, bottom_j, borderI);
	projSpeedPositionCalc(cell.at(n), P_sn, top_j, bottom_j, PROJECTILE);
	projCheckIfChanged(cell.at(n), cell.at(n+1), borderI_prev, borderI);
	projBorderMove(cell.at(n), borderI);
	projParCalc(cell, borderI_prev, borderI, var, debug);
}

