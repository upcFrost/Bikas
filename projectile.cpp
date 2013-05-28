/*
 * projectile.cpp
 *
 *  Created on: May 22, 2013
 *      Author: frost
 */

#include "projectile.h"

void projPCalc(cell2dStatic & cell, double & P_sn,
		int & top_j, int & bottom_j) {
	double count = 0;

	for (int j = 0; j < max_j; j++) {
		if (cell.at(i_sn-1).at(j).type != 18) {
			count++;
			top_j = j;
			P_sn += cell.at(i_sn-1).at(j).P[0];
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
		int top_j, int bottom_j) {
	U_sn.push_back(euler_Usn(P_sn, M_PI*pow(top_j*dr,2) - M_PI*pow(bottom_j*dr,2),
			0, dt, U_sn.back()));
	x_sn.push_back(euler_Xsn(x_sn.back(), U_sn.back()));
	i_sn = floor(x_sn.back() / dx);
}

void projCheckIfChanged(cell2dStatic & cell, cell2dStatic & nextTCell,
		int i_sn_prev) {
	if (i_sn_prev != i_sn) {
		for (int j = 0; j < max_j; j++) {
			/* Return cell to its original shape */
			// For n
			double tempArray[5];
			pre_cell_geometry(tempArray, cell.at(i_sn_prev-1).at(j), i_sn_prev-1, j);
			gasCell * oldCell = &cell.at(i_sn_prev-1).at(j);
			gasCell * oldTCell = &nextTCell.at(i_sn_prev-1).at(j);
			gasCell * newCell = &cell.at(i_sn-1).at(j);

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

			/* Weights */
//			double weightP = 0.4999;
//			double weightRho = 0.4999;
//			double weightE = 0.4999;
//			double weightV = 0.6;
//
//			newCell->P[0] = 2*weightP*oldCell->P[0];
//			newCell->P[1] = 2*weightP*oldCell->P[1];
//			newCell->P[2] = 2*weightP*oldCell->P[2];
//			newCell->P[3] = 2*weightP*oldCell->P[3];
//			newCell->P[4] = 2*weightP*oldCell->P[4];
//			newCell->rho = 2*weightRho*oldCell->rho;
//			newCell->e = 2*weightE*oldCell->e;
//			newCell->Vx[0] = 2*weightV*oldCell->Vx[0];
//			newCell->Vx[1] = 2*weightV*oldCell->Vx[1];
//			newCell->Vx[2] = 2*weightV*oldCell->Vx[2];
//			newCell->Vx[3] = 2*weightV*oldCell->Vx[3];
//			newCell->Vx[4] = 2*weightV*oldCell->Vx[4];
//			newCell->Vr[0] = 2*weightV*oldCell->Vr[0];
//			newCell->Vr[1] = 2*weightV*oldCell->Vr[1];
//			newCell->Vr[2] = 2*weightV*oldCell->Vr[2];
//			newCell->Vr[3] = 2*weightV*oldCell->Vr[3];
//			newCell->Vr[4] = 2*weightV*oldCell->Vr[4];
//			newCell->final_z = 2*weightP*oldCell->final_z;
//			newCell->final_psi = 2*weightP*oldCell->final_psi;
//
//			oldCell->P[0] = 2*(1 - weightP)*oldCell->P[0];
//			oldCell->P[1] = 2*(1 - weightP)*oldCell->P[1];
//			oldCell->P[2] = 2*(1 - weightP)*oldCell->P[2];
//			oldCell->P[3] = 2*(1 - weightP)*oldCell->P[3];
//			oldCell->P[4] = 2*(1 - weightP)*oldCell->P[4];
//			oldCell->rho = 2*(1 - weightRho)*oldCell->rho;
//			oldCell->e = 2*(1 - weightE)*oldCell->e;
//			oldCell->Vx[0] = 2*(1 - weightV)*oldCell->Vx[0];
//			oldCell->Vx[1] = 2*(1 - weightV)*oldCell->Vx[1];
//			oldCell->Vx[2] = 2*(1 - weightV)*oldCell->Vx[2];
//			oldCell->Vx[3] = 2*(1 - weightV)*oldCell->Vx[3];
//			oldCell->Vx[4] = 2*(1 - weightV)*oldCell->Vx[4];
//			oldCell->Vr[0] = 2*(1 - weightV)*oldCell->Vr[0];
//			oldCell->Vr[1] = 2*(1 - weightV)*oldCell->Vr[1];
//			oldCell->Vr[2] = 2*(1 - weightV)*oldCell->Vr[2];
//			oldCell->Vr[3] = 2*(1 - weightV)*oldCell->Vr[3];
//			oldCell->Vr[4] = 2*(1 - weightV)*oldCell->Vr[4];
//			oldCell->final_z = 2*(1 - weightP)*oldCell->final_z;
//			oldCell->final_psi = 2*(1 - weightP)*oldCell->final_psi;
		}
	}
}

void projBorderMove(cell2dStatic & cell) {
	double arrayT[5] = {0};
	for (int j = 0; j < max_j; j++) {
		if (cell.at(i_sn-1).at(j).type != 18) {
			euler_proj_broder(arrayT, j, x_sn.back(), dx, dr);
			// For n
			for (int iter = 0; iter < 5; iter++) {
				cell.at(i_sn-1).at(j).A[iter] = arrayT[iter];
			}
			if (cell.at(i_sn-1).at(j+1).type == 18) cell.at(i_sn-1).at(j).A[4] = 0;
			if (cell.at(i_sn-1).at(j-1).type == 18) cell.at(i_sn-1).at(j).A[3] = 0;

			for (int iter = 0; iter < 5; iter++) {
				if (cell.at(i_sn-1).at(j).A[iter] >= 2) {
					cell.at(i_sn-1).at(j).A[iter] -= 1;
				}
			}
		}
	}
}

void projParCalc(cell2d & cell, int i_sn_prev, int var, bool debug) {
	for (int j = 0; j < max_j; j++) {
		if (cell.at(n).at(i_sn-1).at(j).type != 18) {
			gasCell * curCell = &cell.at(n).at(i_sn-1).at(j);
			double barQi;
			double Qi;
			double newP = 0;

			if (i_sn_prev == i_sn) {
				if (debug) {
					if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
					if (j == 2) printf("Prev A[0] = %16.16f\n",
							cell.at(n-1).at(i_sn-1).at(j).A[0]);
				}
				barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(n-1).at(i_sn-1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			} else {
				if (debug) {
					if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
					if (j == 2) printf("Prev A[0] = %16.16f\n",
							cell.at(n-1).at(i_sn-2).at(j).A[0]);
				}
				barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(n-1).at(i_sn-2).at(j).A[0]-1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
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
				if (j == 2) debug_projectile_par(i_sn, j, borderP, newRho,
						newVx, newE, newP, curCell->final_psi, x_sn.back());
			}

			curCell->P[0] = newP;
			curCell->e = newE;
			curCell->Vx[0] = newVx;
			curCell->rho = newRho;
		}
	}
}

void projCalc(cell2d & cell, int var, bool debug) {
	double P_sn = 0; int top_j = 0; int bottom_j = 0;
	int i_sn_prev = i_sn;
	projPCalc(cell.at(n), P_sn, top_j, bottom_j);
	projSpeedPositionCalc(cell.at(n), P_sn, top_j, bottom_j);
	projCheckIfChanged(cell.at(n), cell.at(n+1), i_sn_prev);
	projBorderMove(cell.at(n));
	projParCalc(cell, i_sn_prev, var, debug);
}

