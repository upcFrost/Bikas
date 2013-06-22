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
			P_sn += cell.at(borderI-1).at(j).P;
		} else {
			if (bottom_j == 0) bottom_j++;
		}
	}
	bottom_j++;

	P_sn /= count;

//
//	gasCell * curCell = &cell.at(i_sn-1).at(10);
//	// Local speed of sound
//	double ai = sqrt(k * curCell->P / curCell->rho);
//	// Pressure at the border
//	P_sn = curCell->P + ai*curCell->rho *
//		(U_sn.back() - curCell->Vx);
}

void projSpeedPositionCalc(cell2dStatic & cell, double P_sn,
		int top_j, int bottom_j, bool PROJECTILE, int & borderI) {

	double S = M_PI*pow(top_j*dr,2) - M_PI*pow(bottom_j*dr,2);

	if (PROJECTILE) {

		U_sn.push_back(euler_Usn(P_sn, S,
			0, dt, U_sn.back(), m_sn));
		x_sn.push_back(euler_Xsn(x_sn.back(), U_sn.back()));
		i_sn = floor(x_sn.back() / dx);
		borderI = i_sn;

	} else {

		double P_pist = 0; int count = 0;
		for (int j = 0; j < max_j; j++) {
			if (cell.at(i_pist+1).at(j).type != 18) {
				count++;
				P_pist += cell.at(i_pist+1).at(j).P;
			}
		}
		P_pist /= count;

		U_pist.push_back(euler_Usn(P_sn, S,
			P_pist*S, dt, U_pist.back(), 0.01));
		if (U_pist.back() < 0) {
			printf("123");
		}
		x_pist.push_back(euler_Xsn(x_pist.back(), U_pist.back()));
		i_pist = floor(x_pist.back() / dx);
		borderI = i_pist;

	}
}

void projCheckIfChanged(cell2d & fullCell,
		cell2dStatic & cell, cell2dStatic & nextTCell,
		int borderI_prev, int borderI) {
	if (borderI_prev != borderI && borderI > borderI_prev) {
//		/* Return cell to its original shape */
//		for (int j = 0; j < max_j; j++) {
//			double tempArray[5];
//			pre_cell_geometry(tempArray, cell.at(borderI_prev-1).at(j), borderI_prev-1, j);
//			for (int n = 0; n < maxN; n++) {
//				fullCell[n][borderI_prev-1][j].A[0] = tempArray[0];
//				fullCell[n][borderI_prev-1][j].A[1] = tempArray[1];
//				fullCell[n][borderI_prev-1][j].A[2] = tempArray[2];
//				fullCell[n][borderI_prev-1][j].A[3] = tempArray[3];
//				fullCell[n][borderI_prev-1][j].A[4] = tempArray[4];
//			}
//		}
		for (int j = 0; j < max_j; j++) {
			gasCell * oldCell = &cell.at(borderI_prev-1).at(j);
			gasCell * newCell = &cell.at(borderI-1).at(j);

			// Copy values to the new cell
			newCell->P = oldCell->P;
			newCell->rho = oldCell->rho;
			newCell->e = oldCell->e;
			newCell->bar_Vx = oldCell->bar_Vx;
			newCell->bar_Vr = oldCell->bar_Vr;
			newCell->Vx = oldCell->Vx;
			newCell->Vr = oldCell->Vr;
			newCell->final_z = oldCell->final_z;
			newCell->final_psi = oldCell->final_psi;

			if (borderI == i_pist) {
				gasCell * empty = &cell.at(borderI).at(j);
				gasCell * nextTEmpty = &nextTCell.at(borderI).at(j);
				empty->P = 0;
				empty->rho = 0;
				empty->e = 0;
				empty->bar_Vx = 0;
				empty->bar_Vr = 0;
				empty->Vx = 0;
				empty->Vr = 0;
				empty->final_z = 0;
				empty->final_psi = 0;
				nextTEmpty->P = 0;
				nextTEmpty->rho = 0;
				nextTEmpty->e = 0;
				nextTEmpty->bar_Vx = 0;
				nextTEmpty->bar_Vr = 0;
				nextTEmpty->Vx = 0;
				nextTEmpty->Vr = 0;
				nextTEmpty->final_z = 0;
				nextTEmpty->final_psi = 0;
			}
		}
	}
}

void projBorderMove(cell2dStatic & cell, int borderI, bool PROJECTILE) {
	double arrayT[5] = {0};
	for (int j = 0; j < max_j; j++) {
		if (cell.at(borderI-1).at(j).type != 18) {
			if (PROJECTILE)
				euler_proj_broder(arrayT, j, x_sn.back(), dx, dr, PROJECTILE);
			else
				euler_proj_broder(arrayT, j, x_pist.back(), dx, dr, PROJECTILE);
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

void projParCalc(cell2d & cell, int borderI_prev, int borderI, int var,
		bool PROJECTILE, bool debug) {
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
							cell.at(prevN).at(borderI-1).at(j).A[0]);
				}
				barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(prevN).at(borderI-1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			} else if (borderI > borderI_prev) {
				if (debug) {
					if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
					if (j == 2) printf("Prev A[0] = %16.16f\n",
							cell.at(prevN).at(borderI-2).at(j).A[0]);
				}
				barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(prevN).at(borderI-2).at(j).A[0]-1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			} else {
				if (debug) {
					if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
					if (j == 2) printf("Prev A[0] = %16.16f\n",
							cell.at(prevN).at(borderI-2).at(j).A[0]);
				}
				barQi = (curCell->A[0]-1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
				Qi = (cell.at(prevN).at(borderI-2).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			}

			// Local speed of sound
			double ai = soundSpeed(curCell->P, curCell->rho, curCell->final_psi, var);
			// Pressure at the border
			double borderP;
			if (PROJECTILE)
				borderP = curCell->P + ai*curCell->rho *
					(U_sn.back() - curCell->Vx);
			else
				borderP = curCell->P + ai*curCell->rho *
					(U_pist.back() - curCell->Vx);
			// Density at the center of the cell
			double newRho = curCell->rho * Qi / barQi;
			// Gas velocity at the center of the cell
			double newVx = curCell->rho / newRho * Qi / barQi * curCell->Vx +
					(borderP - curCell->P) / newRho / barQi *
					dt * M_PI*(2*(j-axis_j)+1)*pow(dr,2);
			// Gas full energy at the center of the cell
			double newE = curCell->e + borderP * (Qi - barQi) / curCell->rho / Qi;
			// Gas pressure at the center of the cell
			switch (var) {
			case IDEAL_GAS:
				newP = (k-1) * (newE - (pow(newVx,2) + pow(curCell->Vr,2))/2) * newRho;
				break;

			case POWDER_EQ:
				newP = (k-1) * (newE - (pow(newVx,2) + pow(curCell->Vr,2))/2 -
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
				if (j == 2) printf("delta P = %10.10f\n",curCell->P - newP);
				if (j == 2) debug_projectile_par(borderI, j, borderP, newRho,
						newVx, newE, newP, curCell->final_psi, x_sn.back());
			}

			if (j == 10 && (newP != newP || newE != newE || newVx != newVx)) {
				printf("123");
			}

			curCell->P = newP;
			curCell->e = newE;
			curCell->Vx = newVx;
			curCell->rho = newRho;
		}
	}
}

void projReturnSizeIfChanged(cell2d & fullCell,
		cell2dStatic & cell, cell2dStatic & nextTCell,
		int borderI_prev, int borderI) {
	if (borderI_prev != borderI && borderI > borderI_prev) {
	/* Return cell to its original shape */
		for (int j = 0; j < max_j; j++) {
			double tempArray[5];
			pre_cell_geometry(tempArray, cell.at(borderI_prev-1).at(j), borderI_prev-1, j);
			for (int n = 0; n < maxN; n++) {
				fullCell[n][borderI_prev-1][j].A[0] = tempArray[0];
				fullCell[n][borderI_prev-1][j].A[1] = tempArray[1];
				fullCell[n][borderI_prev-1][j].A[2] = tempArray[2];
				fullCell[n][borderI_prev-1][j].A[3] = tempArray[3];
				fullCell[n][borderI_prev-1][j].A[4] = tempArray[4];
			}
		}
	}
}

void projCalc(cell2d & cell, int var, int borderI,
		bool PROJECTILE, bool debug) {
	double P_sn = 0; int top_j = 0; int bottom_j = 0;
	int borderI_prev = borderI;
	projPCalc(cell.at(n), P_sn, top_j, bottom_j, borderI);
	projSpeedPositionCalc(cell.at(n), P_sn, top_j, bottom_j, PROJECTILE, borderI);
	projCheckIfChanged(cell, cell.at(n), cell.at(nextN), borderI_prev, borderI);
	projBorderMove(cell.at(n), borderI, PROJECTILE);
	projParCalc(cell, borderI_prev, borderI, var, PROJECTILE, debug);
	projReturnSizeIfChanged(cell, cell.at(n), cell.at(nextN), borderI_prev, borderI);
}

void pistonCalc(cell2d & cell, int borderI_prev, int borderI,
		int var, bool debug) {
	double arrayT[5] = {0};

	for (int j = 0; j < max_j; j++) {
		if (cell.at(n).at(borderI+1).at(j).type != 18) {
			euler_pist_broder(arrayT, j, x_pist.back(), dx, dr);
			// For n
			for (int iter = 0; iter < 5; iter++) {
				cell.at(n).at(borderI+1).at(j).A[iter] = arrayT[iter];
			}
			if (cell.at(n).at(borderI+1).at(j+1).type == 18) cell.at(n).at(borderI+1).at(j).A[4] = 0;
			if (cell.at(n).at(borderI+1).at(j-1).type == 18) cell.at(n).at(borderI+1).at(j).A[3] = 0;

			for (int iter = 0; iter < 5; iter++) {
				if (cell.at(n).at(borderI+1).at(j).A[iter] >= 2) {
					cell.at(n).at(borderI+1).at(j).A[iter] -= 1;
				}
			}
		}
	}

	if (borderI_prev != borderI) {
		for (int j = 0; j < max_j; j++) {
			gasCell * oldCell = &cell.at(n).at(borderI_prev+1).at(j);
			gasCell * newCell = &cell.at(n).at(borderI+1).at(j);

			newCell->P = oldCell->P;
			newCell->rho = oldCell->rho;
			newCell->e = oldCell->e;
			newCell->bar_Vx = oldCell->bar_Vx;
			newCell->bar_Vr = oldCell->bar_Vr;
			newCell->Vx = oldCell->Vx;
			newCell->Vr = oldCell->Vr;
			newCell->final_z = oldCell->final_z;
			newCell->final_psi = oldCell->final_psi;

//			if (borderI == i_pist) {
//				gasCell * empty = &cell.at(borderI).at(j);
//				gasCell * nextTEmpty = cell.at(nextN).at(borderI).at(j);
//				empty->P = 0;
//				empty->P[1] = 0;
//				empty->P[2] = 0;
//				empty->P[3] = 0;
//				empty->P[4] = 0;
//				empty->rho = 0;
//				empty->e = 0;
//				empty->bar_Vx = 0;
//				empty->bar_Vx[1] = 0;
//				empty->bar_Vx[2] = 0;
//				empty->bar_Vx[3] = 0;
//				empty->bar_Vx[4] = 0;
//				empty->bar_Vr = 0;
//				empty->bar_Vr[1] = 0;
//				empty->bar_Vr[2] = 0;
//				empty->bar_Vr[3] = 0;
//				empty->bar_Vr[4] = 0;
//				empty->Vx = 0;
//				empty->Vx[1] = 0;
//				empty->Vx[2] = 0;
//				empty->Vx[3] = 0;
//				empty->Vx[4] = 0;
//				empty->Vr = 0;
//				empty->Vr[1] = 0;
//				empty->Vr[2] = 0;
//				empty->Vr[3] = 0;
//				empty->Vr[4] = 0;
//				empty->final_z = 0;
//				empty->final_psi = 0;
//				nextTEmpty->P = 0;
//				nextTEmpty->P[1] = 0;
//				nextTEmpty->P[2] = 0;
//				nextTEmpty->P[3] = 0;
//				nextTEmpty->P[4] = 0;
//				nextTEmpty->rho = 0;
//				nextTEmpty->e = 0;
//				nextTEmpty->bar_Vx = 0;
//				nextTEmpty->bar_Vr = 0;
//				nextTEmpty->Vx = 0;
//				nextTEmpty->Vx[1] = 0;
//				nextTEmpty->Vx[2] = 0;
//				nextTEmpty->Vx[3] = 0;
//				nextTEmpty->Vx[4] = 0;
//				nextTEmpty->Vr = 0;
//				nextTEmpty->Vr[1] = 0;
//				nextTEmpty->Vr[2] = 0;
//				nextTEmpty->Vr[3] = 0;
//				nextTEmpty->Vr[4] = 0;
//				nextTEmpty->final_z = 0;
//				nextTEmpty->final_psi = 0;
//			}
		}
	}

	for (int j = 0; j < max_j; j++) {
		if (cell.at(n).at(borderI+1).at(j).type != 18) {
			gasCell * curCell = &cell.at(n).at(borderI+1).at(j);
			double barQi;
			double Qi;
			double newP = 0;

			if (debug) {
				if (j == 2) printf("A[0] = %16.16f\n", curCell->A[0]);
				if (j == 2) printf("Prev A[0] = %16.16f\n",
						cell.at(prevN).at(borderI+1).at(j).A[0]);
			}
			barQi = (curCell->A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
			Qi = (cell.at(prevN).at(borderI+1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);

			if (barQi - Qi > 0.1*Qi) Qi = (cell.at(prevN).at(borderI+1).at(j).A[0]+1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
			if (Qi - barQi > 0.1*Qi) barQi = (cell.at(prevN).at(borderI+1).at(j).A[0]+1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);

			// Local speed of sound
			double ai = sqrt(k * curCell->P / curCell->rho);
			// Pressure at the border
			double borderP = curCell->P + ai*curCell->rho *
					(U_pist.back() - curCell->Vx);
			// Density at the center of the cell
			double newRho = curCell->rho * Qi / barQi;
			// Gas velocity at the center of the cell
			double newVx = curCell->rho / newRho * Qi / barQi * curCell->Vx +
					(borderP - curCell->P) / newRho / barQi *
					dt * M_PI*(2*(j-axis_j)+1)*pow(dr,2);
			// Gas full energy at the center of the cell
			double newE = curCell->e + borderP * (Qi - barQi) / curCell->rho / Qi;
			// Gas pressure at the center of the cell
			switch (var) {
				case IDEAL_GAS:
					newP = (k-1) * (newE - (pow(newVx,2) + pow(curCell->Vr,2))/2) * newRho;
					break;

				case POWDER_EQ:
					newP = (k-1) * (newE - (pow(newVx,2) + pow(curCell->Vr,2))/2 -
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

			if (j == 10) {
				std::cout << " A[0] = " << curCell->A[0] << ", borderI = " << borderI
						<< ", borderI_prev = " << borderI_prev << ", U = " << U_pist.back() << std::endl;
			}

			if (j == 10 && (newP > 35 || newP != newP || newE != newE || newVx != newVx)) {
				printf("123");
			}

			curCell->P = newP;
			curCell->e = newE;
			curCell->Vx = newVx;
			curCell->rho = newRho;
		}
	}
}
