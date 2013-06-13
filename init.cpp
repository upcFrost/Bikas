/*
 * init.cpp
 *
 *  Created on: May 20, 2013
 *      Author: frost
 */

#include "init.h"

/* Alpha calculation */
double pre_alpha(gasCell cell, int i, int i_p, double delta_0, double delta) {
    if (i <= i_p)
		//return delta_0/delta;
		return 1;
	else
		return 1;
}

void getGlobalVars(std::ifstream & inputFile, bool havePiston) {
	std::string line;
	std::string tmpLine;

	/* Drop header line */
	std::getline (inputFile, line);

	/* Get initial vars from input file */
	int charNum = 0;
	getline (inputFile, line);
	std::stringstream varsLine(line);
	while(std::getline (varsLine, tmpLine, ',')) {
		switch (charNum) {
		case 0:
			max_i = atoi(tmpLine.c_str());
			break;
		case 1:
			max_j = atoi(tmpLine.c_str());
			break;
		case 2:
			delta = atof(tmpLine.c_str());
			break;
		case 3:
			delta_0 = atof(tmpLine.c_str());
			break;
		case 4:
			dr = atof(tmpLine.c_str());
			break;
		case 5:
			dx = atof(tmpLine.c_str());
			break;
		case 6:
			p0 = atof(tmpLine.c_str());
			break;
		case 7:
			kappa = atof(tmpLine.c_str());
			break;
		case 8:
			lambda = atof(tmpLine.c_str());
			break;
		case 9:
			I_k = atoi(tmpLine.c_str());
			break;
		case 10:
			//P_v = atof(tmpLine.c_str());
			break;
		case 11:
			i_v = atoi(tmpLine.c_str());
			break;
		case 12:
			j_v = atoi(tmpLine.c_str());
			break;
		case 13:
			x_sn.at(0) = atof(tmpLine.c_str());
			break;
		case 14:
			if (havePiston) {
				if (x_pist.size() == 0)
					x_pist.push_back(0);
				x_pist.at(0) = atof(tmpLine.c_str());
			}
			break;
		default:
			break;
		}
		charNum++;
	}
}

void scaleGlobalVars(bool havePiston) {
	x_sn.at(0) = x_sn.at(0) / scaleD;
	if (havePiston) x_pist.at(0) = x_pist.at(0) / scaleD;
	m_sn = m_sn / scaleM;
	P_f = P_f / scaleP;
	P_atm = P_atm / scaleP;
	P_v = P_v / scaleP;
	dx = dx / scaleD;
	dr = dr / scaleD;
	delta_0 = delta_0 / scaleRho;
	delta = delta / scaleRho;
	dt = dt / scaleT;
	I_k = I_k / scaleIK;
	f = f / scaleFF;
}

void initCellVector(cell2d & cell) {
	n = 0;
	cell.resize(1);
	cell[n].resize(max_i);

	for (int i = 0; i < max_i; i++)
	{
	    cell[n][i].resize(max_j);
	}
}

gasCell setEmptyCell(gasCell & cell) {
	gasCell empty = cell;
	empty.P[0] = 0;
	empty.P[1] = 0;
	empty.P[2] = 0;
	empty.P[3] = 0;
	empty.P[4] = 0;
	empty.rho = 0;
	empty.e = 0;
	empty.bar_Vx[0] = 0;
	empty.bar_Vx[1] = 0;
	empty.bar_Vx[2] = 0;
	empty.bar_Vx[3] = 0;
	empty.bar_Vx[4] = 0;
	empty.bar_Vr[0] = 0;
	empty.bar_Vr[1] = 0;
	empty.bar_Vr[2] = 0;
	empty.bar_Vr[3] = 0;
	empty.bar_Vr[4] = 0;
	empty.Vx[0] = 0;
	empty.Vx[1] = 0;
	empty.Vx[2] = 0;
	empty.Vx[3] = 0;
	empty.Vx[4] = 0;
	empty.Vr[0] = 0;
	empty.Vr[1] = 0;
	empty.Vr[2] = 0;
	empty.Vr[3] = 0;
	empty.Vr[4] = 0;
	empty.final_z = 0;
	empty.final_psi = 0;
	return empty;
}

void populateCellVector(std::ifstream & inputFile, cell2d & cell,
		int var, bool havePiston, bool debug) {
	std::string line;
	std::string tmpLine;
	int charNum = 0;
	int i = 0;
	int j = 0;

	/* Skip header line */
	std::getline(inputFile, line);
	while (inputFile.good()) {

		/* Set types, alpha, A and F from input file */
		getline(inputFile, line);
		charNum = 0;
		std::stringstream varsLine(line);
		while(getline (varsLine, tmpLine, ',')) {
			switch (charNum) {
				case 0:
					i = atoi(tmpLine.c_str());
					break;
				case 1:
					j = atoi(tmpLine.c_str());
					break;
				case 2:
					cell.at(n).at(i).at(j).type = atoi(tmpLine.c_str());
					break;
				case 3:
					cell.at(n).at(i).at(j).r_1 = atof(tmpLine.c_str());
					break;
				case 4:
					cell.at(n).at(i).at(j).r_2 = atof(tmpLine.c_str());
					break;
				case 5:
					cell.at(n).at(i).at(j).x_1 = atof(tmpLine.c_str());
					break;
				case 6:
					cell.at(n).at(i).at(j).x_2 = atof(tmpLine.c_str());
					break;
				default:
					break;
			}
			charNum++;
		}

		gasCell * curCell = &cell.at(n).at(i).at(j);

		curCell->alpha = pre_alpha(cell.at(n).at(i).at(j), i, max_i, delta_0, delta);
		double A[5];
		pre_cell_geometry(A, cell.at(n).at(i).at(j), i, j);
		if (curCell->type == 1 || curCell->type == 3 ||	curCell->type == 8 ||
				curCell->type == 10 || curCell->type == 22 )	{
			curCell->weightVector = wightVectorsCalc(cell, i, j, n, debug);
		}
		for (unsigned int iter = 0; iter < 5; iter++)
			curCell->A[iter] = A[iter];
		i_sn = floor(x_sn.at(0)/dx);
		i_sn_0 = i_sn;
		if (havePiston) {
			if (fmod(x_pist.back(),dx) < 0.75*dx) {
				i_pist = floor(x_pist.back() / dx);
				mergedI = false;
			} else {
				i_pist = floor(x_pist.back() / dx) + 1;
				mergedI = true;
			}
			i_pist_0 = i_pist;
			if (U_pist.size() == 0)
				U_pist.push_back(0);
		}
		max_z_0 = max_z;

		/* Set gas init parameters */
		if ((i < i_pist && havePiston) || (i < i_sn && !havePiston)) {
			for (int iter = 0; iter < 5; iter++) curCell->P[iter] = P_v;
			for (int iter = 0; iter < 5; iter++) curCell->Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->bar_Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->bar_Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->dM[iter] = 0;
			curCell->rho = delta_0;
			switch (var) {
			case IDEAL_GAS:
				curCell->final_psi = 1;
				curCell->final_z = 1;
				curCell->e = cell.at(n).at(i).at(j).P[0] / (k-1) / delta_0;
				break;
			case ABEL_DUPRE:
				curCell->final_psi = 1;
				curCell->final_z = 1;
				curCell->e = cell.at(n).at(i).at(j).P[0] / (k-1) / delta_0;
				break;
			case POWDER_EQ:
				curCell->final_psi = (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1);
				curCell->final_z = 2 * curCell->final_psi / (kappa * (1 + sqrt(1 + 4*lambda*curCell->final_psi/kappa)));
				curCell->e = curCell->P[0] / (k-1) * (1/curCell->rho -
					(1 - curCell->final_psi)/delta - alpha_k * curCell->final_psi) +
					f/(k-1)*(1-curCell->final_psi);
				break;
			default:
				break;
			}
		} else if (i < i_sn && havePiston) {
//			curCell->rho = (PISTON_RHO / scaleRho);
//			curCell->rho = 931.779;
			curCell->rho = 300;
//			for (int iter = 0; iter < 5; iter++) curCell->P[iter] = PISTON_B *
//					curCell->rho/PISTON_RHO * (curCell->rho/PISTON_RHO - 1) /
//					pow(PISTON_C - curCell->rho/PISTON_RHO, 2);
			for (int iter = 0; iter < 5; iter++) curCell->P[iter] = P_v;
			for (int iter = 0; iter < 5; iter++) curCell->Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->bar_Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->bar_Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->dM[iter] = 0;
			curCell->e = cell.at(n).at(i).at(j).P[0] / (k-1) / curCell->rho;
//			curCell->e = cell.at(n).at(i).at(j).P[0] / (k-1) / (PISTON_RHO/ scaleRho);
			curCell->final_psi = 1;
			curCell->final_z = 1;
		} else {
			for (int iter = 0; iter < 5; iter++) curCell->P[iter] = P_atm;
			for (int iter = 0; iter < 5; iter++) curCell->Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->dM[iter] = 0;
			curCell->rho = rho_atm; // Air
			curCell->e = P_atm / rho_atm / (k - 1); // Ideal gas
			curCell->final_psi = 1;
			curCell->final_z = 1;
		}

		if (i == i_pist) {
			cell.at(n).at(i_pist).at(j) = setEmptyCell(cell.at(n).at(i_pist).at(j));
		}
	}
}

void init(std::ifstream & inputFile, cell2d & cell, int var,
		bool havePiston, bool debug) {
	getGlobalVars(inputFile, havePiston);
	U_sn.at(0) = 0;
	scaleGlobalVars(havePiston);
	initCellVector(cell);
	populateCellVector(inputFile, cell, var, havePiston, debug);
}

void borderCellsFix(cell2d & cell, bool havePiston) {
	 for (int j = 2; j < max_j-2; j++) {
		double arrayT[5] = {0};
		euler_proj_broder(arrayT, j, x_sn.back(), dx, dr, true);
		// For n
		cell.at(n).at(i_sn-1).at(j).A[0] = arrayT[0];
		cell.at(n).at(i_sn-1).at(j).A[1] = arrayT[1];
		cell.at(n).at(i_sn-1).at(j).A[2] = arrayT[2];
		cell.at(n).at(i_sn-1).at(j).A[3] = arrayT[3];
		cell.at(n).at(i_sn-1).at(j).A[4] = arrayT[4];
		if (j == max_j - 2) cell.at(n).at(i_sn-1).at(j).A[4] = 0;
		if (j == 2) cell.at(n).at(i_sn-1).at(j).A[3] = 0;
	}
	if (havePiston) {
		for (int j = 2; j < max_j-2; j++) {
			double arrayT[5] = {0};
			euler_proj_broder(arrayT, j, x_pist.back(), dx, dr, true);
			// For n
			cell.at(n).at(i_pist-1).at(j).A[0] = arrayT[0];
			cell.at(n).at(i_pist-1).at(j).A[1] = arrayT[1];
			cell.at(n).at(i_pist-1).at(j).A[2] = arrayT[2];
			cell.at(n).at(i_pist-1).at(j).A[3] = arrayT[3];
			cell.at(n).at(i_pist-1).at(j).A[4] = arrayT[4];
			if (cell.at(n).at(i_pist-1).at(j+1).type == 18) cell.at(n).at(i_pist-1).at(j).A[4] = 0;
			if (cell.at(n).at(i_pist-1).at(j-1).type == 18) cell.at(n).at(i_pist-1).at(j).A[3] = 0;
		}
		for (int j = 2; j < max_j-2; j++) {
			double arrayT[5] = {0};
			euler_pist_broder(arrayT, j, x_pist.back(), dx, dr);
			// For n
			cell.at(n).at(i_pist+1).at(j).A[0] = arrayT[0];
			cell.at(n).at(i_pist+1).at(j).A[1] = arrayT[1];
			cell.at(n).at(i_pist+1).at(j).A[2] = arrayT[2];
			cell.at(n).at(i_pist+1).at(j).A[3] = arrayT[3];
			cell.at(n).at(i_pist+1).at(j).A[4] = arrayT[4];
			if (cell.at(n).at(i_pist+1).at(j+1).type == 18) cell.at(n).at(i_pist+1).at(j).A[4] = 0;
			if (cell.at(n).at(i_pist+1).at(j-1).type == 18) cell.at(n).at(i_pist+1).at(j).A[3] = 0;
		}
	}
}
