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

void getGlobalVars(std::ifstream & inputFile) {
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
		default:
			break;
		}
		charNum++;
	}
}

void scaleGlobalVars() {
	x_sn.at(0) = x_sn.at(0) / scaleD;
	m_sn = m_sn / scaleM;
	P_f = P_f / scaleP;
	P_atm = P_atm / scaleP;
	P_v = P_v / scaleP;
	dx = dx / scaleD;
	dr = dr / scaleD;
	delta_0 = delta_0 / scaleR;
	delta = delta / scaleR;
	dt = dt / scaleT;
	I_k = I_k / scaleIK;
	f = f / scaleFF;
}

void initCellVector(cell2d & cell) {
	n = 0;
	cell.resize(1);
	cell[n].resize(max_i);

	for (i = 0; i < max_i; i++)
	{
	    cell[n][i].resize(max_j);
	}
}

void populateCellVector(std::ifstream & inputFile, cell2d & cell) {
	std::string line;
	std::string tmpLine;
	int charNum = 0;

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
			printf("Making weightVector for %d:%d\n",i,j);
			curCell->weightVector = wightVectorsCalc(cell, i, j, n, false);
			printf("weightVector size on %d for %d:%d = %u\n", n,i,j,
				(unsigned int)curCell->weightVector.y.size());
		}
		for (unsigned int iter = 0; iter < 5; iter++)
			curCell->A[iter] = A[iter];
		i_sn = floor(x_sn.at(0)/dx);
		i_sn_0 = i_sn;
		max_z_0 = max_z;

		/* Set gas init parameters */
		if (i < i_sn) {
			for (int iter = 0; iter < 5; iter++) curCell->P[iter] = P_v;
			for (int iter = 0; iter < 5; iter++) curCell->Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->bar_Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->bar_Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->dM[iter] = 0;
			curCell->rho = delta_0;
			curCell->final_psi = (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1);
			curCell->final_z = 2 * curCell->final_psi / (kappa * (1 + sqrt(1 + 4*lambda*curCell->final_psi/kappa)));
			//~ cell.at(n).at(i).at(j).e = cell.at(n).at(i).at(j).P[0] / (k-1) / delta_0; // Ideal gas
			//~ cell.at(n).at(i).at(j).e = cell.at(n).at(i).at(j).P[0] / (k-1) * (1/delta_0 - 1/delta); // BMSTU var
			curCell->e = curCell->P[0] / (k-1) * (1/curCell->rho - (1 - curCell->final_psi)/delta - alpha_k * curCell->final_psi); // From BMSTU P equation
			curCell->psi = (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1);
			curCell->z = 2 * curCell->psi / (kappa * (1 + sqrt(1 + 4*lambda*curCell->psi/kappa)));
		} else {
			for (int iter = 0; iter < 5; iter++) curCell->P[iter] = P_atm;
			for (int iter = 0; iter < 5; iter++) curCell->Vx[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->Vr[iter] = 0;
			for (int iter = 0; iter < 5; iter++) curCell->dM[iter] = 0;
			curCell->rho = rho_atm; // Air
			//~ cell.at(n).at(i).at(j).e = P_atm * 2*M_PI*(j*dr + dr/2)*dx*dr / 0.4; // Used for powder eq
			curCell->e = P_atm / rho_atm / (k - 1); // Ideal gas
			curCell->psi = 1;
			curCell->z = 1;
			curCell->final_psi = 1;
			curCell->final_z = 1;
		}
	}
}

void init(std::ifstream & inputFile, cell2d & cell) {
	getGlobalVars(inputFile);
	U_sn.at(0) = 0;
	scaleGlobalVars();
	initCellVector(cell);
	populateCellVector(inputFile, cell);
}
