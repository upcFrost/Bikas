#include "debug.h"

void debug_projectile_par(int i_sn, int j, double borderP, 
	double newRho, double newVx, double newE, double newP, 
	double final_psi, double x_sn)
{
	std::cout << std::endl << 
		"i_sn = " << i_sn << std::endl <<
		"x_sn = " << x_sn << std::endl <<
		"j = " << j << std::endl <<
		"borderP = " << borderP << std::endl <<
		"newRho = " << newRho << std::endl <<
		"newVx = " << newVx << std::endl <<
		"newE = " << newE << std::endl <<
		"newP = " << newP << std::endl <<
		"Final_psi = " << final_psi << std::endl << std::endl;
}

void debug_type1_cell(int i, int j, double sin_a, double cos_a, 
	double xbegin, double ybegin, double xend, double yend, 
	double r_1, double r_2,
	double center_d, double center_dx, double center_dy,
	Point2D vertices[4])
{
	std::cout << "i = " << i << std::endl <<
			"j = " << j << std::endl <<
			"cos_a = " << cos_a << std::endl <<
			"sin_a = " << sin_a << std::endl <<
			"r_1 = " << r_1 << std::endl <<
			"r_2 = " << r_2 << std::endl <<
			"xbegin = " << xbegin << std::endl <<
			"ybegin = " << ybegin << std::endl <<
			"xend = " << xend << std::endl <<
			"yend = " << yend << std::endl <<
			"point 0 = " << vertices[0].x << ":" << vertices[0].y << std::endl <<
			"point 1 = " << vertices[1].x << ":" << vertices[1].y << std::endl <<
			"point 2 = " << vertices[2].x << ":" << vertices[2].y << std::endl <<
			"point 3 = " << vertices[3].x << ":" << vertices[3].y << std::endl <<
			"center_d = " << center_d << std::endl <<
			"center_dx = " << center_dx << std::endl <<
			"center_dy = " << center_dy << std::endl;
	std::cout << std::endl;
}


void debug_weights(int i, int j, double P, std::vector < WeightPart > weightVector)
{
	std::cout << "i = " << i << ", j = " << j << ", P[0] = " << P << std::endl;
	std::cout << "Weights: " << std::endl;
	for (unsigned int idx2 = 0; idx2 < weightVector.size(); idx2++) {
		std::cout << "Weight: i = " << weightVector[idx2].i << ", j = " << weightVector[idx2].j << ", weight = " << weightVector[idx2].weight << std::endl;
	}
	std::cout << std::endl;
}


void debug_equality_Vx_e(int i_sn, int max_j, int n, cell2d cell)
{
	for (int i = 1; i < i_sn; i++) {
		for (int j = 1; j < max_j-3; j++) {
			if (cell.at(n).at(i).at(j).type != 18) {
				if (cell.at(n+1).at(i).at(j).e != cell.at(n+1).at(i).at(j+1).e) {
					std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(16) << "E is different" << std::endl <<
							"i    = " << i << std::endl <<
							"j    = " << j << std::endl <<
							"bar_E       = " << cell.at(n).at(i).at(j).bar_e << std::endl <<
							"lower bar_E = " << cell.at(n).at(i).at(j-1).bar_e << std::endl <<
							"upper bar_E = " << cell.at(n).at(i).at(j+1).bar_e << std::endl <<
							"E       = " << cell.at(n+1).at(i).at(j).e << std::endl <<
							"lower E = " << cell.at(n+1).at(i).at(j-1).e << std::endl <<
							"upper E = " << cell.at(n+1).at(i).at(j+1).e << std::endl <<
							"P       = " << cell.at(n+1).at(i).at(j).P[0] << std::endl <<
							"lower P = " << cell.at(n+1).at(i).at(j-1).P[0] << std::endl <<
							"upper P = " << cell.at(n+1).at(i).at(j+1).P[0] << std::endl <<
							"bar_Vx       = " << cell.at(n).at(i).at(j).bar_Vx[0] << std::endl <<
							"lower bar_Vx = " << cell.at(n).at(i).at(j-1).bar_Vx[0] << std::endl <<
							"upper bar_Vx = " << cell.at(n).at(i).at(j+1).bar_Vx[0] << std::endl <<
							"Vx       = " << cell.at(n+1).at(i).at(j).Vx[0] << std::endl <<
							"lower Vx = " << cell.at(n+1).at(i).at(j-1).Vx[0] << std::endl <<
							"upper Vx = " << cell.at(n+1).at(i).at(j+1).Vx[0] << std::endl <<
							"rho       = " << cell.at(n+1).at(i).at(j).rho << std::endl <<
							"lower rho = " << cell.at(n+1).at(i).at(j-1).rho << std::endl <<
							"upper rho = " << cell.at(n+1).at(i).at(j+1).rho << std::endl <<
							"prev rho       = " << cell.at(n).at(i).at(j).rho << std::endl <<
							"lower prev rho = " << cell.at(n).at(i).at(j-1).rho << std::endl <<
							"upper prev rho = " << cell.at(n).at(i).at(j+1).rho << std::endl <<
							"dM       = {" << cell.at(n).at(i).at(j).dM[1] << ", " << cell.at(n).at(i).at(j).dM[2] << ", " << cell.at(n).at(i).at(j).dM[3] << ", " << cell.at(n).at(i).at(j).dM[4] << "} " << std::endl <<
							"lower dM = {" << cell.at(n).at(i).at(j-1).dM[1] << ", " << cell.at(n).at(i).at(j-1).dM[2] << ", " << cell.at(n).at(i).at(j-1).dM[3] << ", " << cell.at(n).at(i).at(j-1).dM[4] << "} " << std::endl <<
							"upper dM = {" << cell.at(n).at(i).at(j+1).dM[1] << ", " << cell.at(n).at(i).at(j+1).dM[2] << ", " << cell.at(n).at(i).at(j+1).dM[3] << ", " << cell.at(n).at(i).at(j+1).dM[4] << "} " << std::endl << std::endl;
							for (int k = 1; k < 5; k++)
								if (cell.at(n).at(i).at(j).dM[k] != 0)
									std::cout << "dM[" << k << "] != 0" << std::endl;
							if (cell.at(n+1).at(i).at(j).rho != cell.at(n+1).at(i).at(j+1).rho)
								std::cout << "rho(j) != rho(j+1)" << std::endl;
							if (cell.at(n).at(i).at(j).bar_e != cell.at(n).at(i).at(j+1).bar_e)
								std::cout << "bar_e(j) != bar_e(j+1)" << std::endl;
					getchar();
				}
				if (cell.at(n+1).at(i).at(j).Vx[0] != cell.at(n+1).at(i).at(j+1).Vx[0]) {
					std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(16) << "Vx is different" << std::endl <<
							"i    = " << i << std::endl <<
							"j    = " << j << std::endl <<
							"bar_E       = " << cell.at(n).at(i).at(j).bar_e << std::endl <<
							"lower bar_E = " << cell.at(n).at(i).at(j-1).bar_e << std::endl <<
							"upper bar_E = " << cell.at(n).at(i).at(j+1).bar_e << std::endl <<
							"E       = " << cell.at(n+1).at(i).at(j).e << std::endl <<
							"lower E = " << cell.at(n+1).at(i).at(j-1).e << std::endl <<
							"upper E = " << cell.at(n+1).at(i).at(j+1).e << std::endl <<
							"P       = " << cell.at(n+1).at(i).at(j).P[0] << std::endl <<
							"lower P = " << cell.at(n+1).at(i).at(j-1).P[0] << std::endl <<
							"upper P = " << cell.at(n+1).at(i).at(j+1).P[0] << std::endl <<
							"bar_Vx       = " << cell.at(n).at(i).at(j).bar_Vx[0] << std::endl <<
							"lower bar_Vx = " << cell.at(n).at(i).at(j-1).bar_Vx[0] << std::endl <<
							"upper bar_Vx = " << cell.at(n).at(i).at(j+1).bar_Vx[0] << std::endl <<
							"Vx       = " << cell.at(n+1).at(i).at(j).Vx[0] << std::endl <<
							"lower Vx = " << cell.at(n+1).at(i).at(j-1).Vx[0] << std::endl <<
							"upper Vx = " << cell.at(n+1).at(i).at(j+1).Vx[0] << std::endl <<
							"rho       = " << cell.at(n+1).at(i).at(j).rho << std::endl <<
							"lower rho = " << cell.at(n+1).at(i).at(j-1).rho << std::endl <<
							"upper rho = " << cell.at(n+1).at(i).at(j+1).rho << std::endl <<
							"prev rho       = " << cell.at(n).at(i).at(j).rho << std::endl <<
							"lower prev rho = " << cell.at(n).at(i).at(j-1).rho << std::endl <<
							"upper prev rho = " << cell.at(n).at(i).at(j+1).rho << std::endl <<
							"dM       = {" << cell.at(n).at(i).at(j).dM[1] << ", " << cell.at(n).at(i).at(j).dM[2] << ", " << cell.at(n).at(i).at(j).dM[3] << ", " << cell.at(n).at(i).at(j).dM[4] << "} " << std::endl <<
							"lower dM = {" << cell.at(n).at(i).at(j-1).dM[1] << ", " << cell.at(n).at(i).at(j-1).dM[2] << ", " << cell.at(n).at(i).at(j-1).dM[3] << ", " << cell.at(n).at(i).at(j-1).dM[4] << "} " << std::endl <<
							"upper dM = {" << cell.at(n).at(i).at(j+1).dM[1] << ", " << cell.at(n).at(i).at(j+1).dM[2] << ", " << cell.at(n).at(i).at(j+1).dM[3] << ", " << cell.at(n).at(i).at(j+1).dM[4] << "} " << std::endl << std::endl;
							for (int k = 1; k < 5; k++)
								if (cell.at(n).at(i).at(j).dM[k] != 0)
									std::cout << "dM[" << k << "] != 0" << std::endl;
							if (cell.at(n+1).at(i).at(j).rho != cell.at(n+1).at(i).at(j+1).rho)
								std::cout << "rho(j) != rho(j+1)" << std::endl;
							if (cell.at(n).at(i).at(j).bar_e != cell.at(n).at(i).at(j+1).bar_e)
								std::cout << "bar_e(j) != bar_e(j+1)" << std::endl;
					getchar();
				}
			}
		}
	}
}

void debug_p_output(int n, int nArray, int *i, int *j, cell2d cell)
{
	for (int iter = 0; iter < nArray; iter++) {
		std::cout << "For i = " << i[iter] << ", j =  " << j[iter] << ", P[0] = " << cell.at(n+1).at(i[iter]).at(j[iter]).P[0] << std::endl;
	}
	std::cout << std::endl;
}

void debug_dM_rho_output(int n, int nArray, int *i, int *j, cell2d cell)
{
	for (int iter = 0; iter < nArray; iter++) {
		std::cout << "For i = " << i[iter] << ", j = " << j[iter] << ", dM = {" << cell.at(n).at(i[iter]).at(j[iter]).dM[1] << ", " << cell.at(n).at(i[iter]).at(j[iter]).dM[2] << ", " << cell.at(n).at(i[iter]).at(j[iter]).dM[3] << ", " << cell.at(n).at(i[iter]).at(j[iter]).dM[4] << "}" << std::endl <<
			"rho = " << std::setiosflags(std::ios::fixed) << std::setprecision(16) << cell.at(n+1).at(i[iter]).at(j[iter]).rho << std::endl;
		}
	std::cout << std::endl;
}

void debug_Vx_Vr_P_A_barVx_output(int n, int nArray, int *i, int *j, cell2d cell)
{
	for (int iter = 0; iter < nArray; iter++) {
		std::cout << "For i = " << i[iter] << ", j = " << j[iter] << std::endl <<
			"Vx = {" << cell.at(n).at(i[iter]).at(j[iter]).Vx[1] << ", " << cell.at(n).at(i[iter]).at(j[iter]).Vx[2] << ", " << cell.at(n).at(i[iter]).at(j[iter]).Vx[3] << ", " << cell.at(n).at(i[iter]).at(j[iter]).Vx[4] << "}" << std::endl <<
			"Vr = {" << cell.at(n).at(i[iter]).at(j[iter]).Vr[1] << ", " << cell.at(n).at(i[iter]).at(j[iter]).Vr[2] << ", " << cell.at(n).at(i[iter]).at(j[iter]).Vr[3] << ", " << cell.at(n).at(i[iter]).at(j[iter]).Vr[4] << "}" << std::endl <<
			"P = {" << cell.at(n).at(i[iter]).at(j[iter]).P[0] << ", "  << cell.at(n).at(i[iter]).at(j[iter]).P[1] << ", " << cell.at(n).at(i[iter]).at(j[iter]).P[2] << ", " << cell.at(n).at(i[iter]).at(j[iter]).P[3] << ", " << cell.at(n).at(i[iter]).at(j[iter]).P[4] << "}" << std::endl <<
			"A = {" << cell.at(n).at(i[iter]).at(j[iter]).A[1] << ", " << cell.at(n).at(i[iter]).at(j[iter]).A[2] << ", " << cell.at(n).at(i[iter]).at(j[iter]).A[3] << ", " << cell.at(n).at(i[iter]).at(j[iter]).A[4] << "}" << std::endl <<
			"bar_Vx = " << cell.at(n).at(i[iter]).at(j[iter]).bar_Vx[0] << ", bar_Vr = " << cell.at(n).at(i[iter]).at(j[iter]).bar_Vr[0] << ", bar_e = " << cell.at(n).at(i[iter]).at(j[iter]).bar_e << std::endl << std::endl;
	}
	std::cout << std::endl;
}

void debug_final_output(int n, int nArray, int *i, int *j, cell2d cell) {
	for (int iter = 0; iter < nArray; iter++) {
		printf("E at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n+1).at(i[iter]).at(j[iter]).e);
		printf("rho at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n+1).at(i[iter]).at(j[iter]).rho);
		printf("dM[1] at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n).at(i[iter]).at(j[iter]).dM[1]);
		printf("dM[2] at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n).at(i[iter]).at(j[iter]).dM[2]);
		printf("dM[3] at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n).at(i[iter]).at(j[iter]).dM[3]);
		printf("dM[4] at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n).at(i[iter]).at(j[iter]).dM[4]);
		printf("barVx at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n).at(i[iter]).at(j[iter]).bar_Vx[0]);
		printf("barVr at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n).at(i[iter]).at(j[iter]).bar_Vr[0]);
		printf("Vx at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n+1).at(i[iter]).at(j[iter]).Vx[0]);
		printf("Vr at %d:%d = %16.16f\n",
				i[iter],j[iter],cell.at(n+1).at(i[iter]).at(j[iter]).Vr[0]);
	}
}
