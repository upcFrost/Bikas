/* 
 * File:   main.cpp
 * Author: Frost
 *
 * Created on October 27, 2012, 4:49 PM
 */


#include "main.h"
#include "globals.h"
#include "interp.h"
#include <ctime>

using namespace std;


void finishInit(double e_0) {
	printf("Pre-init complete\n\n");
	cout << "Initial values: " << endl
		<< "P_vsp, Pa: " << P_v << endl
		<< "e_0, Dzh: " << e_0 << endl
		<< "psi_0: " << (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1) << endl
		<< "f, : " << f << endl
		<< "I_k, : " << I_k << endl
		<< "m_sn = " << m_sn << endl
		<< "k = " << k << endl
		<< "kappa = " << kappa << endl
		<< "lambda = " << lambda << endl
		<< "d = " << 2*(max_j-4)*dr << endl
		<< "S = " << M_PI*pow((max_j-4)*dr,2) << endl
		<< "S_km = " << (i_sn_0-1)*dx * 2*M_PI*(max_j-4)*dr << endl
		<< "V_km = " << (i_sn_0-1)*dx * M_PI*pow((max_j-4)*dr,2) << endl
		<< "omega = " << (i_sn_0-1)*dx * M_PI*pow((max_j-4)*dr,2) * delta_0 << endl
		<< "rho_0 = " << delta_0 << endl 
		<< "axis_j = " << axis_j << endl 
		<< endl;
}



int main(int argc, char** argv) {
    cell2d cell;
    string line;
    string tmpLine;
    bool verbose = false;
    double timestep = 0;
    
    t.resize(1);
    x_sn.resize(1);
    U_sn.resize(1);
    /* Open input file and create output files */
    ifstream inputFile ("input.csv");
    ofstream outputGas ("outputGas.csv");
	ofstream outputDyn ("outputDyn.csv");

	ofstream pvdGas ("outputGas.pvd");
	pvdGas.is_open();
	pvdGas << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl <<
	"	<Collection>" << endl;
	
    if (inputFile.is_open() && outputDyn.is_open() && outputGas.is_open())
    {
		/* Drop header line */
		getline (inputFile, line);
		
		/* Get initial vars from input file */
		int charNum = 0;
		getline (inputFile, line);
		stringstream varsLine(line);
		while(getline (varsLine, tmpLine, ',')) {
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
		U_sn.at(0) = 0;
		
		/* Init cell vector */
		n = 0;
		cell.resize(1);
		cell[n].resize(max_i);
		
		for (i = 0; i < max_i; i++) 
		{
		    cell[n][i].resize(max_j);
		}
		
		/* Skip header line */
		getline(inputFile, line);
		
		/* Scaling */
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
		
		/* Init cell array */
		while (inputFile.good()) {
		    
		    /* Set types, alpha, A and F from input file */
		    getline(inputFile, line);
		    charNum = 0;
		    stringstream varsLine(line);
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
			
		    cell.at(n).at(i).at(j).alpha = pre_alpha(cell.at(n).at(i).at(j), i, max_i, delta_0, delta);
		    double A[5];
		    pre_cell_geometry(A, cell.at(n).at(i).at(j), i, j);
		    if (cell.at(n).at(i).at(j).type == 1 || cell.at(n).at(i).at(j).type == 3 ||
					cell.at(n).at(i).at(j).type == 8 || cell.at(n).at(i).at(j).type == 10 )	{
				printf("Making weightVector for %d:%d\n",i,j);
				cell.at(n).at(i).at(j).weightVector = wightVectorsCalc(cell, i, j, n, false);
				printf("weightVector size on %d for %d:%d = %u\n", n,i,j,
					(unsigned int)cell.at(n).at(i).at(j).weightVector.y.size());
			}
		    for (unsigned int iter = 0; iter < 5; iter++) 
				cell.at(n).at(i).at(j).A[iter] = A[iter];
		    i_sn = floor(x_sn.at(0)/dx);
		    i_sn_0 = i_sn;
		    max_z_0 = max_z;
	    
		    /* Set gas init parameters */
		    if (i < i_sn) {
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).P[iter] = P_v;
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).Vx[iter] = 0;
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).Vr[iter] = 0;
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).dM[iter] = 0;
				cell.at(n).at(i).at(j).rho = delta_0;
				cell.at(n).at(i).at(j).final_psi = (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1);
				cell.at(n).at(i).at(j).final_z = 2 * cell.at(n).at(i).at(j).final_psi / (kappa * (1 + sqrt(1 + 4*lambda*cell.at(n).at(i).at(j).final_psi/kappa)));
				//~ cell.at(n).at(i).at(j).e = cell.at(n).at(i).at(j).P[0] / (k-1) / delta_0; // Ideal gas
				//~ cell.at(n).at(i).at(j).e = cell.at(n).at(i).at(j).P[0] / (k-1) * (1/delta_0 - 1/delta); // BMSTU var
				cell.at(n).at(i).at(j).e = cell.at(n).at(i).at(j).P[0] / (k-1) * (1/cell.at(n).at(i).at(j).rho - (1 - cell.at(n).at(i).at(j).final_psi)/delta - alpha_k * cell.at(n).at(i).at(j).final_psi); // From BMSTU P equation
				cell.at(n).at(i).at(j).psi = (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1);
				cell.at(n).at(i).at(j).z = 2 * cell.at(n).at(i).at(j).psi / (kappa * (1 + sqrt(1 + 4*lambda*cell.at(n).at(i).at(j).psi/kappa)));
			} else {
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).P[iter] = P_atm;
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).Vx[iter] = 0;
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).Vr[iter] = 0;
				for (int iter = 0; iter < 5; iter++) cell.at(n).at(i).at(j).dM[iter] = 0;
				cell.at(n).at(i).at(j).rho = rho_atm; // Air
				//~ cell.at(n).at(i).at(j).e = P_atm * 2*M_PI*(j*dr + dr/2)*dx*dr / 0.4; // Used for powder eq
				cell.at(n).at(i).at(j).e = P_atm / rho_atm / (k - 1); // Ideal gas
				cell.at(n).at(i).at(j).psi = 1;
				cell.at(n).at(i).at(j).z = 1;
				
				cell.at(n).at(i).at(j).final_psi = 1;
				cell.at(n).at(i).at(j).final_z = 1;
			}
		}
		
		/** Test - Riemann problem 1 **/
		//~ i_sn = max_i - 10;
		//~ x_sn.back() = i_sn * dx;
				
		/* Prepare output files */
		prepOutputDynCSV(outputDyn);
		if (verbose) {
			outputGas << "t,i,j,P[0],rho,e,Vx[0],Vr[0],bar_Vx[0],bar_Vr[0],bar_e,m,z,psi,dM[1],dM[2],dM[3],dM[4],A[0],A[1],A[2],A[3],A[4],IntE" << endl;
		} else {
			outputGas << "t,x,y,z,P[0],rho,e,Vx[0],Vr[0],z,psi,IntE" << endl;
		}
		
		/* Set A for projectile border - needed if starting with non-full cell */
	    double arrayT[5] = {0};
	    for (j = 2; j < max_j-2; j++) {
		    euler_proj_broder(arrayT, j, x_sn.back(), dx, dr);
			// For n
			cell.at(n).at(i_sn-1).at(j).A[0] = arrayT[0];
			cell.at(n).at(i_sn-1).at(j).A[1] = arrayT[1];
			cell.at(n).at(i_sn-1).at(j).A[2] = arrayT[2];
			cell.at(n).at(i_sn-1).at(j).A[3] = arrayT[3];
			cell.at(n).at(i_sn-1).at(j).A[4] = arrayT[4];
			if (j == max_j - 2) cell.at(n).at(i_sn-1).at(j).A[4] = 0;
			if (j == 2) cell.at(n).at(i_sn-1).at(j).A[3] = 0;
		}
		
		
		/* Add one more position as as ending */
		printf("weightVector size on %d for 99:13 = %u\n", 
			n, (unsigned int)cell.at(n).at(99).at(13).weightVector.y.size());
		cell2dStatic currentCell = cell.back();
		cell.push_back(currentCell);
		printf("weightVector size on %d for 99:13 = %u\n", 
			n+1, (unsigned int)cell.at(n+1).at(99).at(13).weightVector.y.size());

		finishInit(cell.at(n).at(5).at(5).e);
		
		/**************
		 *  Main loop *
		 **************/
		float Ku;
		int iter_count;

		if (argc > 1) {
			for (int argNum = 1; argNum < argc; argNum += 2) {
				if (strcmp (argv[argNum], "-Ku") == 0) {
					Ku = atof(argv[argNum+1]);
				} else if (strcmp (argv[argNum], "-iter") == 0) {
					iter_count = atoi(argv[argNum+1]);
				} else {
					cout << "Invalid parameters.\nValid usage: -Ku [float] -iter [float]";
					exit (EXIT_FAILURE);
				}
			}
		} else {
			try {
				cout << "Max iterations: \n";
				cin >> iter_count;
				if (iter_count <= 0) iter_count = 400;
			} catch (int e) {
				iter_count = 400;
			}
			cout << "Ku: \n";
			try {
				cin >> Ku;
				if (Ku < 0) Ku = 0.05;
			} catch (int e) {
				Ku = 0.05;
			}
		}
		
		clock_t start = clock();
		clock_t stop;
		double speed = 0;
		bool need_out;
		int iteration = 0;
		//~ while (x_sn.at(x_sn.size()-1) < (max_i-8)*dx) {
		while (iteration < iter_count) {
			/** Test - Riemann problem 2 **/
			//~ if (iteration == 10) x_sn.back() += dx * 0.05;
			//~ if (iteration == 11) x_sn.back() += dx * 0.05;
			//~ if (iteration == 12) x_sn.back() += dx * 0.05;
			//~ if (iteration == 13) x_sn.back() += dx * 0.05;
			
			need_out = false;
			if (iteration % 1000 == 0 && iteration > 0) {
				stop = clock();
				speed = (stop - start) / CLOCKS_PER_SEC;
				start = clock();
			}

			cell2dStatic currentCell = cell.back();
			cell.push_back(currentCell);
			n = cell.size() - 2;
			
			/* Projectile position calculation */
			/** TODO: fix P_sn **/
			int i_sn_prev = i_sn;
			double P_sn = cell.at(n).at(i_sn-1).at(4).P[0]; // Just 4 :)
			U_sn.push_back(euler_Usn(P_sn, 3.1415*pow((max_j-4)*dr,2), 0, dt, U_sn.back()));
			x_sn.push_back(euler_Xsn(x_sn.back(), U_sn.back()));
			i_sn = floor(x_sn.back() / dx);
			double arrayT[5];
			bool changed = false;
			if (i_sn_prev != i_sn) {
				changed = true;
				cout << "i_sn changed" << endl;
				for (j = 0; j < max_j; j++) {					
					/* Return cell to its original shape */
					// For n
					double tempArray[5];
					pre_cell_geometry(tempArray, cell.at(n).at(i_sn_prev-1).at(j), i_sn_prev-1, j);
					
					cell.at(n).at(i_sn_prev-1).at(j).A[0] = tempArray[0];
					cell.at(n).at(i_sn_prev-1).at(j).A[1] = tempArray[1];
					cell.at(n).at(i_sn_prev-1).at(j).A[2] = tempArray[2];
					cell.at(n).at(i_sn_prev-1).at(j).A[3] = tempArray[3];
					cell.at(n).at(i_sn_prev-1).at(j).A[4] = tempArray[4];
					// For n+1
					cell.at(n+1).at(i_sn_prev-1).at(j).A[0] = tempArray[0];
					cell.at(n+1).at(i_sn_prev-1).at(j).A[1] = tempArray[1];
					cell.at(n+1).at(i_sn_prev-1).at(j).A[2] = tempArray[2];
					cell.at(n+1).at(i_sn_prev-1).at(j).A[3] = tempArray[3];
					cell.at(n+1).at(i_sn_prev-1).at(j).A[4] = tempArray[4];

					cell.at(n).at(i_sn-1).at(j).P[0] = cell.at(n).at(i_sn_prev-1).at(j).P[0];
					cell.at(n).at(i_sn-1).at(j).P[1] = cell.at(n).at(i_sn_prev-1).at(j).P[1];
					cell.at(n).at(i_sn-1).at(j).P[2] = cell.at(n).at(i_sn_prev-1).at(j).P[2];
					cell.at(n).at(i_sn-1).at(j).P[3] = cell.at(n).at(i_sn_prev-1).at(j).P[3];
					cell.at(n).at(i_sn-1).at(j).P[4] = cell.at(n).at(i_sn_prev-1).at(j).P[4];
					cell.at(n).at(i_sn-1).at(j).rho = cell.at(n).at(i_sn_prev-1).at(j).rho;
					cell.at(n).at(i_sn-1).at(j).e = cell.at(n).at(i_sn_prev-1).at(j).e;
					cell.at(n).at(i_sn-1).at(j).Vx[0] = cell.at(n).at(i_sn_prev-1).at(j).Vx[0];
					cell.at(n).at(i_sn-1).at(j).Vx[1] = cell.at(n).at(i_sn_prev-1).at(j).Vx[1];
					cell.at(n).at(i_sn-1).at(j).Vx[2] = cell.at(n).at(i_sn_prev-1).at(j).Vx[2];
					cell.at(n).at(i_sn-1).at(j).Vx[3] = cell.at(n).at(i_sn_prev-1).at(j).Vx[3];
					cell.at(n).at(i_sn-1).at(j).Vx[4] = cell.at(n).at(i_sn_prev-1).at(j).Vx[4];
					cell.at(n).at(i_sn-1).at(j).Vr[0] = cell.at(n).at(i_sn_prev-1).at(j).Vr[0];
					cell.at(n).at(i_sn-1).at(j).Vr[1] = cell.at(n).at(i_sn_prev-1).at(j).Vr[1];
					cell.at(n).at(i_sn-1).at(j).Vr[2] = cell.at(n).at(i_sn_prev-1).at(j).Vr[2];
					cell.at(n).at(i_sn-1).at(j).Vr[3] = cell.at(n).at(i_sn_prev-1).at(j).Vr[3];
					cell.at(n).at(i_sn-1).at(j).Vr[4] = cell.at(n).at(i_sn_prev-1).at(j).Vr[4];
					cell.at(n).at(i_sn-1).at(j).final_z = cell.at(n).at(i_sn_prev-1).at(j).final_z;
					cell.at(n).at(i_sn-1).at(j).final_psi = cell.at(n).at(i_sn_prev-1).at(j).final_psi;

				}
				//~ for (j = 0; j < max_j; j++) {					
					//~ cell.at(n).at(i_sn).at(j).bar_Vx[0] = -cell.at(n).at(i_sn-1).at(j).bar_Vx[0];
					//~ cell.at(n).at(i_sn).at(j).bar_Vr[0] = cell.at(n).at(i_sn-1).at(j).bar_Vr[0];
					//~ cell.at(n).at(i_sn).at(j).P[0] = cell.at(n).at(i_sn-1).at(j).P[0];
					//~ cell.at(n).at(i_sn).at(j).Vx[0] = -cell.at(n).at(i_sn-1).at(j).Vx[0];
					//~ cell.at(n).at(i_sn).at(j).Vr[0] = cell.at(n).at(i_sn-1).at(j).Vr[0];
					//~ cell.at(n).at(i_sn).at(j).rho = cell.at(n).at(i_sn-1).at(j).rho;
					//~ cell.at(n).at(i_sn).at(j).e = cell.at(n).at(i_sn-1).at(j).e;
					//~ cell.at(n).at(i_sn).at(j).bar_e = cell.at(n).at(i_sn-1).at(j).bar_e;
				//~ }
			}
			
			/* Moving projectile borders */
			for (j = 0; j < max_j; j++) {
				if (cell.at(n).at(i_sn).at(j).type != 18) {
					euler_proj_broder(arrayT, j, x_sn.back(), dx, dr);
					// For n
					for (int iter = 0; iter < 5; iter++) {
						cell.at(n).at(i_sn-1).at(j).A[iter] = arrayT[iter];
						cell.at(n+1).at(i_sn-1).at(j).A[iter] = arrayT[iter];
					}
					if (cell.at(n).at(i_sn-1).at(j+1).type == 18) cell.at(n).at(i_sn-1).at(j).A[4] = 0;
					if (cell.at(n).at(i_sn-1).at(j-1).type == 18) cell.at(n).at(i_sn-1).at(j).A[3] = 0;
					if (cell.at(n).at(i_sn-1).at(j+1).type == 18) cell.at(n+1).at(i_sn-1).at(j).A[4] = 0;
					if (cell.at(n).at(i_sn-1).at(j-1).type == 18) cell.at(n+1).at(i_sn-1).at(j).A[3] = 0;
					
					for (int iter = 0; iter < 5; iter++) {
						if (cell.at(n).at(i_sn-1).at(j).A[iter] >= 2) {
							cell.at(n).at(i_sn-1).at(j).A[iter] -= 1;
						}
					}
				}
			}
			
			// From Zenkin V.A. "Issledovanie gazodin processov v dizelyah", 05.04.02 - Teplovie dvigateli, str 9
			// Adjusting values
			/** TODO: Check borderP (Qi) and newP (ceil) **/
			for (j = 0; j < max_j; j++) {
				if (cell.at(n).at(i_sn-1).at(j).type != 18) {
					double barQi;
					double Qi;
					
					if (!changed) {
						if (j == 2) cout << "A[0] = " << cell.at(n).at(i_sn-1).at(j).A[0] << endl;
						if (j == 2) cout << "Prev A[0] = " << cell.at(n-1).at(i_sn-1).at(j).A[0] << endl;
						barQi = (cell.at(n).at(i_sn-1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
						Qi = (cell.at(n-1).at(i_sn-1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
					} else {
						//~ cell.at(n).at(i_sn-1).at(j).rho = (2*cell.at(n-1).at(i_sn_prev-1).at(j).rho - cell.at(n-2).at(i_sn_prev-1).at(j).rho) * cell.at(n).at(i_sn-1).at(j).A[0]/(cell.at(n-1).at(i_sn_prev-1).at(j).A[0]-1);
						if (j == 2) cout << "A[0] = " << cell.at(n).at(i_sn-1).at(j).A[0] << endl;
						if (j == 2) cout << "Prev A[0] = " << cell.at(n-1).at(i_sn-2).at(j).A[0] << endl;
						barQi = (cell.at(n).at(i_sn-1).at(j).A[0]) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2); // Right side is the full cell volume, so we'll get absolute value
						Qi = (cell.at(n-1).at(i_sn-2).at(j).A[0]-1) * M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
					}
					
					double ai = sqrt(k*cell.at(n).at(i_sn-1).at(j).P[0]/cell.at(n).at(i_sn-1).at(j).rho); // Local speed of sound
					double borderP = cell.at(n).at(i_sn-1).at(j).P[0] + ai*cell.at(n).at(i_sn-1).at(j).rho * (U_sn.back() - cell.at(n).at(i_sn-1).at(j).Vx[0]);
					double newRho = cell.at(n).at(i_sn-1).at(j).rho * Qi / barQi;
					double newVx = cell.at(n).at(i_sn-1).at(j).rho / newRho * Qi / barQi * cell.at(n).at(i_sn-1).at(j).Vx[0] + (borderP - cell.at(n).at(i_sn-1).at(j).P[0]) / newRho / barQi * dt * M_PI*(2*(j-axis_j)+1)*pow(dr,2);
					double newE = cell.at(n).at(i_sn-1).at(j).e + borderP * (Qi - barQi) / newRho / barQi;
					double newP = (k-1) * (newE - pow(newVx,2)/2 - pow(cell.at(n).at(i_sn-1).at(j).Vr[0],2)) / ( 1/newRho - (1 - cell.at(n).at(i_sn-1).at(j).final_psi)/delta - alpha_k * cell.at(n).at(i_sn-1).at(j).final_psi); // BMSTU var
					//~ double newP = (k-1) * (newE - pow(newVx,2)/2) * newRho; // If no powder present

					/** DEBUG **/
					if (j == 2) cout << endl << "Old rho = " << cell.at(n).at(i_sn-1).at(j).rho << endl;
					if (j == 2) cout << "Qi = " << Qi << endl;
					if (j == 2) cout << "barQi = " << barQi << endl;
					if (j == 2) debug_projectile_par(i_sn, j, borderP, newRho, newVx, newE, newP, cell.at(n).at(i_sn-1).at(j).final_psi, x_sn.back());
					
					cell.at(n).at(i_sn-1).at(j).P[0] = newP;
					//~ cell.at(n).at(i_sn-1).at(j).P[2] = borderP;
					cell.at(n).at(i_sn-1).at(j).e = newE;
					cell.at(n).at(i_sn-1).at(j).Vx[0] = newVx;
					cell.at(n).at(i_sn-1).at(j).rho = newRho;
					

				}
			}
			
			//~ for (j = 0; j < max_j; j++) {
				//~ cell.at(n).at(i_sn).at(j).P[0] = cell.at(n).at(i_sn-1).at(j).P[0];
				//~ cell.at(n).at(i_sn).at(j).Vr[0] = cell.at(n).at(i_sn-1).at(j).Vr[0];
				//~ cell.at(n).at(i_sn).at(j).Vx[0] = -cell.at(n).at(i_sn-1).at(j).Vx[0];
				//~ cell.at(n).at(i_sn).at(j).rho = cell.at(n).at(i_sn-1).at(j).rho;
				//~ cell.at(n).at(i_sn).at(j).e = cell.at(n).at(i_sn-1).at(j).e;
				//~ cell.at(n).at(i_sn).at(j).final_z = cell.at(n).at(i_sn-1).at(j).final_z;
				//~ cell.at(n).at(i_sn).at(j).final_psi = cell.at(n).at(i_sn-1).at(j).final_psi;
//~ 
				//~ cell.at(n+1).at(i_sn).at(j).P[0] = cell.at(n+1).at(i_sn-1).at(j).P[0];
				//~ cell.at(n+1).at(i_sn).at(j).Vr[0] = cell.at(n+1).at(i_sn-1).at(j).Vr[0];
				//~ cell.at(n+1).at(i_sn).at(j).Vx[0] = -cell.at(n+1).at(i_sn-1).at(j).Vx[0];
				//~ cell.at(n+1).at(i_sn).at(j).rho = cell.at(n+1).at(i_sn-1).at(j).rho;
				//~ cell.at(n+1).at(i_sn).at(j).e = cell.at(n+1).at(i_sn-1).at(j).e;
			//~ }			
			
			
			
			/* Euler stage */
			/** TODO: change i_sn to max_i **/
			for (i = 1; i < i_sn; i++) {
				for (j = 1; j < max_j-1; j++) {
					if (cell.at(n).at(i).at(j).type != 18) {
						cell.at(n).at(i).at(j).bar_z = euler_z(&cell, &cell.at(n).at(i).at(j), n, i, j);
						cell.at(n).at(i).at(j).bar_psi = euler_psi(&cell, &cell.at(n).at(i).at(j), n, i, j);
						cell.at(n).at(i).at(j).bar_Vx[0] = euler_bar_Vx(&cell,n,i,j,dt,dx,dr);
						cell.at(n).at(i).at(j).bar_Vr[0] = euler_bar_Vr(&cell,n,i,j,dt,dx,dr);
						//~ if (cell.at(n).at(i).at(j).type == 1 || cell.at(n).at(i).at(j).type == 3 ||
						//~ 		cell.at(n).at(i).at(j).type == 8 || cell.at(n).at(i).at(j).type == 10 )	{
						//~ 	printf("Before rotation: %10.10f:%10.10f\n",
						//~ 		cell.at(n).at(i).at(j).bar_Vx[0], cell.at(n).at(i).at(j).bar_Vr[0]);
						//~ 	rotateVectors(cell.at(n).at(i).at(j).bar_Vx[0],
						//~ 		cell.at(n).at(i).at(j).bar_Vr[0], cell.at(n).at(i).at(j).angle);
						//~ 	printf("After rotation: %10.10f:%10.10f\n",
						//~ 		cell.at(n).at(i).at(j).bar_Vx[0], cell.at(n).at(i).at(j).bar_Vr[0]);
						//~ }
						cell.at(n).at(i).at(j).bar_e = euler_bar_e(&cell,n,i,j,dt,dx,dr);
					}
				}
			}
			
			/** DEBUG **/
			if (true) {
				int nArray = 1; int *iArray = new int[nArray]; int *jArray = new int[nArray];
				iArray[0] = 101; jArray[0] = 31;
				debug_Vx_Vr_P_A_barVx_output(n, nArray, iArray, jArray, cell);
			}

			/* Lagrange stage */
			for (i = 1; i < i_sn; i++) {
				for (j = 1; j < max_j-1; j++) {
					if (cell.at(n).at(i).at(j).type != 18) {
						cell.at(n).at(i).at(j).m = lagrange_m(&cell.at(n).at(i).at(j));
						double array[21] = {0};
						lagrange_mass(array, &cell, i, j, n, dx, dr, dt);
						cell.at(n).at(i).at(j).dM[0] = 0;
						cell.at(n).at(i).at(j).D[0] = 0;
						for (int iter = 1; iter < 5; iter++) {
							cell.at(n).at(i).at(j).dM[iter] = array[iter];
							cell.at(n).at(i).at(j).D[iter] = array[4+iter];
						}
						if (cell.at(n).at(i).at(j).A[1] == 0 && cell.at(n).at(i).at(j).dM[1] != 0) cout << "dM1 = " << cell.at(n).at(i).at(j).dM[1] << "\n";
						if (cell.at(n).at(i).at(j).A[2] == 0 && cell.at(n).at(i).at(j).dM[2] != 0) cout << "dM2 = " << cell.at(n).at(i).at(j).dM[2] << "\n";
						if (cell.at(n).at(i).at(j).A[3] == 0 && cell.at(n).at(i).at(j).dM[3] != 0) cout << "dM3 = " << cell.at(n).at(i).at(j).dM[3] << "\n";
						if (cell.at(n).at(i).at(j).A[4] == 0 && cell.at(n).at(i).at(j).dM[4] != 0) cout << "dM4 = " << cell.at(n).at(i).at(j).dM[4] << "\n";
						
						cell.at(n+1).at(i).at(j).rho = lagrange_rho(&cell.at(n).at(i).at(j),&cell.at(n-1).at(i).at(j),i,j,dt,dx,dr);
						if (cell.at(n+1).at(i).at(j).rho < 0) {
							cout << "i_sn = " << i_sn << endl
								<< "P[0] at i_sn = " << cell.at(n).at(i_sn).at(j).P[0] << endl;
							broken_dt = true;
						}
					}
				}
			}
			
			/** DEBUG **/
			if (true) {
				int nArray = 1; int *iArray = new int[nArray]; int *jArray = new int[nArray];
				iArray[0] = 101; jArray[0] = 31;
				debug_dM_rho_output(n, nArray, iArray, jArray, cell);			
			}
			
			
			/* Final stage */
			for (i = 1; i < i_sn; i++) {
				for (j = 1; j < max_j-1; j++) {
					if (cell.at(n).at(i).at(j).type != 18) {
						cell.at(n+1).at(i).at(j).Vx[0] = final_calc_Vx(&cell,i,j,n,dx,dr,dt);
						cell.at(n+1).at(i).at(j).Vr[0] = final_calc_Vr(&cell,i,j,n,dx,dr,dt);
						cell.at(n+1).at(i).at(j).e = final_calc_e(&cell,i,j,n,dx,dr,dt);
						if (cell.at(n+1).at(i).at(j).e < 0) need_out = true;
						
						cell.at(n+1).at(i).at(j).final_z = new_final_z(&cell,i,j,n,dx,dr,dt);
						cell.at(n+1).at(i).at(j).final_psi = new_final_psi(&cell,i,j,n,dx,dr,dt);
						
						// Post-final stage
						cell.at(n+1).at(i).at(j).z = final_calc_z(&cell, &cell.at(n+1).at(i).at(j), n, i, j);
						cell.at(n+1).at(i).at(j).psi = final_calc_psi(&cell, &cell.at(n+1).at(i).at(j), n, i, j);
						cell.at(n+1).at(i).at(j).P[0] = final_calc_p(&cell.at(n).at(i).at(j), &cell.at(n+1).at(i).at(j));
							
						//~ if (changed && i == i_sn-1 && j==6) {
							//~ cout << "Old P[0] = " << cell.at(n).at(i-1).at(j).P[0] << endl
								//~ << "New P[0] = " << cell.at(n+1).at(i).at(j).P[0] << endl
								//~ << "next cell P[0] = " << cell.at(n+1).at(i+1).at(j).P[0] << endl;
							//~ cout << "New rho = " << cell.at(n+1).at(i).at(j).rho << endl
								//~ << "Old rho = " << cell.at(n-1).at(i-1).at(j).rho << endl
								//~ << "New intE = " << (cell.at(n+1).at(i).at(j).e - (pow(cell.at(n+1).at(i).at(j).Vx[0],2)+pow(cell.at(n+1).at(i).at(j).Vr[0],2))/2) << endl
								//~ << "Old intE = " << (cell.at(n-1).at(i-1).at(j).e - (pow(cell.at(n-1).at(i).at(j).Vx[0],2)+pow(cell.at(n-1).at(i).at(j).Vr[0],2))/2) << endl
								//~ << "New final_psi = " << cell.at(n+1).at(i).at(j).final_psi << endl
								//~ << "Old final_psi = " << cell.at(n-1).at(i-1).at(j).final_psi << endl;
						//~ }
					}
				}
			}
			//~ /* P, Vx, Vr and E border conditions for j = 1, max_j - 2, i = 0 and i = i_sn*/
			//~ for (i = 1; i < max_i-1; i++) {
				//~ for (j = 1; j < max_j-1; j++) {
					//~ switch (cell.at(n+1).at(i).at(j).type) {
						//~ // type 1 --> top border closed, left and right partially closed
						//~ case 1:
							//~ cell.at(n+1).at(i).at(j+1).P[0] = 0;
							//~ for (unsigned int idx2 = 0; idx2 < cell.at(n+1).at(i).at(j+1).weightVector.size(); idx2++) {
								//~ cell.at(n+1).at(i).at(j+1).P[0] += cell.at(n+1).at(i).at(j+1).weightVector[idx2].weight * cell.at(n+1).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].i).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].j).P[0];
							//~ }
							//~ cell.at(n+1).at(i).at(j+1).Vx[0] = 0;
							//~ for (unsigned int idx2 = 0; idx2 < cell.at(n+1).at(i).at(j+1).weightVector.size(); idx2++) {
								//~ cell.at(n+1).at(i).at(j+1).Vx[0] += cell.at(n+1).at(i).at(j+1).weightVector[idx2].weight * cell.at(n+1).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].i).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].j).Vx[0];
							//~ }
							//~ cell.at(n+1).at(i).at(j+1).Vr[0] = 0;
							//~ for (unsigned int idx2 = 0; idx2 < cell.at(n+1).at(i).at(j+1).weightVector.size(); idx2++) {
								//~ cell.at(n+1).at(i).at(j+1).Vr[0] += cell.at(n+1).at(i).at(j+1).weightVector[idx2].weight * cell.at(n+1).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].i).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].j).Vr[0];
							//~ }
							//~ cell.at(n+1).at(i).at(j+1).e = 0;
							//~ for (unsigned int idx2 = 0; idx2 < cell.at(n+1).at(i).at(j+1).weightVector.size(); idx2++) {
								//~ cell.at(n+1).at(i).at(j+1).e += cell.at(n+1).at(i).at(j+1).weightVector[idx2].weight * cell.at(n+1).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].i).at(cell.at(n+1).at(i).at(j+1).weightVector[idx2].j).e;
							//~ }
							//~ break;
						//~ 
						//~ // type 15 --> top and left borders closed
						//~ case 15:
							//~ cell.at(n+1).at(i).at(j+1).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i-1).at(j).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i).at(j+1).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i-1).at(j).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i).at(j+1).Vr[0] = -cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i-1).at(j).Vr[0] = -cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i).at(j+1).e = cell.at(n+1).at(i).at(j).e;
							//~ cell.at(n+1).at(i-1).at(j).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
						//~ 
						//~ // type 17 --> left border closed
						//~ case 17:
							//~ cell.at(n+1).at(i-1).at(j).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i-1).at(j).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i-1).at(j).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i-1).at(j).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
						//~ 
						//~ // type 16 --> bottom and left borders closed
						//~ case 16:
							//~ cell.at(n+1).at(i).at(j-1).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i-1).at(j).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i).at(j-1).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i-1).at(j).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i).at(j-1).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i-1).at(j).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i).at(j-1).e = cell.at(n+1).at(i).at(j).e;
							//~ cell.at(n+1).at(i-1).at(j).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
					//~ 
						//~ // type 13 --> top border closed
						//~ case 13:
							//~ cell.at(n+1).at(i).at(j+1).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i).at(j+1).Vx[0] = cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i).at(j+1).Vr[0] = -cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i).at(j+1).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
					//~ 
						//~ // type 14 --> bottom border closed
						//~ case 14:
							//~ cell.at(n+1).at(i).at(j-1).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i).at(j-1).Vx[0] = cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i).at(j-1).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i).at(j-1).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
						//~ 
						//~ // type 19 --> right border closed
						//~ case 19:
							//~ cell.at(n+1).at(i+1).at(j).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i+1).at(j).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i+1).at(j).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i+1).at(j).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
						//~ 
						//~ // type 21 --> bottom and right borders closed
						//~ case 21:
							//~ cell.at(n+1).at(i).at(j-1).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i+1).at(j).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i).at(j-1).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i+1).at(j).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i).at(j-1).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i+1).at(j).Vr[0] = cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i).at(j-1).e = cell.at(n+1).at(i).at(j).e;
							//~ cell.at(n+1).at(i+1).at(j).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
							//~ 
						//~ // type 20 --> top and right borders closed
						//~ case 20:
							//~ cell.at(n+1).at(i).at(j+1).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i+1).at(j).P[0] = cell.at(n+1).at(i).at(j).P[0];
							//~ cell.at(n+1).at(i).at(j+1).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i+1).at(j).Vx[0] = -cell.at(n+1).at(i).at(j).Vx[0];
							//~ cell.at(n+1).at(i).at(j+1).Vr[0] = -cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i+1).at(j).Vr[0] = -cell.at(n+1).at(i).at(j).Vr[0];
							//~ cell.at(n+1).at(i).at(j+1).e = cell.at(n+1).at(i).at(j).e;
							//~ cell.at(n+1).at(i+1).at(j).e = cell.at(n+1).at(i).at(j).e;
							//~ break;
							//~ 
						//~ default:								
							//~ break;
					//~ }
				//~ }
			//~ }
			
			/** DEBUG **/
			if (true) {
				int nArray = 1; int *iArray = new int[nArray]; int *jArray = new int[nArray];
				iArray[0] = 101; jArray[0] = 31;
				debug_p_output(n, nArray, iArray, jArray, cell);
			}
					
			/** Smoothing **/
			//~ /* X axis */
			//~ int pointNum = 4; // Should be even, number of interpolating points
			//~ for (j = 0; j < max_j; j++) {
				//~ for (i = 1+pointNum/2; i < i_sn-1-pointNum/2; i++) {
					//~ if (cell.at(n+1).at(i).at(j).type != 18 && 
							//~ cell.at(n+1).at(i+1).at(j).type != 18 && 
							//~ cell.at(n+1).at(i+2).at(j).type != 18) {
						//~ double xPoints[pointNum]; 
						//~ double VxPoints[pointNum]; double VxPointsLin[pointNum]; double VxPointsQuad[pointNum]; double VxPointsCub[pointNum];
						//~ double VrPoints[pointNum]; double VrPointsLin[pointNum]; double VrPointsQuad[pointNum]; double VrPointsCub[pointNum]; 
						//~ double PPoints[pointNum]; double PPointsLin[pointNum]; double PPointsQuad[pointNum]; double PPointsCub[pointNum]; 
						//~ double ePoints[pointNum]; double ePointsLin[pointNum]; double ePointsQuad[pointNum]; double ePointsCub[pointNum]; 
						//~ for (int iter = 0; iter < pointNum; iter++) {
							//~ int current_i = iter < pointNum/2 ? i - pointNum/2 + iter : i - pointNum/2 + iter + 1;
							//~ xPoints[iter] = current_i * dx;
							//~ VxPoints[iter] = cell.at(n+1).at(current_i).at(j).Vx[0];
							//~ VrPoints[iter] = cell.at(n+1).at(current_i).at(j).Vr[0];
							//~ PPoints[iter] = cell.at(n+1).at(current_i).at(j).P[0];
							//~ ePoints[iter] = cell.at(n+1).at(current_i).at(j).e;
						//~ }
						//~ cubic_nak ( pointNum, xPoints, VxPoints, VxPointsLin, VxPointsQuad, VxPointsCub );
						//~ cubic_nak ( pointNum, xPoints, VrPoints, VrPointsLin, VrPointsQuad, VrPointsCub );
						//~ cubic_nak ( pointNum, xPoints, PPoints, PPointsLin, PPointsQuad, PPointsCub );
						//~ cubic_nak ( pointNum, xPoints, ePoints, ePointsLin, ePointsQuad, ePointsCub );
						//~ 
						//~ cell.at(n+1).at(i).at(j).Vx[0] = spline_eval ( pointNum, xPoints, VxPoints, VxPointsLin, VxPointsQuad, VxPointsCub, i*dx );
						//~ cell.at(n+1).at(i).at(j).Vr[0] = spline_eval ( pointNum, xPoints, VrPoints, VrPointsLin, VrPointsQuad, VrPointsCub, i*dx );
						//~ cell.at(n+1).at(i).at(j).P[0] = spline_eval ( pointNum, xPoints, PPoints, PPointsLin, PPointsQuad, PPointsCub, i*dx );
						//~ cell.at(n+1).at(i).at(j).e = spline_eval ( pointNum, xPoints, ePoints, ePointsLin, ePointsQuad, ePointsCub, i*dx );
					//~ }
				//~ }
			//~ }
			/* Y axis */
			//~ for (i = 0; i < max_i; i++) {
				//~ for (j = 0; j < max_j; j++) {
					//~ switch (cell.at(n+1).at(i).at(j).type) {
						//~ // type 0 --> free cell
						//~ case 0:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
							//~ 
						//~ // type 15 --> top and left borders closed
						//~ case 15:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] - 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
						//~ 
						//~ // type 17 --> left border closed
						//~ case 17:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
						//~ 
						//~ // type 16 --> bottom and left borders closed
						//~ case 16:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
					//~ 
						//~ // type 13 --> top border closed
						//~ case 13:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] - 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
					//~ 
						//~ // type 14 --> bottom border closed
						//~ case 14:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
						//~ 
						//~ // type 19 --> right border closed
						//~ case 19:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
						//~ 
						//~ // type 21 --> bottom and right borders closed
						//~ case 21:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
							//~ 
						//~ // type 20 --> top and right borders closed
						//~ case 20:
							//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] - 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
							//~ break;
							//~ 
						//~ default:								
							//~ break;
					//~ }
				//~ }
			//~ }
			
			/*****************
			 * DEBUG: Vx and E test *
			 *****************/
			//~ debug_equality_Vx_e(i_sn, max_j, n, cell);

			/* Next dt calculation */
			t.push_back(t.at(t.size()-1) + dt);
			vector <double> array;
			vector <double> minimum;
			for (i = 0; i < max_i; i++) {
				array.resize(max_j);
				for (j = 0; j < max_j; j++) {
					array[j] = fmin(dx,dr) / ( sqrt( fabs ( 
							k * cell.at(n+1).at(i).at(j).P[0] / 
							cell.at(n+1).at(i).at(j).rho)) 
							+ 
							sqrt(pow(cell.at(n+1).at(i).at(j).Vx[0],2) + pow(cell.at(n+1).at(i).at(j).Vr[0],2))
						);
					//~ array[j] = 5*pow(10,-1) * dx*pow(dr,2)/(cell.at(n+1).at(i).at(j).Vx[0]*dr + cell.at(n+1).at(i).at(j).Vr[0]*dx);
					//~ if (array[j] < 0) array[j] = pow(10,-6);
				}
				minimum.push_back(*min_element(array.begin(), array.end()));
			}
			
			for (i = 101; i < 102; i++) {
				for (j = 31; j < 32; j++) {
					printf("E at %d:%d = %16.16f\nE at %d:%d = %16.16f\n",
						i,j,cell[n+1][i][j].e,
						i-5,j-3,cell[n+1][i-5][j-3].e);
					printf("rho at %d:%d = %16.16f\nrho at %d:%d = %16.16f\n",
						i,j,cell[n+1][i][j].rho,
						i-5,j-3,cell[n+1][i-5][j-3].rho);
					printf("dM[1] at %d:%d = %16.16f\ndM[1] at %d:%d = %16.16f\n",
						i,j,cell[n][i][j].dM[1],
						i-5,j-3,cell[n][i-5][j-3].dM[1]);
					printf("dM[2] at %d:%d = %16.16f\ndM[2] at %d:%d = %16.16f\n",
						i,j,cell[n][i][j].dM[2],
						i-5,j-3,cell[n][i-5][j-3].dM[2]);
					printf("dM[3] at %d:%d = %16.16f\ndM[3] at %d:%d = %16.16f\n",
						i,j,cell[n][i][j].dM[3],
						i-5,j-3,cell[n][i-5][j-3].dM[3]);
					printf("dM[4] at %d:%d = %16.16f\ndM[4] at %d:%d = %16.16f\n",
						i,j,cell[n][i][j].dM[4],
						i-5,j-3,cell[n][i-5][j-3].dM[4]);
					printf("barVx at %d:%d = %16.16f\nbarVx at %d:%d = %16.16f\n",
						i,j,cell[n][i][j].bar_Vx[0],
						i-5,j-3,cell[n][i-5][j-3].bar_Vx[0]);
					printf("barVr at %d:%d = %16.16f\nbarVr at %d:%d = %16.16f\n",
						i,j,cell[n][i][j].bar_Vr[0],
						i-5,j-3,cell[n][i-5][j-3].bar_Vr[0]);
				}
			}

			double next = Ku * *min_element(minimum.begin(), minimum.end());
			if (next > 1.1*dt) {
				dt = 1.1*dt;
			} else {
				dt = next;
			}
			
			if (dt > 0.5*pow(10,-4)/scaleT) dt = 0.5*pow(10,-4)/scaleT;
			cout << "\r Num: " << iteration << ", dt: " << dt << ", speed: " << speed << " sec / 1000 iter, x_sn: " << x_sn.at(x_sn.size()-1) << ", Ak: " << cell.at(n-1).at(i_sn-1).at(4).A[0]/cell.at(n).at(i_sn-1).at(4).A[0];
			
			
			
			/* Output to file */
			if (fabs(t.back() - timestep) > pow(10,-5) || need_out) {
				timestep = t.back();
				
				/* Dynamics - to csv */
				outputDynCSV(outputDyn, t.back(), i_sn, x_sn.back(), U_sn.back(), cell.at(n+1).at(1).at(4).P[0], cell.at(n+1).at(i_sn-1).at(4).P[0]);
				
				/* Gas dynamics - verbose to csv and non-verbose to pvd/vtp */
				if (verbose) {
					outputCSV(cell, outputGas);
				} else {
					ostringstream os;
					os << "./result/outputGas_" << iteration << ".vtp";
					string filename = os.str();
					
					pvdGas << setiosflags(ios::fixed) << setprecision(10) << 
						"<DataSet timestep=\""<< t.back() << 
						"\" group=\"\" part=\"0\" file=\"" 
						<< filename 
						<< "\"/>" << endl;
					
					OutputPVD(cell, filename);
				}
			}
			
			/* Memory cleaning */
			if (t.size() > 5) {
				cell.erase(cell.begin());
				t.erase(t.begin());
			}
			iteration++;
		}
		
		pvdGas << 	"	</Collection>" << endl << "</VTKFile>" << endl;
		
    };

    return 0;
}

double truncNdigit(double value, int N) {
	return floor((value*pow(10,N))+0.5)/pow(10, N);
}

double Bi(int n, int i, double x) {
	return factorial(n) / (factorial(i)*factorial(n-i)) * pow(x, i) * pow(1-x, n-i);		
}

double Bj(int m, int j, double y) {
	return factorial(m) / (factorial(j)*factorial(m-j)) * pow(y, j) * pow(1-y, m-j);
}

double splineEval(double x, double y, int n, int m, double ** k) {
	double result = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result += Bi(n,i,x) * Bj(m,j,y) * k[i][j];
			//~ cout << "Bi = " << Bi(n,i,x) << endl;
			//~ cout << "Bj = " << Bj(m,j,y) << endl;
			//~ cout << "k  = " << k[i][j] << endl;
		}
	}
	return result;
}

/* Yep, using static table */
unsigned int factorial(unsigned int n) {
  static const unsigned int table[] = {1, 1, 2, 6, 24, 120, 720,
    5040, 40320, 362880, 3628800, 39916800, 479001600};
  if (n >= sizeof table / sizeof *table) // if appropriate, omit test if NDEBUG
    throw "Range error";
  return table[n];
}


std::string base64_encode(unsigned char const* bytes_to_encode, unsigned int in_len) {
  std::string ret;
  int i = 0;
  int j = 0;
  unsigned char char_array_3[3];
  unsigned char char_array_4[4];

  while (in_len--) {
    char_array_3[i++] = *(bytes_to_encode++);
    if (i == 3) {
      char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
      char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
      char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
      char_array_4[3] = char_array_3[2] & 0x3f;

      for(i = 0; (i <4) ; i++)
        ret += base64_chars[char_array_4[i]];
      i = 0;
    }
  }

  if (i)
  {
    for(j = i; j < 3; j++)
      char_array_3[j] = '\0';

    char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
    char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
    char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
    char_array_4[3] = char_array_3[2] & 0x3f;

    for (j = 0; (j < i + 1); j++)
      ret += base64_chars[char_array_4[j]];

    while((i++ < 3))
      ret += '=';

  }

  return ret;

}
