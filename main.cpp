/* 
 * File:   main.cpp
 * Author: Frost
 *
 * Created on October 27, 2012, 4:49 PM
 */


#include "main.h"

using namespace std;


void finishInit(double e_0) {
	double psi_0 = (delta/delta_0 - 1) / (f * delta / P_v + alpha_k * delta - 1);
	double d = 2*(max_j-4)*dr;
	double S = M_PI*pow((max_j-4)*dr,2);
	double S_km = (i_sn_0-1)*dx * 2*M_PI*(max_j-4)*dr;
	double V_km = (i_sn_0-1)*dx * M_PI*pow((max_j-4)*dr,2);
	double omega = (i_sn_0-1)*dx * M_PI*pow((max_j-4)*dr,2) * delta_0;
	printf("Pre-init complete\n\n");
	printf("Initial values:\n"
			"P_vsp, Pa: %10.10f\n"
			"e_0, Dzh: %10.10f\n"
			"psi_0: %10.10f\n"
			"f: %10.10f\n"
			"I_k: %10.10f\n"
			"m_sn: %10.10f\n"
			"k: %10.10f\n"
			"kappa: %10.10f\n"
			"lambda: %10.10f\n"
			"d: %10.10f\n"
			"S: %10.10f\n"
			"S_km: %10.10f\n"
			"V_km: %10.10f\n"
			"omega: %10.10f\n"
			"rho_0: %10.10f\n"
			"axis_j: %d\n\n",
			P_v*scaleP, e_0*scaleE, psi_0, f*scaleFF,
			I_k*scaleIK, m_sn*scaleM,
			k, kappa, lambda, d*scaleD,
			S*scaleD*scaleD, S_km*scaleD*scaleD,
			V_km*scaleD*scaleD*scaleD, omega*scaleM,
			delta_0*scaleRho, axis_j);

	int numThreads = 0;
#pragma omp parallel
	{
//		numThreads = omp_get_num_threads();
	}
	printf("OpenMP loaded, number of threads: %d\n", numThreads);
}



int main(int argc, char** argv) {
    cell2d cell;
    string line;
    string tmpLine;
    bool verbose = false;
    int gasVar = IDEAL_GAS;
    double timestep = 0;
    bool havePiston = false;
    bool debug = false;
    
    t.resize(1);
    x_sn.resize(1);
    U_sn.resize(1);

    string inputFileName;
    cout << "Type input file name and path\n";
    cin >> inputFileName;
    /* Open input file and create output files */
    ifstream inputFile (inputFileName.c_str());
    ofstream outputGas ("outputGas.csv");
	ofstream outputDyn ("outputDyn.csv");

	ofstream pvdGas ("outputGas.pvd");
	pvdGas.is_open();
	pvdGas << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl <<
	"	<Collection>" << endl;
	
#ifdef _WIN32
	system("mkdir ./result");
#else
	system("mkdir ./result");
#endif

    if (inputFile.is_open() && outputDyn.is_open() && outputGas.is_open())
    {
    	printf("Which gas model is used?\n"
    			"0 - Ideal gas, 2 - Powder\n");
    	cin >> gasVar;
//    	while (gasVar != IDEAL_GAS && gasVar != POWDER_EQ) {
//    		printf("\nInvalid input, please try again\n");
//    		scanf("%d", &gasVar);
//    	}

    	int ifPiston = 0;
    	printf("Is plastic piston present?\n"
				"1 - yes, any other digit - no\n");
		cin >> ifPiston;
		if (ifPiston == 1) {
			havePiston = true;
		}

    	int debugOutput = 0;
		printf("\nDo you need debug output? (1 - yes, any other digit - no)\n");
		cin >> debugOutput;
		if (debugOutput == 1) {
			debug = true;
		}

    	init(inputFile, cell, gasVar, havePiston, debug);

		/** Test - Riemann problem 1 **/
    	int riemannTest = 0;
    	printf("\nPreform riemann test? (1 - yes, any other digit - no)\n");
    	cin >> riemannTest;
    	if (riemannTest == 1) {
    		i_sn = max_i - 10;
    		x_sn.back() = i_sn * dx;
    	}

				
		/* Prepare output files */
		prepOutputDynCSV(outputDyn);
		prepOutputGasCSV(outputGas,  verbose);
		
		/* Set A for projectile border - needed if starting with non-full cell */
		borderCellsFix(cell, havePiston);
		
		/* Add one more position as as ending */
		cell2dStatic currentCell = cell.back();
		cell.push_back(currentCell);

		finishInit(cell.at(n).at(5).at(5).e);
		
		/**************
		 *  Main loop *
		 **************/
		float Ku;
		int iter_count;

		printf("Max iterations: \n");
		cin >> iter_count;
		if (iter_count <= 0) iter_count = 400;

		printf("Ku: \n");
		cin >> Ku;
		if (Ku <= 0) Ku = 0.05;
		
		clock_t start = clock();
		clock_t stop;
		double speed = 0;
		bool need_out;
		int iteration = 0;
		while (x_sn.at(x_sn.size()-1) < (max_i-8)*dx) {
//		while (iteration < iter_count) {
			need_out = false;
			if (iteration % 100 == 0 && iteration > 0) {
				stop = clock();
				speed = (stop - start) / CLOCKS_PER_SEC;
				start = clock();
			}

			cell2dStatic currentCell = cell.back();
			cell.push_back(currentCell);
			n = cell.size() - 2;
			
			/* Projectile-related calculation */
			if (havePiston) {
				int i_pist_prev = i_pist;
				projCalc(cell, gasVar, i_pist, false, debug);
				pistonCalc(cell, i_pist_prev, i_pist, IDEAL_GAS, debug);
				projCalc(cell, IDEAL_GAS, i_sn, true, debug);
			} else {
				projCalc(cell, gasVar, i_sn, true, debug);
			}
			
			/* Euler stage */
#pragma omp parallel for num_threads(4) schedule(dynamic,1) collapse(2)
			for (int i = 1; i < i_sn; i++) {
				for (int j = 1; j < max_j-1; j++) {
					if (cell.at(n).at(i).at(j).type != 18) {
						if (havePiston && i == i_pist)
							continue;

						gasCell & curCell = cell.at(n).at(i).at(j);

						curCell.bar_z = euler_z(&cell, &cell.at(n).at(i).at(j), n, i, j);
						curCell.bar_psi = euler_psi(cell.at(n).at(i).at(j), n, i, j);
						curCell.bar_Vx[0] = euler_bar_Vx(cell,n,i,j,dt,dx,dr,FIRST_ORDER);
						curCell.bar_Vr[0] = euler_bar_Vr(cell,n,i,j,dt,dx,dr,FIRST_ORDER);
						curCell.bar_e = euler_bar_e(cell,n,i,j,dt,dx,dr,FIRST_ORDER);
					}
				}
			}

			/* Lagrange stage */
#pragma omp parallel for num_threads(4) schedule(dynamic,1) collapse(2)
			for (int i = 1; i < i_sn; i++) {
				for (int j = 1; j < max_j-1; j++) {
					if (cell.at(n).at(i).at(j).type != 18) {
						if (havePiston && i == i_pist)
							continue;

						gasCell & curCell = cell.at(n).at(i).at(j);
						gasCell & nextTCell = cell.at(n+1).at(i).at(j);

						double array[21] = {0};
						lagrange_mass(array, cell, i, j, n, dx, dr, dt);
						curCell.dM[0] = 0;
						curCell.D[0] = 0;
						for (int iter = 1; iter < 5; iter++) {
							curCell.dM[iter] = fabs(array[iter]);
							curCell.D[iter] = array[4+iter];
						}
						nextTCell.rho = lagrange_rho(&cell.at(n).at(i).at(j),
								&cell.at(n-1).at(i).at(j),i,j,dt,dx,dr);

						if (i == i_pist-1 && j == 10) {
							printf("dM = %10.10f, %10.10f, %10.10f, %10.10f, rho = %10.10f",
									curCell.dM[1],curCell.dM[2],curCell.dM[3],curCell.dM[4],curCell.rho);
							cout << " rho = " <<  curCell.rho << endl;
						}
					}
				}
			}

			/* Final stage */
#pragma omp parallel for num_threads(4) schedule(dynamic,1) collapse(2)
			for (int i = 1; i < i_sn; i++) {
				for (int j = 1; j < max_j-1; j++) {
					if (cell.at(n).at(i).at(j).type != 18) {
						if (havePiston && i == i_pist)
							continue;

						gasCell & nextTCell = cell.at(n+1).at(i).at(j);
						
						nextTCell.Vx[0] = final_calc_Vx(cell,i,j,n,dx,dr,dt);
						nextTCell.Vr[0] = final_calc_Vr(cell,i,j,n,dx,dr,dt);
						nextTCell.e = final_calc_e(cell,i,j,n,dx,dr,dt);
						nextTCell.final_z = new_final_z(cell,i,j,n,dx,dr,dt);
						nextTCell.final_psi = new_final_psi(cell,i,j,n,dx,dr,dt);
						
						// Post-final stage
						if (i < i_pist || !havePiston) {
							nextTCell.P[0] = final_calc_p(&cell.at(n).at(i).at(j), &cell.at(n+1).at(i).at(j),
								gasVar, i);
						} else if (i >= i_pist && i <= i_sn && havePiston) {
//							nextTCell.P[0] = final_calc_p(&cell.at(n).at(i).at(j), &cell.at(n+1).at(i).at(j),
//								PISTON, i);
							nextTCell.P[0] = final_calc_p(&cell.at(n).at(i).at(j), &cell.at(n+1).at(i).at(j),
								IDEAL_GAS, i);
						}
					}
				}
			}
					
			/*****************
			 * DEBUG: Vx and E test *
			 *****************/
			//~ debug_equality_Vx_e(i_sn, max_j, n, cell);

			/* Next dt calculation */
			t.push_back(t.at(t.size()-1) + dt);
			vector <double> array;
			vector <double> minimum;

			// Numerator
			float num = fminf(dx, dr);

			for (int i = 0; i < i_sn; i++) {
				array.resize(max_j);
				for (int j = 0; j < max_j; j++) {
					if (cell.at(n).at(i).at(j).A[0] != 0) {
						gasCell & curCell = cell.at(n+1).at(i).at(j);
						array[j] = fmin(static_cast<double> (dx),static_cast<double> (dr)) /
							(
								sqrt( fabs (k * curCell.P[0] /	curCell.rho))
								+
								sqrt(pow(curCell.Vx[0],2) + pow(curCell.Vr[0],2))
							);
					} else {
						array[j] = 1;
					}
				}
				minimum.push_back(*min_element(array.begin(), array.end()));
			}
			
			if (debug) {
				int debug_I = 87, debug_J = 15;
				if (true && debug_I != 0 && debug_J != 0) {
					int nArray = 3; int *iArray = new int[nArray]; int *jArray = new int[nArray];
					iArray[0] = debug_I; jArray[0] = debug_J;
					iArray[1] = debug_I+1; jArray[1] = debug_J;
					iArray[2] = debug_I+2; jArray[2] = debug_J;
					debug_Vx_Vr_P_A_barVx_output(n, nArray, iArray, jArray, cell);
					debug_dM_rho_output(n, nArray, iArray, jArray, cell);
					debug_p_output(n, nArray, iArray, jArray, cell);
					debug_final_output(n, nArray, iArray, jArray, cell);
				}
			}

			double next = Ku * *min_element(minimum.begin(), minimum.end());
			if (next > 1.1*dt) {
				dt = 1.1*dt;
			} else {
				dt = next;
			}
			
			if (dt > Ku*pow(10.0,-6)/scaleT) dt = Ku*pow(10.0,-6)/scaleT;
			cout << "\r Num: " << iteration << ", dt: " << dt << ", speed: " << speed << " sec / 100 iter, x_sn: " << x_sn.at(x_sn.size()-1) << ", Ak: " << cell.at(n-1).at(i_sn-1).at(4).A[0]/cell.at(n).at(i_sn-1).at(4).A[0];
			
			
			
			/* Output to file */
//			if (iteration % 25 == 0) {
			if (fabs(t.back() - timestep) > pow(10.0,-5)/scaleT) {
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
		getchar();
		getchar();
		getchar();
    };
    return 0;
}

double truncNdigit(double value, int N) {
	return floor((value*pow(10.0,N))+0.5)/pow(10.0, N);
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
