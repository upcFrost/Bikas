#include <vector>
#include <cmath>

bool broken_dt = false;
bool mergedI = false;
int n = 0;
int max_i = 0;
int max_j = 0;
int I_k = 1.69*pow(10,6);
int i_sn_0;
int i_pist_0;
double max_z_0;
double delta = 0;
double delta_0 = 0;
double dx = 0;
double dr = 0;
double dt = pow(10,-8);
double lambda = 0;
double kappa = 1;
float p0 = 0;
float V0 = 0;
float dM0 = 0;
double i_v = 0;
double j_v = 0;
double k = 1.262;
double m_sn = 0.3;
double j_sn = 0;
double i_sn = 0;
double i_pist = 0;
double max_z = 1;
std::vector <double> t (0);
std::vector <double> x_sn (0);
std::vector <double> x_pist (0);
std::vector <double> U_sn (0);
std::vector <double> U_pist (0);
double f = 0.922*pow(10,6);
int P_f = 30000000;
double P_atm = 100000;
double rho_atm = 1.2754;
double P_v = 30000000;

/* STUB */
int Qr = 0;
int m = 0;
int beta = 1;
