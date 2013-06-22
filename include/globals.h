/* 
 * File:   main.h
 * Author: Frost
 *
 * Created on October 27, 2012, 5:47 PM
 */

#ifndef GLOB_H
#define	GLOB_H

#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include "types.h"

unsigned int const P_POS		= 0;
unsigned int const VX_POS		= 1;
unsigned int const VR_POS		= 2;
unsigned int const E_POS		= 3;
unsigned int const RHO_POS		= 4;
unsigned int const BAR_VX_POS	= 5;
unsigned int const BAR_VR_POS	= 6;
unsigned int const BAR_E_POS	= 7;
unsigned int const BAR_Z_POS	= 8;
unsigned int const BAR_PSI_POS	= 9;

#define P_PAR		(0x1 << 0);
#define VX_PAR		(0x1 << 1);
#define VR_PAR		(0x1 << 2);
#define E_PAR		(0x1 << 3);
#define RHO_PAR		(0x1 << 4);
#define BAR_VX_PAR	(0x1 << 5);
#define BAR_VR_PAR	(0x1 << 6);
#define BAR_E_PAR	(0x1 << 7);
#define Z_PAR		(0x1 << 8);
#define PSI_PAR		(0x1 << 9);

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)));


#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

static const int FIRST_ORDER 		= 0;
static const int SECOND_ORDER 		= 1;
static const int FIRST_ORDER_NS		= 2;
static const int SECOND_ORDER_NS	= 3;

static const int IDEAL_GAS			= 0;
static const int ABEL_DUPRE			= 1;
static const int POWDER_EQ			= 2;
static const int PISTON				= 3;


extern int n, prevN, nextN, max_i, max_j, I_k, i_v, j_v, i_sn, i_pist, j_sn, i_sn_0, i_pist_0, P_f;
extern int maxN;
extern double delta, delta_0, dx, dr, dt, p0, lambda, kappa, k, m_sn, max_z,max_z_0, f;
extern double V0, P_atm, rho_atm, P_v;
extern double dM0;
extern bool broken_dt;
extern bool mergedI;
extern std::vector <double> t;
extern std::vector <double> x_sn;
extern std::vector <double> x_pist;
extern std::vector <double> U_sn;
extern std::vector <double> U_pist;

const int axis_j = 0; 

/* Scaling */
static const double scaleD = pow(10.0,	-2);
static const double scaleT = pow(10.0,	-5);
static const double scaleM = pow(10.0,	-6);
static const double scaleV = pow(10.0,	3);
static const double scaleE = pow(10.0,	0);
static const double scaleF = pow(10.0,	2);
static const double scaleP = pow(10.0,	6);
static const double scaleRho = pow(10.0,	0);
static const double scaleDM = pow(10.0,	1);
static const double scaleFF = pow(10.0,	6);
static const double scaleIK = pow(10.0,	1);
static const double scaleGasKAPPA = pow(10.0, -3);
static const double scaleGasMU = pow(10.0, -1);
static const double scaleGasCP = pow(10.0, 6);
static const double scaleGasALPHAK = pow(10.0, 0);

/* Gas constants */
static const double alpha_k = 0.001006 / scaleGasALPHAK;
static const double gasA = -0.666667;
//static const double gasMu = 0.001 / scaleGasMU; // Zenkin, str 44, "Dlya turbulentnih techeniy"
static const double gasMu = 2.85*pow(10.0,-5) / scaleGasMU; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasLambda = pow(10,-3);
//static const double gasLambda = pow(10.0,-2); // Zenkin, str 44, "Dlya turbulentnih techeniy"
//static const double gasKappa = 1.21 / scaleGasKAPPA; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasKappa = 1.21; // Решение задачи по определению эффективности многокамерного дульного тормоза
//static const double gasCp = 1650 / scaleGasCP; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasCp = 1650; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasPr = gasMu * gasCp / gasKappa;
static const double gasB = gasKappa / gasPr;

/* Piston constants */
static const double PISTON_RHO 	= 920;
static const double PISTON_B 	= 1.19 * pow(10.0, 9) / scaleP;
static const double PISTON_C 	= 1.73;

/* STUB */
extern int Qr;
extern int m;
extern int beta;
#endif	/* GLOB_H */


