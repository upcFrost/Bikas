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
#include <math.h>
#include <algorithm>
#include <iomanip>
#include "types.h"


extern int n, i, j, max_i, max_j, I_k, i_v, j_v, i_sn, j_sn, i_sn_0, P_f;
extern double delta, delta_0, dx, dr, dt, p0, lambda, kappa, k, m_sn, max_z,max_z_0, f;
extern double V0;
extern double dM0;
extern bool broken_dt;
extern std::vector <double> t;
extern std::vector <double> x_sn;
extern std::vector <double> U_sn;

const int axis_j = 0; 

/* Scaling */
static const double scaleD = pow(10,	-2);
static const double scaleT = pow(10,	-5);
static const double scaleM = pow(10,	-6);
static const double scaleV = pow(10,	3);
static const double scaleE = pow(10,	0);
static const double scaleF = pow(10,	2);
static const double scaleP = pow(10,	6);
static const double scaleR = pow(10,	0);
static const double scaleDM = pow(10,	1);
static const double scaleFF = pow(10,	6);
static const double scaleIK = pow(10,	1);
static const double scaleGasKAPPA = pow(10, -3);
static const double scaleGasMU = pow(10, -1);
static const double scaleGasCP = pow(10, 6);
static const double scaleGasALPHAK = pow(10, 0);

/* Gas constants */
static const double alpha_k = 0.001006 / scaleGasALPHAK;
static const double gasA = -0.666667;
static const double gasMu = 2.85*pow(10,-5) / scaleGasMU; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasLambda = pow(10,-3);
static const double gasKappa = 1.21 / scaleGasKAPPA; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasCp = 1650 / scaleGasCP; // Решение задачи по определению эффективности многокамерного дульного тормоза
static const double gasPr = gasMu * gasCp / gasKappa;
static const double gasB = gasKappa / gasPr;


/* STUB */
extern int Qr;
extern int m;
extern int beta;
#endif	/* GLOB_H */


