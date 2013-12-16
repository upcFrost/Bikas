/*
 * gasCell.hpp
 *
 *  Created on: 16.12.2013
 *      Author: PBelyaev
 */

#ifndef GASCELL_HPP_
#define GASCELL_HPP_

#include "Line2D.hpp"

class GasCell {
public:
    // Geometry parameters
    int type;
    double r_1, r_2, x_1, x_2;

    // Gas parameters
    double bar_Vx, Vx;
    double bar_Vr, Vr;
    double P;
    double dM[5]; // Array of 5
    int D[5];
    double A[5]; // A[0] == f, A[0] = A_i-1/2, A[1] = A_i+1/2, A[2] = A_j-1/2, A[3] = A_j+1/2
    double rho;
    double e;
    double bar_e;
    double alpha;

    WeightVector weightVector;
    Line2D line;

    double bar_z, final_z;
    double bar_psi, final_psi;

    GasCell();
    void SetGeometry(int, double, double, double, double);
    void SetLine(Line2D);
    void SetBarVx(double);
    void SetVx(double);
    void SetBarVr(double);
    void SetVr(double);
    void SetP(double);
    void SetEmptyCell();
}


#endif /* GASCELL_HPP_ */
