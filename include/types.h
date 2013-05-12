#ifndef TYPES_H
#define	TYPES_H

#include <vector>

typedef struct Line2D { 
	double xbegin; 
	double xend; 
	double ybegin; 
	double yend;  
} Line2D;

typedef struct LineAngle2D { 
	double sin_a; 
	double cos_a; 
	double sin_2a; 
	double cos_2a;  
} LineAngle2D;

typedef struct WeightPart { double weight; int i; int j; } WeightPart;
typedef struct WeightVector { std::vector < WeightPart > x; std::vector < WeightPart > y; } WeightVector;

/* Cell parameters */
struct gasCell {
    // Geometry parameters
    int type;
    double r_1;
    double r_2;
    double x_1;
    double x_2;
    // Gas parameters
    double bar_Vx[5];
    double Vx[5]; // Array of 5
    double bar_Vr[5];
    double Vr[5]; // Array of 5
    double P[5]; // Array of 5
    double dM[5]; // Array of 5
    double dMVx[5]; // Array of 5
    double dMVr[5]; // Array of 5
    double dME[5]; // Array of 5
    int D[5];
    double A[5]; // A[0] == f, A[0] = A_i-1/2, A[1] = A_i+1/2, A[2] = A_j-1/2, A[3] = A_j+1/2
    double rho;
    double e;
    double bar_e;
    double alpha;
    double z;
    double psi;
    double m;
    
    WeightVector weightVector;
    LineAngle2D angle;
    
    double bar_z;
    double bar_psi;
    double final_z;
    double final_psi;
};

typedef std::vector < std::vector < std::vector < gasCell > > > cell2d;
typedef struct Point2D { double x; double y; } Point2D;
typedef struct Int2D { int i; int j; } Int2D;
typedef struct Triangle2D { Point2D point[3]; } Triangle2D;
typedef struct Trapezoid2D { Point2D point[4]; } Trapezoid2D;


#endif	/* TYPES_H */
