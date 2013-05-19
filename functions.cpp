#include "functions.h"

using namespace std;

/* Alpha calculation */
double pre_alpha(gasCell cell, int i, int i_p, double delta_0, double delta) {
    if (i <= i_p)
		//return delta_0/delta; 
		return 1;
	else
		return 1;
}

/* Cell geometry parameter calculation
 * 
 * TODO: Fix full[0], [1] and [2] EVERYWHERE
 * 
 *  */
void pre_cell_geometry(double array[5], gasCell cell, int i, int j) {
    double full[5];
	full[0] = M_PI*(2*(j-axis_j)+1)*pow(dr,2)*dx;
	full[1] = M_PI*(2*(j-axis_j)+1)*pow(dr,2);
	full[2] = M_PI*(2*(j-axis_j)+1)*pow(dr,2);
	full[3] = 2*M_PI*(j-axis_j)*dr*dx;
	full[4] = 2*M_PI*(j-axis_j+1)*dr*dx;
    switch (cell.type) {
	case 0:
		array[0] = 1;
		array[1] = 1;
		array[2] = 1;
		array[3] = 1;
		array[4] = 1;
		break;

	case 1:
	    //~ array[0] = (2*M_PI * (j*dr + (cell.r_1*dr + cell.r_2*dr)/2) * (cell.r_1*dr + cell.r_2*dr)/2 * dx) / full[0] ;
	    //~ array[1] = (M_PI * (2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
	    //~ array[2] = (M_PI * (2*j*dr*cell.r_2*dr + pow(cell.r_2*dr, 2))) / full[2] ;
	    //~ array[3] = 1;
	    //~ array[4] = 0;
	    array[0] = ((j-1)*(cell.r_1+cell.r_2) + cell.r_1*cell.r_2 +
			pow(cell.r_2-cell.r_1, 2)/3) / (2*j-1);
		array[1] = cell.r_1 * (j-1+cell.r_1/2) / (j-0.5);
		array[2] = cell.r_2 * (j-1+cell.r_2/2) / (j-0.5);
		array[3] = 1;
		array[4] = 0;
	    printf("Cell at (%d, %d) with type 1; A[1] = %4.4f, A[2] = %4.4f, A[3] = %4.4f, A[4] = %4.4f\n", i,j,array[1],array[2],array[3],array[4]);
	    break;

	case 2:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI * (j*dr + (cell.r_1*dr + cell.r_2*dr)/2) * (cell.r_1*dr + cell.r_2*dr)/2 * dx) / full[0] ;
	    array[1] = (M_PI * (2*(j+1)*dr*cell.r_1*dr - pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = (M_PI * (2*(j+1)*dr*cell.r_2*dr - pow(cell.r_2*dr, 2))) / full[2] ;
	    array[3] = 0;
	    array[4] = 1;
	    break;

	case 3:
	    //~ array[0] = (2*M_PI * (j*dr + dr/2) * (cell.x_1*dx + cell.x_2*dx)/2 * dr) / full[0] ;
	    //~ array[1] = 1;
	    //~ array[2] = 0;
	    //~ array[3] = (2*M_PI*j*cell.x_1*dx*dr) / full[3] ;
	    //~ array[4] = (2*M_PI*(j+1)*cell.x_2*dx*dr) / full[4] ;
	    array[0] = (j*cell.x_2 + (j-1)*cell.x_1 - (cell.x_2 - cell.x_1)/3) / (2*j-1);
	    array[1] = 1;
	    array[2] = 0;
	    array[3] = cell.x_1;
	    array[4] = cell.x_2;
	    printf("Cell at (%d, %d) with type 3; A[1] = %4.4f, A[2] = %4.4f, A[3] = %4.4f, A[4] = %4.4f\n", i,j,array[1],array[2],array[3],array[4]);
	    break;

	case 4:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI * (j*dr + dr/2) * (cell.x_1*dx + cell.x_2*dx)/2 * dr) / full[0] ;
	    array[1] = 0;
	    array[2] = 1;
	    array[3] = (2*M_PI*j*(dx - cell.x_1*dx)*dr) / full[3] ;
	    array[4] = (2*M_PI*j*(dx - cell.x_2*dx)*dr) / full[4] ;
	    break;

	case 5:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI*(j*dr + cell.r_1*dr + (dr-cell.r_1*dr)/2) * (dr-cell.r_1*dr)*cell.x_2*dx/2) / full[0] ;
	    array[1] = (M_PI*(2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = 1;
	    array[3] = 1;
	    array[4] = (2*M_PI*j*(dx-cell.x_2*dx)*dr) / full[4] ;
	    break;

	case 6:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI*(j*dr + cell.r_1*dr/2) * cell.r_1*dr*cell.x_1*dx / 2) / full[0] ;
	    array[1] = (M_PI*(2*(j+1)*dr*cell.r_1*dr - pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = 1;
	    array[3] = (2*M_PI*j*(dx - cell.x_1*dx)*dr) / full[3] ;
	    array[4] = 1;
	    break;

	case 7:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI*(j*dr + cell.r_1*dr/2)*cell.r_1*dr * (dx - cell.x_1*dx)/2) / full[0] ;
	    array[1] = 1;
	    array[2] = (M_PI*(2*(j+1)*dr*cell.r_2*dr - pow(cell.r_2*dr,2))) / full[2] ;
	    array[3] = (2*M_PI*j*cell.x_1*dx*dr) / full[3] ;
	    array[4] = 1;
	    break;

	case 8:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI*(j*dr + cell.r_1*dr + (dr - cell.r_1*dr)/2) * (dr - cell.r_1*dr)*(dx - cell.x_2*dx)/2) / full[0] ;
	    array[1] = 1;
	    array[2] = (M_PI*(2*j*dr*cell.r_1*dr - pow(cell.r_1*dr,2))) / full[2] ;
	    array[3] = 1;
	    array[4] = (2*M_PI*(j+1)*cell.x_2*dx*dr) / full[4] ;
	    break;

	case 9:
	    array[0] = (2*M_PI*(j*dr + cell.r_1*dr + (dr-cell.r_1*dr)/2) * (dr-cell.r_1*dr)*cell.x_2*dx/2) / full[0] ;
	    array[1] = (2*M_PI*j*pow(dr,2)*dx - M_PI*(2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = 0;
	    array[3] = 0;
	    array[4] = (2*M_PI*j*pow(dr,2)*dx - 2*M_PI*j*(dx-cell.x_2*dx)*dr) / full[4] ;
	    break;

	case 10:
	    //~ array[0] = (2*M_PI*(j*dr + cell.r_1*dr/2) * cell.r_1*dr*cell.x_1*dx / 2) / full[0] ;
	    //~ array[1] = (M_PI * (2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
	    //~ array[2] = 0;
	    //~ array[3] = (2*M_PI*j*cell.x_1*dx*dr) / full[3] ;
	    //~ array[4] = 0;
	    array[0] = cell.x_1*cell.r_1*(j-1 + cell.r_1/3) / (2*j-1);
	    array[1] = cell.r_1 * (j-1 + cell.r_1/2) / (j-0.5);
	    array[2] = 0;
	    array[3] = cell.x_1;
	    array[4] = 0;
	    printf("Cell at (%d, %d) with type 10; A[1] = %4.4f, A[2] = %4.4f, A[3] = %4.4f, A[4] = %4.4f\n", i,j,array[1],array[2],array[3],array[4]);
	    break;

	case 11:
	    array[0] = (2*M_PI*(j*dr + cell.r_1*dr/2)*cell.r_1*dr * (dx - cell.x_1*dx)/2) / full[0] ;
	    array[1] = 0;
	    array[2] = (2*M_PI*j*pow(dr,2)*dx - M_PI*(2*(j+1)*dr*cell.r_2*dr - pow(cell.r_2*dr,2))) / full[2] ;
	    array[3] = (2*M_PI*j*pow(dr,2)*dx - 2*M_PI*j*cell.x_1*dx*dr) / full[3] ;
	    array[4] = 0;
	    break;

	case 12:
	    array[0] = (2*M_PI*(j*dr + cell.r_1*dr + (dr - cell.r_1*dr)/2) * (dr - cell.r_1*dr)*(dx - cell.x_2*dx)/2) / full[0] ;
	    array[1] = 0;
	    array[2] = (2*M_PI*j*pow(dr,2)*dx - M_PI*(2*j*dr*cell.r_1*dr - pow(cell.r_1*dr,2))) / full[2] ;
	    array[3] = 0;
	    array[4] = (2*M_PI*j*pow(dr,2)*dx - 2*M_PI*(j+1)*cell.x_2*dx*dr) / full[4] ;
	    break;

	case 13:
		array[0] = 1;
		array[1] = 1;
		array[2] = 1;
		array[3] = 1;
		array[4] = 0;
		break;

	case 14:
		array[0] = 1;
		array[1] = 1;
		array[2] = 1;
		array[3] = 1;
		array[4] = 1;
		break;

	case 15:
		array[0] = 1;
		array[1] = 0;
		array[2] = 1;
		array[3] = 1;
		array[4] = 0;
		break;

	case 16:
		array[0] = 1;
		array[1] = 0;
		array[2] = 1;
		array[3] = 0;
		array[4] = 1;
		break;

	case 17:
		array[0] = 1;
		array[1] = 0;
		array[2] = 1;
		array[3] = 1;
		array[4] = 1;
		break;

	case 18:
		array[0] = 0;
		array[1] = 0;
		array[2] = 0;
		array[3] = 0;
		array[4] = 0;
		break;

	case 19:
		array[0] = 1;
		array[1] = 1;
		array[2] = 0;
		array[3] = 1;
		array[4] = 1;
		break;

	case 20:
		array[0] = 1;
		array[1] = 1;
		array[2] = 0;
		array[3] = 1;
		array[4] = 0;
		break;

	case 21:
		array[0] = 1;
		array[1] = 1;
		array[2] = 0;
		array[3] = 0;
		array[4] = 1;
		break;

	case 22:
	    array[0] = ((j-1)*(cell.r_1+cell.r_2) + cell.r_1*cell.r_2 +
			pow(cell.r_2-cell.r_1, 2)/3) / (2*j-1);
		array[1] = cell.r_1 * (j-1+cell.r_1/2) / (j-0.5);
		array[2] = 0;
		array[3] = 1;
		array[4] = 0;


	default:
		break;
	}
}


/* Euler stage U_sn calculation */
double euler_Usn(double P_sn, double S, double F, double dt, double U_prev) {

	if (P_sn > P_f || U_prev > 0) {
		//return cell->Vx[0] + (S*cell->P[4] - F)*dt / m_sn;
		return U_prev + (S*P_sn - F)*dt / (m_sn);
	} else {
		return 0;
	}
}

/* Euler stage X_sn calculation */
double euler_Xsn(double x_prev, double U_new) {
	return x_prev + U_new*dt;
}

/* Euler stage projectile border calculation */
void euler_proj_broder(double array[5], int j, double Xsn, double dx, double dr) {

	double full[5];
	full[0] = M_PI*(2*(j-axis_j)+1)*dx*pow(dr,2);
	full[1] = M_PI*(2*(j-axis_j)+1)*pow(dr,2);
	full[2] = M_PI*(2*(j-axis_j)+1)*pow(dr,2);
	full[3] = 2*M_PI*(j-axis_j)*dr*dx;
	full[4] = 2*M_PI*(j-axis_j+1)*dr*dx;

	array[0] = 1 + M_PI*(2*(j-axis_j)+1)*pow(dr,2)*fmod(Xsn, dx) / full[0];
	array[1] = 1;
	array[2] = 0;
	array[3] = 1 + 2*M_PI * (j-axis_j)*dr * fmod(Xsn, dx) / full[3];
	array[4] = 1 + 2*M_PI * (j-axis_j+1)*dr * fmod(Xsn, dx) / full[4];

	if (j == max_j - 2) array[4] = 0;
	if (j == axis_j) array[3] = 0;
}

 /* Vx calculation on euler stage */
double euler_bar_Vx(cell2d& cell, int n, int i, int j,
		double dt, double dx, double dr, int var) {

	gasCell curCell = cell[n][i][j];
	double result = 0;

	unsigned NEEDED_COND = 0;
	NEEDED_COND += P_PAR;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += VX_PAR;
	NEEDED_COND += VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

    //~ double P_i12 = (P_i1 + P)/2 * (1 - (k-1)*(Vx_i1 - Vx)*dt/dx);
	//~ double P_i_12 = (P_i_1 + P)/2 * (1 - (k-1)*(Vx - Vx_i_1)*dt/dx);
	double P_i12 = (brd[P_POS].i1j + brd[P_POS].ij)/2;
	double P_i_12 = (brd[P_POS].i_1j + brd[P_POS].ij)/2;

	switch (var) {
	case FIRST_ORDER:
		/** First order **/
		result = curCell.Vx[0] -
			(P_i12 - P_i_12) / dx * fmax(curCell.A[1], curCell.A[2]) * dt
			/ (curCell.rho * curCell.A[0]);
		break;

	case SECOND_ORDER:
		/** Second order by x, only when [i-2] and [i+2] cells're available **/
		result = curCell.Vx[0] -
			(1.0/12.0*cell[n][i-2][j].P[0] - 2.0/3.0*brd[P_POS].i_1j +
					2.0/3.0*brd[P_POS].i1j - 1.0/12.0*cell[n][i+2][j].P[0])
			/ dx * fmax(curCell.A[1],curCell.A[2]) * dt / (curCell.rho * curCell.A[0]);
		break;

	case FIRST_ORDER_NS:
		/** With Navier-Stocks **/
		result = curCell.Vx[0] - dt  / (curCell.rho * curCell.A[0] * dx) *
			(
				// Pressure
				(P_i12 - P_i_12) * fmax(curCell.A[1], curCell.A[2]) -
				// Viscosity
				gasMu * dx * (
					(gasA+2)/pow(dx,2) * (
						brd[VX_POS].i1j - 2*curCell.Vx[0] + brd[VX_POS].i_1j
					) + (
						brd[VX_POS].ij1 - 2*curCell.Vx[0] + brd[VX_POS].ij_1
					) / pow(dr,2) +
					(gasA+1)/(4*dx*dr) * (
						brd[VR_POS].i1j1 + brd[VR_POS].i_1j_1 - brd[VR_POS].i_1j1 - brd[VR_POS].i1j_1
					)
				)
			);
		break;

	case SECOND_ORDER_NS:
		/** Second order Navier-Stocks, only when [i-2] and [i+2] cells're available **/
		result = curCell.Vx[0] - dt / (curCell.rho * curCell.A[0] * dx) *
			(
				// Pressure
				(1.0/12.0*cell[n][i-2][j].P[0] - 2.0/3.0*brd[P_POS].i_1j +
						2.0/3.0*brd[P_POS].i1j - 1.0/12.0*cell[n][i+2][j].P[0])
						* fmax(curCell.A[1],curCell.A[2]) -
				// Viscosity
				gasMu * dx * (
					(gasA+2)/pow(dx,2) * (
						brd[VX_POS].i1j - 2*curCell.Vx[0] + brd[VX_POS].i_1j
					) + (
						brd[VX_POS].ij1 - 2*curCell.Vx[0] + brd[VX_POS].ij_1
					) / pow(dr,2) +
					(gasA+1)/(4*dx*dr) * (
						brd[VR_POS].i1j1 + brd[VR_POS].i_1j_1 - brd[VR_POS].i_1j1 - brd[VR_POS].i1j_1
					)
				)
			);

		break;

	default:
		break;
	}

	if (i == 107 && j == 15) {
		gasCell cell_106 = cell[n][106][15];
		gasCell cell_107 = cell[n][107][15];
		gasCell cell_108 = cell[n][108][15];
		printf("123");
	}

    if (fabs(result-brd[VX_POS].ij) < pow(10.0,-15)) result = brd[VX_POS].ij;
	return result;
}

/* Vr calculation on euler stage */
double euler_bar_Vr(cell2d& cell, int n, int i, int j,
		double dt, double dx, double dr, int var) {

	gasCell curCell = cell[n][i][j];
	double result = 0;

	unsigned NEEDED_COND = 0;
	NEEDED_COND += P_PAR;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += VX_PAR;
	NEEDED_COND += VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);
	BorderCond borderArray[8];
	borderArray[0] = brd[0];
	borderArray[1] = brd[1];
	borderArray[2] = brd[2];
	borderArray[3] = brd[3];
	borderArray[4] = brd[4];
	borderArray[5] = brd[5];
	borderArray[6] = brd[6];
	borderArray[7] = brd[7];

    //~ double P_j12 = (P_j1 + P)/2 * (1 - (k-1)*(Vr_j1 - Vr)*dt/dr);
	//~ double P_j_12 = (P_j_1 + P)/2 * (1 - (k-1)*(Vr - Vr_j_1)*dt/dr);
	double P_j12 = (brd[P_POS].ij1 + brd[P_POS].ij)/2;
	double P_j_12 = (brd[P_POS].ij_1 + brd[P_POS].ij)/2;

	switch (var) {
	case FIRST_ORDER:
		/** First order **/
		result = curCell.Vr[0] - (P_j12 - P_j_12)
			/ brd[RHO_POS].ij * dt / curCell.A[0] *
			fmax(curCell.A[4],curCell.A[3]) / dr;
		break;

	case SECOND_ORDER:
		/** Second order by r, only when [j-2] and [j+2] cells're available **/
		result = curCell.Vr[0] -
			(1.0/12.0*cell[n][i][j-2].P[0] - 2.0/3.0*brd[P_POS].ij_1
					+ 2.0/3.0*brd[P_POS].ij1 - 1.0/12.0*cell[n][i][j+2].P[0])
			/ curCell.rho * dt / curCell.A[0] *	fmax(curCell.A[4],curCell.A[3]) / dr;
		break;

	case FIRST_ORDER_NS:
		/** With Navier-Stocks **/
		result = curCell.Vr[0] - dt  / (curCell.rho * curCell.A[0] * dr) *
			(
			// Pressure
			(P_j12 - P_j_12)*fmax(curCell.A[3],curCell.A[4]) -
			// Viscosity
			gasMu * dr *(
				(gasA+2)/pow(dr,2) * (
					brd[VR_POS].ij1 - 2*curCell.Vr[0] + brd[VR_POS].ij_1
				) + (
					brd[VR_POS].i1j - 2*curCell.Vr[0] + brd[VR_POS].i_1j
				) / pow(dx,2) +
				(gasA+1)/(4*dx*dr) * (
					brd[VX_POS].i1j1 + brd[VX_POS].i_1j_1 - brd[VX_POS].i_1j1 - brd[VX_POS].i1j_1
				)
			)
		) ;
		break;

	case SECOND_ORDER_NS:
		result = curCell.Vr[0] - dt  / (curCell.rho * curCell.A[0] * dr) *
			(
				// Pressure
				(1.0/12.0*cell[n][i][j-2].P[0] - 2.0/3.0*brd[P_POS].ij_1
					+ 2.0/3.0*brd[P_POS].ij1 - 1.0/12.0*cell[n][i][j+2].P[0])
					* fmax(curCell.A[4],curCell.A[3]) -
				// Viscosity
				gasMu * dr *(
					(gasA+2)/pow(dr,2) * (
						brd[VR_POS].ij1 - 2*curCell.Vr[0] + brd[VR_POS].ij_1
					) + (
						brd[VR_POS].i1j - 2*curCell.Vr[0] + brd[VR_POS].i_1j
					) / pow(dx,2) +
					(gasA+1)/(4*dx*dr) * (
						brd[VX_POS].i1j1 + brd[VX_POS].i_1j_1 - brd[VX_POS].i_1j1 - brd[VX_POS].i1j_1
					)
				)
			);
		break;

	default:
		break;
	}

	/** From my calc **/
    //~ double result = cell.Vr[0] -
        //~ cell.alpha *
        //~ (
            //~ cell.P[4] * cell.A[4] -
            //~ cell.P[3] * cell.A[3]
        //~ )
        //~ * dt / (cell.rho * cell.A[0] * dr)
        //~ -
        //~ cell.alpha * cell.P[0] * dt / cell.rho;

    /** Cleaning noise **/
    //~ if (fabs(curCell.P[4] - curCell.P[3]) < pow(10,-5)/scaleV) result = curCell.Vr[0];

	if (fabs(result-curCell.Vr[0]) < pow(10.0,-15)) result = curCell.Vr[0];
	return result;
}

void rotateVectors(double& Vx, double& Vr, LineAngle2D angle) {
	double V = sqrt(pow(Vx,2)+pow(Vr,2));
	Vx = V*angle.cos_a;
	Vr = V*angle.sin_a;
}

double * smoothSpeed(double * Vx, double * Vr, LineAngle2D angle) {
	double V[3];
	double * result = new double[2];

	for (int i = 0; i < 3; i++)
		V[i] = sqrt(pow(Vx[i],2)+pow(Vr[i],2));
	V[1] = 0.1*V[0] + 0.8*V[1] + 0.1*V[2];
	result[0] = V[1]*angle.cos_a;
	result[1] = V[1]*angle.sin_a;

	return result;
}

/* Energy calculation on euler stage */
double euler_bar_e(cell2d& cell, int n, int i, int j,
		double dt, double dx, double dr, int var) {

	gasCell curCell = cell[n][i][j];
	double result = 0;

	unsigned NEEDED_COND = 0;
	NEEDED_COND += P_PAR;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += E_PAR;
	NEEDED_COND += VX_PAR;
	NEEDED_COND += VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

	double P_i12 = (brd[P_POS].i1j + brd[P_POS].ij)/2;
	double P_i_12 = (brd[P_POS].i_1j + brd[P_POS].ij)/2;
	double P_j12 = (brd[P_POS].ij1 + brd[P_POS].ij)/2;
	double P_j_12 = (brd[P_POS].ij_1 + brd[P_POS].ij)/2;
	double Vx_i12 = (brd[VX_POS].i1j + brd[VX_POS].ij)/2;
	double Vx_i_12 = (brd[VX_POS].i_1j + brd[VX_POS].ij)/2;
	double Vr_j12 = (brd[VR_POS].ij1 + brd[VR_POS].ij)/2;
	double Vr_j_12 = (brd[VR_POS].ij_1 + brd[VR_POS].ij)/2;

	switch (var) {
	case FIRST_ORDER:
	    /** First order **/
		result = curCell.e -
			(
				(P_i12*Vx_i12 - P_i_12*Vx_i_12) / dx * fmax(curCell.A[1],curCell.A[2]) +
				(fabs(j-axis_j)*P_j12*Vr_j12 - fabs(j-axis_j-1)*P_j_12*Vr_j_12) / ((fabs(j-axis_j-0.5))*pow(dr,2)) * fmax(curCell.A[3],curCell.A[4])
			) * dt / (curCell.A[0] * curCell.rho)
		/** If powder present **/
			+ curCell.P[0]/((k-1)*I_k) * (1 + lambda * cell[n][i][j].final_z) * f*kappa*dt
		/** End here **/
			;
		break;

	case SECOND_ORDER:
		break;

	case FIRST_ORDER_NS: {
	    /** With Navier-Stocks **/
	    double I1 = curCell.e - (pow(curCell.Vx[0],2) + pow(curCell.Vr[0],2))/2;
	    double I2 = brd[E_POS].i_1j - (pow(brd[VX_POS].i_1j,2) + pow(brd[VR_POS].i_1j,2))/2;
	    double I3 = brd[E_POS].i1j - (pow(brd[VX_POS].i1j,2) + pow(brd[VR_POS].i1j,2))/2;
	    double I4 = brd[E_POS].ij_1 - (pow(brd[VX_POS].ij_1,2) + pow(brd[VR_POS].ij_1,2))/2;
	    double I5 = brd[E_POS].ij1 - (pow(brd[VX_POS].ij1,2) + pow(brd[VR_POS].ij1,2))/2;
	    result = curCell.e - dt / (curCell.A[0] * curCell.rho) *
			(
				// Pressure
				(P_i12*Vx_i12 - P_i_12*Vx_i_12) / dx * fmax(curCell.A[1],curCell.A[2]) +
				(fabs(j-axis_j)*P_j12*Vr_j12 - fabs(j-axis_j-1)*P_j_12*Vr_j_12) /
					((fabs(j-axis_j-0.5))*pow(dr,2)) * fmax(curCell.A[3],curCell.A[4])
				-
				// Viscosity
				// X axis
				gasMu * (gasA + 2) / (2 * pow(dx,2)) * (
					pow(brd[VX_POS].i1j, 2) + pow(brd[VX_POS].i_1j, 2) - 2*pow(curCell.Vx[0], 2)
				)
				+
				gasA * gasMu / (4*dx*dr) * (
					brd[VX_POS].i1j * (brd[VR_POS].i1j1 - brd[VR_POS].i1j_1) -
					brd[VX_POS].i_1j * (brd[VR_POS].i_1j1 - brd[VR_POS].i_1j_1)
				)
				+
				gasMu / (pow(dx,2)) * (
					pow(brd[VR_POS].i1j, 2) + pow(brd[VR_POS].i_1j, 2) - 2*pow(curCell.Vr[0], 2)
				)
				+
				gasMu / (4*dx*dr) * (
					brd[VR_POS].i1j * (brd[VX_POS].i1j1 - brd[VX_POS].i1j_1) -
					brd[VR_POS].i_1j * (brd[VX_POS].i_1j1 - brd[VX_POS].i_1j_1)
				)
				+
				gasB * gasMu / pow(dx,2) * (I3 + I2 - 2*I1)
				+
				// R axis
				gasMu * (gasA + 2) / (2 * pow(dr,2)) * (
					pow(brd[VR_POS].ij1, 2) + pow(brd[VR_POS].ij_1, 2) - 2*pow(curCell.Vr[0], 2)
				)
				+
				gasA * gasMu / (4*dx*dr) * (
					brd[VR_POS].ij1 * (brd[VX_POS].i1j1 - brd[VX_POS].i_1j1) -
					brd[VR_POS].ij_1 * (brd[VX_POS].i1j_1 - brd[VX_POS].i_1j_1)
				)
				+
				gasMu / (pow(dr,2)) * (
					pow(brd[VX_POS].ij1, 2) + pow(brd[VX_POS].ij_1, 2) - 2*pow(curCell.Vx[0], 2)
				)
				+
				gasMu / (4*dx*dr) * (
					brd[VX_POS].ij1 * (brd[VR_POS].i1j1 - brd[VR_POS].i_1j1) -
					brd[VX_POS].ij_1 * (brd[VR_POS].i1j_1 - brd[VR_POS].i_1j_1)
				)
				+
				gasB * gasMu / pow(dr,2) * (I5 + I4 - 2*I1)
		)
		/** If powder present **/
		+ curCell.P[0]/((k-1)*I_k) * (1 + lambda*curCell.final_z) * f*kappa*dt
		;
	    break;
	}

	default:
		break;
	}

	/** From my calc **/
    //~ double result = cell->e -
        //~ cell->A[0] * cell->P[0] *
        //~ (
            //~ cell->Vx[2] - cell->Vx[1]
        //~ ) * dt /
        //~ (
            //~ cell->A[0] * cell->rho * dx
        //~ ) -
        //~ cell->A[0] * cell->P[0] *
        //~ (
            //~ cell->Vr[4] - cell->Vr[3]
        //~ ) * dt /
        //~ (
            //~ cell->A[0] * cell->rho * dr
        //~ ) -
        //~ cell->alpha * cell->Vx[0] *
        //~ (
            //~ cell->P[2] * cell->A[2] - cell->P[1] * cell->A[1]
        //~ ) * dt /
        //~ (
            //~ cell->A[0] * cell->rho * dx
        //~ ) -
        //~ cell->alpha * cell->Vr[0] *
        //~ (
            //~ cell->P[4]*cell->A[4] - cell->P[3]*cell->A[3]
        //~ ) * dt /
        //~ (
            //~ cell->A[0] * cell->rho * dr
        //~ ) -
        //~ cell->A[0] * cell->P[0] * cell->Vr[0] * dt /
        //~ (
            //~ cell->A[0] * cell->rho * j * dr
        //~ ) +
            //~ Qr * dt /
        //~ (
            //~ cell->A[0] * cell->rho
        //~ );

	/** Catch less-than-zero rezult **/
	//~ if (result < 0) {
		//~ cout << "bar_e < 0" << endl;
		//~ cout << "e = " << (*cell)[n][i][j].e << endl <<
				//~ "P[1-4] = {" << (*cell)[n][i][j].P[1] << "    " << (*cell)[n][i][j].P[2] << "    " << (*cell)[n][i][j].P[3] << "    " << (*cell)[n][i][j].P[4] << "}" << endl <<
				//~ "V[1-4] = {" << (*cell)[n][i][j].Vx[1] << "    " << (*cell)[n][i][j].Vx[2] << "    " << (*cell)[n][i][j].Vr[3] << "    " << (*cell)[n][i][j].Vr[4] << "}" << endl <<
				//~ "A[1-4] = {" << (*cell)[n][i][j].A[1] << "    " << (*cell)[n][i][j].A[2] << "    " << (*cell)[n][i][j].A[3] << "    " << (*cell)[n][i][j].A[4] << "}" << endl;
		//~ cout << "first brackets: " << ((*cell)[n][i][j].P[2]*(*cell)[n][i][j].Vx[2] - (*cell)[n][i][j].P[1]*(*cell)[n][i][j].Vx[1]) / dx *
		//~ fmax((*cell)[n][i][j].A[1],(*cell)[n][i][j].A[2]) << endl;
		//~ cout << "second brackets: " << ((*cell)[n][i][j].P[4]*(*cell)[n][i][j].Vr[4] - (*cell)[n][i][j].P[3]*(*cell)[n][i][j].Vr[3]) / ((fabs(j-axis_j-0.5))*dr) *
		//~ fmax((*cell)[n][i][j].A[3],(*cell)[n][i][j].A[4]) << endl;
		//~ getchar();
	//~ }

	//~ if (result < 0 || result > pow(10,8)) {
		//~ cout << "bar_E = " << result << endl;
		//~ cout << "i = " << i << ", j = " << j << endl;
		//~ cout << "Pressure part = " << 1/(2*dx) * (
				//~ (*cell)[n][i][j].Vx[0] * (
					//~ (*cell)[n][i-1][j].P[0] - (*cell)[n][i+1][j].P[0]
				//~ ) + (*cell)[n][i][j].P[0] * (
					//~ (*cell)[n][i-1][j].Vx[0] - (*cell)[n][i+1][j].Vx[0]
				//~ )
			//~ ) + 1/(2*dr) * (
				//~ (*cell)[n][i][j].Vr[0] * (
					//~ (*cell)[n][i][j-1].P[0] - (*cell)[n][i][j+1].P[0]
				//~ ) + (*cell)[n][i][j].P[0] * (
					//~ (*cell)[n][i][j-1].Vr[0] - (*cell)[n][i][j+1].Vr[0]
				//~ )
			//~ ) << endl;
		//~ cout << "Viscosity part X axis first = " << gasMu * (gasA + 2) / (2 * pow(dx,2)) * (
				//~ pow((*cell)[n][i+1][j].Vx[0], 2) +
				//~ pow((*cell)[n][i-1][j].Vx[0], 2) -
				//~ 2*pow((*cell)[n][i][j].Vx[0], 2)
			//~ ) << endl;
		//~ cout << "Viscosity part X axis second = " << gasA * gasMu / (4*dx*dr) * (
				//~ (*cell)[n][i+1][j].Vx[0] * (
					//~ (*cell)[n][i+1][j+1].Vr[0] - (*cell)[n][i+1][j-1].Vr[0]
				//~ ) - (*cell)[n][i-1][j].Vx[0] * (
					//~ (*cell)[n][i-1][j+1].Vr[0] - (*cell)[n][i-1][j-1].Vr[0]
				//~ )
			//~ ) << endl;
		//~ cout << "Viscosity part X axis third = " << gasMu / (pow(dx,2)) * (
				//~ pow((*cell)[n][i+1][j].Vr[0], 2) +
				//~ pow((*cell)[n][i-1][j].Vr[0], 2) -
				//~ 2*pow((*cell)[n][i][j].Vr[0], 2)
			//~ ) << endl;
		//~ cout << "Viscosity part X axis forth = " << gasMu / (4*dx*dr) * (
				//~ (*cell)[n][i+1][j].Vr[0] * (
					//~ (*cell)[n][i+1][j+1].Vx[0] - (*cell)[n][i+1][j-1].Vx[0]
				//~ ) - (*cell)[n][i-1][j].Vr[0] * (
					//~ (*cell)[n][i-1][j+1].Vx[0] - (*cell)[n][i-1][j-1].Vx[0]
				//~ )
			//~ ) << endl;
		//~ cout << "Viscosity part X axis sixth = " << gasB * gasMu / pow(dx,2) * (
				//~ I3 + I2 - 2*I1
			//~ ) << endl;
		//~ cout << "Viscosity part X axis sixth in brackets = " << (I3 + I2 - 2*I1) << endl;
		//~ cout << "Viscosity part X axis sixth in brackets first = " << I3 << endl;
		//~ cout << "Viscosity part X axis sixth in brackets second = " << I2 << endl;
		//~ cout << "e[i-1] = " << (*cell)[n][i-1][j].e << endl;
		//~ cout << "(Vx[0]^2 + Vr[0]^2)/2 = " << (pow((*cell)[n][i-1][j].Vx[0],2) + pow((*cell)[n][i-1][j].Vr[0],2))/2 << endl;
		//~ cout << "Viscosity part X axis sixth in brackets third = " << 2*I1 << endl;
			//~ // R axis
			//~ cout << "Viscosity part R axis = " << gasMu * (gasA + 2) / (2 * pow(dr,2)) * (
				//~ pow((*cell)[n][i][j+1].Vr[0], 2) +
				//~ pow((*cell)[n][i][j-1].Vr[0], 2) -
				//~ 2*pow((*cell)[n][i][j].Vr[0], 2)
			//~ ) + gasA * gasMu / (4*dx*dr) * (
				//~ (*cell)[n][i][j+1].Vr[0] * (
					//~ (*cell)[n][i+1][j+1].Vx[0] - (*cell)[n][i-1][j+1].Vx[0]
				//~ ) - (*cell)[n][i][j-1].Vr[0] * (
					//~ (*cell)[n][i+1][j-1].Vx[0] - (*cell)[n][i-1][j-1].Vx[0]
				//~ )
			//~ ) + gasMu / (pow(dr,2)) * (
				//~ pow((*cell)[n][i][j+1].Vx[0], 2) +
				//~ pow((*cell)[n][i][j-1].Vx[0], 2) -
				//~ 2*pow((*cell)[n][i][j].Vx[0], 2)
			//~ ) + gasMu / (4*dx*dr) * (
				//~ (*cell)[n][i][j+1].Vx[0] * (
					//~ (*cell)[n][i+1][j+1].Vr[0] - (*cell)[n][i-1][j+1].Vr[0]
				//~ ) - (*cell)[n][i][j-1].Vx[0] * (
					//~ (*cell)[n][i+1][j-1].Vr[0] - (*cell)[n][i-1][j-1].Vr[0]
				//~ )
			//~ ) + gasB * gasMu / pow(dr,2) * (
				//~ I5 + I4 - 2*I1
			//~ ) << endl;
		//~ getchar();
	//~ }

	if (fabs(result-brd[E_POS].ij) < pow(10.0,-15)) result = brd[E_POS].ij;
	if (result == result) { // if not NaN
		return result;
	} else {
		cout << "bar_e is NaN" << endl
			<< "i = " << i << endl
			<< "j = " << j << endl
			<< "Vr[3] = " << cell[n][i][j].Vr[3] << endl
			<< "Vr[4] = " << cell[n][i][j].Vr[4] << endl
			<< "P[3] = " << cell[n][i][j].P[3] << endl
			<< "P[4] = " << cell[n][i][j].P[4] << endl
			<< "Second brackets = " <<((j-axis_j)*cell[n][i][j].P[4]*cell[n][i][j].Vr[4] - (j-axis_j-1)*cell[n][i][j].P[3]*cell[n][i][j].Vr[3]) / ((fabs(j-axis_j-0.5))*dr) *
				fmax(cell[n][i][j].A[3],cell[n][i][j].A[4]) << endl
			<< "(fabs(j-axis_j-0.5))*dr = " << (fabs(j-axis_j-0.5))*dr << endl
			<< "fmax(cell[n][i][j].A[3],cell[n][i][j].A[4]) = " << fmax(cell[n][i][j].A[3],cell[n][i][j].A[4]) << endl << endl;
		return 0;
	}
}

double lagrange_e(gasCell * prevCell, gasCell * cell) {
	double Ak = prevCell->A[0]/cell->A[0] < 1 ? prevCell->A[0]/cell->A[0] : 1;
    //~ double scaled_e = (cell->bar_e - (pow(prevCell->Vx[0],2) + pow(prevCell->Vr[0],2))/2)*Ak +
		//~ (pow(prevCell->Vx[0],2) + pow(prevCell->Vr[0],2))/2;
	double scaled_e = cell->bar_e * Ak;
	return scaled_e;
}

/* Powder -> gas */
double lagrange_m(gasCell * cell) {
	if (cell->psi < 1) {
		return (1 - cell->alpha)*delta / (1 - cell->psi) * (kappa + 2 * kappa * lambda * cell->z) * cell->P[0]/I_k;
	} else {
		return 0;
	}
}

/* Density on lagrange stage */
double lagrange_rho(gasCell * cell, gasCell * prevCell, int i, int j, double dt, double dx, double dr) {
    double result;
    if (fabs(cell->dM[1]) > pow(10.0,-15) ||
    		fabs(cell->dM[2]) > pow(10.0,-15) ||
    		fabs(cell->dM[3]) > pow(10.0,-15) ||
    		fabs(cell->dM[4]) > pow(10.0,-15)) {
	result = cell->rho +
		(
            cell->dM[1] +
            cell->dM[3] -
            cell->dM[2] -
            cell->dM[4]
        ) /
        (cell->A[0] * dx * fabs(j-axis_j-0.5) * pow(dr,2))
        ;
    } else {
    	result = cell->rho;
    }
        //~ if (j == 6 && prevCell->A[0]/cell->A[0] > 1) {
			//~ cout << "K > 1" << endl
				//~ << "j = " << j << endl
				//~ << "rho = " << cell->rho << endl
				//~ << "result = " << result << endl
				//~ << "k = " << k << endl
				//~ << "prev A[0] = " << prevCell->A[0] << endl
				//~ << "A[0] = " << cell->A[0] << endl
				//~ << "dM[1] = " << cell->dM[1] << endl
				//~ << "dM[2] = " << cell->dM[2] << endl
				//~ << "dM[3] = " << cell->dM[3] << endl
				//~ << "dM[4] = " << cell->dM[4] << endl;
			//~ getchar();
		//~ }
        if (result < 0) {
			cout << "rho < 0" << endl
				<< "i = " << i << endl
				<< "j = " << j << endl
				<< "rho = " << cell->rho << endl
				<< "A[0] = " << cell->A[0] << endl
				<< "dM[1] = " << cell->dM[1]*dt << endl
				<< "dM[2] = " << cell->dM[2]*dt << endl
				<< "dM[3] = " << cell->dM[3]*dt << endl
				<< "dM[4] = " << cell->dM[4]*dt << endl
				<< "dx = " << dx << endl
				<< "dr = " << dr << endl
				<< "dt = " << dt << endl << endl;
			getchar();
		}
    return result;
}

/* Lagrangian stage mass calculation */
void lagrange_mass(double array[21], cell2d& cell, int i, int j, int n,
        double dx, double dr, double dt) {

    /****************************************************************
     * Additional zero to have the same length as A					*
     * first empty + dM[4] + dMVx[4] + dMVr[4] + dME[4] + D[4]		*
     * **************************************************************/
    for (int iter = 0; iter < 21; iter++) array[iter] = 0;

	gasCell curCell = cell[n][i][j];

	unsigned NEEDED_COND = 0;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += BAR_VX_PAR;
	NEEDED_COND += BAR_VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

    double Vx_i12 = (brd[BAR_VX_POS].i1j + brd[BAR_VX_POS].ij)/2;
	double Vx_i_12 = (brd[BAR_VX_POS].i_1j + brd[BAR_VX_POS].ij)/2;
	double Vr_j12 = (brd[BAR_VR_POS].ij1 + brd[BAR_VR_POS].ij)/2;
	double Vr_j_12 = (brd[BAR_VR_POS].ij_1 + brd[BAR_VR_POS].ij)/2;

	/** First order **/
	bool ruleVx1 = Vx_i_12 > 0 ? true : false;
	bool ruleVx2 = Vx_i12 > 0 ? true : false;
	bool ruleVr1 = Vr_j_12 > 0 ? true : false;
	bool ruleVr2 = Vr_j12 > 0 ? true : false;

	if (ruleVx1) {
		array[1] = (fabs(j-axis_j-0.5))*cell[n][i-1][j].A[2] * brd[RHO_POS].i_1j * Vx_i_12 * pow(dr,2) * dt;
	} else {
		array[1] = (fabs(j-axis_j-0.5))*curCell.A[1] * brd[RHO_POS].ij * Vx_i_12 * pow(dr,2) * dt;
	}
	if (fabs(array[1]) < pow(10,-14)) array[1] = 0;
	if (ruleVx2) {
		array[2] = (fabs(j-axis_j-0.5))*curCell.A[2] * brd[RHO_POS].ij * Vx_i12 * pow(dr,2) * dt;
	} else {
		array[2] = (fabs(j-axis_j-0.5))*cell[n][i+1][j].A[1] * brd[RHO_POS].i1j * Vx_i12 * pow(dr,2) * dt;
	}
	if (fabs(array[2]) < pow(10,-14)) array[1] = 0;
	if (ruleVr1) {
		array[3] = cell[n][i][j-1].A[4] * brd[RHO_POS].ij_1 * Vr_j_12 * j*dx*dr * dt;
	} else {
		array[3] = curCell.A[3] * brd[RHO_POS].ij * Vr_j_12 * j*dx*dr * dt;
	}
	if (fabs(array[3]) < pow(10,-14)) array[1] = 0;
	if (ruleVr2) {
		array[4] = curCell.A[4] * brd[RHO_POS].ij * Vr_j12 * (j+1)*dx*dr * dt;
	} else {
		array[4] = cell[n][i][j+1].A[3] * brd[RHO_POS].ij1 * Vr_j12 * (j+1)*dx*dr * dt;
	}
	if (fabs(array[4]) < pow(10,-14)) array[1] = 0;
	array[5] = ruleVx1 ? 1 : 0;
	array[6] = ruleVx2 ? 0 : 1;
	array[7] = ruleVr1 ? 1 : 0;
	array[8] = ruleVr2 ? 0 : 1;

	/** Central **/
//	array[1] = brd[RHO_POS].i_1j * brd[BAR_VX_POS].i_1j * (fabs(j-axis_j-0.5)) * pow(dr,2) * dt;
//	array[2] = brd[RHO_POS].i1j * brd[BAR_VX_POS].i1j * (fabs(j-axis_j-0.5)) * pow(dr,2) * dt;
//	array[3] = brd[RHO_POS].ij_1 * brd[BAR_VR_POS].ij_1 * j*dx*dr * dt;
//	array[4] = brd[RHO_POS].ij1 * brd[BAR_VR_POS].ij1 * (j+1)*dx*dr * dt;
//	array[5] = ruleVx1 ? 1 : 0;
//	array[6] = ruleVx2 ? 0 : 1;
//	array[7] = ruleVr1 ? 1 : 0;
//	array[8] = ruleVr2 ? 0 : 1;

	/** Second order **/
    //~ bool ruleVx11 = Vx + Vx_i_1 > 0 ? true : false;
    //~ bool ruleVx12 = Vx - (Vx_i1 - Vx_i_1)/4 > 0 ? true : false;
    //~ bool ruleVx13 = Vx_i_1 + (Vx - Vx_i_2)/4 > 0 ? true : false;
    //~ bool ruleVx21 = Vx + Vx_i1 > 0 ? true : false;
    //~ bool ruleVx22 = Vx + (Vx_i1 - Vx_i_1)/4 > 0 ? true : false;
    //~ bool ruleVx23 = Vx_i1 - (Vx_i2 - Vx)/4 > 0 ? true : false;
    //~ bool ruleVr11 = Vr + Vr_j_1 > 0 ? true : false;
    //~ bool ruleVr12 = Vr + (Vr_j1 - Vr_j_1)/4 > 0 ? true : false;
    //~ bool ruleVr13 = Vr_j_1 - (Vr_j_2 - Vr)/4 > 0 ? true : false;
    //~ bool ruleVr21 = Vr + Vr_j1 > 0 ? true : false;
    //~ bool ruleVr22 = Vr + (Vr_j1 - Vr_j_1)/4 > 0 ? true : false;
    //~ bool ruleVr23 = Vr_j1 - (Vr_j2 - Vr)/4 > 0 ? true : false;
    //~
    //~ // dM(i-1/2,j)
    //~ // From right to left
    //~ if (ruleVx11 == ruleVx12 && !ruleVx11 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i1) > pow(10,-15) || fabs(Vx_i_1) > pow(10,-15))) {
        //~ array[1] = (Vx - (Vx_i1 - Vx_i_1)/4) *
            //~ (rho - (rho_i1 - rho_i_1)/4) * cell[n][i-1][j].A[2] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ }
    //~ // From left to right
    //~ if (ruleVx11 == ruleVx13 && ruleVx11 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i_2) > pow(10,-15) || fabs(Vx_i_1) > pow(10,-15))) {
	//~ array[1] = (Vx_i_1 + (Vx - Vx_i_2)/4) *
            //~ (rho_i_1 + (rho - rho_i_2)/4) * cell[n][i][j].A[1] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ }
    //~ // dM(i+1/2,j)
    //~ // From left to right
    //~ if (ruleVx21 == ruleVx22 && ruleVx21 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i1) > pow(10,-15) || fabs(Vx_i_1) > pow(10,-15))) {
        //~ array[2] = (Vx + (Vx_i1 - Vx_i_1)/4) *
            //~ (rho + (rho_i1 - rho_i_1)/4) * cell[n][i][j].A[2] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ }
    //~ // From right to left
    //~ if (ruleVx21 == ruleVx23 && !ruleVx21 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i1) > pow(10,-15) || fabs(Vx_i2) > pow(10,-15))) {
        //~ array[2] = (Vx_i1 - (Vx_i2 - Vx)/4) *
            //~ (rho_i1 - (rho_i2 - rho)/4) * cell[n][i+1][j].A[1] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ }
    //~ // dM(i,j-1/2)
    //~ if (ruleVr11 == ruleVr13 && ruleVr11 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j_1) > pow(10,-15) || fabs(Vr_j_2) > pow(10,-15))) {
        //~ array[3] = (Vr_j_1 - (Vr_j_2 - Vr)/4) *
            //~ (rho_j_1 - (rho_j_2 - rho)/4) * cell[n][i][j-1].A[4] * dx * dt;
    //~ }
    //~ if (ruleVr11 == ruleVr12 && !ruleVr11 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j1) > pow(10,-15) || fabs(Vr_j_1) > pow(10,-15))) {
	//~ array[3] = (Vr + (Vr_j_1 - Vr_j1)/4) *
            //~ (rho + (rho_j_1 - rho_j1)/4) * cell[n][i][j].A[3] * dx * dt;
    //~ }
    //~ // dM(i,j+1/2)
    //~ if (ruleVr21 == ruleVr22 && ruleVr21 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j1) > pow(10,-15) || fabs(Vr_j_1) > pow(10,-15))) {
            //~ array[4] = (Vr + (Vr_j1 - Vr_j_1)/4) *
                    //~ (rho + (rho_j1 - rho_j_1)/4) * cell[n][i][j].A[4] * dx * dt;
	//~ }
    //~ if (ruleVr21 == ruleVr23 && !ruleVr21 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j1) > pow(10,-15) || fabs(Vr_j2) > pow(10,-15))) {
            //~ array[4] = (Vr_j1 - (Vr_j2 - Vr)/4) *
                    //~ (rho_j1 - (rho_j2 - rho)/4) * cell[n][i][j+1].A[3] * dx * dt;
    //~ }
    //~
    //~ if (cell[n][i][j].A[1] == 0) array[1] = 0;
	//~ if (cell[n][i][j].A[2] == 0) array[2] = 0;
	//~ if (cell[n][i][j].A[3] == 0) array[3] = 0;
	//~ if (cell[n][i][j].A[4] == 0) array[4] = 0;
    //~
    //~ array[5] = array[1] > 0 ? 1 : 0;
    //~ array[6] = array[2] > 0 ? 0 : 1;
    //~ array[7] = array[3] > 0 ? 1 : 0;
    //~ array[8] = array[4] > 0 ? 0 : 1;
}

/* Final stage z calculation */
double final_calc_z (cell2d * previousCell, gasCell * cell, int n, int i, int j) {
	double z1 = ((*previousCell)[n][i-1][j].z + (*previousCell)[n][i][j].z) / 2;
	double z2 = ((*previousCell)[n][i+1][j].z + (*previousCell)[n][i][j].z) / 2;
	double z3 = ((*previousCell)[n][i][j-1].z + (*previousCell)[n][i][j].z) / 2;
	double z4 = ((*previousCell)[n][i][j+1].z + (*previousCell)[n][i][j].z) / 2;
	double result = (*previousCell)[n][i][j].z +
		(
			cell->P[0] / I_k -
			cell->Vx[0] * (z2 - z1) / dx -
			cell->Vr[0] * (z4 - z3) / ((fabs(j-axis_j-0.5))*pow(dr,2))
		) * dt;
	if (result > max_z) {
		return (*previousCell)[n][i][j].z;
	} else if (result < 0) {
		return 0;
	} else {
		return result;
	}
}

/* Final stage psi calculation */
double final_calc_psi (cell2d * previousCell, gasCell * cell, int n, int i, int j) {
	double psi1 = ((*previousCell)[n][i-1][j].psi + (*previousCell)[n][i][j].psi) / 2;
	double psi2 = ((*previousCell)[n][i+1][j].psi + (*previousCell)[n][i][j].psi) / 2;
	double psi3 = ((*previousCell)[n][i][j-1].psi + (*previousCell)[n][i][j].psi) / 2;
	double psi4 = ((*previousCell)[n][i][j+1].psi + (*previousCell)[n][i][j].psi) / 2;
	double result = (*previousCell)[n][i][j].psi +
		(
			(kappa + 2*kappa*lambda*cell->z) * cell->P[0] / I_k -
			cell->Vx[0] * (psi2 - psi1) / dx -
			cell->Vr[0] * (psi4 - psi3) / ((fabs(j-axis_j-0.5))*pow(dr,2))
		) * dt;
	if (result > 1) {
		return 1;
	} else if (result < 0) {
		return 0;
	} else {
		return result;
	}
}


double euler_z (cell2d * previousCell, gasCell * cell, int n, int i, int j) {
	double result = (*previousCell)[n][i][j].final_z +
		cell->P[0] / I_k * dt;
	if (result > max_z) {
		return max_z;
	} else {
		return result;
	}
}
double euler_psi (cell2d * previousCell, gasCell * cell, int n, int i, int j) {
	double result = (*previousCell)[n][i][j].final_psi +
		(kappa + 2*kappa*lambda*cell->bar_z) * cell->P[0] / I_k * dt;
	if (result > 1) {
		return 1;
	} else {
		return result;
	}
}

double new_final_z (cell2d& cell, int i, int j, int n,
        double dx, double dr, double dt) {

	gasCell curCell = cell[n][i][j];

	unsigned NEEDED_COND = 0;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += BAR_VX_PAR;
	NEEDED_COND += BAR_VR_PAR;
	NEEDED_COND += Z_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

    double result = brd[Z_POS].ij * brd[RHO_POS].ij / cell[n+1][i][j].rho
    +
    (
		curCell.D[1] * brd[Z_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[2] * brd[Z_POS].i1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[3] * brd[Z_POS].i_1j_1 * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[4] * brd[Z_POS].i_1j1 * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[1]) * brd[Z_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[2]) * brd[Z_POS].i_1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[3]) * brd[Z_POS].i_1j * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[4]) * brd[Z_POS].i_1j * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
    )
	;
    return result;
}

double new_final_psi (cell2d& cell, int i, int j, int n,
        double dx, double dr, double dt) {

	gasCell curCell = cell[n][i][j];

	unsigned NEEDED_COND = 0;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += BAR_VX_PAR;
	NEEDED_COND += BAR_VR_PAR;
	NEEDED_COND += PSI_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

	double result = brd[PSI_POS].ij * brd[RHO_POS].ij / cell[n+1][i][j].rho
	    +
	    (
			curCell.D[1] * brd[PSI_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    +
		    curCell.D[2] * brd[PSI_POS].i1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    +
		    curCell.D[3] * brd[PSI_POS].i_1j_1 * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    +
		    curCell.D[4] * brd[PSI_POS].i_1j1 * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    -
		    (1-curCell.D[1]) * brd[PSI_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    -
		    (1-curCell.D[2]) * brd[PSI_POS].i_1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    -
		    (1-curCell.D[3]) * brd[PSI_POS].i_1j * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
		    -
		    (1-curCell.D[4]) * brd[PSI_POS].i_1j * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    )
		;
    return result;
}








/* Final stage P calculation */
double final_calc_p(gasCell * prevCell, gasCell * cell) {
	double result = (cell->e - (pow(cell->Vx[0],2)+pow(cell->Vr[0],2))/2 ) * (k-1) /
		(
//			1/cell->rho // No powder present
			1/cell->rho - (1 - cell->final_psi)/delta - alpha_k * cell->final_psi // BMSTU var
//			alpha_k // Abel (Dupre) equation
		)
//		- (pow(cell->Vx[0],2)+pow(cell->Vr[0],2))/2 // From Ershov, UDK 519.6:532.6, #775, 2007, str. 159-173
		;
	if (result < 0) {
		broken_dt = true;
		cout << "P < 0" << endl;
		cout << "alpha_k = " << alpha_k << endl;
		cout << "E = " << cell->e << endl
			<< "Internal energy = " << cell->e  * (k-1) * cell->rho << endl
			<< "rho = " << cell->rho << endl
			<< "final_psi = " << cell->final_psi << endl
			<< "delta = " << delta << endl
			<< "alpha = " << cell->alpha << endl
			<< "m*dt = " << cell->m*dt << endl
			<< "Vx[0] = " << cell->Vx[0] << endl
			<< "Vr[0] = " << cell->Vr[0] << endl;
		getchar();
	}
	return result;
}

/* Final stage Vx calculation */
double final_calc_Vx(cell2d& cell, int i, int j, int n,
        double dx, double dr, double dt) {

	gasCell curCell = cell[n][i][j];

	unsigned NEEDED_COND = 0;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += BAR_VX_PAR;
	NEEDED_COND += BAR_VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

    double result = brd[BAR_VX_POS].ij * brd[RHO_POS].ij / cell[n+1][i][j].rho
    +
    (
		curCell.D[1] * brd[BAR_VX_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[2] * brd[BAR_VX_POS].i1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[3] * brd[BAR_VX_POS].i_1j_1 * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[4] * brd[BAR_VX_POS].i_1j1 * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[1]) * brd[BAR_VX_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[2]) * brd[BAR_VX_POS].i_1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[3]) * brd[BAR_VX_POS].i_1j * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[4]) * brd[BAR_VX_POS].i_1j * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
    )
	;
	if (result != result) {
		printf("Vx is NaN!");
		getchar();
	}
    return result;
}

/* Final stage Vr calculation */
double final_calc_Vr(cell2d& cell, int i, int j, int n,
        double dx, double dr, double dt) {

	gasCell curCell = cell[n][i][j];

	unsigned NEEDED_COND = 0;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += BAR_VX_PAR;
	NEEDED_COND += BAR_VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

    double result = brd[BAR_VR_POS].ij * brd[RHO_POS].ij / cell[n+1][i][j].rho
    +
    (
	    curCell.D[1] * brd[BAR_VR_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[2] * brd[BAR_VR_POS].i1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[3] * brd[BAR_VR_POS].ij_1 * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[4] * brd[BAR_VR_POS].ij1 * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[1]) * brd[BAR_VR_POS].ij * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[2]) * brd[BAR_VR_POS].ij * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[3]) * brd[BAR_VR_POS].ij * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[4]) * brd[BAR_VR_POS].ij * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	)
    ;
    return result;
}


/* Final stage E calculation */
double final_calc_e(cell2d& cell, int i, int j, int n,
        double dx, double dr, double dt) {

	gasCell curCell = cell[n][i][j];

	unsigned NEEDED_COND = 0;
	NEEDED_COND += RHO_PAR;
	NEEDED_COND += BAR_VX_PAR;
	NEEDED_COND += BAR_E_PAR;
	BorderCond brd[10];
	calculateBorder(n, cell[n], NEEDED_COND, brd);

    double result = brd[BAR_E_POS].ij * brd[RHO_POS].ij / cell[n+1][i][j].rho
    +
    (
	    curCell.D[1] * brd[BAR_E_POS].i_1j * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[2] * brd[BAR_E_POS].i1j * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[3] * brd[BAR_E_POS].ij_1 * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    +
	    curCell.D[4] * brd[BAR_E_POS].ij1 * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[1]) * brd[BAR_E_POS].ij * curCell.dM[1] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[2]) * brd[BAR_E_POS].ij * curCell.dM[2] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[3]) * brd[BAR_E_POS].ij * curCell.dM[3] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	    -
	    (1-curCell.D[4]) * brd[BAR_E_POS].ij * curCell.dM[4] / ((fabs(j-axis_j-0.5)) * cell[n+1][i][j].rho * curCell.A[0] * dx * pow(dr,2))
	)
    ;
	if (result < 0) {
		cout << "E < 0" << endl;
		cout << "i = " << i << endl
			<< "j = " << j << endl;
		cout << setiosflags(ios::fixed) << setprecision(16) <<
			"bar_E       = " << cell[n][i][j].bar_e << endl <<
			"lower bar_E = " << cell[n][i][j-1].bar_e << endl <<
			"upper bar_E = " << cell[n][i][j+1].bar_e << endl <<
			"E       = " << cell[n+1][i][j].e << endl <<
			"lower E = " << cell[n+1][i][j-1].e << endl <<
			"upper E = " << cell[n+1][i][j+1].e << endl <<
			"P[0] = " << cell[n][i][j].P[0] << endl <<
			"A[0] = " << cell[n][i][j].A[0] << endl <<
			"bar_Vx       = " << cell[n][i][j].bar_Vx[0] << endl <<
			"lower bar_Vx = " << cell[n][i][j-1].bar_Vx[0] << endl <<
			"upper bar_Vx = " << cell[n][i][j+1].bar_Vx[0] << endl <<
			"rho       = " << cell[n+1][i][j].rho << endl <<
			"lower rho = " << cell[n+1][i][j-1].rho << endl <<
			"upper rho = " << cell[n+1][i][j+1].rho << endl <<
			"prev rho       = " << cell[n][i][j].rho << endl <<
			"lower prev rho = " << cell[n][i][j-1].rho << endl <<
			"upper prev rho = " << cell[n][i][j+1].rho << endl <<
			"D       = {" << cell[n][i][j].D[1] << ", " << cell[n][i][j].D[2] << ", " << cell[n][i][j].D[3] << ", " << cell[n][i][j].D[4] << "} " << endl <<
			"dM       = {" << cell[n][i][j].dM[1] << ", " << cell[n][i][j].dM[2] << ", " << cell[n][i][j].dM[3] << ", " << cell[n][i][j].dM[4] << "} " << endl <<
			"lower dM = {" << cell[n][i][j-1].dM[1] << ", " << cell[n][i][j-1].dM[2] << ", " << cell[n][i][j-1].dM[3] << ", " << cell[n][i][j-1].dM[4] << "} " << endl <<
			"upper dM = {" << cell[n][i][j+1].dM[1] << ", " << cell[n][i][j+1].dM[2] << ", " << cell[n][i][j+1].dM[3] << ", " << cell[n][i][j+1].dM[4] << "} " << endl << endl;
		getchar();
	}

	if (result != result) {
		cout << "final_e is NaN" << endl;
		cout << "i = " << i << endl
			<< "j = " << j << endl;
		cout << setiosflags(ios::fixed) << setprecision(16) <<
			"bar_E       = " << cell[n+1][i][j].bar_e << endl <<
			"lower bar_E = " << cell[n+1][i][j-1].bar_e << endl <<
			"upper bar_E = " << cell[n+1][i][j+1].bar_e << endl <<
			"E       = " << cell[n+1][i][j].e << endl <<
			"lower E = " << cell[n+1][i][j-1].e << endl <<
			"upper E = " << cell[n+1][i][j+1].e << endl <<
			"P[0] = " << cell[n][i][j].P[0] << endl <<
			"A[0] = " << cell[n][i][j].A[0] << endl <<
			"bar_Vx       = " << cell[n+1][i][j].bar_Vx[0] << endl <<
			"lower bar_Vx = " << cell[n+1][i][j-1].bar_Vx[0] << endl <<
			"upper bar_Vx = " << cell[n+1][i][j+1].bar_Vx[0] << endl <<
			"rho       = " << cell[n+1][i][j].rho << endl <<
			"lower rho = " << cell[n+1][i][j-1].rho << endl <<
			"upper rho = " << cell[n+1][i][j+1].rho << endl <<
			"prev rho       = " << cell[n][i][j].rho << endl <<
			"lower prev rho = " << cell[n][i][j-1].rho << endl <<
			"upper prev rho = " << cell[n][i][j+1].rho << endl <<
			"dM       = {" << cell[n+1][i][j].dM[1] << ", " << cell[n+1][i][j].dM[2] << ", " << cell[n+1][i][j].dM[3] << ", " << cell[n+1][i][j].dM[4] << "} " << endl <<
			"lower dM = {" << cell[n+1][i][j-1].dM[1] << ", " << cell[n+1][i][j-1].dM[2] << ", " << cell[n+1][i][j-1].dM[3] << ", " << cell[n+1][i][j-1].dM[4] << "} " << endl <<
			"upper dM = {" << cell[n+1][i][j+1].dM[1] << ", " << cell[n+1][i][j+1].dM[2] << ", " << cell[n+1][i][j+1].dM[3] << ", " << cell[n+1][i][j+1].dM[4] << "} " << endl << endl;
		getchar();
	}

    return result;
}

double smooth_Vr(cell2d * cell) {
	double result;

	unsigned NEEDED_COND = 0;
	NEEDED_COND += VR_PAR;
	BorderCond brd[10];
	calculateBorder(n, (*cell)[n], NEEDED_COND, brd);

	result = 0.05 * brd[VR_POS].i_1j_1 + 0.05 * brd[VR_POS].ij_1 + 0.05 * brd[VR_POS].i1j_1 +
			0.05 * brd[VR_POS].i_1j + 0.6 * brd[VR_POS].ij + 0.05 * brd[VR_POS].i1j +
			0.05 * brd[VR_POS].i_1j1 + 0.05 * brd[VR_POS].ij1 + 0.05 * brd[VR_POS].i1j1;
	return result;
}
