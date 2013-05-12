#include "globals.h"
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

//~ double polygonArea(double *X, double *Y, int points) {
//~ 
  //~ double  area=0. ;
  //~ int     i, j=points-1  ;
//~ 
  //~ for (i=0; i<points; i++) {
    //~ area+=(X[j]+X[i])*(Y[j]-Y[i]); j=i; }
//~ 
//~ return area*.5; }

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
	    printf("Cell at (%d, %d) with type 3; A[1] = %4.4f, A[2] = %4.4f, A[3] = %4.4f, A[4] = %4.4f\n", i,j,array[1],array[2],array[3],array[4]);
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
		
	
	default:
		break;
	}
}

//~ WeightVector wightVectorsCalc(cell2d cell, int i, int j, int n) {
	//~ Point2D vertices[4];
	//~ Int2D vertices_ij[4];
	//~ WeightVector result;
	//~ /** 
	 //~ * Transform points from center of the cell
	 //~ * 
	 //~ * 	[3]	+---+ [2]
	 //~ * 		|	|
	 //~ * 	[0]	+---+ [1]
	 //~ * 
	 //~ * 			|  cos(a)  -sin(a)	|
	 //~ * R(-a) = 	|					|
	 //~ * 			|  sin(a)  cos(a)	|
	 //~ * 
	 //~ * weightCell[0] is i+1 cell (X), [1] is j+1 (Y)
	 //~ * 
	 //~ **/
	//~ int curI = i;
	//~ int curJ = j;
	//~ double ybegin = j+1;
	//~ double yend = j;
	//~ double xbegin = i;
	//~ double xend = i+1;
	//~ double cos_a;
	//~ double sin_a;
	//~ double cos_2a;
	//~ double sin_2a;
	//~ printf("\n\n**************************************************\n\n");
	//~ for (unsigned int weightCell = 0; weightCell < 2; weightCell++) {
		//~ if (weightCell == 0) {
			//~ curI = i+1; 
			//~ curJ = j;
		//~ } else if (weightCell == 1) {
			//~ curI = i; 
			//~ curJ = j+1;
		//~ }
		//~ 
		//~ /* TODO: here i assume we have angle from left top to right bottom */
		//~ switch (cell.at(n).at(i).at(j).type) {
			//~ case 1:
			//~ xbegin = i*dx;
			//~ xend = (i+1)*dx;
			//~ ybegin = j*dr + cell.at(n).at(i).at(j).r_1*dr;
			//~ yend = j*dr + cell.at(n).at(i).at(j).r_2*dr;
			//~ break;
			//~ 
			//~ case 3:
			//~ xbegin = i*dx + cell.at(n).at(i).at(j).x_2*dx;
			//~ xend = i*dx + cell.at(n).at(i).at(j).x_1*dx;
			//~ ybegin = (j+1)*dr;
			//~ yend = j*dr;
			//~ break;
			//~ 
			//~ case 8:
			//~ xbegin = i*dx + cell.at(n).at(i).at(j).x_2*dx;
			//~ xend = (i+1)*dx;
			//~ ybegin = (j+1)*dr;
			//~ yend = j*dr + cell.at(n).at(i).at(j).r_2*dr;
			//~ break;
			//~ 
			//~ case 10:
			//~ xbegin = i*dx;
			//~ xend = i*dx + cell.at(n).at(i).at(j).x_1*dx;
			//~ ybegin = j*dr + cell.at(n).at(i).at(j).r_1*dr;
			//~ yend = j*dr;
			//~ break;
			//~ 
			//~ default:
			//~ break;
		//~ }
		//~ cos_a = (xend - xbegin) / sqrt(pow(yend-ybegin,2)+pow(xend-xbegin,2));
		//~ sin_a = (yend - ybegin) / sqrt(pow(yend-ybegin,2)+pow(xend-xbegin,2));
		//~ cos_2a = pow(cos_a,2) - pow(sin_a,2);
		//~ sin_2a = 2*cos_a*sin_a;
		//~ 
		//~ // (curI+0.5)*dx returns us to the center of the cell
		//~ vertices[0].x = (curI)*dx; 
		//~ vertices[1].x = (curI+1)*dx; 
		//~ vertices[2].x = (curI+1)*dx; 
		//~ vertices[3].x = (curI)*dx; 
		//~ vertices[0].y = (curJ)*dr;
		//~ vertices[1].y = (curJ)*dr;
		//~ vertices[2].y = (curJ+1)*dr;
		//~ vertices[3].y = (curJ+1)*dr;
		//~ 
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ cout << "x_1 = " << cell.at(n).at(i).at(j).x_1 << ", r_1 = " << cell.at(n).at(i).at(j).r_1 << ", x_2 = " << cell.at(n).at(i).at(j).x_2 << ", r_2 = " << cell.at(n).at(i).at(j).r_2 << endl; 
		//~ cout << "Line begin = " << xbegin << ":" << ybegin << ", line end = " << xend << ":" << yend << endl;
		//~ cout << "Line angle = " << asin(sin_a)*180/M_PI << endl;
		//~ cout << "Orig Vertices at i = " << i << ", j = " << j << ", weightCell = " << weightCell << endl;
		//~ cout << "0: " << vertices[0].x << ":" << vertices[0].y << endl;
		//~ cout << "1: " << vertices[1].x << ":" << vertices[1].y << endl;
		//~ cout << "2: " << vertices[2].x << ":" << vertices[2].y << endl;
		//~ cout << "3: " << vertices[3].x << ":" << vertices[3].y << endl;
		//~ getchar();
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ 
		//~ // Get distance from each point to line
		//~ for (unsigned int idx = 0; idx < 4; idx++) {
			//~ double dist = (
				//~ (xend - xbegin) * (ybegin - vertices[idx].y) - 
				//~ (xbegin - vertices[idx].x) * (yend - ybegin)
			    //~ ) / sqrt(pow(xend-xbegin,2)+pow(yend-ybegin,2));
			//~ double diffX = 2 * dist * sin_a;
			//~ double diffY = 2 * dist * cos_a;
			//~ vertices[idx].x -= diffX;
			//~ vertices[idx].y += diffY;
			//~ vertices_ij[idx].i = floor(vertices[idx].x/dx); 
			//~ vertices_ij[idx].j = floor(vertices[idx].y/dr);
			//~ printf("Point %u, distance total: %4.4f, x: %4.4f, y: %4.4f\n", idx,dist,diffX,diffY);
		//~ }
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ cout << "Vertices at i = " << i << ", j = " << j << ", weightCell = " << weightCell << endl;
		//~ cout << "0: " << vertices[0].x << ":" << vertices[0].y << endl;
		//~ cout << "1: " << vertices[1].x << ":" << vertices[1].y << endl;
		//~ cout << "2: " << vertices[2].x << ":" << vertices[2].y << endl;
		//~ cout << "3: " << vertices[3].x << ":" << vertices[3].y << endl;
		//~ cout << "Vertices_ij at i = " << i << ", j = " << j << ", weightCell = " << weightCell << endl;
		//~ cout << "0: " << vertices_ij[0].i << ":" << vertices_ij[0].j << endl;
		//~ cout << "1: " << vertices_ij[1].i << ":" << vertices_ij[1].j << endl;
		//~ cout << "2: " << vertices_ij[2].i << ":" << vertices_ij[2].j << endl;
		//~ cout << "3: " << vertices_ij[3].i << ":" << vertices_ij[3].j << endl;
		//~ getchar();
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /**
		 //~ * Axis intersections and point types
		 //~ * 
		 //~ * Point types:
		 //~ * 0 - vertex
		 //~ * 1 - axis intersection
		 //~ * 2 - internal point 
		 //~ * 
		 //~ **/
		//~ Point2D point;
		//~ vector <Point2D> points;
		//~ vector <int> pointType;
		//~ for (int idx1 = 0; idx1 < 4; idx1++) {
			//~ point.x = vertices[idx1].x; point.y = vertices[idx1].y;
			//~ points.push_back(point);
			//~ pointType.push_back(0);
			//~ if (fabs(fmod(point.x, dx)) < pow(10,-6) && fabs(fmod(point.y, dr)) < pow(10,-6)) continue;
			//~ int idx2 = idx1 != 3 ? idx1 + 1 : 0;
			//~ // One intersection point on Y axis
			//~ if (vertices_ij[idx1].i != vertices_ij[idx2].i && vertices_ij[idx1].j == vertices_ij[idx2].j) {
				//~ double interX = floor(fmax(vertices[idx1].x,vertices[idx2].x)/dx)*dx;
				//~ double interY = vertices[idx2].y + (interX - vertices[idx2].x) * (vertices[idx2].y - vertices[idx1].y) / (vertices[idx2].x - vertices[idx1].x);
				//~ point.x = interX; point.y = interY;
				//~ points.push_back(point);
				//~ pointType.push_back(1);
			//~ } else
			//~ // One intersection point on X axis
			//~ if (vertices_ij[idx1].i == vertices_ij[idx2].i && vertices_ij[idx1].j != vertices_ij[idx2].j) {
				//~ u_int idx_max = vertices[idx1].y > vertices[idx2].y ? idx1 : idx2;
				//~ double interY = floor(vertices[idx_max].y/dr)*dr;
				//~ double interX = vertices[idx2].x + (interY - vertices[idx2].y) * (vertices[idx2].x - vertices[idx1].x) / (vertices[idx2].y - vertices[idx1].y);
				//~ point.x = interX; point.y = interY;
				//~ points.push_back(point);
				//~ pointType.push_back(1);
			//~ } else
			//~ // If has 2 intersection points - first serve closest (because triangle gen is buggy)
			//~ if (vertices_ij[idx1].i != vertices_ij[idx2].i && vertices_ij[idx1].j != vertices_ij[idx2].j) {
				//~ double interX1 = floor(fmax(vertices[idx1].x,vertices[idx2].x)/dx)*dx;
				//~ double interY1 = vertices[idx2].y + (interX1 - vertices[idx2].x) * (vertices[idx2].y - vertices[idx1].y) / (vertices[idx2].x - vertices[idx1].x);
				//~ double interY2 = floor(fmax(vertices[idx1].y,vertices[idx2].y)/dr)*dr;
				//~ double interX2 = vertices[idx2].x + (interY2 - vertices[idx2].y) * (vertices[idx2].x - vertices[idx1].x) / (vertices[idx2].y - vertices[idx1].y);
				//~ if (pow(vertices[idx1].x-interX1,2) + pow(vertices[idx1].y-interY1,2) < pow(vertices[idx1].x-interX2,2) + pow(vertices[idx1].y-interY2,2)) {
					//~ point.x = interX1; point.y = interY1;
					//~ points.push_back(point);
					//~ pointType.push_back(1);
					//~ point.x = interX2; point.y = interY2 ;
					//~ points.push_back(point);
					//~ pointType.push_back(1);
				//~ } else {
					//~ point.x = interX2; point.y = interY2;
					//~ points.push_back(point);
					//~ pointType.push_back(1);
					//~ point.x = interX1; point.y = interY1;
					//~ points.push_back(point);
					//~ pointType.push_back(1);
				//~ }
			//~ }
		//~ }
		//~ // Internal points and global type calculation
		//~ int globalType = -1;
		//~ int max_i_diff = 0;
		//~ int max_j_diff = 0;
		//~ unsigned int max_i_point[2] = {0};
		//~ unsigned int max_j_point[2] = {0};
		//~ for (unsigned int idx2 = 0; idx2 < 4; idx2++) {
			//~ for (unsigned int idx3 = 0; idx3 < 4; idx3++) {
				//~ if (vertices_ij[idx2].i - vertices_ij[idx3].i > max_i_diff) {
					//~ max_i_diff = vertices_ij[idx2].i - vertices_ij[idx3].i;
					//~ max_i_point[0] = idx2;
					//~ max_i_point[1] = idx3;
				//~ }
				//~ if (vertices_ij[idx2].j - vertices_ij[idx3].j > max_j_diff) {
					//~ max_j_diff = vertices_ij[idx2].j - vertices_ij[idx3].j;
					//~ max_j_point[0] = idx2;
					//~ max_j_point[1] = idx3;
				//~ }
			//~ }
		//~ }
		//~ 
		//~ // Cleaning points
		//~ bool changed = true;
		//~ while (changed) {
			//~ changed = false;
			//~ for (unsigned int idx = 0; idx < points.size(); idx++) {
				//~ for (unsigned int idx2 = 0; idx2 < points.size(); idx2++) {
					//~ if (fabs(points.at(idx2).x - points.at(idx).x) < pow(10,-6) && 
					    //~ fabs(points.at(idx2).y - points.at(idx).y) < pow(10,-6) && 
					    //~ pointType.at(idx) != 2 &&
					    //~ idx2 != idx) {
						//~ points.erase(points.begin()+idx2);
						//~ pointType.erase(pointType.begin()+idx2);
						//~ changed = true;
					//~ }
				//~ }
			//~ }
		//~ }
//~ 
		//~ cout << "Cell : " << i << ":" << j << endl;
		//~ cout << "max_i_diff = " << max_i_diff << ", max_i_points : " << max_i_point[0] << " : " << max_i_point[1] << endl;
		//~ cout << "max_j_diff = " << max_j_diff << ", max_j_points : " << max_j_point[0] << " : " << max_j_point[1] <<  endl;
		//~ if (max_i_diff == 1 && max_j_diff == 1) {
			//~ double internalX = floor(fmax(fmax(vertices[0].x,vertices[1].x),fmax(vertices[2].x,vertices[3].x))/dx)*dx;
			//~ double internalY = floor(fmax(fmax(vertices[0].y,vertices[1].y),fmax(vertices[2].y,vertices[3].y))/dr)*dr;
			//~ point.x = internalX; point.y = internalY;
			//~ points.push_back(point);
			//~ pointType.push_back(2);
			//~ globalType = 0;
		//~ } else {
			//~ unsigned int prevPointI[2] = {0};
			//~ unsigned int nextPointI[2] = {0};
			//~ unsigned int prevPointJ[2] = {0};
			//~ unsigned int nextPointJ[2] = {0};
			//~ prevPointI[0] = max_i_point[0] == 0 ? points.size()-1 : max_i_point[0]-1;
			//~ prevPointI[1] = max_i_point[1] == 0 ? points.size()-1 : max_i_point[1]-1;
			//~ nextPointI[0] = max_i_point[0] == points.size()-1 ? 0 : max_i_point[0]+1;
			//~ nextPointI[1] = max_i_point[1] == points.size()-1 ? 0 : max_i_point[1]+1;
			//~ prevPointJ[0] = max_j_point[0] == 0 ? points.size()-1 : max_j_point[0]-1;
			//~ prevPointJ[1] = max_j_point[1] == 0 ? points.size()-1 : max_j_point[1]-1;
			//~ nextPointJ[0] = max_j_point[0] == points.size()-1 ? 0 : max_j_point[0]+1;
			//~ nextPointJ[1] = max_j_point[1] == points.size()-1 ? 0 : max_j_point[1]+1;
			//~ 
			//~ unsigned int max_i_difference = 0;
			//~ unsigned int max_j_difference = 0;
			//~ if (fabs(floor(points.at(nextPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx)) > max_i_difference) max_i_difference = fabs(floor(points.at(nextPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx));
			//~ if (fabs(floor(points.at(prevPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx)) > max_i_difference) max_i_difference = fabs(floor(points.at(prevPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx));
			//~ if (fabs(floor(points.at(nextPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr)) > max_j_difference) max_j_difference = fabs(floor(points.at(nextPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr));
			//~ if (fabs(floor(points.at(prevPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr)) > max_j_difference) max_j_difference = fabs(floor(points.at(prevPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr));
			//~ 
			//~ //if (max_i_difference == 1 && max_j_difference == 1) {
				//~ double internalX = floor(fmax(fmax(vertices[0].x,vertices[1].x),fmax(vertices[2].x,vertices[3].x))/dx)*dx;
				//~ double internalY = floor(fmax(fmax(vertices[0].y,vertices[1].y),fmax(vertices[2].y,vertices[3].y))/dr)*dr;
				//~ point.x = internalX; point.y = internalY;
				//~ points.push_back(point);
				//~ pointType.push_back(2);
				//~ globalType = 1;
			//~ //}
		//~ }
		//~ /** TODO: other types of internal points **/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ for (unsigned int pointNum = 0; pointNum < points.size(); pointNum++) {
			//~ cout << "Point " << pointNum << " at i = " << i << ", j = " << j << " : " << points.at(pointNum).x << "::" << points.at(pointNum).y << ", type = " << pointType.at(pointNum) << endl;
		//~ }
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ /******************************************************************************************/
		//~ // Triangles and weights
		//~ std::vector < Triangle2D > triangleVector;
		//~ std::vector < Trapezoid2D > trapezoidVector;
		//~ double origArea = dx*dr;
		//~ double area;
		//~ 
		//~ 
		//~ // We'll now make a vector of all cells where our polygon is
		//~ vector <Int2D> cells; 
		//~ Int2D newCell;
		//~ unsigned int i_left, i_right, j_top, j_bottom = 0;
		//~ // Left and right
		//~ i_left = vertices_ij[max_i_point[0]].i;
		//~ i_right = vertices_ij[max_i_point[1]].i;
		//~ if (i_left > i_right) {
			//~ unsigned int tmp = i_left;
			//~ i_left = i_right;
			//~ i_right = tmp;
		//~ }
		//~ // Top and bottom
		//~ j_bottom = vertices_ij[max_j_point[0]].j;
		//~ j_top = vertices_ij[max_j_point[1]].j;
		//~ if (j_bottom > j_top) {
			//~ unsigned int tmp = j_bottom;
			//~ j_bottom = j_top;
			//~ j_top = tmp;
		//~ }
		//~ printf("\nWe have the following cells: \n");
		//~ // Make cells
		//~ for (unsigned int idx2 = i_left; idx2 <= i_right; idx2++) {
			//~ for (unsigned int idx3 = j_bottom; idx3 <= j_top; idx3++) {
				//~ newCell.i = idx2; newCell.j = idx3;
				//~ cells.push_back(newCell);
				//~ printf("Cell %u: %d:%d\n", (unsigned int)(cells.size()-1), idx2, idx3);
			//~ }
		//~ }
		//~ printf("\nNow will arrange points in those cells\n");
		//~ // We'll use center of each of those cells to determine how much points we have in each cell
		//~ 
		//~ for (unsigned int idx = 0; idx < cells.size(); idx++) {
			//~ vector <Point2D> pointsInCell;
			//~ vector <int> pointsInCellTypes;
			//~ // If vertex - check for outer triangle or trapezoid
			//~ double delta = pow(10,-10);
			//~ for (unsigned int idx2 = 0; idx2 < points.size(); idx2++) {
				//~ bool rule[4];
				//~ rule[0] = points.at(idx2).x > cells.at(idx).i*dx - delta;
				//~ rule[1] = points.at(idx2).x < (cells.at(idx).i+1)*dx + delta;
				//~ rule[2] = points.at(idx2).y > cells.at(idx).j*dr - delta;
				//~ rule[3] = points.at(idx2).y < (cells.at(idx).j+1)*dr + delta;
				//~ if (rule[0] && rule[1] && rule[2] && rule[3]) {
					//~ pointsInCell.push_back(points.at(idx2));
					//~ pointsInCellTypes.push_back(pointType.at(idx2));
					//~ printf("Cell %u: %4.4f:%4.4f, type %d \n", 
					    //~ idx,points.at(idx2).x,points.at(idx2).y,pointType.at(idx2));
				//~ }
			//~ }
			//~ // If we have internal point - fix ID
			//~ if(std::find(pointsInCellTypes.begin(), pointsInCellTypes.end(), 2) != pointsInCellTypes.end()) {
				//~ printf("We have an internal point!\n");
				//~ unsigned int intIdx = 0;
				//~ Point2D intPoint;
				//~ for (unsigned int idx2 = 0; idx2 < pointsInCellTypes.size(); idx2++) {
					//~ if (pointsInCellTypes.at(idx2) == 2) {
						//~ intIdx = idx2;
						//~ intPoint = pointsInCell.at(idx2);
						//~ break;
					//~ }
				//~ }
				//~ double minHor = 1000, minVer = 1000;
				//~ int minHorIdx = 0, minVerIdx = 0;
				//~ for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
					//~ if (idx2 != intIdx) {
						//~ Point2D point = pointsInCell.at(idx2);
						//~ printf("Testing point %u\n", idx2);
						//~ 
						//~ if (fabs(point.x - intPoint.x) < minHor && fabs(point.y - intPoint.y) < pow(10,-6)
						    //~ && pointsInCellTypes.at(idx2) == 1) {
							//~ printf("hor = %4.4f\n", point.x - intPoint.x);
							//~ minHor = point.x - intPoint.x;
							//~ minHorIdx = idx2;
						//~ }
						//~ 
						//~ if (fabs(point.y - intPoint.y) < minVer && fabs(point.x - intPoint.x) < pow(10,-6)
						    //~ && pointsInCellTypes.at(idx2) == 1) {
							//~ printf("ver = %4.4f\n", point.x - intPoint.x);
							//~ minVer = point.y - intPoint.y;
							//~ minVerIdx = idx2;
						//~ }
					//~ }
				//~ }
				//~ printf("Closest points to internal: hor %4.4f:%4.4f, ver %4.4f:%4.4f\n\n", 
				    //~ pointsInCell.at(minHorIdx).x,pointsInCell.at(minHorIdx).y,pointsInCell.at(minVerIdx).x,pointsInCell.at(minVerIdx).y);
				//~ 
				//~ if (minVerIdx != 0 && minHorIdx != 0) {
				//~ std::vector<Point2D>::iterator it = pointsInCell.begin();
				//~ pointsInCell.erase(it+intIdx);
					//~ if (minVerIdx > minHorIdx) {
						//~ pointsInCell.insert(it+minVerIdx, intPoint);
					//~ } else {
						//~ pointsInCell.insert(it+minHorIdx, intPoint);
					//~ }
				//~ }
				//~ 
				//~ printf("Points in current cell in order: ");
				//~ for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
					//~ printf("%4.4f:%4.4f, ",pointsInCell.at(idx2).x, pointsInCell.at(idx2).y);
				//~ }
				//~ printf("\n");
			//~ }
			//~ printf("PointsArray size: %u\n", (unsigned int)pointsInCell.size());
			//~ double XPointsArray[pointsInCell.size()];
			//~ double YPointsArray[pointsInCell.size()];
			//~ for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
				//~ XPointsArray[idx2] = pointsInCell.at(idx2).x;
				//~ YPointsArray[idx2] = pointsInCell.at(idx2).y;
				//~ printf("Point %4.4f:%4.4f;  ",XPointsArray[idx2],YPointsArray[idx2]);
			//~ }
			//~ area = fabs(polygonArea(XPointsArray, YPointsArray, pointsInCell.size()));
			//~ printf("Polygon area: %8.8f\n", area);
			//~ WeightPart weightPart;
			//~ weightPart.weight = area / origArea;
			//~ weightPart.i = cells.at(idx).i;
			//~ weightPart.j = cells.at(idx).j;
			//~ printf("Weight of this cell: %4.4f\n\n\n", weightPart.weight);
			//~ if (weightCell == 0) {
				//~ cell.at(n).at(i).at(j).weightVector.x.push_back(weightPart);
			//~ } else if (weightCell == 1) {
				//~ cell.at(n).at(i).at(j).weightVector.y.push_back(weightPart);
			//~ }
		//~ }
	//~ }
	//~ printf("Total weight parts at %d in cell %d:%d\nby x: %u\nby y: %u\n", n,i,j,
		//~ (unsigned int) cell.at(n).at(i).at(j).weightVector.x.size(),
		//~ (unsigned int) cell.at(n).at(i).at(j).weightVector.y.size());
	//~ result = cell.at(n).at(i).at(j).weightVector;
	//~ return result;
//~ }

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
double euler_bar_Vx(cell2d * cell, int n, int i, int j, double dt, double dx, double dr) {
	/* Dummy cells calc */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }
    /* I'll use the following naming rule P[i-2] = P_i_2, P[i] = P, P[i+1] = P_i1 */
    double P_i_2 = 0; double P_i_1 = 0; double P = 0; double P_i1 = 0; double P_i2 = 0;
    double Vx_i_2 = 0; double Vx_i_1 = 0; double Vx = 0; double Vx_i1 = 0; double Vx_i2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 1:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    
	    P_i1 = 0;
	    Vx_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			P_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vx[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    
	    P_i1 = 0;
	    Vx_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			P_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vx[0];
	    }
	    break;
	
	// Top and left borders closed
	case 15:
	    P_i_2 = (*cell)[n][i+1][j].P[0]; P_i_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; // Mirror i-2 -> i+1, i-1 -> self
	    Vx_i_2 = (*cell)[n][i+1][j].Vx[0]; Vx_i_1 = (*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	    
	// Left border closed
	case 17:
	    P_i_2 = (*cell)[n][i+1][j].P[0]; P_i_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; // Mirror i-2 -> i+1, i-1 -> self 
	    Vx_i_2 = (*cell)[n][i+1][j].Vx[0]; Vx_i_1 = (*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Bottom and left borders closed
	case 16:
	    P_i_2 = (*cell)[n][i+1][j].P[0]; P_i_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; // Mirror i-2 -> i+1, i-1 -> self 
	    Vx_i_2 = (*cell)[n][i+1][j].Vx[0]; Vx_i_1 = (*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Top border closed
	case 13:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Bottom border closed
	case 14:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
	    break;
	
	// Right border closed
	case 19:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; // Mirror i+2 -> i-1, i+1 -> self
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i][j].Vx[0]; Vx_i2 = (*cell)[n][i-1][j].Vx[0]; 
	    break;
	
	// Top and right borders closed
	case 20:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; // Mirror i+2 -> i-1, i+1 -> self 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i][j].Vx[0]; Vx_i2 = (*cell)[n][i-1][j].Vx[0]; 
	    break;
	
	// Bottom and right borders closed
	case 21:
	    P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; // Mirror i+2 -> i-1, i+1 -> self 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i][j].Vx[0]; Vx_i2 = (*cell)[n][i-1][j].Vx[0]; 
	    break;
	
	default:
	    break;
    }
	
	/** First order **/
	//~ double P_i12 = (P_i1 + P)/2 * (1 - (k-1)*(Vx_i1 - Vx)*dt/dx);
	//~ double P_i_12 = (P_i_1 + P)/2 * (1 - (k-1)*(Vx - Vx_i_1)*dt/dx);
	double P_i12 = (P_i1 + P)/2;
	double P_i_12 = (P_i_1 + P)/2;
	double result = (*cell)[n][i][j].Vx[0] - (P_i12 - P_i_12) / dx * fmax((*cell)[n][i][j].A[1],(*cell)[n][i][j].A[2]) * dt / ((*cell)[n][i][j].rho * (*cell)[n][i][j].A[0]);
	if (fabs(result-(*cell)[n][i][j].Vx[0]) < pow(10,-15)) result = (*cell)[n][i][j].Vx[0];
	
    /** Second order by x **/
    //~ double result = (*cell)[n][i][j].Vx[0] - (1.0/12.0*P_i_2 - 2.0/3.0*P_i_1 + 2.0/3.0*P_i1 - 1.0/12.0*P_i2) / dx * fmax((*cell)[n][i][j].A[1],(*cell)[n][i][j].A[2]) * dt / ((*cell)[n][i][j].rho * (*cell)[n][i][j].A[0]);
	//~ 
	/** With Navier-Stocks **/
    //~ double result = (*cell)[n][i][j].Vx[0] +
        //~ dt  / ((*cell)[n][i][j].rho * (*cell)[n][i][j].A[0] * dx) * 
        //~ (
			//~ // Pressure
			//~ ((*cell)[n][i][j].P[1] - (*cell)[n][i][j].P[2]) * fmax((*cell)[n][i][j].A[1], (*cell)[n][i][j].A[2]) +
			//~ // Viscosity
			//~ gasMu * dx * (
				//~ (gasA+2)/pow(dx,2) * (
					//~ (*cell)[n][i+1][j].Vx[0] - 2*(*cell)[n][i][j].Vx[0] + (*cell)[n][i-1][j].Vx[0]
				//~ ) + (
					//~ (*cell)[n][i][j+1].Vx[0] - 2*(*cell)[n][i][j].Vx[0] + (*cell)[n][i][j-1].Vx[0]
				//~ ) / pow(dr,2) +
				//~ (gasA+1)/(4*dx*dr) * (
					//~ (*cell)[n][i+1][j+1].Vr[0] + (*cell)[n][i-1][j-1].Vr[0] - (*cell)[n][i-1][j+1].Vr[0] - (*cell)[n][i+1][j-1].Vr[0]
				//~ )
			//~ )
        //~ );
    //~ 
    /** DEBUG **/
    if (i == 248 && j == 5) {
	cout << "i = 248, j = 5, P{1,2} = {" << (*cell)[n][i][j].P[1] << ", " << (*cell)[n][i][j].P[2] << "}, Vx{i-1,i,i+1} = " << (*cell)[n][i-1][j].Vx[0] << ", " << (*cell)[n][i][j].Vx[0] << ", " << (*cell)[n][i+1][j].Vx[0] << "}" << endl;
    cout << "result = " << result << endl;
    }
    if (i == 249 && j == 5) {
	cout << "i = 250, j = 5, P{1,2} = {" << (*cell)[n][i][j].P[1] << ", " << (*cell)[n][i][j].P[2] << "}, Vx{i-1,i,i+1} = " << (*cell)[n][i-1][j].Vx[0] << ", " << (*cell)[n][i][j].Vx[0] << ", " << (*cell)[n][i+1][j].Vx[0] << "}" << endl;
	cout << "result = " << result << endl;
    }
    
	return result;
}

/* Vr calculation on euler stage */
double euler_bar_Vr(cell2d * cell, int n, int i, int j, double dt, double dx, double dr) {
	/* Dummy cells calc */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }
    // I'll use the following naming rule P[j-2] = P_j_2, P[j] = P, P[j+1] = P_j1
    double P_j_2 = 0; double P_j_1 = 0; double P = 0; double P_j1 = 0; double P_j2 = 0;
    double Vr_j_2 = 0; double Vr_j_1 = 0; double Vr = 0; double Vr_j1 = 0; double Vr_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 1:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 
	    
	    P_j1 = 0;
	    Vr_j1 = 0;
	    printf("Total weightVector size at cell %d:%d by y: %u\n",i,j,
			(unsigned int) (*cell)[n][i][j].weightVector.y.size());
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			P_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vr[0];
			printf("Getting weight = %10.10f from cell %d:%d\n", weight, weightCell.i, weightCell.j);
	    }
	    printf("\nP_j1 in cell %d:%d = %10.10f\n", i,j,P_j1);
	    printf("(P_j1 + P)/2 - (P_j_1 + P)/2 = %10.10f\n\n",(P_j1 + P)/2 - (P_j_1 + P)/2);
	    
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 
	    
	    P_j1 = 0;
	    Vr_j1 = 0;
	    printf("Total weightVector size at cell %d:%d by y: %u\n",i,j,
			(unsigned int) (*cell)[n][i][j].weightVector.y.size());
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			P_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vr[0];
			printf("Getting weight = %10.10f from cell %d:%d\n", 
				weight, weightCell.i, weightCell.j);
	    }
	    printf("\nP_j1 in cell %d:%d = %10.10f\n", i,j,P_j1);
	    printf("(P_j1 + P)/2 - (P_j_1 + P)/2 = %10.10f\n\n",(P_j1 + P)/2 - (P_j_1 + P)/2);
	    break;
	
	// Top and left borders closed
	case 15:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; // Mirror j+2 -> j-1, j+1 -> self
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 
	    break;
	    
	// Left border closed
	case 17:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Bottom and left borders closed
	case 16:
	    P_j_2 = (*cell)[n][i][j+1].P[0]; P_j_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; // Mirror j-2 -> j+1, j-1 -> self 
	    Vr_j_2 = (*cell)[n][i][j+1].Vr[0]; Vr_j_1 = (*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Top border closed
	case 13:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; // Mirror j+2 -> j-1, j+1 -> self
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 
	    break;
	
	// Bottom border closed
	case 14:
	    P_j_2 = (*cell)[n][i][j+1].P[0]; P_j_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; // Mirror j-2 -> j+1, j-1 -> self 
	    Vr_j_2 = (*cell)[n][i][j+1].Vr[0]; Vr_j_1 = (*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Right border closed
	case 19:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Top and right borders closed
	case 20:
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; // Mirror j+2 -> j-1, j+1 -> self
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 
	    break;
	
	// Bottom and right borders closed
	case 21:
	    P_j_2 = (*cell)[n][i][j+1].P[0]; P_j_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; // Mirror j-2 -> j+1, j-1 -> self 
	    Vr_j_2 = (*cell)[n][i][j+1].Vr[0]; Vr_j_1 = (*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	default:
	    break;
    }
    
    /** First order **/
	//~ double P_j12 = (P_j1 + P)/2 * (1 - (k-1)*(Vr_j1 - Vr)*dt/dr);
	//~ double P_j_12 = (P_j_1 + P)/2 * (1 - (k-1)*(Vr - Vr_j_1)*dt/dr);
	double P_j12 = (P_j1 + P)/2;
	double P_j_12 = (P_j_1 + P)/2;
	double result = (*cell)[n][i][j].Vr[0] - (P_j12 - P_j_12) / (*cell)[n][i][j].rho * dt / (*cell)[n][i][j].A[0] * fmax((*cell)[n][i][j].A[4],(*cell)[n][i][j].A[3]) / dr;
	if (fabs(result-(*cell)[n][i][j].Vr[0]) < pow(10,-15)) result = (*cell)[n][i][j].Vr[0];
	
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
    
    /** Second order by x **/
    //~ double result = (*cell)[n][i][j].Vr[0] - (1.0/12.0*P_j_2 - 2.0/3.0*P_j_1 + 2.0/3.0*P_j1 - 1.0/12.0*P_j2) / (*cell)[n][i][j].rho * dt / (*cell)[n][i][j].A[0] * fmax((*cell)[n][i][j].A[4],(*cell)[n][i][j].A[3]) / dr;
    //~ 
    /** With Navier-Stocks **/
    //~ double result = (*cell)[n][i][j].Vr[0] +
        //~ dt  / ((*cell)[n][i][j].rho * (*cell)[n][i][j].A[0] * dr) * 
        //~ (
			//~ // Pressure
			//~ ((*cell)[n][i][j].P[3] - (*cell)[n][i][j].P[4])*fmax((*cell)[n][i][j].A[3],(*cell)[n][i][j].A[4]) +
			//~ // Viscosity
			//~ gasMu * dr *(
				//~ (gasA+2)/pow(dr,2) * (
					//~ (*cell)[n][i][j+1].Vr[0] - 2*(*cell)[n][i][j].Vr[0] + (*cell)[n][i][j-1].Vr[0]
				//~ ) + (
					//~ (*cell)[n][i+1][j].Vr[0] - 2*(*cell)[n][i][j].Vr[0] + (*cell)[n][i+1][j].Vr[0]
				//~ ) / pow(dx,2) + 
				//~ (gasA+1)/(4*dx*dr) * (
					//~ (*cell)[n][i+1][j+1].Vx[0] + (*cell)[n][i-1][j-1].Vx[0] - (*cell)[n][i-1][j+1].Vx[0] - (*cell)[n][i+1][j-1].Vx[0]
				//~ )				
			//~ )
        //~ ) ;
    
    /** Cleaning noise **/    
    //~ if (fabs((*cell)[n][i][j].P[4] - (*cell)[n][i][j].P[3]) < pow(10,-5)/scaleV) result = (*cell)[n][i][j].Vr[0];
    
	return result;
}

void rotateVectors(double& Vx, double& Vr, LineAngle2D angle) {
	Vx = Vx*angle.cos_a - Vr*angle.sin_a;
	Vr = Vx*angle.sin_a + Vr*angle.cos_a;
}

/* Energy calculation on euler stage */
double euler_bar_e(cell2d * cell, int n, int i, int j, double dt, double dx, double dr) {
	/* Dummy cells calc */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }
    // I'll use the following naming rule P[j-2] = P_j_2, P[j] = P, P[j+1] = P_j1
    double P_i_2 = 0; double P_i_1 = 0; double P = 0; double P_i1 = 0; double P_i2 = 0;
    double P_j_2 = 0; double P_j_1 = 0; double P_j1 = 0; double P_j2 = 0;
    double Vx_i_2 = 0; double Vx_i_1 = 0; double Vx = 0; double Vx_i1 = 0; double Vx_i2 = 0;
    double Vx_j_2 = 0; double Vx_j_1 = 0; double Vx_j1 = 0; double Vx_j2 = 0;
    double Vr_i_2 = 0; double Vr_i_1 = 0; double Vr = 0; double Vr_i1 = 0; double Vr_i2 = 0;
    double Vr_j_2 = 0; double Vr_j_1 = 0; double Vr_j1 = 0; double Vr_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+1].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
		Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0]; 
		Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 1:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
		P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; 
		Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
		Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j].Vx[0]; Vx_j2 = (*cell)[n][i][j-1].Vx[0]; 
		Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 

		P_j1 = 0;
		Vx_j1 = 0;
		Vr_j1 = 0;
		for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			P_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vr[0];
		}
		break;
	
	// Partial cell, only for first-order calc
	case 3:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; 
		P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; 
		Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i][j].Vx[0]; Vx_i2 = (*cell)[n][i-1][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i][j].Vr[0]; Vr_i2 = (*cell)[n][i-1][j].Vr[0]; 
		Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0]; 
		Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
		
		P_i1 = 0;
		Vx_i1 = 0;
		Vr_i1 = 0;
		for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			P_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vr[0];
		}
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
		P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+1].P[0]; 
		Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
		Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0]; 
		Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; 
		P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; 
		Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i][j].Vx[0]; Vx_i2 = (*cell)[n][i-1][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i][j].Vr[0]; Vr_i2 = (*cell)[n][i-1][j].Vr[0]; 
		Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j].Vx[0]; Vx_j2 = (*cell)[n][i][j-1].Vx[0]; 
		Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j].Vr[0]; Vr_j2 = (*cell)[n][i][j-1].Vr[0]; 
		
		P_i1 = 0;
		Vx_i1 = 0;
		Vr_i1 = 0;
		for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			P_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vr[0];
		}
		
		P_j1 = 0;
		Vx_j1 = 0;
		Vr_j1 = 0;
		for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			P_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].P[0];
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].Vr[0];
		}
	    break;
	
	
	// Top and left borders closed
	case 15:
		P_i_2 = (*cell)[n][i+1][j].P[0]; P_i_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; // Mirror i-2 -> i+1, i-1 -> self
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; // Mirror j+2 -> j-1, j+1 -> self
	    Vx_i_2 = -(*cell)[n][i+1][j].Vx[0]; Vx_i_1 = -(*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i+1][j].Vr[0]; Vr_i_1 = (*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j].Vx[0]; Vx_j2 = (*cell)[n][i][j-1].Vx[0];
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = -(*cell)[n][i][j].Vr[0]; Vr_j2 = -(*cell)[n][i][j-1].Vr[0];
	    break;
	    
	// Left border closed
	case 17:
		P_i_2 = (*cell)[n][i+1][j].P[0]; P_i_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; // Mirror i-2 -> i+1, i-1 -> self 
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+1].P[0]; 
	    Vx_i_2 = -(*cell)[n][i+1][j].Vx[0]; Vx_i_1 = -(*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i+1][j].Vr[0]; Vr_i_1 = (*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Bottom and left borders closed
	case 16:
		P_i_2 = (*cell)[n][i+1][j].P[0]; P_i_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; // Mirror i-2 -> i+1, i-1 -> self 
	    P_j_2 = (*cell)[n][i][j+1].P[0]; P_j_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; // Mirror j-2 -> j+1, j-1 -> self 
	    Vx_i_2 = -(*cell)[n][i+1][j].Vx[0]; Vx_i_1 = -(*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i+1][j].Vr[0]; Vr_i_1 = (*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j+1].Vx[0]; Vx_j_1 = (*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0];
	    Vr_j_2 = -(*cell)[n][i][j+1].Vr[0]; Vr_j_1 = -(*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0];
	    break;
	
	// Top border closed
	case 13:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; // Mirror j+2 -> j-1, j+1 -> self
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j].Vx[0]; Vx_j2 = (*cell)[n][i][j-1].Vx[0];
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = -(*cell)[n][i][j].Vr[0]; Vr_j2 = -(*cell)[n][i][j-1].Vr[0];
	    break;
	
	// Bottom border closed
	case 14:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i+1][j].P[0]; P_i2 = (*cell)[n][i+2][j].P[0]; 
	    P_j_2 = (*cell)[n][i][j+1].P[0]; P_j_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; // Mirror j-2 -> j+1, j-1 -> self 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = (*cell)[n][i+1][j].Vx[0]; Vx_i2 = (*cell)[n][i+2][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i+1][j].Vr[0]; Vr_i2 = (*cell)[n][i+2][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j+1].Vx[0]; Vx_j_1 = (*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0];
	    Vr_j_2 = -(*cell)[n][i][j+1].Vr[0]; Vr_j_1 = -(*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0];
	    break;
	
	// Right border closed
	case 19:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; // Mirror i+2 -> i-1, i+1 -> self
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+1].P[0]; 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = -(*cell)[n][i][j].Vx[0]; Vx_i2 = -(*cell)[n][i-1][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i][j].Vr[0]; Vr_i2 = (*cell)[n][i-1][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0]; 
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0]; 
	    break;
	
	// Top and right borders closed
	case 20:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; // Mirror i+2 -> i-1, i+1 -> self 
	    P_j_2 = (*cell)[n][i][j-2].P[0]; P_j_1 = (*cell)[n][i][j-1].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j].P[0]; P_j2 = (*cell)[n][i][j-1].P[0]; // Mirror j+2 -> j-1, j+1 -> self
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = -(*cell)[n][i][j].Vx[0]; Vx_i2 = -(*cell)[n][i-1][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i][j].Vr[0]; Vr_i2 = (*cell)[n][i-1][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j-2].Vx[0]; Vx_j_1 = (*cell)[n][i][j-1].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j].Vx[0]; Vx_j2 = (*cell)[n][i][j-1].Vx[0];
	    Vr_j_2 = (*cell)[n][i][j-2].Vr[0]; Vr_j_1 = (*cell)[n][i][j-1].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = -(*cell)[n][i][j].Vr[0]; Vr_j2 = -(*cell)[n][i][j-1].Vr[0];
	    break;
	
	// Bottom and right borders closed
	case 21:
		P_i_2 = (*cell)[n][i-2][j].P[0]; P_i_1 = (*cell)[n][i-1][j].P[0]; P = (*cell)[n][i][j].P[0]; P_i1 = (*cell)[n][i][j].P[0]; P_i2 = (*cell)[n][i-1][j].P[0]; // Mirror i+2 -> i-1, i+1 -> self 
	    P_j_2 = (*cell)[n][i][j+1].P[0]; P_j_1 = (*cell)[n][i][j].P[0]; P = (*cell)[n][i][j].P[0]; P_j1 = (*cell)[n][i][j+1].P[0]; P_j2 = (*cell)[n][i][j+2].P[0]; // Mirror j-2 -> j+1, j-1 -> self 
	    Vx_i_2 = (*cell)[n][i-2][j].Vx[0]; Vx_i_1 = (*cell)[n][i-1][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_i1 = -(*cell)[n][i][j].Vx[0]; Vx_i2 = -(*cell)[n][i-1][j].Vx[0]; 
		Vr_i_2 = (*cell)[n][i-2][j].Vr[0]; Vr_i_1 = (*cell)[n][i-1][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_i1 = (*cell)[n][i][j].Vr[0]; Vr_i2 = (*cell)[n][i-1][j].Vr[0]; 
	    Vx_j_2 = (*cell)[n][i][j+1].Vx[0]; Vx_j_1 = (*cell)[n][i][j].Vx[0]; Vx = (*cell)[n][i][j].Vx[0]; Vx_j1 = (*cell)[n][i][j+1].Vx[0]; Vx_j2 = (*cell)[n][i][j+2].Vx[0];
	    Vr_j_2 = -(*cell)[n][i][j+1].Vr[0]; Vr_j_1 = -(*cell)[n][i][j].Vr[0]; Vr = (*cell)[n][i][j].Vr[0]; Vr_j1 = (*cell)[n][i][j+1].Vr[0]; Vr_j2 = (*cell)[n][i][j+2].Vr[0];
	    break;
	
	default:
	    break;
    }
	    
    /** First order **/
	double P_i12 = (P_i1 + P)/2;
	double P_i_12 = (P_i_1 + P)/2;
	double P_j12 = (P_j1 + P)/2;
	double P_j_12 = (P_j_1 + P)/2;
	double result = (*cell)[n][i][j].e - 
		(
			(P_i12*(Vx_i1 + Vx)/2 - P_i_12*(Vx_i_1 + Vx)/2) / dx * fmax((*cell)[n][i][j].A[1],(*cell)[n][i][j].A[2]) +
			(fabs(j-axis_j)*P_j12*(Vr_j1 + Vr)/2 - fabs(j-axis_j-1)*P_j_12*(Vr_j_1 + Vr)/2) / ((fabs(j-axis_j-0.5))*dr) *	fmax((*cell)[n][i][j].A[3],(*cell)[n][i][j].A[4])
		) * dt / ((*cell)[n][i][j].A[0] * (*cell)[n][i][j].rho)
	/** without fmax **/
		//~ (
		//~ ((*cell)[n][i][j].P[2]*(*cell)[n][i][j].Vx[2]*(*cell)[n][i][j].A[2] - (*cell)[n][i][j].P[1]*(*cell)[n][i][j].Vx[1]*(*cell)[n][i][j].A[1]) / dx +
		//~ (j*(*cell)[n][i][j].P[4]*(*cell)[n][i][j].Vr[4]*(*cell)[n][i][j].A[4] - (j-1)*(*cell)[n][i][j].P[3]*(*cell)[n][i][j].Vr[3]*(*cell)[n][i][j].A[3]) / ((j-1/2)*dr)
		//~ ) * dt / ((*cell)[n][i][j].A[0] * (*cell)[n][i][j].rho)
	/** If powder present **/
		+ (*cell)[n][i][j].P[0]/((k-1)*I_k) * (1 + lambda*(*cell)[n][i][j].final_z) * f*kappa*dt
	/** End here **/
		;
		
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

    /** With Navier-Stocks **/
    //~ // Gas internal energy
    //~ double I1 = (*cell)[n][i][j].e - (pow((*cell)[n][i][j].Vx[0],2) + pow((*cell)[n][i][j].Vr[0],2))/2;
    //~ double I2 = (*cell)[n][i-1][j].e - (pow((*cell)[n][i-1][j].Vx[0],2) + pow((*cell)[n][i-1][j].Vr[0],2))/2;
    //~ double I3 = (*cell)[n][i+1][j].e - (pow((*cell)[n][i+1][j].Vx[0],2) + pow((*cell)[n][i+1][j].Vr[0],2))/2;
    //~ double I4 = (*cell)[n][i][j-1].e - (pow((*cell)[n][i][j-1].Vx[0],2) + pow((*cell)[n][i][j-1].Vr[0],2))/2;
    //~ double I5 = (*cell)[n][i][j+1].e - (pow((*cell)[n][i][j+1].Vx[0],2) + pow((*cell)[n][i][j+1].Vr[0],2))/2;
    //~ 
    //~ double result = (*cell)[n][i][j].e + dt / ((*cell)[n][i][j].rho) *
		//~ (
			//~ // Pressure
			//~ 1/(2*dx) * (
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
			//~ ) + 
			//~ // Viscosity
			//~ // X axis
			//~ gasMu * (gasA + 2) / (2 * pow(dx,2)) * (
				//~ pow((*cell)[n][i+1][j].Vx[0], 2) +
				//~ pow((*cell)[n][i-1][j].Vx[0], 2) -
				//~ 2*pow((*cell)[n][i][j].Vx[0], 2)
			//~ ) + gasA * gasMu / (4*dx*dr) * (
				//~ (*cell)[n][i+1][j].Vx[0] * (
					//~ (*cell)[n][i+1][j+1].Vr[0] - (*cell)[n][i+1][j-1].Vr[0]
				//~ ) - (*cell)[n][i-1][j].Vx[0] * (
					//~ (*cell)[n][i-1][j+1].Vr[0] - (*cell)[n][i-1][j-1].Vr[0]
				//~ )
			//~ ) + gasMu / (pow(dx,2)) * (
				//~ pow((*cell)[n][i+1][j].Vr[0], 2) +
				//~ pow((*cell)[n][i-1][j].Vr[0], 2) -
				//~ 2*pow((*cell)[n][i][j].Vr[0], 2)
			//~ ) + gasMu / (4*dx*dr) * (
				//~ (*cell)[n][i+1][j].Vr[0] * (
					//~ (*cell)[n][i+1][j+1].Vx[0] - (*cell)[n][i+1][j-1].Vx[0]
				//~ ) - (*cell)[n][i-1][j].Vr[0] * (
					//~ (*cell)[n][i-1][j+1].Vx[0] - (*cell)[n][i-1][j-1].Vx[0]
				//~ )
			//~ ) + gasB * gasMu / pow(dx,2) * (
				//~ I3 + I2 - 2*I1
			//~ ) +
			//~ // R axis
			//~ gasMu * (gasA + 2) / (2 * pow(dr,2)) * (
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
			//~ )
	//~ )
	/** If powder present **/
	//~ + (*cell)[n][i][j].P[0]/((k-1)*I_k) * (1 + lambda*(*cell)[n][i][j].final_z) * f*kappa*dt
	//~ ;
	
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
	if (fabs(result-(*cell)[n][i][j].e) < pow(10,-15)) result = (*cell)[n][i][j].e;
	if (result == result) { // if not NaN
		return result;
	} else {
		cout << "bar_e is NaN" << endl
			<< "i = " << i << endl
			<< "j = " << j << endl
			<< "Vr[3] = " << (*cell)[n][i][j].Vr[3] << endl
			<< "Vr[4] = " << (*cell)[n][i][j].Vr[4] << endl
			<< "P[3] = " << (*cell)[n][i][j].P[3] << endl
			<< "P[4] = " << (*cell)[n][i][j].P[4] << endl
			<< "Second brackets = " <<((j-axis_j)*(*cell)[n][i][j].P[4]*(*cell)[n][i][j].Vr[4] - (j-axis_j-1)*(*cell)[n][i][j].P[3]*(*cell)[n][i][j].Vr[3]) / ((fabs(j-axis_j-0.5))*dr) *
				fmax((*cell)[n][i][j].A[3],(*cell)[n][i][j].A[4]) << endl
			<< "(fabs(j-axis_j-0.5))*dr = " << (fabs(j-axis_j-0.5))*dr << endl
			<< "fmax((*cell)[n][i][j].A[3],(*cell)[n][i][j].A[4]) = " << fmax((*cell)[n][i][j].A[3],(*cell)[n][i][j].A[4]) << endl << endl;
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
    if (fabs(cell->dM[1]) > pow(10,-15) || fabs(cell->dM[2]) > pow(10,-15) || fabs(cell->dM[3]) > pow(10,-15) || fabs(cell->dM[4]) > pow(10,-15)) {
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
        //+
        //cell->m*dt/(cell->A[0]);
    return result;
}

/* Lagrangian stage mass calculation */
void lagrange_mass(double array[21], cell2d * cell, int i, int j, int n,
        double dx, double dr, double dt) {
			
    /****************************************************************
     * Additional zero to have the same length as A					*
     * first empty + dM[4] + dMVx[4] + dMVr[4] + dME[4] + D[4]		*
     * **************************************************************/
    for (int iter = 0; iter < 21; iter++) array[iter] = 0;
    
    /* Small hack to calc i_sn correctly */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }

    /* I'll use the following naming rule Vx_i1 = Vx[i+1], Vx_i_1 = Vx[i-1] */
    double Vx_i_2 = 0;	double Vr_i_2 = 0;	double rho_i_2 = 0;
    double Vx_i_1 = 0;	double Vr_i_1 = 0;	double rho_i_1 = 0;
    double Vx = 0;	double Vr = 0;		double rho = 0;
    double Vx_i1 = 0;	double Vr_i1 = 0;	double rho_i1 = 0;
    double Vx_i2 = 0;	double Vr_i2 = 0;	double rho_i2 = 0;
    double Vx_j_2 = 0;	double Vr_j_2 = 0;	double rho_j_2 = 0;
    double Vx_j_1 = 0;	double Vr_j_1 = 0;	double rho_j_1 = 0;
    double Vx_j1 = 0;	double Vr_j1 = 0;	double rho_j1 = 0;
    double Vx_j2 = 0;	double Vr_j2 = 0;	double rho_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Partial cell, only for first-order calc
	case 1:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho;
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho;

	    rho_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho;
	    Vx_i2 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;

	    rho_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho;
	    Vx_i2 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho;
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho;
		
	    rho_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	
	    rho_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	
	// Top and left borders closed
	case 15:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	    
	// Left border closed
	case 17:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Bottom and left borders closed
	case 16:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Top border closed
	case 13:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom border closed
	case 14:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Right border closed
	case 19:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0]; 	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Top and right borders closed
	case 20:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom and right borders closed
	case 21:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	default:
	    break;
    }
	
	/** First order **/
	bool ruleVx1 = Vx + Vx_i_1 > 0 ? true : false;
	bool ruleVx2 = Vx + Vx_i1 > 0 ? true : false;
	bool ruleVr1 = Vr + Vr_j_1 > 0 ? true : false;
	bool ruleVr2 = Vr + Vr_j1 > 0 ? true : false;
	
	if (ruleVx1) {
			array[1] = (fabs(j-axis_j-0.5))*(*cell)[n][i-1][j].A[2] * rho_i_1 * (Vx + Vx_i_1)/2 * pow(dr,2) * dt;
	} else {
			array[1] = (fabs(j-axis_j-0.5))*(*cell)[n][i][j].A[1] * rho * (Vx + Vx_i_1)/2 * pow(dr,2) * dt;
		}
	if (ruleVx2) {
			array[2] = (fabs(j-axis_j-0.5))*(*cell)[n][i][j].A[2] * rho * (Vx + Vx_i1)/2 * pow(dr,2) * dt;
		} else {
			array[2] = (fabs(j-axis_j-0.5))*(*cell)[n][i+1][j].A[1] * rho_i1 * (Vx + Vx_i1)/2 * pow(dr,2) * dt;
		}
	if (ruleVr1) {
			array[3] = (*cell)[n][i][j-1].A[4] * rho_j_1 * (Vr + Vr_j_1)/2 * dx * dt;
		} else {
			array[3] = (*cell)[n][i][j].A[3] * rho * (Vr + Vr_j_1)/2 * dx * dt;
		}
	if (ruleVr2) {
			array[4] = (*cell)[n][i][j].A[4] * rho * (Vr + Vr_j1)/2 * dx * dt;
		} else {
			array[4] = (*cell)[n][i][j+1].A[3] * rho_j1 * (Vr + Vr_j1)/2 * dx * dt;
		}
	array[5] = ruleVx1 ? 1 : 0;
	array[6] = ruleVx2 ? 0 : 1;
	array[7] = ruleVr1 ? 1 : 0;
	array[8] = ruleVr2 ? 0 : 1;
	
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
            //~ (rho - (rho_i1 - rho_i_1)/4) * (*cell)[n][i-1][j].A[2] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ } 
    //~ // From left to right
    //~ if (ruleVx11 == ruleVx13 && ruleVx11 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i_2) > pow(10,-15) || fabs(Vx_i_1) > pow(10,-15))) {
	//~ array[1] = (Vx_i_1 + (Vx - Vx_i_2)/4) *
            //~ (rho_i_1 + (rho - rho_i_2)/4) * (*cell)[n][i][j].A[1] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ }
    //~ // dM(i+1/2,j)
    //~ // From left to right
    //~ if (ruleVx21 == ruleVx22 && ruleVx21 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i1) > pow(10,-15) || fabs(Vx_i_1) > pow(10,-15))) {
        //~ array[2] = (Vx + (Vx_i1 - Vx_i_1)/4) *
            //~ (rho + (rho_i1 - rho_i_1)/4) * (*cell)[n][i][j].A[2] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ } 
    //~ // From right to left
    //~ if (ruleVx21 == ruleVx23 && !ruleVx21 && (fabs(Vx) > pow(10,-15) || fabs(Vx_i1) > pow(10,-15) || fabs(Vx_i2) > pow(10,-15))) {
        //~ array[2] = (Vx_i1 - (Vx_i2 - Vx)/4) *
            //~ (rho_i1 - (rho_i2 - rho)/4) * (*cell)[n][i+1][j].A[1] * fabs(j-axis_j-0.5)*pow(dr,2) * dt;
    //~ }
    //~ // dM(i,j-1/2)
    //~ if (ruleVr11 == ruleVr13 && ruleVr11 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j_1) > pow(10,-15) || fabs(Vr_j_2) > pow(10,-15))) {
        //~ array[3] = (Vr_j_1 - (Vr_j_2 - Vr)/4) *
            //~ (rho_j_1 - (rho_j_2 - rho)/4) * (*cell)[n][i][j-1].A[4] * dx * dt;
    //~ }
    //~ if (ruleVr11 == ruleVr12 && !ruleVr11 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j1) > pow(10,-15) || fabs(Vr_j_1) > pow(10,-15))) {
	//~ array[3] = (Vr + (Vr_j_1 - Vr_j1)/4) *
            //~ (rho + (rho_j_1 - rho_j1)/4) * (*cell)[n][i][j].A[3] * dx * dt;
    //~ }
    //~ // dM(i,j+1/2)
    //~ if (ruleVr21 == ruleVr22 && ruleVr21 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j1) > pow(10,-15) || fabs(Vr_j_1) > pow(10,-15))) {
            //~ array[4] = (Vr + (Vr_j1 - Vr_j_1)/4) *
                    //~ (rho + (rho_j1 - rho_j_1)/4) * (*cell)[n][i][j].A[4] * dx * dt;
	//~ }
    //~ if (ruleVr21 == ruleVr23 && !ruleVr21 && (fabs(Vr) > pow(10,-15) || fabs(Vr_j1) > pow(10,-15) || fabs(Vr_j2) > pow(10,-15))) {
            //~ array[4] = (Vr_j1 - (Vr_j2 - Vr)/4) *
                    //~ (rho_j1 - (rho_j2 - rho)/4) * (*cell)[n][i][j+1].A[3] * dx * dt;
    //~ }
    //~ 
    //~ if ((*cell)[n][i][j].A[1] == 0) array[1] = 0;
	//~ if ((*cell)[n][i][j].A[2] == 0) array[2] = 0;
	//~ if ((*cell)[n][i][j].A[3] == 0) array[3] = 0;
	//~ if ((*cell)[n][i][j].A[4] == 0) array[4] = 0;
    //~ 
    //~ array[5] = array[1] > 0 ? 1 : 0;
    //~ array[6] = array[2] > 0 ? 0 : 1;
    //~ array[7] = array[3] > 0 ? 1 : 0;
    //~ array[8] = array[4] > 0 ? 0 : 1;
    
    /** DEBUG **/
    if (i == 248 && j == 5) {
	cout << "i = 248, j = 5, cell type = " << type << ", Vx{i-1,i,i+1} = {" << Vx_i_1 << ", " << Vx << ", " << Vx_i1 << "}" << endl;
    }
    if (i == 249 && j == 5) {
	cout << "i = 249, j = 5, cell type = " << type << ", Vx{i-1,i,i+1} = {" << Vx_i_1 << ", " << Vx << ", " << Vx_i1 << "}" << endl;
    }
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

double new_final_z (cell2d * cell, int i, int j, int n,
        double dx, double dr, double dt) {
		/* Small hack to calc i_sn correctly */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }

    /* I'll use the following naming rule Vx_i1 = Vx[i+1], Vx_i_1 = Vx[i-1] */
    double Vx_i_2 = 0;	double Vr_i_2 = 0;	double bar_z_i_2 = 0;
    double Vx_i_1 = 0;	double Vr_i_1 = 0;	double bar_z_i_1 = 0;
    double Vx = 0;	double Vr = 0;		double bar_z = 0;
    double Vx_i1 = 0;	double Vr_i1 = 0;	double bar_z_i1 = 0;
    double Vx_i2 = 0;	double Vr_i2 = 0;	double bar_z_i2 = 0;
    double Vx_j_2 = 0;	double Vr_j_2 = 0;	double bar_z_j_2 = 0;
    double Vx_j_1 = 0;	double Vr_j_1 = 0;	double bar_z_j_1 = 0;
    double Vx_j1 = 0;	double Vr_j1 = 0;	double bar_z_j1 = 0;
    double Vx_j2 = 0;	double Vr_j2 = 0;	double bar_z_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	// Partial cell, only for first-order calc
	case 1:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;

	    bar_z_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
		bar_z_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_z;
		Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;

	    bar_z_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
		bar_z_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_z;
		Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
		
	    bar_z_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
		bar_z_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_z;
		Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	
	    bar_z_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
		bar_z_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_z;
		Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	
	// Top and left borders closed
	case 15:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i+1][j].bar_z; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i][j].bar_z; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j].bar_z; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j-1].bar_z; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	    
	// Left border closed
	case 17:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i+1][j].bar_z; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i][j].bar_z; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	// Bottom and left borders closed
	case 16:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i+1][j].bar_z; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i][j].bar_z; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_z_j_2 = (*cell)[n][i][j+1].bar_z; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j].bar_z; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	// Top border closed
	case 13:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j].bar_z; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j-1].bar_z; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom border closed
	case 14:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i+1][j].bar_z;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i+2][j].bar_z;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_z_j_2 = (*cell)[n][i][j+1].bar_z; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j].bar_z; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	// Right border closed
	case 19:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0]; 	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i][j].bar_z; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i-1][j].bar_z; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	// Top and right borders closed
	case 20:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i][j].bar_z; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i-1][j].bar_z; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_z_j_2 = (*cell)[n][i][j-2].bar_z;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j-1].bar_z;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j].bar_z; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j-1].bar_z; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom and right borders closed
	case 21:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_z_i_2 = (*cell)[n][i-2][j].bar_z;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i_1 = (*cell)[n][i-1][j].bar_z;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_z = (*cell)[n][i][j].bar_z;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_z_i1 = (*cell)[n][i][j].bar_z; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_z_i2 = (*cell)[n][i-1][j].bar_z; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_z_j_2 = (*cell)[n][i][j+1].bar_z; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_z_j_1 = (*cell)[n][i][j].bar_z; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_z_j1 = (*cell)[n][i][j+1].bar_z;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_z_j2 = (*cell)[n][i][j+2].bar_z;
	    break;
	
	default:
	    break;
    }
	double result = (
        (*cell)[n][i][j].D[1] * bar_z_i_1 * (*cell)[n][i][j].dM[1] +
        (*cell)[n][i][j].D[2] * bar_z_i1 * (*cell)[n][i][j].dM[2] +
        (*cell)[n][i][j].D[3] * bar_z_j_1 * (*cell)[n][i][j].dM[3] +
        (*cell)[n][i][j].D[4] * bar_z_j1 * (*cell)[n][i][j].dM[4] +
        bar_z * (
			(*cell)[n][i][j].rho * (*cell)[n][i][j].A[0] * dx * (fabs(j-axis_j-0.5))*pow(dr,2) -
            (1-(*cell)[n][i][j].D[1]) * (*cell)[n][i][j].dM[1] -
            (1-(*cell)[n][i][j].D[2]) * (*cell)[n][i][j].dM[2] -
            (1-(*cell)[n][i][j].D[3]) * (*cell)[n][i][j].dM[3] -
            (1-(*cell)[n][i][j].D[4]) * (*cell)[n][i][j].dM[4]
        )
    ) / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2));
    return result;
}

double new_final_psi (cell2d * cell, int i, int j, int n,
        double dx, double dr, double dt) {
		/* Small hack to calc i_sn correctly */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }

    /* I'll use the following naming rule Vx_i1 = Vx[i+1], Vx_i_1 = Vx[i-1] */
    double Vx_i_2 = 0;	double Vr_i_2 = 0;	double bar_psi_i_2 = 0;
    double Vx_i_1 = 0;	double Vr_i_1 = 0;	double bar_psi_i_1 = 0;
    double Vx = 0;	double Vr = 0;		double bar_psi = 0;
    double Vx_i1 = 0;	double Vr_i1 = 0;	double bar_psi_i1 = 0;
    double Vx_i2 = 0;	double Vr_i2 = 0;	double bar_psi_i2 = 0;
    double Vx_j_2 = 0;	double Vr_j_2 = 0;	double bar_psi_j_2 = 0;
    double Vx_j_1 = 0;	double Vr_j_1 = 0;	double bar_psi_j_1 = 0;
    double Vx_j1 = 0;	double Vr_j1 = 0;	double bar_psi_j1 = 0;
    double Vx_j2 = 0;	double Vr_j2 = 0;	double bar_psi_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	
	// Partial cell, only for first-order calc
	case 1:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;

	    bar_psi_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
		bar_psi_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_psi;
		Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;

	    bar_psi_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
		bar_psi_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_psi;
		Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
		
	    bar_psi_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
		bar_psi_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_psi;
		Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	
	    bar_psi_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
		Int2D weightCell;
		double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
		weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
		weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
		bar_psi_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_psi;
		Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
		Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	
	
	// Top and left borders closed
	case 15:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i+1][j].bar_psi; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i][j].bar_psi; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j].bar_psi; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j-1].bar_psi; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	    
	// Left border closed
	case 17:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i+1][j].bar_psi; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i][j].bar_psi; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	// Bottom and left borders closed
	case 16:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i+1][j].bar_psi; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i][j].bar_psi; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_psi_j_2 = (*cell)[n][i][j+1].bar_psi; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j].bar_psi; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	// Top border closed
	case 13:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j].bar_psi; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j-1].bar_psi; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom border closed
	case 14:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i+1][j].bar_psi;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i+2][j].bar_psi;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_psi_j_2 = (*cell)[n][i][j+1].bar_psi; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j].bar_psi; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	// Right border closed
	case 19:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i][j].bar_psi; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i-1][j].bar_psi; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	// Top and right borders closed
	case 20:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i][j].bar_psi; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i-1][j].bar_psi; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_psi_j_2 = (*cell)[n][i][j-2].bar_psi;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j-1].bar_psi;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j].bar_psi; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j-1].bar_psi; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom and right borders closed
	case 21:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_psi_i_2 = (*cell)[n][i-2][j].bar_psi;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i_1 = (*cell)[n][i-1][j].bar_psi;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_psi = (*cell)[n][i][j].bar_psi;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_psi_i1 = (*cell)[n][i][j].bar_psi; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_psi_i2 = (*cell)[n][i-1][j].bar_psi; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_psi_j_2 = (*cell)[n][i][j+1].bar_psi; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_psi_j_1 = (*cell)[n][i][j].bar_psi; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_psi_j1 = (*cell)[n][i][j+1].bar_psi;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_psi_j2 = (*cell)[n][i][j+2].bar_psi;
	    break;
	
	default:
	    break;
    }
	double result = (
        (*cell)[n][i][j].D[1] * bar_psi_i_1 * (*cell)[n][i][j].dM[1] +
        (*cell)[n][i][j].D[2] * bar_psi_i1 * (*cell)[n][i][j].dM[2] +
        (*cell)[n][i][j].D[3] * bar_psi_j_1 * (*cell)[n][i][j].dM[3] +
        (*cell)[n][i][j].D[4] * bar_psi_j1 * (*cell)[n][i][j].dM[4] +
        bar_psi * (
        	(*cell)[n][i][j].rho * (*cell)[n][i][j].A[0] * dx * (fabs(j-axis_j-0.5))*pow(dr,2) -
            (1-(*cell)[n][i][j].D[1]) * (*cell)[n][i][j].dM[1] -
            (1-(*cell)[n][i][j].D[2]) * (*cell)[n][i][j].dM[2] -
            (1-(*cell)[n][i][j].D[3]) * (*cell)[n][i][j].dM[3] -
            (1-(*cell)[n][i][j].D[4]) * (*cell)[n][i][j].dM[4]
        )
    ) / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2));
    return result;
}








/* Final stage P calculation */
double final_calc_p(gasCell * prevCell, gasCell * cell) {
	double result = (cell->e - (pow(cell->Vx[0],2)+pow(cell->Vr[0],2))/2 ) * (k-1) /
		(
			//~ 1/cell->rho // No powder present
			1/cell->rho - (1 - cell->final_psi)/delta - alpha_k * cell->final_psi // BMSTU var
			//~ alpha_k // Abel (Dupre) equation
		)
		//~ - (pow(cell->Vx[0],2)+pow(cell->Vr[0],2))/2 // From Ershov, UDK 519.6:532.6, #775, 2007, str. 159-173
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
double final_calc_Vx(cell2d * cell, int i, int j, int n,
        double dx, double dr, double dt) {
	/* Small hack to calc i_sn correctly */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }

    /* I'll use the following naming rule Vx_i1 = Vx[i+1], Vx_i_1 = Vx[i-1] */
    double Vx_i_2 = 0;	double Vr_i_2 = 0;	double rho_i_2 = 0;
    double Vx_i_1 = 0;	double Vr_i_1 = 0;	double rho_i_1 = 0;
    double Vx = 0;	double Vr = 0;		double rho = 0;
    double Vx_i1 = 0;	double Vr_i1 = 0;	double rho_i1 = 0;
    double Vx_i2 = 0;	double Vr_i2 = 0;	double rho_i2 = 0;
    double Vx_j_2 = 0;	double Vr_j_2 = 0;	double rho_j_2 = 0;
    double Vx_j_1 = 0;	double Vr_j_1 = 0;	double rho_j_1 = 0;
    double Vx_j1 = 0;	double Vr_j1 = 0;	double rho_j1 = 0;
    double Vx_j2 = 0;	double Vr_j2 = 0;	double rho_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	
	// Partial cell, only for first-order calc
	case 1:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;

	    rho_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;

	    rho_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
		
	    rho_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	
	    rho_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	
	
	
	// Top and left borders closed
	case 15:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	    
	// Left border closed
	case 17:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Bottom and left borders closed
	case 16:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Top border closed
	case 13:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom border closed
	case 14:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Right border closed
	case 19:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0]; 	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Top and right borders closed
	case 20:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom and right borders closed
	case 21:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	default:
	    break;
    }
    /* With m */
    double result = (*cell)[n][i][j].bar_Vx[0] * (*cell)[n][i][j].rho / (*cell)[n+1][i][j].rho 
    + 
    (
	    (*cell)[n][i][j].D[1] * Vx_i_1 * (*cell)[n][i][j].dM[1] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[2] * Vx_i1 * (*cell)[n][i][j].dM[2] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[3] * Vx_j_1 * (*cell)[n][i][j].dM[3] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[4] * Vx_j1 * (*cell)[n][i][j].dM[4] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[1]) * Vx * (*cell)[n][i][j].dM[1] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[2]) * Vx * (*cell)[n][i][j].dM[2] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[3]) * Vx * (*cell)[n][i][j].dM[3] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[4]) * Vx * (*cell)[n][i][j].dM[4] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
    )
	;
    return result;
}

/* Final stage Vr calculation */
double final_calc_Vr(cell2d * cell, int i, int j, int n,
        double dx, double dr, double dt) {
	/* Small hack to calc i_sn correctly */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }

    /* I'll use the following naming rule Vx_i1 = Vx[i+1], Vx_i_1 = Vx[i-1] */
    double Vx_i_2 = 0;	double Vr_i_2 = 0;	double rho_i_2 = 0;
    double Vx_i_1 = 0;	double Vr_i_1 = 0;	double rho_i_1 = 0;
    double Vx = 0;	double Vr = 0;		double rho = 0;
    double Vx_i1 = 0;	double Vr_i1 = 0;	double rho_i1 = 0;
    double Vx_i2 = 0;	double Vr_i2 = 0;	double rho_i2 = 0;
    double Vx_j_2 = 0;	double Vr_j_2 = 0;	double rho_j_2 = 0;
    double Vx_j_1 = 0;	double Vr_j_1 = 0;	double rho_j_1 = 0;
    double Vx_j1 = 0;	double Vr_j1 = 0;	double rho_j1 = 0;
    double Vx_j2 = 0;	double Vr_j2 = 0;	double rho_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	    
	// Partial cell, only for first-order calc
	case 1:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;

	    rho_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;

	    rho_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
		
	    rho_i1 = 0;
	    Vx_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	
	    rho_j1 = 0;
	    Vx_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vx_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vx[0];
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	
	
	// Top and left borders closed
	case 15:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	    
	// Left border closed
	case 17:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Bottom and left borders closed
	case 16:
	    Vx_i_2 = -(*cell)[n][i+1][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i+1][j].rho; // Mirror -2 -> +1, Vx = -bar_Vx[i+1]
	    Vx_i_1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vx = -bar_Vx[i]
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Top border closed
	case 13:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom border closed
	case 14:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = (*cell)[n][i+1][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i+1][j].rho;
	    Vx_i2 = (*cell)[n][i+2][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i+2][j].rho;
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Right border closed
	case 19:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0]; 	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	// Top and right borders closed
	case 20:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0];		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j-2].bar_Vx[0];	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	rho_j_2 = (*cell)[n][i][j-2].rho;
	    Vx_j_1 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j-1].rho;
	    Vx_j1 = (*cell)[n][i][j].bar_Vx[0];		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    Vx_j2 = (*cell)[n][i][j-1].bar_Vx[0];	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j-1].rho; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom and right borders closed
	case 21:
	    Vx_i_2 = (*cell)[n][i-2][j].bar_Vx[0];	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	rho_i_2 = (*cell)[n][i-2][j].rho;
	    Vx_i_1 = (*cell)[n][i-1][j].bar_Vx[0];	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i_1 = (*cell)[n][i-1][j].rho;
	    Vx = (*cell)[n][i][j].bar_Vx[0]; 		Vr = (*cell)[n][i][j].bar_Vr[0];	rho = (*cell)[n][i][j].rho;
	    Vx_i1 = -(*cell)[n][i][j].bar_Vx[0];	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	rho_i1 = (*cell)[n][i][j].rho; // Mirror +1 -> self, Vx = -bar_Vx[i]
	    Vx_i2 = -(*cell)[n][i-1][j].bar_Vx[0];	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	rho_i2 = (*cell)[n][i-1][j].rho; // Mirror +2 -> -1, Vx = -bar_Vx[i-1]
	    Vx_j_2 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	rho_j_2 = (*cell)[n][i][j+1].rho; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    Vx_j_1 = (*cell)[n][i][j].bar_Vx[0];	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	rho_j_1 = (*cell)[n][i][j].rho; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    Vx_j1 = (*cell)[n][i][j+1].bar_Vx[0];	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	rho_j1 = (*cell)[n][i][j+1].rho;
	    Vx_j2 = (*cell)[n][i][j+2].bar_Vx[0];	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	rho_j2 = (*cell)[n][i][j+2].rho;
	    break;
	
	default:
	    break;
    }
    double result = (*cell)[n][i][j].bar_Vr[0] * (*cell)[n][i][j].rho / (*cell)[n+1][i][j].rho 
    + 
    (
	    (*cell)[n][i][j].D[1] * Vr_i_1 * (*cell)[n][i][j].dM[1] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[2] * Vr_i1 * (*cell)[n][i][j].dM[2] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[3] * Vr_j_1 * (*cell)[n][i][j].dM[3] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[4] * Vr_j1 * (*cell)[n][i][j].dM[4] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[1]) * Vr * (*cell)[n][i][j].dM[1] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[2]) * Vr * (*cell)[n][i][j].dM[2] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[3]) * Vr * (*cell)[n][i][j].dM[3] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[4]) * Vr * (*cell)[n][i][j].dM[4] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	)
    ;
    return result;
}


/* Final stage E calculation */
double final_calc_e(cell2d * cell, int i, int j, int n,
        double dx, double dr, double dt) {
	/* Small hack to calc i_sn correctly */
    int type = (*cell)[n][i][j].type;
    if (i == i_sn - 1) {
		if (type == 13) {
		    type = 20;
		} else if (type == 14) {
		    type = 21;
		} else if (type == 0) {
		    type = 19;
		}
    }

    /* I'll use the following naming rule Vx_i1 = Vx[i+1], Vx_i_1 = Vx[i-1] */
    double rho_i_2 = 0;	double Vr_i_2 = 0;	double bar_e_i_2 = 0;
    double rho_i_1 = 0;	double Vr_i_1 = 0;	double bar_e_i_1 = 0;
    double rho = 0;	double Vr = 0;		double bar_e = 0;
    double rho_i1 = 0;	double Vr_i1 = 0;	double bar_e_i1 = 0;
    double rho_i2 = 0;	double Vr_i2 = 0;	double bar_e_i2 = 0;
    double rho_j_2 = 0;	double Vr_j_2 = 0;	double bar_e_j_2 = 0;
    double rho_j_1 = 0;	double Vr_j_1 = 0;	double bar_e_j_1 = 0;
    double rho_j1 = 0;	double Vr_j1 = 0;	double bar_e_j1 = 0;
    double rho_j2 = 0;	double Vr_j2 = 0;	double bar_e_j2 = 0;
    
    switch (type) {
	// Free cell
	case 0:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho;		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	// Partial cell, only for first-order calc
	case 1:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho;		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;

	    bar_e_j1 = 0;
	    rho_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			bar_e_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_e;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 3:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho;		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;

	    bar_e_i1 = 0;
	    rho_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			bar_e_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_e;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	// Partial cell, only for first-order calc
	case 8:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho;		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	// Partial cell, only for first-order calc
	case 10:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho;		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
		
	    bar_e_i1 = 0;
	    rho_i1 = 0;
	    Vr_i1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.x.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.x.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.x.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.x.at(idx).j;
			bar_e_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_e;
			rho_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vr_i1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	
	    bar_e_j1 = 0;
	    rho_j1 = 0;
	    Vr_j1 = 0;
	    for (unsigned int idx = 0; idx < (*cell)[n][i][j].weightVector.y.size(); idx++) {
			Int2D weightCell;
			double weight = (*cell)[n][i][j].weightVector.y.at(idx).weight;
			weightCell.i = (*cell)[n][i][j].weightVector.y.at(idx).i;
			weightCell.j = (*cell)[n][i][j].weightVector.y.at(idx).j;
			bar_e_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_e;
			rho_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].rho;
			Vr_j1 += weight*(*cell)[n][weightCell.i][weightCell.j].bar_Vr[0];
	    }
	    break;
	
	
	
	// Top and left borders closed
	case 15:
	    rho_i_2 = -(*cell)[n][i+1][j].rho;	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i+1][j].bar_e; // Mirror -2 -> +1, rho = -bar_rho[i+1]
	    rho_i_1 = -(*cell)[n][i][j].rho;	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i][j].bar_e; // Mirror -1 -> self, rho = -bar_rho[i]
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j].rho;		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j].bar_e; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    rho_j2 = (*cell)[n][i][j-1].rho;	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j-1].bar_e; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	    
	// Left border closed
	case 17:
	    rho_i_2 = -(*cell)[n][i+1][j].rho;	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i+1][j].bar_e; // Mirror -2 -> +1, rho = -bar_rho[i+1]
	    rho_i_1 = -(*cell)[n][i][j].rho;	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i][j].bar_e; // Mirror -1 -> self, rho = -bar_rho[i]
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	// Bottom and left borders closed
	case 16:
	    rho_i_2 = -(*cell)[n][i+1][j].rho;	Vr_i_2 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i+1][j].bar_e; // Mirror -2 -> +1, rho = -bar_rho[i+1]
	    rho_i_1 = -(*cell)[n][i][j].rho;	Vr_i_1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i][j].bar_e; // Mirror -1 -> self, rho = -bar_rho[i]
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j+1].rho;	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_e_j_2 = (*cell)[n][i][j+1].bar_e; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    rho_j_1 = (*cell)[n][i][j].rho;	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j].bar_e; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	// Top border closed
	case 13:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j].rho;		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j].bar_e; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    rho_j2 = (*cell)[n][i][j-1].rho;	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j-1].bar_e; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom border closed
	case 14:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = (*cell)[n][i+1][j].rho;	Vr_i1 = (*cell)[n][i+1][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i+1][j].bar_e;
	    rho_i2 = (*cell)[n][i+2][j].rho;	Vr_i2 = (*cell)[n][i+2][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i+2][j].bar_e;
	    rho_j_2 = (*cell)[n][i][j+1].rho;	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_e_j_2 = (*cell)[n][i][j+1].bar_e; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    rho_j_1 = (*cell)[n][i][j].rho;	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j].bar_e; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	// Right border closed
	case 19:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0]; 	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = -(*cell)[n][i][j].rho;	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i][j].bar_e; // Mirror +1 -> self, rho = -bar_rho[i]
	    rho_i2 = -(*cell)[n][i-1][j].rho;	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i-1][j].bar_e; // Mirror +2 -> -1, rho = -bar_rho[i-1]
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	// Top and right borders closed
	case 20:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho;		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = -(*cell)[n][i][j].rho;	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i][j].bar_e; // Mirror +1 -> self, rho = -bar_rho[i]
	    rho_i2 = -(*cell)[n][i-1][j].rho;	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i-1][j].bar_e; // Mirror +2 -> -1, rho = -bar_rho[i-1]
	    rho_j_2 = (*cell)[n][i][j-2].rho;	Vr_j_2 = (*cell)[n][i][j-2].bar_Vr[0]; 	bar_e_j_2 = (*cell)[n][i][j-2].bar_e;
	    rho_j_1 = (*cell)[n][i][j-1].rho;	Vr_j_1 = (*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j-1].bar_e;
	    rho_j1 = (*cell)[n][i][j].rho;		Vr_j1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j].bar_e; // Mirror +1 -> self, Vr = -bar_Vr[j]
	    rho_j2 = (*cell)[n][i][j-1].rho;	Vr_j2 = -(*cell)[n][i][j-1].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j-1].bar_e; // Mirror +2 -> -1, Vr = -bar_Vr[j-1]
	    break;
	
	// Bottom and right borders closed
	case 21:
	    rho_i_2 = (*cell)[n][i-2][j].rho;	Vr_i_2 = (*cell)[n][i-2][j].bar_Vr[0]; 	bar_e_i_2 = (*cell)[n][i-2][j].bar_e;
	    rho_i_1 = (*cell)[n][i-1][j].rho;	Vr_i_1 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i_1 = (*cell)[n][i-1][j].bar_e;
	    rho = (*cell)[n][i][j].rho; 		Vr = (*cell)[n][i][j].bar_Vr[0];	bar_e = (*cell)[n][i][j].bar_e;
	    rho_i1 = -(*cell)[n][i][j].rho;	Vr_i1 = (*cell)[n][i][j].bar_Vr[0]; 	bar_e_i1 = (*cell)[n][i][j].bar_e; // Mirror +1 -> self, rho = -bar_rho[i]
	    rho_i2 = -(*cell)[n][i-1][j].rho;	Vr_i2 = (*cell)[n][i-1][j].bar_Vr[0]; 	bar_e_i2 = (*cell)[n][i-1][j].bar_e; // Mirror +2 -> -1, rho = -bar_rho[i-1]
	    rho_j_2 = (*cell)[n][i][j+1].rho;	Vr_j_2 = -(*cell)[n][i][j+1].bar_Vr[0];	bar_e_j_2 = (*cell)[n][i][j+1].bar_e; // Mirror -2 -> +1, Vr = -bar_Vr[j+1]
	    rho_j_1 = (*cell)[n][i][j].rho;	Vr_j_1 = -(*cell)[n][i][j].bar_Vr[0]; 	bar_e_j_1 = (*cell)[n][i][j].bar_e; // Mirror -1 -> self, Vr = -bar_Vr[j]
	    rho_j1 = (*cell)[n][i][j+1].rho;	Vr_j1 = (*cell)[n][i][j+1].bar_Vr[0]; 	bar_e_j1 = (*cell)[n][i][j+1].bar_e;
	    rho_j2 = (*cell)[n][i][j+2].rho;	Vr_j2 = (*cell)[n][i][j+2].bar_Vr[0]; 	bar_e_j2 = (*cell)[n][i][j+2].bar_e;
	    break;
	
	default:
	    break;
    }
    double result = (*cell)[n][i][j].bar_e * (*cell)[n][i][j].rho / (*cell)[n+1][i][j].rho 
    + 
    (
	    (*cell)[n][i][j].D[1] * bar_e_i_1 * (*cell)[n][i][j].dM[1] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[2] * bar_e_i1 * (*cell)[n][i][j].dM[2] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[3] * bar_e_j_1 * (*cell)[n][i][j].dM[3] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    +
	    (*cell)[n][i][j].D[4] * bar_e_j1 * (*cell)[n][i][j].dM[4] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[1]) * bar_e * (*cell)[n][i][j].dM[1] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[2]) * bar_e * (*cell)[n][i][j].dM[2] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[3]) * bar_e * (*cell)[n][i][j].dM[3] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	    -
	    (1-(*cell)[n][i][j].D[4]) * bar_e * (*cell)[n][i][j].dM[4] / ((fabs(j-axis_j-0.5)) * (*cell)[n+1][i][j].rho * (*cell)[n][i][j].A[0] * dx * pow(dr,2))
	)
    ;
	if (result < 0) {
		cout << "E < 0" << endl;
		cout << "i = " << i << endl
			<< "j = " << j << endl;
		cout << setiosflags(ios::fixed) << setprecision(16) << 
			"bar_E       = " << (*cell)[n][i][j].bar_e << endl <<
			"lower bar_E = " << (*cell)[n][i][j-1].bar_e << endl <<
			"upper bar_E = " << (*cell)[n][i][j+1].bar_e << endl <<
			"E       = " << (*cell)[n+1][i][j].e << endl <<
			"lower E = " << (*cell)[n+1][i][j-1].e << endl <<
			"upper E = " << (*cell)[n+1][i][j+1].e << endl <<
			"P[0] = " << (*cell)[n][i][j].P[0] << endl <<
			"A[0] = " << (*cell)[n][i][j].A[0] << endl <<
			"bar_Vx       = " << (*cell)[n][i][j].bar_Vx[0] << endl <<
			"lower bar_Vx = " << (*cell)[n][i][j-1].bar_Vx[0] << endl <<
			"upper bar_Vx = " << (*cell)[n][i][j+1].bar_Vx[0] << endl <<
			"rho       = " << (*cell)[n+1][i][j].rho << endl <<
			"lower rho = " << (*cell)[n+1][i][j-1].rho << endl <<
			"upper rho = " << (*cell)[n+1][i][j+1].rho << endl <<
			"prev rho       = " << (*cell)[n][i][j].rho << endl <<
			"lower prev rho = " << (*cell)[n][i][j-1].rho << endl <<
			"upper prev rho = " << (*cell)[n][i][j+1].rho << endl <<
			"D       = {" << (*cell)[n][i][j].D[1] << ", " << (*cell)[n][i][j].D[2] << ", " << (*cell)[n][i][j].D[3] << ", " << (*cell)[n][i][j].D[4] << "} " << endl <<
			"dM       = {" << (*cell)[n][i][j].dM[1] << ", " << (*cell)[n][i][j].dM[2] << ", " << (*cell)[n][i][j].dM[3] << ", " << (*cell)[n][i][j].dM[4] << "} " << endl <<
			"lower dM = {" << (*cell)[n][i][j-1].dM[1] << ", " << (*cell)[n][i][j-1].dM[2] << ", " << (*cell)[n][i][j-1].dM[3] << ", " << (*cell)[n][i][j-1].dM[4] << "} " << endl <<
			"upper dM = {" << (*cell)[n][i][j+1].dM[1] << ", " << (*cell)[n][i][j+1].dM[2] << ", " << (*cell)[n][i][j+1].dM[3] << ", " << (*cell)[n][i][j+1].dM[4] << "} " << endl << endl;
		getchar();
	}
	
	if (result != result) {
		cout << "final_e is NaN" << endl;
		cout << "i = " << i << endl
			<< "j = " << j << endl;
		cout << setiosflags(ios::fixed) << setprecision(16) << 
			"bar_E       = " << (*cell)[n+1][i][j].bar_e << endl <<
			"lower bar_E = " << (*cell)[n+1][i][j-1].bar_e << endl <<
			"upper bar_E = " << (*cell)[n+1][i][j+1].bar_e << endl <<
			"E       = " << (*cell)[n+1][i][j].e << endl <<
			"lower E = " << (*cell)[n+1][i][j-1].e << endl <<
			"upper E = " << (*cell)[n+1][i][j+1].e << endl <<
			"P[0] = " << (*cell)[n][i][j].P[0] << endl <<
			"A[0] = " << (*cell)[n][i][j].A[0] << endl <<
			"bar_Vx       = " << (*cell)[n+1][i][j].bar_Vx[0] << endl <<
			"lower bar_Vx = " << (*cell)[n+1][i][j-1].bar_Vx[0] << endl <<
			"upper bar_Vx = " << (*cell)[n+1][i][j+1].bar_Vx[0] << endl <<
			"rho       = " << (*cell)[n+1][i][j].rho << endl <<
			"lower rho = " << (*cell)[n+1][i][j-1].rho << endl <<
			"upper rho = " << (*cell)[n+1][i][j+1].rho << endl <<
			"prev rho       = " << (*cell)[n][i][j].rho << endl <<
			"lower prev rho = " << (*cell)[n][i][j-1].rho << endl <<
			"upper prev rho = " << (*cell)[n][i][j+1].rho << endl <<
			"dM       = {" << (*cell)[n+1][i][j].dM[1] << ", " << (*cell)[n+1][i][j].dM[2] << ", " << (*cell)[n+1][i][j].dM[3] << ", " << (*cell)[n+1][i][j].dM[4] << "} " << endl <<
			"lower dM = {" << (*cell)[n+1][i][j-1].dM[1] << ", " << (*cell)[n+1][i][j-1].dM[2] << ", " << (*cell)[n+1][i][j-1].dM[3] << ", " << (*cell)[n+1][i][j-1].dM[4] << "} " << endl <<
			"upper dM = {" << (*cell)[n+1][i][j+1].dM[1] << ", " << (*cell)[n+1][i][j+1].dM[2] << ", " << (*cell)[n+1][i][j+1].dM[3] << ", " << (*cell)[n+1][i][j+1].dM[4] << "} " << endl << endl;
		getchar();
	}
	
    return result;
}
