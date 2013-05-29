#include "border.h"


/* Cell geometry parameter calculation
 *
 * TODO: Fix full[0], [1] and [2] EVERYWHERE
 *
 *  */
void pre_cell_geometry(double array[5], gasCell cell, int i, int j) {
    double full[5];
	full[0] = M_PI*(2*(j-axis_j)+0.5)*pow(dr,2)*dx;
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
	    array[0] = (2*M_PI * (j*dr + (cell.r_1*dr + cell.r_2*dr)/2) * (cell.r_1*dr + cell.r_2*dr)/2 * dx) / full[0] ;
	    array[1] = (M_PI * (2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = (M_PI * (2*j*dr*cell.r_2*dr + pow(cell.r_2*dr, 2))) / full[2] ;
	    array[3] = 1;
	    array[4] = 0;
//	    array[0] = ((j-1)*(cell.r_1+cell.r_2) + cell.r_1*cell.r_2 +
//			pow(cell.r_2-cell.r_1, 2)/3) / (2*j-1);
//		array[1] = cell.r_1 * (j-1+cell.r_1/2) / (j-0.5);
//		array[2] = cell.r_2 * (j-1+cell.r_2/2) / (j-0.5);
//		array[3] = 1;
//		array[4] = 0;
	    break;

	case 2:
	    array[0] = (2*M_PI*(j*dr + dr/2)*dx*dr - 2*M_PI * (j*dr + (cell.r_1*dr + cell.r_2*dr)/2) * (cell.r_1*dr + cell.r_2*dr)/2 * dx) / full[0] ;
	    array[1] = (M_PI * (2*(j+1)*dr*cell.r_1*dr - pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = (M_PI * (2*(j+1)*dr*cell.r_2*dr - pow(cell.r_2*dr, 2))) / full[2] ;
	    array[3] = 0;
	    array[4] = 1;
	    break;

	case 3:
	    array[0] = (2*M_PI * (j*dr + dr/2) * (cell.x_1*dx + cell.x_2*dx)/2 * dr) / full[0] ;
	    array[1] = 1;
	    array[2] = 0;
	    array[3] = (2*M_PI*j*cell.x_1*dx*dr) / full[3] ;
	    array[4] = (2*M_PI*(j+1)*cell.x_2*dx*dr) / full[4] ;
//	    array[0] = (j*cell.x_2 + (j-1)*cell.x_1 - (cell.x_2 - cell.x_1)/3) / (2*j-1);
//	    array[1] = 1;
//	    array[2] = 0;
//	    array[3] = cell.x_1;
//	    array[4] = cell.x_2;
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
	    array[0] = (2*M_PI*(j*dr + cell.r_1*dr/2) * cell.r_1*dr*cell.x_1*dx / 2) / full[0] ;
	    array[1] = (M_PI * (2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
	    array[2] = 0;
	    array[3] = (2*M_PI*j*cell.x_1*dx*dr) / full[3];
	    array[4] = 0;
//	    array[0] = cell.x_1*cell.r_1*(j-1 + cell.r_1/3) / (2*j-1);
//	    array[1] = cell.r_1 * (j-1 + cell.r_1/2) / (j-0.5);
//	    array[2] = 0;
//	    array[3] = cell.x_1;
//	    array[4] = 0;
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
		array[3] = 0;
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
		array[0] = (2*M_PI*(j*dr + cell.r_1*dr/2) * cell.r_1*dr*cell.x_1*dx / 2) / full[0] ;
		array[1] = (M_PI * (2*j*dr*cell.r_1*dr + pow(cell.r_1*dr, 2))) / full[1] ;
		array[2] = 0;
		array[3] = (2*M_PI*j*cell.x_1*dx*dr) / full[3];
		array[4] = 0;
//	    array[0] = ((j-1)*(cell.r_1+cell.r_2) + cell.r_1*cell.r_2 +
//			pow(cell.r_2-cell.r_1, 2)/3) / (2*j-1);
//		array[1] = cell.r_1 * (j-1+cell.r_1/2) / (j-0.5);
//		array[2] = 0;
//		array[3] = 1;
//		array[4] = 0;
		break;


	default:
		break;
	}
}




double polygonArea(double *X, double *Y, int points) {

  double  area=0. ;
  int     i, j=points-1  ;

  for (i=0; i<points; i++) {
    area+=(X[j]+X[i])*(Y[j]-Y[i]); j=i; }

return area*.5; }



void getLineAngle(cell2d cell, int i, int j, int n, Line2D& line,
	LineAngle2D& angle, bool debug) {
	
	line.ybegin = j+1;
	line.yend = j;
	line.xbegin = i;
	line.xend = i+1;
	
	/* TODO: here i assume we have angle from left top to right bottom */
	switch (cell.at(n).at(i).at(j).type) {
		case 1:
		line.xbegin = i*dx;
		line.xend = (i+1)*dx;
		line.ybegin = j*dr + cell.at(n).at(i).at(j).r_1*dr;
		line.yend = j*dr + cell.at(n).at(i).at(j).r_2*dr;
		break;
		
		case 3:
		line.xbegin = i*dx + cell.at(n).at(i).at(j).x_2*dx;
		line.xend = i*dx + cell.at(n).at(i).at(j).x_1*dx;
		line.ybegin = (j+1)*dr;
		line.yend = j*dr;
		break;
		
		case 8:
		line.xbegin = i*dx + cell.at(n).at(i).at(j).x_2*dx;
		line.xend = (i+1)*dx;
		line.ybegin = (j+1)*dr;
		line.yend = j*dr + cell.at(n).at(i).at(j).r_2*dr;
		break;
		
		case 10:
		line.xbegin = i*dx;
		line.xend = i*dx + cell.at(n).at(i).at(j).x_1*dx;
		line.ybegin = j*dr + cell.at(n).at(i).at(j).r_1*dr;
		line.yend = j*dr;
		break;
		
		case 22:
		line.xbegin = i*dx;
		line.xend = (i+1)*dx;
		line.ybegin = j*dr + cell.at(n).at(i).at(j).r_1*dr;
		line.yend = j*dr + cell.at(n).at(i).at(j).r_2*dr;
		break;

		default:
		break;
	}
	angle.cos_a = (line.xend - line.xbegin) / sqrt(pow(line.yend-line.ybegin,2)+pow(line.xend-line.xbegin,2));
	angle.sin_a = (line.yend - line.ybegin) / sqrt(pow(line.yend-line.ybegin,2)+pow(line.xend-line.xbegin,2));
	angle.cos_2a = pow(angle.cos_a,2) - pow(angle.sin_a,2);
	angle.sin_2a = 2*angle.cos_a*angle.sin_a;
	
	if (debug) {
		printf("x_1 = %4.4f, r_1 = %4.4f, x_2 = %4.4f, r_2 = %4.4f\n",
			cell.at(n).at(i).at(j).x_1,cell.at(n).at(i).at(j).r_1,
			cell.at(n).at(i).at(j).x_2,cell.at(n).at(i).at(j).r_2);
		printf("Line begin = %4.4f:%4.4f, line end = %4.4f:%4.4f, line angle = %4.4f\n",
			line.xbegin,line.ybegin,  line.xend,line.yend,  
			asin(angle.sin_a)*180/M_PI);
		
		getchar();
	}
}



void setVertices(int weightCell, Point2D vertices[4], bool debug, int i, int j) {
	int curI = i;
	int curJ = j;
	if (weightCell == 0) {
		curI = i+1; 
		curJ = j;
	} else if (weightCell == 1) {
		curI = i; 
		curJ = j+1;
	} else if (weightCell == 2) {
		curI = i+1;
		curJ = j+1;
	}
	// (curI+0.5)*dx returns us to the center of the cell
	vertices[0].x = (curI)*dx; 
	vertices[1].x = (curI+1)*dx; 
	vertices[2].x = (curI+1)*dx; 
	vertices[3].x = (curI)*dx; 
	vertices[0].y = (curJ)*dr;
	vertices[1].y = (curJ)*dr;
	vertices[2].y = (curJ+1)*dr;
	vertices[3].y = (curJ+1)*dr;
	
	if (debug) {
		printf("Orig vertices at i = %d, j = %d, weightCell = %d\n",
			i,j,weightCell);
		for (unsigned int idx2 = 0; idx2 < 4; idx2++) 
			printf("%d: %4.4f:%4.4f\n", idx2, vertices[idx2].x, vertices[idx2].y);
		
		getchar();
	}
}

void getMirrorVerts(Point2D vertices[4], Int2D vertices_ij[4], 
	Line2D line, LineAngle2D angle, bool debug, int i, int j) {
		
	for (unsigned int idx = 0; idx < 4; idx++) {
		double dist = (
			(line.xend - line.xbegin) * (line.ybegin - vertices[idx].y) - 
			(line.xbegin - vertices[idx].x) * (line.yend - line.ybegin)
		    ) / sqrt(pow(line.xend-line.xbegin,2)+pow(line.yend-line.ybegin,2));
		double diffX = 2 * dist * angle.sin_a;
		double diffY = 2 * dist * angle.cos_a;
		vertices[idx].x -= diffX;
		vertices[idx].y += diffY;
		vertices_ij[idx].i = floor(vertices[idx].x/dx); 
		vertices_ij[idx].j = floor(vertices[idx].y/dr);
		
		if (debug) {
			printf("Point %u, distance total: %4.4f, x: %4.4f, y: %4.4f\n", idx,dist,diffX,diffY);
		}
	}
	
	if (debug) {
		printf("Vertices at i = %d, j = %d\n", i,j);
			for (unsigned int idx2 = 0; idx2 < 4; idx2++) 
				printf("%d: %4.4f:%4.4f\n", idx2, vertices[idx2].x, vertices[idx2].y);
		
		printf("Vertices_ij at i = %d, j = %d\n", i,j);
			for (unsigned int idx2 = 0; idx2 < 4; idx2++) 
				printf("%d: %d:%d\n", idx2, vertices_ij[idx2].i, vertices_ij[idx2].j);
				
		getchar();
	}
}

bool onCross(double x, double y) {
	double deltaX = pow(10,-6)*dx;
	double deltaY = pow(10,-6)*dx;
	double delta = pow(10,-5)*dx;
	return (fmod(x+deltaX,dx) < delta && fmod(y+deltaY,dr) < delta);
}

std::vector <TPoint2D> getExtPoints(Point2D vertices[4],
		Int2D vertices_ij[4], bool debug) {
	/**
	 * Axis intersections and point types
	 * 
	 * Point types:
	 * 0 - vertex
	 * 1 - axis intersection
	 * 2 - internal point 
	 * 
	 **/
	std::vector <TPoint2D> result;
	for (int idx1 = 0; idx1 < 4; idx1++) {
		TPoint2D point;
		point.x = vertices[idx1].x; point.y = vertices[idx1].y;
		point.type = 0;
		result.push_back(point);

		if (debug) {
			std::string cross = onCross(point.x, point.y) ? "true" : "false";
			printf("onCross returned %s for point %6.6f:%6.6f\n",
					cross.c_str(), point.x, point.y);
		}

		if (onCross(point.x,point.y)) {
			continue;
		}
		if (fabs(fmod(point.x, dx)) < pow(10,-6) && fabs(fmod(point.y, dr)) < pow(10,-6)) continue;
		int idx2 = idx1 != 3 ? idx1 + 1 : 0;
		// One intersection point on Y axis
		if (vertices_ij[idx1].i != vertices_ij[idx2].i && vertices_ij[idx1].j == vertices_ij[idx2].j) {
			double interX = floor(fmax(vertices[idx1].x,vertices[idx2].x)/dx)*dx;
			double interY = vertices[idx2].y + (interX - vertices[idx2].x) * (vertices[idx2].y - vertices[idx1].y) / (vertices[idx2].x - vertices[idx1].x);
			point.x = interX; point.y = interY; point.type = 1;
			if (!onCross(point.x,point.y)) {
				result.push_back(point);
			}
		} else
		// One intersection point on X axis
		if (vertices_ij[idx1].i == vertices_ij[idx2].i && vertices_ij[idx1].j != vertices_ij[idx2].j) {
			unsigned int idx_max = vertices[idx1].y > vertices[idx2].y ? idx1 : idx2;
			double interY = floor(vertices[idx_max].y/dr)*dr;
			double interX = vertices[idx2].x + (interY - vertices[idx2].y) * (vertices[idx2].x - vertices[idx1].x) / (vertices[idx2].y - vertices[idx1].y);
			point.x = interX; point.y = interY; point.type = 1;
			if (!onCross(point.x,point.y)) {
				result.push_back(point);
			}
		} else
		// If has 2 intersection result - first serve the closest one (because triangle gen is buggy)
		if (vertices_ij[idx1].i != vertices_ij[idx2].i && vertices_ij[idx1].j != vertices_ij[idx2].j) {
			double interX1 = floor(fmax(vertices[idx1].x,vertices[idx2].x)/dx)*dx;
			double interY1 = vertices[idx2].y + (interX1 - vertices[idx2].x) * (vertices[idx2].y - vertices[idx1].y) / (vertices[idx2].x - vertices[idx1].x);
			double interY2 = floor(fmax(vertices[idx1].y,vertices[idx2].y)/dr)*dr;
			double interX2 = vertices[idx2].x + (interY2 - vertices[idx2].y) * (vertices[idx2].x - vertices[idx1].x) / (vertices[idx2].y - vertices[idx1].y);
			if (pow(vertices[idx1].x-interX1,2) + pow(vertices[idx1].y-interY1,2) < pow(vertices[idx1].x-interX2,2) + pow(vertices[idx1].y-interY2,2)) {
				point.x = interX1; point.y = interY1; point.type = 1;
				if (!onCross(point.x,point.y)) {
					result.push_back(point);
				}
				point.x = interX2; point.y = interY2 ; point.type = 1;
				if (!onCross(point.x,point.y)) {
					result.push_back(point);
				}
			} else {
				point.x = interX2; point.y = interY2; point.type = 1;
				if (!onCross(point.x,point.y)) {
					result.push_back(point);
				}
				point.x = interX1; point.y = interY1; point.type = 1;
				if (!onCross(point.x,point.y)) {
					result.push_back(point);
				}
			}
		}
	}
	
	return result;
}



Int2D getMaxDifference(unsigned int max_i_point[2], unsigned int max_j_point[2],
		Int2D vertices_ij[4], bool debug, int i, int j) {
			
	Int2D max_diff;
	max_diff.i = -1;
	max_diff.j = -1;
	for (unsigned int idx2 = 0; idx2 < 4; idx2++) {
		for (unsigned int idx3 = 0; idx3 < 4; idx3++) {
			if (vertices_ij[idx2].i - vertices_ij[idx3].i > max_diff.i) {
				max_diff.i = vertices_ij[idx2].i - vertices_ij[idx3].i;
				max_i_point[0] = idx2;
				max_i_point[1] = idx3;
			}
			if (vertices_ij[idx2].j - vertices_ij[idx3].j > max_diff.j) {
				max_diff.j = vertices_ij[idx2].j - vertices_ij[idx3].j;
				max_j_point[0] = idx2;
				max_j_point[1] = idx3;
			}
		}
	}
	
	if (debug) {
		printf("Cell %d:%d\n", i,j);
		printf("max_diff.i = %d, max_i_points : %d, %d\n",
			max_diff.i,max_i_point[0],max_i_point[1]);
		printf("max_diff.j = %d, max_j_points : %d, %d\n",
			max_diff.j,max_j_point[0],max_j_point[1]);
	}
	
	return max_diff;
}



std::vector <TPoint2D> clearDoubles(std::vector <TPoint2D> points) {
	//~ bool changed = true;
	//~ while (changed) {
		//~ changed = false;
		//~ for (unsigned int idx = 0; idx < points.size(); idx++) {
			//~ for (unsigned int idx2 = 0; idx2 < points.size(); idx2++) {
				//~ if (fabs(points.at(idx2).x - points.at(idx).x) < pow(10,-6) &&
				    //~ fabs(points.at(idx2).y - points.at(idx).y) < pow(10,-6) &&
				    //~ points.at(idx).type != 2 &&
				    //~ idx2 != idx) {
					//~ points.erase(points.begin()+idx2);
					//~ changed = true;
				//~ }
			//~ }
		//~ }
	//~ }
	//~
	return points;
}

double getDiff(std::vector <TPoint2D> points, unsigned int idx1,
		unsigned int idx2, bool byX) {
	if (byX) {
		return floor(points.at(idx1).x/dx) - floor(points.at(idx2).x/dx);
	} else {
		return floor(points.at(idx1).y/dr) - floor(points.at(idx2).y/dr);
	}
}

/*
 * nvert: Number of vertices in the polygon. Whether to repeat the first vertex at the end.
 * vertx, verty: Arrays containing the x- and y-coordinates of the polygon's vertices.
 * testx, testy: X- and y-coordinate of the test point.
 */

int pnpoly(int nvert, float *vertx, float *verty, float testx, float testy)
{
  int i, j, c = 0;
  for (i = 0, j = nvert-1; i < nvert; j = i++) {
    if ( ((verty[i]>testy) != (verty[j]>testy)) &&
     (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) )
       c = !c;
  }
  return c;
}


std::vector <TPoint2D> getIntPoints(std::vector <TPoint2D> points,
		Point2D vertices[4], Int2D max_diff, unsigned int max_i_point[2],
		unsigned int max_j_point[2], bool debug) {
			
	TPoint2D point;
	
	if (max_diff.i == 1 && max_diff.j == 1) {
		double internalX = floor(fmax(fmax(vertices[0].x,vertices[1].x),fmax(vertices[2].x,vertices[3].x))/dx)*dx;
		double internalY = floor(fmax(fmax(vertices[0].y,vertices[1].y),fmax(vertices[2].y,vertices[3].y))/dr)*dr;
		point.x = internalX; point.y = internalY; point.type = 2;
		points.push_back(point);
	} else {
		int i0 = floor(points.at(max_i_point[0]).x/dx);
		int i1 = floor(points.at(max_i_point[1]).x/dx);
		int j0 = floor(points.at(max_j_point[0]).y/dr);
		int j1 = floor(points.at(max_j_point[1]).y/dr);
		int minI = i0 < i1 ? i0 : i1;
		int maxI = i0 > i1 ? i0 : i1;
		int minJ = j0 < j1 ? j0 : j1;
		int maxJ= j0 > j1 ? j0 : j1;
		float vertx[4]; float verty[4]; int nvert = 4;
		for (int idx = 0; idx < 4; idx++) {
			vertx[idx] = vertices[idx].x;
			verty[idx] = vertices[idx].y;
		}
		for (int i = minI-1; i <= maxI+1; i++) {
			for (int j = minJ-1; j <= maxJ+1; j++) {
				int ans = pnpoly(nvert, vertx, verty, i*dx, j*dr);
				if (ans == 1) {
					if (debug)
						printf("pnpoly returned %d on point %d:%d\n", ans, i, j);
					point.x = i*dx; point.y = j*dr; point.type = 2;
					points.push_back(point);
				}
			}
		}
	}
	
	return points;
}



std::vector <Int2D> makeCells(Int2D vertices_ij[4],
		unsigned int max_i_point[2], unsigned int max_j_point[2], bool debug) {
	
	std::vector <Int2D> cells; 
	
	unsigned int i_left, i_right, j_top, j_bottom = 0;
	// Left and right
	i_left = vertices_ij[max_i_point[0]].i;
	i_right = vertices_ij[max_i_point[1]].i;
	if (i_left > i_right) {
		unsigned int tmp = i_left;
		i_left = i_right;
		i_right = tmp;
	}
	// Top and bottom
	j_bottom = vertices_ij[max_j_point[0]].j;
	j_top = vertices_ij[max_j_point[1]].j;
	if (j_bottom > j_top) {
		unsigned int tmp = j_bottom;
		j_bottom = j_top;
		j_top = tmp;
	}
	if (debug)
		printf("\nWe have the following cells: \n");
	// Make cells
	for (unsigned int idx2 = i_left; idx2 <= i_right; idx2++) {
		for (unsigned int idx3 = j_bottom; idx3 <= j_top; idx3++) {
			Int2D newCell;
			newCell.i = idx2; newCell.j = idx3;
			cells.push_back(newCell);
			if (debug)
				printf("Cell %u: %d:%d\n", (unsigned int)(cells.size()-1), idx2, idx3);
		}
	}
	
	return cells;
}



std::vector <TPoint2D> getPointsInCell(unsigned int idx,
		std::vector <TPoint2D> points, std::vector <Int2D> cells, 
		bool debug) {
	
	std::vector <TPoint2D> pointsInCell;
	
	double delta = pow(10,-10);
	for (unsigned int idx2 = 0; idx2 < points.size(); idx2++) {
		bool rule[4];
		rule[0] = points.at(idx2).x > cells.at(idx).i*dx - delta;
		rule[1] = points.at(idx2).x < (cells.at(idx).i+1)*dx + delta;
		rule[2] = points.at(idx2).y > cells.at(idx).j*dr - delta;
		rule[3] = points.at(idx2).y < (cells.at(idx).j+1)*dr + delta;
		if (rule[0] && rule[1] && rule[2] && rule[3]) {
			pointsInCell.push_back(points.at(idx2));
			if (debug)
				printf("Cell %u: %4.4f:%4.4f, type %d \n",
						idx,points.at(idx2).x,points.at(idx2).y,
						points.at(idx2).type);
		}
	}
	
	return pointsInCell;
}


static bool IsOnSegment(double xi, double yi, double xj, double yj,
                        double xk, double yk) {
  return (xi <= xk || xj <= xk) && (xk <= xi || xk <= xj) &&
         (yi <= yk || yj <= yk) && (yk <= yi || yk <= yj);
}

static char ComputeDirection(double xi, double yi, double xj, double yj,
                             double xk, double yk) {
  double a = (xk - xi) * (yj - yi);
  double b = (xj - xi) * (yk - yi);
  return a < b ? -1 : a > b ? 1 : 0;
}

/** Do line segments (x1, y1)--(x2, y2) and (x3, y3)--(x4, y4) intersect? */
bool DoLineSegmentsIntersect(double x1, double y1, double x2, double y2,
                             double x3, double y3, double x4, double y4) {
  char d1 = ComputeDirection(x3, y3, x4, y4, x1, y1);
  char d2 = ComputeDirection(x3, y3, x4, y4, x2, y2);
  char d3 = ComputeDirection(x1, y1, x2, y2, x3, y3);
  char d4 = ComputeDirection(x1, y1, x2, y2, x4, y4);
  return (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
          ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) ||
         (d1 == 0 && IsOnSegment(x3, y3, x4, y4, x1, y1)) ||
         (d2 == 0 && IsOnSegment(x3, y3, x4, y4, x2, y2)) ||
         (d3 == 0 && IsOnSegment(x1, y1, x2, y2, x3, y3)) ||
         (d4 == 0 && IsOnSegment(x1, y1, x2, y2, x4, y4));
}

std::vector <TPoint2D> fixIntPointID(std::vector <TPoint2D> pointsInCell,
		bool debug) {
	
	if (debug)
		printf("We have an internal point!\n");

	unsigned int intIdx = 100;
	TPoint2D intPoint;
	for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
		if (pointsInCell.at(idx2).type == 2) {
			intIdx = idx2;
			intPoint = pointsInCell.at(idx2);
			break;
		}
	}
	// If none found - return
	if (intIdx == 100)
		return pointsInCell;

	for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
		Point2D l1[2];
		l1[0].x = pointsInCell.at(intIdx).x;
		l1[0].y = pointsInCell.at(intIdx).y;
		l1[1].x = pointsInCell.at(idx2).x;
		l1[1].y = pointsInCell.at(idx2).y;
		bool intersect = false;
		for (unsigned int idx3 = 0; idx3 < pointsInCell.size(); idx3++) {
			for (unsigned int idx4 = 0; idx4 < pointsInCell.size(); idx4++) {
				if (idx3 != idx4 && idx3 != intIdx && idx3 != idx2
						&& idx4 != intIdx && idx4 != idx2) {
					Point2D l2[2];
					l2[0].x = pointsInCell.at(idx3).x;
					l2[0].y = pointsInCell.at(idx3).y;
					l2[1].x = pointsInCell.at(idx4).x;
					l2[1].y = pointsInCell.at(idx4).y;
					intersect = DoLineSegmentsIntersect(l1[0].x,l1[0].y,l1[1].x,l1[1].y,
							l2[0].x,l2[0].y,l2[1].x,l2[1].y);
				}
			}
		}
		if (!intersect && idx2 != 0) {
			if (debug)
				printf("point will be placed between %u and %u\n", idx2, idx2+1);
			std::vector<TPoint2D>::iterator it = pointsInCell.begin();
			pointsInCell.erase(it+intIdx);
			pointsInCell.insert(it+idx2+1, intPoint);
			break;
		}
	}
	
	return pointsInCell;
}


Vector2dVector triangulateCell(std::vector <TPoint2D> points, bool debug) {
	Vector2dVector a;
	Vector2dVector result;

	for (unsigned int idx = 0; idx < points.size(); idx++) {
		a.push_back( Vector2d(points.at(idx).x, points.at(idx).y) );
	}
	if (debug) {
		printf("Size of point vector %u\n", (unsigned int) a.size());
		printf("Area from triangulate = %6.6f\n", Triangulate::Area(a));
	}
	Triangulate::Process(a, result);
	return result;
}

double triangleArea(double dX0, double dY0, double dX1, double dY1, double dX2, double dY2)
{
    double dArea = ((dX1 - dX0)*(dY2 - dY0) - (dX2 - dX0)*(dY1 - dY0))/2.0;
    return (dArea > 0.0) ? dArea : -dArea;
}
 /**
     *  Returns a convex hull given an unordered array of points.
     */
//    public static function convexHull(data:Array):Array
//    {
//        return findHull( order(data) );
//    }
    /**
     *  Orders an array of points counterclockwise.
     */
std::vector <TPoint2D> fixPointOrder(std::vector <TPoint2D> points, bool debug) {
        // first run through all the points and find the upper left
	TPoint2D p = points.front();
	int n = points.size();
	for (int i = 1; i < n; i++) {
		if (points.at(i).y > p.y)
		{
			p = points.at(i);
		}
		else if (points.at(i).y == p.y && points.at(i).x < p.x)
		{
			p = points.at(i);
		}
	}
	if (debug)
		printf("Top left point is at %6.6f:%6.6f\n", p.x, p.y);
	// next find all the cotangents of the angles made by the point P and the
	// other points
	std::vector <CPoint2D> sorted;
	// we need arrays for positive and negative values, because Array.sort
	// will put sort the negatives backwards.
	std::vector <CPoint2D> pos;
	std::vector <CPoint2D> neg;
	// add points back in order
	for (int i = 0; i < n; i++)
	{
		double a = points.at(i).x - p.x - 0.00001*dx;
		double b = points.at(i).y - p.y - 0.00001*dr;
		double cot = b/a;
		if (cot < 0) {
			CPoint2D cPoint;
			cPoint.point = points.at(i);
			cPoint.cot = cot;
			neg.push_back(cPoint);
		} else {
			CPoint2D cPoint;
			cPoint.point = points.at(i);
			cPoint.cot = cot;
			pos.push_back(cPoint);
		}
	}
	// sort the arrays
	std::sort(pos.begin(), pos.end(), lesserCot());
	std::sort(neg.begin(), neg.end(), lesserCot());
	neg.insert(neg.end(), pos.begin(), pos.end());
	sorted = neg;

	std::vector <TPoint2D> ordered;
	ordered.push_back(p);
	for (int i = 0; i < n; i++)
	{
		if (p.x == sorted.at(i).point.x && p.y == sorted.at(i).point.y)
			continue;
		ordered.push_back(sorted.at(i).point);
	}
	return ordered;
}
/**
 *
 */
double direction(TPoint2D p1, TPoint2D p2, TPoint2D p3) {
	// > 0  is right turn
	// == 0 is collinear
	// < 0  is left turn
	// we only want right turns, usually we want right turns, but
	// flash's grid is flipped on y.
	return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}
/**
 *  Given an array of points ordered counterclockwise, findHull will
 *  filter the points and return an array containing the vertices of a
 *  convex polygon that envelopes those points.
 */
std::vector <TPoint2D> findHull(std::vector <TPoint2D> points) {
	int n = points.size();
	std::vector <TPoint2D> hull;
	hull.push_back(points.at(0)); // add the pivot
	hull.push_back(points.at(1)); // makes first vector

	if (n > 3) {
		for (int i = 2; i < n; i++)	{
			while (direction(hull.at(hull.size() - 2), hull.at(hull.size() - 1), points.at(i)) >= 0)
				hull.pop_back();
			hull.push_back(points.at(i));
		}
	}

	return hull;
}


WeightVector wightVectorsCalc(cell2d& cell, int i, int j, int n, bool debug) {
	/** 
	 * Transform points from center of the cell
	 * 
	 * 	[3]	+---+ [2]
	 * 		|	|
	 * 	[0]	+---+ [1]
	 * 
	 * 			|  cos(a)  -sin(a)	|
	 * R(-a) = 	|					|
	 * 			|  sin(a)  cos(a)	|
	 * 
	 * weightCell[0] is i+1 cell (X), [1] is j+1 (Y), [2] is i+1,j+1
	 * 
	 **/
	if (debug)
		printf("\n\n**************************************************\n\n");
	
	Point2D vertices[4];
	Int2D vertices_ij[4];
	Line2D line;
	LineAngle2D angle;
	WeightVector result;
	std::vector <TPoint2D> points;
	unsigned int max_i_point[2] = {0};
	unsigned int max_j_point[2] = {0};
	Int2D max_diff;
	double origArea = dx*dr;
	std::vector <Int2D> cells;
	
	for (unsigned int weightCell = 0; weightCell < 3; weightCell++) {
		// First, we'll set original cell's vertices
		setVertices(weightCell, vertices, debug, i, j);
		
		// Determining line begin-end points and angle
		getLineAngle(cell, i, j, n, line, angle, debug);
		cell.at(n).at(i).at(j).angle = angle;
		
		// Getting mirrored over line vertices array
		getMirrorVerts(vertices, vertices_ij, line, angle, debug, i, j);
		
		// Getting all points of intersection between mirrored cell and grid
		points = getExtPoints(vertices, vertices_ij, debug);
		
		// Getting maximum i and j difference between mirrored vertices
		max_diff = getMaxDifference(max_i_point, max_j_point,
			vertices_ij, debug, i, j);
			
		// Cleaning points from doubles (occures sometimes)
		points = clearDoubles(points);
		
		// Getting internal points (grid's own intersections)
		points = getIntPoints(points, vertices, max_diff, 
			max_i_point, max_j_point, debug);
			
		// Some debug info about our points
		if (debug) {
			for (unsigned int pointNum = 0; pointNum < points.size(); pointNum++) {
				printf("Point %u at i = %d, j = %d: %4.4f:%4.4f, type = %d\n",
					pointNum, i, j, 
					points.at(pointNum).x, points.at(pointNum).y,
					points.at(pointNum).type);
			}
		}
		
		// We'll now make a vector of all cells where our polygon is
		cells = makeCells(vertices_ij, max_i_point, max_j_point, debug); 
		
		// Debug for making a good text ;)
		if (debug)
			printf("\nNow will arrange points in those cells\n");
		
		// Main loop where each cell's weight is calculated
		double totalWeight = 0;
		for (unsigned int idx = 0; idx < cells.size(); idx++) {
			// We'll use center of each of those cells to determine how much points we have in each cell
			std::vector <TPoint2D> pointsInCell = getPointsInCell(idx, 
				points, cells, debug);
			if (pointsInCell.size() < 3)
				continue;
			// Reorder IDs
			pointsInCell = fixPointOrder(pointsInCell, debug);
			if (debug) {
				printf("test points in order: \n");
				for (unsigned int testIdx = 0; testIdx < pointsInCell.size(); testIdx++) {
					printf("%10.10f:%10.10f\n",pointsInCell.at(testIdx).x,
							pointsInCell.at(testIdx).y);
				}
			getchar();
			}
			// Triangulate
			Vector2dVector triangles = triangulateCell(pointsInCell, debug);
			int tcount = triangles.size()/3;
			// Now we'll use the result in polygonArea function
			double area = 0;
			for (int i=0; i<tcount; i++)
			{
			  const Vector2d &p1 = triangles[i*3+0];
			  const Vector2d &p2 = triangles[i*3+1];
			  const Vector2d &p3 = triangles[i*3+2];
			  if (debug)
				  printf("Triangle %d => (%6.6f,%6.6f) (%6.6f,%6.6f) (%6.6f,%6.6f)\n",
						  i+1,p1.GetX(),p1.GetY(),p2.GetX(),p2.GetY(),p3.GetX(),p3.GetY());
			  area += fabs(triangleArea(p1.GetX(), p1.GetY(), p2.GetX(), p2.GetY(),
					  p3.GetX(), p3.GetY()));
			}
			if (debug) printf("Polygon area: %8.8f\n", area);
			// Populating weightPart and pushing it to the cell's weightVector
			WeightPart weightPart;
			weightPart.weight = area / origArea;
			weightPart.i = cells.at(idx).i;
			weightPart.j = cells.at(idx).j;
			if (debug)
				printf("Weight of this cell: %4.4f\n\n\n",weightPart.weight);
			if (weightCell == 0) {
				result.x.push_back(weightPart);
			} else if (weightCell == 1) {
				result.y.push_back(weightPart);
			} else if (weightCell == 2) {
				result.xy.push_back(weightPart);
			}
			totalWeight += weightPart.weight;
		}
//		 Scaling to 1
//		double scale = 1.0/totalWeight;
//		if (weightCell == 0) {
//			for (unsigned int idx = 0; idx < result.x.size(); idx++) {
//				result.x.at(idx).weight *= scale;
//			}
//		} else if (weightCell == 1) {
//			for (unsigned int idx = 0; idx < result.y.size(); idx++) {
//				result.y.at(idx).weight *= scale;
//			}
//		} else if (weightCell == 2) {
//			for (unsigned int idx = 0; idx < result.xy.size(); idx++) {
//				result.xy.at(idx).weight *= scale;
//			}
//		}
	}

	printf("Weights for cell %d:%d\n",i,j);
	printf("For x+1\n");
	for (unsigned int idx = 0; idx < result.x.size(); idx++) {
		printf("i = %d, j = %d, area = %10.10f, weight = %10.10f\n",
				result.x.at(idx).i,result.x.at(idx).j,result.x.at(idx).weight*origArea,
				result.x.at(idx).weight);
	}
	printf("For y+1\n");
	for (unsigned int idx = 0; idx < result.y.size(); idx++) {
		printf("i = %d, j = %d, area = %10.10f, weight = %10.10f\n",
			result.y.at(idx).i,result.y.at(idx).j,result.y.at(idx).weight*origArea,
			result.y.at(idx).weight);
	}

	if (debug) printf("Total weight parts at %d in cell %d:%d\nby x: %u\nby y: %u\nby xy: %u\n", n,i,j,
		(unsigned int) result.x.size(),
		(unsigned int) result.y.size(),
		(unsigned int) result.xy.size());

	return result;
}



void calculateBorder(int n, cell2dStatic& cell, unsigned long ctrl,
		BorderCond result[10], int i, int j) {

	/* Dummy cells calc */
	int type = cell[i][j].type;
	if (i == i_sn - 1) {
		if (type == 13) {
			type = 20;
		} else if (type == 14) {
			type = 21;
		} else if (type == 0) {
			type = 19;
		}
	}
	if (i == i_sn) {
		if (type == 13) {
			type = 20;
		} else if (type == 14) {
			type = 21;
		} else if (type == 0) {
			type = 19;
		}
	}

	/**
	 * ********Naming rule*********
	 *
	 *			i-1		i		i+1
	 * 	j+1		6	-	7	-	8
	 * 	j		3	-	4	-	5
	 * 	j-1		0	-	1	-	2
	 *
	 * 	*******Control array*******
	 *	Shift	-	parameter
	 *	0		-	P
	 *	1		-	Vx
	 *	2		-	Vr
	 *	3		-	E
	 *	4		-	rho
	 *	5		-	bar_Vx
	 *	6		-	bar_Vr
	 *	7		-	bar_e
	 *
	 */

	bool isSet_P 		= 	CHECK_BIT(ctrl, 0);
	bool isSet_Vx 		= 	CHECK_BIT(ctrl, 1);
	bool isSet_Vr 		= 	CHECK_BIT(ctrl, 2);
	bool isSet_E 		= 	CHECK_BIT(ctrl, 3);
	bool isSet_rho 		= 	CHECK_BIT(ctrl, 4);
	bool isSet_barVx 	= 	CHECK_BIT(ctrl, 5);
	bool isSet_barVr 	= 	CHECK_BIT(ctrl, 6);
	bool isSet_barE 	= 	CHECK_BIT(ctrl, 7);
	bool isSet_z 		= 	CHECK_BIT(ctrl, 8);
	bool isSet_psi 		= 	CHECK_BIT(ctrl, 9);

	bool no_I_1			=	cell[i-1][j+1].type == 18 ? true : false;
	bool no_I1			=	cell[i+1][j+1].type == 18 ? true : false;

	switch (type) {
		// Free cell
		case 0:
			if (isSet_P) {
				result[0].i_1j1 = no_I_1 ? (cell[i][j+1].P[0]+cell[i-1][j].P[0])/2 : cell[i-1][j+1].P[0];
				result[0].ij1 	= cell[i][j+1].P[0];
				result[0].i1j1 	= no_I1 ? (cell[i][j+1].P[0]+cell[i+1][j].P[0])/2 : cell[i+1][j+1].P[0];
				result[0].i_1j 	= cell[i-1][j].P[0];
				result[0].ij 	= cell[i][j].P[0];
				result[0].i1j 	= cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0];
				result[0].ij_1 	= cell[i][j-1].P[0];
				result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = no_I_1 ? (cell[i][j+1].Vx[0]+cell[i-1][j].Vx[0])/2 : cell[i-1][j+1].Vx[0];
				result[1].ij1 = cell[i][j+1].Vx[0];
				result[1].i1j1 = no_I1 ? (cell[i][j+1].Vx[0]+cell[i+1][j].Vx[0])/2 : cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0];
				result[1].ij = cell[i][j].Vx[0];
				result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];
				result[1].ij_1 = cell[i][j-1].Vx[0];
				result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = no_I_1 ? (cell[i][j+1].Vr[0]+cell[i-1][j].Vr[0])/2 : cell[i-1][j+1].Vr[0];
				result[2].ij1 = cell[i][j+1].Vr[0];
				result[2].i1j1 = no_I1 ? (cell[i][j+1].Vr[0]+cell[i+1][j].Vr[0])/2 : cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0];
				result[2].ij = cell[i][j].Vr[0];
				result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];
				result[2].ij_1 = cell[i][j-1].Vr[0];
				result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = no_I_1 ? (cell[i][j+1].e+cell[i-1][j].e)/2 : cell[i-1][j+1].e;
				result[3].ij1 = cell[i][j+1].e;
				result[3].i1j1 = no_I1 ? (cell[i][j+1].e+cell[i+1][j].e)/2 : cell[i+1][j+1].e;
				result[3].i_1j = cell[i-1][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = no_I_1 ? (cell[i][j+1].rho+cell[i-1][j].rho)/2 : cell[i-1][j+1].rho;
				result[4].ij1 = cell[i][j+1].rho;
				result[4].i1j1 = no_I1 ? (cell[i][j+1].rho+cell[i+1][j].rho)/2 : cell[i+1][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = no_I_1 ? (cell[i][j+1].bar_Vx[0]+cell[i-1][j].bar_Vx[0])/2 : cell[i-1][j+1].bar_Vx[0];
				result[5].ij1 = cell[i][j+1].bar_Vx[0];
				result[5].i1j1 = no_I1 ? (cell[i][j+1].bar_Vx[0]+cell[i+1][j].bar_Vx[0])/2 : cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];
				result[5].ij = cell[i][j].bar_Vx[0];
				result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 = cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = no_I_1 ? (cell[i][j+1].bar_Vr[0]+cell[i-1][j].bar_Vr[0])/2 : cell[i-1][j+1].bar_Vr[0];
				result[6].ij1 = cell[i][j+1].bar_Vr[0];
				result[6].i1j1 = no_I1 ? (cell[i][j+1].bar_Vr[0]+cell[i+1][j].bar_Vr[0])/2 : cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];
				result[6].ij = cell[i][j].bar_Vr[0];
				result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 = cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = no_I_1 ? (cell[i][j+1].bar_e+cell[i-1][j].bar_e)/2 : cell[i-1][j+1].bar_e;
				result[7].ij1 = cell[i][j+1].bar_e;
				result[7].i1j1 = no_I1 ? (cell[i][j+1].bar_e+cell[i+1][j].bar_e)/2 : cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;
				result[7].ij = cell[i][j].bar_e;
				result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;
				result[7].ij_1 = cell[i][j-1].bar_e;
				result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = no_I_1 ? (cell[i][j+1].bar_z+cell[i-1][j].bar_z)/2 : cell[i-1][j+1].bar_z;
				result[8].ij1 = cell[i][j+1].bar_z;
				result[8].i1j1 = no_I1 ? (cell[i][j+1].bar_z+cell[i+1][j].bar_z)/2 : cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;
				result[8].ij = cell[i][j].bar_z;
				result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;
				result[8].ij_1 = cell[i][j-1].bar_z;
				result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = no_I_1 ? (cell[i][j+1].bar_psi+cell[i-1][j].bar_psi)/2 : cell[i-1][j+1].bar_psi;
				result[9].ij1 = cell[i][j+1].bar_psi;
				result[9].i1j1 = no_I1 ? (cell[i][j+1].bar_psi+cell[i+1][j].bar_psi)/2 : cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;
				result[9].ij = cell[i][j].bar_psi;
				result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;
				result[9].ij_1 = cell[i][j-1].bar_psi;
				result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}
			break;

		case 1:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j].P[0];
				result[0].ij1 = cell[i][j].P[0];
				result[0].i1j1 = cell[i+1][j].P[0];
				result[0].i_1j = cell[i-1][j].P[0];
				result[0].ij = cell[i][j].P[0];
				result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0];
				result[0].ij_1 = cell[i][j-1].P[0];
				result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j].Vx[0];
				result[1].ij1 = cell[i][j].Vx[0];
				result[1].i1j1 = cell[i+1][j].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0];
				result[1].ij = cell[i][j].Vx[0];
				result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];
				result[1].ij_1 = cell[i][j-1].Vx[0];
				result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = -cell[i-1][j].Vr[0];
				result[2].ij1 = -cell[i][j].Vr[0];
				result[2].i1j1 = -cell[i+1][j].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0];
				result[2].ij = cell[i][j].Vr[0];
				result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];
				result[2].ij_1 = cell[i][j-1].Vr[0];
				result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j].e;
				result[3].ij1 = cell[i][j].e;
				result[3].i1j1 = cell[i+1][j].e;
				result[3].i_1j = cell[i-1][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j].rho;
				result[4].ij1 = cell[i][j].rho;
				result[4].i1j1 = cell[i+1][j].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j].bar_Vx[0];
				result[5].ij1 = cell[i][j].bar_Vx[0];
				result[5].i1j1 = cell[i+1][j].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];
				result[5].ij = cell[i][j].bar_Vx[0];
				result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 = cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = -cell[i-1][j].bar_Vr[0];
				result[6].ij1 = -cell[i][j].bar_Vr[0];
				result[6].i1j1 = -cell[i+1][j].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];
				result[6].ij = cell[i][j].bar_Vr[0];
				result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 = cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j].bar_e;
				result[7].ij1 = cell[i][j].bar_e;
				result[7].i1j1 = cell[i+1][j].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;
				result[7].ij = cell[i][j].bar_e;
				result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;
				result[7].ij_1 = cell[i][j-1].bar_e;
				result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j].bar_z;
				result[8].ij1 = cell[i][j].bar_z;
				result[8].i1j1 = cell[i+1][j].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;
				result[8].ij = cell[i][j].bar_z;
				result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;
				result[8].ij_1 = cell[i][j-1].bar_z;
				result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j].bar_psi;
				result[9].ij1 = cell[i][j].bar_psi;
				result[9].i1j1 = cell[i+1][j].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;
				result[9].ij = cell[i][j].bar_psi;
				result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;
				result[9].ij_1 = cell[i][j-1].bar_psi;
				result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].ij1 = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.y.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.y.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.y.at(idx).i;
				weightCell.j = cell[i][j].weightVector.y.at(idx).j;
				if (isSet_P) {
					result[0].ij1 += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].ij1 += weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].ij1 -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].ij1 += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].ij1 += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].ij1 += weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].ij1 -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].ij1 += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].ij1 += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].ij1 += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].i1j1 = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.xy.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.xy.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.xy.at(idx).i;
				weightCell.j = cell[i][j].weightVector.xy.at(idx).j;
				if (isSet_P) {
					result[0].i1j1 += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].i1j1 += weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].i1j1 -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].i1j1 += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].i1j1 += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].i1j1 -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}
			break;

		case 3:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j+1].P[0];
				result[0].ij1 		= cell[i][j+1].P[0];
				result[0].i1j1 		= cell[i+1][j+1].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i+1][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= -cell[i][j+1].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= -cell[i][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= -cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0];
				result[2].ij1 = cell[i][j+1].Vr[0];
				result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0];
				result[2].ij = cell[i][j].Vr[0];
				result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];
				result[2].ij_1 = cell[i][j-1].Vr[0];
				result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= -cell[i][j+1].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= -cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}
			break;

		case 8:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j+1].P[0]; 	result[0].ij1 = cell[i][j+1].P[0];	result[0].i1j1 = cell[i+1][j+1].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];	result[1].i1j1 = cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];	result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];	result[5].ij_1 = cell[i][j-1].bar_Vx[0];	result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}
			break;


		case 10:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j].P[0];
				result[0].ij1 		= cell[i][j].P[0];
				result[0].i1j1 		= cell[i][j].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j].Vx[0];
				result[1].ij1 		= cell[i][j].Vx[0];
				result[1].i1j1 		= -cell[i][j].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= -cell[i][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= -cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= -cell[i-1][j].Vr[0];
				result[2].ij1 		= -cell[i][j].Vr[0];
				result[2].i1j1 		= -cell[i][j].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i][j].Vr[0];
				result[2].i_1j_1 	= cell[i-1][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j].e;
				result[3].ij1 = cell[i][j].e;
				result[3].i1j1 = cell[i][j].e;
				result[3].i_1j = cell[i-1][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j].rho;
				result[4].ij1 = cell[i][j].rho;
				result[4].i1j1 = cell[i][j].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j].bar_Vx[0];
				result[5].ij1 		= cell[i][j].bar_Vx[0];
				result[5].i1j1 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= -cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= -cell[i-1][j].bar_Vr[0];
				result[6].ij1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j1 		= -cell[i][j].bar_Vr[0];
				result[6].i_1j 		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i][j].bar_Vr[0];
				result[6].i_1j_1 	= cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 		= cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 	= cell[i][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j].bar_e;	result[7].ij1 = cell[i][j].bar_e;	result[7].i1j1 = cell[i][j].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j].bar_z;	result[8].ij1 = cell[i][j].bar_z;	result[8].i1j1 = cell[i][j].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j].bar_psi;	result[9].ij1 = cell[i][j].bar_psi;	result[9].i1j1 = cell[i][j].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i][j-1].bar_psi;
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].i1j = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.x.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.x.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.x.at(idx).i;
				weightCell.j = cell[i][j].weightVector.x.at(idx).j;
				if (isSet_P) {
					result[0].i1j += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].i1j -= weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].i1j -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].i1j += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].i1j += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].i1j -= weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].i1j -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].i1j += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].i1j += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].i1j += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].ij1 = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.y.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.y.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.y.at(idx).i;
				weightCell.j = cell[i][j].weightVector.y.at(idx).j;
				if (isSet_P) {
					result[0].ij1 += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].ij1 -= weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].ij1 -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].ij1 += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].ij1 += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].ij1 -= weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].ij1 -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].ij1 += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].ij1 += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].ij1 += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].i1j1 = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.xy.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.xy.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.xy.at(idx).i;
				weightCell.j = cell[i][j].weightVector.xy.at(idx).j;
				if (isSet_P) {
					result[0].i1j1 += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].i1j1 -= weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].i1j1 -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].i1j1 += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].i1j1 += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].i1j1 -= weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].i1j1 -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}
			break;

		// Top border closed
		case 13:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j].P[0];
				result[0].ij1 		= cell[i][j].P[0];
				result[0].i1j1 		= cell[i+1][j].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i+1][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j].Vx[0];
				result[1].ij1 		= cell[i][j].Vx[0];
				result[1].i1j1 		= cell[i+1][j].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= cell[i+1][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= -cell[i-1][j].Vr[0];
				result[2].ij1 		= -cell[i][j].Vr[0];
				result[2].i1j1 		= -cell[i+1][j].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i+1][j].Vr[0];
				result[2].i_1j_1 	= cell[i-1][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 	= cell[i-1][j].e;
				result[3].ij1 		= cell[i][j].e;
				result[3].i1j1 		= cell[i+1][j].e;
				result[3].i_1j 		= cell[i-1][j].e;
				result[3].ij 		= cell[i][j].e;
				result[3].i1j 		= cell[i+1][j].e;
				result[3].i_1j_1 	= cell[i-1][j-1].e;
				result[3].ij_1 		= cell[i][j-1].e;
				result[3].i1j_1 	= cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j].rho;
				result[4].ij1 = cell[i][j].rho;
				result[4].i1j1 = cell[i+1][j].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j].bar_Vx[0];
				result[5].ij1 		= cell[i][j].bar_Vx[0];
				result[5].i1j1 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= -cell[i-1][j].bar_Vr[0];
				result[6].ij1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j1 		= -cell[i+1][j].bar_Vr[0];
				result[6].i_1j 		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 	= cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 		= cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 	= cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j].bar_e;
				result[7].ij1 = cell[i][j].bar_e;
				result[7].i1j1 = cell[i+1][j].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;
				result[7].ij = cell[i][j].bar_e;
				result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;
				result[7].ij_1 = cell[i][j-1].bar_e;
				result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = no_I_1 ? cell[i-1][j].bar_z : cell[i-1][j+1].bar_z;
				result[8].ij1 = cell[i][j].bar_z;
				result[8].i1j1 = no_I1 ? cell[i+1][j].bar_z : cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;
				result[8].ij = cell[i][j].bar_z;
				result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;
				result[8].ij_1 = cell[i][j-1].bar_z;
				result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = no_I_1 ? cell[i-1][j].bar_psi : cell[i-1][j+1].bar_psi;
				result[9].ij1 = cell[i][j].bar_psi;
				result[9].i1j1 = no_I1 ? cell[i+1][j].bar_psi : cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;
				result[9].ij = cell[i][j].bar_psi;
				result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;
				result[9].ij_1 = cell[i][j-1].bar_psi;
				result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}
			break;

		// Bottom border closed
		case 14:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j+1].P[0];
				result[0].ij1 		= cell[i][j+1].P[0];
				result[0].i1j1 		= cell[i+1][j+1].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i+1][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j].P[0];
				result[0].ij_1 		= cell[i][j].P[0];
				result[0].i1j_1 	= cell[i+1][j].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= cell[i+1][j+1].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= cell[i+1][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j].Vx[0];
				result[1].ij_1 		= cell[i][j].Vx[0];
				result[1].i1j_1 	= cell[i+1][j].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= cell[i-1][j+1].Vr[0];
				result[2].ij1 		= cell[i][j+1].Vr[0];
				result[2].i1j1 		= cell[i+1][j+1].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i+1][j].Vr[0];
				result[2].i_1j_1 	= -cell[i-1][j].Vr[0];
				result[2].ij_1 		= -cell[i][j].Vr[0];
				result[2].i1j_1 	= -cell[i+1][j].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 	= cell[i-1][j+1].e;
				result[3].ij1 		= cell[i][j+1].e;
				result[3].i1j1 		= cell[i+1][j+1].e;
				result[3].i_1j 		= cell[i-1][j].e;
				result[3].ij 		= cell[i][j].e;
				result[3].i1j 		= cell[i+1][j].e;
				result[3].i_1j_1 	= cell[i-1][j].e;
				result[3].ij_1 		= cell[i][j].e;
				result[3].i1j_1 	= cell[i+1][j].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;
				result[4].ij1 = cell[i][j+1].rho;
				result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j].rho;
				result[4].ij_1 = cell[i][j].rho;
				result[4].i1j_1 = cell[i+1][j].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j].bar_Vx[0];
				result[5].ij_1 		= cell[i][j].bar_Vx[0];
				result[5].i1j_1 	= cell[i+1][j].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= cell[i-1][j+1].bar_Vr[0];
				result[6].ij1 		= cell[i][j+1].bar_Vr[0];
				result[6].i1j1 		= cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j 		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 	= -cell[i-1][j].bar_Vr[0];
				result[6].ij_1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j_1 	= -cell[i+1][j].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;
				result[7].ij1 = cell[i][j+1].bar_e;
				result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;
				result[7].ij = cell[i][j].bar_e;
				result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j].bar_e;
				result[7].ij_1 = cell[i][j].bar_e;
				result[7].i1j_1 = cell[i+1][j].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j].bar_z;	result[8].ij_1 = cell[i][j].bar_z;	result[8].i1j_1 = cell[i+1][j].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j].bar_psi;	result[9].ij_1 = cell[i][j].bar_psi;	result[9].i1j_1 = cell[i+1][j].bar_psi;
			}
			break;

		// Top and left borders closed
		case 15:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i][j].P[0];
				result[0].ij1 		= cell[i][j].P[0];
				result[0].i1j1 		= cell[i+1][j].P[0];
				result[0].i_1j 		= cell[i][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i+1][j].P[0];
				result[0].i_1j_1 	= cell[i][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= -cell[i][j].Vx[0];
				result[1].ij1 		= cell[i][j].Vx[0];
				result[1].i1j1 		= cell[i+1][j].Vx[0];
				result[1].i_1j 		= -cell[i][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= cell[i+1][j].Vx[0];
				result[1].i_1j_1 	= -cell[i][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= -cell[i][j].Vr[0];
				result[2].ij1 		= -cell[i][j].Vr[0];
				result[2].i1j1 		= -cell[i+1][j].Vr[0];
				result[2].i_1j 		= cell[i][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i+1][j].Vr[0];
				result[2].i_1j_1 	= cell[i][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i][j].e;
				result[3].ij1 = cell[i][j].e;
				result[3].i1j1 = cell[i+1][j].e;
				result[3].i_1j = cell[i][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i][j].rho;
				result[4].ij1 = cell[i][j].rho;
				result[4].i1j1 = cell[i+1][j].rho;
				result[4].i_1j = cell[i][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= -cell[i][j].bar_Vx[0];
				result[5].ij1 		= cell[i][j].bar_Vx[0];
				result[5].i1j1 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j 		= -cell[i][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 	= -cell[i][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= -cell[i][j].bar_Vr[0];
				result[6].ij1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j1 		= -cell[i+1][j].bar_Vr[0];
				result[6].i_1j 		= cell[i][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 	= cell[i][j-1].bar_Vr[0];
				result[6].ij_1 		= cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 	= cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i][j].bar_e;	result[7].ij1 = cell[i][j].bar_e;	result[7].i1j1 = cell[i+1][j].bar_e;
				result[7].i_1j = cell[i][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i][j].bar_z;	result[8].ij1 = cell[i][j].bar_z;	result[8].i1j1 = cell[i+1][j].bar_z;
				result[8].i_1j = cell[i][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i][j].bar_psi;	result[9].ij1 = cell[i][j].bar_psi;	result[9].i1j1 = cell[i+1][j].bar_psi;
				result[9].i_1j = cell[i][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}
			break;

		// Bottom and left borders closed
		case 16:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i][j+1].P[0];
				result[0].ij1 		= cell[i][j+1].P[0];
				result[0].i1j1 		= cell[i+1][j+1].P[0];
				result[0].i_1j 		= cell[i][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i+1][j].P[0];
				result[0].i_1j_1 	= cell[i][j].P[0];
				result[0].ij_1 		= cell[i][j].P[0];
				result[0].i1j_1 	= cell[i+1][j].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= -cell[i][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= cell[i+1][j+1].Vx[0];
				result[1].i_1j 		= -cell[i][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= cell[i+1][j].Vx[0];
				result[1].i_1j_1 	= -cell[i][j].Vx[0];
				result[1].ij_1 		= cell[i][j].Vx[0];
				result[1].i1j_1 	= cell[i+1][j].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= cell[i][j+1].Vr[0];
				result[2].ij1 		= cell[i][j+1].Vr[0];
				result[2].i1j1 		= cell[i+1][j+1].Vr[0];
				result[2].i_1j 		= cell[i][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i+1][j].Vr[0];
				result[2].i_1j_1 	= -cell[i][j].Vr[0];
				result[2].ij_1 		= -cell[i][j].Vr[0];
				result[2].i1j_1 	= -cell[i+1][j].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i][j+1].e;
				result[3].ij1 = cell[i][j+1].e;
				result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i][j].e;
				result[3].ij_1 = cell[i][j].e;
				result[3].i1j_1 = cell[i+1][j].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i][j+1].rho;
				result[4].ij1 = cell[i][j+1].rho;
				result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i][j].rho;
				result[4].ij_1 = cell[i][j].rho;
				result[4].i1j_1 = cell[i+1][j].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= -cell[i][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j 		= -cell[i][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 	= -cell[i][j].bar_Vx[0];
				result[5].ij_1 		= cell[i][j].bar_Vx[0];
				result[5].i1j_1 	= cell[i+1][j].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= cell[i][j+1].bar_Vr[0];
				result[6].ij1 		= cell[i][j+1].bar_Vr[0];
				result[6].i1j1 		= cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j 		= cell[i][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 	= -cell[i][j].bar_Vr[0];
				result[6].ij_1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j_1 	= -cell[i+1][j].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i][j].bar_e;	result[7].ij_1 = cell[i][j].bar_e;	result[7].i1j_1 = cell[i+1][j].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i][j].bar_z;	result[8].ij_1 = cell[i][j].bar_z;	result[8].i1j_1 = cell[i+1][j].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i][j].bar_psi;	result[9].ij_1 = cell[i][j].bar_psi;	result[9].i1j_1 = cell[i+1][j].bar_psi;
			}
			break;

		// Left border closed
		case 17:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i][j+1].P[0];
				result[0].ij1 		= cell[i][j+1].P[0];
				result[0].i1j1 		= cell[i+1][j+1].P[0];
				result[0].i_1j 		= cell[i][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i+1][j].P[0];
				result[0].i_1j_1 	= cell[i][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= -cell[i][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= cell[i+1][j+1].Vx[0];
				result[1].i_1j 		= -cell[i][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= cell[i+1][j].Vx[0];
				result[1].i_1j_1 	= -cell[i][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= cell[i][j+1].Vr[0];
				result[2].ij1 		= cell[i][j+1].Vr[0];
				result[2].i1j1 		= cell[i+1][j+1].Vr[0];
				result[2].i_1j 		= cell[i][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i+1][j].Vr[0];
				result[2].i_1j_1 	= cell[i][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i][j+1].e;
				result[3].ij1 = cell[i][j+1].e;
				result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i][j+1].rho;
				result[4].ij1 = cell[i][j+1].rho;
				result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= -cell[i][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j 		= -cell[i][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 	= -cell[i][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i][j+1].bar_Vr[0];
				result[6].ij1 = cell[i][j+1].bar_Vr[0];
				result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i][j].bar_Vr[0];
				result[6].ij = cell[i][j].bar_Vr[0];
				result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i][j-1].bar_Vr[0];
				result[6].ij_1 = cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i][j+1].bar_e;
				result[7].ij1 = cell[i][j+1].bar_e;
				result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i][j].bar_e;
				result[7].ij = cell[i][j].bar_e;
				result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i][j-1].bar_e;
				result[7].ij_1 = cell[i][j-1].bar_e;
				result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i+1][j+1].bar_z;
				result[8].i_1j = cell[i][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i+1][j].bar_z;
				result[8].i_1j_1 = cell[i][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i+1][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i+1][j+1].bar_psi;
				result[9].i_1j = cell[i][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i+1][j].bar_psi;
				result[9].i_1j_1 = cell[i][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i+1][j-1].bar_psi;
			}
			break;

		// Right border closed
		case 19:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j+1].P[0];
				result[0].ij1 = cell[i][j+1].P[0];
				result[0].i1j1 = cell[i][j+1].P[0];
				result[0].i_1j = cell[i-1][j].P[0];
				result[0].ij = cell[i][j].P[0];
				result[0].i1j = cell[i][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0];
				result[0].ij_1 = cell[i][j-1].P[0];
				result[0].i1j_1 = cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= -cell[i][j+1].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= -cell[i][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= -cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= cell[i-1][j+1].Vr[0];
				result[2].ij1 		= cell[i][j+1].Vr[0];
				result[2].i1j1 		= cell[i][j+1].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i][j].Vr[0];
				result[2].i_1j_1 	= cell[i-1][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e;
				result[3].ij1 = cell[i][j+1].e;
				result[3].i1j1 = cell[i][j+1].e;
				result[3].i_1j = cell[i-1][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;
				result[4].ij1 = cell[i][j+1].rho;
				result[4].i1j1 = cell[i][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= -cell[i][j+1].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= -cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= cell[i-1][j+1].bar_Vr[0];
				result[6].ij1 		= cell[i][j+1].bar_Vr[0];
				result[6].i1j1 		= cell[i][j+1].bar_Vr[0];
				result[6].i_1j 		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i][j].bar_Vr[0];
				result[6].i_1j_1 	= cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 		= cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 	= cell[i][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i][j-1].bar_psi;
			}
			break;

		// Top and right borders closed
		case 20:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j].P[0];
				result[0].ij1 		= cell[i][j].P[0];
				result[0].i1j1 		= cell[i][j].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j].Vx[0];
				result[1].ij1 		= cell[i][j].Vx[0];
				result[1].i1j1 		= -cell[i][j].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= -cell[i][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= -cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= -cell[i-1][j].Vr[0];
				result[2].ij1 		= -cell[i][j].Vr[0];
				result[2].i1j1 		= -cell[i][j].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i][j].Vr[0];
				result[2].i_1j_1 	= cell[i-1][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j].e;
				result[3].ij1 = cell[i][j].e;
				result[3].i1j1 = cell[i][j].e;
				result[3].i_1j = cell[i-1][j].e;
				result[3].ij = cell[i][j].e;
				result[3].i1j = cell[i][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;
				result[3].ij_1 = cell[i][j-1].e;
				result[3].i1j_1 = cell[i][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j].rho;
				result[4].ij1 = cell[i][j].rho;
				result[4].i1j1 = cell[i][j].rho;
				result[4].i_1j = cell[i-1][j].rho;
				result[4].ij = cell[i][j].rho;
				result[4].i1j = cell[i][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;
				result[4].ij_1 = cell[i][j-1].rho;
				result[4].i1j_1 = cell[i][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j].bar_Vx[0];
				result[5].ij1 		= cell[i][j].bar_Vx[0];
				result[5].i1j1 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= -cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= -cell[i-1][j].bar_Vr[0];
				result[6].ij1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j1 		= -cell[i][j].bar_Vr[0];
				result[6].i_1j 		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i][j].bar_Vr[0];
				result[6].i_1j_1 	= cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 		= cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 	= cell[i][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j].bar_e;
				result[7].ij1 = cell[i][j].bar_e;
				result[7].i1j1 = cell[i][j].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;
				result[7].ij = cell[i][j].bar_e;
				result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;
				result[7].ij_1 = cell[i][j-1].bar_e;
				result[7].i1j_1 = cell[i][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j].bar_z;	result[8].ij1 = cell[i][j].bar_z;	result[8].i1j1 = cell[i][j].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j].bar_psi;	result[9].ij1 = cell[i][j].bar_psi;	result[9].i1j1 = cell[i][j].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i][j-1].bar_psi;
			}
			break;

		// Bottom and right borders closed
		case 21:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j+1].P[0];
				result[0].ij1 		= cell[i][j+1].P[0];
				result[0].i1j1 		= cell[i][j+1].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j].P[0];
				result[0].ij_1 		= cell[i][j].P[0];
				result[0].i1j_1 	= cell[i][j].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= -cell[i][j+1].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= -cell[i][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j].Vx[0];
				result[1].ij_1 		= cell[i][j].Vx[0];
				result[1].i1j_1 	= -cell[i][j].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= cell[i-1][j+1].Vr[0];
				result[2].ij1 		= cell[i][j+1].Vr[0];
				result[2].i1j1 		= cell[i][j+1].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i][j].Vr[0];
				result[2].i_1j_1 	= -cell[i-1][j].Vr[0];
				result[2].ij_1 		= -cell[i][j].Vr[0];
				result[2].i1j_1 	= -cell[i][j].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 	= cell[i-1][j+1].e;
				result[3].ij1 		= cell[i][j+1].e;
				result[3].i1j1 		= cell[i][j+1].e;
				result[3].i_1j 		= cell[i-1][j].e;
				result[3].ij 		= cell[i][j].e;
				result[3].i1j 		= cell[i][j].e;
				result[3].i_1j_1 	= cell[i-1][j].e;
				result[3].ij_1 		= cell[i][j].e;
				result[3].i1j_1 	= cell[i][j].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 	= cell[i-1][j+1].rho;
				result[4].ij1 		= cell[i][j+1].rho;
				result[4].i1j1 		= cell[i][j+1].rho;
				result[4].i_1j 		= cell[i-1][j].rho;
				result[4].ij 		= cell[i][j].rho;
				result[4].i1j 		= cell[i][j].rho;
				result[4].i_1j_1 	= cell[i-1][j].rho;
				result[4].ij_1 		= cell[i][j].rho;
				result[4].i1j_1 	= cell[i][j].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= -cell[i][j+1].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j].bar_Vx[0];
				result[5].ij_1 		= cell[i][j].bar_Vx[0];
				result[5].i1j_1 	= -cell[i][j].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= cell[i-1][j+1].bar_Vr[0];
				result[6].ij1 		= cell[i][j+1].bar_Vr[0];
				result[6].i1j1		= cell[i][j+1].bar_Vr[0];
				result[6].i_1j		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i][j].bar_Vr[0];
				result[6].i_1j_1 	= -cell[i-1][j].bar_Vr[0];
				result[6].ij_1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j_1 	= -cell[i][j].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 	= cell[i-1][j+1].bar_e;
				result[7].ij1 		= cell[i][j+1].bar_e;
				result[7].i1j1 		= cell[i][j+1].bar_e;
				result[7].i_1j 		= cell[i-1][j].bar_e;
				result[7].ij 		= cell[i][j].bar_e;
				result[7].i1j 		= cell[i][j].bar_e;
				result[7].i_1j_1 	= cell[i-1][j].bar_e;
				result[7].ij_1 		= cell[i][j].bar_e;
				result[7].i1j_1 	= cell[i][j].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j].bar_z;	result[8].ij_1 = cell[i][j].bar_z;	result[8].i1j_1 = cell[i][j].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j].bar_psi;	result[9].ij_1 = cell[i][j].bar_psi;	result[9].i1j_1 = cell[i][j].bar_psi;
			}
			break;

		case 22:
			if (isSet_P) {
				result[0].i_1j1 	= cell[i-1][j+1].P[0];
				result[0].ij1 		= cell[i][j+1].P[0];
				result[0].i1j1 		= cell[i][j].P[0];
				result[0].i_1j 		= cell[i-1][j].P[0];
				result[0].ij 		= cell[i][j].P[0];
				result[0].i1j 		= cell[i][j].P[0];
				result[0].i_1j_1 	= cell[i-1][j-1].P[0];
				result[0].ij_1 		= cell[i][j-1].P[0];
				result[0].i1j_1 	= cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 	= cell[i-1][j+1].Vx[0];
				result[1].ij1 		= cell[i][j+1].Vx[0];
				result[1].i1j1 		= -cell[i][j].Vx[0];
				result[1].i_1j 		= cell[i-1][j].Vx[0];
				result[1].ij 		= cell[i][j].Vx[0];
				result[1].i1j 		= -cell[i][j].Vx[0];
				result[1].i_1j_1 	= cell[i-1][j-1].Vx[0];
				result[1].ij_1 		= cell[i][j-1].Vx[0];
				result[1].i1j_1 	= -cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 	= -cell[i-1][j].Vr[0];
				result[2].ij1 		= -cell[i][j].Vr[0];
				result[2].i1j1 		= -cell[i][j].Vr[0];
				result[2].i_1j 		= cell[i-1][j].Vr[0];
				result[2].ij 		= cell[i][j].Vr[0];
				result[2].i1j 		= cell[i][j].Vr[0];
				result[2].i_1j_1 	= cell[i-1][j-1].Vr[0];
				result[2].ij_1 		= cell[i][j-1].Vr[0];
				result[2].i1j_1 	= cell[i][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 	= cell[i-1][j+1].e;
				result[3].ij1 		= cell[i][j+1].e;
				result[3].i1j1 		= cell[i][j].e;
				result[3].i_1j 		= cell[i-1][j].e;
				result[3].ij 		= cell[i][j].e;
				result[3].i1j 		= cell[i][j].e;
				result[3].i_1j_1 	= cell[i-1][j-1].e;
				result[3].ij_1 		= cell[i][j-1].e;
				result[3].i1j_1 	= cell[i][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 	= cell[i-1][j+1].rho;
				result[4].ij1 		= cell[i][j+1].rho;
				result[4].i1j1 		= cell[i][j].rho;
				result[4].i_1j 		= cell[i-1][j].rho;
				result[4].ij 		= cell[i][j].rho;
				result[4].i1j 		= cell[i][j].rho;
				result[4].i_1j_1 	= cell[i-1][j-1].rho;
				result[4].ij_1 		= cell[i][j-1].rho;
				result[4].i1j_1 	= cell[i][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 	= cell[i-1][j+1].bar_Vx[0];
				result[5].ij1 		= cell[i][j+1].bar_Vx[0];
				result[5].i1j1 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j 		= cell[i-1][j].bar_Vx[0];
				result[5].ij 		= cell[i][j].bar_Vx[0];
				result[5].i1j 		= -cell[i][j].bar_Vx[0];
				result[5].i_1j_1 	= cell[i-1][j-1].bar_Vx[0];
				result[5].ij_1 		= cell[i][j-1].bar_Vx[0];
				result[5].i1j_1 	= -cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 	= -cell[i-1][j].bar_Vr[0];
				result[6].ij1 		= -cell[i][j].bar_Vr[0];
				result[6].i1j1 		= -cell[i][j].bar_Vr[0];
				result[6].i_1j 		= cell[i-1][j].bar_Vr[0];
				result[6].ij 		= cell[i][j].bar_Vr[0];
				result[6].i1j 		= cell[i][j].bar_Vr[0];
				result[6].i_1j_1 	= cell[i-1][j-1].bar_Vr[0];
				result[6].ij_1 		= cell[i][j-1].bar_Vr[0];
				result[6].i1j_1 	= cell[i][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i][j-1].bar_e;
			}
			if (isSet_z) {
				result[8].i_1j1 = cell[i-1][j+1].bar_z;	result[8].ij1 = cell[i][j+1].bar_z;	result[8].i1j1 = cell[i][j+1].bar_z;
				result[8].i_1j = cell[i-1][j].bar_z;		result[8].ij = cell[i][j].bar_z; 	result[8].i1j = cell[i][j].bar_z;
				result[8].i_1j_1 = cell[i-1][j-1].bar_z;	result[8].ij_1 = cell[i][j-1].bar_z;	result[8].i1j_1 = cell[i][j-1].bar_z;
			}
			if (isSet_psi) {
				result[9].i_1j1 = cell[i-1][j+1].bar_psi;	result[9].ij1 = cell[i][j+1].bar_psi;	result[9].i1j1 = cell[i][j+1].bar_psi;
				result[9].i_1j = cell[i-1][j].bar_psi;		result[9].ij = cell[i][j].bar_psi; 	result[9].i1j = cell[i][j].bar_psi;
				result[9].i_1j_1 = cell[i-1][j-1].bar_psi;	result[9].ij_1 = cell[i][j-1].bar_psi;	result[9].i1j_1 = cell[i][j-1].bar_psi;
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].i1j = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.x.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.x.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.x.at(idx).i;
				weightCell.j = cell[i][j].weightVector.x.at(idx).j;
				if (isSet_P) {
					result[0].i1j += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].i1j -= weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].i1j -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].i1j += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].i1j += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].i1j -= weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].i1j -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].i1j += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].i1j += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].i1j += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].ij1 = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.y.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.y.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.y.at(idx).i;
				weightCell.j = cell[i][j].weightVector.y.at(idx).j;
				if (isSet_P) {
					result[0].ij1 += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].ij1 -= weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].ij1 -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].ij1 += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].ij1 += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].ij1 -= weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].ij1 -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].ij1 += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].ij1 += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].ij1 += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}

			for (unsigned int idx = 0; idx < 10; idx++) {
				result[idx].i1j1 = 0;
			}
			for (unsigned int idx = 0; idx < cell[i][j].weightVector.xy.size(); idx++) {
				Int2D weightCell;
				double weight = cell[i][j].weightVector.xy.at(idx).weight;
				weightCell.i = cell[i][j].weightVector.xy.at(idx).i;
				weightCell.j = cell[i][j].weightVector.xy.at(idx).j;
				if (isSet_P) {
					result[0].i1j1 += weight*cell[weightCell.i][weightCell.j].P[0];
				}
				if (isSet_Vx) {
					result[1].i1j1 -= weight*cell[weightCell.i][weightCell.j].Vx[0];
				}
				if (isSet_Vr) {
					result[2].i1j1 -= weight*cell[weightCell.i][weightCell.j].Vr[0];
				}
				if (isSet_E) {
					result[3].i1j1 += weight*cell[weightCell.i][weightCell.j].e;
				}
				if (isSet_rho) {
					result[4].i1j1 += weight*cell[weightCell.i][weightCell.j].rho;
				}
				if (isSet_barVx) {
					result[5].i1j1 -= weight*cell[weightCell.i][weightCell.j].bar_Vx[0];
				}
				if (isSet_barVr) {
					result[6].i1j1 -= weight*cell[weightCell.i][weightCell.j].bar_Vr[0];
				}
				if (isSet_barE) {
					result[7].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_e;
				}
				if (isSet_z) {
					result[8].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_z;
				}
				if (isSet_psi) {
					result[9].i1j1 += weight*cell[weightCell.i][weightCell.j].bar_psi;
				}
			}
			break;

		default:
			break;
	}
}
