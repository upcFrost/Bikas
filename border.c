#include "border.h"

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



void setVertices(int weightCell, Point2D vertices[4], bool debug) {
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
	Line2D line, LineAngle2D angle, bool debug) {
		
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



std::vector <TPoint2D> getExtPoints(Point2D vertices[4], Int2D vertices_ij[4]) {
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
		if (fabs(fmod(point.x, dx)) < pow(10,-6) && fabs(fmod(point.y, dr)) < pow(10,-6)) continue;
		int idx2 = idx1 != 3 ? idx1 + 1 : 0;
		// One intersection point on Y axis
		if (vertices_ij[idx1].i != vertices_ij[idx2].i && vertices_ij[idx1].j == vertices_ij[idx2].j) {
			double interX = floor(fmax(vertices[idx1].x,vertices[idx2].x)/dx)*dx;
			double interY = vertices[idx2].y + (interX - vertices[idx2].x) * (vertices[idx2].y - vertices[idx1].y) / (vertices[idx2].x - vertices[idx1].x);
			point.x = interX; point.y = interY; point.type = 1;
			result.push_back(point);
		} else
		// One intersection point on X axis
		if (vertices_ij[idx1].i == vertices_ij[idx2].i && vertices_ij[idx1].j != vertices_ij[idx2].j) {
			unsigned int idx_max = vertices[idx1].y > vertices[idx2].y ? idx1 : idx2;
			double interY = floor(vertices[idx_max].y/dr)*dr;
			double interX = vertices[idx2].x + (interY - vertices[idx2].y) * (vertices[idx2].x - vertices[idx1].x) / (vertices[idx2].y - vertices[idx1].y);
			point.x = interX; point.y = interY; point.type = 1;
			result.push_back(point);
		} else
		// If has 2 intersection result - first serve closest (because triangle gen is buggy)
		if (vertices_ij[idx1].i != vertices_ij[idx2].i && vertices_ij[idx1].j != vertices_ij[idx2].j) {
			double interX1 = floor(fmax(vertices[idx1].x,vertices[idx2].x)/dx)*dx;
			double interY1 = vertices[idx2].y + (interX1 - vertices[idx2].x) * (vertices[idx2].y - vertices[idx1].y) / (vertices[idx2].x - vertices[idx1].x);
			double interY2 = floor(fmax(vertices[idx1].y,vertices[idx2].y)/dr)*dr;
			double interX2 = vertices[idx2].x + (interY2 - vertices[idx2].y) * (vertices[idx2].x - vertices[idx1].x) / (vertices[idx2].y - vertices[idx1].y);
			if (pow(vertices[idx1].x-interX1,2) + pow(vertices[idx1].y-interY1,2) < pow(vertices[idx1].x-interX2,2) + pow(vertices[idx1].y-interY2,2)) {
				point.x = interX1; point.y = interY1; point.type = 1;
				result.push_back(point);
				point.x = interX2; point.y = interY2 ; point.type = 1;
				result.push_back(point);
			} else {
				point.x = interX2; point.y = interY2; point.type = 1;
				result.push_back(point);
				point.x = interX1; point.y = interY1; point.type = 1;
				result.push_back(point);
			}
		}
	}
	
	return result;
}



Int2D getMaxDifference(unsigned int max_i_point[2], unsigned int max_j_point[2],
		Int2D vertices_ij[4], bool debug) {
			
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



std::vector <TPoint2D> getIntPoints(std::vector <TPoint2D> points,
		Point2D vertices[4], Int2D max_diff,
		unsigned int max_i_point[2], unsigned int max_j_point[2]) {
			
	TPoint2D point;
	
	if (max_diff.i == 1 && max_diff.j == 1) {
		double internalX = floor(fmax(fmax(vertices[0].x,vertices[1].x),fmax(vertices[2].x,vertices[3].x))/dx)*dx;
		double internalY = floor(fmax(fmax(vertices[0].y,vertices[1].y),fmax(vertices[2].y,vertices[3].y))/dr)*dr;
		point.x = internalX; point.y = internalY; point.type = 2;
		points.push_back(point);
	} else {
		unsigned int prevPointI[2] = {0};
		unsigned int nextPointI[2] = {0};
		unsigned int prevPointJ[2] = {0};
		unsigned int nextPointJ[2] = {0};
		prevPointI[0] = max_i_point[0] == 0 ? points.size()-1 : max_i_point[0]-1;
		prevPointI[1] = max_i_point[1] == 0 ? points.size()-1 : max_i_point[1]-1;
		nextPointI[0] = max_i_point[0] == points.size()-1 ? 0 : max_i_point[0]+1;
		nextPointI[1] = max_i_point[1] == points.size()-1 ? 0 : max_i_point[1]+1;
		prevPointJ[0] = max_j_point[0] == 0 ? points.size()-1 : max_j_point[0]-1;
		prevPointJ[1] = max_j_point[1] == 0 ? points.size()-1 : max_j_point[1]-1;
		nextPointJ[0] = max_j_point[0] == points.size()-1 ? 0 : max_j_point[0]+1;
		nextPointJ[1] = max_j_point[1] == points.size()-1 ? 0 : max_j_point[1]+1;
		
		unsigned int max_i_difference = 0;
		unsigned int max_j_difference = 0;
		if (fabs(floor(points.at(nextPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx)) > max_i_difference) 
				max_i_difference = fabs(floor(points.at(nextPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx));
		if (fabs(floor(points.at(prevPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx)) > max_i_difference) 
				max_i_difference = fabs(floor(points.at(prevPointI[0]).x/dx) - floor(points.at(max_i_point[1]).x/dx));
		if (fabs(floor(points.at(nextPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr)) > max_j_difference) 
				max_j_difference = fabs(floor(points.at(nextPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr));
		if (fabs(floor(points.at(prevPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr)) > max_j_difference) 
				max_j_difference = fabs(floor(points.at(prevPointJ[0]).y/dr) - floor(points.at(max_j_point[1]).y/dr));
		
		//if (max_i_difference == 1 && max_j_difference == 1) {
			double internalX = floor(fmax(fmax(vertices[0].x,vertices[1].x),fmax(vertices[2].x,vertices[3].x))/dx)*dx;
			double internalY = floor(fmax(fmax(vertices[0].y,vertices[1].y),fmax(vertices[2].y,vertices[3].y))/dr)*dr;
			point.x = internalX; point.y = internalY; point.type = 2;
			points.push_back(point);
		//}
		
		//~ cout << "max_i_difference = " << max_i_difference << endl;
		//~ cout << "max_j_difference = " << max_j_difference << endl;
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
			if (debug) printf("Cell %u: %4.4f:%4.4f, type %d \n", 
			    idx,points.at(idx2).x,points.at(idx2).y,points.at(idx2).type);
		}
	}
	
	return pointsInCell;
}


std::vector <TPoint2D> fixIntPointID(std::vector <TPoint2D> pointsInCell,
		bool debug) {
	
	if (debug) printf("We have an internal point!\n");
	unsigned int intIdx = 0;
	TPoint2D intPoint;
	for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
		if (pointsInCell.at(idx2).type == 2) {
			intIdx = idx2;
			intPoint = pointsInCell.at(idx2);
			break;
		}
	}
	double minHor = 1000, minVer = 1000;
	int minHorIdx = 0, minVerIdx = 0;
	for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
		if (idx2 != intIdx) {
			TPoint2D point = pointsInCell.at(idx2);
			if (debug) printf("Testing point %u\n", idx2);
			
			if (fabs(point.x - intPoint.x) < minHor && fabs(point.y - intPoint.y) < pow(10,-6)
			    && pointsInCell.at(idx2).type == 1) {
				if (debug) printf("hor = %4.4f\n", point.x - intPoint.x);
				minHor = point.x - intPoint.x;
				minHorIdx = idx2;
			}
			
			if (fabs(point.y - intPoint.y) < minVer && fabs(point.x - intPoint.x) < pow(10,-6)
			    && pointsInCell.at(idx2).type == 1) {
				if (debug) printf("ver = %4.4f\n", point.x - intPoint.x);
				minVer = point.y - intPoint.y;
				minVerIdx = idx2;
			}
		}
	}
	if (debug) 
		printf("Closest points to internal: hor %4.4f:%4.4f, ver %4.4f:%4.4f\n\n", 
			pointsInCell.at(minHorIdx).x,pointsInCell.at(minHorIdx).y,
			pointsInCell.at(minVerIdx).x,pointsInCell.at(minVerIdx).y);
	
	if (minVerIdx != 0 && minHorIdx != 0) {
	std::vector<TPoint2D>::iterator it = pointsInCell.begin();
	pointsInCell.erase(it+intIdx);
		if (minVerIdx > minHorIdx) {
			pointsInCell.insert(it+minVerIdx, intPoint);
		} else {
			pointsInCell.insert(it+minHorIdx, intPoint);
		}
	}
	
	if (debug) printf("Points in current cell in order: ");
	for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
		if (debug) printf("%4.4f:%4.4f, ",
			pointsInCell.at(idx2).x, pointsInCell.at(idx2).y);
	}
	if (debug) printf("\n");
	
	return pointsInCell;
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
	double area;
	std::vector <Int2D> cells;
	
	for (unsigned int weightCell = 0; weightCell < 3; weightCell++) {
		// First, we'll set original cell's vertices
		setVertices(weightCell, vertices, debug);
		
		// Determining line begin-end points and angle
		getLineAngle(cell, i, j, n, line, angle, debug);
		cell.at(n).at(i).at(j).angle = angle;
		
		// Getting mirrored over line vertices array
		getMirrorVerts(vertices, vertices_ij, line, angle, debug);
		
		// Getting all points of intersection between mirrored cell and grid
		points = getExtPoints(vertices, vertices_ij);
		
		// Getting maximum i and j difference between mirrored vertices
		max_diff = getMaxDifference(max_i_point, max_j_point,
			vertices_ij, debug);
			
		// Cleaning points from doubles (occures sometimes)
		points = clearDoubles(points);
		
		// Getting internal points (grid's own intersections)
		points = getIntPoints(points, vertices, max_diff, 
			max_i_point, max_j_point);
			
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
		if (debug) printf("\nNow will arrange points in those cells\n");
		
		// Main loop where each cell's weight is calculated
		double totalWeight = 0;
		for (unsigned int idx = 0; idx < cells.size(); idx++) {
			// We'll use center of each of those cells to determine how much points we have in each cell
			std::vector <TPoint2D> pointsInCell = getPointsInCell(idx, 
				points, cells, debug);
			// If we have internal point - fix ID
			if (std::find_if(pointsInCell.begin(), pointsInCell.end(), 
					find_type(2)) != pointsInCell.end()) {
				pointsInCell = fixIntPointID(pointsInCell, debug);
			}
			// Some nice debug
			if (debug) printf("PointsArray size: %u\n", 
				(unsigned int)pointsInCell.size());
			// Now we'll use the following arrays in polygonArea function
			double XPointsArray[pointsInCell.size()];
			double YPointsArray[pointsInCell.size()];
			for (unsigned int idx2 = 0; idx2 < pointsInCell.size(); idx2++) {
				XPointsArray[idx2] = pointsInCell.at(idx2).x;
				YPointsArray[idx2] = pointsInCell.at(idx2).y;
				if (debug) printf("Point %4.4f:%4.4f;  ",
					XPointsArray[idx2],YPointsArray[idx2]);
			}
			area = fabs(polygonArea(XPointsArray, YPointsArray, 
				pointsInCell.size()));
			if (debug) printf("Polygon area: %8.8f\n", area);
			// Populating weightPart and pushing it to the cell's weightVector
			WeightPart weightPart;
			weightPart.weight = area / origArea;
			weightPart.i = cells.at(idx).i;
			weightPart.j = cells.at(idx).j;
			if (debug) printf("Weight of this cell: %4.4f\n\n\n", 
				weightPart.weight);
			if (weightCell == 0) {
				result.x.push_back(weightPart);
			} else if (weightCell == 1) {
				result.y.push_back(weightPart);
			} else if (weightCell == 2) {
				result.xy.push_back(weightPart);
			}
			totalWeight += weightPart.weight;
		}
		// Scaling to 1
		double scale = 1.0/totalWeight;
		if (weightCell == 0) {
			for (unsigned int idx = 0; idx < result.x.size(); idx++) {
				result.x.at(idx).weight *= scale;
			}
		} else if (weightCell == 1) {
			for (unsigned int idx = 0; idx < result.y.size(); idx++) {
				result.y.at(idx).weight *= scale;
			}
		} else if (weightCell == 2) {
			for (unsigned int idx = 0; idx < result.xy.size(); idx++) {
				result.xy.at(idx).weight *= scale;
			}
		}
	}

	if (debug) printf("Total weight parts at %d in cell %d:%d\nby x: %u\nby y: %u\nby xy: %u\n", n,i,j,
		(unsigned int) result.x.size(),
		(unsigned int) result.y.size(),
		(unsigned int) result.xy.size());

	return result;
}



void calculateBorder(int n, cell2dStatic& cell, unsigned long ctrl, BorderCond result[8]) {

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

	switch (type) {
		// Free cell
		case 0:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j+1].P[0]; 		result[0].ij1 = cell[i][j+1].P[0];		result[0].i1j1 = cell[i+1][j+1].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];		result[1].i1j1 = cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 		result[1].ij = cell[i][j].Vx[0]; 		result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];		result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 		result[2].ij = cell[i][j].Vr[0]; 		result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e; 		result[3].ij1 = cell[i][j+1].e;			result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i-1][j].e; 			result[3].ij = cell[i][j].e; 			result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;		result[4].ij1 = cell[i][j+1].rho;		result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;			result[4].ij = cell[i][j].rho; 			result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;		result[4].ij_1 = cell[i][j-1].rho;		result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];	result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];result[5].ij_1 = cell[i][j-1].bar_Vx[0];result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];	result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];result[6].ij_1 = cell[i][j-1].bar_Vr[0];result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;		result[7].ij1 = cell[i][j+1].bar_e;		result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 		result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			break;

		case 1:
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

			for (unsigned int idx = 0; idx < 8; idx++) {
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
			}

			for (unsigned int idx = 0; idx < 8; idx++) {
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
			}
			break;

		case 3:
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
			break;


		case 10:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j].P[0]; 		result[0].ij1 = cell[i][j].P[0];		result[0].i1j1 = cell[i][j].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];		result[1].i1j1 = cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 		result[1].ij = cell[i][j].Vx[0]; 		result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];		result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 		result[2].ij = cell[i][j].Vr[0]; 		result[2].i1j = cell[i+1][j].Vr[0];
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

			for (unsigned int idx = 0; idx < 8; idx++) {
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
			}

			for (unsigned int idx = 0; idx < 8; idx++) {
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
			}

			for (unsigned int idx = 0; idx < 8; idx++) {
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
			}
			break;

		// Top border closed
		case 13:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j].P[0]; 	result[0].ij1 = cell[i][j].P[0];	result[0].i1j1 = cell[i+1][j].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j].Vx[0]; 	result[1].ij1 = cell[i][j].Vx[0];	result[1].i1j1 = cell[i+1][j].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j].Vr[0]; 	result[2].ij1 = cell[i][j].Vr[0];	result[2].i1j1 = cell[i+1][j].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j].e; 		result[3].ij1 = cell[i][j].e;		result[3].i1j1 = cell[i+1][j].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j].rho;	result[4].ij1 = cell[i][j].rho;	result[4].i1j1 = cell[i+1][j].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j].bar_Vx[0];	result[5].ij1 = cell[i][j].bar_Vx[0];	result[5].i1j1 = cell[i+1][j].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];	result[5].ij_1 = cell[i][j-1].bar_Vx[0];	result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j].bar_Vr[0];	result[6].ij1 = cell[i][j].bar_Vr[0];	result[6].i1j1 = cell[i+1][j].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j].bar_e;	result[7].ij1 = cell[i][j].bar_e;	result[7].i1j1 = cell[i+1][j].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			break;

		// Bottom border closed
		case 14:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j+1].P[0]; 	result[0].ij1 = cell[i][j+1].P[0];	result[0].i1j1 = cell[i+1][j+1].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i-1][j].P[0]; 	result[0].ij_1 = cell[i][j].P[0]; 	result[0].i1j_1 = cell[i+1][j].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];	result[1].i1j1 = cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j].Vx[0];	result[1].ij_1 = cell[i][j].Vx[0];	result[1].i1j_1 = cell[i+1][j].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];	result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j].Vr[0];	result[2].ij_1 = cell[i][j].Vr[0];	result[2].i1j_1 = cell[i+1][j].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i-1][j].e;		result[3].ij_1 = cell[i][j].e;		result[3].i1j_1 = cell[i+1][j].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i-1][j].rho;	result[4].ij_1 = cell[i][j].rho;	result[4].i1j_1 = cell[i+1][j].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j].bar_Vx[0];	result[5].ij_1 = cell[i][j].bar_Vx[0];	result[5].i1j_1 = cell[i+1][j].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j].bar_Vr[0];	result[6].ij_1 = cell[i][j].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j].bar_e;	result[7].ij_1 = cell[i][j].bar_e;	result[7].i1j_1 = cell[i+1][j].bar_e;
			}
			break;

		// Top and left borders closed
		case 15:
			if (isSet_P) {
				result[0].i_1j1 = cell[i][j].P[0]; 	result[0].ij1 = cell[i][j].P[0];	result[0].i1j1 = cell[i+1][j].P[0];
				result[0].i_1j = cell[i][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i][j].Vx[0]; 	result[1].ij1 = cell[i][j].Vx[0];	result[1].i1j1 = cell[i+1][j].Vx[0];
				result[1].i_1j = cell[i][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i][j].Vr[0]; 	result[2].ij1 = cell[i][j].Vr[0];	result[2].i1j1 = cell[i+1][j].Vr[0];
				result[2].i_1j = cell[i][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i][j].e; 		result[3].ij1 = cell[i][j].e;		result[3].i1j1 = cell[i+1][j].e;
				result[3].i_1j = cell[i][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i][j].rho;	result[4].ij1 = cell[i][j].rho;	result[4].i1j1 = cell[i+1][j].rho;
				result[4].i_1j = cell[i][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i][j].bar_Vx[0];	result[5].ij1 = cell[i][j].bar_Vx[0];	result[5].i1j1 = cell[i+1][j].bar_Vx[0];
				result[5].i_1j = cell[i][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i][j-1].bar_Vx[0];	result[5].ij_1 = cell[i][j-1].bar_Vx[0];	result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i][j].bar_Vr[0];	result[6].ij1 = cell[i][j].bar_Vr[0];	result[6].i1j1 = cell[i+1][j].bar_Vr[0];
				result[6].i_1j = cell[i][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i][j].bar_e;	result[7].ij1 = cell[i][j].bar_e;	result[7].i1j1 = cell[i+1][j].bar_e;
				result[7].i_1j = cell[i][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			break;

		// Bottom and left borders closed
		case 16:
			if (isSet_P) {
				result[0].i_1j1 = cell[i][j+1].P[0]; 	result[0].ij1 = cell[i][j+1].P[0];	result[0].i1j1 = cell[i+1][j+1].P[0];
				result[0].i_1j = cell[i][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i][j].P[0]; 	result[0].ij_1 = cell[i][j].P[0]; 	result[0].i1j_1 = cell[i+1][j].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];	result[1].i1j1 = cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i][j].Vx[0];	result[1].ij_1 = cell[i][j].Vx[0];	result[1].i1j_1 = cell[i+1][j].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];	result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i][j].Vr[0];	result[2].ij_1 = cell[i][j].Vr[0];	result[2].i1j_1 = cell[i+1][j].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i][j].e;		result[3].ij_1 = cell[i][j].e;		result[3].i1j_1 = cell[i+1][j].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i][j].rho;	result[4].ij_1 = cell[i][j].rho;	result[4].i1j_1 = cell[i+1][j].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j = cell[i][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i][j].bar_Vx[0];	result[5].ij_1 = cell[i][j].bar_Vx[0];	result[5].i1j_1 = cell[i+1][j].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i][j].bar_Vr[0];	result[6].ij_1 = cell[i][j].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i][j].bar_e;	result[7].ij_1 = cell[i][j].bar_e;	result[7].i1j_1 = cell[i+1][j].bar_e;
			}
			break;

		// Left border closed
		case 17:
			if (isSet_P) {
				result[0].i_1j1 = cell[i][j+1].P[0]; 	result[0].ij1 = cell[i][j+1].P[0];	result[0].i1j1 = cell[i+1][j+1].P[0];
				result[0].i_1j = cell[i][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i+1][j].P[0];
				result[0].i_1j_1 = cell[i][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i+1][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];	result[1].i1j1 = cell[i+1][j+1].Vx[0];
				result[1].i_1j = cell[i][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i+1][j].Vx[0];
				result[1].i_1j_1 = cell[i][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i+1][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];	result[2].i1j1 = cell[i+1][j+1].Vr[0];
				result[2].i_1j = cell[i][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i+1][j].Vr[0];
				result[2].i_1j_1 = cell[i][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i+1][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i+1][j+1].e;
				result[3].i_1j = cell[i][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i+1][j].e;
				result[3].i_1j_1 = cell[i][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i+1][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i+1][j+1].rho;
				result[4].i_1j = cell[i][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i+1][j].rho;
				result[4].i_1j_1 = cell[i][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i+1][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i+1][j+1].bar_Vx[0];
				result[5].i_1j = cell[i][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i+1][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i][j-1].bar_Vx[0];	result[5].ij_1 = cell[i][j-1].bar_Vx[0];	result[5].i1j_1 = cell[i+1][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i+1][j+1].bar_Vr[0];
				result[6].i_1j = cell[i][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i+1][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i+1][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i+1][j+1].bar_e;
				result[7].i_1j = cell[i][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i+1][j].bar_e;
				result[7].i_1j_1 = cell[i][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i+1][j-1].bar_e;
			}
			break;

		// Right border closed
		case 19:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j+1].P[0]; 	result[0].ij1 = cell[i][j+1].P[0];	result[0].i1j1 = cell[i][j+1].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];	result[1].i1j1 = cell[i][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];	result[2].i1j1 = cell[i][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i][j+1].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i][j+1].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];	result[5].ij_1 = cell[i][j-1].bar_Vx[0];	result[5].i1j_1 = cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i][j-1].bar_e;
			}
			break;

		// Top and right borders closed
		case 20:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j].P[0]; 	result[0].ij1 = cell[i][j].P[0];	result[0].i1j1 = cell[i][j].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i][j].P[0];
				result[0].i_1j_1 = cell[i-1][j-1].P[0]; 	result[0].ij_1 = cell[i][j-1].P[0]; 	result[0].i1j_1 = cell[i][j-1].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j].Vx[0]; 	result[1].ij1 = cell[i][j].Vx[0];	result[1].i1j1 = cell[i][j].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j-1].Vx[0];	result[1].ij_1 = cell[i][j-1].Vx[0];	result[1].i1j_1 = cell[i][j-1].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j].Vr[0]; 	result[2].ij1 = cell[i][j].Vr[0];	result[2].i1j1 = cell[i][j].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j-1].Vr[0];	result[2].ij_1 = cell[i][j-1].Vr[0];	result[2].i1j_1 = cell[i][j-1].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j].e; 		result[3].ij1 = cell[i][j].e;		result[3].i1j1 = cell[i][j].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i][j].e;
				result[3].i_1j_1 = cell[i-1][j-1].e;		result[3].ij_1 = cell[i][j-1].e;		result[3].i1j_1 = cell[i][j-1].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j].rho;	result[4].ij1 = cell[i][j].rho;	result[4].i1j1 = cell[i][j].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i][j].rho;
				result[4].i_1j_1 = cell[i-1][j-1].rho;	result[4].ij_1 = cell[i][j-1].rho;	result[4].i1j_1 = cell[i][j-1].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j].bar_Vx[0];	result[5].ij1 = cell[i][j].bar_Vx[0];	result[5].i1j1 = cell[i][j].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j-1].bar_Vx[0];	result[5].ij_1 = cell[i][j-1].bar_Vx[0];	result[5].i1j_1 = cell[i][j-1].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j].bar_Vr[0];	result[6].ij1 = cell[i][j].bar_Vr[0];	result[6].i1j1 = cell[i][j].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j-1].bar_Vr[0];	result[6].ij_1 = cell[i][j-1].bar_Vr[0];	result[6].i1j_1 = cell[i][j-1].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j].bar_e;	result[7].ij1 = cell[i][j].bar_e;	result[7].i1j1 = cell[i][j].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j-1].bar_e;	result[7].ij_1 = cell[i][j-1].bar_e;	result[7].i1j_1 = cell[i][j-1].bar_e;
			}
			break;

		// Bottom and right borders closed
		case 21:
			if (isSet_P) {
				result[0].i_1j1 = cell[i-1][j+1].P[0]; 	result[0].ij1 = cell[i][j+1].P[0];	result[0].i1j1 = cell[i][j+1].P[0];
				result[0].i_1j = cell[i-1][j].P[0]; 		result[0].ij = cell[i][j].P[0]; 		result[0].i1j = cell[i][j].P[0];
				result[0].i_1j_1 = cell[i-1][j].P[0]; 	result[0].ij_1 = cell[i][j].P[0]; 	result[0].i1j_1 = cell[i][j].P[0];
			}
			if (isSet_Vx) {
				result[1].i_1j1 = cell[i-1][j+1].Vx[0]; 	result[1].ij1 = cell[i][j+1].Vx[0];	result[1].i1j1 = cell[i][j+1].Vx[0];
				result[1].i_1j = cell[i-1][j].Vx[0]; 	result[1].ij = cell[i][j].Vx[0]; 	result[1].i1j = cell[i][j].Vx[0];
				result[1].i_1j_1 = cell[i-1][j].Vx[0];	result[1].ij_1 = cell[i][j].Vx[0];	result[1].i1j_1 = cell[i][j].Vx[0];
			}
			if (isSet_Vr) {
				result[2].i_1j1 = cell[i-1][j+1].Vr[0]; 	result[2].ij1 = cell[i][j+1].Vr[0];	result[2].i1j1 = cell[i][j+1].Vr[0];
				result[2].i_1j = cell[i-1][j].Vr[0]; 	result[2].ij = cell[i][j].Vr[0]; 	result[2].i1j = cell[i][j].Vr[0];
				result[2].i_1j_1 = cell[i-1][j].Vr[0];	result[2].ij_1 = cell[i][j].Vr[0];	result[2].i1j_1 = cell[i][j].Vr[0];
			}
			if (isSet_E) {
				result[3].i_1j1 = cell[i-1][j+1].e; 		result[3].ij1 = cell[i][j+1].e;		result[3].i1j1 = cell[i][j+1].e;
				result[3].i_1j = cell[i-1][j].e; 		result[3].ij = cell[i][j].e; 		result[3].i1j = cell[i][j].e;
				result[3].i_1j_1 = cell[i-1][j].e;		result[3].ij_1 = cell[i][j].e;		result[3].i1j_1 = cell[i][j].e;
			}
			if (isSet_rho) {
				result[4].i_1j1 = cell[i-1][j+1].rho;	result[4].ij1 = cell[i][j+1].rho;	result[4].i1j1 = cell[i][j+1].rho;
				result[4].i_1j = cell[i-1][j].rho;		result[4].ij = cell[i][j].rho; 		result[4].i1j = cell[i][j].rho;
				result[4].i_1j_1 = cell[i-1][j].rho;	result[4].ij_1 = cell[i][j].rho;	result[4].i1j_1 = cell[i][j].rho;
			}
			if (isSet_barVx) {
				result[5].i_1j1 = cell[i-1][j+1].bar_Vx[0];	result[5].ij1 = cell[i][j+1].bar_Vx[0];	result[5].i1j1 = cell[i][j+1].bar_Vx[0];
				result[5].i_1j = cell[i-1][j].bar_Vx[0];		result[5].ij = cell[i][j].bar_Vx[0]; 	result[5].i1j = cell[i][j].bar_Vx[0];
				result[5].i_1j_1 = cell[i-1][j].bar_Vx[0];	result[5].ij_1 = cell[i][j].bar_Vx[0];	result[5].i1j_1 = cell[i][j].bar_Vx[0];
			}
			if (isSet_barVr) {
				result[6].i_1j1 = cell[i-1][j+1].bar_Vr[0];	result[6].ij1 = cell[i][j+1].bar_Vr[0];	result[6].i1j1 = cell[i][j+1].bar_Vr[0];
				result[6].i_1j = cell[i-1][j].bar_Vr[0];		result[6].ij = cell[i][j].bar_Vr[0]; 	result[6].i1j = cell[i][j].bar_Vr[0];
				result[6].i_1j_1 = cell[i-1][j].bar_Vr[0];	result[6].ij_1 = cell[i][j].bar_Vr[0];	result[6].i1j_1 = cell[i][j].bar_Vr[0];
			}
			if (isSet_barE) {
				result[7].i_1j1 = cell[i-1][j+1].bar_e;	result[7].ij1 = cell[i][j+1].bar_e;	result[7].i1j1 = cell[i][j+1].bar_e;
				result[7].i_1j = cell[i-1][j].bar_e;		result[7].ij = cell[i][j].bar_e; 	result[7].i1j = cell[i][j].bar_e;
				result[7].i_1j_1 = cell[i-1][j].bar_e;	result[7].ij_1 = cell[i][j].bar_e;	result[7].i1j_1 = cell[i][j].bar_e;
			}
			break;
	}
}
