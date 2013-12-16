/*
 * Line2D.hpp
 *
 *  Created on: 16.12.2013
 *      Author: PBelyaev
 */

#ifndef LINE2D_HPP_
#define LINE2D_HPP_

#include <cmath>

class Line2D {
public:
	double xbegin, xend, ybegin, yend;
	double sin_a, cos_a, sin_2a, cos_2a;
	Line2D ();
	Line2D (double, double, double, double);
	void setXBegin(double);
	void setXEnd(double);
	void setYBegin(double);
	void setYEnd(double);
	void setX(double, double);
	void setY(double, double);

	void setLine(double, double, double, double);
private:
	void setSinA();
	void setCosA();
	void setSin2A();
	void setCos2A();
};

#endif /* LINE2D_HPP_ */
