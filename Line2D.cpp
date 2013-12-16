/*
 * Line2D.cpp
 *
 *  Created on: 16.12.2013
 *      Author: PBelyaev
 */

#include "include/Line2D.hpp"

Line2D::Line2D() {

}

Line2D::Line2D(double x1, double x2, double y1, double y2) {
	setLine(x1, x2, y1, y2);
	setSinA();
	setCosA();
	setSin2A();
	setCos2A();
}

void Line2D::setXBegin(double x) {
	xbegin = x;
}

void Line2D::setXEnd(double x) {
	xend = x;
}

void Line2D::setYBegin(double y) {
	ybegin = y;
}

void Line2D::setYEnd(double y) {
	yend = y;
}

void Line2D::setX(double x1, double x2) {
	xbegin = x1;
	xend = x2;
}

void Line2D::setY(double y1, double y2) {
	ybegin = y1;
	yend = y2;
}

void Line2D::setLine(double x1, double x2, double y1, double y2) {
	xbegin = x1;
	xend = x2;
	ybegin = y1;
	yend = y2;
}

void Line2D::setSinA() {
	sin_a = (yend - ybegin)
			/ std::sqrt(pow(yend - ybegin, 2) + std::pow(xend - xbegin, 2));
}

void Line2D::setCosA() {
	cos_a = (xend - xbegin)
			/ std::sqrt(pow(yend - ybegin, 2) + std::pow(xend - xbegin, 2));
}

void Line2D::setSin2A() {
	sin_2a = 2 * cos_a * sin_a;
}

void Line2D::setCos2A() {
	cos_2a = std::pow(cos_a, 2) - std::pow(sin_a, 2);
}

