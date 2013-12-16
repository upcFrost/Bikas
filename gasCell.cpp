/*
 * gasCell.cpp
 *
 *  Created on: 16.12.2013
 *      Author: PBelyaev
 */

GasCell::GasCell() {}

void GasCell::SetGeometry(int type, double r_1, double r_2, double x_1, double x_2) {
	this->type = type;
	this->x_1 = x_1;
	this->x_2 = x_2;
	this->r_1 = r_1;
	this->r_2 = r_2;
}

void GasCell::SetLine(Line2D line) {
	this->line = line;
}

void GasCell::SetBarVx(double bar_Vx) {
	this->bar_Vx = bar_Vx;
}

void GasCell::SetVx(double Vx) {
	this->Vx = Vx;
}

void GasCell::SetBarVr(double bar_Vr) {
	this->bar_Vr = bar_Vr;
}

void GasCell::SetVr(double Vr) {
	this->Vr = Vr;
}

void GasCell::SetP(double P) {
	this->P = P;
}

void GasCell::SetEmptyCell() {
	P = 0;
	rho = 0;
	e = 0;
	bar_Vx = 0;
	Vx = 0;
	Vr = 0;
	final_z = 0;
	final_psi = 0;
}
