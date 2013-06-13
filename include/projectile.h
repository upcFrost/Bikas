/*
 * projectile.h
 *
 *  Created on: May 22, 2013
 *      Author: frost
 */

#ifndef PROJECTILE_H_
#define PROJECTILE_H_

#include <cmath>
#include "globals.h"
#include "border.h"
#include "functions.h"
#include "debug.h"

void projCalc(cell2d & cell, int var, int borderI,
		bool PROJECTILE, bool debug);

void pistonCalc(cell2d & cell, int borderI_prev,
		int borderI, int gasVar, bool debug);

#endif /* PROJECTILE_H_ */
