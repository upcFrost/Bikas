/*
 * init.h
 *
 *  Created on: May 20, 2013
 *      Author: frost
 */

#ifndef INIT_H_
#define INIT_H_

#include "globals.h"
#include "border.h"
#include "math_functions.h"
#include "Line2D.hpp"

void init(std::ifstream & inputFile, cell2d & cell, int var,
		bool havePiston, bool debug);

void borderCellsFix(cell2d & cell, bool havePiston);


#endif /* INIT_H_ */
