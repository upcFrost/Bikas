#ifndef BORDER_H
#define	BORDER_H

#include "types.h"
#include "globals.h"
#include <vector>
#include <stdbool.h>

typedef struct TPoint2D { double x; double y; int type; } TPoint2D;

struct find_type : std::unary_function<TPoint2D, bool> {
    int type;
    find_type(int type):type(type) { }
    bool operator()(TPoint2D const& p) const {
        return p.type == type;
    }
};

WeightVector wightVectorsCalc(cell2d& cell, int i, int j, int n, bool debug);
void calculateBorder(int n, cell2dStatic& cell, unsigned long ctrl,
		BorderCond result[8], int i, int j);
void pre_cell_geometry(double array[5], gasCell cell, int i, int j);

#endif	/* BORDER_H */
