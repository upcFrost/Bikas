#ifndef BORDER_H
#define	BORDER_H

#include "types.h"
#include "globals.h"
#include "triangulate.h"
#include <vector>
#include <stdbool.h>

typedef struct TPoint2D { double x; double y; int type; } TPoint2D;

typedef struct CPoint2D { TPoint2D point; double cot; } CPoint2D;

struct lesserCot
{
    inline bool operator() (const CPoint2D& cPoint1, const CPoint2D& cPoint2)
    {
        return (cPoint1.cot > cPoint2.cot);
    }
};

struct find_type : std::unary_function<TPoint2D, bool> {
    int type;
    find_type(int type):type(type) { }
    bool operator()(TPoint2D const& p) const {
        return p.type == type;
    }
};

WeightVector wightVectorsCalc(cell2d& cell, int i, int j, int n, bool debug);
void calculateBorder(int n, const cell2dStatic& cell, unsigned long ctrl,
		BorderCond (&result)[10], int i, int j);
void pre_cell_geometry(double array[5], gasCell cell, int i, int j);

#endif	/* BORDER_H */
