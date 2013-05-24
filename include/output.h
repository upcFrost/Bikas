#ifndef OUTPUT_H
#define	OUTPUT_H


#include "globals.h"
#include "types.h"
#include <vector>

#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

void OutputPVD(cell2d cell, std::string filename);
void outputCSV(cell2d cell, std::ofstream & outputGas);
void prepOutputDynCSV(std::ofstream & outputDyn);
void prepOutputGasCSV(std::ofstream & outputGas, bool verbose);
void outputDynCSV(std::ofstream & outputDyn, double t, int i_sn, double x_sn, double U_sn, double boltP, double projP);

#endif	/* OUTPUT_H */
