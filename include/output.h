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

typedef struct linkID {
	vtkIdType id[5];
	int i;
	int j;
} linkID;

class Output {
public:
	Output(cell2dStatic & cell);

	void MakePVDOutput(cell2d & cell, std::string & filename);
	void MakeCSVOutput(cell2d & cell, std::ofstream & outputGas);
	void MakeDynOutput(std::ofstream & outputDyn, double t, int i_sn,
			double x_sn, double U_sn, double boltP, double projP);

private:
	// Points array
	vtkSmartPointer<vtkPoints> points;
	// Vertices array
	vtkSmartPointer<vtkCellArray> verts;
	// Cells array
	vtkSmartPointer<vtkCellArray> strips;
	// General data array
	vtkSmartPointer<vtkPolyData> polydata;
	// Data arrays
	vtkSmartPointer<vtkDoubleArray> outputP;
	vtkSmartPointer<vtkDoubleArray> outputE;
	vtkSmartPointer<vtkDoubleArray> outputZ;
	vtkSmartPointer<vtkDoubleArray> outputRho;
	vtkSmartPointer<vtkDoubleArray> outputV;
	// ID-point conversion vector
	std::vector<linkID> idVector;

	void initPVDArrays(cell2dStatic & cell);
};

//void OutputPVD(cell2d & cell, std::string & filename);
void outputCSV(cell2d & cell, std::ofstream & outputGas);
void prepOutputDynCSV(std::ofstream & outputDyn);
void prepOutputGasCSV(std::ofstream & outputGas, bool verbose);
void outputDynCSV(std::ofstream & outputDyn, double t, int i_sn, double x_sn, double U_sn, double boltP, double projP);

#endif	/* OUTPUT_H */
