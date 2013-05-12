#include <stdio.h>
#include "output.h"

void OutputPVD(cell2d cell, std::string filename) {
	// Points array
	vtkSmartPointer<vtkPoints> points = 
	    vtkSmartPointer<vtkPoints>::New();
	// General polydata array
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	// Data arrays
	vtkSmartPointer<vtkDoubleArray> outputP = vtkSmartPointer<vtkDoubleArray>::New();
	outputP->SetName("P");
	vtkSmartPointer<vtkDoubleArray> outputE = vtkSmartPointer<vtkDoubleArray>::New();
	outputE->SetName("E");
	vtkSmartPointer<vtkDoubleArray> outputZ = vtkSmartPointer<vtkDoubleArray>::New();
	outputZ->SetName("Z");
	vtkSmartPointer<vtkDoubleArray> outputRho = vtkSmartPointer<vtkDoubleArray>::New();
	outputRho->SetName("Rho");
	vtkSmartPointer<vtkDoubleArray> outputV = vtkSmartPointer<vtkDoubleArray>::New();
	outputV->SetName("V");
	outputV->SetNumberOfComponents(3);
	outputV->SetComponentName(0, "Vx");
	outputV->SetComponentName(1, "Vr");
	outputV->SetComponentName(2, "Vz");
	// Vertices array
	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	// Cells array
	vtkSmartPointer<vtkCellArray> strips = vtkSmartPointer<vtkCellArray>::New();
	
	for (int i = 0; i < max_i-1; i++) {
		for (int j = 0; j < max_j-1; j++) {
			if (cell.at(n).at(i).at(j).type != 18) {
				vtkIdType IDs[5];
				bool allCorners = false;
				
				//~ if (c_0.r_2 != 0 && c_0.x_2 != 0) {
					//~ strips->InsertNextCell(5);
				//~ } else {
					//~ strips->InsertNextCell(4);
				//~ }
				
				gasCell c_0 = cell.at(n).at(i).at(j);
				gasCell c_i1 = cell.at(n).at(i+1).at(j).type != 18 ? 
					cell.at(n).at(i+1).at(j) : cell.at(n).at(i).at(j);
				gasCell c_i_1 = cell.at(n).at(i-1).at(j).type != 18 ? 
					cell.at(n).at(i-1).at(j) : cell.at(n).at(i).at(j);
				gasCell c_j1 = cell.at(n).at(i).at(j+1).type != 18 ? 
					cell.at(n).at(i).at(j+1) : cell.at(n).at(i).at(j);
				gasCell c_j_1 = cell.at(n).at(i).at(j-1).type != 18 ? 
					cell.at(n).at(i).at(j-1) : cell.at(n).at(i).at(j);
				gasCell c_i1j1 = cell.at(n).at(i+1).at(j+1).type != 18 ?
					cell.at(n).at(i+1).at(j+1) : c_i1;
				gasCell c_i1j_1 = cell.at(n).at(i+1).at(j-1).type != 18 ?
					cell.at(n).at(i+1).at(j-1) : c_i1;
				gasCell c_i_1j1 = cell.at(n).at(i-1).at(j+1).type != 18 ?
					cell.at(n).at(i-1).at(j+1) : c_i_1;
				gasCell c_i_1j_1 = cell.at(n).at(i-1).at(j-1).type != 18 ?
					cell.at(n).at(i-1).at(j-1) : c_i_1;
				
				IDs[0] = points->InsertNextPoint (i*dx, j*dr, 0 );
				outputP->InsertNextValue ( (c_0.P[0]+c_i_1j_1.P[0]+c_i_1.P[0]+c_j_1.P[0])/4 );
				outputE->InsertNextValue ( (c_0.e+c_i_1j_1.e+c_i_1.e+c_j_1.e)/4 );
				outputRho->InsertNextValue ( (c_0.rho+c_i_1j_1.rho+c_i_1.rho+c_j_1.rho)/4 );
				outputZ->InsertNextValue ( (c_0.final_z+c_i_1j_1.final_z+c_i_1.final_z+c_j_1.final_z)/4 );
				outputV->InsertNextTuple3 ( (c_0.Vx[0]+c_i_1j_1.Vx[0]+c_i_1.Vx[0]+c_j_1.Vx[0])/4, 
					(c_0.Vr[0]+c_i_1j_1.Vr[0]+c_i_1.Vr[0]+c_j_1.Vr[0])/4, 0 );
				//~ strips->InsertCellPoint(IDs[0]);
				
				IDs[1] = points->InsertNextPoint ((i+c_0.x_1)*dx, j*dr, 0 );
				outputP->InsertNextValue ( (c_0.P[0]+c_i1j_1.P[0]+c_i1.P[0]+c_j_1.P[0])/4 );
				outputE->InsertNextValue ( (c_0.e+c_i1j_1.e+c_i1.e+c_j_1.e)/4 );
				outputRho->InsertNextValue ( (c_0.rho+c_i1j_1.rho+c_i1.rho+c_j_1.rho)/4 );
				outputZ->InsertNextValue ( (c_0.final_z+c_i1j_1.final_z+c_i1.final_z+c_j_1.final_z)/4 );
				outputV->InsertNextTuple3 ( (c_0.Vx[0]+c_i1j_1.Vx[0]+c_i1.Vx[0]+c_j_1.Vx[0])/4, 
					(c_0.Vr[0]+c_i1j_1.Vr[0]+c_i1.Vr[0]+c_j_1.Vr[0])/4, 0 );
				//~ strips->InsertCellPoint(IDs[1]);
				
				IDs[2] = points->InsertNextPoint (i*dx, (j+c_0.r_1)*dr, 0 );
				outputP->InsertNextValue ( (c_0.P[0]+c_i_1j1.P[0]+c_i_1.P[0]+c_j1.P[0])/4 );
				outputE->InsertNextValue ( (c_0.e+c_i_1j1.e+c_i_1.e+c_j1.e)/4 );
				outputRho->InsertNextValue ( (c_0.rho+c_i_1j1.rho+c_i_1.rho+c_j1.rho)/4 );
				outputZ->InsertNextValue ( (c_0.final_z+c_i_1j1.final_z+c_i_1.final_z+c_j1.final_z)/4 );
				outputV->InsertNextTuple3 ( (c_0.Vx[0]+c_i_1j1.Vx[0]+c_i_1.Vx[0]+c_j1.Vx[0])/4, 
					(c_0.Vr[0]+c_i_1j1.Vr[0]+c_i_1.Vr[0]+c_j1.Vr[0])/4, 0 );
				//~ strips->InsertCellPoint(IDs[2]);
				
				if (c_0.r_2 != 0 && c_0.x_2 != 0) {
					IDs[3] = points->InsertNextPoint ((i+c_0.x_2)*dx, (j+c_0.r_2)*dr, 0 );
					outputP->InsertNextValue ( (c_0.P[0]+c_i1j1.P[0]+c_i1.P[0]+c_j1.P[0])/4 );
					outputE->InsertNextValue ( (c_0.e+c_i1j1.e+c_i1.e+c_j1.e)/4 );
					outputRho->InsertNextValue ( (c_0.rho+c_i1j1.rho+c_i1.rho+c_j1.rho)/4 );
					outputZ->InsertNextValue ( (c_0.final_z+c_i1j1.final_z+c_i1.final_z+c_j1.final_z)/4 );
					outputV->InsertNextTuple3 ( (c_0.Vx[0]+c_i1j1.Vx[0]+c_i1.Vx[0]+c_j1.Vx[0])/4, 
						(c_0.Vr[0]+c_i1j1.Vr[0]+c_i1.Vr[0]+c_j1.Vr[0])/4, 0 );
					//~ strips->InsertCellPoint(IDs[3]);
					allCorners = true;
				}
				
				IDs[4] = points->InsertNextPoint ((i+(c_0.x_1+c_0.x_2)/4)*dx, 
					(j+(c_0.r_1+c_0.r_2)/4)*dr, 0 );
				outputP->InsertNextValue ( c_0.P[0] );
				outputE->InsertNextValue ( c_0.e );
				outputRho->InsertNextValue ( c_0.rho );
				outputZ->InsertNextValue ( c_0.final_z );
				outputV->InsertNextTuple3 ( c_0.Vx[0], c_0.Vr[0], 0 );
				//~ strips->InsertCellPoint(IDs[5]);
				
				strips->InsertNextCell(3);
				strips->InsertCellPoint(IDs[0]);
				strips->InsertCellPoint(IDs[1]);
				strips->InsertCellPoint(IDs[4]);
				strips->InsertNextCell(3);
				strips->InsertCellPoint(IDs[0]);
				strips->InsertCellPoint(IDs[2]);
				strips->InsertCellPoint(IDs[4]);
				if (allCorners) {
					strips->InsertNextCell(3);
					strips->InsertCellPoint(IDs[1]);
					strips->InsertCellPoint(IDs[3]);
					strips->InsertCellPoint(IDs[4]);
					strips->InsertNextCell(3);
					strips->InsertCellPoint(IDs[2]);
					strips->InsertCellPoint(IDs[3]);
					strips->InsertCellPoint(IDs[4]);
				} else {
					strips->InsertNextCell(3);
					strips->InsertCellPoint(IDs[1]);
					strips->InsertCellPoint(IDs[2]);
					strips->InsertCellPoint(IDs[4]);
				}
			}
		}
	}
	
	// Create a polydata object and add the points to it.
	polydata->SetPoints(points);
	polydata->SetStrips(strips);
	polydata->GetPointData()->AddArray(outputP);
	polydata->GetPointData()->AddArray(outputE);
	polydata->GetPointData()->AddArray(outputRho);
	polydata->GetPointData()->SetVectors(outputV);
	polydata->GetPointData()->AddArray(outputZ);
	 
	// Write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
		writer->SetInput(polydata);
	#else
		writer->SetInputData(polydata);
	#endif
	
	// Optional - set the mode. The default is binary.
	writer->SetDataModeToBinary();
	//writer->SetDataModeToAscii();
	writer->Write();
}

void outputCSV(cell2d cell, std::ofstream & outputGas) {
	for (i = 0; i < max_i; i++) {
		for (j = 0; j < max_j; j++) {
			outputGas << setiosflags(ios::fixed) << setprecision(10)
				<< t.at(n) << "," 
				<< i << "," 
				<< j << "," 
				<< cell.at(n+1).at(i).at(j).P[0] << "," 
				<< cell.at(n+1).at(i).at(j).rho << "," 
				<< cell.at(n+1).at(i).at(j).e << ","
				<< cell.at(n+1).at(i).at(j).Vx[0] << ","
				<< cell.at(n+1).at(i).at(j).Vr[0] << ","
				<< cell.at(n).at(i).at(j).bar_Vx[0] << ","
				<< cell.at(n).at(i).at(j).bar_Vr[0] << ","
				<< cell.at(n).at(i).at(j).bar_e << ","
				<< cell.at(n).at(i).at(j).m << ","
				<< cell.at(n+1).at(i).at(j).z << ","
				<< cell.at(n+1).at(i).at(j).psi << ","
				<< cell.at(n).at(i).at(j).dM[1] << ","
				<< cell.at(n).at(i).at(j).dM[2] << ","
				<< cell.at(n).at(i).at(j).dM[3] << ","
				<< cell.at(n).at(i).at(j).dM[4] << ","
				<< cell.at(n).at(i).at(j).A[0] << ","
				<< cell.at(n).at(i).at(j).A[1] << ","
				<< cell.at(n).at(i).at(j).A[2] << ","
				<< cell.at(n).at(i).at(j).A[3] << ","
				<< cell.at(n).at(i).at(j).A[4] << ","
				<< cell.at(n+1).at(i).at(j).e - (pow(cell.at(n+1).at(i).at(j).Vx[0],2)+pow(cell.at(n+1).at(i).at(j).Vr[0],2))/2 << endl;
		}
	}
}

void prepOutputDynCSV(std::ofstream & outputDyn) {
	outputDyn << "t,i_sn,x_sn,U_sn,P_0,P_sn" << endl;
}

void outputDynCSV(std::ofstream & outputDyn, double t, int i_sn, double x_sn, double U_sn, double boltP, double projP) {
	outputDyn << setiosflags(ios::fixed) << setprecision(10) 
		<< t << ","
		<< i_sn << ","
		<< x_sn << ","
		<< U_sn << ","
		<< boltP << ","
		<< projP << endl;
}
