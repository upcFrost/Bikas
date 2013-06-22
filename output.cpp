#include <stdio.h>
#include "output.h"

OutPVD::OutPVD(cell2dStatic & cell) {
	// Points array
	points = vtkSmartPointer<vtkPoints>::New();
	// General polydata array
	polydata = vtkSmartPointer<vtkPolyData>::New();
	// Data arrays
	outputP = vtkSmartPointer<vtkDoubleArray>::New();
	outputE = vtkSmartPointer<vtkDoubleArray>::New();
	outputZ = vtkSmartPointer<vtkDoubleArray>::New();
	outputRho = vtkSmartPointer<vtkDoubleArray>::New();
	outputV = vtkSmartPointer<vtkDoubleArray>::New();
	// Vertices array
	verts = vtkSmartPointer<vtkCellArray>::New();
	// Cells array
	strips = vtkSmartPointer<vtkCellArray>::New();

	outputP->SetName("P");
	outputE->SetName("E");
	outputRho->SetName("Rho");
	outputV->SetName("V");
	outputV->SetNumberOfComponents(3);
	outputV->SetComponentName(0, "Vx");
	outputV->SetComponentName(1, "Vr");
	outputV->SetComponentName(2, "Vz");
	outputZ->SetName("Z");

	OutPVD::initArrays(cell);
}

void OutPVD::initArrays(cell2dStatic & cell) {
	for (int i = 0; i < max_i-1; i++) {
		for (int j = 0; j < max_j-1; j++) {
			if (cell[i][j].type != 18) {
				linkID IDs;
				bool allCorners = false;

				gasCell c_0 = cell[i][j];
				gasCell c_i1 = cell[i+1][j].type != 18 ? cell[i+1][j] : cell[i][j];
				gasCell c_i_1 = cell[i-1][j].type != 18 ? cell[i-1][j] : cell[i][j];
				gasCell c_j1 = cell[i][j+1].type != 18 ? cell[i][j+1] : cell[i][j];
				gasCell c_j_1 = cell[i][j-1].type != 18 ? cell[i][j-1] : cell[i][j];
				gasCell c_i1j1 = cell[i+1][j+1].type != 18 ? cell[i+1][j+1] : c_i1;
				gasCell c_i1j_1 = cell[i+1][j-1].type != 18 ? cell[i+1][j-1] : c_i1;
				gasCell c_i_1j1 = cell[i-1][j+1].type != 18 ? cell[i-1][j+1] : c_i_1;
				gasCell c_i_1j_1 = cell[i-1][j-1].type != 18 ? cell[i-1][j-1] : c_i_1;

				// Detect IDs...
				IDs.i = i;
				IDs.j = j;
				IDs.id[0] = points->InsertNextPoint (i*dx*scaleD, j*dr*scaleD, 0 );
				IDs.id[1] = points->InsertNextPoint ((i+c_0.x_1)*dx*scaleD, j*dr*scaleD, 0 );
				IDs.id[2] = points->InsertNextPoint (i*dx*scaleD, (j+c_0.r_1)*dr*scaleD, 0 );
				if (c_0.r_2 != 0 && c_0.x_2 != 0) {
					IDs.id[3] = points->InsertNextPoint ((i+c_0.x_2)*dx*scaleD, (j+c_0.r_2)*dr*scaleD, 0 );
					allCorners = true;
				} else {
					IDs.id[3] = -1;
				}
				IDs.id[4] = points->InsertNextPoint ((i+(c_0.x_1+c_0.x_2)/4)*dx*scaleD,
					(j+(c_0.r_1+c_0.r_2)/4)*dr*scaleD, 0 );
				// ...and assign them to the variable
				OutPVD::idVector.push_back(IDs);

				// Insert first values into arrays to resize them
				// TODO: Fix it. Should use Resize function
				outputP->InsertNextValue ( (c_0.P+c_i_1j_1.P+c_i_1.P+c_j_1.P)/4 * scaleP);
				outputE->InsertNextValue ( (c_0.e+c_i_1j_1.e+c_i_1.e+c_j_1.e)/4 * scaleE );
				outputRho->InsertNextValue ( (c_0.rho+c_i_1j_1.rho+c_i_1.rho+c_j_1.rho)/4 * scaleRho);
				outputZ->InsertNextValue ( (c_0.final_z+c_i_1j_1.final_z+c_i_1.final_z+c_j_1.final_z)/4 );
				outputV->InsertNextTuple3 ( (c_0.Vx+c_i_1j_1.Vx+c_i_1.Vx+c_j_1.Vx)/4  * scaleV,
					(c_0.Vr+c_i_1j_1.Vr+c_i_1.Vr+c_j_1.Vr)/4 * scaleV, 0 );

				outputP->InsertNextValue ( (c_0.P+c_i1j_1.P+c_i1.P+c_j_1.P)/4 * scaleP );
				outputE->InsertNextValue ( (c_0.e+c_i1j_1.e+c_i1.e+c_j_1.e)/4 * scaleE );
				outputRho->InsertNextValue ( (c_0.rho+c_i1j_1.rho+c_i1.rho+c_j_1.rho)/4 * scaleRho );
				outputZ->InsertNextValue ( (c_0.final_z+c_i1j_1.final_z+c_i1.final_z+c_j_1.final_z)/4 );
				outputV->InsertNextTuple3 ( (c_0.Vx+c_i1j_1.Vx+c_i1.Vx+c_j_1.Vx)/4 * scaleV,
					(c_0.Vr+c_i1j_1.Vr+c_i1.Vr+c_j_1.Vr)/4 * scaleV, 0 );
				//~ strips->InsertCellPoint(IDs[1]);

				outputP->InsertNextValue ( (c_0.P+c_i_1j1.P+c_i_1.P+c_j1.P)/4 * scaleP );
				outputE->InsertNextValue ( (c_0.e+c_i_1j1.e+c_i_1.e+c_j1.e)/4 * scaleE );
				outputRho->InsertNextValue ( (c_0.rho+c_i_1j1.rho+c_i_1.rho+c_j1.rho)/4 * scaleRho );
				outputZ->InsertNextValue ( (c_0.final_z+c_i_1j1.final_z+c_i_1.final_z+c_j1.final_z)/4 );
				outputV->InsertNextTuple3 ( (c_0.Vx+c_i_1j1.Vx+c_i_1.Vx+c_j1.Vx)/4 * scaleV,
					(c_0.Vr+c_i_1j1.Vr+c_i_1.Vr+c_j1.Vr)/4 * scaleV, 0 );
				//~ strips->InsertCellPoint(IDs[2]);

				if (c_0.r_2 != 0 && c_0.x_2 != 0) {
					outputP->InsertNextValue ( (c_0.P+c_i1j1.P+c_i1.P+c_j1.P)/4 * scaleP );
					outputE->InsertNextValue ( (c_0.e+c_i1j1.e+c_i1.e+c_j1.e)/4 * scaleE );
					outputRho->InsertNextValue ( (c_0.rho+c_i1j1.rho+c_i1.rho+c_j1.rho)/4 * scaleRho );
					outputZ->InsertNextValue ( (c_0.final_z+c_i1j1.final_z+c_i1.final_z+c_j1.final_z)/4 );
					outputV->InsertNextTuple3 ( (c_0.Vx+c_i1j1.Vx+c_i1.Vx+c_j1.Vx)/4 * scaleV,
						(c_0.Vr+c_i1j1.Vr+c_i1.Vr+c_j1.Vr)/4 * scaleV, 0 );
					//~ strips->InsertCellPoint(IDs[3]);
					allCorners = true;
				}

				outputP->InsertNextValue ( c_0.P * scaleP );
				outputE->InsertNextValue ( c_0.e * scaleE );
				outputRho->InsertNextValue ( c_0.rho * scaleRho );
				outputZ->InsertNextValue ( c_0.final_z );
				outputV->InsertNextTuple3 ( c_0.Vx * scaleV, c_0.Vr * scaleV, 0 );
				//~ strips->InsertCellPoint(IDs[5]);

				// And finally - populate cell vector
				strips->InsertNextCell(3);
				strips->InsertCellPoint(IDs.id[0]);
				strips->InsertCellPoint(IDs.id[1]);
				strips->InsertCellPoint(IDs.id[4]);
				strips->InsertNextCell(3);
				strips->InsertCellPoint(IDs.id[0]);
				strips->InsertCellPoint(IDs.id[2]);
				strips->InsertCellPoint(IDs.id[4]);
				if (allCorners) {
					strips->InsertNextCell(3);
					strips->InsertCellPoint(IDs.id[1]);
					strips->InsertCellPoint(IDs.id[3]);
					strips->InsertCellPoint(IDs.id[4]);
					strips->InsertNextCell(3);
					strips->InsertCellPoint(IDs.id[2]);
					strips->InsertCellPoint(IDs.id[3]);
					strips->InsertCellPoint(IDs.id[4]);
				} else {
					strips->InsertNextCell(3);
					strips->InsertCellPoint(IDs.id[1]);
					strips->InsertCellPoint(IDs.id[2]);
					strips->InsertCellPoint(IDs.id[4]);
				}
			}
		}
	}

	// Now we'll point polydata array to data arrays
	polydata->SetPoints(points);
	polydata->SetStrips(strips);
	polydata->GetPointData()->AddArray(outputP);
	polydata->GetPointData()->AddArray(outputE);
	polydata->GetPointData()->AddArray(outputRho);
	polydata->GetPointData()->SetVectors(outputV);
	polydata->GetPointData()->AddArray(outputZ);
}

void OutPVD::MakeOutput(cell2d & cell, std::string & filename) {
#pragma omp parallel for num_threads(4) schedule(dynamic)
	for (unsigned int idx = 0; idx < idVector.size(); idx++) {
		int i = idVector[idx].i;
		int j = idVector[idx].j;

		gasCell & c_0 = cell[n][i][j];
		gasCell & c_i1 = cell[n][i+1][j].type != 18 ? cell[n][i+1][j] : cell[n][i][j];
		gasCell & c_i_1 = cell[n][i-1][j].type != 18 ? cell[n][i-1][j] : cell[n][i][j];
		gasCell & c_j1 = cell[n][i][j+1].type != 18 ? cell[n][i][j+1] : cell[n][i][j];
		gasCell & c_j_1 = cell[n][i][j-1].type != 18 ? cell[n][i][j-1] : cell[n][i][j];
		gasCell & c_i1j1 = cell[n][i+1][j+1].type != 18 ? cell[n][i+1][j+1] : c_i1;
		gasCell & c_i1j_1 = cell[n][i+1][j-1].type != 18 ? cell[n][i+1][j-1] : c_i1;
		gasCell & c_i_1j1 = cell[n][i-1][j+1].type != 18 ? cell[n][i-1][j+1] : c_i_1;
		gasCell & c_i_1j_1 = cell[n][i-1][j-1].type != 18 ? cell[n][i-1][j-1] : c_i_1;

		outputP->SetValue(idVector[idx].id[0], (c_0.P+c_i_1j_1.P+c_i_1.P+c_j_1.P)/4 * scaleP);
		outputE->SetValue(idVector[idx].id[0],  (c_0.e+c_i_1j_1.e+c_i_1.e+c_j_1.e)/4 * scaleE );
		outputRho->SetValue(idVector[idx].id[0],  (c_0.rho+c_i_1j_1.rho+c_i_1.rho+c_j_1.rho)/4 * scaleRho);
		outputZ->SetValue(idVector[idx].id[0],  (c_0.final_z+c_i_1j_1.final_z+c_i_1.final_z+c_j_1.final_z)/4 );
		outputV->SetTuple3(idVector[idx].id[0],  (c_0.Vx+c_i_1j_1.Vx+c_i_1.Vx+c_j_1.Vx)/4  * scaleV,
			(c_0.Vr+c_i_1j_1.Vr+c_i_1.Vr+c_j_1.Vr)/4 * scaleV, 0 );

		outputP->SetValue (idVector[idx].id[1], (c_0.P+c_i1j_1.P+c_i1.P+c_j_1.P)/4 * scaleP );
		outputE->SetValue (idVector[idx].id[1], (c_0.e+c_i1j_1.e+c_i1.e+c_j_1.e)/4 * scaleE );
		outputRho->SetValue (idVector[idx].id[1], (c_0.rho+c_i1j_1.rho+c_i1.rho+c_j_1.rho)/4 * scaleRho );
		outputZ->SetValue (idVector[idx].id[1], (c_0.final_z+c_i1j_1.final_z+c_i1.final_z+c_j_1.final_z)/4 );
		outputV->SetTuple3 (idVector[idx].id[1], (c_0.Vx+c_i1j_1.Vx+c_i1.Vx+c_j_1.Vx)/4 * scaleV,
			(c_0.Vr+c_i1j_1.Vr+c_i1.Vr+c_j_1.Vr)/4 * scaleV, 0 );

		outputP->SetValue (idVector[idx].id[2],  (c_0.P+c_i_1j1.P+c_i_1.P+c_j1.P)/4 * scaleP );
		outputE->SetValue (idVector[idx].id[2],  (c_0.e+c_i_1j1.e+c_i_1.e+c_j1.e)/4 * scaleE );
		outputRho->SetValue (idVector[idx].id[2],  (c_0.rho+c_i_1j1.rho+c_i_1.rho+c_j1.rho)/4 * scaleRho );
		outputZ->SetValue (idVector[idx].id[2],  (c_0.final_z+c_i_1j1.final_z+c_i_1.final_z+c_j1.final_z)/4 );
		outputV->SetTuple3 (idVector[idx].id[2],  (c_0.Vx+c_i_1j1.Vx+c_i_1.Vx+c_j1.Vx)/4 * scaleV,
			(c_0.Vr+c_i_1j1.Vr+c_i_1.Vr+c_j1.Vr)/4 * scaleV, 0 );

		if (idVector[idx].id[3] != -1) {
			outputP->SetValue (idVector[idx].id[3], (c_0.P+c_i1j1.P+c_i1.P+c_j1.P)/4 * scaleP );
			outputE->SetValue (idVector[idx].id[3], (c_0.e+c_i1j1.e+c_i1.e+c_j1.e)/4 * scaleE );
			outputRho->SetValue (idVector[idx].id[3], (c_0.rho+c_i1j1.rho+c_i1.rho+c_j1.rho)/4 * scaleRho );
			outputZ->SetValue (idVector[idx].id[3], (c_0.final_z+c_i1j1.final_z+c_i1.final_z+c_j1.final_z)/4 );
			outputV->SetTuple3 (idVector[idx].id[3], (c_0.Vx+c_i1j1.Vx+c_i1.Vx+c_j1.Vx)/4 * scaleV,
				(c_0.Vr+c_i1j1.Vr+c_i1.Vr+c_j1.Vr)/4 * scaleV, 0 );
		}

		outputP->SetValue (idVector[idx].id[4], c_0.P * scaleP );
		outputE->SetValue (idVector[idx].id[4], c_0.e * scaleE );
		outputRho->SetValue (idVector[idx].id[4], c_0.rho * scaleRho );
		outputZ->SetValue (idVector[idx].id[4], c_0.final_z );
		outputV->SetTuple3 (idVector[idx].id[4],
				c_0.Vx * scaleV, c_0.Vr * scaleV, 0 );
	}

	// Write the file
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetFileName(filename.c_str());
	#if VTK_MAJOR_VERSION <= 5
		writer->SetInput(polydata);
	#else
		writer->SetInputData(polydata);
	#endif

	// Set the mode
	writer->SetDataModeToBinary();
	writer->Write();
}

//void OutputPVD(cell2d & cell, std::string & filename) {
//	// Points array
//	vtkSmartPointer<vtkPoints> points =
//	    vtkSmartPointer<vtkPoints>::New();
//	// General polydata array
//	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
//	// Data arrays
//	vtkSmartPointer<vtkDoubleArray> outputP = vtkSmartPointer<vtkDoubleArray>::New();
//	outputP->SetName("P");
//	vtkSmartPointer<vtkDoubleArray> outputE = vtkSmartPointer<vtkDoubleArray>::New();
//	outputE->SetName("E");
//	vtkSmartPointer<vtkDoubleArray> outputZ = vtkSmartPointer<vtkDoubleArray>::New();
//	outputZ->SetName("Z");
//	vtkSmartPointer<vtkDoubleArray> outputRho = vtkSmartPointer<vtkDoubleArray>::New();
//	outputRho->SetName("Rho");
//	vtkSmartPointer<vtkDoubleArray> outputV = vtkSmartPointer<vtkDoubleArray>::New();
//	outputV->SetName("V");
//	outputV->SetNumberOfComponents(3);
//	outputV->SetComponentName(0, "Vx");
//	outputV->SetComponentName(1, "Vr");
//	outputV->SetComponentName(2, "Vz");
//	// Vertices array
//	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
//	// Cells array
//	vtkSmartPointer<vtkCellArray> strips = vtkSmartPointer<vtkCellArray>::New();
//
//	for (int i = 0; i < max_i-1; i++) {
//		for (int j = 0; j < max_j-1; j++) {
//			if (cell.at(n).at(i).at(j).type != 18) {
//				vtkIdType IDs[5] = {0};
//				bool allCorners = false;
//
//				//~ if (c_0.r_2 != 0 && c_0.x_2 != 0) {
//					//~ strips->InsertNextCell(5);
//				//~ } else {
//					//~ strips->InsertNextCell(4);
//				//~ }
//				gasCell c_0 = cell[n][i][j];
//				gasCell c_i1 = cell[n][i+1][j].type != 18 ?
//					cell[n][i+1][j] : cell[n][i][j];
//				gasCell c_i_1 = cell[n][i-1][j].type != 18 ?
//					cell[n][i-1][j] : cell[n][i][j];
//				gasCell c_j1 = cell[n][i][j+1].type != 18 ?
//					cell[n][i][j+1] : cell[n][i][j];
//				gasCell c_j_1 = cell[n][i][j-1].type != 18 ?
//					cell[n][i][j-1] : cell[n][i][j];
//				gasCell c_i1j1 = cell[n][i+1][j+1].type != 18 ?
//					cell[n][i+1][j+1] : c_i1;
//				gasCell c_i1j_1 = cell[n][i+1][j-1].type != 18 ?
//					cell[n][i+1][j-1] : c_i1;
//				gasCell c_i_1j1 = cell[n][i-1][j+1].type != 18 ?
//					cell[n][i-1][j+1] : c_i_1;
//				gasCell c_i_1j_1 = cell[n][i-1][j-1].type != 18 ?
//					cell[n][i-1][j-1] : c_i_1;
//
//				IDs[0] = points->InsertNextPoint (i*dx*scaleD, j*dr*scaleD, 0 );
//				outputP->InsertNextValue ( (c_0.P+c_i_1j_1.P+c_i_1.P+c_j_1.P)/4 * scaleP);
//				outputE->InsertNextValue ( (c_0.e+c_i_1j_1.e+c_i_1.e+c_j_1.e)/4 * scaleE );
//				outputRho->InsertNextValue ( (c_0.rho+c_i_1j_1.rho+c_i_1.rho+c_j_1.rho)/4 * scaleRho);
//				outputZ->InsertNextValue ( (c_0.final_z+c_i_1j_1.final_z+c_i_1.final_z+c_j_1.final_z)/4 );
//				outputV->InsertNextTuple3 ( (c_0.Vx+c_i_1j_1.Vx+c_i_1.Vx+c_j_1.Vx)/4  * scaleV,
//					(c_0.Vr+c_i_1j_1.Vr+c_i_1.Vr+c_j_1.Vr)/4 * scaleV, 0 );
//				//~ strips->InsertCellPoint(IDs[0]);
//
//				IDs[1] = points->InsertNextPoint ((i+c_0.x_1)*dx*scaleD, j*dr*scaleD, 0 );
//				outputP->InsertNextValue ( (c_0.P+c_i1j_1.P+c_i1.P+c_j_1.P)/4 * scaleP );
//				outputE->InsertNextValue ( (c_0.e+c_i1j_1.e+c_i1.e+c_j_1.e)/4 * scaleE );
//				outputRho->InsertNextValue ( (c_0.rho+c_i1j_1.rho+c_i1.rho+c_j_1.rho)/4 * scaleRho );
//				outputZ->InsertNextValue ( (c_0.final_z+c_i1j_1.final_z+c_i1.final_z+c_j_1.final_z)/4 );
//				outputV->InsertNextTuple3 ( (c_0.Vx+c_i1j_1.Vx+c_i1.Vx+c_j_1.Vx)/4 * scaleV,
//					(c_0.Vr+c_i1j_1.Vr+c_i1.Vr+c_j_1.Vr)/4 * scaleV, 0 );
//				//~ strips->InsertCellPoint(IDs[1]);
//
//				IDs[2] = points->InsertNextPoint (i*dx*scaleD, (j+c_0.r_1)*dr*scaleD, 0 );
//				outputP->InsertNextValue ( (c_0.P+c_i_1j1.P+c_i_1.P+c_j1.P)/4 * scaleP );
//				outputE->InsertNextValue ( (c_0.e+c_i_1j1.e+c_i_1.e+c_j1.e)/4 * scaleE );
//				outputRho->InsertNextValue ( (c_0.rho+c_i_1j1.rho+c_i_1.rho+c_j1.rho)/4 * scaleRho );
//				outputZ->InsertNextValue ( (c_0.final_z+c_i_1j1.final_z+c_i_1.final_z+c_j1.final_z)/4 );
//				outputV->InsertNextTuple3 ( (c_0.Vx+c_i_1j1.Vx+c_i_1.Vx+c_j1.Vx)/4 * scaleV,
//					(c_0.Vr+c_i_1j1.Vr+c_i_1.Vr+c_j1.Vr)/4 * scaleV, 0 );
//				//~ strips->InsertCellPoint(IDs[2]);
//
//				if (c_0.r_2 != 0 && c_0.x_2 != 0) {
//					IDs[3] = points->InsertNextPoint ((i+c_0.x_2)*dx*scaleD, (j+c_0.r_2)*dr*scaleD, 0 );
//					outputP->InsertNextValue ( (c_0.P+c_i1j1.P+c_i1.P+c_j1.P)/4 * scaleP );
//					outputE->InsertNextValue ( (c_0.e+c_i1j1.e+c_i1.e+c_j1.e)/4 * scaleE );
//					outputRho->InsertNextValue ( (c_0.rho+c_i1j1.rho+c_i1.rho+c_j1.rho)/4 * scaleRho );
//					outputZ->InsertNextValue ( (c_0.final_z+c_i1j1.final_z+c_i1.final_z+c_j1.final_z)/4 );
//					outputV->InsertNextTuple3 ( (c_0.Vx+c_i1j1.Vx+c_i1.Vx+c_j1.Vx)/4 * scaleV,
//						(c_0.Vr+c_i1j1.Vr+c_i1.Vr+c_j1.Vr)/4 * scaleV, 0 );
//					//~ strips->InsertCellPoint(IDs[3]);
//					allCorners = true;
//				}
//
//				IDs[4] = points->InsertNextPoint ((i+(c_0.x_1+c_0.x_2)/4)*dx*scaleD,
//					(j+(c_0.r_1+c_0.r_2)/4)*dr*scaleD, 0 );
//				outputP->InsertNextValue ( c_0.P * scaleP );
//				outputE->InsertNextValue ( c_0.e * scaleE );
//				outputRho->InsertNextValue ( c_0.rho * scaleRho );
//				outputZ->InsertNextValue ( c_0.final_z );
//				outputV->InsertNextTuple3 ( c_0.Vx * scaleV, c_0.Vr * scaleV, 0 );
//				//~ strips->InsertCellPoint(IDs[5]);
//
//				strips->InsertNextCell(3);
//				strips->InsertCellPoint(IDs[0]);
//				strips->InsertCellPoint(IDs[1]);
//				strips->InsertCellPoint(IDs[4]);
//				strips->InsertNextCell(3);
//				strips->InsertCellPoint(IDs[0]);
//				strips->InsertCellPoint(IDs[2]);
//				strips->InsertCellPoint(IDs[4]);
//				if (allCorners) {
//					strips->InsertNextCell(3);
//					strips->InsertCellPoint(IDs[1]);
//					strips->InsertCellPoint(IDs[3]);
//					strips->InsertCellPoint(IDs[4]);
//					strips->InsertNextCell(3);
//					strips->InsertCellPoint(IDs[2]);
//					strips->InsertCellPoint(IDs[3]);
//					strips->InsertCellPoint(IDs[4]);
//				} else {
//					strips->InsertNextCell(3);
//					strips->InsertCellPoint(IDs[1]);
//					strips->InsertCellPoint(IDs[2]);
//					strips->InsertCellPoint(IDs[4]);
//				}
//			}
//		}
//	}
//
//	// Create a polydata object and add the points to it.
//	polydata->SetPoints(points);
//	polydata->SetStrips(strips);
//	polydata->GetPointData()->AddArray(outputP);
//	polydata->GetPointData()->AddArray(outputE);
//	polydata->GetPointData()->AddArray(outputRho);
//	polydata->GetPointData()->SetVectors(outputV);
//	polydata->GetPointData()->AddArray(outputZ);
//
//	// Write the file
//	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//	writer->SetFileName(filename.c_str());
//	#if VTK_MAJOR_VERSION <= 5
//		writer->SetInput(polydata);
//	#else
//		writer->SetInputData(polydata);
//	#endif
//
//	// Optional - set the mode. The default is binary.
//	writer->SetDataModeToBinary();
////	writer->SetDataModeToAscii();
//	writer->Write();
//}

void outputCSV(cell2d & cell, std::ofstream & outputGas) {
	for (int i = 0; i < max_i; i++) {
		for (int j = 0; j < max_j; j++) {
			outputGas << std::setiosflags(std::ios::fixed) << std::setprecision(10)
				<< t.at(n) << "," 
				<< i << "," 
				<< j << "," 
				<< cell.at(nextN).at(i).at(j).P << ","
				<< cell.at(nextN).at(i).at(j).rho << ","
				<< cell.at(nextN).at(i).at(j).e << ","
				<< cell.at(nextN).at(i).at(j).Vx << ","
				<< cell.at(nextN).at(i).at(j).Vr << ","
				<< cell.at(n).at(i).at(j).bar_Vx << ","
				<< cell.at(n).at(i).at(j).bar_Vr << ","
				<< cell.at(n).at(i).at(j).bar_e << ","
				<< cell.at(nextN).at(i).at(j).final_z << ","
				<< cell.at(nextN).at(i).at(j).final_psi << ","
				<< cell.at(n).at(i).at(j).dM[1] << ","
				<< cell.at(n).at(i).at(j).dM[2] << ","
				<< cell.at(n).at(i).at(j).dM[3] << ","
				<< cell.at(n).at(i).at(j).dM[4] << ","
				<< cell.at(n).at(i).at(j).A[0] << ","
				<< cell.at(n).at(i).at(j).A[1] << ","
				<< cell.at(n).at(i).at(j).A[2] << ","
				<< cell.at(n).at(i).at(j).A[3] << ","
				<< cell.at(n).at(i).at(j).A[4] << ","
				<< cell.at(nextN).at(i).at(j).e - (pow(cell.at(nextN).at(i).at(j).Vx,2)+pow(cell.at(nextN).at(i).at(j).Vr,2))/2 << std::endl;
		}
	}
}

void prepOutputDynCSV(std::ofstream & outputDyn) {
	outputDyn << "t,i_sn,x_sn,U_sn,P_0,P_sn" << std::endl;
}

void prepOutputGasCSV(std::ofstream & outputGas, bool verbose) {
	if (verbose) {
		outputGas << "t,i,j,P,rho,e,Vx,Vr,bar_Vx,bar_Vr,bar_e,m,z,psi,dM[1],dM[2],dM[3],dM[4],A[0],A[1],A[2],A[3],A[4],IntE" << std::endl;
	} else {
		outputGas << "t,x,y,z,P,rho,e,Vx,Vr,z,psi,IntE" << std::endl;
	}
}

void outputDynCSV(std::ofstream & outputDyn, double t, int i_sn, double x_sn, double U_sn, double boltP, double projP) {
	outputDyn << std::setiosflags(std::ios::fixed) << std::setprecision(10)
		<< t * scaleT << ","
		<< i_sn << ","
		<< x_sn * scaleD << ","
		<< U_sn * scaleV << ","
		<< boltP * scaleP << ","
		<< projP * scaleP << std::endl;
}
