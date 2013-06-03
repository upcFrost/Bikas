/*
 * unused.cpp
 *
 *  Created on: Jun 3, 2013
 *      Author: frost
 *
 *      Here the unused code is stored. Sometimes it becomes useful,
 *      so i don't want to delete it.
 */

	/** Smoothing **/
	/* X axis */
//	int pointNum = 4; // Should be even, number of interpolating points
//	for (j = 0; j < max_j; j++) {
//		for (i = i_sn-10; i < i_sn-3; i++) {
//			if (cell.at(n+1).at(i).at(j).type != 18 &&
//					cell.at(n+1).at(i+1).at(j).type != 18 &&
//					cell.at(n+1).at(i+2).at(j).type != 18 &&
//					cell.at(n+1).at(i-1).at(j).type != 18 &&
//					cell.at(n+1).at(i-2).at(j).type != 18) {
//				double xPoints[pointNum];
//				double VxPoints[pointNum]; double VxPointsLin[pointNum]; double VxPointsQuad[pointNum]; double VxPointsCub[pointNum];
//				double VrPoints[pointNum]; double VrPointsLin[pointNum]; double VrPointsQuad[pointNum]; double VrPointsCub[pointNum];
//				double PPoints[pointNum]; double PPointsLin[pointNum]; double PPointsQuad[pointNum]; double PPointsCub[pointNum];
//				double ePoints[pointNum]; double ePointsLin[pointNum]; double ePointsQuad[pointNum]; double ePointsCub[pointNum];
//				double rhoPoints[pointNum]; double rhoPointsLin[pointNum]; double rhoPointsQuad[pointNum]; double rhoPointsCub[pointNum];
//				for (int iter = 0; iter < pointNum; iter++) {
//					int current_i = iter < pointNum/2 ? i - pointNum/2 + iter : i - pointNum/2 + iter + 1;
//					xPoints[iter] = current_i * dx;
//					VxPoints[iter] = cell.at(n+1).at(current_i).at(j).Vx[0];
//					VrPoints[iter] = cell.at(n+1).at(current_i).at(j).Vr[0];
//					PPoints[iter] = cell.at(n+1).at(current_i).at(j).P[0];
//					ePoints[iter] = cell.at(n+1).at(current_i).at(j).e;
//					rhoPoints[iter] = cell.at(n+1).at(current_i).at(j).rho;
//				}
//				cubic_nak ( pointNum, xPoints, VxPoints, VxPointsLin, VxPointsQuad, VxPointsCub );
//				cubic_nak ( pointNum, xPoints, VrPoints, VrPointsLin, VrPointsQuad, VrPointsCub );
//				cubic_nak ( pointNum, xPoints, PPoints, PPointsLin, PPointsQuad, PPointsCub );
//				cubic_nak ( pointNum, xPoints, ePoints, ePointsLin, ePointsQuad, ePointsCub );
//				cubic_nak ( pointNum, xPoints, rhoPoints, rhoPointsLin, rhoPointsQuad, rhoPointsCub );
//
//				cell.at(n+1).at(i).at(j).Vx[0] = spline_eval ( pointNum, xPoints, VxPoints, VxPointsLin, VxPointsQuad, VxPointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).Vr[0] = spline_eval ( pointNum, xPoints, VrPoints, VrPointsLin, VrPointsQuad, VrPointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).P[0] = spline_eval ( pointNum, xPoints, PPoints, PPointsLin, PPointsQuad, PPointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).e = spline_eval ( pointNum, xPoints, ePoints, ePointsLin, ePointsQuad, ePointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).rho = spline_eval ( pointNum, xPoints, rhoPoints, rhoPointsLin, rhoPointsQuad, rhoPointsCub, i*dx );
//			}
//		}
//	}
//	int pointNum = 4; // Should be even, number of interpolating points
//	for (j = 0; j < max_j; j++) {
//		for (i = 1+pointNum/2; i < i_sn-1-pointNum/2; i++) {
//			if (cell.at(n+1).at(i).at(j).type != 18 &&
//					cell.at(n+1).at(i+1).at(j).type != 18 &&
//					cell.at(n+1).at(i+2).at(j).type != 18 &&
//					cell.at(n+1).at(i-1).at(j).type != 18 &&
//					cell.at(n+1).at(i-2).at(j).type != 18) {
//				double xPoints[pointNum];
//				double VxPoints[pointNum]; double VxPointsLin[pointNum]; double VxPointsQuad[pointNum]; double VxPointsCub[pointNum];
//				double VrPoints[pointNum]; double VrPointsLin[pointNum]; double VrPointsQuad[pointNum]; double VrPointsCub[pointNum];
//				double PPoints[pointNum]; double PPointsLin[pointNum]; double PPointsQuad[pointNum]; double PPointsCub[pointNum];
//				double ePoints[pointNum]; double ePointsLin[pointNum]; double ePointsQuad[pointNum]; double ePointsCub[pointNum];
//				for (int iter = 0; iter < pointNum; iter++) {
//					int current_i = iter < pointNum/2 ? i - pointNum/2 + iter : i - pointNum/2 + iter + 1;
//					xPoints[iter] = current_i * dx;
//					VxPoints[iter] = cell.at(n+1).at(current_i).at(j).Vx[0];
//					VrPoints[iter] = cell.at(n+1).at(current_i).at(j).Vr[0];
//					PPoints[iter] = cell.at(n+1).at(current_i).at(j).P[0];
//					ePoints[iter] = cell.at(n+1).at(current_i).at(j).e;
//				}
//				cubic_nak ( pointNum, xPoints, VxPoints, VxPointsLin, VxPointsQuad, VxPointsCub );
//				cubic_nak ( pointNum, xPoints, VrPoints, VrPointsLin, VrPointsQuad, VrPointsCub );
//				cubic_nak ( pointNum, xPoints, PPoints, PPointsLin, PPointsQuad, PPointsCub );
//				cubic_nak ( pointNum, xPoints, ePoints, ePointsLin, ePointsQuad, ePointsCub );
//
//				cell.at(n+1).at(i).at(j).Vx[0] = spline_eval ( pointNum, xPoints, VxPoints, VxPointsLin, VxPointsQuad, VxPointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).Vr[0] = spline_eval ( pointNum, xPoints, VrPoints, VrPointsLin, VrPointsQuad, VrPointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).P[0] = spline_eval ( pointNum, xPoints, PPoints, PPointsLin, PPointsQuad, PPointsCub, i*dx );
//				cell.at(n+1).at(i).at(j).e = spline_eval ( pointNum, xPoints, ePoints, ePointsLin, ePointsQuad, ePointsCub, i*dx );
//			}
//		}
//	}
	/* Y axis */
	//~ for (i = 0; i < max_i; i++) {
		//~ for (j = 0; j < max_j; j++) {
			//~ switch (cell.at(n+1).at(i).at(j).type) {
				//~ // type 0 --> free cell
				//~ case 0:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
					//~
				//~ // type 15 --> top and left borders closed
				//~ case 15:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] - 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
				//~
				//~ // type 17 --> left border closed
				//~ case 17:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
				//~
				//~ // type 16 --> bottom and left borders closed
				//~ case 16:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
			//~
				//~ // type 13 --> top border closed
				//~ case 13:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] - 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
			//~
				//~ // type 14 --> bottom border closed
				//~ case 14:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
				//~
				//~ // type 19 --> right border closed
				//~ case 19:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
				//~
				//~ // type 21 --> bottom and right borders closed
				//~ case 21:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] + 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
					//~
				//~ // type 20 --> top and right borders closed
				//~ case 20:
					//~ cell.at(n+1).at(i).at(j).Vr[0] = 0.1*cell.at(n+1).at(i).at(j-1).Vr[0] + 0.8*cell.at(n+1).at(i).at(j).Vr[0] - 0.1*cell.at(n+1).at(i).at(j+1).Vr[0];
					//~ break;
					//~
				//~ default:
					//~ break;
			//~ }
		//~ }
	//~ }




