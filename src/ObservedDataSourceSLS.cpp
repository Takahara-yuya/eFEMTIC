//-------------------------------------------------------------------------------------------------------
// This file is part of eFEMTIC
//
// Copyright (C) 2025 Zuwei Huang
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.
//-------------------------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "ObservedDataSourceSLS.h"
#include "ObservedDataSource.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

//-------------------------------------------------------
//----- Definition of the class ObservedDataSourceSLS -----
//-------------------------------------------------------

// Constructer
ObservedDataSourceSLS::ObservedDataSourceSLS() :
	ObservedDataSource(),
	x_start(0.0),
	y_start(0.0),
	z_start(0.0),
	x_end(0.0),
	y_end(0.0),
	z_end(0.0),
	lengthAfterpartioning(0.0),
	m_sourceDirection(1.0),
	currentIntensity(1.0)
{
}

// Destructer
ObservedDataSourceSLS::~ObservedDataSourceSLS() {

}

// Read data from input file
void ObservedDataSourceSLS::inputObservedData(std::ifstream& inFile) {

	const double convKilometerToMeter = 1000.0;
	//
	dipoleParameter.numOfEquivalentDipoleInThisSource = 0;
	dipoleParameter.sourceTypeOfAllEquivalentDipoleSource = ELECTRIC_DIPOLE_SOURCE;
	int nED = maxNumOfEquivalentDipoles;
	//
	dipoleParameter.elementIncludingEquivalentDipoleInThisSource = new int[nED];
	dipoleParameter.sourceID = new int[nED];
	dipoleParameter.midPoint = new CommonParameters::locationXYZ[nED];
	dipoleParameter.parameterOfEquivalentDipoleInThisSource = new double* [nED];
	dipoleParameter.directionOfEquivalentDipoles = new double* [nED];
	for (int i = 0; i < nED; i++) {
		dipoleParameter.parameterOfEquivalentDipoleInThisSource[i] = new double[2];
		dipoleParameter.directionOfEquivalentDipoles[i] = new double[3];
	}
	//start reading file
	std::string tmpchar;
	inFile >> tmpchar;
	//maybe check here
	inFile >> m_sourceID;
	inFile >> m_sourceDirection >> currentIntensity;
	inFile >> x_start >> y_start >> z_start;
	inFile >> x_end >> y_end >> z_end;
	x_start *= convKilometerToMeter; y_start *= convKilometerToMeter; z_start *= convKilometerToMeter;
	x_end *= convKilometerToMeter; y_end *= convKilometerToMeter; z_end *= convKilometerToMeter;

	calculateDirecctionVector();

	//Partion the current line !!
	partionTheCurrentLine();

	OutputFiles::m_logFile << "# Partitioning of the long wire is completed." << std::endl;
	if (lengthAfterpartioning != 0) {
		lengthAfterpartioning = 0;
	}
	for (int i = 0; i < dipoleParameter.numOfEquivalentDipoleInThisSource; i++) {
		lengthAfterpartioning += dipoleParameter.parameterOfEquivalentDipoleInThisSource[i][1];
	}
}

void ObservedDataSourceSLS::calculateDirecctionVector() {
	//calculate direction
	double length(0);
	length = sqrt(pow(x_start - x_end, 2) + pow(y_start - y_end, 2) + pow(z_start - z_end, 2));
	if (length == 0) {
		OutputFiles::m_logFile << "Source length in 0 !!" << std::endl;
		exit(1);
	}
	directionX = (x_end - x_start) / length * m_sourceDirection;
	directionY = (y_end - y_start) / length * m_sourceDirection;
	directionZ = (z_end - z_start) / length * m_sourceDirection;
	//check if the direction vector unit
	double units = std::sqrt(directionX * directionX + directionY * directionY + directionZ * directionZ);
	if (units != 1) {
		OutputFiles::m_logFile << "Direction vector don't unit in source ID: " << dipoleParameter.sourceID << std::endl;
		exit(1);
	}
}

void ObservedDataSourceSLS::partionTheCurrentLine() {
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	double overlapLength, parallelLength, crossLength;
	Point3D overlapMidpoint, parallelMidpoint, crossMidpoint, TetraCrossStart, TetraCrossEnd;
	int nElem, nFace, faceStatus[4], iDipole;
	bool doNotFindInThisElement(false);
	bool isEquivalentDipoleInFace(false);
	bool isEquivalentDipoleInEdge(false);
	std::vector<int> IDsOfdoNotFindInThisElement;
	std::vector<Point3D> overlapEdgeNodeCoordinateFirst, overlapEdgeNodeCoordinateSecond;
	if (meshType == MeshData::TETRA) {
		const MeshDataTetraElement* const pMeshDataTetraElement = (AnalysisControl::getInstance())->getPointerOfMeshDataTetraElement();
		nElem = pMeshDataTetraElement->getNumElemTotal();
		Point3D segStart, segEnd;
		segStart.x = x_start; segEnd.x = x_end;
		segStart.y = y_start; segEnd.y = y_end;
		segStart.z = z_start; segEnd.z = z_end;
		for (int iElem = 0; iElem < nElem; iElem++) {
			//jump the source in the air !! < FOR NOW, MAY BE DELETED IN FUTURE? >
			if (log10(pResistivityBlock->getResistivityValuesFromElemID(iElem)) > 7.0) {
				continue;
			}
			//jump the share face source
			if (find(IDsOfdoNotFindInThisElement.begin(), IDsOfdoNotFindInThisElement.end(), iElem) != IDsOfdoNotFindInThisElement.end()) {
				continue;
			}
			//jump the share edge source
			Point3D tetA, tetB, tetC, tetD, intersection[4], result1[3], result2[3];
			int nodeID;
			nodeID = pMeshDataTetraElement->getNodesOfElements(iElem, 0);
			tetA.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
			tetA.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
			tetA.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
			nodeID = pMeshDataTetraElement->getNodesOfElements(iElem, 1);
			tetB.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
			tetB.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
			tetB.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
			nodeID = pMeshDataTetraElement->getNodesOfElements(iElem, 2);
			tetC.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
			tetC.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
			tetC.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
			nodeID = pMeshDataTetraElement->getNodesOfElements(iElem, 3);
			tetD.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
			tetD.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
			tetD.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
			if (isOverlapedEdgeHasBeenCalculatedInThisTetrahedron(tetA, tetB, tetC, tetD, overlapEdgeNodeCoordinateFirst, overlapEdgeNodeCoordinateSecond)) {
				continue;
			}
			//1.0 test if start and end located in this tetrahedron !!
			bool isSegStartInThisTetrahedron, isSegEndInThisTetrahedron;
			isSegStartInThisTetrahedron = isPointInTetrahedron(segStart, tetA, tetB, tetC, tetD, 1e-10);
			isSegEndInThisTetrahedron = isPointInTetrahedron(segEnd, tetA, tetB, tetC, tetD, 1e-10);
			//
			if (isSegStartInThisTetrahedron && isSegEndInThisTetrahedron) {//all the segment is in the tetrahedron
				dipoleParameter.numOfEquivalentDipoleInThisSource++;
				iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
				if (iDipole > maxNumOfEquivalentDipoles - 1) {
					OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
					exit(1);
				}
				crossLength = distance(segStart, segEnd);
				dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
				dipoleParameter.midPoint[iDipole].X = 0.5 * (segStart.x + segEnd.x);
				dipoleParameter.midPoint[iDipole].Y = 0.5 * (segStart.y + segEnd.y);
				dipoleParameter.midPoint[iDipole].Z = 0.5 * (segStart.z + segEnd.z);
				dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
				dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
				dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
				dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
				dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
				break;
			}
			else if (isSegStartInThisTetrahedron || isSegEndInThisTetrahedron) {//one of the node is in the Tetrahedron
				/*OutputFiles::m_logFile << tetA.x << " " << tetA.y << " " << tetA.z << std::endl;
				OutputFiles::m_logFile << tetB.x << " " << tetB.y << " " << tetB.z << std::endl;
				OutputFiles::m_logFile << tetC.x << " " << tetC.y << " " << tetC.z << std::endl;
				OutputFiles::m_logFile << tetD.x << " " << tetD.y << " " << tetD.z << std::endl;
				OutputFiles::m_logFile << segStart.x << " " << segStart.y << " " << segStart.z << std::endl;
				OutputFiles::m_logFile << segEnd.x << " " << segEnd.y << " " << segEnd.z << std::endl;*/
				//initialize the original face state
				for (int iFace = 0; iFace < 4; iFace++) {
					faceStatus[iFace] = 4;
					assignmentPoint3D(intersection[iFace], defaultValue, defaultValue, defaultValue);
				}
				//
				Point3D crossStart[4], crossEnd[4];
				isEquivalentDipoleInEdge = false;
				isEquivalentDipoleInFace = false;
				for (int iFace = 0; iFace < 4; iFace++) {
					if (isEquivalentDipoleInEdge || isEquivalentDipoleInFace) {
						isEquivalentDipoleInEdge = false;
						isEquivalentDipoleInFace = false;
						break;
					}
					nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 0);
					tetA.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
					tetA.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
					tetA.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
					nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 1);
					tetB.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
					tetB.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
					tetB.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
					nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 2);
					tetC.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
					tetC.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
					tetC.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
					//
					/*OutputFiles::m_logFile << iFace << std::endl;
					OutputFiles::m_logFile << tetA.x << " " << tetA.y << " " << tetA.z << std::endl;
					OutputFiles::m_logFile << tetB.x << " " << tetB.y << " " << tetB.z << std::endl;
					OutputFiles::m_logFile << tetC.x << " " << tetC.y << " " << tetC.z << std::endl;
					OutputFiles::m_logFile << segStart.x << " " << segStart.y << " " << segStart.z << std::endl;
					OutputFiles::m_logFile << segEnd.x << " " << segEnd.y << " " << segEnd.z << std::endl;*/
					//1.0 chek if the segment parallel to the face
					if (isSegmentParallelToPlane(segStart, segEnd, tetA, tetB, tetC)) {// segment parallel to this face
						// check the status with the line and segment
						int status1 = checkSegmentsIntersection(tetA, tetB, segStart, segEnd, result1[0], result2[0]);
						int status2 = checkSegmentsIntersection(tetA, tetC, segStart, segEnd, result1[1], result2[1]);
						int status3 = checkSegmentsIntersection(tetB, tetC, segStart, segEnd, result1[2], result2[2]);
						//
						if (status1 == -1 && status2 == -1 && status3 == -1) {// no cross
							faceStatus[iFace] = 4;
						}
						else if (status1 == 1 || status2 == 1 || status3 == 1) {// overlap with the edge0
							faceStatus[iFace] = 1;
							isEquivalentDipoleInEdge = true;
							if (status1 == 1) {
								overlapEdgeNodeCoordinateFirst.push_back(tetA);
								overlapEdgeNodeCoordinateSecond.push_back(tetB);
								overlapLength = distance(result1[0], result2[0]);
								overlapMidpoint.x = 0.5 * (result1[0].x + result2[0].x);
								overlapMidpoint.y = 0.5 * (result1[0].y + result2[0].y);
								overlapMidpoint.z = 0.5 * (result1[0].z + result2[0].z);
							}
							else if (status2 == 1) {
								overlapEdgeNodeCoordinateFirst.push_back(tetA);
								overlapEdgeNodeCoordinateSecond.push_back(tetC);
								overlapLength = distance(result1[1], result2[1]);
								overlapMidpoint.x = 0.5 * (result1[1].x + result2[1].x);
								overlapMidpoint.y = 0.5 * (result1[1].y + result2[1].y);
								overlapMidpoint.z = 0.5 * (result1[1].z + result2[1].z);
							}
							else if (status3 == 1) {
								overlapEdgeNodeCoordinateFirst.push_back(tetB);
								overlapEdgeNodeCoordinateSecond.push_back(tetC);
								overlapLength = distance(result1[2], result2[2]);
								overlapMidpoint.x = 0.5 * (result1[2].x + result2[2].x);
								overlapMidpoint.y = 0.5 * (result1[2].y + result2[2].y);
								overlapMidpoint.z = 0.5 * (result1[2].z + result2[2].z);
							}
						}
						else {//layer on the face, calculate the intersection point with the segment
							faceStatus[iFace] = 2;
							isEquivalentDipoleInFace = true;
							//
							if (isSegStartInThisTetrahedron) {
								crossEnd[iFace] = segStart;
							}
							else {
								crossEnd[iFace] = segEnd;
							}
							if (checkIfSegmentCrossOnePointOfThisFace(result1, crossEnd[iFace])) {
								if (isSegStartInThisTetrahedron) {
									crossStart[iFace] = segStart;
								}
								else {
									crossStart[iFace] = segEnd;
								}
								parallelLength = distance(crossStart[iFace], crossEnd[iFace]);
								parallelMidpoint.x = 0.5 * (crossStart[iFace].x + crossEnd[iFace].x);
								parallelMidpoint.y = 0.5 * (crossStart[iFace].y + crossEnd[iFace].y);
								parallelMidpoint.z = 0.5 * (crossStart[iFace].z + crossEnd[iFace].z);
							}
							else {
								OutputFiles::m_logFile << "Insection points do not equal to 1 in one face" << std::endl;
								exit(1);
							}
						}
					}
					else {//check if segment cross this face 
						bool isSegmentCrossThisFace;
						isSegmentCrossThisFace = segmentTriangleIntersection(segStart, segEnd, tetA, tetB, tetC, intersection[iFace]);
						if (isSegmentCrossThisFace) {
							faceStatus[iFace] = 3;
						}
						else {
							faceStatus[iFace] = 4;
						}
					}
					//4 status in one face:
					//1. overlap with th edge.--
					//2. layer on the face.--
					//3. cross the face in one point (intersection).
					//4. do not cross the face.
				}
				//already loop in the four faces!!
				int numOfInsectionFace(0);
				for (int iFace = 0; iFace < 4; iFace++) {
					if (faceStatus[iFace] == 1) {
						if (overlapLength > 0) {
							//only need to find eqdipole in this face !!
							doNotFindInThisElement = true;
							dipoleParameter.numOfEquivalentDipoleInThisSource++;
							iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
							if (iDipole > maxNumOfEquivalentDipoles - 1) {
								OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
								exit(1);
							}
							dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
							dipoleParameter.midPoint[iDipole].X = overlapMidpoint.x;
							dipoleParameter.midPoint[iDipole].Y = overlapMidpoint.y;
							dipoleParameter.midPoint[iDipole].Z = overlapMidpoint.z;
							dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
							dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = overlapLength;
							dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
							dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
							dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
						}
						break;
					}
					else if (faceStatus[iFace] == 2) {
						if (parallelLength > 0) {
							//only need to find eqdipole in this face !!
							doNotFindInThisElement = true;
							int neiNum = 1;//this face!!
							for (int i = 0; i < neiNum; i++) {
								IDsOfdoNotFindInThisElement.push_back(pMeshDataTetraElement->getIDOfNeighborElement(iElem, iFace));
								//	OutputFiles::m_logFile << pMeshDataTetraElement->getIDOfNeighborElement(iElem, iFace) << std::endl;
							}
							dipoleParameter.numOfEquivalentDipoleInThisSource++;
							iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
							if (iDipole > maxNumOfEquivalentDipoles - 1) {
								OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
								exit(1);
							}
							dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
							dipoleParameter.midPoint[iDipole].X = parallelMidpoint.x;
							dipoleParameter.midPoint[iDipole].Y = parallelMidpoint.y;
							dipoleParameter.midPoint[iDipole].Z = parallelMidpoint.z;
							dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
							dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = parallelLength;
							dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
							dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
							dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
						}
						//OutputFiles::m_logFile << parallelMidpoint.x << " " << parallelMidpoint.y << " " << parallelMidpoint.z << " " << parallelLength << std::endl;
						break;
					}
					else if (faceStatus[iFace] == 3) {
						//find two insection point!!
						numOfInsectionFace++;
					}
					else if (faceStatus[iFace] == 4) {
						continue;
					}
					else {
						OutputFiles::m_logFile << "Face status error in partioning the dipoles !!" << std::endl;
						exit(1);
					}
				}
				//cross the tetrahedron without parallel to any face and have at least one intersection point !!
				if (numOfInsectionFace > 0 && !doNotFindInThisElement) {
					if (isSegStartInThisTetrahedron) {
						TetraCrossStart = segStart;
					}
					else {
						TetraCrossStart = segEnd;
					}
					if (checkIfSegmentCrossOnePointOfThisTetrahedron(intersection, TetraCrossStart, TetraCrossEnd)) {
						dipoleParameter.numOfEquivalentDipoleInThisSource++;
						iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
						if (iDipole > maxNumOfEquivalentDipoles - 1) {
							OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
							exit(1);
						}
						crossLength = distance(TetraCrossStart, TetraCrossEnd);
						dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
						dipoleParameter.midPoint[iDipole].X = 0.5 * (TetraCrossStart.x + TetraCrossEnd.x);
						dipoleParameter.midPoint[iDipole].Y = 0.5 * (TetraCrossStart.y + TetraCrossEnd.y);
						dipoleParameter.midPoint[iDipole].Z = 0.5 * (TetraCrossStart.z + TetraCrossEnd.z);
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
						doNotFindInThisElement = false;
					}
				}
				//
			}
			else { //segment cross this tetrahedron
				//initialize the original face state
				for (int iFace = 0; iFace < 4; iFace++) {
					faceStatus[iFace] = 4;
					assignmentPoint3D(intersection[iFace], defaultValue, defaultValue, defaultValue);
				}
				//start loop
				Point3D crossStart[4], crossEnd[4];
				isEquivalentDipoleInEdge = false;
				isEquivalentDipoleInFace = false;
				for (int iFace = 0; iFace < 4; iFace++) {
					if (isEquivalentDipoleInEdge || isEquivalentDipoleInFace) {
						isEquivalentDipoleInEdge = false;
						isEquivalentDipoleInFace = false;
						break;
					}
					nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 0);
					tetA.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
					tetA.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
					tetA.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
					nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 1);
					tetB.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
					tetB.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
					tetB.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
					nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 2);
					tetC.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
					tetC.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
					tetC.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
					//1.0 chek if the segment parallel to the face
					if (isSegmentParallelToPlane(segStart, segEnd, tetA, tetB, tetC)) {// segment parallel to this face
						// break;
						// check the status with the line and segment
						int status1 = checkSegmentsIntersection(tetA, tetB, segStart, segEnd, result1[0], result2[0]);
						int status2 = checkSegmentsIntersection(tetA, tetC, segStart, segEnd, result1[1], result2[1]);
						int status3 = checkSegmentsIntersection(tetB, tetC, segStart, segEnd, result1[2], result2[2]);
						//
						if (status1 == -1 && status2 == -1 && status3 == -1) {// no cross 
							faceStatus[iFace] = 4;
						}
						else if (status1 == 1 || status2 == 1 || status3 == 1) {// overlap with the edge0
							isEquivalentDipoleInEdge = true;
							faceStatus[iFace] = 1;
							if (status1 == 1) {
								overlapEdgeNodeCoordinateFirst.push_back(tetA);
								overlapEdgeNodeCoordinateSecond.push_back(tetB);
								overlapLength = distance(result1[0], result2[0]);
								overlapMidpoint.x = 0.5 * (result1[0].x + result2[0].x);
								overlapMidpoint.y = 0.5 * (result1[0].y + result2[0].y);
								overlapMidpoint.z = 0.5 * (result1[0].z + result2[0].z);
							}
							else if (status2 == 1) {
								overlapEdgeNodeCoordinateFirst.push_back(tetA);
								overlapEdgeNodeCoordinateSecond.push_back(tetC);
								overlapLength = distance(result1[1], result2[1]);
								overlapMidpoint.x = 0.5 * (result1[1].x + result2[1].x);
								overlapMidpoint.y = 0.5 * (result1[1].y + result2[1].y);
								overlapMidpoint.z = 0.5 * (result1[1].z + result2[1].z);
							}
							else if (status3 == 1) {
								overlapEdgeNodeCoordinateFirst.push_back(tetB);
								overlapEdgeNodeCoordinateSecond.push_back(tetC);
								overlapLength = distance(result1[2], result2[2]);
								overlapMidpoint.x = 0.5 * (result1[2].x + result2[2].x);
								overlapMidpoint.y = 0.5 * (result1[2].y + result2[2].y);
								overlapMidpoint.z = 0.5 * (result1[2].z + result2[2].z);
							}
						}
						else {//layer on the face, calculate the intersection point with the segment
							isEquivalentDipoleInFace = true;
							faceStatus[iFace] = 2;
							//
							if (checkIfSegmentCrossTwoPointOfThisFace(result1, crossStart[iFace], crossEnd[iFace])) {
								parallelLength = distance(crossStart[iFace], crossEnd[iFace]);
								parallelMidpoint.x = 0.5 * (crossStart[iFace].x + crossEnd[iFace].x);
								parallelMidpoint.y = 0.5 * (crossStart[iFace].y + crossEnd[iFace].y);
								parallelMidpoint.z = 0.5 * (crossStart[iFace].z + crossEnd[iFace].z);
							}
							else {
								OutputFiles::m_logFile << "Insection points do not equal to 2 in one face" << std::endl;
								exit(1);
							}
						}
					}
					else {//check if segment cross this face 
						bool isSegmentCrossThisFace;
						isSegmentCrossThisFace = segmentTriangleIntersection(segStart, segEnd, tetA, tetB, tetC, intersection[iFace]);
						if (isSegmentCrossThisFace) {
							faceStatus[iFace] = 3;
						}
						else {
							faceStatus[iFace] = 4;
						}
					}
					//4 status in one face:
					//1. overlap with th edge.--
					//2. layer on the face.--
					//3. cross the face in one point (intersection).
					//4. do not cross the face.
				}
				//already loop in the four faces!!
				int numOfInsectionFace(0);
				for (int iFace = 0; iFace < 4; iFace++) {
					//OutputFiles::m_logFile << faceStatus[iFace] << std::endl;
					if (faceStatus[iFace] == 1) {
						//only need to find eqdipole in this face !!
						doNotFindInThisElement = true;
						dipoleParameter.numOfEquivalentDipoleInThisSource++;
						iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
						if (iDipole > maxNumOfEquivalentDipoles - 1) {
							OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
							exit(1);
						}
						dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
						dipoleParameter.midPoint[iDipole].X = overlapMidpoint.x;
						dipoleParameter.midPoint[iDipole].Y = overlapMidpoint.y;
						dipoleParameter.midPoint[iDipole].Z = overlapMidpoint.z;
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = overlapLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
						break;
					}
					else if (faceStatus[iFace] == 2) {
						//only need to find eqdipole in this face !!
						doNotFindInThisElement = true;
						int neiNum = pMeshDataTetraElement->getNumNeighborElement();
						for (int i = 0; i < neiNum; i++) {
							IDsOfdoNotFindInThisElement.push_back(pMeshDataTetraElement->getIDOfNeighborElement(iElem, iFace));
						}
						dipoleParameter.numOfEquivalentDipoleInThisSource++;
						iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
						if (iDipole > maxNumOfEquivalentDipoles - 1) {
							OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
							exit(1);
						}
						//
						dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
						dipoleParameter.midPoint[iDipole].X = parallelMidpoint.x;
						dipoleParameter.midPoint[iDipole].Y = parallelMidpoint.y;
						dipoleParameter.midPoint[iDipole].Z = parallelMidpoint.z;
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = parallelLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
						break;
					}
					else if (faceStatus[iFace] == 3) {
						//find two insection point!!
						nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 0);
						tetA.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
						tetA.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
						tetA.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
						nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 1);
						tetB.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
						tetB.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
						tetB.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
						nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, 2);
						tetC.x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
						tetC.y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
						tetC.z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
						numOfInsectionFace++;
					}
					else if (faceStatus[iFace] == 4) {
						continue;
					}
					else {
						OutputFiles::m_logFile << "Face status error in partioning the dipoles !!" << std::endl;
						exit(1);
					}
				}
				//cross the tetrahedron without parallel to any face and have two intersection points !!
				if (numOfInsectionFace > 1 && !doNotFindInThisElement) {
					if (checkIfSegmentCrossTwoPointOfThisTetrahedron(intersection, TetraCrossStart, TetraCrossEnd)) {
						dipoleParameter.numOfEquivalentDipoleInThisSource++;
						iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
						if (iDipole > maxNumOfEquivalentDipoles - 1) {
							OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
							exit(1);
						}
						crossLength = distance(TetraCrossStart, TetraCrossEnd);
						dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
						dipoleParameter.midPoint[iDipole].X = 0.5 * (TetraCrossStart.x + TetraCrossEnd.x);
						dipoleParameter.midPoint[iDipole].Y = 0.5 * (TetraCrossStart.y + TetraCrossEnd.y);
						dipoleParameter.midPoint[iDipole].Z = 0.5 * (TetraCrossStart.z + TetraCrossEnd.z);
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity;
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ;
						doNotFindInThisElement = false;
					}
				}
			}
		}
	}
}

double ObservedDataSourceSLS::getTotalLengthOfEquivalentDipoles() {

	return lengthAfterpartioning;

}