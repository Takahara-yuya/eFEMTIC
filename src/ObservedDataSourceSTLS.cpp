//-------------------------------------------------------------------------------------------------------
// This module provides CSEM (Controlled Source Electromagnetics) modeling capabilities within FEMTIC.
// Implemented by Zuwei Huang, 2025.
//
// For support, bug reports, or questions, please contact the author at: hzw1498218560@tongji.edu.cn
//-------------------------------------------------------------------------------------------------------
#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
#include <iomanip>

#include "ObservedDataSourceSTLS.h"
#include "ObservedDataSource.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

//-------------------------------------------------------
//----- Definition of the class ObservedDataSourceSTLS -----
//-------------------------------------------------------

// Constructer
ObservedDataSourceSTLS::ObservedDataSourceSTLS() :
	ObservedDataSource(),
	x_start(NULL),
	y_start(NULL),
	z_start(NULL),
	x_end(NULL),
	y_end(NULL),
	z_end(NULL),
	lengthAfterPartitioning(0.0),
	m_sourceDirection(NULL),
	currentIntensity(NULL),
	directionX(NULL),
	directionY(NULL),
	directionZ(NULL),
	m_findCurrentLineOnEarthSurface(false)
{
}

// Destructer
ObservedDataSourceSTLS::~ObservedDataSourceSTLS() {

	if (x_start != NULL) {
		delete[] x_start;
		x_start = NULL;
	}

	if (y_start != NULL) {
		delete[] y_start;
		y_start = NULL;
	}

	if (z_start != NULL) {
		delete[] z_start;
		z_start = NULL;
	}

	if (x_end != NULL) {
		delete[] x_end;
		x_end = NULL;
	}

	if (y_end != NULL) {
		delete[] y_end;
		y_end = NULL;
	}

	if (z_end != NULL) {
		delete[] z_end;
		z_end = NULL;
	}

	if (m_sourceDirection != NULL) {
		delete[] m_sourceDirection;
		m_sourceDirection = NULL;
	}

	if (currentIntensity != NULL) {
		delete[] currentIntensity;
		currentIntensity = NULL;
	}

	if (directionX != NULL) {
		delete[] directionX;
		directionX = NULL;
	}

	if (directionY != NULL) {
		delete[] directionY;
		directionY = NULL;
	}

	if (directionZ != NULL) {
		delete[] directionZ;
		directionZ = NULL;
	}

}

// Read data from input file
void ObservedDataSourceSTLS::inputObservedData(std::ifstream& inFile) {

	m_findCurrentLineOnEarthSurface = (AnalysisControl::getInstance())->getIsSourcePartitionedOnSurface();
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

	//allocate memories
	inFile >> m_sourceID;
	inFile >> m_numOfSegments;
	m_sourceDirection = new int[m_numOfSegments];
	currentIntensity = new double[m_numOfSegments];
	directionX = new double[m_numOfSegments];
	directionY = new double[m_numOfSegments];
	directionZ = new double[m_numOfSegments];
	x_start = new double[m_numOfSegments];
	y_start = new double[m_numOfSegments];
	z_start = new double[m_numOfSegments];
	x_end = new double[m_numOfSegments];
	y_end = new double[m_numOfSegments];
	z_end = new double[m_numOfSegments];
	//start partition
	for (int iSeg = 0; iSeg < m_numOfSegments; iSeg++) {

		inFile >> m_sourceDirection[iSeg] >> currentIntensity[iSeg];
		inFile >> x_start[iSeg] >> y_start[iSeg] >> z_start[iSeg];
		inFile >> x_end[iSeg] >> y_end[iSeg] >> z_end[iSeg];

		x_start[iSeg] *= convKilometerToMeter; y_start[iSeg] *= convKilometerToMeter; z_start[iSeg] *= convKilometerToMeter;
		x_end[iSeg] *= convKilometerToMeter; y_end[iSeg] *= convKilometerToMeter; z_end[iSeg] *= convKilometerToMeter;

		if (!m_findCurrentLineOnEarthSurface)calculateDirecctionVector(iSeg);

		//Partition the current line !!
		partitionTheCurrentLine(iSeg);

	}

	OutputFiles::m_logFile << "# Partitioning of the long wire is completed." << std::endl;
	if (lengthAfterPartitioning != 0) {
		lengthAfterPartitioning = 0;
	}

	//OutputFiles::m_logFile << "--------------------" << std::endl;
	for (int i = 0; i < dipoleParameter.numOfEquivalentDipoleInThisSource; i++) {
		lengthAfterPartitioning += dipoleParameter.parameterOfEquivalentDipoleInThisSource[i][1];
	}

}

void ObservedDataSourceSTLS::calculateDirecctionVector(int iSeg) {
	//calculate direction
	double length(0);
	length = sqrt(pow(x_start[iSeg] - x_end[iSeg], 2) + pow(y_start[iSeg] - y_end[iSeg], 2) + pow(z_start[iSeg] - z_end[iSeg], 2));
	if (length == 0) {
		OutputFiles::m_logFile << "Source length is 0 !!" << std::endl;
		OutputFiles::m_logFile << x_start[iSeg] << " " << y_start[iSeg] << " " << z_start[iSeg] << std::endl;
		OutputFiles::m_logFile << x_end[iSeg] << " " << y_end[iSeg] << " " << z_end[iSeg] << std::endl;
		exit(1);
	}
	directionX[iSeg] = (x_end[iSeg] - x_start[iSeg]) / length * m_sourceDirection[iSeg];
	directionY[iSeg] = (y_end[iSeg] - y_start[iSeg]) / length * m_sourceDirection[iSeg];
	directionZ[iSeg] = (z_end[iSeg] - z_start[iSeg]) / length * m_sourceDirection[iSeg];
	//check if the direction vector unit
	double units = std::sqrt(directionX[iSeg] * directionX[iSeg] + directionY[iSeg] * directionY[iSeg] + directionZ[iSeg] * directionZ[iSeg]);
	if (abs(units - 1) > 0.001) {
		OutputFiles::m_logFile << "Direction vector don't unit in source ID: " << m_sourceID << std::endl;
		exit(1);
	}
}

void ObservedDataSourceSTLS::partitionTheCurrentLine(int iSeg) {
	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	//
	if (meshType == MeshData::TETRA) {
		if (m_findCurrentLineOnEarthSurface) {
			findEquivalentDipolesInTetraMeshOnEarthSurface(iSeg);
		}
		else {
			findEquivalentDipolesInTetraMeshConsideringZCoordiante(iSeg);
		}
	}
	else if (meshType == MeshData::NONCONFORMING_HEXA) {
		findEquivalentDipolesInDHexaMeshOnEarthSurface(iSeg);
	}
}

double ObservedDataSourceSTLS::getTotalLengthOfEquivalentDipoles() {

	return lengthAfterPartitioning;

}

void ObservedDataSourceSTLS::findEquivalentDipolesInTetraMeshConsideringZCoordiante(int iSeg) {
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	double overlapLength, parallelLength, crossLength;
	Point3D overlapMidpoint, parallelMidpoint, crossMidpoint, TetraCrossStart, TetraCrossEnd;
	int nElem, nFace, faceStatus[4], iDipole;
	bool doNotFindInThisElement(false);
	bool isEquivalentDipoleInFace(false);
	bool isEquivalentDipoleInEdge(false);
	std::vector<int> IDsOfdoNotFindInThisElement;
	std::vector<Point3D> overlapEdgeNodeCoordinateFirst, overlapEdgeNodeCoordinateSecond;
	const MeshDataTetraElement* const pMeshDataTetraElement = (AnalysisControl::getInstance())->getPointerOfMeshDataTetraElement();
	nElem = pMeshDataTetraElement->getNumElemTotal();
	Point3D segStart, segEnd;
	segStart.x = x_start[iSeg]; segEnd.x = x_end[iSeg];
	segStart.y = y_start[iSeg]; segEnd.y = y_end[iSeg];
	segStart.z = z_start[iSeg]; segEnd.z = z_end[iSeg];
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
			dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
			dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
			dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
			dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
			dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
			break;
		}
		else if (isSegStartInThisTetrahedron || isSegEndInThisTetrahedron) {//one of the node is in the Tetrahedron
			//OutputFiles::m_logFile << iElem << " hahaha" << std::endl;
			//OutputFiles::m_logFile << tetA.x << " " << tetA.y << " " << tetA.z << std::endl;
			//OutputFiles::m_logFile << tetB.x << " " << tetB.y << " " << tetB.z << std::endl;
			//OutputFiles::m_logFile << tetC.x << " " << tetC.y << " " << tetC.z << std::endl;
			//OutputFiles::m_logFile << tetD.x << " " << tetD.y << " " << tetD.z << std::endl;
			//OutputFiles::m_logFile << segStart.x << " " << segStart.y << " " << segStart.z << std::endl;
			//OutputFiles::m_logFile << segEnd.x << " " << segEnd.y << " " << segEnd.z << std::endl;
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
				for (int i = 0; i < 3; i++) {
					result1[i] = defaultValue;
					result2[i] = defaultValue;
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
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = overlapLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
					}
					break;
				}
				else if (faceStatus[iFace] == 2) {
					//only need to find eqdipole in this face !!
					if (parallelLength > 0) {
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
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = parallelLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
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
					OutputFiles::m_logFile << "Face status error in Partitioning the dipoles !!" << std::endl;
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
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
					doNotFindInThisElement = false;
				}
			}
			//
		}
		else { //segment cross this tetrahedron
			/*	OutputFiles::m_logFile << tetA.x << " " << tetA.y << " " << tetA.z << std::endl;
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
			//start loop
			Point3D crossStart[4], crossEnd[4];
			isEquivalentDipoleInEdge = false;
			isEquivalentDipoleInFace = false;
			for (int iFace = 0; iFace < 4; iFace++) {
				for (int i = 0; i < 3; i++) {
					result1[i] = defaultValue;
					result2[i] = defaultValue;
				}
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
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = overlapLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
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
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = parallelLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
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
					OutputFiles::m_logFile << "Face status error in Partitioning the dipoles !!" << std::endl;
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
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
					doNotFindInThisElement = false;
				}
			}
		}
	}
}

void ObservedDataSourceSTLS::findEquivalentDipolesInDHexaMeshOnEarthSurface(int iSeg) {
	int upperNode2Edge[4][2];
	upperNode2Edge[0][0] = 0;
	upperNode2Edge[0][1] = 1;
	upperNode2Edge[1][0] = 1;
	upperNode2Edge[1][1] = 2;
	upperNode2Edge[2][0] = 2;
	upperNode2Edge[2][1] = 3;
	upperNode2Edge[3][0] = 3;
	upperNode2Edge[3][1] = 0;
	int shearEdgeFace[4];
	shearEdgeFace[0] = 2;
	shearEdgeFace[1] = 1;
	shearEdgeFace[2] = 3;
	shearEdgeFace[3] = 0;
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	double crossLength;
	int nElem, iDipole, edgeStatus[4];
	bool isEquivalentDipoleInEdge(false);
	bool isThisElemOnEarthSurface(true);
	std::vector<int> IDsOfdoNotFindInThisElement;
	std::vector<Point3D> p(4);
	const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
	nElem = ptrMeshDataNonConformingHexaElement->getNumElemTotal();
	Point3D segStart, segEnd;
	segStart.x = x_start[iSeg]; segEnd.x = x_end[iSeg];
	segStart.y = y_start[iSeg]; segEnd.y = y_end[iSeg];
	segStart.z = z_start[iSeg]; segEnd.z = z_end[iSeg];
	for (int iElem = 0; iElem < nElem; iElem++) {
		//	for (int iElem = 170889; iElem < 170891; iElem++) {
				//jump the source in the air !! < FOR NOW, MAY BE DELETED IN FUTURE? >
		if (log10(pResistivityBlock->getResistivityValuesFromElemID(iElem)) > 7.0) {
			continue;
		}
		//jump the share face source
		if (find(IDsOfdoNotFindInThisElement.begin(), IDsOfdoNotFindInThisElement.end(), iElem) != IDsOfdoNotFindInThisElement.end()) {
			continue;
		}
		//
		int boundaryTest = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, 4);//upper face
		if (boundaryTest < 1) {
			continue;//do not find on XY-mins boundary
		}
		for (int iNeigh = 0; iNeigh < boundaryTest; iNeigh++) {
			int upperElem = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, 4, iNeigh);
			if (log10(pResistivityBlock->getResistivityValuesFromElemID(upperElem)) < 7.0) {
				isThisElemOnEarthSurface = false;
				break;
			}
		}
		//
		if (!isThisElemOnEarthSurface) {
			isThisElemOnEarthSurface = true;
			continue;
		}
		//if point in this Hexa?
		std::vector<quadrilateralFace> faces(6);
		int nodeID;
		for (int iFace = 0; iFace < 6; iFace++) {
			for (int iNode = 0; iNode < 4; iNode++) {
				nodeID = ptrMeshDataNonConformingHexaElement->getNodeIDGlobalFromElementAndFace(iElem, iFace, iNode);
				faces[iFace].p[iNode].x = ptrMeshDataNonConformingHexaElement->getXCoordinatesOfNodes(nodeID);
				faces[iFace].p[iNode].y = ptrMeshDataNonConformingHexaElement->getYCoordinatesOfNodes(nodeID);
				faces[iFace].p[iNode].z = ptrMeshDataNonConformingHexaElement->getZCoordinatesOfNodes(nodeID);
			}
		}
		/*OutputFiles::m_logFile << iElem << " " << faces[4].p[0].x << " " << faces[4].p[0].y << " " << faces[4].p[0].z << " " << std::endl;
		OutputFiles::m_logFile << iElem << " " << faces[4].p[1].x << " " << faces[4].p[1].y << " " << faces[4].p[1].z << " " << std::endl;
		OutputFiles::m_logFile << iElem << " " << faces[4].p[2].x << " " << faces[4].p[2].y << " " << faces[4].p[2].z << " " << std::endl;
		OutputFiles::m_logFile << iElem << " " << faces[4].p[3].x << " " << faces[4].p[3].y << " " << faces[4].p[3].z << " " << std::endl;*/
		//
		bool isSegStartInThisTetrahedron, isSegEndInThisTetrahedron;
		isSegStartInThisTetrahedron = isPointOnThisHexaXYPlanes(faces[4].p, segStart, 1e-10);
		isSegEndInThisTetrahedron = isPointOnThisHexaXYPlanes(faces[4].p, segEnd, 1e-10);
		//OutputFiles::m_logFile << isSegStartInThisTetrahedron << " " << isSegEndInThisTetrahedron << std::endl;
		//starting find cross point on surface!
		if (isSegStartInThisTetrahedron && isSegEndInThisTetrahedron) {
			dipoleParameter.numOfEquivalentDipoleInThisSource++;
			iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
			if (iDipole > maxNumOfEquivalentDipoles - 1) {
				OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
				exit(1);
			}
			segStart.z = getZOnQuadPlane(faces[4].p[0], faces[4].p[1], faces[4].p[2], faces[4].p[3], segStart.x, segStart.y);
			segEnd.z = getZOnQuadPlane(faces[4].p[0], faces[4].p[1], faces[4].p[2], faces[4].p[3], segEnd.x, segEnd.y);
			crossLength = distance(segStart, segEnd);
			dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
			dipoleParameter.midPoint[iDipole].X = 0.5 * (segStart.x + segEnd.x);
			dipoleParameter.midPoint[iDipole].Y = 0.5 * (segStart.y + segEnd.y);
			dipoleParameter.midPoint[iDipole].Z = 0.5 * (segStart.z + segEnd.z);
			dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
			dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
			dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
			dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
			dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
			/*OutputFiles::m_logFile << iElem << " " << faces[4].p[0].x << " " << faces[4].p[0].y <<" "<< faces[4].p[0].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;
			OutputFiles::m_logFile << iElem << " " << faces[4].p[1].x << " " << faces[4].p[1].y <<" "<< faces[4].p[1].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;
			OutputFiles::m_logFile << iElem << " " << faces[4].p[2].x << " " << faces[4].p[2].y <<" "<< faces[4].p[2].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;
			OutputFiles::m_logFile << iElem << " " << faces[4].p[3].x << " " << faces[4].p[3].y <<" "<< faces[4].p[3].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;*/
			break;
		}
		else if (isSegStartInThisTetrahedron || isSegEndInThisTetrahedron) {
			for (int iEdge = 0; iEdge < 4; iEdge++) {
				edgeStatus[iEdge] = -1;
			}
			Point3D crossStart[4], crossEnd[4];
			isEquivalentDipoleInEdge = false;
			for (int iEdge = 0; iEdge < 4; iEdge++) {
				edgeStatus[iEdge] = checkSegmentIntersectionInXYPlane(segStart, segEnd, faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossStart[iEdge], crossEnd[iEdge]);
			}
			for (int iEdge = 0; iEdge < 4; iEdge++) {
				crossLength = distance(crossStart[iEdge], crossEnd[iEdge]);
				//OutputFiles::m_logFile << edgeStatus[iEdge] << " " << faces[4].p[upperNode2Edge[iEdge][0]].x << " " << faces[4].p[upperNode2Edge[iEdge][0]].y << " " << faces[4].p[upperNode2Edge[iEdge][1]].x << " " << faces[4].p[upperNode2Edge[iEdge][1]].y<< " " << segStart.x << " " << segStart.y << " " << segEnd.x << " " << segEnd.y  << std::endl;
				//	OutputFiles::m_logFile<<iElem<<" " << edgeStatus[iEdge] << " " << crossStart[iEdge].x << " " << crossStart[iEdge].y << " " << crossEnd[iEdge].x << " " << crossEnd[iEdge].y << " " << segStart.x << " " << segStart.y << " " << segEnd.x << " " << segEnd.y << std::endl;
				if (edgeStatus[iEdge] == 0 && crossLength > 1e-10) {
					int numOfNeighbourElement(-1);
					numOfNeighbourElement = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, shearEdgeFace[iEdge]);
					if (numOfNeighbourElement <= 1) {
						isEquivalentDipoleInEdge = true;
						crossStart[iEdge].z = getZOnLineSegment(faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossStart[iEdge].x, crossStart[iEdge].y);
						crossEnd[iEdge].z = getZOnLineSegment(faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossEnd[iEdge].x, crossEnd[iEdge].y);
						//calculate parameter
						dipoleParameter.numOfEquivalentDipoleInThisSource++;
						iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
						if (iDipole > maxNumOfEquivalentDipoles - 1) {
							OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
							exit(1);
						}
						crossLength = distance(crossStart[iEdge], crossEnd[iEdge]);
						dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
						dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + crossEnd[iEdge].x);
						dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + crossEnd[iEdge].y);
						dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + crossEnd[iEdge].z);
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
						int IDOfNeighbourElement(-1);
						if (numOfNeighbourElement > 0) {
							IDOfNeighbourElement = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, shearEdgeFace[iEdge], 0);
							IDsOfdoNotFindInThisElement.push_back(IDOfNeighbourElement);
						}
						/*	OutputFiles::m_logFile << numOfNeighbourElement << std::endl;
							OutputFiles::m_logFile << IDOfNeighbourElement << std::endl;
							OutputFiles::m_logFile << iElem << " " << faces[4].p[0].x << " " << faces[4].p[0].y << " " << faces[4].p[0].z << " " << segStart.x << " " << segStart.y << std::endl;
							OutputFiles::m_logFile << iElem << " " << faces[4].p[1].x << " " << faces[4].p[1].y << " " << faces[4].p[1].z << " " << segStart.x << " " << segStart.y << std::endl;
							OutputFiles::m_logFile << iElem << " " << faces[4].p[2].x << " " << faces[4].p[2].y << " " << faces[4].p[2].z << " " << crossStart[iEdge].x << " " << crossStart[iEdge].y << std::endl;
							OutputFiles::m_logFile << iElem << " " << faces[4].p[3].x << " " << faces[4].p[3].y << " " << faces[4].p[3].z << " " << crossEnd[iEdge].x << " " << crossEnd[iEdge].y << std::endl;
							exit(1);*/
						break;
					}
				}
			}
			if (!isEquivalentDipoleInEdge) {
				for (int iEdge = 0; iEdge < 4; iEdge++) {
					if (edgeStatus[iEdge] == 1) {
						//find the only cross point 
						if ((isSegStartInThisTetrahedron && distanceXYPlane(segStart, crossStart[iEdge]) > 1e-10) || (isSegEndInThisTetrahedron && distanceXYPlane(segEnd, crossStart[iEdge]) > 1e-10)) {
							if (isSegStartInThisTetrahedron) {
								segStart.z = getZOnQuadPlane(faces[4].p[0], faces[4].p[1], faces[4].p[2], faces[4].p[3], segStart.x, segStart.y);
							}
							else {
								segEnd.z = getZOnQuadPlane(faces[4].p[0], faces[4].p[1], faces[4].p[2], faces[4].p[3], segEnd.x, segEnd.y);
							}
							crossStart[iEdge].z = getZOnLineSegment(faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossStart[iEdge].x, crossStart[iEdge].y);
							//calculate the parameter
							dipoleParameter.numOfEquivalentDipoleInThisSource++;
							iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
							if (iDipole > maxNumOfEquivalentDipoles - 1) {
								OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
								exit(1);
							}
							if (isSegStartInThisTetrahedron) {
								crossLength = distance(crossStart[iEdge], segStart);
								dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
								dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + segStart.x);
								dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + segStart.y);
								dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + segStart.z);
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
								dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
							}
							else {
								crossLength = distance(crossStart[iEdge], segEnd);
								dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
								dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + segEnd.x);
								dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + segEnd.y);
								dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + segEnd.z);
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
								dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
							}
							/*	OutputFiles::m_logFile << iElem << " " << faces[4].p[0].x << " " << faces[4].p[0].y <<" "<< faces[4].p[0].z<<" "<< crossStart[iEdge].x <<" "<< crossStart[iEdge].y << std::endl;
								OutputFiles::m_logFile << iElem << " " << faces[4].p[1].x << " " << faces[4].p[1].y <<" "<< faces[4].p[1].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;
								OutputFiles::m_logFile << iElem << " " << faces[4].p[2].x << " " << faces[4].p[2].y <<" "<< faces[4].p[2].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;
								OutputFiles::m_logFile << iElem << " " << faces[4].p[3].x << " " << faces[4].p[3].y <<" "<< faces[4].p[3].z<<" "<< segStart.x<<" "<< segStart.y << std::endl;
								exit(1);*/
							break;
						}
					}
				}
			}
			//
		}
		else {//both start and end is out of the face
			for (int iEdge = 0; iEdge < 4; iEdge++) {
				edgeStatus[iEdge] = -1;
			}
			Point3D crossStart[4], crossEnd[4];
			isEquivalentDipoleInEdge = false;
			for (int iEdge = 0; iEdge < 4; iEdge++) {
				edgeStatus[iEdge] = checkSegmentIntersectionInXYPlane(segStart, segEnd, faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossStart[iEdge], crossEnd[iEdge]);
			}
			/*	for (int iEdge = 0; iEdge < 4; iEdge++) {
					OutputFiles::m_logFile << edgeStatus[iEdge] << std::endl;
				}
				OutputFiles::m_logFile << iElem << " " << faces[4].p[0].x << " " << faces[4].p[0].y << " " << faces[4].p[0].z << " " << segStart.x << " " << segStart.y << std::endl;
				OutputFiles::m_logFile << iElem << " " << faces[4].p[1].x << " " << faces[4].p[1].y << " " << faces[4].p[1].z << " " << segStart.x << " " << segStart.y << std::endl;
				OutputFiles::m_logFile << iElem << " " << faces[4].p[2].x << " " << faces[4].p[2].y << " " << faces[4].p[2].z << " " << segStart.x << " " << segStart.y << std::endl;
				OutputFiles::m_logFile << iElem << " " << faces[4].p[3].x << " " << faces[4].p[3].y << " " << faces[4].p[3].z << " " << segStart.x << " " << segStart.y << std::endl;
				exit(1); */
				//check status
			int numOfInterSectionPoints(0);
			for (int iEdge = 0; iEdge < 4; iEdge++) {
				if (edgeStatus[iEdge] == 0) {
					int numOfNeighbourElement(-1);
					numOfNeighbourElement = ptrMeshDataNonConformingHexaElement->getNumNeighborElement(iElem, shearEdgeFace[iEdge]);
					if (numOfNeighbourElement <= 1) {
						isEquivalentDipoleInEdge = true;
						crossStart[iEdge].z = getZOnLineSegment(faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossStart[iEdge].x, crossStart[iEdge].y);
						crossEnd[iEdge].z = getZOnLineSegment(faces[4].p[upperNode2Edge[iEdge][0]], faces[4].p[upperNode2Edge[iEdge][1]], crossEnd[iEdge].x, crossEnd[iEdge].y);
						//calculate parameter
						dipoleParameter.numOfEquivalentDipoleInThisSource++;
						iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
						if (iDipole > maxNumOfEquivalentDipoles - 1) {
							OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
							exit(1);
						}
						crossLength = distance(crossStart[iEdge], crossEnd[iEdge]);
						dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
						dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + crossEnd[iEdge].x);
						dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + crossEnd[iEdge].y);
						dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + crossEnd[iEdge].z);
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
						dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
						dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
						dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
						if (numOfNeighbourElement > 0) {
							int IDOfNeighbourElement = ptrMeshDataNonConformingHexaElement->getIDOfNeighborElement(iElem, shearEdgeFace[iEdge], 0);
							IDsOfdoNotFindInThisElement.push_back(IDOfNeighbourElement);
						}
						break;
					}
				}
				if (edgeStatus[iEdge] > 0)numOfInterSectionPoints++;
			}
			if (!isEquivalentDipoleInEdge && numOfInterSectionPoints > 1) {
				Point3D triangleCrossStart, triangleCrossEnd;
				if (checkIfSegmentCrossTwoPointOfHexaXYPlane(crossStart, triangleCrossStart, triangleCrossEnd)) {
					triangleCrossStart.z = getZOnQuadPlane(faces[4].p[0], faces[4].p[1], faces[4].p[2], faces[4].p[3], triangleCrossStart.x, triangleCrossStart.y);
					triangleCrossEnd.z = getZOnQuadPlane(faces[4].p[0], faces[4].p[1], faces[4].p[2], faces[4].p[3], triangleCrossEnd.x, triangleCrossEnd.y);
					//calculate parameter
					dipoleParameter.numOfEquivalentDipoleInThisSource++;
					iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
					if (iDipole > maxNumOfEquivalentDipoles - 1) {
						OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
						exit(1);
					}
					crossLength = distance(triangleCrossStart, triangleCrossEnd);
					dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
					dipoleParameter.midPoint[iDipole].X = 0.5 * (triangleCrossStart.x + triangleCrossEnd.x);
					dipoleParameter.midPoint[iDipole].Y = 0.5 * (triangleCrossStart.y + triangleCrossEnd.y);
					dipoleParameter.midPoint[iDipole].Z = 0.5 * (triangleCrossStart.z + triangleCrossEnd.z);
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
				}
			}
		}
	}
}

void ObservedDataSourceSTLS::findEquivalentDipolesInTetraMeshOnEarthSurface(int iSeg) {
	int nEdge = 3;
	int upperNode2Edge[3][2];
	upperNode2Edge[0][0] = 0;
	upperNode2Edge[0][1] = 1;
	upperNode2Edge[1][0] = 1;
	upperNode2Edge[1][1] = 2;
	upperNode2Edge[2][0] = 2;
	upperNode2Edge[2][1] = 0;
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	double crossLength;
	std::vector<Point3D> nodesInDuplicate;
	int nElemOnSurface, elemID, faceID, iDipole;
	int* edgeStatus;
	edgeStatus = new int[nEdge];
	bool isEquivalentDipoleInEdge(false);
	bool isThisElemOnEarthSurface(true);
	const MeshDataTetraElement* const pMeshDataTetraElement = (AnalysisControl::getInstance())->getPointerOfMeshDataTetraElement();

	nElemOnSurface = pMeshDataTetraElement->getNumElemOnLandSurface();
	Point3D segStart, segEnd;
	segStart.x = x_start[iSeg]; segEnd.x = x_end[iSeg];
	segStart.y = y_start[iSeg]; segEnd.y = y_end[iSeg];
	segStart.z = z_start[iSeg]; segEnd.z = z_end[iSeg];
	for (int iElem = 0; iElem < nElemOnSurface; iElem++) {
		//
		std::tie(elemID, faceID) = pMeshDataTetraElement->getElemIDAndFaceIDOnSurface(iElem);
		//if point in this Tetra?
		triangleFace faceOnSurface;
		int nodeID;
		int numOfNodesInDuplicate = 0;
		for (int iNode = 0; iNode < 3; iNode++) {
			nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(elemID, faceID, iNode);
			faceOnSurface.p[iNode].x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
			faceOnSurface.p[iNode].y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
			faceOnSurface.p[iNode].z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
			//check duplicate nodes
			if (checkIfPointExistsInVector(faceOnSurface.p[iNode], nodesInDuplicate)) {
				numOfNodesInDuplicate++;
			}
		}
		//check duplicate nodes
		if (numOfNodesInDuplicate > 1) {
			continue;
		}
		bool isSegStartInThisTetrahedron, isSegEndInThisTetrahedron;
		isSegStartInThisTetrahedron = isPointInTriangleXY(faceOnSurface.p, segStart, 1e-10);
		isSegEndInThisTetrahedron = isPointInTriangleXY(faceOnSurface.p, segEnd, 1e-10);
		//starting find cross point on surface!
		//first check both point in the tetrahedron
		if (isSegStartInThisTetrahedron && isSegEndInThisTetrahedron) {
			dipoleParameter.numOfEquivalentDipoleInThisSource++;
			iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
			if (iDipole > maxNumOfEquivalentDipoles - 1) {
				OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
				exit(1);
			}
			segStart.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], segStart.x, segStart.y);
			segEnd.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], segEnd.x, segEnd.y);
			calculateDirecctionVectorOnSurface(segStart, segEnd, iSeg);
			crossLength = distance(segStart, segEnd);
			dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
			dipoleParameter.midPoint[iDipole].X = 0.5 * (segStart.x + segEnd.x);
			dipoleParameter.midPoint[iDipole].Y = 0.5 * (segStart.y + segEnd.y);
			dipoleParameter.midPoint[iDipole].Z = 0.5 * (segStart.z + segEnd.z);
			dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
			dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
			dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
			dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
			dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
			break;
		}
		//second check one point in the tetrahedron
		else if (isSegStartInThisTetrahedron || isSegEndInThisTetrahedron) {

			for (int iEdge = 0; iEdge < 3; iEdge++) {
				edgeStatus[iEdge] = -1;
			}
			Point3D crossStart[3], crossEnd[3];
			isEquivalentDipoleInEdge = false;
			for (int iEdge = 0; iEdge < 3; iEdge++) {
				edgeStatus[iEdge] = checkSegmentIntersectionInXYPlane(segStart, segEnd, faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossStart[iEdge], crossEnd[iEdge]);
				//OutputFiles::m_logFile << edgeStatus[iEdge] << " " << faceOnSurface.p[upperNode2Edge[iEdge][0]].x << " " << faceOnSurface.p[upperNode2Edge[iEdge][0]].y << " " << faceOnSurface.p[upperNode2Edge[iEdge][1]].x << " " << faceOnSurface.p[upperNode2Edge[iEdge][1]].y << " " << crossStart[iEdge].x << " " << crossStart[iEdge].y << " " << crossEnd[iEdge].x << " " << crossEnd[iEdge].y << std::endl;
			}
			//if (isSegStartInThisTetrahedron) {
			//	for (int iNode = 0; iNode < 3; iNode++) {
			//		OutputFiles::m_logFile << elemID << std::endl;
			//		OutputFiles::m_logFile << iNode << " " << faceOnSurface.p[iNode].x << " " << faceOnSurface.p[iNode].y << " " << faceOnSurface.p[iNode].z << std::endl;
			//	}
			//	OutputFiles::m_logFile << segStart.x << " " << segStart.y << " " << segStart.z << std::endl;
			//	//OutputFiles::m_logFile << segEnd.x << " " << segEnd.y << " " << segEnd.z << std::endl;
			//	exit(1);
			//}
			for (int iEdge = 0; iEdge < 3; iEdge++) {
				crossLength = distance(crossStart[iEdge], crossEnd[iEdge]);
				if (edgeStatus[iEdge] == 0 && crossLength > 1e-10) {
					isEquivalentDipoleInEdge = true;
					crossStart[iEdge].z = getZOnLineSegment(faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossStart[iEdge].x, crossStart[iEdge].y);
					crossEnd[iEdge].z = getZOnLineSegment(faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossEnd[iEdge].x, crossEnd[iEdge].y);
					//calculate parameter
					dipoleParameter.numOfEquivalentDipoleInThisSource++;
					iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
					if (iDipole > maxNumOfEquivalentDipoles - 1) {
						OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
						exit(1);
					}
					crossLength = distance(crossStart[iEdge], crossEnd[iEdge]);
					calculateDirecctionVectorOnSurface(segStart, segEnd, iSeg);
					dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
					dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + crossEnd[iEdge].x);
					dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + crossEnd[iEdge].y);
					dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + crossEnd[iEdge].z);
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
					//
					nodesInDuplicate.push_back(faceOnSurface.p[upperNode2Edge[iEdge][0]]);
					nodesInDuplicate.push_back(faceOnSurface.p[upperNode2Edge[iEdge][1]]);
					break;
				}
			}
			if (!isEquivalentDipoleInEdge) {
				for (int iEdge = 0; iEdge < 3; iEdge++) {
					if (edgeStatus[iEdge] == 1) {
						//find the only cross point 
						if ((isSegStartInThisTetrahedron && distanceXYPlane(segStart, crossStart[iEdge]) > 1e-10) || (isSegEndInThisTetrahedron && distanceXYPlane(segEnd, crossStart[iEdge]) > 1e-10)) {
							if (isSegStartInThisTetrahedron) {
								segStart.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], segStart.x, segStart.y);
							}
							else {
								segEnd.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], segEnd.x, segEnd.y);
							}
							crossStart[iEdge].z = getZOnLineSegment(faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossStart[iEdge].x, crossStart[iEdge].y);
							//calculate the parameter
							dipoleParameter.numOfEquivalentDipoleInThisSource++;
							iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
							if (iDipole > maxNumOfEquivalentDipoles - 1) {
								OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
								exit(1);
							}
							if (isSegStartInThisTetrahedron) {
								crossLength = distance(crossStart[iEdge], segStart);
								calculateDirecctionVectorOnSurface(segStart, segEnd, iSeg);
								dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
								dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + segStart.x);
								dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + segStart.y);
								dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + segStart.z);
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
								dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
							}
							else {
								crossLength = distance(crossStart[iEdge], segEnd);
								calculateDirecctionVectorOnSurface(segStart, segEnd, iSeg);
								dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
								dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + segEnd.x);
								dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + segEnd.y);
								dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + segEnd.z);
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
								dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
								dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
								dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
							}
							break;
						}
					}
				}
			}
			//
		}
		//third check both point out of the tetrahedron
		else {
			for (int iEdge = 0; iEdge < 3; iEdge++) {
				edgeStatus[iEdge] = -1;
			}
			Point3D crossStart[3], crossEnd[3];
			isEquivalentDipoleInEdge = false;
			for (int iEdge = 0; iEdge < 3; iEdge++) {
				edgeStatus[iEdge] = checkSegmentIntersectionInXYPlane(segStart, segEnd, faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossStart[iEdge], crossEnd[iEdge]);
			}
			//check status
			int numOfInterSectionPoints(0);
			for (int iEdge = 0; iEdge < 3; iEdge++) {
				if (edgeStatus[iEdge] == 0) {
					isEquivalentDipoleInEdge = true;
					crossStart[iEdge].z = getZOnLineSegment(faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossStart[iEdge].x, crossStart[iEdge].y);
					crossEnd[iEdge].z = getZOnLineSegment(faceOnSurface.p[upperNode2Edge[iEdge][0]], faceOnSurface.p[upperNode2Edge[iEdge][1]], crossEnd[iEdge].x, crossEnd[iEdge].y);
					//calculate parameter
					dipoleParameter.numOfEquivalentDipoleInThisSource++;
					iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
					if (iDipole > maxNumOfEquivalentDipoles - 1) {
						OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
						exit(1);
					}
					crossLength = distance(crossStart[iEdge], crossEnd[iEdge]);
					calculateDirecctionVectorOnSurface(segStart, segEnd, iSeg);
					dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
					dipoleParameter.midPoint[iDipole].X = 0.5 * (crossStart[iEdge].x + crossEnd[iEdge].x);
					dipoleParameter.midPoint[iDipole].Y = 0.5 * (crossStart[iEdge].y + crossEnd[iEdge].y);
					dipoleParameter.midPoint[iDipole].Z = 0.5 * (crossStart[iEdge].z + crossEnd[iEdge].z);
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
					nodesInDuplicate.push_back(faceOnSurface.p[upperNode2Edge[iEdge][0]]);
					nodesInDuplicate.push_back(faceOnSurface.p[upperNode2Edge[iEdge][1]]);
					break;
				}
				if (edgeStatus[iEdge] > 0)numOfInterSectionPoints++;
			}
			if (!isEquivalentDipoleInEdge && numOfInterSectionPoints > 1) {
				Point3D triangleCrossStart, triangleCrossEnd;
				if (checkIfSegmentCrossTwoPointOfThisFace(crossStart, triangleCrossStart, triangleCrossEnd)) {
					triangleCrossStart.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], triangleCrossStart.x, triangleCrossStart.y);
					triangleCrossEnd.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], triangleCrossEnd.x, triangleCrossEnd.y);
					//calculate parameter
					dipoleParameter.numOfEquivalentDipoleInThisSource++;
					iDipole = dipoleParameter.numOfEquivalentDipoleInThisSource - 1;
					if (iDipole > maxNumOfEquivalentDipoles - 1) {
						OutputFiles::m_logFile << "Number of equivalent dipoles exceeds the max number of toralence ." << std::endl;
						exit(1);
					}
					crossLength = distance(triangleCrossStart, triangleCrossEnd);
					calculateDirecctionVectorOnSurface(segStart, segEnd, iSeg);
					dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
					dipoleParameter.midPoint[iDipole].X = 0.5 * (triangleCrossStart.x + triangleCrossEnd.x);
					dipoleParameter.midPoint[iDipole].Y = 0.5 * (triangleCrossStart.y + triangleCrossEnd.y);
					dipoleParameter.midPoint[iDipole].Z = 0.5 * (triangleCrossStart.z + triangleCrossEnd.z);
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][0] = currentIntensity[iSeg];
					dipoleParameter.parameterOfEquivalentDipoleInThisSource[iDipole][1] = crossLength;
					dipoleParameter.directionOfEquivalentDipoles[iDipole][0] = directionX[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][1] = directionY[iSeg];
					dipoleParameter.directionOfEquivalentDipoles[iDipole][2] = directionZ[iSeg];
				}
			}
		}
	}
}

void ObservedDataSourceSTLS::calculateDirecctionVectorOnSurface(Point3D segStart, Point3D segEnd, int iSeg) {
	//calculate direction
	double length(0);
	length = sqrt(pow(segStart.x - segEnd.x, 2) + pow(segStart.y - segEnd.y, 2) + pow(segStart.z - segEnd.z, 2));
	if (length == 0) {
		directionX[iSeg] = 0;
		directionY[iSeg] = 0;
		directionZ[iSeg] = 0;
	}
	else {
		directionX[iSeg] = -(segStart.x - segEnd.x) / length * m_sourceDirection[iSeg];
		directionY[iSeg] = -(segStart.y - segEnd.y) / length * m_sourceDirection[iSeg];
		directionZ[iSeg] = -(segStart.z - segEnd.z) / length * m_sourceDirection[iSeg];
	}
	//check if the direction vector unit
	double units = std::sqrt(directionX[iSeg] * directionX[iSeg] + directionY[iSeg] * directionY[iSeg] + directionZ[iSeg] * directionZ[iSeg]);
	if (abs(units - 1) > 0.001) {
		OutputFiles::m_logFile << "Direction vector don't unit in source ID: " << m_sourceID << std::endl;
		OutputFiles::m_logFile << segStart.x << " " << segStart.y << " " << segStart.z << std::endl;
		OutputFiles::m_logFile << segEnd.x << " " << segEnd.y << " " << segEnd.z << std::endl;
		OutputFiles::m_logFile << length << " " << m_sourceDirection[iSeg] << " " << units << std::endl;
		exit(1);
	}
}