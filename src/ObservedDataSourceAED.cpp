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

#include "ObservedDataSourceAED.h"
#include "ObservedDataSource.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

//-------------------------------------------------------
//----- Definition of the class ObservedDataSourceAED -----
//-------------------------------------------------------

// Constructer
ObservedDataSourceAED::ObservedDataSourceAED() :
	ObservedDataSource(),
	m_sourceDirection(NULL),
	dip(NULL),
	azimuth(NULL),
	m_findDipoleSourceOnEarthSurface(false)
{
}

// Destructer
ObservedDataSourceAED::~ObservedDataSourceAED() {

	if (m_sourceDirection != NULL) {
		delete[] m_sourceDirection;
		m_sourceDirection = NULL;
	}

	if (dip != NULL) {
		delete[] dip;
		dip = NULL;
	}

	if (azimuth != NULL) {
		delete[] azimuth;
		azimuth = NULL;
	}

}

// Read data from input file
void ObservedDataSourceAED::inputObservedData(std::ifstream& inFile) {

	m_findDipoleSourceOnEarthSurface = (AnalysisControl::getInstance())->getIsSourcePartitionedOnSurface();
	const double convKilometerToMeter = 1000.0;
	//single dipole source
	int nED(0);

	//start reading file
	std::string tmpchar;
	inFile >> tmpchar;
	//maybe check here
	inFile >> m_sourceID;
	inFile >> dipoleParameter.numOfEquivalentDipoleInThisSource;
	nED = dipoleParameter.numOfEquivalentDipoleInThisSource;
	//
	dipoleParameter.elementIncludingEquivalentDipoleInThisSource = new int[nED];
	dipoleParameter.sourceID = new int[nED];
	dipoleParameter.midPoint = new CommonParameters::locationXYZ[nED];
	dipoleParameter.sourceTypeOfAllEquivalentDipoleSource = ELECTRIC_DIPOLE_SOURCE;
	dipoleParameter.parameterOfEquivalentDipoleInThisSource = new double* [nED];
	dipoleParameter.directionOfEquivalentDipoles = new double* [nED];
	for (int i = 0; i < nED; i++) {
		dipoleParameter.parameterOfEquivalentDipoleInThisSource[i] = new double[2];
		dipoleParameter.directionOfEquivalentDipoles[i] = new double[3];
	}
	m_sourceDirection = new int[nED];
	azimuth = new double[nED];
	dip = new double[nED];

	for (int i = 0; i < nED; i++) {
		inFile >> m_sourceDirection[i] >> dipoleParameter.parameterOfEquivalentDipoleInThisSource[i][0] >> dipoleParameter.parameterOfEquivalentDipoleInThisSource[i][1];
		//read midpoint of dipole
		inFile >> dipoleParameter.midPoint[i].X >> dipoleParameter.midPoint[i].Y >> dipoleParameter.midPoint[i].Z;
		//read azimuth and dipole
		inFile >> azimuth[i] >> dip[i];
		dipoleParameter.sourceID[i] = i;
		calculateDirecctionVector(i);

		//Convert kilometer to meter
		dipoleParameter.midPoint[i].X *= convKilometerToMeter;
		dipoleParameter.midPoint[i].Y *= convKilometerToMeter;
		dipoleParameter.midPoint[i].Z *= convKilometerToMeter;

		if (!m_findDipoleSourceOnEarthSurface) {
			findElementIncludingDipoleSource(i);
		}
		else {
			findElementIncludingDipoleSourceOnEarthSurface(i);
		}
	}
}

void ObservedDataSourceAED::calculateDirecctionVector(int iED) {
	//calculate direction
	dipoleParameter.directionOfEquivalentDipoles[iED][0] = static_cast<double>(m_sourceDirection[iED]) * std::cos(degreesToRadians(azimuth[iED])) * std::cos(degreesToRadians(dip[iED]));

	dipoleParameter.directionOfEquivalentDipoles[iED][1] = static_cast<double>(m_sourceDirection[iED]) * std::sin(degreesToRadians(azimuth[iED])) * std::cos(degreesToRadians(dip[iED]));

	dipoleParameter.directionOfEquivalentDipoles[iED][2] = static_cast<double>(m_sourceDirection[iED]) * std::sin(degreesToRadians(dip[iED]));

	//check if the direction vector unit
	double units = std::sqrt(dipoleParameter.directionOfEquivalentDipoles[iED][0] * dipoleParameter.directionOfEquivalentDipoles[iED][0]
		+ dipoleParameter.directionOfEquivalentDipoles[iED][1] * dipoleParameter.directionOfEquivalentDipoles[iED][1]
		+ dipoleParameter.directionOfEquivalentDipoles[iED][2] * dipoleParameter.directionOfEquivalentDipoles[iED][2]);
	if (units != 1) {
		OutputFiles::m_logFile << "Direction vector don't unit in source ID: " << dipoleParameter.sourceID[iED] << std::endl;
		exit(1);
	}
}

void ObservedDataSourceAED::findElementIncludingDipoleSource(int iDipole) {
	//initialize parameters
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	const bool modLoc = (AnalysisControl::getInstance())->getIsSourceLocMovedToCenter();
	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	int nElem;
	Point3D dipoleMidpoint;
	CommonParameters::locationXYZ CenterCoord;
	//
	assignmentPoint3D(dipoleMidpoint, dipoleParameter.midPoint[iDipole].X, dipoleParameter.midPoint[iDipole].Y, dipoleParameter.midPoint[iDipole].Z);
	if (meshType == MeshData::TETRA) {
		const MeshDataTetraElement* const pMeshDataTetraElement = (AnalysisControl::getInstance())->getPointerOfMeshDataTetraElement();
		nElem = pMeshDataTetraElement->getNumElemTotal();
		dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = -1;
		for (int iElem = 0; iElem < nElem; iElem++) {
			//jump the source in the air !! < FOR NOW, MAY BE DELETED IN FUTURE? >
			if (log10(pResistivityBlock->getResistivityValuesFromElemID(iElem)) > 7.0) {
				continue;
			}
			//
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
			bool isDipoleInThisTetrahedron;
			isDipoleInThisTetrahedron = isPointInTetrahedron(dipoleMidpoint, tetA, tetB, tetC, tetD, 1e-10);
			//
			if (isDipoleInThisTetrahedron) {//all the segment is in the tetrahedron
				dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
				if (modLoc) {
					CenterCoord = pMeshDataTetraElement->getCenterCoord(iElem);
					dipoleParameter.midPoint[iDipole].X = CenterCoord.X;
					dipoleParameter.midPoint[iDipole].Y = CenterCoord.Y;
					dipoleParameter.midPoint[iDipole].Z = CenterCoord.Z;
				}
				break;
			}
		}
		if (dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] < 0) {
			OutputFiles::m_logFile << "Find dipole in element failed !!" << std::endl;
		}
	}
	else if (meshType == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		nElem = ptrMeshDataNonConformingHexaElement->getNumElemTotal();
		dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = -1;
		for (int iElem = 0; iElem < nElem; iElem++) {
			//jump the source in the air !! < FOR NOW, MAY BE DELETED IN FUTURE? >
			if (log10(pResistivityBlock->getResistivityValuesFromElemID(iElem)) > 7.0) {
				continue;
			}
			//
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
			bool isDipoleInThisNonConformingHexa;
			isDipoleInThisNonConformingHexa = isPointInNonConformingHexa(faces, dipoleMidpoint, 1e-10);
			//
			if (isDipoleInThisNonConformingHexa) {//all the segment is in the tetrahedron
				dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
				if (modLoc) {
					CenterCoord = ptrMeshDataNonConformingHexaElement->getCenterCoord(iElem);
					dipoleParameter.midPoint[iDipole].X = CenterCoord.X;
					dipoleParameter.midPoint[iDipole].Y = CenterCoord.Y;
					dipoleParameter.midPoint[iDipole].Z = CenterCoord.Z;
				}
				break;
			}
		}
		if (dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] < 0) {
			OutputFiles::m_logFile << "Find dipole in element failed !!" << std::endl;
		}
	}
}

void ObservedDataSourceAED::findElementIncludingDipoleSourceOnEarthSurface(int iDipole) {
	//initialize parameters
	const ResistivityBlock* const pResistivityBlock = ResistivityBlock::getInstance();
	const bool modLoc = (AnalysisControl::getInstance())->getIsSourceLocMovedToCenter();
	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	int nElem;
	Point3D dipoleMidpoint;
	CommonParameters::locationXYZ CenterCoord;
	//
	assignmentPoint3D(dipoleMidpoint, dipoleParameter.midPoint[iDipole].X, dipoleParameter.midPoint[iDipole].Y, dipoleParameter.midPoint[iDipole].Z);
	if (meshType == MeshData::TETRA) {
		const MeshDataTetraElement* const pMeshDataTetraElement = (AnalysisControl::getInstance())->getPointerOfMeshDataTetraElement();
		nElem = pMeshDataTetraElement->getNumElemTotal();
		dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = -1;
		int nElemOnSurface = pMeshDataTetraElement->getNumElemOnLandSurface();
		int elemID, faceID;
		for (int iElem = 0; iElem < nElemOnSurface; iElem++) {
			std::tie(elemID, faceID) = pMeshDataTetraElement->getElemIDAndFaceIDOnSurface(iElem);
			//
			triangleFace faceOnSurface;
			int nodeID;
			for (int iNode = 0; iNode < 3; iNode++) {
				nodeID = pMeshDataTetraElement->getNodeIDGlobalFromElementAndFace(elemID, faceID, iNode);
				faceOnSurface.p[iNode].x = pMeshDataTetraElement->getXCoordinatesOfNodes(nodeID);
				faceOnSurface.p[iNode].y = pMeshDataTetraElement->getYCoordinatesOfNodes(nodeID);
				faceOnSurface.p[iNode].z = pMeshDataTetraElement->getZCoordinatesOfNodes(nodeID);
				//check duplicate nodes
				//OutputFiles::m_logFile << faceOnSurface.p[iNode].x << " " << faceOnSurface.p[iNode].y << " " << faceOnSurface.p[iNode].z << std::endl;
			}
		////	OutputFiles::m_logFile << dipoleMidpoint.x << " " << dipoleMidpoint.y << " " << dipoleMidpoint.z << std::endl;
			bool isDipoleInThisTetrahedron = isPointInTriangleXY(faceOnSurface.p, dipoleMidpoint, 1e-10);
		//	//
			if (isDipoleInThisTetrahedron) {//all the segment is in the tetrahedron
				dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = elemID;
				if (modLoc) {
					CenterCoord = pMeshDataTetraElement->getCenterCoord(elemID);
					dipoleParameter.midPoint[iDipole].X = CenterCoord.X;
					dipoleParameter.midPoint[iDipole].Y = CenterCoord.Y;
					dipoleParameter.midPoint[iDipole].Z = CenterCoord.Z;
				}
				else {
					dipoleMidpoint.z = getZOnTrianglePlane(faceOnSurface.p[0], faceOnSurface.p[1], faceOnSurface.p[2], dipoleMidpoint.x, dipoleMidpoint.y);
					dipoleParameter.midPoint[iDipole].Z = dipoleMidpoint.z;
				}
				break;
			}
		}
		if (dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] < 0) {
			OutputFiles::m_logFile << "Find dipole on surface element failed !!" << std::endl;
		}
	}
	else if (meshType == MeshData::NONCONFORMING_HEXA) {
		const MeshDataNonConformingHexaElement* const ptrMeshDataNonConformingHexaElement = (AnalysisControl::getInstance())->getPointerOfMeshDataNonConformingHexaElement();
		nElem = ptrMeshDataNonConformingHexaElement->getNumElemTotal();
		dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = -1;
		for (int iElem = 0; iElem < nElem; iElem++) {
			//jump the source in the air !! < FOR NOW, MAY BE DELETED IN FUTURE? >
			if (log10(pResistivityBlock->getResistivityValuesFromElemID(iElem)) > 7.0) {
				continue;
			}
			//
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
			bool isDipoleInThisNonConformingHexa;
			isDipoleInThisNonConformingHexa = isPointInNonConformingHexa(faces, dipoleMidpoint, 1e-10);
			//
			if (isDipoleInThisNonConformingHexa) {//all the segment is in the tetrahedron
				dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] = iElem;
				if (modLoc) {
					CenterCoord = ptrMeshDataNonConformingHexaElement->getCenterCoord(iElem);
					dipoleParameter.midPoint[iDipole].X = CenterCoord.X;
					dipoleParameter.midPoint[iDipole].Y = CenterCoord.Y;
					dipoleParameter.midPoint[iDipole].Z = CenterCoord.Z;
				}
				break;
			}
		}
		if (dipoleParameter.elementIncludingEquivalentDipoleInThisSource[iDipole] < 0) {
			OutputFiles::m_logFile << "Find dipole in element failed !!" << std::endl;
		}
	}
}
