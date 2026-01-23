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

#include "ObservedDataStationCSEMAmplitudeAndPhase.h"
#include "ObservedDataStationPoint.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationCSEMAmplitudeAndPhase::ObservedDataStationCSEMAmplitudeAndPhase() :
	ObservedDataStationPoint(),
	m_dataTypeToBeUsed(NULL),
	m_ExCalculated(NULL),
	m_EyCalculated(NULL),
	m_HxCalculatedCSEM(NULL),
	m_HyCalculatedCSEM(NULL),
	m_HzCalculated(NULL),
	m_ExCalculatedAmplitude(NULL),
	m_EyCalculatedAmplitude(NULL),
	m_HxCalculatedAmplitude(NULL),
	m_HyCalculatedAmplitude(NULL),
	m_HzCalculatedAmplitude(NULL),
	m_ExCalculatedPhase(NULL),
	m_EyCalculatedPhase(NULL),
	m_HxCalculatedPhase(NULL),
	m_HyCalculatedPhase(NULL),
	m_HzCalculatedPhase(NULL),
	m_ExObservedAmplitude(NULL),
	m_EyObservedAmplitude(NULL),
	m_HxObservedAmplitude(NULL),
	m_HyObservedAmplitude(NULL),
	m_HzObservedAmplitude(NULL),
	m_ExObservedPhase(NULL),
	m_EyObservedPhase(NULL),
	m_HxObservedPhase(NULL),
	m_HyObservedPhase(NULL),
	m_HzObservedPhase(NULL),
	m_ExObservedAmplitudeSD(NULL),
	m_EyObservedAmplitudeSD(NULL),
	m_HxObservedAmplitudeSD(NULL),
	m_HyObservedAmplitudeSD(NULL),
	m_HzObservedAmplitudeSD(NULL),
	m_ExObservedPhaseSD(NULL),
	m_EyObservedPhaseSD(NULL),
	m_HxObservedPhaseSD(NULL),
	m_HyObservedPhaseSD(NULL),
	m_HzObservedPhaseSD(NULL),
	m_ExResidualAmplitude(NULL),
	m_EyResidualAmplitude(NULL),
	m_HxResidualAmplitude(NULL),
	m_HyResidualAmplitude(NULL),
	m_HzResidualAmplitude(NULL),
	m_ExResidualPhase(NULL),
	m_EyResidualPhase(NULL),
	m_HxResidualPhase(NULL),
	m_HyResidualPhase(NULL),
	m_HzResidualPhase(NULL),
	m_rhsVectorIDOfEx(0),
	m_rhsVectorIDOfEy(0),
	m_nTypeOfData(0),
	m_dataIDOfExAmplitude(NULL),
	m_dataIDOfEyAmplitude(NULL),
	m_dataIDOfHxAmplitude(NULL),
	m_dataIDOfHyAmplitude(NULL),
	m_dataIDOfHzAmplitude(NULL),
	m_dataIDOfExPhase(NULL),
	m_dataIDOfEyPhase(NULL),
	m_dataIDOfHxPhase(NULL),
	m_dataIDOfHyPhase(NULL),
	m_dataIDOfHzPhase(NULL),
	m_fixDistortionMatrix(true),
	m_useExData(false),
	m_useEyData(false),
	m_useHxData(false),
	m_useHyData(false),
	m_useHzData(false),
	m_typeOfElectricField(AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD),
	m_arrayDistortionMatrixDifferences(NULL),
	m_arrayGainsAndRotations(NULL)
{

	if ((AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::DISTORTION_TYPE_UNDEFINED) {
		OutputFiles::m_logFile << "Error : Type of distortion must be defined before instantiation of ObservedDataStationMT !!" << std::endl;
		exit(1);
	}
	else if ((AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
		m_arrayDistortionMatrixDifferences = new DistortionMatrixDifferences;
		for (int i = 0; i < 4; ++i) {
			m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[i] = 0.0;
			m_arrayDistortionMatrixDifferences->distortionMatrixDifference[i] = 0.0;
			m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[i] = 0.0;
			m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[i] = -1;
		}
	}
	else if ((AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS ||
		(AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY) {
		m_arrayGainsAndRotations = new GainsAndRotations;
		for (int i = 0; i < 4; ++i) {
			m_arrayGainsAndRotations->gainsAndRotationsPre[i] = 0.0;
			m_arrayGainsAndRotations->gainsAndRotations[i] = 0.0;
			m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[i] = 0.0;
			m_arrayGainsAndRotations->IDsOfGainsAndRotations[i] = -1;
		}
	}

}

// Destructer
ObservedDataStationCSEMAmplitudeAndPhase::~ObservedDataStationCSEMAmplitudeAndPhase() {

	if (m_ExCalculated != NULL) {
		delete[] m_ExCalculated;
		m_ExCalculated = NULL;
	}

	if (m_EyCalculated != NULL) {
		delete[] m_EyCalculated;
		m_EyCalculated = NULL;
	}

	if (m_HxCalculatedCSEM != NULL) {
		delete[] m_HxCalculatedCSEM;
		m_HxCalculatedCSEM = NULL;
	}

	if (m_HyCalculatedCSEM != NULL) {
		delete[] m_HyCalculatedCSEM;
		m_HyCalculatedCSEM = NULL;
	}

	if (m_HzCalculated != NULL) {
		delete[] m_HzCalculated;
		m_HzCalculated = NULL;
	}

	if (m_ExCalculatedAmplitude != NULL) {
		delete[] m_ExCalculatedAmplitude;
		m_ExCalculatedAmplitude = NULL;
	}
	if (m_EyCalculatedAmplitude != NULL) {
		delete[] m_EyCalculatedAmplitude;
		m_EyCalculatedAmplitude = NULL;
	}
	if (m_HxCalculatedAmplitude != NULL) {
		delete[] m_HxCalculatedAmplitude;
		m_HxCalculatedAmplitude = NULL;
	}
	if (m_HyCalculatedAmplitude != NULL) {
		delete[] m_HyCalculatedAmplitude;
		m_HyCalculatedAmplitude = NULL;
	}
	if (m_HzCalculatedAmplitude != NULL) {
		delete[] m_HzCalculatedAmplitude;
		m_HzCalculatedAmplitude = NULL;
	}
	if (m_ExCalculatedPhase != NULL) {
		delete[] m_ExCalculatedPhase;
		m_ExCalculatedPhase = NULL;
	}
	if (m_EyCalculatedPhase != NULL) {
		delete[] m_EyCalculatedPhase;
		m_EyCalculatedPhase = NULL;
	}
	if (m_HxCalculatedPhase != NULL) {
		delete[] m_HxCalculatedPhase;
		m_HxCalculatedPhase = NULL;
	}
	if (m_HyCalculatedPhase != NULL) {
		delete[] m_HyCalculatedPhase;
		m_HyCalculatedPhase = NULL;
	}
	if (m_HzCalculatedPhase != NULL) {
		delete[] m_HzCalculatedPhase;
		m_HzCalculatedPhase = NULL;
	}

	if (m_ExObservedAmplitude != NULL) {
		delete[] m_ExObservedAmplitude;
		m_ExObservedAmplitude = NULL;
	}
	if (m_EyObservedAmplitude != NULL) {
		delete[] m_EyObservedAmplitude;
		m_EyObservedAmplitude = NULL;
	}
	if (m_HxObservedAmplitude != NULL) {
		delete[] m_HxObservedAmplitude;
		m_HxObservedAmplitude = NULL;
	}
	if (m_HyObservedAmplitude != NULL) {
		delete[] m_HyObservedAmplitude;
		m_HyObservedAmplitude = NULL;
	}
	if (m_HzObservedAmplitude != NULL) {
		delete[] m_HzObservedAmplitude;
		m_HzObservedAmplitude = NULL;
	}
	if (m_ExObservedPhase != NULL) {
		delete[] m_ExObservedPhase;
		m_ExObservedPhase = NULL;
	}
	if (m_EyObservedPhase != NULL) {
		delete[] m_EyObservedPhase;
		m_EyObservedPhase = NULL;
	}
	if (m_HxObservedPhase != NULL) {
		delete[] m_HxObservedPhase;
		m_HxObservedPhase = NULL;
	}
	if (m_HyObservedPhase != NULL) {
		delete[] m_HyObservedPhase;
		m_HyObservedPhase = NULL;
	}
	if (m_HzObservedPhase != NULL) {
		delete[] m_HzObservedPhase;
		m_HzObservedPhase = NULL;
	}
	if (m_ExObservedAmplitudeSD != NULL) {
		delete[] m_ExObservedAmplitudeSD;
		m_ExObservedAmplitudeSD = NULL;
	}
	if (m_EyObservedAmplitudeSD != NULL) {
		delete[] m_EyObservedAmplitudeSD;
		m_EyObservedAmplitudeSD = NULL;
	}
	if (m_HxObservedAmplitudeSD != NULL) {
		delete[] m_HxObservedAmplitudeSD;
		m_HxObservedAmplitudeSD = NULL;
	}
	if (m_HyObservedAmplitudeSD != NULL) {
		delete[] m_HyObservedAmplitudeSD;
		m_HyObservedAmplitudeSD = NULL;
	}
	if (m_HzObservedAmplitudeSD != NULL) {
		delete[] m_HzObservedAmplitudeSD;
		m_HzObservedAmplitudeSD = NULL;
	}
	if (m_ExObservedPhaseSD != NULL) {
		delete[] m_ExObservedPhaseSD;
		m_ExObservedPhaseSD = NULL;
	}
	if (m_EyObservedPhaseSD != NULL) {
		delete[] m_EyObservedPhaseSD;
		m_EyObservedPhaseSD = NULL;
	}
	if (m_HxObservedPhaseSD != NULL) {
		delete[] m_HxObservedPhaseSD;
		m_HxObservedPhaseSD = NULL;
	}
	if (m_HyObservedPhaseSD != NULL) {
		delete[] m_HyObservedPhaseSD;
		m_HyObservedPhaseSD = NULL;
	}
	if (m_HzObservedPhaseSD != NULL) {
		delete[] m_HzObservedPhaseSD;
		m_HzObservedPhaseSD = NULL;
	}
	if (m_ExResidualAmplitude != NULL) {
		delete[] m_ExResidualAmplitude;
		m_ExResidualAmplitude = NULL;
	}
	if (m_EyResidualAmplitude != NULL) {
		delete[] m_EyResidualAmplitude;
		m_EyResidualAmplitude = NULL;
	}
	if (m_HxResidualAmplitude != NULL) {
		delete[] m_HxResidualAmplitude;
		m_HxResidualAmplitude = NULL;
	}
	if (m_HyResidualAmplitude != NULL) {
		delete[] m_HyResidualAmplitude;
		m_HyResidualAmplitude = NULL;
	}
	if (m_HzResidualAmplitude != NULL) {
		delete[] m_HzResidualAmplitude;
		m_HzResidualAmplitude = NULL;
	}
	if (m_ExResidualPhase != NULL) {
		delete[] m_ExResidualPhase;
		m_ExResidualPhase = NULL;
	}
	if (m_EyResidualPhase != NULL) {
		delete[] m_EyResidualPhase;
		m_EyResidualPhase = NULL;
	}
	if (m_HxResidualPhase != NULL) {
		delete[] m_HxResidualPhase;
		m_HxResidualPhase = NULL;
	}
	if (m_HyResidualPhase != NULL) {
		delete[] m_HyResidualPhase;
		m_HyResidualPhase = NULL;
	}
	if (m_HzResidualPhase != NULL) {
		delete[] m_HzResidualPhase;
		m_HzResidualPhase = NULL;
	}

	if (m_dataIDOfExAmplitude != NULL) {
		delete[] m_dataIDOfExAmplitude;
		m_dataIDOfExAmplitude = NULL;
	}
	if (m_dataIDOfEyAmplitude != NULL) {
		delete[] m_dataIDOfEyAmplitude;
		m_dataIDOfEyAmplitude = NULL;
	}
	if (m_dataIDOfHxAmplitude != NULL) {
		delete[] m_dataIDOfHxAmplitude;
		m_dataIDOfHxAmplitude = NULL;
	}
	if (m_dataIDOfHyAmplitude != NULL) {
		delete[] m_dataIDOfHyAmplitude;
		m_dataIDOfHyAmplitude = NULL;
	}
	if (m_dataIDOfHzAmplitude != NULL) {
		delete[] m_dataIDOfHzAmplitude;
		m_dataIDOfHzAmplitude = NULL;
	}
	if (m_dataIDOfExPhase != NULL) {
		delete[] m_dataIDOfExPhase;
		m_dataIDOfExPhase = NULL;
	}
	if (m_dataIDOfEyPhase != NULL) {
		delete[] m_dataIDOfEyPhase;
		m_dataIDOfEyPhase = NULL;
	}
	if (m_dataIDOfHxPhase != NULL) {
		delete[] m_dataIDOfHxPhase;
		m_dataIDOfHxPhase = NULL;
	}
	if (m_dataIDOfHyPhase != NULL) {
		delete[] m_dataIDOfHyPhase;
		m_dataIDOfHyPhase = NULL;
	}
	if (m_dataIDOfHzPhase != NULL) {
		delete[] m_dataIDOfHzPhase;
		m_dataIDOfHzPhase = NULL;
	}

	if (m_dataTypeToBeUsed != NULL) {
		delete[] m_dataTypeToBeUsed;
		m_dataTypeToBeUsed = NULL;
	}

	if (m_arrayDistortionMatrixDifferences != NULL) {
		delete m_arrayDistortionMatrixDifferences;
		m_arrayDistortionMatrixDifferences = NULL;
	}

	if (m_arrayGainsAndRotations != NULL) {
		delete m_arrayGainsAndRotations;
		m_arrayGainsAndRotations = NULL;
	}

}

// Read data from input file
void ObservedDataStationCSEMAmplitudeAndPhase::inputObservedData(std::ifstream& inFile) {

	inFile >> m_stationID;

	inFile >> m_IDOfMagneticFieldStation;

	OutputFiles::m_logFile << "# " << std::setw(15) << std::left << m_sourceID << std::setw(15) << std::left << m_stationID << std::setw(18) << std::left << m_IDOfMagneticFieldStation;

	const AnalysisControl* const pAnalysisControl = AnalysisControl::getInstance();
	// Type of owner element
	int ownerType(-1);
	std::string ownerElemType;
	if (pAnalysisControl->isTypeOfOwnerElementSetIndivisually()) {
		// Owner element type of each site is specified
		inFile >> ownerType;
	}
	else {
		ownerType = pAnalysisControl->getTypeOfOwnerElement();
	}
	switch (ownerType) {
	case AnalysisControl::USE_LOWER_ELEMENT:
		m_useUpperElementForInterpolationOfEMField = false;
		ownerElemType = "Lower";
		break;
	case AnalysisControl::USE_UPPER_ELEMENT:
		m_useUpperElementForInterpolationOfEMField = true;
		ownerElemType = "Upper";
		break;
	default:
		OutputFiles::m_logFile << std::endl << "Error : Unknown type of owner element : " << ownerType << std::endl;
		exit(1);
		break;
	}
	OutputFiles::m_logFile << std::setw(15) << std::left << ownerElemType;

	// Type of electric field
	if (pAnalysisControl->isTypeOfElectricFieldSetIndivisually()) {
		// Electric field type of each site is specified
		inFile >> m_typeOfElectricField;
	}
	else {
		m_typeOfElectricField = pAnalysisControl->getTypeOfElectricField();
	}
	if (pAnalysisControl->getTypeOfMesh() == MeshData::HEXA) {
		if (m_typeOfElectricField != AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD) {
			OutputFiles::m_logFile << std::endl << "Warning : Horizontal electric field must be used for hexahedral mesh." << std::endl;
		}
		m_typeOfElectricField = AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD;
	}
	std::string elecType;
	switch (m_typeOfElectricField) {
	case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
		elecType = "Tangential";
		break;
	case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
		elecType = "Horizontal";
		break;
	default:
		OutputFiles::m_logFile << std::endl << "Error : Unknown type of the electric field : " << m_typeOfElectricField << std::endl;
		exit(1);
		break;
	}
	OutputFiles::m_logFile << std::left << elecType << std::endl;

	double dbuf(0.0);
	inFile >> dbuf;
	m_location.X = dbuf * CommonParameters::convKilometerToMeter;
	inFile >> dbuf;
	m_location.Y = dbuf * CommonParameters::convKilometerToMeter;

	//read type of data
	inFile >> m_nTypeOfData;
	m_dataTypeToBeUsed = new int[m_nTypeOfData];
	for (int i = 0; i < m_nTypeOfData; i++) {
		inFile >> m_dataTypeToBeUsed[i];
	}
	if (isDataTypeDuplicated(m_dataTypeToBeUsed)) {
		OutputFiles::m_logFile << "CSEM data type duplicated in station: " << m_stationID << " !!" << std::endl;
		exit(1);
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, EX_AMPPHA)) {
		m_useExData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, EY_AMPPHA)) {
		m_useEyData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, HX_AMPPHA)) {
		m_useHxData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, HY_AMPPHA)) {
		m_useHyData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, HZ_AMPPHA)) {
		m_useHzData = true;
	}

	//read frequency
	inFile >> m_numOfFrequency;

	const int nFreq = m_numOfFrequency;

	if (nFreq > 0) {
		m_freq = new double[nFreq];
		m_ExObservedAmplitude = new double[nFreq];
		m_ExObservedPhase = new double[nFreq];
		m_ExObservedAmplitudeSD = new double[nFreq];
		m_ExObservedPhaseSD = new double[nFreq];
		m_EyObservedAmplitude = new double[nFreq];
		m_EyObservedPhase = new double[nFreq];
		m_EyObservedAmplitudeSD = new double[nFreq];
		m_EyObservedPhaseSD = new double[nFreq];
		m_HxObservedAmplitude = new double[nFreq];
		m_HxObservedPhase = new double[nFreq];
		m_HxObservedAmplitudeSD = new double[nFreq];
		m_HxObservedPhaseSD = new double[nFreq];
		m_HyObservedAmplitude = new double[nFreq];
		m_HyObservedPhase = new double[nFreq];
		m_HyObservedAmplitudeSD = new double[nFreq];
		m_HyObservedPhaseSD = new double[nFreq];
		m_HzObservedAmplitude = new double[nFreq];
		m_HzObservedPhase = new double[nFreq];
		m_HzObservedAmplitudeSD = new double[nFreq];
		m_HzObservedPhaseSD = new double[nFreq];
		//temporary variables to read data

		for (int i = 0; i < nFreq; ++i) {
			m_ExObservedAmplitude[i] = defaultValue;
			m_ExObservedPhase[i] = defaultValue;
			m_ExObservedAmplitudeSD[i] = defaultValue;
			m_ExObservedPhaseSD[i] = defaultValue;
			m_EyObservedAmplitude[i] = defaultValue;
			m_EyObservedPhase[i] = defaultValue;
			m_EyObservedAmplitudeSD[i] = defaultValue;
			m_EyObservedPhaseSD[i] = defaultValue;
			m_HxObservedAmplitude[i] = defaultValue;
			m_HxObservedPhase[i] = defaultValue;
			m_HxObservedAmplitudeSD[i] = defaultValue;
			m_HxObservedPhaseSD[i] = defaultValue;
			m_HyObservedAmplitude[i] = defaultValue;
			m_HyObservedPhase[i] = defaultValue;
			m_HyObservedAmplitudeSD[i] = defaultValue;
			m_HyObservedPhaseSD[i] = defaultValue;
			m_HzObservedAmplitude[i] = defaultValue;
			m_HzObservedPhase[i] = defaultValue;
			m_HzObservedAmplitudeSD[i] = defaultValue;
			m_HzObservedPhaseSD[i] = defaultValue;
		}

		for (int i = 0; i < nFreq; ++i) {
			inFile >> m_freq[i];
			//read data
			for (int iD = 0; iD < m_nTypeOfData; iD++) {
				double dbuf1(0.0);
				double dbuf2(0.0);
				inFile >> dbuf1 >> dbuf2;
				switch (m_dataTypeToBeUsed[iD]) {
				case EX_AMPPHA:
					m_ExObservedAmplitude[i] = dbuf1;
					m_ExObservedPhase[i] = dbuf2;
					break;
				case EY_AMPPHA:
					m_EyObservedAmplitude[i] = dbuf1;
					m_EyObservedPhase[i] = dbuf2;
					break;
				case HX_AMPPHA:
					m_HxObservedAmplitude[i] = dbuf1;
					m_HxObservedPhase[i] = dbuf2;
					break;
				case HY_AMPPHA:
					m_HyObservedAmplitude[i] = dbuf1;
					m_HyObservedPhase[i] = dbuf2;
					break;
				case HZ_AMPPHA:
					m_HzObservedAmplitude[i] = dbuf1;
					m_HzObservedPhase[i] = dbuf2;
					break;
				}
			}
			//read error
			for (int iD = 0; iD < m_nTypeOfData; iD++) {
				switch (m_dataTypeToBeUsed[iD]) {
				case EX_AMPPHA:
					inFile >> m_ExObservedAmplitudeSD[i];
					inFile >> m_ExObservedPhaseSD[i];
					break;
				case EY_AMPPHA:
					inFile >> m_EyObservedAmplitudeSD[i];
					inFile >> m_EyObservedPhaseSD[i];
					break;
				case HX_AMPPHA:
					inFile >> m_HxObservedAmplitudeSD[i];
					inFile >> m_HxObservedPhaseSD[i];
					break;
				case HY_AMPPHA:
					inFile >> m_HyObservedAmplitudeSD[i];
					inFile >> m_HyObservedPhaseSD[i];
					break;
				case HZ_AMPPHA:
					inFile >> m_HzObservedAmplitudeSD[i];
					inFile >> m_HzObservedPhaseSD[i];
					break;
				}
			}
		}
	}

}

// Calulate electric field
void ObservedDataStationCSEMAmplitudeAndPhase::calculateElectricField(const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy, const int ifre) {

	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	if (meshType == MeshData::TETRA) {// Tetra mesh
		switch (getTypeOfElectricField()) {
		case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
			m_ExCalculated[ifre] = ptrForward3D->calcValueElectricFieldXDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3);
			m_EyCalculated[ifre] = ptrForward3D->calcValueElectricFieldYDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3);
			break;
		case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
			m_ExCalculated[ifre] = ptrForward3D->calcValueElectricFieldTangentialX(m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2);
			m_EyCalculated[ifre] = ptrForward3D->calcValueElectricFieldTangentialY(m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2);
			break;
		default:
			OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
			exit(1);
			break;
		}
	}
	else if (meshType == MeshData::HEXA) {// Hexa mesh
		//do not support!!
		OutputFiles::m_logFile << "Error : CSEM module does not support the BRICK mesh!" << std::endl;
		exit(1);
		/*	m_ExCalculated[ifre] = ptrForward3D->calcValueElectricFieldXDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);
			m_EyCalculated[ifre] = ptrForward3D->calcValueElectricFieldYDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);*/
	}
	else if (meshType == MeshData::NONCONFORMING_HEXA) {
		switch (getTypeOfElectricField()) {
		case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
			m_ExCalculated[ifre] = ptrForward3D->calcValueElectricFieldXDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);
			m_EyCalculated[ifre] = ptrForward3D->calcValueElectricFieldYDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);
			break;
		case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
			m_ExCalculated[ifre] = ptrForward3D->calcValueElectricFieldTangentialXFromAllEdges(m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);
			m_EyCalculated[ifre] = ptrForward3D->calcValueElectricFieldTangentialYFromAllEdges(m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);
			break;
		default:
			OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
			exit(1);
			break;
		}
	}
	else {
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

	// For inversion
	if (m_useExData) {
		m_rhsVectorIDOfEx = rhsVectorIDOfEx;
	}
	if (m_useEyData) {
		m_rhsVectorIDOfEy = rhsVectorIDOfEy;
	}

}

// Calulate vertical magnetic field
void ObservedDataStationCSEMAmplitudeAndPhase::calculateVerticalMagneticField(const Forward3D* const ptrForward3D, const int rhsVectorIDOfHz, const int ifre) {

	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	if (meshType == MeshData::TETRA) {// Tetra mesh
		m_HzCalculated[ifre] = ptrForward3D->calcValueMagneticFieldZDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3);
	}
	else if (meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA) {
		m_HzCalculated[ifre] = ptrForward3D->calcValueMagneticFieldZDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z);
	}
	else {
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

	// For inversion
	if (m_useHzData) {
		m_rhsVectorIDOfHz = rhsVectorIDOfHz;
	}
}

// Calulate Impedance tensor
void ObservedDataStationCSEMAmplitudeAndPhase::calculateDataTensor(const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount) {

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE(freq);

	if (freqIDThisPEInSta < 0) {// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[freqIDThisPEInSta];

	m_HxCalculatedCSEM[freqIDThisPEInSta] = ptrStationOfMagneticField->getHxCalculated(0);
	m_HyCalculatedCSEM[freqIDThisPEInSta] = ptrStationOfMagneticField->getHyCalculated(0);

	m_ExCalculatedAmplitude[freqIDThisPEInSta] = log10(std::abs(m_ExCalculated[freqIDGlobalInSta]));
	m_ExCalculatedPhase[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2(m_ExCalculated[freqIDGlobalInSta].imag(), m_ExCalculated[freqIDGlobalInSta].real());
	m_EyCalculatedAmplitude[freqIDThisPEInSta] = log10(std::abs(m_EyCalculated[freqIDGlobalInSta]));
	m_EyCalculatedPhase[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2(m_EyCalculated[freqIDGlobalInSta].imag(), m_EyCalculated[freqIDGlobalInSta].real());
	m_HxCalculatedAmplitude[freqIDThisPEInSta] = log10(std::abs(m_HxCalculatedCSEM[freqIDGlobalInSta]));
	m_HxCalculatedPhase[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2(m_HxCalculatedCSEM[freqIDGlobalInSta].imag(), m_HxCalculatedCSEM[freqIDGlobalInSta].real());
	m_HyCalculatedAmplitude[freqIDThisPEInSta] = log10(std::abs(m_HyCalculatedCSEM[freqIDGlobalInSta]));
	m_HyCalculatedPhase[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2(m_HyCalculatedCSEM[freqIDGlobalInSta].imag(), m_HyCalculatedCSEM[freqIDGlobalInSta].real());
	m_HzCalculatedAmplitude[freqIDThisPEInSta] = log10(std::abs(m_HzCalculated[freqIDGlobalInSta]));
	m_HzCalculatedPhase[freqIDThisPEInSta] = CommonParameters::rad2deg * atan2(m_HzCalculated[freqIDGlobalInSta].imag(), m_HzCalculated[freqIDGlobalInSta].real());

	if (m_useExData) {
		m_ExResidualAmplitude[freqIDThisPEInSta] = (m_ExObservedAmplitude[freqIDGlobalInSta] - m_ExCalculatedAmplitude[freqIDThisPEInSta]) / m_ExObservedAmplitudeSD[freqIDGlobalInSta];
		m_ExResidualPhase[freqIDThisPEInSta] = (m_ExObservedPhase[freqIDGlobalInSta] - m_ExCalculatedPhase[freqIDThisPEInSta]) / m_ExObservedPhaseSD[freqIDGlobalInSta];
		m_dataIDOfExAmplitude[freqIDThisPEInSta] = icount++;
		m_dataIDOfExPhase[freqIDThisPEInSta] = icount++;
	}

	if (m_useEyData) {
		m_EyResidualAmplitude[freqIDThisPEInSta] = (m_EyObservedAmplitude[freqIDGlobalInSta] - m_EyCalculatedAmplitude[freqIDThisPEInSta]) / m_EyObservedAmplitudeSD[freqIDGlobalInSta];
		m_EyResidualPhase[freqIDThisPEInSta] = (m_EyObservedPhase[freqIDGlobalInSta] - m_EyCalculatedPhase[freqIDThisPEInSta]) / m_EyObservedPhaseSD[freqIDGlobalInSta];
		m_dataIDOfEyAmplitude[freqIDThisPEInSta] = icount++;
		m_dataIDOfEyPhase[freqIDThisPEInSta] = icount++;
	}

	if (m_useHxData) {
		m_HxResidualAmplitude[freqIDThisPEInSta] = (m_HxObservedAmplitude[freqIDGlobalInSta] - m_HxCalculatedAmplitude[freqIDThisPEInSta]) / m_HxObservedAmplitudeSD[freqIDGlobalInSta];
		m_HxResidualPhase[freqIDThisPEInSta] = (m_HxObservedPhase[freqIDGlobalInSta] - m_HxCalculatedPhase[freqIDThisPEInSta]) / m_HxObservedPhaseSD[freqIDGlobalInSta];
		m_dataIDOfHxAmplitude[freqIDThisPEInSta] = icount++;
		m_dataIDOfHxPhase[freqIDThisPEInSta] = icount++;
	}

	if (m_useHyData) {
		m_HyResidualAmplitude[freqIDThisPEInSta] = (m_HyObservedAmplitude[freqIDGlobalInSta] - m_HyCalculatedAmplitude[freqIDThisPEInSta]) / m_HyObservedAmplitudeSD[freqIDGlobalInSta];
		m_HyResidualPhase[freqIDThisPEInSta] = (m_HyObservedPhase[freqIDGlobalInSta] - m_HyCalculatedPhase[freqIDThisPEInSta]) / m_HyObservedPhaseSD[freqIDGlobalInSta];
		m_dataIDOfHyAmplitude[freqIDThisPEInSta] = icount++;
		m_dataIDOfHyPhase[freqIDThisPEInSta] = icount++;
	}

	if (m_useHzData) {
		m_HzResidualAmplitude[freqIDThisPEInSta] = (m_HzObservedAmplitude[freqIDGlobalInSta] - m_HzCalculatedAmplitude[freqIDThisPEInSta]) / m_HzObservedAmplitudeSD[freqIDGlobalInSta];
		m_HzResidualPhase[freqIDThisPEInSta] = (m_HzObservedPhase[freqIDGlobalInSta] - m_HzCalculatedPhase[freqIDThisPEInSta]) / m_HzObservedPhaseSD[freqIDGlobalInSta];
		m_dataIDOfHzAmplitude[freqIDThisPEInSta] = icount++;
		m_dataIDOfHzPhase[freqIDThisPEInSta] = icount++;
	}

	//std::complex<double> cxx(1.0, 0.0);
	//std::complex<double> cxy(0.0, 0.0);
	//std::complex<double> cyx(0.0, 0.0);
	//std::complex<double> cyy(1.0, 0.0);
	//const int typeOfDistortion = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//if (typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);
	//	cxx = std::complex<double>(m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXX] + 1.0, 0.0);
	//	cxy = std::complex<double>(m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXY], 0.0);
	//	cyx = std::complex<double>(m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYX], 0.0);
	//	cyy = std::complex<double>(m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYY] + 1.0, 0.0);
	//}
	//else if (typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	const double gX = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN]);
	//	const double gY = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN]);
	//	const double betaX = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_ROTATION];
	//	const double betaY = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_ROTATION];

	//	cxx = std::complex<double>(gX * cos(betaX), 0.0);
	//	cxy = std::complex<double>(-gY * sin(betaY), 0.0);
	//	cyx = std::complex<double>(gX * sin(betaX), 0.0);
	//	cyy = std::complex<double>(gY * cos(betaY), 0.0);
	//}
	//else if (typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	const double gX = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN]);
	//	const double gY = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN]);

	//	cxx = std::complex<double>(gX, 0.0);
	//	cxy = std::complex<double>(0.0, 0.0);
	//	cyx = std::complex<double>(0.0, 0.0);
	//	cyy = std::complex<double>(gY, 0.0);
	//}


#ifdef _DEBUG_WRITE
	std::cout << "freqIDThisPEInSta Zxx Zxy Zyx Zyy : " << freqIDThisPEInSta << " " << m_ZxxCalculated[freqIDThisPEInSta] << " " << m_ZxyCalculated[freqIDThisPEInSta] << " " << m_ZyxCalculated[freqIDThisPEInSta] << " " << m_ZyyCalculated[freqIDThisPEInSta] << std::endl;
#endif

	// For inversion


#ifdef _DEBUG_WRITE
	std::cout << "m_dataIDOfZxx[freqIDThisPEInSta].realPart : " << m_dataIDOfZxx[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZxx[freqIDThisPEInSta].imagPart : " << m_dataIDOfZxx[freqIDThisPEInSta].imagPart << std::endl;
	std::cout << "m_dataIDOfZxy[freqIDThisPEInSta].realPart : " << m_dataIDOfZxy[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZxy[freqIDThisPEInSta].imagPart : " << m_dataIDOfZxy[freqIDThisPEInSta].imagPart << std::endl;
	std::cout << "m_dataIDOfZyx[freqIDThisPEInSta].realPart : " << m_dataIDOfZyx[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZyx[freqIDThisPEInSta].imagPart : " << m_dataIDOfZyx[freqIDThisPEInSta].imagPart << std::endl;
	std::cout << "m_dataIDOfZyy[freqIDThisPEInSta].realPart : " << m_dataIDOfZyy[freqIDThisPEInSta].realPart << std::endl;
	std::cout << "m_dataIDOfZyy[freqIDThisPEInSta].imagPart : " << m_dataIDOfZyy[freqIDThisPEInSta].imagPart << std::endl;
#endif

}

// Initialize electric field
void ObservedDataStationCSEMAmplitudeAndPhase::initializeElectricField(const int freq) {

	m_ExCalculated[freq] = std::complex<double>(0.0, 0.0);
	m_EyCalculated[freq] = std::complex<double>(0.0, 0.0);

}

// Initialize vertical magnetic field
void ObservedDataStationCSEMAmplitudeAndPhase::initializeVerticalMagneticField(const int freq) {

	m_HzCalculated[freq] = std::complex<double>(0.0, 0.0);

}

// Initialize Impedance tensor and errors
void ObservedDataStationCSEMAmplitudeAndPhase::initializeDataTensorsAndErrors() {

	for (int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i) {
		m_ExResidualAmplitude[i] = 0.0;
		m_EyResidualAmplitude[i] = 0.0;
		m_HxResidualAmplitude[i] = 0.0;
		m_HyResidualAmplitude[i] = 0.0;
		m_HzResidualAmplitude[i] = 0.0;
		m_ExResidualPhase[i] = 0.0;
		m_EyResidualPhase[i] = 0.0;
		m_HxResidualPhase[i] = 0.0;
		m_HyResidualPhase[i] = 0.0;
		m_HzResidualPhase[i] = 0.0;

		m_ExCalculatedAmplitude[i] = 0.0;
		m_EyCalculatedAmplitude[i] = 0.0;
		m_HxCalculatedAmplitude[i] = 0.0;
		m_HyCalculatedAmplitude[i] = 0.0;
		m_HzCalculatedAmplitude[i] = 0.0;
		m_ExCalculatedPhase[i] = 0.0;
		m_EyCalculatedPhase[i] = 0.0;
		m_HxCalculatedPhase[i] = 0.0;
		m_HyCalculatedPhase[i] = 0.0;
		m_HzCalculatedPhase[i] = 0.0;

	}

}

// Allocate memory for the calculated Impedance and errors
void ObservedDataStationCSEMAmplitudeAndPhase::allocateMemoryForCalculatedValues() {

	if (m_ExCalculated != NULL) {
		delete[] m_ExCalculated;
		m_ExCalculated = NULL;
	}

	if (m_EyCalculated != NULL) {
		delete[] m_EyCalculated;
		m_EyCalculated = NULL;
	}

	if (m_HxCalculatedCSEM != NULL) {
		delete[] m_HxCalculatedCSEM;
		m_HxCalculatedCSEM = NULL;
	}

	if (m_HyCalculatedCSEM != NULL) {
		delete[] m_HyCalculatedCSEM;
		m_HyCalculatedCSEM = NULL;
	}

	if (m_ExCalculatedAmplitude != NULL) {
		delete[] m_ExCalculatedAmplitude;
		m_ExCalculatedAmplitude = NULL;
	}
	if (m_EyCalculatedAmplitude != NULL) {
		delete[] m_EyCalculatedAmplitude;
		m_EyCalculatedAmplitude = NULL;
	}
	if (m_HxCalculatedAmplitude != NULL) {
		delete[] m_HxCalculatedAmplitude;
		m_HxCalculatedAmplitude = NULL;
	}
	if (m_HyCalculatedAmplitude != NULL) {
		delete[] m_HyCalculatedAmplitude;
		m_HyCalculatedAmplitude = NULL;
	}
	if (m_HzCalculatedAmplitude != NULL) {
		delete[] m_HzCalculatedAmplitude;
		m_HzCalculatedAmplitude = NULL;
	}
	if (m_ExCalculatedPhase != NULL) {
		delete[] m_ExCalculatedPhase;
		m_ExCalculatedPhase = NULL;
	}
	if (m_EyCalculatedPhase != NULL) {
		delete[] m_EyCalculatedPhase;
		m_EyCalculatedPhase = NULL;
	}
	if (m_HxCalculatedPhase != NULL) {
		delete[] m_HxCalculatedPhase;
		m_HxCalculatedPhase = NULL;
	}
	if (m_HyCalculatedPhase != NULL) {
		delete[] m_HyCalculatedPhase;
		m_HyCalculatedPhase = NULL;
	}
	if (m_HzCalculatedPhase != NULL) {
		delete[] m_HzCalculatedPhase;
		m_HzCalculatedPhase = NULL;
	}

	if (m_ExResidualAmplitude != NULL) {
		delete[] m_ExResidualAmplitude;
		m_ExResidualAmplitude = NULL;
	}
	if (m_EyResidualAmplitude != NULL) {
		delete[] m_EyResidualAmplitude;
		m_EyResidualAmplitude = NULL;
	}
	if (m_HxResidualAmplitude != NULL) {
		delete[] m_HxResidualAmplitude;
		m_HxResidualAmplitude = NULL;
	}
	if (m_HyResidualAmplitude != NULL) {
		delete[] m_HyResidualAmplitude;
		m_HyResidualAmplitude = NULL;
	}
	if (m_HzResidualAmplitude != NULL) {
		delete[] m_HzResidualAmplitude;
		m_HzResidualAmplitude = NULL;
	}
	if (m_ExResidualPhase != NULL) {
		delete[] m_ExResidualPhase;
		m_ExResidualPhase = NULL;
	}
	if (m_EyResidualPhase != NULL) {
		delete[] m_EyResidualPhase;
		m_EyResidualPhase = NULL;
	}
	if (m_HxResidualPhase != NULL) {
		delete[] m_HxResidualPhase;
		m_HxResidualPhase = NULL;
	}
	if (m_HyResidualPhase != NULL) {
		delete[] m_HyResidualPhase;
		m_HyResidualPhase = NULL;
	}
	if (m_HzResidualPhase != NULL) {
		delete[] m_HzResidualPhase;
		m_HzResidualPhase = NULL;
	}
	if (m_dataIDOfExAmplitude != NULL) {
		delete[] m_dataIDOfExAmplitude;
		m_dataIDOfExAmplitude = NULL;
	}
	if (m_dataIDOfEyAmplitude != NULL) {
		delete[] m_dataIDOfEyAmplitude;
		m_dataIDOfEyAmplitude = NULL;
	}
	if (m_dataIDOfHxAmplitude != NULL) {
		delete[] m_dataIDOfHxAmplitude;
		m_dataIDOfHxAmplitude = NULL;
	}
	if (m_dataIDOfHyAmplitude != NULL) {
		delete[] m_dataIDOfHyAmplitude;
		m_dataIDOfHyAmplitude = NULL;
	}
	if (m_dataIDOfHzAmplitude != NULL) {
		delete[] m_dataIDOfHzAmplitude;
		m_dataIDOfHzAmplitude = NULL;
	}
	if (m_dataIDOfExPhase != NULL) {
		delete[] m_dataIDOfExPhase;
		m_dataIDOfExPhase = NULL;
	}
	if (m_dataIDOfEyPhase != NULL) {
		delete[] m_dataIDOfEyPhase;
		m_dataIDOfEyPhase = NULL;
	}
	if (m_dataIDOfHxPhase != NULL) {
		delete[] m_dataIDOfHxPhase;
		m_dataIDOfHxPhase = NULL;
	}
	if (m_dataIDOfHyPhase != NULL) {
		delete[] m_dataIDOfHyPhase;
		m_dataIDOfHyPhase = NULL;
	}
	if (m_dataIDOfHzPhase != NULL) {
		delete[] m_dataIDOfHzPhase;
		m_dataIDOfHzPhase = NULL;
	}

	if (m_numOfFreqCalculatedByThisStaAndPE > 0) {

		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
		m_ExCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_EyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_HxCalculatedCSEM = new  std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_HyCalculatedCSEM = new  std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_HzCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];

		m_ExCalculatedAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_EyCalculatedAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HxCalculatedAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HyCalculatedAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HzCalculatedAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_ExCalculatedPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_EyCalculatedPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HxCalculatedPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HyCalculatedPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HzCalculatedPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_ExResidualAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_EyResidualAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HxResidualAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HyResidualAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HzResidualAmplitude = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_ExResidualPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_EyResidualPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HxResidualPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HyResidualPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_HzResidualPhase = new double[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfExAmplitude = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfEyAmplitude = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHxAmplitude = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHyAmplitude = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHzAmplitude = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfExPhase = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfEyPhase = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHxPhase = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHyPhase = new int[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHzPhase = new int[m_numOfFreqCalculatedByThisStaAndPE];

		// Initialize
		for (int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i) {
			m_ExCalculated[i] = std::complex<double>(0.0, 0.0);
			m_EyCalculated[i] = std::complex<double>(0.0, 0.0);
			m_HxCalculatedCSEM[i] = std::complex<double>(0.0, 0.0);
			m_HyCalculatedCSEM[i] = std::complex<double>(0.0, 0.0);
			m_HzCalculated[i] = std::complex<double>(0.0, 0.0);
			m_ExCalculatedAmplitude[i] = 0.0;
			m_EyCalculatedAmplitude[i] = 0.0;
			m_HxCalculatedAmplitude[i] = 0.0;
			m_HyCalculatedAmplitude[i] = 0.0;
			m_HzCalculatedAmplitude[i] = 0.0;
			m_ExCalculatedPhase[i] = 0.0;
			m_EyCalculatedPhase[i] = 0.0;
			m_HxCalculatedPhase[i] = 0.0;
			m_HyCalculatedPhase[i] = 0.0;
			m_HzCalculatedPhase[i] = 0.0;
			m_ExResidualAmplitude[i] = 0.0;
			m_EyResidualAmplitude[i] = 0.0;
			m_HxResidualAmplitude[i] = 0.0;
			m_HyResidualAmplitude[i] = 0.0;
			m_HzResidualAmplitude[i] = 0.0;
			m_ExResidualPhase[i] = 0.0;
			m_EyResidualPhase[i] = 0.0;
			m_HxResidualPhase[i] = 0.0;
			m_HyResidualPhase[i] = 0.0;
			m_HzResidualPhase[i] = 0.0;
			m_dataIDOfExAmplitude[i] = -1;
			m_dataIDOfEyAmplitude[i] = -1;
			m_dataIDOfHxAmplitude[i] = -1;
			m_dataIDOfHyAmplitude[i] = -1;
			m_dataIDOfHzAmplitude[i] = -1;
			m_dataIDOfExPhase[i] = -1;
			m_dataIDOfEyPhase[i] = -1;
			m_dataIDOfHxPhase[i] = -1;
			m_dataIDOfHyPhase[i] = -1;
			m_dataIDOfHzPhase[i] = -1;
		}

	}

}

// Output calculated Impedance tensors
void ObservedDataStationCSEMAmplitudeAndPhase::outputCalculatedValues() const {

	int icount(0);
	for (std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr) {
		fprintf(OutputFiles::m_csvFile, "%10d, %10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_sourceID, m_stationID, m_freq[*itr], m_ExCalculatedAmplitude[icount], m_ExCalculatedPhase[icount], m_EyCalculatedAmplitude[icount], m_EyCalculatedPhase[icount],
			m_HxCalculatedAmplitude[icount], m_HxCalculatedPhase[icount], m_HyCalculatedAmplitude[icount], m_HyCalculatedPhase[icount], m_HzCalculatedAmplitude[icount], m_HzCalculatedPhase[icount]);
		fprintf(OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,	%15.6e, %15.6e,",
			m_ExResidualAmplitude[icount], m_ExResidualPhase[icount], m_EyResidualAmplitude[icount], m_EyResidualPhase[icount],
			m_HxResidualAmplitude[icount], m_HxResidualPhase[icount], m_HyResidualAmplitude[icount], m_HyResidualPhase[icount], m_HzResidualAmplitude[icount], m_HzResidualPhase[icount]);
		fprintf(OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_ExObservedAmplitude[*itr], m_ExObservedPhase[*itr], m_EyObservedAmplitude[*itr], m_EyObservedPhase[*itr],
			m_HxObservedAmplitude[*itr], m_HxObservedPhase[*itr], m_HyObservedAmplitude[*itr], m_HyObservedPhase[*itr], m_HzObservedAmplitude[*itr], m_HzObservedPhase[*itr]);
		fprintf(OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n",
			m_ExObservedAmplitudeSD[*itr], m_ExObservedPhaseSD[*itr], m_EyObservedAmplitudeSD[*itr], m_EyObservedPhaseSD[*itr],
			m_HxObservedAmplitudeSD[*itr], m_HxObservedPhaseSD[*itr], m_HyObservedAmplitudeSD[*itr], m_HyObservedPhaseSD[*itr], m_HzObservedAmplitudeSD[*itr], m_HzObservedPhaseSD[*itr]);
		++icount;
	}
}

// Calulate interpolator vector of electric field
void ObservedDataStationCSEMAmplitudeAndPhase::calcInterpolatorVectorOfElectricField(Forward3D* const ptrForward3D) {

	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	if (meshType == MeshData::TETRA) {// Tetra mesh
		switch (getTypeOfElectricField()) {
		case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
			if (m_useExData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldXDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfEx);
			}
			if (m_useEyData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldYDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfEy);
			}
			break;
		case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
			if (m_useExData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialX(m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2, m_rhsVectorIDOfEx);
			}
			if (m_useEyData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialY(m_elementIncludingStation, m_faceIncludingStation, m_areaCoordinateValues.coord1, m_areaCoordinateValues.coord2, m_rhsVectorIDOfEy);
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
			exit(1);
			break;
		}
	}
	else if (meshType == MeshData::HEXA) {// Hexa mesh
		if (m_useExData) {
			ptrForward3D->calcInterpolatorVectorOfElectricFieldXDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx);
		}
		if (m_useEyData) {
			ptrForward3D->calcInterpolatorVectorOfElectricFieldYDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEy);
		}
	}
	else if (meshType == MeshData::NONCONFORMING_HEXA) {
		switch (getTypeOfElectricField()) {
		case AnalysisControl::USE_HORIZONTAL_ELECTRIC_FIELD:
			if (m_useExData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldXDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx);
			}
			if (m_useEyData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldYDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEy);
			}
			break;
		case AnalysisControl::USE_TANGENTIAL_ELECTRIC_FIELD:
			if (m_useExData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialXFromAllEdges(m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEx);
			}
			if (m_useEyData) {
				ptrForward3D->calcInterpolatorVectorOfElectricFieldTangentialYFromAllEdges(m_elementIncludingStation, m_faceIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfEy);
			}
			break;
		default:
			OutputFiles::m_logFile << "Error : Unknown type of the electric field." << std::endl;
			exit(1);
			break;
		}
	}
	else {
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

void ObservedDataStationCSEMAmplitudeAndPhase::calcInterpolatorVectorOfHorizontalMagneticField(Forward3D* const ptrForward3D) {

	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	if (meshType == MeshData::TETRA) {// Tetra mesh
		if (m_useHxData) {
			ptrForward3D->calcInterpolatorVectorOfMagneticFieldXDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfHx);
		}
		if (m_useHyData) {
			ptrForward3D->calcInterpolatorVectorOfMagneticFieldYDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfHy);
		}
	}
	else if (meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA) {
		if (m_useHxData) {
			ptrForward3D->calcInterpolatorVectorOfMagneticFieldXDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfHx);
		}
		if (m_useHyData) {
			ptrForward3D->calcInterpolatorVectorOfMagneticFieldYDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfHy);
		}
	}
	else {
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Calulate interpolator vector of vertical magnetic field
void ObservedDataStationCSEMAmplitudeAndPhase::calcInterpolatorVectorOfVerticalMagneticField(Forward3D* const ptrForward3D) {

	const int meshType = (AnalysisControl::getInstance())->getTypeOfMesh();
	if (meshType == MeshData::TETRA) {// Tetra mesh
		if (m_useHzData) {
			ptrForward3D->calcInterpolatorVectorOfMagneticFieldZDirection(m_elementIncludingStation, m_volumeCoordinateValues.coord1, m_volumeCoordinateValues.coord2, m_volumeCoordinateValues.coord3, m_rhsVectorIDOfHz);
		}
	}
	else if (meshType == MeshData::HEXA || meshType == MeshData::NONCONFORMING_HEXA) {
		if (m_useHzData) {
			ptrForward3D->calcInterpolatorVectorOfMagneticFieldZDirection(m_elementIncludingStation, m_localCoordinateValues.X, m_localCoordinateValues.Y, m_localCoordinateValues.Z, m_rhsVectorIDOfHz);
		}
	}
	else {
		OutputFiles::m_logFile << "Error : Wrong type of mesh : " << meshType << std::endl;
		exit(1);
	}

}

// Calulate sensitivity matrix of Impedance tensors
void ObservedDataStationCSEMAmplitudeAndPhase::calculateSensitivityMatrix(const double freq, const int nModel,
	const ObservedDataStationPoint* const ptrStationOfMagneticField,
	const std::complex<double>* const derivativesOfEMField,
	double* const sensitivityMatrix, const bool forceSDToOne) const {

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE(freq);

	if (freqIDThisPEInSta < 0) {// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	//sourceID
	int iSource = m_iSource;

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[freqIDThisPEInSta];

	const std::complex<double> Ex = m_ExCalculated[freqIDThisPEInSta];
	const std::complex<double> Ey = m_EyCalculated[freqIDThisPEInSta];
	const std::complex<double> Hx = m_HxCalculatedCSEM[freqIDThisPEInSta];
	const std::complex<double> Hy = m_HyCalculatedCSEM[freqIDThisPEInSta];
	const std::complex<double> Hz = m_HzCalculated[freqIDThisPEInSta];

	const long long rhsVectorIDOfHx = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHx());
	const long long rhsVectorIDOfHy = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHy());
	const long long nBlkNotFixed = static_cast<long long>((ResistivityBlock::getInstance())->getNumResistivityBlockNotFixed());
	const double ln10 = log(10.0);
	const double epsE = std::max(std::norm(Ex), std::norm(Ey)) * 1.0e-20;
	const double epsH = std::max({ std::norm(Hx), std::norm(Hy), std::norm(Hz) }) * 1.0e-20;

	for (long long imdl = 0; imdl < nBlkNotFixed; ++imdl) {

		const double dExRe = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEx + imdl].real();
		const double dEyRe = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEy + imdl].real();
		const double dHxRe = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHx + imdl].real();
		const double dHyRe = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHy + imdl].real();
		const double dHzRe = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHz + imdl].real();
		const double dExIm = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEx + imdl].imag();
		const double dEyIm = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEy + imdl].imag();
		const double dHxIm = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHx + imdl].imag();
		const double dHyIm = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHy + imdl].imag();
		const double dHzIm = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHz + imdl].imag();

		if (m_useExData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfExAmplitude[freqIDThisPEInSta]) + imdl]
				= 1.0 / ln10 / std::max(std::norm(Ex), epsE) * (Ex.real() * dExRe + Ex.imag() * dExIm) / m_ExObservedAmplitudeSD[freqIDGlobalInSta];
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfExPhase[freqIDThisPEInSta]) + imdl]
				= CommonParameters::rad2deg / std::max(std::norm(Ex), epsE) * (Ex.real() * dExIm - Ex.imag() * dExRe) / m_ExObservedPhaseSD[freqIDGlobalInSta];
		}
		if (m_useEyData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfEyAmplitude[freqIDThisPEInSta]) + imdl]
				= 1.0 / ln10 / std::max(std::norm(Ey), epsE) * (Ey.real() * dEyRe + Ey.imag() * dEyIm) / m_EyObservedAmplitudeSD[freqIDGlobalInSta];
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfEyPhase[freqIDThisPEInSta]) + imdl]
				= CommonParameters::rad2deg / std::max(std::norm(Ey), epsE) * (Ey.real() * dEyIm - Ey.imag() * dEyRe) / m_EyObservedPhaseSD[freqIDGlobalInSta];
		}
		if (m_useHxData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHxAmplitude[freqIDThisPEInSta]) + imdl]
				= 1.0 / ln10 / std::max(std::norm(Hx), epsH) * (Hx.real() * dHxRe + Hx.imag() * dHxIm) / m_HxObservedAmplitudeSD[freqIDGlobalInSta];
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHxPhase[freqIDThisPEInSta]) + imdl]
				= CommonParameters::rad2deg / std::max(std::norm(Hx), epsH) * (Hx.real() * dHxIm - Hx.imag() * dHxRe) / m_HxObservedPhaseSD[freqIDGlobalInSta];
		}
		if (m_useHyData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHyAmplitude[freqIDThisPEInSta]) + imdl]
				= 1.0 / ln10 / std::max(std::norm(Hy), epsH) * (Hy.real() * dHyRe + Hy.imag() * dHyIm) / m_HyObservedAmplitudeSD[freqIDGlobalInSta];
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHyPhase[freqIDThisPEInSta]) + imdl]
				= CommonParameters::rad2deg / std::max(std::norm(Hy), epsH) * (Hy.real() * dHyIm - Hy.imag() * dHyRe) / m_HyObservedPhaseSD[freqIDGlobalInSta];
		}
		if (m_useHzData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHzAmplitude[freqIDThisPEInSta]) + imdl]
				= 1.0 / ln10 / std::max(std::norm(Hz), epsH) * (Hz.real() * dHzRe + Hz.imag() * dHzIm) / m_HzObservedAmplitudeSD[freqIDGlobalInSta];
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHzPhase[freqIDThisPEInSta]) + imdl]
				= CommonParameters::rad2deg / std::max(std::norm(Hz), epsH) * (Hz.real() * dHzIm - Hz.imag() * dHzRe) / m_HzObservedPhaseSD[freqIDGlobalInSta];
		}
	}

	//const std::complex<double> HxCalculated[2] = { ptrStationOfMagneticField->getHxCalculated(0), ptrStationOfMagneticField->getHxCalculated(1) };
	//const std::complex<double> HyCalculated[2] = { ptrStationOfMagneticField->getHyCalculated(0), ptrStationOfMagneticField->getHyCalculated(1) };

	//const std::complex<double> det = HxCalculated[0] * HyCalculated[1] - HxCalculated[1] * HyCalculated[0];
	//const std::complex<double> divDet = std::complex<double>(1.0, 0.0) / det;
	//const std::complex<double> divDet2 = divDet * divDet;

	//const long long rhsVectorIDOfHx = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHx());
	//const long long rhsVectorIDOfHy = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHy());
	//const long long nBlkNotFixed = static_cast<long long>((ResistivityBlock::getInstance())->getNumResistivityBlockNotFixed());

	////	//----- debug >>>>>
	////#ifdef _DEBUG_WRITE
	////	std::cout << "divDet " << divDet << std::endl;
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Ex imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*m_rhsVectorIDOfEx + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Ey imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*m_rhsVectorIDOfEy + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Hx imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*rhsVectorIDOfHx   + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Hy imdl derivatives[0] " << imdl << " " << derivativesOfEMFieldExPol[ nBlkNotFixed*rhsVectorIDOfHy   + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Ex imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*m_rhsVectorIDOfEx + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Ey imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*m_rhsVectorIDOfEy + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Hx imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*rhsVectorIDOfHx   + imdl ] << std::endl;
	////	}
	////	for( int imdl = 0; imdl < nBlkNotFixed; ++imdl ){
	////		std::cout << "Hy imdl derivatives[1] " << imdl << " " << derivativesOfEMFieldEyPol[ nBlkNotFixed*rhsVectorIDOfHy   + imdl ] << std::endl;
	////	}
	////#endif
	////	//----- debug <<<<<

	//	//const double baseOfStaticShift = ( ObservedData::getInstance() )->getBaseOfStaticShift();

	//	//const double factor( pow( 10.0, m_staticShiftFactor ) );

	//double cxx = 1.0;
	//double cxy = 0.0;
	//double cyx = 0.0;
	//double cyy = 1.0;
	//if (!m_fixDistortionMatrix) { // Distortion matrix is not fixed
	//	const int typeOfDistortion = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//	if (typeOfDistortion == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//		assert(m_arrayDistortionMatrixDifferences != NULL);
	//		cxx = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXX] + 1.0;
	//		cxy = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CXY];
	//		cyx = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYX];
	//		cyy = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[ObservedDataStationMT::COMPONENT_ID_CYY] + 1.0;
	//	}
	//	else if (typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//		assert(m_arrayGainsAndRotations != NULL);
	//		const double gX = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN]);
	//		const double gY = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN]);
	//		const double betaX = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_ROTATION];
	//		const double betaY = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_ROTATION];
	//		cxx = gX * cos(betaX);
	//		cxy = -gY * sin(betaY);
	//		cyx = gX * sin(betaX);
	//		cyy = gY * cos(betaY);
	//	}
	//	else if (typeOfDistortion == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//		assert(m_arrayGainsAndRotations != NULL);
	//		const double gX = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN]);
	//		const double gY = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN]);
	//		cxx = gX;
	//		cxy = 0.0;
	//		cyx = 0.0;
	//		cyy = gY;
	//	}
	//}

	//CommonParameters::DoubleComplexValues ZxxSD = m_ZxxSD[freqIDGlobalInSta];
	//CommonParameters::DoubleComplexValues ZxySD = m_ZxySD[freqIDGlobalInSta];
	//CommonParameters::DoubleComplexValues ZyxSD = m_ZyxSD[freqIDGlobalInSta];
	//CommonParameters::DoubleComplexValues ZyySD = m_ZyySD[freqIDGlobalInSta];
	//if (forceSDToOne) {
	//	// force erro to one
	//	ZxxSD.realPart = 1.0;
	//	ZxySD.realPart = 1.0;
	//	ZyxSD.realPart = 1.0;
	//	ZyySD.realPart = 1.0;
	//	ZxxSD.imagPart = 1.0;
	//	ZxySD.imagPart = 1.0;
	//	ZyxSD.imagPart = 1.0;
	//	ZyySD.imagPart = 1.0;
	//}

	//for (long long imdl = 0; imdl < nBlkNotFixed; ++imdl) {

	//	const std::complex<double> work1 = derivativesOfEMFieldExPol[nBlkNotFixed * rhsVectorIDOfHx + imdl] * HyCalculated[1]
	//		+ derivativesOfEMFieldEyPol[nBlkNotFixed * rhsVectorIDOfHy + imdl] * HxCalculated[0]
	//		- derivativesOfEMFieldEyPol[nBlkNotFixed * rhsVectorIDOfHx + imdl] * HyCalculated[0]
	//		- derivativesOfEMFieldExPol[nBlkNotFixed * rhsVectorIDOfHy + imdl] * HxCalculated[1];

	//	// dZxx/dm
	//	const std::complex<double> workXX1 = derivativesOfEMFieldExPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl] * HyCalculated[1]
	//		+ derivativesOfEMFieldEyPol[nBlkNotFixed * rhsVectorIDOfHy + imdl] * m_ExCalculated[0]
	//		- derivativesOfEMFieldEyPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl] * HyCalculated[0]
	//		- derivativesOfEMFieldExPol[nBlkNotFixed * rhsVectorIDOfHy + imdl] * m_ExCalculated[1];

	//	const std::complex<double> workXX2 = m_ExCalculated[0] * HyCalculated[1] - m_ExCalculated[1] * HyCalculated[0];

	//	const double dZxxRealUndist = (workXX1 * divDet - work1 * workXX2 * divDet2).real();
	//	const double dZxxImagUndist = (workXX1 * divDet - work1 * workXX2 * divDet2).imag();

	//	// dZxy/dm
	//	const std::complex<double> workXY1 = derivativesOfEMFieldEyPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl] * HxCalculated[0]
	//		+ derivativesOfEMFieldExPol[nBlkNotFixed * rhsVectorIDOfHx + imdl] * m_ExCalculated[1]
	//		- derivativesOfEMFieldExPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEx) + imdl] * HxCalculated[1]
	//		- derivativesOfEMFieldEyPol[nBlkNotFixed * rhsVectorIDOfHx + imdl] * m_ExCalculated[0];

	//	const std::complex<double> workXY2 = m_ExCalculated[1] * HxCalculated[0] - m_ExCalculated[0] * HxCalculated[1];

	//	const double dZxyRealUndist = (workXY1 * divDet - work1 * workXY2 * divDet2).real();
	//	const double dZxyImagUndist = (workXY1 * divDet - work1 * workXY2 * divDet2).imag();

	//	// dZyx/dm
	//	const std::complex<double> workYX1 = derivativesOfEMFieldExPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl] * HyCalculated[1]
	//		+ derivativesOfEMFieldEyPol[nBlkNotFixed * rhsVectorIDOfHy + imdl] * m_EyCalculated[0]
	//		- derivativesOfEMFieldEyPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl] * HyCalculated[0]
	//		- derivativesOfEMFieldExPol[nBlkNotFixed * rhsVectorIDOfHy + imdl] * m_EyCalculated[1];

	//	const std::complex<double> workYX2 = m_EyCalculated[0] * HyCalculated[1] - m_EyCalculated[1] * HyCalculated[0];

	//	const double dZyxRealUndist = (workYX1 * divDet - work1 * workYX2 * divDet2).real();
	//	const double dZyxImagUndist = (workYX1 * divDet - work1 * workYX2 * divDet2).imag();

	//	// dZyy/dm
	//	const std::complex<double> workYY1 = derivativesOfEMFieldEyPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl] * HxCalculated[0]
	//		+ derivativesOfEMFieldExPol[nBlkNotFixed * rhsVectorIDOfHx + imdl] * m_EyCalculated[1]
	//		- derivativesOfEMFieldExPol[nBlkNotFixed * static_cast<long long>(m_rhsVectorIDOfEy) + imdl] * HxCalculated[1]
	//		- derivativesOfEMFieldEyPol[nBlkNotFixed * rhsVectorIDOfHx + imdl] * m_EyCalculated[0];

	//	const std::complex<double> workYY2 = m_EyCalculated[1] * HxCalculated[0] - m_EyCalculated[0] * HxCalculated[1];

	//	const double dZyyRealUndist = (workYY1 * divDet - work1 * workYY2 * divDet2).real();
	//	const double dZyyImagUndist = (workYY1 * divDet - work1 * workYY2 * divDet2).imag();

	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + imdl] = (cxx * dZxxRealUndist + cxy * dZyxRealUndist) / ZxxSD.realPart;
	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + imdl] = (cxx * dZxyRealUndist + cxy * dZyyRealUndist) / ZxySD.realPart;
	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + imdl] = (cyx * dZxxRealUndist + cyy * dZyxRealUndist) / ZyxSD.realPart;
	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + imdl] = (cyx * dZxyRealUndist + cyy * dZyyRealUndist) / ZyySD.realPart;

	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + imdl] = (cxx * dZxxImagUndist + cxy * dZyxImagUndist) / ZxxSD.imagPart;
	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + imdl] = (cxx * dZxyImagUndist + cxy * dZyyImagUndist) / ZxySD.imagPart;
	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + imdl] = (cyx * dZxxImagUndist + cyy * dZyxImagUndist) / ZyxSD.imagPart;
	//	sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + imdl] = (cyx * dZxyImagUndist + cyy * dZyyImagUndist) / ZyySD.imagPart;

	//	if (!m_fixDistortionMatrix) { // Distortion matrix is not fixed

	//		const std::complex<double> ZxxUndist = (m_ExCalculated[0] * HyCalculated[1] - m_ExCalculated[1] * HyCalculated[0]) / det;
	//		const std::complex<double> ZxyUndist = (m_ExCalculated[1] * HxCalculated[0] - m_ExCalculated[0] * HxCalculated[1]) / det;
	//		const std::complex<double> ZyxUndist = (m_EyCalculated[0] * HyCalculated[1] - m_EyCalculated[1] * HyCalculated[0]) / det;
	//		const std::complex<double> ZyyUndist = (m_EyCalculated[1] * HxCalculated[0] - m_EyCalculated[0] * HxCalculated[1]) / det;

	//		if ((AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//			const long long ID_Cxx = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CXX]);
	//			const long long ID_Cxy = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CXY]);
	//			const long long ID_Cyx = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CYX]);
	//			const long long ID_Cyy = static_cast<long long>(m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[COMPONENT_ID_CYY]);
	//			if (ID_Cxx >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx] = ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx] = ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx] = ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx] = ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx] = 0.0;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx] = 0.0;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxx] = 0.0;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxx] = 0.0;// For Im(Zyy)
	//			}
	//			if (ID_Cxy >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy] = ZyxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy] = ZyxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy] = ZyyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy] = ZyyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy] = 0.0;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy] = 0.0;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cxy] = 0.0;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cxy] = 0.0;// For Im(Zyy)
	//			}
	//			if (ID_Cyx >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx] = 0.0;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx] = 0.0;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx] = 0.0;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx] = 0.0;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx] = ZxxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx] = ZxxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyx] = ZxyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyx] = ZxyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//			if (ID_Cyy >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy] = 0.0;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy] = 0.0;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy] = 0.0;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy] = 0.0;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy] = ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy] = ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_Cyy] = ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_Cyy] = ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//		}
	//		else if ((AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {

	//			const double ln10 = log(10.0);

	//			const double gX = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN]);
	//			const double gY = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN]);
	//			const double betaX = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_ROTATION];
	//			const double betaY = m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_ROTATION];

	//			const long long ID_GainX = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_GAIN]);
	//			const long long ID_GainY = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_GAIN]);
	//			const long long ID_RotX = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_ROTATION]);
	//			const long long ID_RotY = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_ROTATION]);

	//			if (ID_GainX >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * cos(betaX) * ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * cos(betaX) * ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * cos(betaX) * ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * cos(betaX) * ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * sin(betaX) * ZxxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * sin(betaX) * ZxxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * sin(betaX) * ZxyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * sin(betaX) * ZxyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//			if (ID_GainY >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = -gY * ln10 * sin(betaY) * ZyxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = -gY * ln10 * sin(betaY) * ZyxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = -gY * ln10 * sin(betaY) * ZyyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = -gY * ln10 * sin(betaY) * ZyyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * cos(betaY) * ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * cos(betaY) * ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * cos(betaY) * ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * cos(betaY) * ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//			if (ID_RotX >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX] = -gX * sin(betaX) * ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX] = -gX * sin(betaX) * ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX] = -gX * sin(betaX) * ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX] = -gX * sin(betaX) * ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX] = gX * cos(betaX) * ZxxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX] = gX * cos(betaX) * ZxxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotX] = gX * cos(betaX) * ZxyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotX] = gX * cos(betaX) * ZxyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//			if (ID_RotY >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY] = -gY * cos(betaY) * ZyxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY] = -gY * cos(betaY) * ZyxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY] = -gY * cos(betaY) * ZyyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY] = -gY * cos(betaY) * ZyyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY] = -gY * sin(betaY) * ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY] = -gY * sin(betaY) * ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_RotY] = -gY * sin(betaY) * ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_RotY] = -gY * sin(betaY) * ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//		}
	//		else if ((AnalysisControl::getInstance())->getTypeOfDistortion() == AnalysisControl::ESTIMATE_GAINS_ONLY) {

	//			const double ln10 = log(10.0);

	//			const double gX = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EX_GAIN]);
	//			const double gY = pow(10.0, m_arrayGainsAndRotations->gainsAndRotations[ObservedDataStationMT::EY_GAIN]);

	//			const long long ID_GainX = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_GAIN]);
	//			const long long ID_GainY = static_cast<long long>(m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_GAIN]);
	//			if (ID_GainX >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * ZxxUndist.real() / ZxxSD.realPart;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * ZxxUndist.imag() / ZxxSD.imagPart;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * ZxyUndist.real() / ZxySD.realPart;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = gX * ln10 * ZxyUndist.imag() / ZxySD.imagPart;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = 0.0;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = 0.0;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainX] = 0.0;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainX] = 0.0;// For Im(Zyy)
	//			}
	//			if (ID_GainY >= 0) {
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = 0.0;// For Re(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = 0.0;// For Im(Zxx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = 0.0;// For Re(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZxy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = 0.0;// For Im(Zxy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * ZyxUndist.real() / ZyxSD.realPart;// For Re(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyx[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * ZyxUndist.imag() / ZyxSD.imagPart;// For Im(Zyx)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].realPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * ZyyUndist.real() / ZyySD.realPart;// For Re(Zyy)
	//				sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfZyy[freqIDThisPEInSta].imagPart) + nBlkNotFixed + ID_GainY] = gY * ln10 * ZyyUndist.imag() / ZyySD.imagPart;// For Im(Zyy)
	//			}
	//		}
	//	}
	//}

}

// Calculate data vector of this PE
void ObservedDataStationCSEMAmplitudeAndPhase::calculateResidualVectorOfDataThisPE(const double freq, const int offset, double* vector) const {

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE(freq);

	if (freqIDThisPEInSta < 0) {// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	//?
	if (m_useExData) {
		vector[offset + m_dataIDOfExAmplitude[freqIDThisPEInSta]] = m_ExResidualAmplitude[freqIDThisPEInSta];
		vector[offset + m_dataIDOfExPhase[freqIDThisPEInSta]] = m_ExResidualPhase[freqIDThisPEInSta];
	}
	if (m_useEyData) {
		vector[offset + m_dataIDOfEyAmplitude[freqIDThisPEInSta]] = m_EyResidualAmplitude[freqIDThisPEInSta];
		vector[offset + m_dataIDOfEyPhase[freqIDThisPEInSta]] = m_EyResidualPhase[freqIDThisPEInSta];
	}
	if (m_useHxData) {
		vector[offset + m_dataIDOfHxAmplitude[freqIDThisPEInSta]] = m_HxResidualAmplitude[freqIDThisPEInSta];
		vector[offset + m_dataIDOfHxPhase[freqIDThisPEInSta]] = m_HxResidualPhase[freqIDThisPEInSta];
	}
	if (m_useHyData) {
		vector[offset + m_dataIDOfHyAmplitude[freqIDThisPEInSta]] = m_HyResidualAmplitude[freqIDThisPEInSta];
		vector[offset + m_dataIDOfHyPhase[freqIDThisPEInSta]] = m_HyResidualPhase[freqIDThisPEInSta];
	}
	if (m_useHzData) {
		vector[offset + m_dataIDOfHzAmplitude[freqIDThisPEInSta]] = m_HzResidualAmplitude[freqIDThisPEInSta];
		vector[offset + m_dataIDOfHzPhase[freqIDThisPEInSta]] = m_HzResidualPhase[freqIDThisPEInSta];
	}

}

// Calulate L2 norm of misfit
double ObservedDataStationCSEMAmplitudeAndPhase::calculateErrorSumOfSquaresThisPE() const {

	double misfit(0.0);

	for (int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq) {
		if (m_useExData) {
			misfit += pow(m_ExResidualAmplitude[ifreq], 2);
			misfit += pow(m_ExResidualPhase[ifreq], 2);
		}
		if (m_useEyData) {
			misfit += pow(m_EyResidualAmplitude[ifreq], 2);
			misfit += pow(m_EyResidualPhase[ifreq], 2);
		}
		if (m_useHxData) {
			misfit += pow(m_HxResidualAmplitude[ifreq], 2);
			misfit += pow(m_HxResidualPhase[ifreq], 2);
		}
		if (m_useHyData) {
			misfit += pow(m_HyResidualAmplitude[ifreq], 2);
			misfit += pow(m_HyResidualPhase[ifreq], 2);
		}
		if (m_useHzData) {
			misfit += pow(m_HzResidualAmplitude[ifreq], 2);
			misfit += pow(m_HzResidualPhase[ifreq], 2);
		}
	}

	return misfit;

}

// Copy current distortion parameters to previous ones
void ObservedDataStationCSEMAmplitudeAndPhase::copyDistortionParamsCurToPre(const int iComp) {

	//const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);
	//	assert(iComp >= COMPONENT_ID_CXX);
	//	assert(iComp <= COMPONENT_ID_CYY);
	//	m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp] = m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp];
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp >= EX_GAIN);
	//	assert(iComp <= EY_ROTATION);
	//	m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = m_arrayGainsAndRotations->gainsAndRotations[iComp];
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = m_arrayGainsAndRotations->gainsAndRotations[iComp];
	//}
	//else {
	//	OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
	//	exit(1);
	//}

}

/// Get flag specifing whether distortion matrix are fixed or not
bool ObservedDataStationCSEMAmplitudeAndPhase::doesFixDistortionMatrix() const {

	return m_fixDistortionMatrix;

}

// Get type of the electric field used to calculate response functions
int ObservedDataStationCSEMAmplitudeAndPhase::getTypeOfElectricField() const {

	return m_typeOfElectricField;

}

// Set flag specifing whether distortion matrix are fixed or not
void ObservedDataStationCSEMAmplitudeAndPhase::setFixDistortionMatrix(const bool doesFix) {

	m_fixDistortionMatrix = doesFix;

}

// Set type of the electric field used to calculate response functions
void ObservedDataStationCSEMAmplitudeAndPhase::setTypeOfElectricField(const int type) {

	m_typeOfElectricField = type;

}

// Set distortion parameters of previous iteration
void ObservedDataStationCSEMAmplitudeAndPhase::setDistortionParamsPre(const int iComp, const double val) {

	/*const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
		assert(m_arrayDistortionMatrixDifferences != NULL);
		assert(iComp >= COMPONENT_ID_CXX);
		assert(iComp <= COMPONENT_ID_CYY);
		m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp] = val;
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp >= EX_GAIN);
		assert(iComp <= EY_ROTATION);
		m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = val;
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp == EX_GAIN || iComp == EY_GAIN);
		m_arrayGainsAndRotations->gainsAndRotationsPre[iComp] = val;
	}
	else {
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}*/


}

// Set distortion parameters
void ObservedDataStationCSEMAmplitudeAndPhase::setDistortionParams(const int iComp, const double val) {

	//const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);
	//	assert(iComp >= COMPONENT_ID_CXX);
	//	assert(iComp <= COMPONENT_ID_CYY);
	//	m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp] = val;
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp >= EX_GAIN);
	//	assert(iComp <= EY_ROTATION);
	//	m_arrayGainsAndRotations->gainsAndRotations[iComp] = val;
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp == EX_GAIN || iComp == EY_GAIN);
	//	m_arrayGainsAndRotations->gainsAndRotations[iComp] = val;
	//}
	//else {
	//	OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
	//	exit(1);
	//}

}

// Set ID of distortion parameters
void ObservedDataStationCSEMAmplitudeAndPhase::setIDOfDistortionParams(const int iComp, const int ID) {

	//const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);
	//	assert(iComp >= COMPONENT_ID_CXX);
	//	assert(iComp <= COMPONENT_ID_CYY);
	//	if (doesFixDistortionMatrix()) {
	//		m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp] = -1;
	//		OutputFiles::m_logFile << "Warning : ID of distortion matrix must be -1 if distortion matrix is fixed." << std::endl;
	//	}
	//	else {// Distortion matrix is not fixed
	//		m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp] = ID;
	//	}

	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp >= EX_GAIN);
	//	assert(iComp <= EY_ROTATION);
	//	if (doesFixDistortionMatrix()) {
	//		m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = -1;
	//		OutputFiles::m_logFile << "Warning : ID of distortion matrix must be -1 if distortion matrix is fixed." << std::endl;
	//	}
	//	else {// Distortion matrix is not fixed
	//		m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = ID;
	//	}
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp == EX_GAIN || iComp == EY_GAIN);
	//	if (doesFixDistortionMatrix()) {
	//		m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = -1;
	//		OutputFiles::m_logFile << "Warning : ID of distortion matrix must be -1 if distortion matrix is fixed." << std::endl;
	//	}
	//	else {// Distortion matrix is not fixed
	//		m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] = ID;
	//	}
	//}
	//else {
	//	OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
	//	exit(1);
	//}

}

// Set full updated value of distortion parameters
void ObservedDataStationCSEMAmplitudeAndPhase::setDistortionParamsUpdatedFull(const int iComp, const double val) {

	//const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);
	//	assert(iComp >= COMPONENT_ID_CXX);
	//	assert(iComp <= COMPONENT_ID_CYY);
	//	if (doesFixDistortionMatrix()) {// Not change distortion matrix
	//		m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp] = m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp];
	//	}
	//	else {
	//		m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp] = val;
	//	}
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp >= EX_GAIN);
	//	assert(iComp <= EY_ROTATION);
	//	if (doesFixDistortionMatrix()) {// Not change distortion matrix
	//		m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	//	}
	//	else {
	//		if (iComp == ObservedDataStationMT::EX_ROTATION || iComp == ObservedDataStationMT::EY_ROTATION) {
	//			// Rotation angle is bounded in from -pi/2 to pi/2 radians
	//			if (val > 0.5 * CommonParameters::PI) {
	//				m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = 0.5 * CommonParameters::PI;
	//			}
	//			else if (val < -0.5 * CommonParameters::PI) {
	//				m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = -0.5 * CommonParameters::PI;
	//			}
	//			else {
	//				m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = val;
	//			}
	//		}
	//		else {
	//			m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = val;
	//		}
	//	}

	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp == EX_GAIN || iComp == EY_GAIN);
	//	if (doesFixDistortionMatrix()) {// Not change distortion matrix
	//		m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	//	}
	//	else {
	//		m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] = val;
	//	}
	//}
	//else {
	//	OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
	//	exit(1);
	//}

}

// Update distortion parameters
void ObservedDataStationCSEMAmplitudeAndPhase::updateDistortionParams(const double dampingFactor) {

	//if (doesFixDistortionMatrix()) {// distortion matrix if fixed
	//	return;
	//}

	//const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();

	//if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);

	//	for (int iComp = 0; iComp < 4; ++iComp) {
	//		if (m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp] >= 0) {
	//			m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp] =
	//				dampingFactor * m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp] + (1.0 - dampingFactor) * m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp];
	//		}
	//	}

	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);

	//	for (int iComp = 0; iComp < 4; ++iComp) {
	//		if (m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp] >= 0) {
	//			if (iComp == ObservedDataStationMT::EX_ROTATION || iComp == ObservedDataStationMT::EY_ROTATION) {
	//				// Rotation angle is bounded in from -pi/2 to pi/2 radians
	//				// y = tan(beta)
	//				// beta = Arctan(y)
	//				// dy/dbeta = 1 / cos(beta)^2
	//				const double yUpdatedFull = tan(m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp]);
	//				const double yPre = tan(m_arrayGainsAndRotations->gainsAndRotationsPre[iComp]);
	//				const double yNew = dampingFactor * yUpdatedFull + (1.0 - dampingFactor) * yPre;
	//				const double betaNew = atan(yNew);
	//				m_arrayGainsAndRotations->gainsAndRotations[iComp] = betaNew;
	//			}
	//			else {
	//				m_arrayGainsAndRotations->gainsAndRotations[iComp] =
	//					dampingFactor * m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp] + (1.0 - dampingFactor) * m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	//			}
	//		}
	//	}

	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	if (m_arrayGainsAndRotations->IDsOfGainsAndRotations[EX_GAIN] >= 0) {
	//		m_arrayGainsAndRotations->gainsAndRotations[EX_GAIN] =
	//			dampingFactor * m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[EX_GAIN] + (1.0 - dampingFactor) * m_arrayGainsAndRotations->gainsAndRotationsPre[EX_GAIN];
	//	}
	//	if (m_arrayGainsAndRotations->IDsOfGainsAndRotations[EY_GAIN] >= 0) {
	//		m_arrayGainsAndRotations->gainsAndRotations[EY_GAIN] =
	//			dampingFactor * m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[EY_GAIN] + (1.0 - dampingFactor) * m_arrayGainsAndRotations->gainsAndRotationsPre[EY_GAIN];
	//	}
	//}
	//else {
	//	OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
	//	exit(1);
	//}

}

// Get distortion parameters of previous iteration
double ObservedDataStationCSEMAmplitudeAndPhase::getDistortionParamsPre(const int iComp) const {

	//const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	//if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
	//	assert(m_arrayDistortionMatrixDifferences != NULL);
	//	assert(iComp >= COMPONENT_ID_CXX);
	//	assert(iComp <= COMPONENT_ID_CYY);

	//	return m_arrayDistortionMatrixDifferences->distortionMatrixDifferencePre[iComp];
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp >= EX_GAIN);
	//	assert(iComp <= EY_ROTATION);

	//	return m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	//}
	//else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
	//	assert(m_arrayGainsAndRotations != NULL);
	//	assert(iComp == EX_GAIN || iComp == EY_GAIN);

	//	return m_arrayGainsAndRotations->gainsAndRotationsPre[iComp];
	//}
	//else {
	//	OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
	//	exit(1);
	//}

}

// Get distortion parameters
double ObservedDataStationCSEMAmplitudeAndPhase::getDistortionParams(const int iComp) const {

	/*const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
		assert(m_arrayDistortionMatrixDifferences != NULL);
		assert(iComp >= COMPONENT_ID_CXX);
		assert(iComp <= COMPONENT_ID_CYY);

		return m_arrayDistortionMatrixDifferences->distortionMatrixDifference[iComp];
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp >= EX_GAIN);
		assert(iComp <= EY_ROTATION);

		return m_arrayGainsAndRotations->gainsAndRotations[iComp];
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp == EX_GAIN || iComp == EY_GAIN);

		return m_arrayGainsAndRotations->gainsAndRotations[iComp];
	}
	else {
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}*/

}

// Get ID of distortion parameters
int ObservedDataStationCSEMAmplitudeAndPhase::getIDOfDistortionParams(const int iComp) const {

	/*const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
		assert(m_arrayDistortionMatrixDifferences != NULL);
		assert(iComp >= COMPONENT_ID_CXX);
		assert(iComp <= COMPONENT_ID_CYY);

		return m_arrayDistortionMatrixDifferences->IDsOfDistortionMatrixDifference[iComp];
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp >= EX_GAIN);
		assert(iComp <= EY_ROTATION);

		return m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp];
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
		assert(m_arrayGainsAndRotations != NULL);

		return m_arrayGainsAndRotations->IDsOfGainsAndRotations[iComp];
	}
	else {
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}*/

}

// Get full updated value of distortion parameters
double ObservedDataStationCSEMAmplitudeAndPhase::getDistortionParamsUpdatedFull(const int iComp) const {

	/*const int type = (AnalysisControl::getInstance())->getTypeOfDistortion();
	if (type == AnalysisControl::ESTIMATE_DISTORTION_MATRIX_DIFFERENCE) {
		assert(m_arrayDistortionMatrixDifferences != NULL);
		assert(iComp >= COMPONENT_ID_CXX);
		assert(iComp <= COMPONENT_ID_CYY);
		assert(iComp >= COMPONENT_ID_CXX);
		assert(iComp <= COMPONENT_ID_CYY);

		return m_arrayDistortionMatrixDifferences->distortionMatrixDifferenceUpdatedFull[iComp];
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_AND_ROTATIONS) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp >= EX_GAIN);
		assert(iComp <= EY_ROTATION);

		return m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp];
	}
	else if (type == AnalysisControl::ESTIMATE_GAINS_ONLY) {
		assert(m_arrayGainsAndRotations != NULL);
		assert(iComp == EX_GAIN || iComp == EY_GAIN);

		return m_arrayGainsAndRotations->gainsAndRotationsUpdatedFull[iComp];
	}
	else {
		OutputFiles::m_logFile << "Error : Type of galvanic distortion is wrong. : " << type << std::endl;
		exit(1);
	}*/

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isDataTypeUsed(int* dataType, int target) {

	for (int i = 0; i < m_nTypeOfData; i++) {
		if (dataType[i] == target) {
			return true;
		}
	}

	return false;

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isDataTypeDuplicated(int* dataType) {

	for (int i = 0; i < m_nTypeOfData; i++) {
		for (int j = 0; j < m_nTypeOfData; j++) {
			if ((dataType[i] == dataType[j]) && (i != j)) {
				return true;
			}
		}
	}

	return false;

}

void ObservedDataStationCSEMAmplitudeAndPhase::initializeSourceID(const int sourceID, const int iSource) {

	m_sourceID = sourceID;
	m_iSource = iSource;

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isThisStationUseExData() const {

	if (m_useExData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isThisStationUseEyData() const {

	if (m_useEyData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isThisStationUseHxData() const {

	if (m_useHxData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isThisStationUseHyData() const {

	if (m_useHyData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEMAmplitudeAndPhase::isThisStationUseHzData() const {

	if (m_useHzData) {
		return true;
	}
	else {
		return false;
	}

}

int ObservedDataStationCSEMAmplitudeAndPhase::getiSourceOfThisStation() const {

	return m_iSource;

}
