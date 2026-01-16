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

#include "ObservedDataStationCSEM.h"
#include "ObservedDataStationPoint.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

// Constructer
ObservedDataStationCSEM::ObservedDataStationCSEM() :
	ObservedDataStationPoint(),
	m_dataTypeToBeUsed(NULL),
	m_ExObserved(NULL),
	m_EyObserved(NULL),
	m_HxObserved(NULL),
	m_HyObserved(NULL),
	m_HzObserved(NULL),
	m_ExSD(NULL),
	m_EySD(NULL),
	m_HxSD(NULL),
	m_HySD(NULL),
	m_HzSD(NULL),
	m_ExCalculated(NULL),
	m_EyCalculated(NULL),
	m_HxCalculatedCSEM(NULL),
	m_HyCalculatedCSEM(NULL),
	m_HzCalculated(NULL),
	m_ExResidual(NULL),
	m_EyResidual(NULL),
	m_HxResidual(NULL),
	m_HyResidual(NULL),
	m_HzResidual(NULL),
	m_rhsVectorIDOfEx(0),
	m_rhsVectorIDOfEy(0),
	m_nTypeOfData(0),
	m_dataIDOfEx(NULL),
	m_dataIDOfEy(NULL),
	m_dataIDOfHx(NULL),
	m_dataIDOfHy(NULL),
	m_dataIDOfHz(NULL),
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
ObservedDataStationCSEM::~ObservedDataStationCSEM() {

	if (m_ExObserved != NULL) {
		delete[] m_ExObserved;
		m_ExObserved = NULL;
	}

	if (m_EyObserved != NULL) {
		delete[] m_EyObserved;
		m_EyObserved = NULL;
	}

	if (m_HxObserved != NULL) {
		delete[] m_HxObserved;
		m_HxObserved = NULL;
	}

	if (m_HyObserved != NULL) {
		delete[] m_HyObserved;
		m_HyObserved = NULL;
	}

	if (m_HzObserved != NULL) {
		delete[] m_HzObserved;
		m_HzObserved = NULL;
	}

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

	if (m_ExSD != NULL) {
		delete[] m_ExSD;
		m_ExSD = NULL;
	}

	if (m_EySD != NULL) {
		delete[] m_EySD;
		m_EySD = NULL;
	}

	if (m_HxSD != NULL) {
		delete[] m_HxSD;
		m_HxSD = NULL;
	}

	if (m_HySD != NULL) {
		delete[] m_HySD;
		m_HySD = NULL;
	}

	if (m_HzSD != NULL) {
		delete[] m_HzSD;
		m_HzSD = NULL;
	}

	if (m_ExResidual != NULL) {
		delete[] m_ExResidual;
		m_ExResidual = NULL;
	}

	if (m_EyResidual != NULL) {
		delete[] m_EyResidual;
		m_EyResidual = NULL;
	}

	if (m_HxResidual != NULL) {
		delete[] m_HxResidual;
		m_HxResidual = NULL;
	}

	if (m_HyResidual != NULL) {
		delete[] m_HyResidual;
		m_HyResidual = NULL;
	}

	if (m_HzResidual != NULL) {
		delete[] m_HzResidual;
		m_HzResidual = NULL;
	}

	if (m_dataIDOfEx != NULL) {
		delete[] m_dataIDOfEx;
		m_dataIDOfEx = NULL;
	}

	if (m_dataIDOfEy != NULL) {
		delete[] m_dataIDOfEy;
		m_dataIDOfEy = NULL;
	}

	if (m_dataIDOfHx != NULL) {
		delete[] m_dataIDOfHx;
		m_dataIDOfHx = NULL;
	}

	if (m_dataIDOfHy != NULL) {
		delete[] m_dataIDOfHy;
		m_dataIDOfHy = NULL;
	}

	if (m_dataIDOfHz != NULL) {
		delete[] m_dataIDOfHz;
		m_dataIDOfHz = NULL;
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
void ObservedDataStationCSEM::inputObservedData(std::ifstream& inFile) {

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
	if (isDataTypeUsed(m_dataTypeToBeUsed, EX_COMPLEX)) {
		m_useExData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, EY_COMPLEX)) {
		m_useEyData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, HX_COMPLEX)) {
		m_useHxData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, HY_COMPLEX)) {
		m_useHyData = true;
	}
	if (isDataTypeUsed(m_dataTypeToBeUsed, HZ_COMPLEX)) {
		m_useHzData = true;
	}

	//read frequency
	inFile >> m_numOfFrequency;

	const int nFreq = m_numOfFrequency;

	if (nFreq > 0) {
		m_freq = new double[nFreq];
		//	if (m_useExData) {
		m_ExObserved = new std::complex<double>[nFreq];
		m_ExSD = new CommonParameters::DoubleComplexValues[nFreq];
		//		}
		//		if (m_useEyData) {
		m_EyObserved = new std::complex<double>[nFreq];
		m_EySD = new CommonParameters::DoubleComplexValues[nFreq];
		//		}
		//		if (m_useHxData) {
		m_HxObserved = new std::complex<double>[nFreq];
		m_HxSD = new CommonParameters::DoubleComplexValues[nFreq];
		//		}
		//		if (m_useHyData) {
		m_HyObserved = new std::complex<double>[nFreq];
		m_HySD = new CommonParameters::DoubleComplexValues[nFreq];
		//		}
		//		if (m_useHzData) {
		m_HzObserved = new std::complex<double>[nFreq];
		m_HzSD = new CommonParameters::DoubleComplexValues[nFreq];
		//		}
		for (int i = 0; i < nFreq; ++i) {
			m_ExObserved[i] = complexDefaultValue;
			m_ExSD[i].realPart = defaultValue;
			m_ExSD[i].imagPart = defaultValue;
			m_EyObserved[i] = complexDefaultValue;
			m_EySD[i].realPart = defaultValue;
			m_EySD[i].imagPart = defaultValue;
			m_HxObserved[i] = complexDefaultValue;
			m_HxSD[i].realPart = defaultValue;
			m_HxSD[i].imagPart = defaultValue;
			m_HyObserved[i] = complexDefaultValue;
			m_HySD[i].realPart = defaultValue;
			m_HySD[i].imagPart = defaultValue;
			m_HzObserved[i] = complexDefaultValue;
			m_HzSD[i].realPart = defaultValue;
			m_HzSD[i].imagPart = defaultValue;
		}

		for (int i = 0; i < nFreq; ++i) {
			inFile >> m_freq[i];
			//read data
			for (int iD = 0; iD < m_nTypeOfData; iD++) {
				double dbuf1(0.0);
				double dbuf2(0.0);
				inFile >> dbuf1 >> dbuf2;
				switch (m_dataTypeToBeUsed[iD]) {
				case EX_COMPLEX:
					m_ExObserved[i] = std::complex<double>(dbuf1, dbuf2);
					break;
				case EY_COMPLEX:
					m_EyObserved[i] = std::complex<double>(dbuf1, dbuf2);
					break;
				case HX_COMPLEX:
					m_HxObserved[i] = std::complex<double>(dbuf1, dbuf2);
					break;
				case HY_COMPLEX:
					m_HyObserved[i] = std::complex<double>(dbuf1, dbuf2);
					break;
				case HZ_COMPLEX:
					m_HzObserved[i] = std::complex<double>(dbuf1, dbuf2);
					break;
				}
			}
			//read error
			for (int iD = 0; iD < m_nTypeOfData; iD++) {
				switch (m_dataTypeToBeUsed[iD]) {
				case EX_COMPLEX:
					inFile >> m_ExSD[i].realPart;
					inFile >> m_ExSD[i].imagPart;
					break;
				case EY_COMPLEX:
					inFile >> m_EySD[i].realPart;
					inFile >> m_EySD[i].imagPart;
					break;
				case HX_COMPLEX:
					inFile >> m_HxSD[i].realPart;
					inFile >> m_HxSD[i].imagPart;
					break;
				case HY_COMPLEX:
					inFile >> m_HySD[i].realPart;
					inFile >> m_HySD[i].imagPart;
					break;
				case HZ_COMPLEX:
					inFile >> m_HzSD[i].realPart;
					inFile >> m_HzSD[i].imagPart;
					break;
				}
			}
		}
	}

#ifdef _DEBUG_WRITE
	std::cout << " MT " << m_stationID << " " << m_IDOfMagneticFieldStation << std::endl;
	std::cout << m_location.X << " " << m_location.Y << std::endl;
	std::cout << m_numOfFrequency << std::endl;
	for (int i = 0; i < m_numOfFrequency; ++i) {
		std::cout << m_freq[i] << " "
			<< m_ExObserved[i] << " "
			<< m_EyObserved[i] << " "
			<< m_HxObserved[i] << " "
			<< m_HyObserved[i] << " "
			<< m_HzObserved[i] << " "
			<< m_ExSD[i].realPart << " "
			<< m_ExSD[i].imagPart << " "
			<< m_EySD[i].realPart << " "
			<< m_EySD[i].imagPart << " "
			<< m_HxSD[i].realPart << " "
			<< m_HxSD[i].imagPart << " "
			<< m_HySD[i].realPart << " "
			<< m_HySD[i].imagPart << " "
			<< m_HzSD[i].realPart << " "
			<< m_HzSD[i].imagPart << std::endl;
	}
#endif

}

// Calulate electric field
void ObservedDataStationCSEM::calculateElectricField(const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy, const int ifre) {

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
void ObservedDataStationCSEM::calculateVerticalMagneticField(const Forward3D* const ptrForward3D, const int rhsVectorIDOfHz, const int ifre) {

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
void ObservedDataStationCSEM::calculateDataTensor(const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount) {

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE(freq);

	if (freqIDThisPEInSta < 0) {// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	const int freqIDGlobalInSta = m_freqIDsAmongThisStationCalculatedByThisPE[freqIDThisPEInSta];

	m_HxCalculatedCSEM[freqIDThisPEInSta] = ptrStationOfMagneticField->getHxCalculated(0);
	m_HyCalculatedCSEM[freqIDThisPEInSta] = ptrStationOfMagneticField->getHyCalculated(0);

	if (m_useExData) {
		m_ExResidual[freqIDThisPEInSta].realPart = (m_ExObserved[freqIDGlobalInSta].real() - m_ExCalculated[freqIDThisPEInSta].real()) / m_ExSD[freqIDGlobalInSta].realPart;
		m_ExResidual[freqIDThisPEInSta].imagPart = (m_ExObserved[freqIDGlobalInSta].imag() - m_ExCalculated[freqIDThisPEInSta].imag()) / m_ExSD[freqIDGlobalInSta].imagPart;
		m_dataIDOfEx[freqIDThisPEInSta].realPart = icount++;
		m_dataIDOfEx[freqIDThisPEInSta].imagPart = icount++;
	}
	if (m_useEyData) {
		m_EyResidual[freqIDThisPEInSta].realPart = (m_EyObserved[freqIDGlobalInSta].real() - m_EyCalculated[freqIDThisPEInSta].real()) / m_EySD[freqIDGlobalInSta].realPart;
		m_EyResidual[freqIDThisPEInSta].imagPart = (m_EyObserved[freqIDGlobalInSta].imag() - m_EyCalculated[freqIDThisPEInSta].imag()) / m_EySD[freqIDGlobalInSta].imagPart;
		m_dataIDOfEy[freqIDThisPEInSta].realPart = icount++;
		m_dataIDOfEy[freqIDThisPEInSta].imagPart = icount++;
	}
	if (m_useHxData) {
		m_HxResidual[freqIDThisPEInSta].realPart = (m_HxObserved[freqIDGlobalInSta].real() - m_HxCalculatedCSEM[freqIDThisPEInSta].real()) / m_HxSD[freqIDGlobalInSta].realPart;
		m_HxResidual[freqIDThisPEInSta].imagPart = (m_HxObserved[freqIDGlobalInSta].imag() - m_HxCalculatedCSEM[freqIDThisPEInSta].imag()) / m_HxSD[freqIDGlobalInSta].imagPart;
		m_dataIDOfHx[freqIDThisPEInSta].realPart = icount++;
		m_dataIDOfHx[freqIDThisPEInSta].imagPart = icount++;
	}
	if (m_useHyData) {
		m_HyResidual[freqIDThisPEInSta].realPart = (m_HyObserved[freqIDGlobalInSta].real() - m_HyCalculatedCSEM[freqIDThisPEInSta].real()) / m_HySD[freqIDGlobalInSta].realPart;
		m_HyResidual[freqIDThisPEInSta].imagPart = (m_HyObserved[freqIDGlobalInSta].imag() - m_HyCalculatedCSEM[freqIDThisPEInSta].imag()) / m_HySD[freqIDGlobalInSta].imagPart;
		m_dataIDOfHy[freqIDThisPEInSta].realPart = icount++;
		m_dataIDOfHy[freqIDThisPEInSta].imagPart = icount++;
	}
	if (m_useHzData) {
		m_HzResidual[freqIDThisPEInSta].realPart = (m_HzObserved[freqIDGlobalInSta].real() - m_HzCalculated[freqIDThisPEInSta].real()) / m_HzSD[freqIDGlobalInSta].realPart;
		m_HzResidual[freqIDThisPEInSta].imagPart = (m_HzObserved[freqIDGlobalInSta].imag() - m_HzCalculated[freqIDThisPEInSta].imag()) / m_HzSD[freqIDGlobalInSta].imagPart;
		m_dataIDOfHz[freqIDThisPEInSta].realPart = icount++;
		m_dataIDOfHz[freqIDThisPEInSta].imagPart = icount++;
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
void ObservedDataStationCSEM::initializeElectricField(const int freq) {

	m_ExCalculated[freq] = std::complex<double>(0.0, 0.0);
	m_EyCalculated[freq] = std::complex<double>(0.0, 0.0);

}

// Initialize vertical magnetic field
void ObservedDataStationCSEM::initializeVerticalMagneticField(const int freq) {

	m_HzCalculated[freq] = std::complex<double>(0.0, 0.0);

}

// Initialize Impedance tensor and errors
void ObservedDataStationCSEM::initializeDataTensorsAndErrors() {

	//for( int i = 0; i < m_numOfFrequency; ++i ){
	for (int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i) {
		m_ExResidual[i].realPart = 0.0;
		m_ExResidual[i].imagPart = 0.0;
		m_EyResidual[i].realPart = 0.0;
		m_EyResidual[i].imagPart = 0.0;
		m_HxResidual[i].realPart = 0.0;
		m_HxResidual[i].imagPart = 0.0;
		m_HyResidual[i].realPart = 0.0;
		m_HyResidual[i].imagPart = 0.0;
		m_HzResidual[i].imagPart = 0.0;
		m_HzResidual[i].imagPart = 0.0;

		m_ExCalculated[i] = std::complex<double>(0.0, 0.0);
		m_EyCalculated[i] = std::complex<double>(0.0, 0.0);
		m_HxCalculatedCSEM[i] = std::complex<double>(0.0, 0.0);
		m_HyCalculatedCSEM[i] = std::complex<double>(0.0, 0.0);
		m_HzCalculated[i] = std::complex<double>(0.0, 0.0);

	}

}

// Allocate memory for the calculated Impedance and errors
void ObservedDataStationCSEM::allocateMemoryForCalculatedValues() {

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

	if (m_HzResidual != NULL) {
		delete[] m_HzResidual;
		m_HzResidual = NULL;
	}

	if (m_dataIDOfEx != NULL) {
		delete[] m_dataIDOfEx;
		m_dataIDOfEx = NULL;
	}

	if (m_dataIDOfEy != NULL) {
		delete[] m_dataIDOfEy;
		m_dataIDOfEy = NULL;
	}

	if (m_dataIDOfHx != NULL) {
		delete[] m_dataIDOfHx;
		m_dataIDOfHx = NULL;
	}

	if (m_dataIDOfHy != NULL) {
		delete[] m_dataIDOfHy;
		m_dataIDOfHy = NULL;
	}

	if (m_dataIDOfHz != NULL) {
		delete[] m_dataIDOfHz;
		m_dataIDOfHz = NULL;
	}

	if (m_numOfFreqCalculatedByThisStaAndPE > 0) {

		// If total number of frequencies calculated by this PE is not more than zero, do not allocate memory
	//	if (m_useExData) {
		m_ExCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_ExResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfEx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		//	}
		//	if (m_useEyData) {
		m_EyCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_EyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfEy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		//	}
		//	if (m_useHxData) {
		m_HxCalculatedCSEM = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_HxResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHx = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		//	}
		//	if (m_useHyData) {
		m_HyCalculatedCSEM = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_HyResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHy = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		//	}
		//	if (m_useHzData) {
		m_HzCalculated = new std::complex<double>[m_numOfFreqCalculatedByThisStaAndPE];
		m_HzResidual = new CommonParameters::DoubleComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		m_dataIDOfHz = new CommonParameters::InitComplexValues[m_numOfFreqCalculatedByThisStaAndPE];
		//	}
			// Initialize
		for (int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i) {
			m_ExCalculated[i] = std::complex<double>(0.0, 0.0);
			m_EyCalculated[i] = std::complex<double>(0.0, 0.0);
			m_HxCalculatedCSEM[i] = std::complex<double>(0.0, 0.0);
			m_HyCalculatedCSEM[i] = std::complex<double>(0.0, 0.0);
			m_HzCalculated[i] = std::complex<double>(0.0, 0.0);

			m_ExResidual[i].realPart = 0.0;
			m_ExResidual[i].imagPart = 0.0;
			m_EyResidual[i].realPart = 0.0;
			m_EyResidual[i].imagPart = 0.0;
			m_HxResidual[i].realPart = 0.0;
			m_HxResidual[i].imagPart = 0.0;
			m_HyResidual[i].realPart = 0.0;
			m_HyResidual[i].imagPart = 0.0;
			m_HzResidual[i].realPart = 0.0;
			m_HzResidual[i].imagPart = 0.0;
			m_dataIDOfEx[i].realPart = -1;
			m_dataIDOfEx[i].imagPart = -1;
			m_dataIDOfEy[i].realPart = -1;
			m_dataIDOfEy[i].imagPart = -1;
			m_dataIDOfHx[i].realPart = -1;
			m_dataIDOfHx[i].imagPart = -1;
			m_dataIDOfHy[i].realPart = -1;
			m_dataIDOfHy[i].imagPart = -1;
			m_dataIDOfHz[i].realPart = -1;
			m_dataIDOfHz[i].imagPart = -1;
		}

	}

}

// Output calculated Impedance tensors
void ObservedDataStationCSEM::outputCalculatedValues() const {

	int icount(0);
	for (std::vector<int>::const_iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr) {
		fprintf(OutputFiles::m_csvFile, "%10d, %10d, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_sourceID, m_stationID, m_freq[*itr], m_ExCalculated[icount].real(), m_ExCalculated[icount].imag(), m_EyCalculated[icount].real(), m_EyCalculated[icount].imag(),
			m_HxCalculatedCSEM[icount].real(), m_HxCalculatedCSEM[icount].imag(), m_HyCalculatedCSEM[icount].real(), m_HyCalculatedCSEM[icount].imag(), m_HzCalculated[icount].real(), m_HzCalculated[icount].imag());
		fprintf(OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,	%15.6e, %15.6e,",
			m_ExResidual[icount].realPart, m_ExResidual[icount].imagPart, m_EyResidual[icount].realPart, m_EyResidual[icount].imagPart,
			m_HxResidual[icount].realPart, m_HxResidual[icount].imagPart, m_HyResidual[icount].realPart, m_HyResidual[icount].imagPart, m_HzResidual[icount].realPart, m_HzResidual[icount].imagPart);
		fprintf(OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,",
			m_ExObserved[*itr].real(), m_ExObserved[*itr].imag(), m_EyObserved[*itr].real(), m_EyObserved[*itr].imag(),
			m_HxObserved[*itr].real(), m_HxObserved[*itr].imag(), m_HyObserved[*itr].real(), m_HyObserved[*itr].imag(), m_HzObserved[*itr].real(), m_HzObserved[*itr].imag());
		fprintf(OutputFiles::m_csvFile, " %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e, %15.6e,\n",
			m_ExSD[*itr].realPart, m_ExSD[*itr].imagPart, m_EySD[*itr].realPart, m_EySD[*itr].imagPart,
			m_HxSD[*itr].realPart, m_HxSD[*itr].imagPart, m_HySD[*itr].realPart, m_HySD[*itr].imagPart, m_HzSD[*itr].realPart, m_HzSD[*itr].imagPart);
		++icount;
	}
}

// Calulate interpolator vector of electric field
void ObservedDataStationCSEM::calcInterpolatorVectorOfElectricField(Forward3D* const ptrForward3D) {

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

void ObservedDataStationCSEM::calcInterpolatorVectorOfHorizontalMagneticField(Forward3D* const ptrForward3D) {

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
void ObservedDataStationCSEM::calcInterpolatorVectorOfVerticalMagneticField(Forward3D* const ptrForward3D) {

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
void ObservedDataStationCSEM::calculateSensitivityMatrix(const double freq, const int nModel,
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

	const long long rhsVectorIDOfHx = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHx());
	const long long rhsVectorIDOfHy = static_cast<long long>(ptrStationOfMagneticField->getRhsVectorIDOfHy());
	const long long nBlkNotFixed = static_cast<long long>((ResistivityBlock::getInstance())->getNumResistivityBlockNotFixed());

	for (long long imdl = 0; imdl < nBlkNotFixed; ++imdl) {
		if (m_useExData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfEx[freqIDThisPEInSta].realPart) + imdl] = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEx + imdl].real() / m_ExSD[freqIDGlobalInSta].realPart;
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfEx[freqIDThisPEInSta].imagPart) + imdl] = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEx + imdl].imag() / m_ExSD[freqIDGlobalInSta].imagPart;
		}
		if (m_useEyData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfEy[freqIDThisPEInSta].realPart) + imdl] = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEy + imdl].real() / m_EySD[freqIDGlobalInSta].realPart;
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfEy[freqIDThisPEInSta].imagPart) + imdl] = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfEy + imdl].imag() / m_EySD[freqIDGlobalInSta].imagPart;
		}
		if (m_useHxData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHx[freqIDThisPEInSta].realPart) + imdl] = derivativesOfEMField[nBlkNotFixed * rhsVectorIDOfHx + imdl].real() / m_HxSD[freqIDGlobalInSta].realPart;
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHx[freqIDThisPEInSta].imagPart) + imdl] = derivativesOfEMField[nBlkNotFixed * rhsVectorIDOfHx + imdl].imag() / m_HxSD[freqIDGlobalInSta].imagPart;
		}
		if (m_useHyData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHy[freqIDThisPEInSta].realPart) + imdl] = derivativesOfEMField[nBlkNotFixed * rhsVectorIDOfHy + imdl].real() / m_HySD[freqIDGlobalInSta].realPart;
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHy[freqIDThisPEInSta].imagPart) + imdl] = derivativesOfEMField[nBlkNotFixed * rhsVectorIDOfHy + imdl].imag() / m_HySD[freqIDGlobalInSta].imagPart;
		}
		if (m_useHzData) {
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHz[freqIDThisPEInSta].realPart) + imdl] = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHz + imdl].real() / m_HzSD[freqIDGlobalInSta].realPart;
			sensitivityMatrix[static_cast<long long>(nModel) * static_cast<long long>(m_dataIDOfHz[freqIDThisPEInSta].imagPart) + imdl] = derivativesOfEMField[nBlkNotFixed * m_rhsVectorIDOfHz + imdl].imag() / m_HzSD[freqIDGlobalInSta].imagPart;
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
void ObservedDataStationCSEM::calculateResidualVectorOfDataThisPE(const double freq, const int offset, double* vector) const {

	const int freqIDThisPEInSta = getFreqIDsAmongThisPE(freq);

	if (freqIDThisPEInSta < 0) {// Specified frequency is not the ones calculated by this PE in this station 
		return;
	}

	//?
	if (m_useExData) {
		vector[offset + m_dataIDOfEx[freqIDThisPEInSta].realPart] = m_ExResidual[freqIDThisPEInSta].realPart;
		vector[offset + m_dataIDOfEx[freqIDThisPEInSta].imagPart] = m_ExResidual[freqIDThisPEInSta].imagPart;
	}
	if (m_useEyData) {
		vector[offset + m_dataIDOfEy[freqIDThisPEInSta].realPart] = m_EyResidual[freqIDThisPEInSta].realPart;
		vector[offset + m_dataIDOfEy[freqIDThisPEInSta].imagPart] = m_EyResidual[freqIDThisPEInSta].imagPart;
	}
	if (m_useHxData) {
		vector[offset + m_dataIDOfHx[freqIDThisPEInSta].realPart] = m_HxResidual[freqIDThisPEInSta].realPart;
		vector[offset + m_dataIDOfHx[freqIDThisPEInSta].imagPart] = m_HxResidual[freqIDThisPEInSta].imagPart;
	}
	if (m_useHyData) {
		vector[offset + m_dataIDOfHy[freqIDThisPEInSta].realPart] = m_HyResidual[freqIDThisPEInSta].realPart;
		vector[offset + m_dataIDOfHy[freqIDThisPEInSta].imagPart] = m_HyResidual[freqIDThisPEInSta].imagPart;
	}
	if (m_useHzData) {
		vector[offset + m_dataIDOfHz[freqIDThisPEInSta].realPart] = m_HzResidual[freqIDThisPEInSta].realPart;
		vector[offset + m_dataIDOfHz[freqIDThisPEInSta].imagPart] = m_HzResidual[freqIDThisPEInSta].imagPart;
	}

}

// Calulate L2 norm of misfit
double ObservedDataStationCSEM::calculateErrorSumOfSquaresThisPE() const {

	double misfit(0.0);

	for (int ifreq = 0; ifreq < m_numOfFreqCalculatedByThisStaAndPE; ++ifreq) {
		if (m_useExData) {
			misfit += pow(m_ExResidual[ifreq].realPart, 2);
			misfit += pow(m_ExResidual[ifreq].imagPart, 2);
		}
		if (m_useEyData) {
			misfit += pow(m_EyResidual[ifreq].realPart, 2);
			misfit += pow(m_EyResidual[ifreq].imagPart, 2);
		}
		if (m_useHxData) {
			misfit += pow(m_HxResidual[ifreq].realPart, 2);
			misfit += pow(m_HxResidual[ifreq].imagPart, 2);
		}
		if (m_useHyData) {
			misfit += pow(m_HyResidual[ifreq].realPart, 2);
			misfit += pow(m_HyResidual[ifreq].imagPart, 2);
		}
		if (m_useHzData) {
			misfit += pow(m_HzResidual[ifreq].imagPart, 2);
			misfit += pow(m_HzResidual[ifreq].imagPart, 2);
		}
	}

	return misfit;

}

// Copy current distortion parameters to previous ones
void ObservedDataStationCSEM::copyDistortionParamsCurToPre(const int iComp) {

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
bool ObservedDataStationCSEM::doesFixDistortionMatrix() const {

	return m_fixDistortionMatrix;

}

// Get type of the electric field used to calculate response functions
int ObservedDataStationCSEM::getTypeOfElectricField() const {

	return m_typeOfElectricField;

}

// Set flag specifing whether distortion matrix are fixed or not
void ObservedDataStationCSEM::setFixDistortionMatrix(const bool doesFix) {

	m_fixDistortionMatrix = doesFix;

}

// Set type of the electric field used to calculate response functions
void ObservedDataStationCSEM::setTypeOfElectricField(const int type) {

	m_typeOfElectricField = type;

}

// Set distortion parameters of previous iteration
void ObservedDataStationCSEM::setDistortionParamsPre(const int iComp, const double val) {

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
void ObservedDataStationCSEM::setDistortionParams(const int iComp, const double val) {

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
void ObservedDataStationCSEM::setIDOfDistortionParams(const int iComp, const int ID) {

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
void ObservedDataStationCSEM::setDistortionParamsUpdatedFull(const int iComp, const double val) {

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
void ObservedDataStationCSEM::updateDistortionParams(const double dampingFactor) {

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
double ObservedDataStationCSEM::getDistortionParamsPre(const int iComp) const {

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
double ObservedDataStationCSEM::getDistortionParams(const int iComp) const {

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
int ObservedDataStationCSEM::getIDOfDistortionParams(const int iComp) const {

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
double ObservedDataStationCSEM::getDistortionParamsUpdatedFull(const int iComp) const {

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

bool ObservedDataStationCSEM::isDataTypeUsed(int* dataType, int target) {

	for (int i = 0; i < m_nTypeOfData; i++) {
		if (dataType[i] == target) {
			return true;
		}
	}

	return false;

}

bool ObservedDataStationCSEM::isDataTypeDuplicated(int* dataType) {

	for (int i = 0; i < m_nTypeOfData; i++) {
		for (int j = 0; j < m_nTypeOfData; j++) {
			if ((dataType[i] == dataType[j]) && (i != j)) {
				return true;
			}
		}
	}

	return false;

}

void ObservedDataStationCSEM::initializeSourceID(const int sourceID, const int iSource) {

	m_sourceID = sourceID;
	m_iSource = iSource;

}

bool ObservedDataStationCSEM::isThisStationUseExData() const {

	if (m_useExData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEM::isThisStationUseEyData() const {

	if (m_useEyData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEM::isThisStationUseHxData() const {

	if (m_useHxData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEM::isThisStationUseHyData() const {

	if (m_useHyData) {
		return true;
	}
	else {
		return false;
	}

}

bool ObservedDataStationCSEM::isThisStationUseHzData() const {

	if (m_useHzData) {
		return true;
	}
	else {
		return false;
	}

}

int ObservedDataStationCSEM::getiSourceOfThisStation() const {

	return m_iSource;

}
