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
#pragma once
#ifndef DBLDEF_OBSERVED_DATA_CSEM
#define DBLDEF_OBSERVED_DATA_CSEM

#include <vector>
#include <complex>

#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of MT station
class ObservedDataStationCSEM : public ObservedDataStationPoint {
public:
	enum {
		EX_COMPLEX = 1,
		EY_COMPLEX = 2,
		HX_COMPLEX = 3,
		HY_COMPLEX = 4,
		HZ_COMPLEX = 5
	};

	// Constructer
	explicit ObservedDataStationCSEM();

	// Destructer
	virtual ~ObservedDataStationCSEM();

	// Read data from input file
	virtual void inputObservedData(std::ifstream& inFile);

	// Calulate electric field
	void calculateElectricField(const Forward3D* const ptrForward3D, const int rhsVectorIDOfEx, const int rhsVectorIDOfEy, const int ifre);

	// Calulate vertical magnetic field
	void calculateVerticalMagneticField(const Forward3D* const ptrForward3D, const int rhsVectorIDOfHz, const int ifre);

	// Calulate Data tensor
	void calculateDataTensor(const double freq, const ObservedDataStationPoint* const ptrStationOfMagneticField, int& icount);

	// Initialize electric field
	void initializeElectricField(int freq);

	// Initialize vertical magnetic field
	void initializeVerticalMagneticField(const int freq);

	// Initialize Impedance tensor and errors
	void initializeDataTensorsAndErrors();

	// Allocate memory for the calculated Impedance tensors and errors
	virtual void allocateMemoryForCalculatedValues();

	// Output calculated Impedance tensors
	virtual void outputCalculatedValues() const;

	// Calulate interpolator vector of electric field
	void calcInterpolatorVectorOfElectricField(Forward3D* const ptrForward3D);

	void calcInterpolatorVectorOfVerticalMagneticField(Forward3D* const ptrForward3D);

	void calcInterpolatorVectorOfHorizontalMagneticField(Forward3D* const ptrForward3D) override;

	// Calulate sensitivity matrix
	virtual void calculateSensitivityMatrix(const double freq, const int nModel,
		const ObservedDataStationPoint* const ptrStationOfMagneticField,
		const std::complex<double>* const derivativesOfEMField,
		double* const sensitivityMatrix, const bool forceSDToOne = false) const;

	// Calculate data vector of this PE
	virtual void calculateResidualVectorOfDataThisPE(const double freq, const int offset, double* vector) const;

	// Calulate sum of square of misfit
	virtual double calculateErrorSumOfSquaresThisPE() const;

	// Copy current distortion parameters to previous ones
	void copyDistortionParamsCurToPre(const int iComp);

	// Get flag specifing whether distortion matrix are fixed or not
	bool doesFixDistortionMatrix() const;

	// Get type of the electric field used to calculate response functions
	int getTypeOfElectricField() const;

	// Set flag specifing whether distortion matrix are fixed or not
	void setFixDistortionMatrix(const bool doesFix);

	// Set type of the electric field used to calculate response functions
	void setTypeOfElectricField(const int type);

	// Set distortion parameters of previous iteration
	void setDistortionParamsPre(const int iComp, const double val);

	// Set distortion parameters 
	void setDistortionParams(const int iComp, const double val);

	// Set ID of distortion parameters
	void setIDOfDistortionParams(const int iComp, const int ID);

	// Set full updated value of distortion parameters 
	void setDistortionParamsUpdatedFull(const int iComp, const double val);

	// Update distortion parameters
	void updateDistortionParams(const double dampingFactor);

	// Get distortion parameters of previous iteration
	double getDistortionParamsPre(const int iComp) const;

	// Get distortion parameters
	double getDistortionParams(const int iComp) const;

	// Get ID of distortion parameters
	int getIDOfDistortionParams(const int iComp) const;

	// Get iSource of station
	int getiSourceOfThisStation() const;

	// Get full updated value of distortion parameters
	double getDistortionParamsUpdatedFull(const int iComp) const;

	void initializeSourceID(const int sourceID, const int iSource);

	bool isThisStationUseExData() const;

	bool isThisStationUseEyData() const;

	bool isThisStationUseHxData() const;

	bool isThisStationUseHyData() const;

	bool isThisStationUseHzData() const;

protected:
	struct DistortionMatrixDifferences {
		// Previous distortion matrix difference
		double distortionMatrixDifferencePre[4];

		// Distortion matrix difference
		double distortionMatrixDifference[4];

		// Distorsion matrix difference obtained by inversion which is the ones fully updated ( damping factor = 1 )
		double distortionMatrixDifferenceUpdatedFull[4];

		// Component ID of distortion matrix difference. This ID is -1 if distortion matrix is fixed or this station has no data
		int IDsOfDistortionMatrixDifference[4];
	};

	struct GainsAndRotations {
		// Previous gains and rotations of distortion matrix
		double gainsAndRotationsPre[4];

		// Gains and rotations of distortion matrix
		double gainsAndRotations[4];

		// Gains and rotations of distortion matrix obtained by inversion which is the ones fully updated ( damping factor = 1 )
		double gainsAndRotationsUpdatedFull[4];

		// Component ID of gains and rotations of distortion matrix. This ID is -1 if distortion matrix is fixed or this station has no data
		int IDsOfGainsAndRotations[4];
	};

	bool m_useExData;
	bool m_useEyData;
	bool m_useHxData;
	bool m_useHyData;
	bool m_useHzData;

	std::complex<double>* m_ExObserved;
	std::complex<double>* m_EyObserved;
	std::complex<double>* m_HxObserved;
	std::complex<double>* m_HyObserved;
	std::complex<double>* m_HzObserved;

	CommonParameters::DoubleComplexValues* m_ExSD;
	CommonParameters::DoubleComplexValues* m_EySD;
	CommonParameters::DoubleComplexValues* m_HxSD;
	CommonParameters::DoubleComplexValues* m_HySD;
	CommonParameters::DoubleComplexValues* m_HzSD;

	std::complex<double>* m_ExCalculated;
	std::complex<double>* m_EyCalculated;
	std::complex<double>* m_HxCalculatedCSEM;
	std::complex<double>* m_HyCalculatedCSEM;
	std::complex<double>* m_HzCalculated;

	//std::complex<double> ExCalculated;
	//std::complex<double> EyCalculated;
	//std::complex<double> HxCalculated;
	//std::complex<double> HyCalculated;
	//std::complex<double> HzCalculated;

	CommonParameters::DoubleComplexValues* m_ExResidual;
	CommonParameters::DoubleComplexValues* m_EyResidual;
	CommonParameters::DoubleComplexValues* m_HxResidual;
	CommonParameters::DoubleComplexValues* m_HyResidual;
	CommonParameters::DoubleComplexValues* m_HzResidual;

	CommonParameters::InitComplexValues* m_dataIDOfEx;
	CommonParameters::InitComplexValues* m_dataIDOfEy;
	CommonParameters::InitComplexValues* m_dataIDOfHx;
	CommonParameters::InitComplexValues* m_dataIDOfHy;
	CommonParameters::InitComplexValues* m_dataIDOfHz;

	// Arrays of distortion matrix differences
	DistortionMatrixDifferences* m_arrayDistortionMatrixDifferences;

	// Arrays of gains and rotations of distortion matrix
	GainsAndRotations* m_arrayGainsAndRotations;

private:

	const double defaultValue = 999999.0;

	const std::complex<double> complexDefaultValue= std::complex<double>(999999.0, 999999.0);

	int m_sourceID, m_iSource;

	int m_nTypeOfData;

	int* m_dataTypeToBeUsed;

	int m_rhsVectorIDOfEx;
	int m_rhsVectorIDOfEy;
	int m_rhsVectorIDOfHz;

	// Flag specifing fix distortion matrix or not
	bool m_fixDistortionMatrix;

	// Type of the electric field used to calculated response functions
	int m_typeOfElectricField;

	bool isDataTypeUsed(int* dataType, int target);

	bool isDataTypeDuplicated(int* dataType);


	// We need to allow the copy in public since there are different source of CSEM need to be considered !!
	// Copy constructer
	ObservedDataStationCSEM(const ObservedDataStationCSEM& rhs);

	// Copy assignment operator
	ObservedDataStationCSEM& operator=(const ObservedDataStationCSEM& rhs);

};

#endif
