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
#ifndef DBLDEF_OBSERVED_DATA_CSEMAmplitudeandPhase
#define DBLDEF_OBSERVED_DATA_CSEMAmplitudeandPhase

#include <vector>
#include <complex>

#include "ObservedDataStationPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of MT station
class ObservedDataStationCSEMAmplitudeAndPhase : public ObservedDataStationPoint {
public:
	enum {
		EX_AMPPHA = 1,
		EY_AMPPHA = 2,
		HX_AMPPHA = 3,
		HY_AMPPHA = 4,
		HZ_AMPPHA = 5
	};

	// Constructer
	explicit ObservedDataStationCSEMAmplitudeAndPhase();

	// Destructer
	virtual ~ObservedDataStationCSEMAmplitudeAndPhase();

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

	std::complex<double>* m_ExCalculated;
	std::complex<double>* m_EyCalculated;
	std::complex<double>* m_HxCalculatedCSEM;
	std::complex<double>* m_HyCalculatedCSEM;
	std::complex<double>* m_HzCalculated;

	int* m_dataIDOfExAmplitude;
	int* m_dataIDOfEyAmplitude;
	int* m_dataIDOfHxAmplitude;
	int* m_dataIDOfHyAmplitude;
	int* m_dataIDOfHzAmplitude;
	int* m_dataIDOfExPhase;
	int* m_dataIDOfEyPhase;
	int* m_dataIDOfHxPhase;
	int* m_dataIDOfHyPhase;
	int* m_dataIDOfHzPhase;


	//
	double* m_ExCalculatedAmplitude;
	double* m_EyCalculatedAmplitude;
	double* m_HxCalculatedAmplitude;
	double* m_HyCalculatedAmplitude;
	double* m_HzCalculatedAmplitude;
	double* m_ExCalculatedPhase;
	double* m_EyCalculatedPhase;
	double* m_HxCalculatedPhase;
	double* m_HyCalculatedPhase;
	double* m_HzCalculatedPhase;
	double* m_ExObservedAmplitude;
	double* m_EyObservedAmplitude;
	double* m_HxObservedAmplitude;
	double* m_HyObservedAmplitude;
	double* m_HzObservedAmplitude;
	double* m_ExObservedPhase;
	double* m_EyObservedPhase;
	double* m_HxObservedPhase;
	double* m_HyObservedPhase;
	double* m_HzObservedPhase;
	double* m_ExObservedAmplitudeSD;
	double* m_EyObservedAmplitudeSD;
	double* m_HxObservedAmplitudeSD;
	double* m_HyObservedAmplitudeSD;
	double* m_HzObservedAmplitudeSD;
	double* m_ExObservedPhaseSD;
	double* m_EyObservedPhaseSD;
	double* m_HxObservedPhaseSD;
	double* m_HyObservedPhaseSD;
	double* m_HzObservedPhaseSD;
	double* m_ExResidualAmplitude;
	double* m_EyResidualAmplitude;
	double* m_HxResidualAmplitude;
	double* m_HyResidualAmplitude;
	double* m_HzResidualAmplitude;
	double* m_ExResidualPhase;
	double* m_EyResidualPhase;
	double* m_HxResidualPhase;
	double* m_HyResidualPhase;
	double* m_HzResidualPhase;

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
	ObservedDataStationCSEMAmplitudeAndPhase(const ObservedDataStationCSEMAmplitudeAndPhase& rhs);

	// Copy assignment operator
	ObservedDataStationCSEMAmplitudeAndPhase& operator=(const ObservedDataStationCSEMAmplitudeAndPhase& rhs);

};

#endif
