//-------------------------------------------------------------------------------------------------------
// The MIT License (MIT)
//
// Copyright (c) 2021 Yoshiya Usui
// 
// Modified by Zuwei Huang (2025): Implemented and integrated the CSEM module into the FEMTIC framework.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------
// Modifications Copyright (c) 2025 Zuwei Huang
// This file is based on original work by Yoshiya Usui, with modifications.
// The modified version is released under the GNU General Public License version 3 (GPLv3).
//
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
//
// Note: This file is dual-licensed:
// 1. Original code: MIT License
// 2. Modifications: GPLv3 License
// Overall distribution is under the terms of GPLv3.
//-------------------------------------------------------------------------------------------------------
#ifndef DBLDEF_OBSERVED_DATA
#define DBLDEF_OBSERVED_DATA

#include <vector>
#include <complex>

#include "ObservedDataSource.h"
#include "ObservedDataSourceAED.h"
#include "ObservedDataSourceAMD.h"
#include "ObservedDataSourceSLS.h"
#include "ObservedDataSourceSTLS.h"
#include "ObservedDataStationPoint.h"
#include "ObservedDataStationMT.h"
#include "ObservedDataStationCSEM.h"
#include "ObservedDataStationApparentResistivityAndPhase.h"
#include "ObservedDataStationHTF.h"
#include "ObservedDataStationVTF.h"
#include "ObservedDataStationPT.h"
#include "ObservedDataStationNMT.h"
#include "ObservedDataStationNMT2.h"
#include "ObservedDataStationNMT2ApparentResistivityAndPhase.h"
#include "AdditinalOutputPoint.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Class of observed data
class ObservedData {

public:
	// Return the the instance of the class
	static ObservedData* getInstance();

	// Type IDs of observed data
	enum typeID {
		MT = 0,
		HTF,
		VTF,
		PT,
		NMT,
		NMT2,
		APP_RES_AND_PHS,
		NMT2_APP_RES_AND_PHS,
		CSEM,
	};

	// Input mesh data from "observe.dat"
	void inputObservedData();

	// Input data of distortion matrix from "distortion_iterX.dat"
	void inputDistortionMatrixData();

	// Get total number of different frquencies
	int getTotalNumberOfDifferenetFrequencies() const;

	// Get frequency value
	double getFrequencyValue(const int num) const;

	// Find element including each station
	void findElementIncludingEachStation();

	// Calulate EM field of all stations
	void calculateEMFieldOfAllStations(const Forward3D* const ptrForward3D, const double freq, const int iSourceOriPol, const int ifreq);

	// Calulate response function of all stations
	void calculateResponseFunctionOfAllStations(const int freqIDAmongThisPE);

	// Initialize response functions and errors of all stations
	void initializeResponseFunctionsAndErrorsOfAllStations();

	// Allocate memory for the calculated values and errors of all stations after setting up frequencies calculated by this PE, at which observed value exists
	void allocateMemoryForCalculatedValuesOfAllStations();

	// Output calculated values of all stations
	void outputCalculatedValuesOfAllStations() const;

	// Calulate right-hand side vectors consisting of interpolator vectors
	void calculateInterpolatorVectors(Forward3D* const ptrForward3D, const double freq,const int iSourceOriPol);

	// Calulate sensitivity matrix
	void calculateSensitivityMatrix(std::complex<double>** const derivatives, const double freq, const int nModel, double* const sensitivityMatrix) const;

	// Calculate data vector of this PE
	void calculateResidualVectorOfDataThisPE(double* const vector) const;

	// Get number of interpolator vectors
	int getNumInterpolatorVectors(const int is) const;

	// Get total number of observed data
	int getNumObservedDataThisPETotal() const;

	// Get number of observed data
	int getNumObservedDataThisPE(const int ifreq) const;

	// Get number of observed data accumulated
	int getNumObservedDataThisPEAccumulated(const int num) const;

	// Increment number of interpolator vectors
	int incrementNumInterpolatorVectors();

	// Multiply inputed matrix by scaled data values and output new vector
	void multiplyMatrixByScaledData(const int nRow, const double* const matIn, const double* const vecIn, double* vecOut) const;

	// Calculate frequencies treated by thie PE
	void calcFrequenciesCalculatedByThisPE();

	// Get Total number of frequencies calculated by this PE
	int getNumOfFrequenciesCalculatedByThisPE() const;

	// Get IDs of Frequencies calculated by this PE
	int getIDsOfFrequenciesCalculatedByThisPE(const int num) const;

	// Get values of Frequencies calculated by this PE
	double getValuesOfFrequenciesCalculatedByThisPE(const int num) const;

	// Calculate error sum of squares
	double calculateErrorSumOfSquaresThisPE() const;

	// Calculate sum of squares of distortion matrix complexity
	double calculateSumSquareOfDistortionMatrixComplexity() const;

	// Calculate sum of squares of distortion matrix gains
	double calculateSumSquareOfDistortionMatrixGains() const;

	// Calculate sum of squares of distortion matrix rotations
	double calculateSumSquareOfDistortionMatrixRotations() const;

	// Get total number of the distortion parameters whose value is not fixed
	int getNumDistortionParamsNotFixed() const;

	// Get total number of the equivalent dipoles
	int getNumOfTotalSources() const;

	// Get element ID of the equivalent dipoles
	CommonParameters::source getSourceParameters(int m_sourceID) const;

	// Collect all equivalent dipoles
	void collectEquivalentDipolesFromAllTypesSource();

	// Get types of distortion parameters whose value is not fixed
	int getTypesOfDistortionParamsNotFixed(const int iParamsNotFixed) const;

	// Copy the distortion parameters not fixed to vector
	void copyDistortionParamsNotFixedToVector(double* vector) const;

	// Copy previous distortion parameters not fixed to vector
	void copyDistortionParamsNotFixedPreToVector(double* vector) const;

	// Copy current distortion parameters to previous values
	void copyDistortionParamsCurToPre() const;

	// Calculate full updated values of distortion parameters
	void calcDistortionParamsUpdatedFullFromIncrements(const double* const distortionParamIncre);

	// Update distortion parameters
	void updateDistortionParams();

	// Output distortion parametersto file
	void outputDistortionParams(const int iterNum) const;

	// Output information of locations of observed stations to vtk file
	void outputLocationsOfObservedStationsToVtk() const;

	// Output induction arrow to vtk file
	void outputInductionArrowToVtk(const int iterNum) const;

	// Find and return frequency IDs from frequency value
	int getFreqIDs(const double freq) const;

	// Find local ID of source and station
	int changeLocalSourceIDArrangeToGlobalSourceIDArrange(int sourceID, int iStation) const;

	// Find local ID of source and station
	int localSourceStationNumToGlobalSourceStationNum() const;

	// Find local ID of source and station
	int getGlobalStationIDFromSourceIDAndLocalID(int iSource, int iStation) const;

	// return current source ID
	int getSourceID(int iSource) const;

	bool isSourceAssignedToThisPE(int sourceID, double freq) const;

private:
	// Constructer
	ObservedData();

	// Destructer
	~ObservedData();

	// Copy constructer
	ObservedData(const ObservedData& rhs);

	// Copy assignment operator
	ObservedData& operator=(const ObservedData& rhs);

	// Number of all stations
	int m_numAllStations;

	// Number of usual CSEM stations
	int* m_numStationsCSEMOfDifferentSource;

	int m_numStationsCSEM;

	int* m_sourceIDOfStationPoint;

	std::vector<int> m_sourceID;

	// Number of usual CSEM stations
	int* m_sourceIDStationNum;

	// Number of usual MT stations
	int m_numStationsMT;

	// Number of usual apparent resistivity and phase stations
	int m_numStationsApparentResistivityAndPhase;

	// Number of HTF stations
	int m_numStationsHTF;

	// Number of VTF stations
	int m_numStationsVTF;

	// Number of PT stations
	int m_numStationsPT;

	// Number of NMT (line) stations
	int m_numStationsNMT;

	// Number of NMT (triangle area) stations
	int m_numStationsNMT2;

	// Number of NMT (triangle area) stations using apparent resistivity and phase stations
	int m_numStationsNMT2ApparentResistivityAndPhase;

	// Number of additional output points
	int m_numAdditinalOutputPoint;

	// Number of interpolator vectors
	int m_numInterpolatorVectors;

	// Number of interpolator vectors
	int* m_numCSEMInterpolatorVectors;

	// Number of Arbitrary Electrical Dipoles
	int m_numArbitraryElectricalDipole;

	// Number of Arbitrary Magnetic Dipoles
	int m_numArbitraryMagneticDipole;

	// Number of segmented line source
	int m_numSegmentedLineSource;

	// Number of segmented line source
	int m_numStaightLineSource;

	int m_numOfTotalSources;

	// Number of all equivalent source
	int m_numOfTotalEquivalentDipoles;

	// Map : key = station ID, value = Types of stations
	std::map<int, int> m_stationID2type;

	// Map : key = station ID, value = IDs among each station type of all stations
	std::map<int, int> m_stationID2IDAmongEachStationType;

	// All frequencies for which observed data exists
	std::vector<double> m_frequencyAll;

	// Total number of frequencies calculated by this PE
	int m_numOfFrequenciesCalculatedByThisPE;

	// IDs of frequencies calculated by this PE
	int* m_IDsOfFrequenciesCalculatedByThisPE;

	// Values of Frequencies calculated by this PE
	double* m_valuesOfFrequenciesCalculatedByThisPE;

	// Observed source of horizontal electrical dipoles
	ObservedDataSourceAED* m_observedSourceAED;

	// Observed source of horizontal electrical dipoles
	ObservedDataSourceAMD* m_observedSourceAMD;

	// Observed source of segmented line source
	ObservedDataSourceSLS* m_observedSourceSLS;

	// Observed source of straight line source
	ObservedDataSourceSTLS* m_observedSourceSTLS;

	// Parameters of controlled source
	CommonParameters::source* m_source;

	// Observed data of CSEM station
	ObservedDataStationCSEM* m_observedStationCSEM;

	// Observed data of MT station
	ObservedDataStationMT* m_observedStationMT;

	// Observed data of apparent resistivity and phase stations
	ObservedDataStationApparentResistivityAndPhase* m_observedStationApparentResistivityAndPhase;

	// Observed data of HTF station
	ObservedDataStationHTF* m_observedStationHTF;

	// Observed data of VTF station
	ObservedDataStationVTF* m_observedStationVTF;

	// Observed data of PT station
	ObservedDataStationPT* m_observedStationPT;

	// Observed data of NMT station ( line )
	ObservedDataStationNMT* m_observedStationNMT;

	// Observed data of NMT station ( triangle area )
	ObservedDataStationNMT2* m_observedStationNMT2;

	// Observed data of NMT station ( triangle area ) using apparent resistivity and phase
	ObservedDataStationNMT2ApparentResistivityAndPhase* m_observedStationNMT2ApparentResistivityAndPhase;

	// Additional output points
	AdditinalOutputPoint* m_additinalOutputPoint;

	// Number of data
	int* m_numObservedDataThisPEAccumulated;

	// Types of distortion parameters whose value is not fixed
	std::vector<int> m_typesOfDistortionParamsNotFixed;

	// Add new frequency to the array of all frequencies for which observed data exists after checking 
	// wheter inputed frequency have not been contained
	void checkAndAddNewFrequency(const double freq);

	// Check whether magnetic field value of the specified station is used to calculate response functions of the specified frequency
	bool isUseMagneticFieldOfTheStation(const double freq, const int staID);

	// Get instance of the station where magnetic field is observed
	ObservedDataStationPoint* getInstanceOfMagneticStation(const int IDOfMagSta) const;

	// Calculate total number of components of distortion parameters whose value is not fixed
	void calcNumDistortionParamsNotFixed();

};

#endif
