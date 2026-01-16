//-------------------------------------------------------------------------------------------------------
// This module provides CSEM (Controlled Source Electromagnetics) modeling capabilities within FEMTIC.
// Implemented by Zuwei Huang, 2025.
//
// For support, bug reports, or questions, please contact the author at: hzw1498218560@tongji.edu.cn
//-------------------------------------------------------------------------------------------------------
#pragma once
#ifndef DBLDEF_OBSERVED_DATA_SOURCE_AMD
#define DBLDEF_OBSERVED_DATA_SOURCE_AMD

#include <vector>
#include <complex>

#include "ObservedDataSource.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of MT station
class ObservedDataSourceAMD : public ObservedDataSource {
public:

	// Constructer
	explicit ObservedDataSourceAMD();

	// Destructer
	virtual ~ObservedDataSourceAMD();

	// Read data from input file
	virtual void inputObservedData(std::ifstream& inFile);

protected:
	
	// Define the azimuth angle in the direction of due north
	double *azimuth;
	
	// Define the inclination angle with respect to the horizontal direction
	double *dip;

	int *m_sourceDirection;

private:

	// Copy constructer
	ObservedDataSourceAMD(const ObservedDataSourceAMD& rhs);

	// Copy assignment operator
	ObservedDataSourceAMD& operator=(const ObservedDataSourceAMD& rhs);

	// Calculate direction vector 
	void calculateDirecctionVector(int nED);

	void findElementIncludingDipoleSource(int iDipole);

};

#endif
