//-------------------------------------------------------------------------------------------------------
// This module provides CSEM (Controlled Source Electromagnetics) modeling capabilities within FEMTIC.
// Implemented by Zuwei Huang, 2025.
//
// For support, bug reports, or questions, please contact the author at: hzw1498218560@tongji.edu.cn
//-------------------------------------------------------------------------------------------------------
#pragma once
#ifndef DBLDEF_OBSERVED_DATA_SOURCE_AED
#define DBLDEF_OBSERVED_DATA_SOURCE_AED

#include <vector>
#include <complex>

#include "ObservedDataSource.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of MT station
class ObservedDataSourceAED : public ObservedDataSource {
public:

	// Constructer
	explicit ObservedDataSourceAED();

	// Destructer
	virtual ~ObservedDataSourceAED();

	// Read data from input file
	virtual void inputObservedData(std::ifstream& inFile);

protected:
	
	// Define the azimuth angle in the direction of due north
	double *azimuth;
	
	// Define the inclination angle with respect to the horizontal direction
	double *dip;

	int *m_sourceDirection; //-1 or 1

	bool m_findDipoleSourceOnEarthSurface;

private:

	// Copy constructer
	ObservedDataSourceAED(const ObservedDataSourceAED& rhs);

	// Copy assignment operator
	ObservedDataSourceAED& operator=(const ObservedDataSourceAED& rhs);

	// Calculate direction vector 
	void calculateDirecctionVector(int iED);

	void findElementIncludingDipoleSource(int iDipole);

	void findElementIncludingDipoleSourceOnEarthSurface(int iDipole);

};

#endif
