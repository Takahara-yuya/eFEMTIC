//-------------------------------------------------------------------------------------------------------
// This module provides CSEM (Controlled Source Electromagnetics) modeling capabilities within FEMTIC.
// Implemented by Zuwei Huang, 2025.
//
// For support, bug reports, or questions, please contact the author at: hzw1498218560@tongji.edu.cn
//-------------------------------------------------------------------------------------------------------
#pragma once
#ifndef DBLDEF_OBSERVED_DATA_SOURCE_SLS
#define DBLDEF_OBSERVED_DATA_SOURCE_SLS

#include <vector>
#include <complex>

#include "ObservedDataSource.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of MT station
class ObservedDataSourceSLS : public ObservedDataSource {
public:

	// Constructer
	explicit ObservedDataSourceSLS();

	// Destructer
	virtual ~ObservedDataSourceSLS();

	// Read data from input file
	virtual void inputObservedData(std::ifstream& inFile);

	double getTotalLengthOfEquivalentDipoles();

protected:

	// Define the azimuth angle in the direction of due north
	double x_start, y_start, z_start;
	double x_end, y_end, z_end;

	int m_sourceDirection;

	double currentIntensity;

	std::vector<double> eqdx, eqdy, eqdz;
	std::vector<int> eqdie;

	double directionX, directionY, directionZ;

	double lengthAfterpartioning;

private:

	void calculateDirecctionVector();

	void partionTheCurrentLine();

	// Copy constructer
	ObservedDataSourceSLS(const ObservedDataSourceSLS& rhs);

	// Copy assignment operator
	ObservedDataSourceSLS& operator=(const ObservedDataSourceSLS& rhs);

};

#endif
