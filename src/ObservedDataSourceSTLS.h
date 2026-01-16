//-------------------------------------------------------------------------------------------------------
// This module provides CSEM (Controlled Source Electromagnetics) modeling capabilities within FEMTIC.
// Implemented by Zuwei Huang, 2025.
//
// For support, bug reports, or questions, please contact the author at: hzw1498218560@tongji.edu.cn
//-------------------------------------------------------------------------------------------------------
#pragma once
#ifndef DBLDEF_OBSERVED_DATA_SOURCE_STLS
#define DBLDEF_OBSERVED_DATA_SOURCE_STLS

#include <vector>
#include <complex>

#include "ObservedDataSource.h"
#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of station
class ObservedDataSourceSTLS : public ObservedDataSource {
public:

	// Constructer
	explicit ObservedDataSourceSTLS();

	// Destructer
	virtual ~ObservedDataSourceSTLS();

	// Read data from input file
	virtual void inputObservedData(std::ifstream& inFile);

	double getTotalLengthOfEquivalentDipoles();

protected:

	// Define the azimuth angle in the direction of due north
	double *x_start, *y_start, *z_start;
	double *x_end, *y_end, *z_end;

	int *m_sourceDirection;

	int m_numOfSegments;

	double *currentIntensity;

	//std::vector<double> eqdx, eqdy, eqdz;
	//std::vector<int> eqdie;

	double *directionX, *directionY, *directionZ;

	double lengthAfterPartitioning;

	bool m_findCurrentLineOnEarthSurface;

private:

	void calculateDirecctionVector(int iSeg);

	void partitionTheCurrentLine(int iSeg);

	void findEquivalentDipolesInTetraMeshConsideringZCoordiante(int iSeg);

	void findEquivalentDipolesInDHexaMeshOnEarthSurface(int iSeg);
	
	void findEquivalentDipolesInTetraMeshOnEarthSurface(int iSeg);

	void calculateDirecctionVectorOnSurface(Point3D segStart, Point3D segEnd, int iSeg);

	// Copy constructer
	ObservedDataSourceSTLS(const ObservedDataSourceSTLS& rhs);

	// Copy assignment operator
	ObservedDataSourceSTLS& operator=(const ObservedDataSourceSTLS& rhs);

};

#endif
