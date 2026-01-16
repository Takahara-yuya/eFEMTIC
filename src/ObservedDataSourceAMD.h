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
