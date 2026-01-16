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
