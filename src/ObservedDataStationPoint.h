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
#ifndef DBLDEF_OBSERVED_DATA_STATION_POINT
#define DBLDEF_OBSERVED_DATA_STATION_POINT

#include <vector>
#include <complex>

#include "ObservedDataStation.h"

// Observed data of point station
class ObservedDataStationPoint: public ObservedDataStation{
	public:
		// Constructer
		explicit ObservedDataStationPoint();

		// Destructer
		~ObservedDataStationPoint();

		// Find element including a point
		void findElementIncludingStation();

		// Get caluculated value of Hx
		std::complex<double> getHxCalculated( const int iPol ) const;

		// Get caluculated value of Hy
		std::complex<double> getHyCalculated( const int iPol ) const;

		// Get right-hand side vector ID of Hx
		int getRhsVectorIDOfHx() const;

		// Get right-hand side vector ID of Hy
		int getRhsVectorIDOfHy() const;

		// Calulate horizontal magnetic field
		void calculateHorizontalMagneticField( const Forward3D* const ptrForward3D, const int rhsVectorIDOfHx, const int rhsVectorIDOfHy, bool useHx, bool useHy);

		// Initialize horizontal magnetic field
		void initializeHorizontalMagneticField( const int iPol );

		// Calulate interpolator vector of horizontal magnetic field
		virtual void calcInterpolatorVectorOfHorizontalMagneticField( Forward3D* const ptrForward3D );

		// Get location of the point
		const CommonParameters::locationXY& getLocationOfPoint() const;

		// Get Z coordinate of the point
		double getZCoordOfPoint() const;

		// Get flag specifing whether the EM field is interpolated from the values of the upper element
		bool useUpperElementForInterpolationOfEMField() const;

		// Set flag specifing whether the EM field is interpolated from the values of the upper element
		void setFlagUseUpperElementForInterpolationOfEMField( const bool useUpperElem );

	protected:
		// Location of the station
		CommonParameters::locationXY m_location;

		// Element including the station
		int m_elementIncludingStation;

		// Face including the station
		int m_faceIncludingStation;

		// Flag specifing whether the EM field is interpolated from the values of the upper element
		bool m_useUpperElementForInterpolationOfEMField;

		// local coordinate values of the location of the station
		CommonParameters::locationXYZ m_localCoordinateValues;

		// Volume coordinate values of the location of start point of line segment for tetrahedral mesh
		CommonParameters::VolumeCoords m_volumeCoordinateValues;

		// Area coordinate values of the location
		CommonParameters::AreaCoords m_areaCoordinateValues;

		std::complex<double> m_HxCalculated[2];
		std::complex<double> m_HyCalculated[2];

		int m_rhsVectorIDOfHx;
		int m_rhsVectorIDOfHy;

	private:
		// Copy constructer
		ObservedDataStationPoint(const ObservedDataStationPoint& rhs);
		
		// Copy assignment operator
		ObservedDataStationPoint& operator=(const ObservedDataStationPoint& rhs);

};

#endif
