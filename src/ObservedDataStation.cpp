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
#include <iostream>
#include "ObservedDataStation.h"


//-------------------------------------------------------
//----- Definition of the class ObservedDataStation -----
//-------------------------------------------------------

// Constructer
ObservedDataStation::ObservedDataStation():
	m_stationID(0),
	m_IDOfMagneticFieldStation(0),
	m_numOfFrequency(0),
	m_freq(NULL),
	m_numOfFreqCalculatedByThisStaAndPE(0)
{
}

// Destructer
ObservedDataStation::~ObservedDataStation(){

	if( m_freq != NULL){
		delete[] m_freq;
		m_freq = NULL;
	}

}

// Get ID of station
int ObservedDataStation::getStationID() const{
	return m_stationID;
}

// Get total number of frequencies
int ObservedDataStation::getTotalNumberOfFrequency() const{
	return m_numOfFrequency;
}

// Get frequencies at which observed value exists
double ObservedDataStation::getFrequencyValues( const int num ) const{
	return m_freq[num];
}

// Find and return frequency IDs ( consecutive number in this station )  from frequency value.
// If the specified frequency is not included, return -1.
int ObservedDataStation::getFreqIDs( const double freq ) const{

	const double EPS = 1.0E-8;

	for( int i = 0; i < m_numOfFrequency; ++i ){
		if( std::fabs( freq - m_freq[i] ) < EPS ){
			return i;
		}
	}

	return -1;
}

// Set up frequencies calculated by this PE, at which observed value exists
void ObservedDataStation::setupFrequenciesCalculatedByThisPE( const int nFreqCalculatedByThisPE, const double* freqCalculatedByThisPE ){
	
	int icount(0);
	for( int i = 0; i < nFreqCalculatedByThisPE; ++i ){
		const int ifreq = getFreqIDs( freqCalculatedByThisPE[i] );
		if( ifreq >= 0 ){
			m_freqIDsAmongThisStationCalculatedByThisPE.push_back( ifreq );
		}
	}

	m_numOfFreqCalculatedByThisStaAndPE = static_cast<int>( m_freqIDsAmongThisStationCalculatedByThisPE.size() );

#ifdef _DEBUG_WRITE
	std::cout << "m_numOfFreqCalculatedByThisStaAndPE = " << m_numOfFreqCalculatedByThisStaAndPE << std::endl;
	std::cout << "m_freqIDsAmongThisStationCalculatedByThisPE : " << std::endl;
	for( std::vector<int>::iterator itr = m_freqIDsAmongThisStationCalculatedByThisPE.begin(); itr != m_freqIDsAmongThisStationCalculatedByThisPE.end(); ++itr ){
		std::cout << *itr << " " << m_freq[*itr] << std::endl;
	}
#endif

}

// Find and return frequency IDs among the ones calculated by this PE ( consecutive number in this PE and station ) from frequency value.
// If the specified frequency is not included, return -1.
int ObservedDataStation::getFreqIDsAmongThisPE( const double freq ) const{

	const double EPS = 1.0E-8;

	for( int i = 0; i < m_numOfFreqCalculatedByThisStaAndPE; ++i ){
		if( std::fabs( freq - m_freq[ m_freqIDsAmongThisStationCalculatedByThisPE[i] ] ) < EPS ){
			return i;
		}
	}

	return -1;

}

// Get ID of the station where magnetic field is observed
int ObservedDataStation::getIDOfMagneticFieldStation() const{
	return m_IDOfMagneticFieldStation;
}

