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
#ifndef DBLDEF_ROUGHENING_MATRIX
#define DBLDEF_ROUGHENING_MATRIX

#include "DoubleSparseMatrix.h"
#include "DoubleSparseSquareSymmetricMatrix.h"

class RougheningMatrix : public DoubleSparseMatrix {

public:

	//Default Constructer
	explicit RougheningMatrix();

	//Constructer
	explicit RougheningMatrix( const int nrows, const int ncols, const int nRhs = 1 );

	//Destructer
	virtual ~RougheningMatrix();
	
	//Input component of matrix.
	virtual void setStructureByTripletFormat( const int row, const int col );
	
	// Set matrix structure ( locations of non-zero components ) and add values by triplet format
	virtual void setStructureAndAddValueByTripletFormat( const int row, const int col, const double val );
	
	//Add non-zero value to matrix
	virtual void addNonZeroValues( const int row, const int col, const double val );
	
	//Make [R]T[R] matrix, where [R] is a constraining matrix
	void makeRTRMatrix( DoubleSparseSquareSymmetricMatrix& RTRMatrix, const double smallValueOnDiagonals = 0.0 ) const;

	// Calculate model roughness
	double calcModelRoughness( const double* modelVec ) const;

	// Calculate vector of model roughness
	void calcVectorOfModelRoughness( const double* modelVec, double* roughnessVec ) const;

	// Postmultiply diagonal matrix
	void postmultiplyDiagonalMatrix( const double* diagMatrix );

	// Output roughneing matrix
	void outputRougheningMatrix( const std::string& fileName ) const;

private:
	//Copy constructer
	RougheningMatrix(const RougheningMatrix &matrix );

	// Assignment operator
	RougheningMatrix& operator=(const RougheningMatrix& rhs);

};

#endif
