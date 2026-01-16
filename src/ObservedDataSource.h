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
#ifndef DBLDEF_OBSERVED_DATA_SOURCE
#define DBLDEF_OBSERVED_DATA_SOURCE

#include <vector>
#include <complex>
#include <algorithm>

#include "Forward3D.h"
#include "CommonParameters.h"
#include "MeshDataTetraElement.h"

// Observed data of each station
class ObservedDataSource {
public:
	// Constructer
	explicit ObservedDataSource();

	// Destructer
	~ObservedDataSource();

	struct Point3D {
		double x, y, z;

		//
		Point3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

		//
		Point3D operator-(const Point3D& other) const {
			return { x - other.x, y - other.y, z - other.z };
		}

		//
		Point3D operator+(const Point3D& other) const {
			return { x + other.x, y + other.y, z + other.z };
		}

		//
		Point3D operator*(double scalar) const {
			return { x * scalar, y * scalar, z * scalar };
		}

		//
		bool isApprox(const Point3D& other, double epsilon = 1e-10) const {
			double dx = x - other.x;
			double dy = y - other.y;
			double dz = z - other.z;
			return dx * dx + dy * dy + dz * dz < epsilon * epsilon;
		}

		// 
		double distanceSquared(const Point3D& other) const {
			double dx = x - other.x;
			double dy = y - other.y;
			double dz = z - other.z;
			return dx * dx + dy * dy + dz * dz;
		}

		//
		double distanceSquaredXYPlane(const Point3D& other) const {
			double dx = x - other.x;
			double dy = y - other.y;
			return dx * dx + dy * dy;
		}
	};

	Point3D normalized(const Point3D& v, double epsilon = 1e-10) {
		double len = length(v);
		if (len < epsilon) return v;
		return { v.x / len, v.y / len, v.z / len };
	}


	double distance(const Point3D& a, const Point3D& b) {
		return std::sqrt(a.distanceSquared(b));
	}

	double distanceXYPlane(const Point3D& a, const Point3D& b) {
		return std::sqrt(a.distanceSquaredXYPlane(b));
	}

	//
	Point3D lerp(const Point3D& a, const Point3D& b, double t) {
		return {
			a.x + t * (b.x - a.x),
			a.y + t * (b.y - a.y),
			a.z + t * (b.z - a.z)
		};
	}

	// 
	Point3D cross(const Point3D& a, const Point3D& b) {
		return {
			a.y * b.z - a.z * b.y,
			a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x
		};
	}

	// 
	double crossXY(const Point3D& a, const Point3D& b) {
		return {
			a.x * b.y - a.y * b.x
		};
	}

	double dot(const Point3D& a, const Point3D& b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}

	double dotXY(const Point3D& a, const Point3D& b) {
		return a.x * b.x + a.y * b.y;
	}

	bool pointEqualXY(const Point3D& a, const Point3D& b) {
		return (std::fabs(a.x - b.x) < 1e-10) && (std::fabs(a.y - b.y) < 1e-10);
	}

	double length(const Point3D& v) {
		return std::sqrt(dot(v, v));
	}

	int orientation(Point3D a, Point3D b, Point3D c) {
		double val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);
		if (val > 1e-9) return 1;   // 
		if (val < -1e-9) return 2;  // 
		return 0;                   // 
	}

	bool isOnSegment(Point3D a, Point3D b, Point3D p) {
		return (p.x <= std::max(a.x, b.x) && p.x >= std::min(a.x, b.x) &&
			p.y <= std::max(a.y, b.y) && p.y >= std::min(a.y, b.y));
	}

	struct quadrilateralFace {
		Point3D p[4];
	};

	struct triangleFace {
		Point3D p[3];
	};

	enum typeOfSource {
		MAGNETIC_DIPOLE_SOURCE = 0,
		ELECTRIC_DIPOLE_SOURCE,
	};

	const double defaultValue = -999999999.;

	const int maxNumOfEquivalentDipoles = 99999.;

	void findElementIncludingDipoleSource();

	void assignmentPoint3D(Point3D& A, double xValue, double yValue, double zValue);

	CommonParameters::source dipoleParameter;

	bool segmentTriangleIntersection(const Point3D& p0, const Point3D& p1, const Point3D& v0, const Point3D& v1, const Point3D& v2, Point3D& intersection, double epsilon = 1e-10);

	bool isSegmentParallelToPlane(const Point3D& p0, const Point3D& p1,
		const Point3D& v0, const Point3D& v1, const Point3D& v2, double epsilon = 1e-10);

	bool isPointOnPlane(const Point3D& point,
		const Point3D& v0, const Point3D& v1, const Point3D& v2, double epsilon = 1e-10);

	bool isPointOnSegment(const Point3D& p, const Point3D& a, const Point3D& b, double epsilon = 1e-10);

	//
	bool checkIfPointExistsInVector(const Point3D& p, std::vector<Point3D> vp, double epsilon = 1e-10);

	int checkSegmentsIntersection(const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D, Point3D& outPoint1, Point3D& outPoint2);

	bool isPointInTetrahedron(const Point3D& P, const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D, double epsilon = 1e-10);
	//
	bool checkIfSegmentCrossOnePointOfThisFace(const Point3D result[3], Point3D& crossEnd, double epsilon = 1e-10);

	bool checkIfSegmentCrossTwoPointOfThisFace(const Point3D result[3], Point3D& crossStart, Point3D& crossEnd, double epsilon = 1e-10);

	bool checkIfSegmentCrossTwoPointOfHexaXYPlane(const Point3D result[4], Point3D& crossStart, Point3D& crossEnd, double epsilon = 1e-10);

	bool checkIfSegmentCrossOnePointOfThisTetrahedron(const Point3D intersection[4], Point3D& terminalPoint, Point3D& crossEnd, double epsilon = 1e-10);

	bool checkIfSegmentCrossTwoPointOfThisTetrahedron(const Point3D intersection[4], Point3D& crossStart, Point3D& crossEnd, double epsilon = 1e-10);

	bool isPointInTriangleXY(const Point3D node[3], const Point3D& p, double epsilon = 1e-10);

	bool isPointOnTriangleEdge(const Point3D node[3], const Point3D& p, double epsilon = 1e-10);

	bool isOverlapedEdgeHasBeenCalculatedInThisTetrahedron(const Point3D tetA, const Point3D tetB, const Point3D tetC, const Point3D tetD, std::vector<Point3D> overlapEdgeNodeCoordinateFirst, std::vector<Point3D> overlapEdgeNodeCoordinateSecond, double epsilon = 1e-10);

	bool isPointInNonConformingHexa(const std::vector<quadrilateralFace>& faces, const Point3D& pt, double epsilon = 1e-10);

	bool isSegmentParallelInXYPlane(const Point3D& p1, const Point3D& p2, const Point3D& q1, const Point3D& q2, double epsilon = 1e-10);

	int checkSegmentIntersectionInXYPlane(Point3D p1, Point3D p2, Point3D q1, Point3D q2, Point3D& result1, Point3D& result2);

	double getZOnLineSegment(const Point3D& p1, const Point3D& p2, double x, double y);

	double getZOnQuadPlane(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4, double x, double y);

	double getZOnTrianglePlane(const Point3D& p1, const Point3D& p2, const Point3D& p3, double x, double y);

	bool isPointOnThisHexaXYPlanes(const Point3D node[4], Point3D& p, double epsilon = 1e-10);

	bool isPointOnThisTatraSurfacePlane(const Point3D node[3], Point3D& p, double epsilon = 1e-10);
	//
	double findMax(const std::vector<double>& numbers);

	double findMin(const std::vector<double>& numbers);

	int getSourceID() const;

protected:
	//Source ID
	int m_sourceID;

	// Transform degrees to Radians
	double degreesToRadians(double degrees);

private:
	// Copy constructer
	ObservedDataSource(const ObservedDataSource& rhs);

	// Copy assignment operator
	ObservedDataSource& operator=(const ObservedDataSource& rhs);

};

#endif
