//-------------------------------------------------------------------------------------------------------
// This module provides CSEM (Controlled Source Electromagnetics) modeling capabilities within FEMTIC.
// Implemented by Zuwei Huang, 2025.
//
// For support, bug reports, or questions, please contact the author at: hzw1498218560@tongji.edu.cn
//-------------------------------------------------------------------------------------------------------
#include <iostream>
#include <optional>
#include "ObservedDataSource.h"
#include "OutputFiles.h"
#include "CommonParameters.h"
#include "ResistivityBlock.h"

//-------------------------------------------------------
//----- Definition of the class ObservedDataSource -----
//-------------------------------------------------------

// Constructer
ObservedDataSource::ObservedDataSource() :
	m_sourceID(-1)
{

	dipoleParameter.directionOfEquivalentDipoles = NULL;
	dipoleParameter.elementIncludingEquivalentDipoleInThisSource = NULL;
	dipoleParameter.midPoint = NULL;
	dipoleParameter.numOfEquivalentDipoleInThisSource = 0;
	dipoleParameter.sourceID = NULL;
	dipoleParameter.parameterOfEquivalentDipoleInThisSource = NULL;
	dipoleParameter.sourceTypeOfAllEquivalentDipoleSource = ELECTRIC_DIPOLE_SOURCE;

}

// Destructer
ObservedDataSource::~ObservedDataSource() {

	if (dipoleParameter.directionOfEquivalentDipoles != NULL) {
		delete[] dipoleParameter.directionOfEquivalentDipoles;
		dipoleParameter.directionOfEquivalentDipoles = NULL;
	}

	if (dipoleParameter.elementIncludingEquivalentDipoleInThisSource != NULL) {
		delete[] dipoleParameter.elementIncludingEquivalentDipoleInThisSource;
		dipoleParameter.elementIncludingEquivalentDipoleInThisSource = NULL;
	}

	if (dipoleParameter.midPoint != NULL) {
		delete[] dipoleParameter.midPoint;
		dipoleParameter.midPoint = NULL;
	}

	if (dipoleParameter.sourceID != NULL) {
		delete[] dipoleParameter.sourceID;
		dipoleParameter.sourceID = NULL;
	}

	if (dipoleParameter.parameterOfEquivalentDipoleInThisSource != NULL) {
		for (int i = 0; i < dipoleParameter.numOfEquivalentDipoleInThisSource; i++) {
			delete[] dipoleParameter.parameterOfEquivalentDipoleInThisSource[i];
		}
		delete[] dipoleParameter.parameterOfEquivalentDipoleInThisSource;
		dipoleParameter.parameterOfEquivalentDipoleInThisSource = NULL;
	}

}

double ObservedDataSource::degreesToRadians(double degrees) {

	return degrees * (3.1415926535897932384626433832795 / 180.0);

}

void ObservedDataSource::assignmentPoint3D(Point3D& A, double xValue, double yValue, double zValue) {

	A.x = xValue;
	A.y = yValue;
	A.z = zValue;

}

bool ObservedDataSource::segmentTriangleIntersection(const Point3D& p0, const Point3D& p1, const Point3D& v0, const Point3D& v1, const Point3D& v2, Point3D& intersection, double epsilon) {

	//using Möller-Trumbore algorithm
	Point3D segDir = p1 - p0;
	Point3D edge1 = v1 - v0;
	Point3D edge2 = v2 - v0;

	Point3D pvec = cross(segDir, edge2);
	double det = dot(edge1, pvec);

	if (std::abs(det) < epsilon) return false;

	double invDet = 1.0 / det;
	Point3D tvec = p0 - v0;

	double u = dot(tvec, pvec) * invDet;
	if (u < -epsilon || u > 1.0 + epsilon) {
		assignmentPoint3D(intersection, defaultValue, defaultValue, defaultValue);
		return false;
	}

	Point3D qvec = cross(tvec, edge1);
	double v = dot(segDir, qvec) * invDet;
	if (v < -epsilon || u + v > 1.0 + epsilon) {
		assignmentPoint3D(intersection, defaultValue, defaultValue, defaultValue);
		return false;
	}

	double t = dot(edge2, qvec) * invDet;
	if (t < -epsilon || t > 1.0 + epsilon) {
		assignmentPoint3D(intersection, defaultValue, defaultValue, defaultValue);
		return false;
	}

	intersection = p0 + segDir * t;
	return true;

}

// Determine whether the line segment is parallel to the plane where the triangle is located
bool ObservedDataSource::isSegmentParallelToPlane(const Point3D& p0, const Point3D& p1,
	const Point3D& v0, const Point3D& v1, const Point3D& v2, double epsilon) {
	// Calculate the normal vector of the plane
	Point3D edge1 = v1 - v0;
	Point3D edge2 = v2 - v0;
	Point3D normal = cross(edge1, edge2);

	// Check if the triangle has degraded (with an area of zero)
	double normalLengthSq = dot(normal, normal);
	if (normalLengthSq < epsilon * epsilon) {
		// Triangle degeneration, unable to define meaningful planes
		return false;
	}

	// Calculate the direction vector of the line segment
	Point3D dir = p1 - p0;

	// Calculate the dot product of line segment vectors and normal vectors
	double dotProduct = dot(normal, dir);

	// If the absolute value of the dot product is close to zero, then the line segment is parallel to the plane
	return std::abs(dotProduct) < epsilon;
}

bool ObservedDataSource::isPointOnPlane(const Point3D& point,
	const Point3D& v0, const Point3D& v1, const Point3D& v2, double epsilon) {

	Point3D edge1 = v1 - v0;
	Point3D edge2 = v2 - v0;
	Point3D normal = cross(edge1, edge2);


	Point3D vec = point - v0;
	double distance = dot(normal, vec);


	return std::abs(distance) < epsilon;
}

// Determine whether the point is on the line segment
bool ObservedDataSource::isPointOnSegment(const Point3D& p, const Point3D& a, const Point3D& b, double epsilon) {
	Point3D ap = p - a;
	Point3D ab = b - a;

	Point3D cross_prod = cross(ap, ab);
	double cross_length_sq = dot(cross_prod, cross_prod);
	if (cross_length_sq > epsilon) {
		return false;
	}


	double ab_sq = dot(ab, ab);
	if (ab_sq < epsilon) { // A line segment is a point
		return ap.isApprox(Point3D(0, 0, 0), epsilon);
	}

	double t = dot(ap, ab) / ab_sq;
	return t >= -epsilon && t <= 1.0 + epsilon;
}

// Main function: Determine the relationship between line segments and return intersection information through references
int ObservedDataSource::checkSegmentsIntersection(const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D, Point3D& outPoint1, Point3D& outPoint2) {
	double epsilon = 1e-10;

	// 0.Complete equality check
	if ((A.isApprox(C, epsilon) || A.isApprox(D, epsilon)) && (B.isApprox(C, epsilon) || B.isApprox(D, epsilon))) {
		outPoint1 = A;
		outPoint2 = B;
		return 1;
	}

	// 1. Endpoint coincidence check
	//
	// 2. Calculate direction vector
	Point3D AB = B - A;
	Point3D CD = D - C;
	Point3D AC = C - A;

	// 3. Check coplanarity (mixed product)
	Point3D cross_ab_ac = cross(AB, AC);
	double triple_scalar = dot(cross_ab_ac, CD);
	if (std::abs(triple_scalar) > epsilon) {
		outPoint1.x = defaultValue;
		outPoint1.y = defaultValue;
		outPoint1.z = defaultValue;
		outPoint2.x = defaultValue;
		outPoint2.y = defaultValue;
		outPoint2.z = defaultValue;
		return -1; // Line segments are not coplanar or intersect
	}

	// 4. Check if they are collinear
	Point3D cross_ab_cd = cross(AB, CD);
	double cross_length_sq = dot(cross_ab_cd, cross_ab_cd);
	if (cross_length_sq < epsilon * epsilon) {

		std::vector<Point3D> overlapPoints;

		if (isPointOnSegment(C, A, B, epsilon)) overlapPoints.push_back(C);
		if (isPointOnSegment(D, A, B, epsilon)) overlapPoints.push_back(D);
		if (isPointOnSegment(A, C, D, epsilon)) overlapPoints.push_back(A);
		if (isPointOnSegment(B, C, D, epsilon)) overlapPoints.push_back(B);

		if (overlapPoints.empty()) {
			return -1;
		}

		for (size_t i = 0; i < overlapPoints.size(); i++) {
			for (size_t j = i + 1; j < overlapPoints.size(); j++) {
				if (overlapPoints[i].isApprox(overlapPoints[j], epsilon)) {
					overlapPoints.erase(overlapPoints.begin() + j);
					j--;
				}
			}
		}

		if (overlapPoints.size() == 1) {
			outPoint1 = overlapPoints[0];
			return 0;
		}

		Point3D direction = normalized(AB, epsilon);
		auto project = [&](const Point3D& p) {
			Point3D vec = p - A;
			return direction.x * vec.x + direction.y * vec.y + direction.z * vec.z;
			};

		std::sort(overlapPoints.begin(), overlapPoints.end(),
			[&](const Point3D& a, const Point3D& b) {
				return project(a) < project(b);
			});

		outPoint1 = overlapPoints.front();
		outPoint2 = overlapPoints.back();
		return 1;
	}

	// 5.Not collinear, but in the same plane, solve the parametric equation
	double t = 0, s = 0;
	double det_xy = AB.x * CD.y - AB.y * CD.x;
	double det_xz = AB.x * CD.z - AB.z * CD.x;
	double det_yz = AB.y * CD.z - AB.z * CD.y;

	double abs_xy = std::abs(det_xy);
	double abs_xz = std::abs(det_xz);
	double abs_yz = std::abs(det_yz);

	if (abs_xy >= abs_xz && abs_xy >= abs_yz) {
		t = (CD.y * AC.x - CD.x * AC.y) / det_xy;
		s = (AB.y * AC.x - AB.x * AC.y) / det_xy;
	}
	else if (abs_xz >= abs_xy && abs_xz >= abs_yz) {
		t = (CD.z * AC.x - CD.x * AC.z) / det_xz;
		s = (AB.z * AC.x - AB.x * AC.z) / det_xz;
	}
	else {
		t = (CD.z * AC.y - CD.y * AC.z) / det_yz;
		s = (AB.z * AC.y - AB.y * AC.z) / det_yz;
	}

	if (t >= -epsilon && t <= 1.0 + epsilon &&
		s >= -epsilon && s <= 1.0 + epsilon) {

		outPoint1 = A + AB * t;
		return 0;
	}

	assignmentPoint3D(outPoint1, defaultValue, defaultValue, defaultValue);
	assignmentPoint3D(outPoint2, defaultValue, defaultValue, defaultValue);
	return -1;
}

bool ObservedDataSource::isPointInTetrahedron(const Point3D& P, const Point3D& A, const Point3D& B, const Point3D& C, const Point3D& D, double epsilon) {

	Point3D AB = B - A;
	Point3D AC = C - A;
	Point3D AD = D - A;
	Point3D cross_AB_AC = cross(AB, AC);
	double V = std::abs(dot(cross_AB_AC, AD)) / 6.0;

	// V_A: P, B, C, D
	Point3D PB = B - P;
	Point3D PC = C - P;
	Point3D PD = D - P;
	Point3D cross_PB_PC = cross(PB, PC);
	double VA = std::abs(dot(cross_PB_PC, PD)) / 6.0;

	// V_B: A, P, C, D
	Point3D PA = A - P;
	Point3D cross_PA_PC = cross(PA, PC);
	double VB = std::abs(dot(cross_PA_PC, PD)) / 6.0;

	// V_C: A, B, P, D
	Point3D cross_PA_PB = cross(PA, PB);
	double VC = std::abs(dot(cross_PA_PB, PD)) / 6.0;

	// V_D: A, B, C, P
	Point3D cross_PA_PB_2 = cross(PA, PB);
	double VD = std::abs(dot(cross_PA_PB_2, PC)) / 6.0;

	double total_volume = VA + VB + VC + VD;
	if (std::abs(total_volume - V) > 1e-10) {
		return false;
	}

	double alpha = VA / V;
	double beta = VB / V;
	double gamma = VC / V;
	double delta = VD / V;

	if (alpha >= 0 && beta >= 0 && gamma >= 0 && delta >= 0 &&
		std::abs(alpha + beta + gamma + delta - 1.0) < 1e-10) {
		return true;
	}

	return false;
	//////////////////////
	//
	//Point3D AB = B - A;
	//Point3D AC = C - A;
	//Point3D AD = D - A;
	//Point3D AP = P - A;

	//
	//double volume = dot(AB, cross(AC, AD));

	//
	//if (std::abs(volume) < epsilon) return false;

	//
	//double v0 = dot(AP, cross(AC, AD));  // 
	//double v1 = dot(AB, cross(AP, AD));  // 
	//double v2 = dot(AB, cross(AC, AP));  // 
	//double v3 = dot(AP, cross(AB, AC));  //

	//
	//if (volume > 0) {
	//	return (v0 >= -epsilon) && (v1 >= -epsilon) &&
	//		(v2 >= -epsilon) && (v3 >= -epsilon);
	//}
	//else {
	//	return (v0 <= epsilon) && (v1 <= epsilon) &&
	//		(v2 <= epsilon) && (v3 <= epsilon);
	//}
	///////////////////
	//bool xin(false);
	//bool yin(false);
	//bool zin(false);
	////
	//std::vector<double> x, y, z;
	//x.resize(4); y.resize(4); z.resize(4);
	//x[0] = A.x; x[1] = B.x; x[2] = C.x; x[3] = D.x;
	//y[0] = A.y; y[1] = B.y; y[2] = C.y; y[3] = D.y;
	//z[0] = A.z; z[1] = B.z; z[2] = C.z; z[3] = D.z;
	////
	//double xmax = findMax(x);
	//double xmin = findMin(x);
	//double ymax = findMax(y);
	//double ymin = findMin(y);
	//double zmax = findMax(z);
	//double zmin = findMin(z);
	//for (int i = 0; i < 4; i++) {
	//	if (P.x >= xmin && P.x <= xmax) {
	//		xin = true;
	//	}
	//	if (P.y >= ymin && P.y <= ymax) {
	//		yin = true;
	//	}
	//	if (P.z >= zmin && P.z <= zmax) {
	//		zin = true;
	//	}
	//}
	//if (xin && yin && zin) {
	//	return true;
	//}
	//else {
	//	return false;
	//}
}

bool ObservedDataSource::checkIfSegmentCrossTwoPointOfThisFace(const Point3D result[3], Point3D& crossStart, Point3D& crossEnd, double epsilon) {
	//check if the cross point = 2
	Point3D result1[2];
	int count(0);
	bool duplicate(false);
	for (int i = 0; i < 3; i++) {
		if (result[i].x != defaultValue) {
			for (int j = 0; j < i; j++) {
				if (distance(result[i], result[j]) < epsilon)duplicate = true;
			}
			if (!duplicate) {
				count++;
				result1[count - 1] = result[i];
				duplicate = false;
			}
		}
	}
	if (count != 2) {
		return false;
	}
	else
	{
		crossStart = result1[0];
		crossEnd = result1[1];
		return true;
	}
}

bool ObservedDataSource::checkIfSegmentCrossTwoPointOfHexaXYPlane(const Point3D result[4], Point3D& crossStart, Point3D& crossEnd, double epsilon) {
	//check if the cross point = 2
	Point3D result1[2];
	int count(0);
	bool duplicate(false);
	for (int i = 0; i < 4; i++) {
		if (result[i].x != defaultValue) {
			for (int j = 0; j < i; j++) {
				if (distanceXYPlane(result[i], result[j]) < epsilon)duplicate = true;
			}
			if (!duplicate) {
				count++;
				result1[count - 1] = result[i];
				duplicate = false;
			}
		}
	}
	if (count != 2) {
		return false;
	}
	else
	{
		crossStart = result1[0];
		crossEnd = result1[1];
		return true;
	}
}

bool ObservedDataSource::checkIfSegmentCrossOnePointOfThisFace(const Point3D result[3], Point3D& crossEnd, double epsilon) {
	//check if the cross point = 1
	Point3D result1, terminal;
	terminal = crossEnd;
	int count(0), count1(0);
	bool duplicate(false);
	for (int i = 0; i < 3; i++) {
		if (result[i].x == defaultValue)count++;
		if (distance(result[i], terminal) < epsilon) {
			result1 = crossEnd;
			count1++;
		}
	}
	if (count == 2 && count1 == 1) {
		crossEnd = result1;
		return true;
	}
	else {
		count = 0;
		for (int i = 0; i < 3; i++) {
			if (result[i].x != defaultValue && distance(result[i], terminal) > epsilon) {
				for (int j = 0; j < i; j++) {
					if (distance(result[i], result[j]) < epsilon)duplicate = true;
				}
				if (!duplicate) {
					count++;
					result1 = result[i];
					duplicate = false;
				}
			}
		}
		if (count != 1) {
			return false;
		}
		else
		{
			crossEnd = result1;
			return true;
		}
	}
}

bool ObservedDataSource::checkIfSegmentCrossTwoPointOfThisTetrahedron(const Point3D intersection[4], Point3D& crossStart, Point3D& crossEnd, double epsilon) {
	//check if the cross point = 2
	Point3D result1[2];
	int count(0);
	bool crossThroughOnePoint(false);
	for (int i = 0; i < 4; i++) {
		if (intersection[i].x != defaultValue) {
			for (int j = 0; j < i; j++) {
				if (distance(intersection[j], intersection[i]) < epsilon)crossThroughOnePoint = true;
			}
			if (!crossThroughOnePoint) {
				count++;
				result1[count - 1] = intersection[i];
				crossThroughOnePoint = false;
			}
		}
	}
	if (count != 2) {
		return false;
	}
	else
	{
		crossStart = result1[0];
		crossEnd = result1[1];
		return true;
	}
}

bool ObservedDataSource::checkIfSegmentCrossOnePointOfThisTetrahedron(const Point3D intersection[4], Point3D& terminalPoint, Point3D& crossEnd, double epsilon) {
	//check if the cross point = 1
	Point3D result1;
	int count(0);
	bool crossThroughOnePoint(false);
	for (int i = 0; i < 4; i++) {
		if (intersection[i].x != defaultValue) {
			for (int j = 0; j < i; j++) {
				if (distance(intersection[j], intersection[i]) < epsilon)crossThroughOnePoint = true;
				if (distance(terminalPoint, intersection[i]) < epsilon)crossThroughOnePoint = true;
			}
			if (!crossThroughOnePoint) {
				count++;
				result1 = intersection[i];
				crossThroughOnePoint = false;
			}
		}
	}
	if (count != 1) {
		return false;
	}
	else
	{
		crossEnd = result1;
		return true;
	}
}

bool ObservedDataSource::isOverlapedEdgeHasBeenCalculatedInThisTetrahedron(const Point3D tetA, const Point3D tetB, const Point3D tetC, const Point3D tetD, std::vector<Point3D> overlapEdgeNodeCoordinateFirst, std::vector<Point3D> overlapEdgeNodeCoordinateSecond, double epsilon) {
	int numOfOverlapedEdge(0);
	int nEdge = overlapEdgeNodeCoordinateFirst.size();
	for (int iEdge = 0; iEdge < nEdge; iEdge++) {
		if (distance(tetA, overlapEdgeNodeCoordinateFirst[iEdge]) < epsilon || distance(tetA, overlapEdgeNodeCoordinateSecond[iEdge]) < epsilon) {
			numOfOverlapedEdge++;
		}
		if (distance(tetB, overlapEdgeNodeCoordinateFirst[iEdge]) < epsilon || distance(tetB, overlapEdgeNodeCoordinateSecond[iEdge]) < epsilon) {
			numOfOverlapedEdge++;
		}
		if (distance(tetC, overlapEdgeNodeCoordinateFirst[iEdge]) < epsilon || distance(tetC, overlapEdgeNodeCoordinateSecond[iEdge]) < epsilon) {
			numOfOverlapedEdge++;
		}
		if (distance(tetD, overlapEdgeNodeCoordinateFirst[iEdge]) < epsilon || distance(tetD, overlapEdgeNodeCoordinateSecond[iEdge]) < epsilon) {
			numOfOverlapedEdge++;
		}
		if (numOfOverlapedEdge > 1) {
			return true;
		}
		numOfOverlapedEdge = 0;
	}
	return false;
}

double ObservedDataSource::findMax(const std::vector<double>& numbers) {
	if (numbers.empty()) {
		throw std::invalid_argument("Vector is empty!");
	}
	return *max_element(numbers.begin(), numbers.end());
}

double ObservedDataSource::findMin(const std::vector<double>& numbers) {
	if (numbers.empty()) {
		throw std::invalid_argument("Vector is empty!");
	}
	return *min_element(numbers.begin(), numbers.end());
}

int ObservedDataSource::getSourceID() const {

	return m_sourceID;

}

// For DHEXA

bool ObservedDataSource::isPointInNonConformingHexa(const std::vector<quadrilateralFace>& faces, const Point3D& pt, double epsilon) {
	if (faces.size() != 6) return false;

	Point3D center(0, 0, 0);
	int count = 0;
	for (const auto& face : faces) {
		for (int i = 0; i < 4; i++) {
			center.x += face.p[i].x;
			center.y += face.p[i].y;
			center.z += face.p[i].z;
			count++;
		}
	}
	center.x /= count;
	center.y /= count;
	center.z /= count;

	int sign = 0;
	for (const auto& face : faces) {
		Point3D v1 = face.p[1] - face.p[0];
		Point3D v2 = face.p[2] - face.p[0];
		Point3D normal = cross(v1, v2);

		if (dot(normal, face.p[0] - center) < 0) {
			normal = normal * -1.0;
		}

		Point3D vec = pt - face.p[0];
		double d = dot(normal, vec);

		if (std::fabs(d) < epsilon) continue;

		int currSign = (d > 0) ? 1 : -1;
		if (sign == 0) {
			sign = currSign;
		}
		else if (sign != currSign) {
			return false;
		}
	}
	return true;
}

bool ObservedDataSource::isSegmentParallelInXYPlane(const Point3D& p1, const Point3D& p2, const Point3D& q1, const Point3D& q2, double epsilon) {
	double x1 = p1.x, y1 = p1.y;
	double x2 = p2.x, y2 = p2.y;
	double x3 = q1.x, y3 = q1.y;
	double x4 = q2.x, y4 = q2.y;

	double denom = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);

	if (std::fabs(denom) < epsilon) {
		return true;
	}

	return false;
}

int ObservedDataSource::checkSegmentIntersectionInXYPlane(Point3D p1, Point3D p2, Point3D q1, Point3D q2, Point3D& result1, Point3D& result2) {

	double epsilon = 1e-10;

	assignmentPoint3D(result1, defaultValue, defaultValue, defaultValue);
	assignmentPoint3D(result2, defaultValue, defaultValue, defaultValue);

	Point3D r = p2 - p1;
	Point3D s = q2 - q1;
	Point3D qp = q1 - p1;

	double rxs = crossXY(r, s);
	double qpxr = crossXY(qp, r);

	//  Non-parallel scenario: the only intersection point (possibly outside the segment)
	if (std::fabs(rxs) > epsilon) {
		double t = crossXY(qp, s) / rxs; // p + t*r = intersection
		double u = crossXY(qp, r) / rxs; // q + u*s = intersection

		if (t >= -epsilon && t <= 1.0 + epsilon && u >= -epsilon && u <= 1.0 + epsilon) {
			Point3D inter = p1 + r * t;
			result1 = inter;
			return 1;
		}
		else {
			return -1; // The lines intersect, but the line segments do not
		}
	}

	// r x s == 0 -> Parallel or collinear
	if (std::fabs(qpxr) > epsilon) {
		// Parallel but not collinear
		return -1;
	}

	// Collinear: Project q1 and q2 onto the parameter t of p1->p2 (t is [0,1] on segment p)
	double r_dot_r = dotXY(r, r);
	// If segment p degenerates into a point
	if (r_dot_r < epsilon) {
		// p
		if (pointEqualXY(p1, q1) || pointEqualXY(p1, q2)) {
			result1 = p1;
			return 1;
		}
		// Check if p1 is on segment q
		double s_dot_s = dot(s, s);
		if (s_dot_s < epsilon) return false; // q also degenerates into a point
		double u = dot(p1 - q1, s) / s_dot_s;
		if (u >= -epsilon && u <= 1.0 + epsilon) {
			result1 = p1;
			return 1;
		}
		return -1;
	}

	// Calculate the parameters t0 and t1 of q1 and q2 on p1->p2
	double t0 = dot(q1 - p1, r) / r_dot_r;
	double t1 = dot(q2 - p1, r) / r_dot_r;

	double ta = std::min(t0, t1);
	double tb = std::max(t0, t1);

	double segA = std::max(0.0, ta);
	double segB = std::min(1.0, tb);

	if (segA > segB + epsilon) {
		// No overlapping intervals
		return -1;
	}

	if (std::fabs(segA - segB) < epsilon) {
		// Overlap on a single point
		Point3D pt = p1 + r * segA;
		result1 = pt;
		return 1;
	}

	// Overlap into an interval, return the endpoint of the interval (from small to large in the direction of p1->p2)
	Point3D a = p1 + r * segA;
	Point3D b = p1 + r * segB;
	result1 = a;
	result2 = b;
	return 0;

}

double ObservedDataSource::getZOnLineSegment(const Point3D& p1, const Point3D& p2, double x, double y) {

	double t;
	if (fabs(p2.x - p1.x) > 1e-9) {
		t = (x - p1.x) / (p2.x - p1.x);
	}
	else if (fabs(p2.y - p1.y) > 1e-9) {
		t = (y - p1.y) / (p2.y - p1.y);
	}
	else {
		return p1.z;
	}

	return p1.z + t * (p2.z - p1.z);

}

double ObservedDataSource::getZOnQuadPlane(const Point3D& p1, const Point3D& p2, const Point3D& p3, const Point3D& p4, double x, double y) {
	//only use three points to verify a specific plane
	Point3D v1 = p2 - p1;
	Point3D v2 = p3 - p1;
	Point3D n = cross(v1, v2);

	double A = n.x;
	double B = n.y;
	double C = n.z;

	// Ax + By + Cz + D = 0
	double D = -(A * p1.x + B * p1.y + C * p1.z);

	if (fabs(C) < 1e-12) {
		OutputFiles::m_logFile << " getZOnQuadPlane function error ! " << std::endl;
		exit(1);
	}

	return (-(A * x + B * y + D)) / C;
}

double ObservedDataSource::getZOnTrianglePlane(const Point3D& p1, const Point3D& p2, const Point3D& p3, double x, double y) {
	//only use three points to verify a specific plane
	Point3D v1 = p2 - p1;
	Point3D v2 = p3 - p1;
	Point3D n = cross(v1, v2);

	double A = n.x;
	double B = n.y;
	double C = n.z;

	// Ax + By + Cz + D = 0
	double D = -(A * p1.x + B * p1.y + C * p1.z);

	if (fabs(C) < 1e-12) {
		OutputFiles::m_logFile << " getZOnTrianglePlane function error ! " << std::endl;
		exit(1);
	}

	return (-(A * x + B * y + D)) / C;
}

bool ObservedDataSource::isPointOnThisHexaXYPlanes(const Point3D node[4], Point3D& p, double epsilon) {
	double xmin, xmax, ymin, ymax;
	xmin = node[0].x;
	xmax = node[0].x;
	ymin = node[0].y;
	ymax = node[0].y;
	for (int i = 1; i < 4; i++) {
		if (node[i].x < xmin) {
			xmin = node[i].x;
		}
		if (node[i].y < ymin) {
			ymin = node[i].y;
		}
		if (node[i].x > xmax) {
			xmax = node[i].x;
		}
		if (node[i].y > ymax) {
			ymax = node[i].y;
		}
	}

	if (p.x >= xmin && p.x <= xmax && p.y >= ymin && p.y <= ymax) {
		return true;
	}
	else {
		return false;
	}
}

bool ObservedDataSource::isPointOnThisTatraSurfacePlane(const Point3D node[3], Point3D& p, double epsilon) {

	for (int i = 0; i < 3; ++i) {
		if (p.distanceSquared(node[i]) < epsilon * epsilon) {
			return true; 
		}
	}

	Point3D v1 = node[1] - node[0];
	Point3D v2 = node[2] - node[0];
	Point3D normal = cross(v1, v2);
	
	double normalLengthSq = normal.distanceSquared(Point3D(0, 0, 0));
	if (normalLengthSq < epsilon * epsilon) {
		return isPointInTriangleXY(node, p, epsilon);
	}

	Point3D v0 = p - node[0];
	double planeTest = dot(v0, normal);

	if (fabs(planeTest) > epsilon) {
		return false;
	}

	double totalArea = length(normal);  

	
	Point3D cross1 = cross(node[1] - node[0], p - node[0]);
	Point3D cross2 = cross(node[2] - node[1], p - node[1]);
	Point3D cross3 = cross(node[0] - node[2], p - node[2]);

	double area1 = length(cross1);
	double area2 = length(cross2);
	double area3 = length(cross3);

	
	double areaSum = area1 + area2 + area3;

	if (fabs(areaSum - totalArea) < epsilon) {
	
		double dot1 = dot(cross1, normal);
		double dot2 = dot(cross2, normal);
		double dot3 = dot(cross3, normal);

		if (dot1 > -epsilon && dot2 > -epsilon && dot3 > -epsilon) {
			return true;
		}
	}


	if (isPointOnTriangleEdge(node, p, epsilon)) {
		return true;
	}

	return false;
}

bool ObservedDataSource::checkIfPointExistsInVector(const Point3D& p, std::vector<Point3D> vp, double epsilon) {
	for (const auto& point : vp) {
		if (distance(point, p) < epsilon) {
			return true;
		}
	}
	return false;
}

bool ObservedDataSource::isPointOnTriangleEdge(const Point3D node[3], const Point3D& p, double epsilon) {
	for (int i = 0; i < 3; ++i) {
		int j = (i + 1) % 3;

		Point3D edge = node[j] - node[i];
		Point3D toPoint = p - node[i];

		double edgeLengthSq = edge.distanceSquared(Point3D(0, 0, 0));
		if (edgeLengthSq < epsilon * epsilon) {
			continue;
		}

		double t = dot(toPoint, edge) / edgeLengthSq;

		if (t >= -epsilon && t <= 1.0 + epsilon) {
			Point3D projection = node[i] + edge * t;

			if (p.distanceSquared(projection) < epsilon * epsilon) {
				return true;
			}
		}
	}

	return false;
}


bool ObservedDataSource::isPointInTriangleXY(const Point3D node[3], const Point3D& p, double epsilon) {

	Point3D v0 = node[2] - node[0];
	Point3D v1 = node[1] - node[0];
	Point3D v2 = p - node[0];

	double dot00 = dotXY(v0, v0);
	double dot01 = dotXY(v0, v1);
	double dot02 = dotXY(v0, v2);
	double dot11 = dotXY(v1, v1);
	double dot12 = dotXY(v1, v2);

	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1.0 + epsilon);
}
