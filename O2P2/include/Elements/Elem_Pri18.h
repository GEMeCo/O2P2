// ================================================================================================
// 
// This file is part of O2P2, an object oriented environment for the positional FEM
//
// Copyright(C) 2023 GEMeCO - All Rights Reserved.
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
//
// Solid element, with quadratic interpolation functions, prism shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

namespace O2P2 {
	namespace Geom {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Pri18
			  *
			  * @brief Prismatic quadratic element with 18 nodes.
			  * @details Solid element, with quadratic interpolation functions, prism shape.
			  * @image html Elem_Pri18.png height=300
			  */
			class Elem_Pri18 : public ElemSolid
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Pri18() = delete;

			public:
				/** Constructor for prism quadratic elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Pri18(std::shared_ptr<O2P2::Geom::Material>& Material)
					: ElemSolid(Material) { }

				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;

					msg << "2 2 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " " << this->mv_Conect[2]->mv_index + add << " "
						<< this->mv_Conect[3]->mv_index + add << " " << this->mv_Conect[4]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << this->mv_Conect[12]->mv_index + add << " " << this->mv_Conect[13]->mv_index + add << " " << this->mv_Conect[14]->mv_index + add << " "
						<< this->mv_Conect[15]->mv_index + add << " " << this->mv_Conect[16]->mv_index + add << " " << this->mv_Conect[17]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 2 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " " << this->mv_Conect[2]->mv_index + add << " "
						<< this->mv_Conect[6]->mv_index + add << " " << this->mv_Conect[7]->mv_index + add << " " << this->mv_Conect[8]->mv_index + add << " "
						<< this->mv_Conect[12]->mv_index + add << " " << this->mv_Conect[13]->mv_index + add << " " << this->mv_Conect[14]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 2 " << this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[4]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " "
						<< this->mv_Conect[8]->mv_index + add << " " << this->mv_Conect[10]->mv_index + add << " " << this->mv_Conect[11]->mv_index + add << " "
						<< this->mv_Conect[14]->mv_index + add << " " << this->mv_Conect[16]->mv_index + add << " " << this->mv_Conect[17]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 2 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " "
						<< this->mv_Conect[6]->mv_index + add << " " << this->mv_Conect[9]->mv_index + add << " " << this->mv_Conect[11]->mv_index + add << " "
						<< this->mv_Conect[12]->mv_index + add << " " << this->mv_Conect[15]->mv_index + add << " " << this->mv_Conect[17]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";

					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;

					msg << "2 2 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (4 + add) << " "
						<< (5 + add) << " " << (6 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << (13 + add) << " " << (14 + add) << " " << (15 + add) << " " << (16 + add) << " "
						<< (17 + add) << " " << (18 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 2 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (7 + add) << " " << (8 + add) << " "
						<< (9 + add) << " " << (13 + add) << " " << (14 + add) << " " << (15 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 2 " << (3 + add) << " " << (5 + add) << " " << (6 + add) << " " << (9 + add) << " " << (11 + add) << " "
						<< (12 + add) << " " << (15 + add) << " " << (17 + add) << " " << (18 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 2 " << (1 + add) << " " << (4 + add) << " " << (6 + add) << " " << (7 + add) << " " << (10 + add) << " "
						<< (12 + add) << " " << (13 + add) << " " << (16 + add) << " " << (18 + add) << " " << this->mv_Mat->mv_index << "\n";

					return msg.str();
				}

				// Evaluates shape function in the point.
				std::vector<double> getShapeFcOnPoint(const double* Point) override;

				// Evaluates the derivative of shape function in the point.
				std::vector<double> getShapeDerivOnPoint(const double* Point) override;

				// Return a vector with values on the integration points currently known in the element' nodes.
				std::vector<double> getValueOnIPs(const double* value) override;

				// Returns a pointer to the first element of the shape functions (with size [nIP][mv_numNodes]).
				double const* getShapeFc() const override { return &mv_Psi[0][0]; }

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][mv_numNodes][mv_ElDim]).
				double const* getShapeDerivative() const override { return &mv_DPsi[0][0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return &mv_weight[0]; }

				// Returns the number of nodes of current element.
				const int getNumNodes() const override { return mv_numNodes; }

				// Returns the number of faces of current element.
				const int getNumFaces() const override { return mv_numFaces; }

				// Returns the number of integration points of current element.
				const int getNumIP() const override { return mv_numIP; }

				// Verifies dimensionless coordinates from input - if it is immersed on the element.
				inline bool evaluateXsi(const double* xsi) override {
					std::array<double, mv_ElDim + 1> new_xsi = {};

					for (int i = 0; i < mv_ElDim - 1; ++i) {
						new_xsi.at(i) = *(xsi + i);
						new_xsi.at(mv_ElDim - 1) -= *(xsi + i);
					}
					new_xsi.at(mv_ElDim - 1) += 1.;
					new_xsi.at(mv_ElDim) = *(xsi + mv_ElDim);

					if (*std::max_element(new_xsi.begin(), new_xsi.end()) < 1.000001 && *std::min_element(new_xsi.begin(), new_xsi.end() - 1) > -0.000001 && new_xsi.at(mv_ElDim) > -1.000001) return true;
					return false;
				}

			private:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			private:
				/** @brief Number of Nodes */
				static const int mv_numNodes{ 18 };

				/** @brief Number of Integration Points */
				static const int mv_numIP{ 18 };

				/** @brief Number of Faces */
				static const int mv_numFaces{ 5 };

				/** @brief Weights for numerical integration */
				static const double mv_weight[mv_numIP];

				/** @brief Shape functions */
				static const double mv_Psi[mv_numIP][mv_numNodes];

				/** @brief Shape functions derivative */
				static const double mv_DPsi[mv_numIP][mv_numNodes][mv_ElDim];

				/** @brief Integration points */
				static const double mv_xsi[mv_numIP][mv_ElDim];
			};
		} // End of Elem Namespace
	} // End of Geom Namespace
} // End of O2P2 Namespace


// ================================================================================================
//
// Implementation of Member Function: getShapeFcOnPoint
// Shape functions evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Geom::Elem::Elem_Pri18::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(18);

	mi_Psi.at(0) = 0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. - Point[2]) * Point[2];
	mi_Psi.at(1) = 2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. - Point[2]) * Point[2];
	mi_Psi.at(2) = -0.5 * Point[0] * (2. * Point[0] - 1.) * (1. - Point[2]) * Point[2];
	mi_Psi.at(3) = 2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. - Point[2]) * Point[2];
	mi_Psi.at(4) = -2. * Point[0] * Point[1] * Point[2] * (1. - Point[2]);
	mi_Psi.at(5) = -0.5 * Point[1] * (2. * Point[1] - 1.) * (1. - Point[2]) * Point[2];
	mi_Psi.at(6) = (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (-1. + Point[2] * Point[2]);
	mi_Psi.at(7) = 4. * Point[0] * (-1. + Point[0] + Point[1]) * (-1. + Point[2] * Point[2]);
	mi_Psi.at(8) = Point[0] * (2. * Point[0] - 1.) * (1. - Point[2] * Point[2]);
	mi_Psi.at(9)  = 4. * Point[1] * (-1. + Point[0] + Point[1]) * (-1. + Point[2] * Point[2]);
	mi_Psi.at(10) = -4. * Point[0] * Point[1] * (-1. + Point[2] * Point[2]);
	mi_Psi.at(11) = Point[1] * (2. * Point[1] - 1.) * (1. - Point[2] * Point[2]);
	mi_Psi.at(12) = -0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. + Point[2]) * Point[2];
	mi_Psi.at(13) = -2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. + Point[2]) * Point[2];
	mi_Psi.at(14) = 0.5 * Point[0] * (2. * Point[0] - 1.) * (1. + Point[2]) * Point[2];
	mi_Psi.at(15) = -2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. + Point[2]) * Point[2];
	mi_Psi.at(16) = 2. * Point[0] * Point[1] * (1. + Point[2]) * Point[2];
	mi_Psi.at(17) = 0.5 * Point[1] * (2. * Point[1] - 1.) * (1. + Point[2]) * Point[2];

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Geom::Elem::Elem_Pri18::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(18 * 3);

	mi_DPsi.at( 0) = 0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. - Point[2]) * Point[2];
	mi_DPsi.at( 1) = (-2. + 4. * Point[0] + 2. * Point[1]) * (1. - Point[2]) * Point[2];
	mi_DPsi.at( 2) = (0.5 - 2. * Point[0]) * (1. - Point[2]) * Point[2];
	mi_DPsi.at( 3) = 2. * Point[1] * (1. - Point[2]) * Point[2];
	mi_DPsi.at( 4) = -2. * Point[1] * Point[2] * (1. - Point[2]);
	mi_DPsi.at( 5) = 0.;
	mi_DPsi.at( 6) = (3. - 4. * Point[0] - 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	mi_DPsi.at( 7) = (-4. + 8. * Point[0] + 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	mi_DPsi.at( 8) = (1. - 4. * Point[0]) * (-1. + Point[2] * Point[2]);
	mi_DPsi.at( 9) =  4. * Point[1] * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(10) = -4. * Point[1] * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(11) = 0.;
	mi_DPsi.at(12) = -0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. + Point[2]) * Point[2];
	mi_DPsi.at(13) = (2. - 4. * Point[0] - 2. * Point[1]) * (1. + Point[2]) * Point[2];
	mi_DPsi.at(14) = (-0.5 + 2. * Point[0]) * (1. + Point[2]) * Point[2];
	mi_DPsi.at(15) = -2. * Point[1] * (1. + Point[2]) * Point[2];
	mi_DPsi.at(16) = 2. * Point[1] * Point[2] * (1. + Point[2]);
	mi_DPsi.at(17) = 0.;

	mi_DPsi.at(18) = 0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. - Point[2]) * Point[2];
	mi_DPsi.at(19) = 2. * Point[0] * (1. - Point[2]) * Point[2];
	mi_DPsi.at(20) = 0.;
	mi_DPsi.at(21) = (-2. + 2. * Point[0] + 4. * Point[1]) * (1. - Point[2]) * Point[2];
	mi_DPsi.at(22) = -2. * Point[0] * Point[2] * (1. - Point[2]);
	mi_DPsi.at(23) = (0.5 - 2. * Point[1]) * (1. - Point[2]) * Point[2];
	mi_DPsi.at(24) = (3. - 4. * Point[0] - 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(25) = 4. * Point[0] * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(26) = 0.;
	mi_DPsi.at(27) =  (-4. + 4. * Point[0] + 8. * Point[1]) * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(28) = -4. * Point[0] * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(29) = (1. - 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	mi_DPsi.at(30) = -0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. + Point[2]) * Point[2];
	mi_DPsi.at(31) = -2. * Point[0] * (1. + Point[2]) * Point[2];
	mi_DPsi.at(32) = 0.;
	mi_DPsi.at(33) = (2. - 2. * Point[0] - 4. * Point[1]) * (1. + Point[2]) * Point[2];
	mi_DPsi.at(34) = 2. * Point[0] * Point[2] * (1. + Point[2]);
	mi_DPsi.at(35) = (-0.5 + 2. * Point[1]) * (1. + Point[2]) * Point[2];

	mi_DPsi.at(36) = 0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. - 2. * Point[2]);
	mi_DPsi.at(37) = 2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[2]);
	mi_DPsi.at(38) = -0.5 * Point[0] * (2. * Point[0] - 1.) * (1. - 2. * Point[2]);
	mi_DPsi.at(39) = 2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[2]);
	mi_DPsi.at(40) = -2. * Point[0] * Point[1] * (1. - 2. * Point[2]);
	mi_DPsi.at(41) = -0.5 * Point[1] * (2. * Point[1] - 1.) * (1. - 2. * Point[2]);
	mi_DPsi.at(42) = (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (2. * Point[2]);
	mi_DPsi.at(43) = 4. * Point[0] * (-1. + Point[0] + Point[1]) * (2. * Point[2]);
	mi_DPsi.at(44) = Point[0] * (2. * Point[0] - 1.) * (-2. * Point[2]);
	mi_DPsi.at(45) =  4. * Point[1] * (-1. + Point[0] + Point[1]) * (2. * Point[2]);
	mi_DPsi.at(46) = -4. * Point[0] * Point[1] * (2. * Point[2]);
	mi_DPsi.at(47) = Point[1] * (2. * Point[1] - 1.) * (-2. * Point[2]);
	mi_DPsi.at(48) = -0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. + 2. * Point[2]);
	mi_DPsi.at(49) = -2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. + 2. * Point[2]);
	mi_DPsi.at(50) = 0.5 * Point[0] * (2. * Point[0] - 1.) * (1. + 2. * Point[2]);
	mi_DPsi.at(51) = -2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. + 2. * Point[2]);
	mi_DPsi.at(52) = 2. * Point[0] * Point[1] * (1. + 2. * Point[2]);
	mi_DPsi.at(53) = 0.5 * Point[1] * (2. * Point[1] - 1.) * (1. + 2. * Point[2]);

	return mi_DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Geom::Elem::Elem_Pri18::setGeomProperties() {

	const int nVertices = 6;
	const int mi_Dim = mv_Conect.at(0)->getDIM();	// Dimensionality of vector space (2D or 3D)

	mv_Centroid = std::make_unique<double[]>(mi_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Geom::Node*, nVertices> vertices;
	vertices[0] = mv_Conect[0].get();
	vertices[1] = mv_Conect[2].get();
	vertices[2] = mv_Conect[5].get();
	vertices[3] = mv_Conect[12].get();
	vertices[4] = mv_Conect[14].get();
	vertices[5] = mv_Conect[17].get();

	// Memory requested by make_unique is not empty
	for (int i = 0; i < mi_Dim; i++) mv_Centroid[i] = 0.;

	for (auto& node : vertices) {
		for (int i = 0; i < mi_Dim; i++) mv_Centroid[i] += node->getInitPos()[i];
	}

	// Finishing up
	for (int i = 0; i < mi_Dim; i++) mv_Centroid[i] /= nVertices;

	// Distance from centroid to vertices
	double dist[nVertices] = {};
	int i = 0;

	for (auto& node : vertices) {
		for (int j = 0; j < mi_Dim; j++) {
			dist[i] += (mv_Centroid[j] - node->getInitPos()[j]) * (mv_Centroid[j] - node->getInitPos()[j]);
		}
		dist[i] = std::sqrt(dist[i]);
		i++;
	}

	// Since centroid is not the circumcenter, the radius is related to the minimum bounding circle
	mv_Radius = *std::max_element(dist, dist + nVertices);
};


// ================================================================================================
//
// Implementation of Member Function: getValueOnIPs
// Return the values on the integration points currently known in the element' nodes
// 
// ================================================================================================
inline std::vector<double> O2P2::Geom::Elem::Elem_Pri18::getValueOnIPs(const double* value) {

	// return value
	std::vector<double> mi_valueOnIp(mv_numIP, 0.);

	for (int i = 0; i < mv_numIP; i++) {
		for (int j = 0; j < this->mv_numNodes; j++) {
			mi_valueOnIp.at(i) += value[i] * mv_Psi[i][j];
		}
	}

	return mi_valueOnIp;
};


// ================================================================================================
//
// Integration Points
//
// ================================================================================================
inline const double O2P2::Geom::Elem::Elem_Pri18::mv_xsi[mv_numIP][mv_ElDim] =
	{ { 0.816847572980459, 0.091576213509771, -0.774596669241483 },
	  { 0.091576213509771, 0.816847572980459, -0.774596669241483 },
	  { 0.091576213509771, 0.091576213509771, -0.774596669241483 },
	  { 0.816847572980459, 0.091576213509771,  0.774596669241483 },
	  { 0.091576213509771, 0.816847572980459,  0.774596669241483 },
	  { 0.091576213509771, 0.091576213509771,  0.774596669241483 },
	  { 0.108103018168070, 0.445948490915965, -0.774596669241483 },
	  { 0.445948490915965, 0.108103018168070, -0.774596669241483 },
	  { 0.445948490915965, 0.445948490915965, -0.774596669241483 },
	  { 0.816847572980459, 0.091576213509771,  0. },
	  { 0.091576213509771, 0.816847572980459,  0. },
	  { 0.091576213509771, 0.091576213509771,  0. },
	  { 0.108103018168070, 0.445948490915965,  0.774596669241483 },
	  { 0.445948490915965, 0.108103018168070,  0.774596669241483 },
	  { 0.445948490915965, 0.445948490915965,  0.774596669241483 },
	  { 0.108103018168070, 0.445948490915965,  0. },
	  { 0.445948490915965, 0.108103018168070,  0. },
	  { 0.445948490915965, 0.445948490915965,  0. } };

// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
inline const double O2P2::Geom::Elem::Elem_Pri18::mv_weight[mv_numIP] = { 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.06205044157722527777778, 0.06205044157722527777778, 0.06205044157722527777778, 0.04886744162458755555556, 0.04886744162458755555556, 0.04886744162458755555556, 0.06205044157722527777778, 0.06205044157722527777778, 0.06205044157722527777778, 0.09928070652356044444444, 0.09928070652356044444444, 0.09928070652356044444444 };


// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double O2P2::Geom::Elem::Elem_Pri18::mv_Psi[mv_numIP][mv_numNodes] = {
	{ 0.5 * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2],
	  2. * mv_xsi[0][0] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2],
	  -0.5 * mv_xsi[0][0] * (2. * mv_xsi[0][0] - 1.) * (1. - mv_xsi[0][2]) * mv_xsi[0][2],
	  2. * mv_xsi[0][1] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2],
	  -2. * mv_xsi[0][0] * mv_xsi[0][1] * mv_xsi[0][2] * (1. - mv_xsi[0][2]),
	  -0.5 * mv_xsi[0][1] * (2. * mv_xsi[0][1] - 1.) * (1. - mv_xsi[0][2]) * mv_xsi[0][2],
	  (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]),
	  4. * mv_xsi[0][0] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]),
	  mv_xsi[0][0] * (2. * mv_xsi[0][0] - 1.) * (1. - mv_xsi[0][2] * mv_xsi[0][2]),
	  4. * mv_xsi[0][1] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]),
	  -4. * mv_xsi[0][0] * mv_xsi[0][1] * (-1. + mv_xsi[0][2] * mv_xsi[0][2]),
	  mv_xsi[0][1] * (2. * mv_xsi[0][1] - 1.) * (1. - mv_xsi[0][2] * mv_xsi[0][2]),
	  -0.5 * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2],
	  -2. * mv_xsi[0][0] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2],
	  0.5 * mv_xsi[0][0] * (2. * mv_xsi[0][0] - 1.) * (1. + mv_xsi[0][2]) * mv_xsi[0][2],
	  -2. * mv_xsi[0][1] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2],
	  2. * mv_xsi[0][0] * mv_xsi[0][1] * (1. + mv_xsi[0][2]) * mv_xsi[0][2],
	  0.5 * mv_xsi[0][1] * (2. * mv_xsi[0][1] - 1.) * (1. + mv_xsi[0][2]) * mv_xsi[0][2] },

	{ 0.5 * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2],
	  2. * mv_xsi[1][0] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2],
	  -0.5 * mv_xsi[1][0] * (2. * mv_xsi[1][0] - 1.) * (1. - mv_xsi[1][2]) * mv_xsi[1][2],
	  2. * mv_xsi[1][1] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2],
	  -2. * mv_xsi[1][0] * mv_xsi[1][1] * mv_xsi[1][2] * (1. - mv_xsi[1][2]),
	  -0.5 * mv_xsi[1][1] * (2. * mv_xsi[1][1] - 1.) * (1. - mv_xsi[1][2]) * mv_xsi[1][2],
	  (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]),
	  4. * mv_xsi[1][0] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]),
	  mv_xsi[1][0] * (2. * mv_xsi[1][0] - 1.) * (1. - mv_xsi[1][2] * mv_xsi[1][2]),
	  4. * mv_xsi[1][1] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]),
	  -4. * mv_xsi[1][0] * mv_xsi[1][1] * (-1. + mv_xsi[1][2] * mv_xsi[1][2]),
	  mv_xsi[1][1] * (2. * mv_xsi[1][1] - 1.) * (1. - mv_xsi[1][2] * mv_xsi[1][2]),
	  -0.5 * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2],
	  -2. * mv_xsi[1][0] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2],
	  0.5 * mv_xsi[1][0] * (2. * mv_xsi[1][0] - 1.) * (1. + mv_xsi[1][2]) * mv_xsi[1][2],
	  -2. * mv_xsi[1][1] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2],
	  2. * mv_xsi[1][0] * mv_xsi[1][1] * (1. + mv_xsi[1][2]) * mv_xsi[1][2],
	  0.5 * mv_xsi[1][1] * (2. * mv_xsi[1][1] - 1.) * (1. + mv_xsi[1][2]) * mv_xsi[1][2] },

	{ 0.5 * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2],
	  2. * mv_xsi[2][0] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2],
	  -0.5 * mv_xsi[2][0] * (2. * mv_xsi[2][0] - 1.) * (1. - mv_xsi[2][2]) * mv_xsi[2][2],
	  2. * mv_xsi[2][1] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2],
	  -2. * mv_xsi[2][0] * mv_xsi[2][1] * mv_xsi[2][2] * (1. - mv_xsi[2][2]),
	  -0.5 * mv_xsi[2][1] * (2. * mv_xsi[2][1] - 1.) * (1. - mv_xsi[2][2]) * mv_xsi[2][2],
	  (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]),
	  4. * mv_xsi[2][0] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]),
	  mv_xsi[2][0] * (2. * mv_xsi[2][0] - 1.) * (1. - mv_xsi[2][2] * mv_xsi[2][2]),
	  4. * mv_xsi[2][1] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]),
	  -4. * mv_xsi[2][0] * mv_xsi[2][1] * (-1. + mv_xsi[2][2] * mv_xsi[2][2]),
	  mv_xsi[2][1] * (2. * mv_xsi[2][1] - 1.) * (1. - mv_xsi[2][2] * mv_xsi[2][2]),
	  -0.5 * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2],
	  -2. * mv_xsi[2][0] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2],
	  0.5 * mv_xsi[2][0] * (2. * mv_xsi[2][0] - 1.) * (1. + mv_xsi[2][2]) * mv_xsi[2][2],
	  -2. * mv_xsi[2][1] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2],
	  2. * mv_xsi[2][0] * mv_xsi[2][1] * (1. + mv_xsi[2][2]) * mv_xsi[2][2],
	  0.5 * mv_xsi[2][1] * (2. * mv_xsi[2][1] - 1.) * (1. + mv_xsi[2][2]) * mv_xsi[2][2] },

	{ 0.5 * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2],
	  2. * mv_xsi[3][0] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2],
	  -0.5 * mv_xsi[3][0] * (2. * mv_xsi[3][0] - 1.) * (1. - mv_xsi[3][2]) * mv_xsi[3][2],
	  2. * mv_xsi[3][1] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2],
	  -2. * mv_xsi[3][0] * mv_xsi[3][1] * mv_xsi[3][2] * (1. - mv_xsi[3][2]),
	  -0.5 * mv_xsi[3][1] * (2. * mv_xsi[3][1] - 1.) * (1. - mv_xsi[3][2]) * mv_xsi[3][2],
	  (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]),
	  4. * mv_xsi[3][0] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]),
	  mv_xsi[3][0] * (2. * mv_xsi[3][0] - 1.) * (1. - mv_xsi[3][2] * mv_xsi[3][2]),
	  4. * mv_xsi[3][1] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]),
	  -4. * mv_xsi[3][0] * mv_xsi[3][1] * (-1. + mv_xsi[3][2] * mv_xsi[3][2]),
	  mv_xsi[3][1] * (2. * mv_xsi[3][1] - 1.) * (1. - mv_xsi[3][2] * mv_xsi[3][2]),
	  -0.5 * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2],
	  -2. * mv_xsi[3][0] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2],
	  0.5 * mv_xsi[3][0] * (2. * mv_xsi[3][0] - 1.) * (1. + mv_xsi[3][2]) * mv_xsi[3][2],
	  -2. * mv_xsi[3][1] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2],
	  2. * mv_xsi[3][0] * mv_xsi[3][1] * (1. + mv_xsi[3][2]) * mv_xsi[3][2],
	  0.5 * mv_xsi[3][1] * (2. * mv_xsi[3][1] - 1.) * (1. + mv_xsi[3][2]) * mv_xsi[3][2] },

	{ 0.5 * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2],
	  2. * mv_xsi[4][0] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2],
	  -0.5 * mv_xsi[4][0] * (2. * mv_xsi[4][0] - 1.) * (1. - mv_xsi[4][2]) * mv_xsi[4][2],
	  2. * mv_xsi[4][1] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2],
	  -2. * mv_xsi[4][0] * mv_xsi[4][1] * mv_xsi[4][2] * (1. - mv_xsi[4][2]),
	  -0.5 * mv_xsi[4][1] * (2. * mv_xsi[4][1] - 1.) * (1. - mv_xsi[4][2]) * mv_xsi[4][2],
	  (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]),
	  4. * mv_xsi[4][0] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]),
	  mv_xsi[4][0] * (2. * mv_xsi[4][0] - 1.) * (1. - mv_xsi[4][2] * mv_xsi[4][2]),
	  4. * mv_xsi[4][1] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]),
	  -4. * mv_xsi[4][0] * mv_xsi[4][1] * (-1. + mv_xsi[4][2] * mv_xsi[4][2]),
	  mv_xsi[4][1] * (2. * mv_xsi[4][1] - 1.) * (1. - mv_xsi[4][2] * mv_xsi[4][2]),
	  -0.5 * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2],
	  -2. * mv_xsi[4][0] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2],
	  0.5 * mv_xsi[4][0] * (2. * mv_xsi[4][0] - 1.) * (1. + mv_xsi[4][2]) * mv_xsi[4][2],
	  -2. * mv_xsi[4][1] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2],
	  2. * mv_xsi[4][0] * mv_xsi[4][1] * (1. + mv_xsi[4][2]) * mv_xsi[4][2],
	  0.5 * mv_xsi[4][1] * (2. * mv_xsi[4][1] - 1.) * (1. + mv_xsi[4][2]) * mv_xsi[4][2] },

	{ 0.5 * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2],
	  2. * mv_xsi[5][0] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2],
	  -0.5 * mv_xsi[5][0] * (2. * mv_xsi[5][0] - 1.) * (1. - mv_xsi[5][2]) * mv_xsi[5][2],
	  2. * mv_xsi[5][1] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2],
	  -2. * mv_xsi[5][0] * mv_xsi[5][1] * mv_xsi[5][2] * (1. - mv_xsi[5][2]),
	  -0.5 * mv_xsi[5][1] * (2. * mv_xsi[5][1] - 1.) * (1. - mv_xsi[5][2]) * mv_xsi[5][2],
	  (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]),
	  4. * mv_xsi[5][0] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]),
	  mv_xsi[5][0] * (2. * mv_xsi[5][0] - 1.) * (1. - mv_xsi[5][2] * mv_xsi[5][2]),
	  4. * mv_xsi[5][1] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]),
	  -4. * mv_xsi[5][0] * mv_xsi[5][1] * (-1. + mv_xsi[5][2] * mv_xsi[5][2]),
	  mv_xsi[5][1] * (2. * mv_xsi[5][1] - 1.) * (1. - mv_xsi[5][2] * mv_xsi[5][2]),
	  -0.5 * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2],
	  -2. * mv_xsi[5][0] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2],
	  0.5 * mv_xsi[5][0] * (2. * mv_xsi[5][0] - 1.) * (1. + mv_xsi[5][2]) * mv_xsi[5][2],
	  -2. * mv_xsi[5][1] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2],
	  2. * mv_xsi[5][0] * mv_xsi[5][1] * (1. + mv_xsi[5][2]) * mv_xsi[5][2],
	  0.5 * mv_xsi[5][1] * (2. * mv_xsi[5][1] - 1.) * (1. + mv_xsi[5][2]) * mv_xsi[5][2] },

	{ 0.5 * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2],
	  2. * mv_xsi[6][0] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2],
	  -0.5 * mv_xsi[6][0] * (2. * mv_xsi[6][0] - 1.) * (1. - mv_xsi[6][2]) * mv_xsi[6][2],
	  2. * mv_xsi[6][1] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2],
	  -2. * mv_xsi[6][0] * mv_xsi[6][1] * mv_xsi[6][2] * (1. - mv_xsi[6][2]),
	  -0.5 * mv_xsi[6][1] * (2. * mv_xsi[6][1] - 1.) * (1. - mv_xsi[6][2]) * mv_xsi[6][2],
	  (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]),
	  4. * mv_xsi[6][0] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]),
	  mv_xsi[6][0] * (2. * mv_xsi[6][0] - 1.) * (1. - mv_xsi[6][2] * mv_xsi[6][2]),
	  4. * mv_xsi[6][1] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]),
	  -4. * mv_xsi[6][0] * mv_xsi[6][1] * (-1. + mv_xsi[6][2] * mv_xsi[6][2]),
	  mv_xsi[6][1] * (2. * mv_xsi[6][1] - 1.) * (1. - mv_xsi[6][2] * mv_xsi[6][2]),
	  -0.5 * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2],
	  -2. * mv_xsi[6][0] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2],
	  0.5 * mv_xsi[6][0] * (2. * mv_xsi[6][0] - 1.) * (1. + mv_xsi[6][2]) * mv_xsi[6][2],
	  -2. * mv_xsi[6][1] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2],
	  2. * mv_xsi[6][0] * mv_xsi[6][1] * (1. + mv_xsi[6][2]) * mv_xsi[6][2],
	  0.5 * mv_xsi[6][1] * (2. * mv_xsi[6][1] - 1.) * (1. + mv_xsi[6][2]) * mv_xsi[6][2] },

	{ 0.5 * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2],
	  2. * mv_xsi[7][0] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2],
	  -0.5 * mv_xsi[7][0] * (2. * mv_xsi[7][0] - 1.) * (1. - mv_xsi[7][2]) * mv_xsi[7][2],
	  2. * mv_xsi[7][1] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2],
	  -2. * mv_xsi[7][0] * mv_xsi[7][1] * mv_xsi[7][2] * (1. - mv_xsi[7][2]),
	  -0.5 * mv_xsi[7][1] * (2. * mv_xsi[7][1] - 1.) * (1. - mv_xsi[7][2]) * mv_xsi[7][2],
	  (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]),
	  4. * mv_xsi[7][0] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]),
	  mv_xsi[7][0] * (2. * mv_xsi[7][0] - 1.) * (1. - mv_xsi[7][2] * mv_xsi[7][2]),
	  4. * mv_xsi[7][1] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]),
	  -4. * mv_xsi[7][0] * mv_xsi[7][1] * (-1. + mv_xsi[7][2] * mv_xsi[7][2]),
	  mv_xsi[7][1] * (2. * mv_xsi[7][1] - 1.) * (1. - mv_xsi[7][2] * mv_xsi[7][2]),
	  -0.5 * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2],
	  -2. * mv_xsi[7][0] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2],
	  0.5 * mv_xsi[7][0] * (2. * mv_xsi[7][0] - 1.) * (1. + mv_xsi[7][2]) * mv_xsi[7][2],
	  -2. * mv_xsi[7][1] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2],
	  2. * mv_xsi[7][0] * mv_xsi[7][1] * (1. + mv_xsi[7][2]) * mv_xsi[7][2],
	  0.5 * mv_xsi[7][1] * (2. * mv_xsi[7][1] - 1.) * (1. + mv_xsi[7][2]) * mv_xsi[7][2] },

	{ 0.5 * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2],
	  2. * mv_xsi[8][0] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2],
	  -0.5 * mv_xsi[8][0] * (2. * mv_xsi[8][0] - 1.) * (1. - mv_xsi[8][2]) * mv_xsi[8][2],
	  2. * mv_xsi[8][1] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2],
	  -2. * mv_xsi[8][0] * mv_xsi[8][1] * mv_xsi[8][2] * (1. - mv_xsi[8][2]),
	  -0.5 * mv_xsi[8][1] * (2. * mv_xsi[8][1] - 1.) * (1. - mv_xsi[8][2]) * mv_xsi[8][2],
	  (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]),
	  4. * mv_xsi[8][0] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]),
	  mv_xsi[8][0] * (2. * mv_xsi[8][0] - 1.) * (1. - mv_xsi[8][2] * mv_xsi[8][2]),
	  4. * mv_xsi[8][1] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]),
	  -4. * mv_xsi[8][0] * mv_xsi[8][1] * (-1. + mv_xsi[8][2] * mv_xsi[8][2]),
	  mv_xsi[8][1] * (2. * mv_xsi[8][1] - 1.) * (1. - mv_xsi[8][2] * mv_xsi[8][2]),
	  -0.5 * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2],
	  -2. * mv_xsi[8][0] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2],
	  0.5 * mv_xsi[8][0] * (2. * mv_xsi[8][0] - 1.) * (1. + mv_xsi[8][2]) * mv_xsi[8][2],
	  -2. * mv_xsi[8][1] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2],
	  2. * mv_xsi[8][0] * mv_xsi[8][1] * (1. + mv_xsi[8][2]) * mv_xsi[8][2],
	  0.5 * mv_xsi[8][1] * (2. * mv_xsi[8][1] - 1.) * (1. + mv_xsi[8][2]) * mv_xsi[8][2] },

	{ 0.5 * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2],
	  2. * mv_xsi[9][0] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2],
	  -0.5 * mv_xsi[9][0] * (2. * mv_xsi[9][0] - 1.) * (1. - mv_xsi[9][2]) * mv_xsi[9][2],
	  2. * mv_xsi[9][1] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2],
	  -2. * mv_xsi[9][0] * mv_xsi[9][1] * mv_xsi[9][2] * (1. - mv_xsi[9][2]),
	  -0.5 * mv_xsi[9][1] * (2. * mv_xsi[9][1] - 1.) * (1. - mv_xsi[9][2]) * mv_xsi[9][2],
	  (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]),
	  4. * mv_xsi[9][0] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]),
	  mv_xsi[9][0] * (2. * mv_xsi[9][0] - 1.) * (1. - mv_xsi[9][2] * mv_xsi[9][2]),
	  4. * mv_xsi[9][1] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]),
	  -4. * mv_xsi[9][0] * mv_xsi[9][1] * (-1. + mv_xsi[9][2] * mv_xsi[9][2]),
	  mv_xsi[9][1] * (2. * mv_xsi[9][1] - 1.) * (1. - mv_xsi[9][2] * mv_xsi[9][2]),
	  -0.5 * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2],
	  -2. * mv_xsi[9][0] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2],
	  0.5 * mv_xsi[9][0] * (2. * mv_xsi[9][0] - 1.) * (1. + mv_xsi[9][2]) * mv_xsi[9][2],
	  -2. * mv_xsi[9][1] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2],
	  2. * mv_xsi[9][0] * mv_xsi[9][1] * (1. + mv_xsi[9][2]) * mv_xsi[9][2],
	  0.5 * mv_xsi[9][1] * (2. * mv_xsi[9][1] - 1.) * (1. + mv_xsi[9][2]) * mv_xsi[9][2] },

	{ 0.5 * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2],
	  2. * mv_xsi[10][0] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2],
	  -0.5 * mv_xsi[10][0] * (2. * mv_xsi[10][0] - 1.) * (1. - mv_xsi[10][2]) * mv_xsi[10][2],
	  2. * mv_xsi[10][1] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2],
	  -2. * mv_xsi[10][0] * mv_xsi[10][1] * mv_xsi[10][2] * (1. - mv_xsi[10][2]),
	  -0.5 * mv_xsi[10][1] * (2. * mv_xsi[10][1] - 1.) * (1. - mv_xsi[10][2]) * mv_xsi[10][2],
	  (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]),
	  4. * mv_xsi[10][0] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]),
	  mv_xsi[10][0] * (2. * mv_xsi[10][0] - 1.) * (1. - mv_xsi[10][2] * mv_xsi[10][2]),
	  4. * mv_xsi[10][1] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]),
	  -4. * mv_xsi[10][0] * mv_xsi[10][1] * (-1. + mv_xsi[10][2] * mv_xsi[10][2]),
	  mv_xsi[10][1] * (2. * mv_xsi[10][1] - 1.) * (1. - mv_xsi[10][2] * mv_xsi[10][2]),
	  -0.5 * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2],
	  -2. * mv_xsi[10][0] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2],
	  0.5 * mv_xsi[10][0] * (2. * mv_xsi[10][0] - 1.) * (1. + mv_xsi[10][2]) * mv_xsi[10][2],
	  -2. * mv_xsi[10][1] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2],
	  2. * mv_xsi[10][0] * mv_xsi[10][1] * (1. + mv_xsi[10][2]) * mv_xsi[10][2],
	  0.5 * mv_xsi[10][1] * (2. * mv_xsi[10][1] - 1.) * (1. + mv_xsi[10][2]) * mv_xsi[10][2] },

	{ 0.5 * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2],
	  2. * mv_xsi[11][0] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2],
	  -0.5 * mv_xsi[11][0] * (2. * mv_xsi[11][0] - 1.) * (1. - mv_xsi[11][2]) * mv_xsi[11][2],
	  2. * mv_xsi[11][1] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2],
	  -2. * mv_xsi[11][0] * mv_xsi[11][1] * mv_xsi[11][2] * (1. - mv_xsi[11][2]),
	  -0.5 * mv_xsi[11][1] * (2. * mv_xsi[11][1] - 1.) * (1. - mv_xsi[11][2]) * mv_xsi[11][2],
	  (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]),
	  4. * mv_xsi[11][0] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]),
	  mv_xsi[11][0] * (2. * mv_xsi[11][0] - 1.) * (1. - mv_xsi[11][2] * mv_xsi[11][2]),
	  4. * mv_xsi[11][1] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]),
	  -4. * mv_xsi[11][0] * mv_xsi[11][1] * (-1. + mv_xsi[11][2] * mv_xsi[11][2]),
	  mv_xsi[11][1] * (2. * mv_xsi[11][1] - 1.) * (1. - mv_xsi[11][2] * mv_xsi[11][2]),
	  -0.5 * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2],
	  -2. * mv_xsi[11][0] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2],
	  0.5 * mv_xsi[11][0] * (2. * mv_xsi[11][0] - 1.) * (1. + mv_xsi[11][2]) * mv_xsi[11][2],
	  -2. * mv_xsi[11][1] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2],
	  2. * mv_xsi[11][0] * mv_xsi[11][1] * (1. + mv_xsi[11][2]) * mv_xsi[11][2],
	  0.5 * mv_xsi[11][1] * (2. * mv_xsi[11][1] - 1.) * (1. + mv_xsi[11][2]) * mv_xsi[11][2] },

	{ 0.5 * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2],
	  2. * mv_xsi[12][0] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2],
	  -0.5 * mv_xsi[12][0] * (2. * mv_xsi[12][0] - 1.) * (1. - mv_xsi[12][2]) * mv_xsi[12][2],
	  2. * mv_xsi[12][1] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2],
	  -2. * mv_xsi[12][0] * mv_xsi[12][1] * mv_xsi[12][2] * (1. - mv_xsi[12][2]),
	  -0.5 * mv_xsi[12][1] * (2. * mv_xsi[12][1] - 1.) * (1. - mv_xsi[12][2]) * mv_xsi[12][2],
	  (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]),
	  4. * mv_xsi[12][0] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]),
	  mv_xsi[12][0] * (2. * mv_xsi[12][0] - 1.) * (1. - mv_xsi[12][2] * mv_xsi[12][2]),
	  4. * mv_xsi[12][1] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]),
	  -4. * mv_xsi[12][0] * mv_xsi[12][1] * (-1. + mv_xsi[12][2] * mv_xsi[12][2]),
	  mv_xsi[12][1] * (2. * mv_xsi[12][1] - 1.) * (1. - mv_xsi[12][2] * mv_xsi[12][2]),
	  -0.5 * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2],
	  -2. * mv_xsi[12][0] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2],
	  0.5 * mv_xsi[12][0] * (2. * mv_xsi[12][0] - 1.) * (1. + mv_xsi[12][2]) * mv_xsi[12][2],
	  -2. * mv_xsi[12][1] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2],
	  2. * mv_xsi[12][0] * mv_xsi[12][1] * (1. + mv_xsi[12][2]) * mv_xsi[12][2],
	  0.5 * mv_xsi[12][1] * (2. * mv_xsi[12][1] - 1.) * (1. + mv_xsi[12][2]) * mv_xsi[12][2] },

	{ 0.5 * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2],
	  2. * mv_xsi[13][0] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2],
	  -0.5 * mv_xsi[13][0] * (2. * mv_xsi[13][0] - 1.) * (1. - mv_xsi[13][2]) * mv_xsi[13][2],
	  2. * mv_xsi[13][1] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2],
	  -2. * mv_xsi[13][0] * mv_xsi[13][1] * mv_xsi[13][2] * (1. - mv_xsi[13][2]),
	  -0.5 * mv_xsi[13][1] * (2. * mv_xsi[13][1] - 1.) * (1. - mv_xsi[13][2]) * mv_xsi[13][2],
	  (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]),
	  4. * mv_xsi[13][0] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]),
	  mv_xsi[13][0] * (2. * mv_xsi[13][0] - 1.) * (1. - mv_xsi[13][2] * mv_xsi[13][2]),
	  4. * mv_xsi[13][1] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]),
	  -4. * mv_xsi[13][0] * mv_xsi[13][1] * (-1. + mv_xsi[13][2] * mv_xsi[13][2]),
	  mv_xsi[13][1] * (2. * mv_xsi[13][1] - 1.) * (1. - mv_xsi[13][2] * mv_xsi[13][2]),
	  -0.5 * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2],
	  -2. * mv_xsi[13][0] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2],
	  0.5 * mv_xsi[13][0] * (2. * mv_xsi[13][0] - 1.) * (1. + mv_xsi[13][2]) * mv_xsi[13][2],
	  -2. * mv_xsi[13][1] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2],
	  2. * mv_xsi[13][0] * mv_xsi[13][1] * (1. + mv_xsi[13][2]) * mv_xsi[13][2],
	  0.5 * mv_xsi[13][1] * (2. * mv_xsi[13][1] - 1.) * (1. + mv_xsi[13][2]) * mv_xsi[13][2] },

	{ 0.5 * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2],
	  2. * mv_xsi[14][0] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2],
	  -0.5 * mv_xsi[14][0] * (2. * mv_xsi[14][0] - 1.) * (1. - mv_xsi[14][2]) * mv_xsi[14][2],
	  2. * mv_xsi[14][1] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2],
	  -2. * mv_xsi[14][0] * mv_xsi[14][1] * mv_xsi[14][2] * (1. - mv_xsi[14][2]),
	  -0.5 * mv_xsi[14][1] * (2. * mv_xsi[14][1] - 1.) * (1. - mv_xsi[14][2]) * mv_xsi[14][2],
	  (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]),
	  4. * mv_xsi[14][0] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]),
	  mv_xsi[14][0] * (2. * mv_xsi[14][0] - 1.) * (1. - mv_xsi[14][2] * mv_xsi[14][2]),
	  4. * mv_xsi[14][1] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]),
	  -4. * mv_xsi[14][0] * mv_xsi[14][1] * (-1. + mv_xsi[14][2] * mv_xsi[14][2]),
	  mv_xsi[14][1] * (2. * mv_xsi[14][1] - 1.) * (1. - mv_xsi[14][2] * mv_xsi[14][2]),
	  -0.5 * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2],
	  -2. * mv_xsi[14][0] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2],
	  0.5 * mv_xsi[14][0] * (2. * mv_xsi[14][0] - 1.) * (1. + mv_xsi[14][2]) * mv_xsi[14][2],
	  -2. * mv_xsi[14][1] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2],
	  2. * mv_xsi[14][0] * mv_xsi[14][1] * (1. + mv_xsi[14][2]) * mv_xsi[14][2],
	  0.5 * mv_xsi[14][1] * (2. * mv_xsi[14][1] - 1.) * (1. + mv_xsi[14][2]) * mv_xsi[14][2] },

	{ 0.5 * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2],
	  2. * mv_xsi[15][0] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2],
	  -0.5 * mv_xsi[15][0] * (2. * mv_xsi[15][0] - 1.) * (1. - mv_xsi[15][2]) * mv_xsi[15][2],
	  2. * mv_xsi[15][1] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2],
	  -2. * mv_xsi[15][0] * mv_xsi[15][1] * mv_xsi[15][2] * (1. - mv_xsi[15][2]),
	  -0.5 * mv_xsi[15][1] * (2. * mv_xsi[15][1] - 1.) * (1. - mv_xsi[15][2]) * mv_xsi[15][2],
	  (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]),
	  4. * mv_xsi[15][0] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]),
	  mv_xsi[15][0] * (2. * mv_xsi[15][0] - 1.) * (1. - mv_xsi[15][2] * mv_xsi[15][2]),
	  4. * mv_xsi[15][1] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]),
	  -4. * mv_xsi[15][0] * mv_xsi[15][1] * (-1. + mv_xsi[15][2] * mv_xsi[15][2]),
	  mv_xsi[15][1] * (2. * mv_xsi[15][1] - 1.) * (1. - mv_xsi[15][2] * mv_xsi[15][2]),
	  -0.5 * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2],
	  -2. * mv_xsi[15][0] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2],
	  0.5 * mv_xsi[15][0] * (2. * mv_xsi[15][0] - 1.) * (1. + mv_xsi[15][2]) * mv_xsi[15][2],
	  -2. * mv_xsi[15][1] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2],
	  2. * mv_xsi[15][0] * mv_xsi[15][1] * (1. + mv_xsi[15][2]) * mv_xsi[15][2],
	  0.5 * mv_xsi[15][1] * (2. * mv_xsi[15][1] - 1.) * (1. + mv_xsi[15][2]) * mv_xsi[15][2] },

	{ 0.5 * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2],
	  2. * mv_xsi[16][0] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2],
	  -0.5 * mv_xsi[16][0] * (2. * mv_xsi[16][0] - 1.) * (1. - mv_xsi[16][2]) * mv_xsi[16][2],
	  2. * mv_xsi[16][1] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2],
	  -2. * mv_xsi[16][0] * mv_xsi[16][1] * mv_xsi[16][2] * (1. - mv_xsi[16][2]),
	  -0.5 * mv_xsi[16][1] * (2. * mv_xsi[16][1] - 1.) * (1. - mv_xsi[16][2]) * mv_xsi[16][2],
	  (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]),
	  4. * mv_xsi[16][0] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]),
	  mv_xsi[16][0] * (2. * mv_xsi[16][0] - 1.) * (1. - mv_xsi[16][2] * mv_xsi[16][2]),
	  4. * mv_xsi[16][1] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]),
	  -4. * mv_xsi[16][0] * mv_xsi[16][1] * (-1. + mv_xsi[16][2] * mv_xsi[16][2]),
	  mv_xsi[16][1] * (2. * mv_xsi[16][1] - 1.) * (1. - mv_xsi[16][2] * mv_xsi[16][2]),
	  -0.5 * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2],
	  -2. * mv_xsi[16][0] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2],
	  0.5 * mv_xsi[16][0] * (2. * mv_xsi[16][0] - 1.) * (1. + mv_xsi[16][2]) * mv_xsi[16][2],
	  -2. * mv_xsi[16][1] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2],
	  2. * mv_xsi[16][0] * mv_xsi[16][1] * (1. + mv_xsi[16][2]) * mv_xsi[16][2],
	  0.5 * mv_xsi[16][1] * (2. * mv_xsi[16][1] - 1.) * (1. + mv_xsi[16][2]) * mv_xsi[16][2] },

	{ 0.5 * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2],
	  2. * mv_xsi[17][0] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2],
	  -0.5 * mv_xsi[17][0] * (2. * mv_xsi[17][0] - 1.) * (1. - mv_xsi[17][2]) * mv_xsi[17][2],
	  2. * mv_xsi[17][1] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2],
	  -2. * mv_xsi[17][0] * mv_xsi[17][1] * mv_xsi[17][2] * (1. - mv_xsi[17][2]),
	  -0.5 * mv_xsi[17][1] * (2. * mv_xsi[17][1] - 1.) * (1. - mv_xsi[17][2]) * mv_xsi[17][2],
	  (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]),
	  4. * mv_xsi[17][0] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]),
	  mv_xsi[17][0] * (2. * mv_xsi[17][0] - 1.) * (1. - mv_xsi[17][2] * mv_xsi[17][2]),
	  4. * mv_xsi[17][1] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]),
	  -4. * mv_xsi[17][0] * mv_xsi[17][1] * (-1. + mv_xsi[17][2] * mv_xsi[17][2]),
	  mv_xsi[17][1] * (2. * mv_xsi[17][1] - 1.) * (1. - mv_xsi[17][2] * mv_xsi[17][2]),
	  -0.5 * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2],
	  -2. * mv_xsi[17][0] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2],
	  0.5 * mv_xsi[17][0] * (2. * mv_xsi[17][0] - 1.) * (1. + mv_xsi[17][2]) * mv_xsi[17][2],
	  -2. * mv_xsi[17][1] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2],
	  2. * mv_xsi[17][0] * mv_xsi[17][1] * (1. + mv_xsi[17][2]) * mv_xsi[17][2],
	  0.5 * mv_xsi[17][1] * (2. * mv_xsi[17][1] - 1.) * (1. + mv_xsi[17][2]) * mv_xsi[17][2] } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
inline const double O2P2::Geom::Elem::Elem_Pri18::mv_DPsi[mv_numIP][mv_numNodes][mv_ElDim] = {
	{ { 0.5 * (3. - 4. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2], 0.5 * (3. - 4. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2], 0.5 * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][2]) },
	  { (-2. + 4. * mv_xsi[0][0] + 2. * mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2], 2. * mv_xsi[0][0] * (1. - mv_xsi[0][2]) * mv_xsi[0][2], 2. * mv_xsi[0][0] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][2]) },
	  { (0.5 - 2. * mv_xsi[0][0]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2], 0., -0.5 * mv_xsi[0][0] * (2. * mv_xsi[0][0] - 1.) * (1. - 2. * mv_xsi[0][2]) },
	  { 2. * mv_xsi[0][1] * (1. - mv_xsi[0][2]) * mv_xsi[0][2], (-2. + 2. * mv_xsi[0][0] + 4. * mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2], 2. * mv_xsi[0][1] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][2]) },
	  { -2. * mv_xsi[0][1] * mv_xsi[0][2] * (1. - mv_xsi[0][2]), -2. * mv_xsi[0][0] * mv_xsi[0][2] * (1. - mv_xsi[0][2]), -2. * mv_xsi[0][0] * mv_xsi[0][1] * (1. - 2. * mv_xsi[0][2]) },
	  { 0., (0.5 - 2. * mv_xsi[0][1]) * (1. - mv_xsi[0][2]) * mv_xsi[0][2], -0.5 * mv_xsi[0][1] * (2. * mv_xsi[0][1] - 1.) * (1. - 2. * mv_xsi[0][2]) },
	  { (3. - 4. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), (3. - 4. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (2. * mv_xsi[0][2]) },
	  { (-4. + 8. * mv_xsi[0][0] + 4. * mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), 4. * mv_xsi[0][0] * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), 4. * mv_xsi[0][0] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (2. * mv_xsi[0][2]) },
	  { (1. - 4. * mv_xsi[0][0]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), 0., mv_xsi[0][0] * (2. * mv_xsi[0][0] - 1.) * (-2. * mv_xsi[0][2]) },
	  { 4. * mv_xsi[0][1] * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), (-4. + 4. * mv_xsi[0][0] + 8. * mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), 4. * mv_xsi[0][1] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (2. * mv_xsi[0][2]) },
	  { -4. * mv_xsi[0][1] * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), -4. * mv_xsi[0][0] * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), -4. * mv_xsi[0][0] * mv_xsi[0][1] * (2. * mv_xsi[0][2]) },
	  { 0., (1. - 4. * mv_xsi[0][1]) * (-1. + mv_xsi[0][2] * mv_xsi[0][2]), mv_xsi[0][1] * (2. * mv_xsi[0][1] - 1.) * (-2. * mv_xsi[0][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2], -0.5 * (3. - 4. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2], -0.5 * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. - 2. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (1. + 2. * mv_xsi[0][2]) },
	  { (2. - 4. * mv_xsi[0][0] - 2. * mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2], -2. * mv_xsi[0][0] * (1. + mv_xsi[0][2]) * mv_xsi[0][2], -2. * mv_xsi[0][0] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. + 2. * mv_xsi[0][2]) },
	  { (-0.5 + 2. * mv_xsi[0][0]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2], 0., 0.5 * mv_xsi[0][0] * (2. * mv_xsi[0][0] - 1.) * (1. + 2. * mv_xsi[0][2]) },
	  { -2. * mv_xsi[0][1] * (1. + mv_xsi[0][2]) * mv_xsi[0][2], (2. - 2. * mv_xsi[0][0] - 4. * mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2], -2. * mv_xsi[0][1] * (-1. + mv_xsi[0][0] + mv_xsi[0][1]) * (1. + 2. * mv_xsi[0][2]) },
	  { 2. * mv_xsi[0][1] * mv_xsi[0][2] * (1. + mv_xsi[0][2]), 2. * mv_xsi[0][0] * mv_xsi[0][2] * (1. + mv_xsi[0][2]), 2. * mv_xsi[0][0] * mv_xsi[0][1] * (1. + 2. * mv_xsi[0][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[0][1]) * (1. + mv_xsi[0][2]) * mv_xsi[0][2], 0.5 * mv_xsi[0][1] * (2. * mv_xsi[0][1] - 1.) * (1. + 2. * mv_xsi[0][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2], 0.5 * (3. - 4. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2], 0.5 * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][2]) },
	  { (-2. + 4. * mv_xsi[1][0] + 2. * mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2], 2. * mv_xsi[1][0] * (1. - mv_xsi[1][2]) * mv_xsi[1][2], 2. * mv_xsi[1][0] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][2]) },
	  { (0.5 - 2. * mv_xsi[1][0]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2], 0., -0.5 * mv_xsi[1][0] * (2. * mv_xsi[1][0] - 1.) * (1. - 2. * mv_xsi[1][2]) },
	  { 2. * mv_xsi[1][1] * (1. - mv_xsi[1][2]) * mv_xsi[1][2], (-2. + 2. * mv_xsi[1][0] + 4. * mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2], 2. * mv_xsi[1][1] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][2]) },
	  { -2. * mv_xsi[1][1] * mv_xsi[1][2] * (1. - mv_xsi[1][2]), -2. * mv_xsi[1][0] * mv_xsi[1][2] * (1. - mv_xsi[1][2]), -2. * mv_xsi[1][0] * mv_xsi[1][1] * (1. - 2. * mv_xsi[1][2]) },
	  { 0., (0.5 - 2. * mv_xsi[1][1]) * (1. - mv_xsi[1][2]) * mv_xsi[1][2], -0.5 * mv_xsi[1][1] * (2. * mv_xsi[1][1] - 1.) * (1. - 2. * mv_xsi[1][2]) },
	  { (3. - 4. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), (3. - 4. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (2. * mv_xsi[1][2]) },
	  { (-4. + 8. * mv_xsi[1][0] + 4. * mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), 4. * mv_xsi[1][0] * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), 4. * mv_xsi[1][0] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (2. * mv_xsi[1][2]) },
	  { (1. - 4. * mv_xsi[1][0]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), 0., mv_xsi[1][0] * (2. * mv_xsi[1][0] - 1.) * (-2. * mv_xsi[1][2]) },
	  { 4. * mv_xsi[1][1] * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), (-4. + 4. * mv_xsi[1][0] + 8. * mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), 4. * mv_xsi[1][1] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (2. * mv_xsi[1][2]) },
	  { -4. * mv_xsi[1][1] * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), -4. * mv_xsi[1][0] * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), -4. * mv_xsi[1][0] * mv_xsi[1][1] * (2. * mv_xsi[1][2]) },
	  { 0., (1. - 4. * mv_xsi[1][1]) * (-1. + mv_xsi[1][2] * mv_xsi[1][2]), mv_xsi[1][1] * (2. * mv_xsi[1][1] - 1.) * (-2. * mv_xsi[1][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2], -0.5 * (3. - 4. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2], -0.5 * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. - 2. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (1. + 2. * mv_xsi[1][2]) },
	  { (2. - 4. * mv_xsi[1][0] - 2. * mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2], -2. * mv_xsi[1][0] * (1. + mv_xsi[1][2]) * mv_xsi[1][2], -2. * mv_xsi[1][0] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. + 2. * mv_xsi[1][2]) },
	  { (-0.5 + 2. * mv_xsi[1][0]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2], 0., 0.5 * mv_xsi[1][0] * (2. * mv_xsi[1][0] - 1.) * (1. + 2. * mv_xsi[1][2]) },
	  { -2. * mv_xsi[1][1] * (1. + mv_xsi[1][2]) * mv_xsi[1][2], (2. - 2. * mv_xsi[1][0] - 4. * mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2], -2. * mv_xsi[1][1] * (-1. + mv_xsi[1][0] + mv_xsi[1][1]) * (1. + 2. * mv_xsi[1][2]) },
	  { 2. * mv_xsi[1][1] * mv_xsi[1][2] * (1. + mv_xsi[1][2]), 2. * mv_xsi[1][0] * mv_xsi[1][2] * (1. + mv_xsi[1][2]), 2. * mv_xsi[1][0] * mv_xsi[1][1] * (1. + 2. * mv_xsi[1][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[1][1]) * (1. + mv_xsi[1][2]) * mv_xsi[1][2], 0.5 * mv_xsi[1][1] * (2. * mv_xsi[1][1] - 1.) * (1. + 2. * mv_xsi[1][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2], 0.5 * (3. - 4. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2], 0.5 * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][2]) },
	  { (-2. + 4. * mv_xsi[2][0] + 2. * mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2], 2. * mv_xsi[2][0] * (1. - mv_xsi[2][2]) * mv_xsi[2][2], 2. * mv_xsi[2][0] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][2]) },
	  { (0.5 - 2. * mv_xsi[2][0]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2], 0., -0.5 * mv_xsi[2][0] * (2. * mv_xsi[2][0] - 1.) * (1. - 2. * mv_xsi[2][2]) },
	  { 2. * mv_xsi[2][1] * (1. - mv_xsi[2][2]) * mv_xsi[2][2], (-2. + 2. * mv_xsi[2][0] + 4. * mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2], 2. * mv_xsi[2][1] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][2]) },
	  { -2. * mv_xsi[2][1] * mv_xsi[2][2] * (1. - mv_xsi[2][2]), -2. * mv_xsi[2][0] * mv_xsi[2][2] * (1. - mv_xsi[2][2]), -2. * mv_xsi[2][0] * mv_xsi[2][1] * (1. - 2. * mv_xsi[2][2]) },
	  { 0., (0.5 - 2. * mv_xsi[2][1]) * (1. - mv_xsi[2][2]) * mv_xsi[2][2], -0.5 * mv_xsi[2][1] * (2. * mv_xsi[2][1] - 1.) * (1. - 2. * mv_xsi[2][2]) },
	  { (3. - 4. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), (3. - 4. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (2. * mv_xsi[2][2]) },
	  { (-4. + 8. * mv_xsi[2][0] + 4. * mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), 4. * mv_xsi[2][0] * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), 4. * mv_xsi[2][0] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (2. * mv_xsi[2][2]) },
	  { (1. - 4. * mv_xsi[2][0]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), 0., mv_xsi[2][0] * (2. * mv_xsi[2][0] - 1.) * (-2. * mv_xsi[2][2]) },
	  { 4. * mv_xsi[2][1] * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), (-4. + 4. * mv_xsi[2][0] + 8. * mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), 4. * mv_xsi[2][1] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (2. * mv_xsi[2][2]) },
	  { -4. * mv_xsi[2][1] * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), -4. * mv_xsi[2][0] * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), -4. * mv_xsi[2][0] * mv_xsi[2][1] * (2. * mv_xsi[2][2]) },
	  { 0., (1. - 4. * mv_xsi[2][1]) * (-1. + mv_xsi[2][2] * mv_xsi[2][2]), mv_xsi[2][1] * (2. * mv_xsi[2][1] - 1.) * (-2. * mv_xsi[2][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2], -0.5 * (3. - 4. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2], -0.5 * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. - 2. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (1. + 2. * mv_xsi[2][2]) },
	  { (2. - 4. * mv_xsi[2][0] - 2. * mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2], -2. * mv_xsi[2][0] * (1. + mv_xsi[2][2]) * mv_xsi[2][2], -2. * mv_xsi[2][0] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. + 2. * mv_xsi[2][2]) },
	  { (-0.5 + 2. * mv_xsi[2][0]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2], 0., 0.5 * mv_xsi[2][0] * (2. * mv_xsi[2][0] - 1.) * (1. + 2. * mv_xsi[2][2]) },
	  { -2. * mv_xsi[2][1] * (1. + mv_xsi[2][2]) * mv_xsi[2][2], (2. - 2. * mv_xsi[2][0] - 4. * mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2], -2. * mv_xsi[2][1] * (-1. + mv_xsi[2][0] + mv_xsi[2][1]) * (1. + 2. * mv_xsi[2][2]) },
	  { 2. * mv_xsi[2][1] * mv_xsi[2][2] * (1. + mv_xsi[2][2]), 2. * mv_xsi[2][0] * mv_xsi[2][2] * (1. + mv_xsi[2][2]), 2. * mv_xsi[2][0] * mv_xsi[2][1] * (1. + 2. * mv_xsi[2][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[2][1]) * (1. + mv_xsi[2][2]) * mv_xsi[2][2], 0.5 * mv_xsi[2][1] * (2. * mv_xsi[2][1] - 1.) * (1. + 2. * mv_xsi[2][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2], 0.5 * (3. - 4. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2], 0.5 * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][2]) },
	  { (-2. + 4. * mv_xsi[3][0] + 2. * mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2], 2. * mv_xsi[3][0] * (1. - mv_xsi[3][2]) * mv_xsi[3][2], 2. * mv_xsi[3][0] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][2]) },
	  { (0.5 - 2. * mv_xsi[3][0]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2], 0., -0.5 * mv_xsi[3][0] * (2. * mv_xsi[3][0] - 1.) * (1. - 2. * mv_xsi[3][2]) },
	  { 2. * mv_xsi[3][1] * (1. - mv_xsi[3][2]) * mv_xsi[3][2], (-2. + 2. * mv_xsi[3][0] + 4. * mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2], 2. * mv_xsi[3][1] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][2]) },
	  { -2. * mv_xsi[3][1] * mv_xsi[3][2] * (1. - mv_xsi[3][2]), -2. * mv_xsi[3][0] * mv_xsi[3][2] * (1. - mv_xsi[3][2]), -2. * mv_xsi[3][0] * mv_xsi[3][1] * (1. - 2. * mv_xsi[3][2]) },
	  { 0., (0.5 - 2. * mv_xsi[3][1]) * (1. - mv_xsi[3][2]) * mv_xsi[3][2], -0.5 * mv_xsi[3][1] * (2. * mv_xsi[3][1] - 1.) * (1. - 2. * mv_xsi[3][2]) },
	  { (3. - 4. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), (3. - 4. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (2. * mv_xsi[3][2]) },
	  { (-4. + 8. * mv_xsi[3][0] + 4. * mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), 4. * mv_xsi[3][0] * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), 4. * mv_xsi[3][0] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (2. * mv_xsi[3][2]) },
	  { (1. - 4. * mv_xsi[3][0]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), 0., mv_xsi[3][0] * (2. * mv_xsi[3][0] - 1.) * (-2. * mv_xsi[3][2]) },
	  { 4. * mv_xsi[3][1] * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), (-4. + 4. * mv_xsi[3][0] + 8. * mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), 4. * mv_xsi[3][1] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (2. * mv_xsi[3][2]) },
	  { -4. * mv_xsi[3][1] * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), -4. * mv_xsi[3][0] * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), -4. * mv_xsi[3][0] * mv_xsi[3][1] * (2. * mv_xsi[3][2]) },
	  { 0., (1. - 4. * mv_xsi[3][1]) * (-1. + mv_xsi[3][2] * mv_xsi[3][2]), mv_xsi[3][1] * (2. * mv_xsi[3][1] - 1.) * (-2. * mv_xsi[3][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2], -0.5 * (3. - 4. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2], -0.5 * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. - 2. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (1. + 2. * mv_xsi[3][2]) },
	  { (2. - 4. * mv_xsi[3][0] - 2. * mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2], -2. * mv_xsi[3][0] * (1. + mv_xsi[3][2]) * mv_xsi[3][2], -2. * mv_xsi[3][0] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. + 2. * mv_xsi[3][2]) },
	  { (-0.5 + 2. * mv_xsi[3][0]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2], 0., 0.5 * mv_xsi[3][0] * (2. * mv_xsi[3][0] - 1.) * (1. + 2. * mv_xsi[3][2]) },
	  { -2. * mv_xsi[3][1] * (1. + mv_xsi[3][2]) * mv_xsi[3][2], (2. - 2. * mv_xsi[3][0] - 4. * mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2], -2. * mv_xsi[3][1] * (-1. + mv_xsi[3][0] + mv_xsi[3][1]) * (1. + 2. * mv_xsi[3][2]) },
	  { 2. * mv_xsi[3][1] * mv_xsi[3][2] * (1. + mv_xsi[3][2]), 2. * mv_xsi[3][0] * mv_xsi[3][2] * (1. + mv_xsi[3][2]), 2. * mv_xsi[3][0] * mv_xsi[3][1] * (1. + 2. * mv_xsi[3][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[3][1]) * (1. + mv_xsi[3][2]) * mv_xsi[3][2], 0.5 * mv_xsi[3][1] * (2. * mv_xsi[3][1] - 1.) * (1. + 2. * mv_xsi[3][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2], 0.5 * (3. - 4. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2], 0.5 * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][2]) },
	  { (-2. + 4. * mv_xsi[4][0] + 2. * mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2], 2. * mv_xsi[4][0] * (1. - mv_xsi[4][2]) * mv_xsi[4][2], 2. * mv_xsi[4][0] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][2]) },
	  { (0.5 - 2. * mv_xsi[4][0]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2], 0., -0.5 * mv_xsi[4][0] * (2. * mv_xsi[4][0] - 1.) * (1. - 2. * mv_xsi[4][2]) },
	  { 2. * mv_xsi[4][1] * (1. - mv_xsi[4][2]) * mv_xsi[4][2], (-2. + 2. * mv_xsi[4][0] + 4. * mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2], 2. * mv_xsi[4][1] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][2]) },
	  { -2. * mv_xsi[4][1] * mv_xsi[4][2] * (1. - mv_xsi[4][2]), -2. * mv_xsi[4][0] * mv_xsi[4][2] * (1. - mv_xsi[4][2]), -2. * mv_xsi[4][0] * mv_xsi[4][1] * (1. - 2. * mv_xsi[4][2]) },
	  { 0., (0.5 - 2. * mv_xsi[4][1]) * (1. - mv_xsi[4][2]) * mv_xsi[4][2], -0.5 * mv_xsi[4][1] * (2. * mv_xsi[4][1] - 1.) * (1. - 2. * mv_xsi[4][2]) },
	  { (3. - 4. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), (3. - 4. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (2. * mv_xsi[4][2]) },
	  { (-4. + 8. * mv_xsi[4][0] + 4. * mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), 4. * mv_xsi[4][0] * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), 4. * mv_xsi[4][0] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (2. * mv_xsi[4][2]) },
	  { (1. - 4. * mv_xsi[4][0]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), 0., mv_xsi[4][0] * (2. * mv_xsi[4][0] - 1.) * (-2. * mv_xsi[4][2]) },
	  { 4. * mv_xsi[4][1] * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), (-4. + 4. * mv_xsi[4][0] + 8. * mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), 4. * mv_xsi[4][1] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (2. * mv_xsi[4][2]) },
	  { -4. * mv_xsi[4][1] * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), -4. * mv_xsi[4][0] * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), -4. * mv_xsi[4][0] * mv_xsi[4][1] * (2. * mv_xsi[4][2]) },
	  { 0., (1. - 4. * mv_xsi[4][1]) * (-1. + mv_xsi[4][2] * mv_xsi[4][2]), mv_xsi[4][1] * (2. * mv_xsi[4][1] - 1.) * (-2. * mv_xsi[4][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2], -0.5 * (3. - 4. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2], -0.5 * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. - 2. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (1. + 2. * mv_xsi[4][2]) },
	  { (2. - 4. * mv_xsi[4][0] - 2. * mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2], -2. * mv_xsi[4][0] * (1. + mv_xsi[4][2]) * mv_xsi[4][2], -2. * mv_xsi[4][0] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. + 2. * mv_xsi[4][2]) },
	  { (-0.5 + 2. * mv_xsi[4][0]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2], 0., 0.5 * mv_xsi[4][0] * (2. * mv_xsi[4][0] - 1.) * (1. + 2. * mv_xsi[4][2]) },
	  { -2. * mv_xsi[4][1] * (1. + mv_xsi[4][2]) * mv_xsi[4][2], (2. - 2. * mv_xsi[4][0] - 4. * mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2], -2. * mv_xsi[4][1] * (-1. + mv_xsi[4][0] + mv_xsi[4][1]) * (1. + 2. * mv_xsi[4][2]) },
	  { 2. * mv_xsi[4][1] * mv_xsi[4][2] * (1. + mv_xsi[4][2]), 2. * mv_xsi[4][0] * mv_xsi[4][2] * (1. + mv_xsi[4][2]), 2. * mv_xsi[4][0] * mv_xsi[4][1] * (1. + 2. * mv_xsi[4][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[4][1]) * (1. + mv_xsi[4][2]) * mv_xsi[4][2], 0.5 * mv_xsi[4][1] * (2. * mv_xsi[4][1] - 1.) * (1. + 2. * mv_xsi[4][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2], 0.5 * (3. - 4. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2], 0.5 * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][2]) },
	  { (-2. + 4. * mv_xsi[5][0] + 2. * mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2], 2. * mv_xsi[5][0] * (1. - mv_xsi[5][2]) * mv_xsi[5][2], 2. * mv_xsi[5][0] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][2]) },
	  { (0.5 - 2. * mv_xsi[5][0]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2], 0., -0.5 * mv_xsi[5][0] * (2. * mv_xsi[5][0] - 1.) * (1. - 2. * mv_xsi[5][2]) },
	  { 2. * mv_xsi[5][1] * (1. - mv_xsi[5][2]) * mv_xsi[5][2], (-2. + 2. * mv_xsi[5][0] + 4. * mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2], 2. * mv_xsi[5][1] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][2]) },
	  { -2. * mv_xsi[5][1] * mv_xsi[5][2] * (1. - mv_xsi[5][2]), -2. * mv_xsi[5][0] * mv_xsi[5][2] * (1. - mv_xsi[5][2]), -2. * mv_xsi[5][0] * mv_xsi[5][1] * (1. - 2. * mv_xsi[5][2]) },
	  { 0., (0.5 - 2. * mv_xsi[5][1]) * (1. - mv_xsi[5][2]) * mv_xsi[5][2], -0.5 * mv_xsi[5][1] * (2. * mv_xsi[5][1] - 1.) * (1. - 2. * mv_xsi[5][2]) },
	  { (3. - 4. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), (3. - 4. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (2. * mv_xsi[5][2]) },
	  { (-4. + 8. * mv_xsi[5][0] + 4. * mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), 4. * mv_xsi[5][0] * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), 4. * mv_xsi[5][0] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (2. * mv_xsi[5][2]) },
	  { (1. - 4. * mv_xsi[5][0]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), 0., mv_xsi[5][0] * (2. * mv_xsi[5][0] - 1.) * (-2. * mv_xsi[5][2]) },
	  { 4. * mv_xsi[5][1] * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), (-4. + 4. * mv_xsi[5][0] + 8. * mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), 4. * mv_xsi[5][1] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (2. * mv_xsi[5][2]) },
	  { -4. * mv_xsi[5][1] * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), -4. * mv_xsi[5][0] * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), -4. * mv_xsi[5][0] * mv_xsi[5][1] * (2. * mv_xsi[5][2]) },
	  { 0., (1. - 4. * mv_xsi[5][1]) * (-1. + mv_xsi[5][2] * mv_xsi[5][2]), mv_xsi[5][1] * (2. * mv_xsi[5][1] - 1.) * (-2. * mv_xsi[5][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2], -0.5 * (3. - 4. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2], -0.5 * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. - 2. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (1. + 2. * mv_xsi[5][2]) },
	  { (2. - 4. * mv_xsi[5][0] - 2. * mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2], -2. * mv_xsi[5][0] * (1. + mv_xsi[5][2]) * mv_xsi[5][2], -2. * mv_xsi[5][0] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. + 2. * mv_xsi[5][2]) },
	  { (-0.5 + 2. * mv_xsi[5][0]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2], 0., 0.5 * mv_xsi[5][0] * (2. * mv_xsi[5][0] - 1.) * (1. + 2. * mv_xsi[5][2]) },
	  { -2. * mv_xsi[5][1] * (1. + mv_xsi[5][2]) * mv_xsi[5][2], (2. - 2. * mv_xsi[5][0] - 4. * mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2], -2. * mv_xsi[5][1] * (-1. + mv_xsi[5][0] + mv_xsi[5][1]) * (1. + 2. * mv_xsi[5][2]) },
	  { 2. * mv_xsi[5][1] * mv_xsi[5][2] * (1. + mv_xsi[5][2]), 2. * mv_xsi[5][0] * mv_xsi[5][2] * (1. + mv_xsi[5][2]), 2. * mv_xsi[5][0] * mv_xsi[5][1] * (1. + 2. * mv_xsi[5][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[5][1]) * (1. + mv_xsi[5][2]) * mv_xsi[5][2], 0.5 * mv_xsi[5][1] * (2. * mv_xsi[5][1] - 1.) * (1. + 2. * mv_xsi[5][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2], 0.5 * (3. - 4. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2], 0.5 * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][2]) },
	  { (-2. + 4. * mv_xsi[6][0] + 2. * mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2], 2. * mv_xsi[6][0] * (1. - mv_xsi[6][2]) * mv_xsi[6][2], 2. * mv_xsi[6][0] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][2]) },
	  { (0.5 - 2. * mv_xsi[6][0]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2], 0., -0.5 * mv_xsi[6][0] * (2. * mv_xsi[6][0] - 1.) * (1. - 2. * mv_xsi[6][2]) },
	  { 2. * mv_xsi[6][1] * (1. - mv_xsi[6][2]) * mv_xsi[6][2], (-2. + 2. * mv_xsi[6][0] + 4. * mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2], 2. * mv_xsi[6][1] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][2]) },
	  { -2. * mv_xsi[6][1] * mv_xsi[6][2] * (1. - mv_xsi[6][2]), -2. * mv_xsi[6][0] * mv_xsi[6][2] * (1. - mv_xsi[6][2]), -2. * mv_xsi[6][0] * mv_xsi[6][1] * (1. - 2. * mv_xsi[6][2]) },
	  { 0., (0.5 - 2. * mv_xsi[6][1]) * (1. - mv_xsi[6][2]) * mv_xsi[6][2], -0.5 * mv_xsi[6][1] * (2. * mv_xsi[6][1] - 1.) * (1. - 2. * mv_xsi[6][2]) },
	  { (3. - 4. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), (3. - 4. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (2. * mv_xsi[6][2]) },
	  { (-4. + 8. * mv_xsi[6][0] + 4. * mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), 4. * mv_xsi[6][0] * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), 4. * mv_xsi[6][0] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (2. * mv_xsi[6][2]) },
	  { (1. - 4. * mv_xsi[6][0]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), 0., mv_xsi[6][0] * (2. * mv_xsi[6][0] - 1.) * (-2. * mv_xsi[6][2]) },
	  { 4. * mv_xsi[6][1] * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), (-4. + 4. * mv_xsi[6][0] + 8. * mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), 4. * mv_xsi[6][1] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (2. * mv_xsi[6][2]) },
	  { -4. * mv_xsi[6][1] * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), -4. * mv_xsi[6][0] * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), -4. * mv_xsi[6][0] * mv_xsi[6][1] * (2. * mv_xsi[6][2]) },
	  { 0., (1. - 4. * mv_xsi[6][1]) * (-1. + mv_xsi[6][2] * mv_xsi[6][2]), mv_xsi[6][1] * (2. * mv_xsi[6][1] - 1.) * (-2. * mv_xsi[6][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2], -0.5 * (3. - 4. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2], -0.5 * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. - 2. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (1. + 2. * mv_xsi[6][2]) },
	  { (2. - 4. * mv_xsi[6][0] - 2. * mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2], -2. * mv_xsi[6][0] * (1. + mv_xsi[6][2]) * mv_xsi[6][2], -2. * mv_xsi[6][0] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. + 2. * mv_xsi[6][2]) },
	  { (-0.5 + 2. * mv_xsi[6][0]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2], 0., 0.5 * mv_xsi[6][0] * (2. * mv_xsi[6][0] - 1.) * (1. + 2. * mv_xsi[6][2]) },
	  { -2. * mv_xsi[6][1] * (1. + mv_xsi[6][2]) * mv_xsi[6][2], (2. - 2. * mv_xsi[6][0] - 4. * mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2], -2. * mv_xsi[6][1] * (-1. + mv_xsi[6][0] + mv_xsi[6][1]) * (1. + 2. * mv_xsi[6][2]) },
	  { 2. * mv_xsi[6][1] * mv_xsi[6][2] * (1. + mv_xsi[6][2]), 2. * mv_xsi[6][0] * mv_xsi[6][2] * (1. + mv_xsi[6][2]), 2. * mv_xsi[6][0] * mv_xsi[6][1] * (1. + 2. * mv_xsi[6][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[6][1]) * (1. + mv_xsi[6][2]) * mv_xsi[6][2], 0.5 * mv_xsi[6][1] * (2. * mv_xsi[6][1] - 1.) * (1. + 2. * mv_xsi[6][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2], 0.5 * (3. - 4. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2], 0.5 * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][2]) },
	  { (-2. + 4. * mv_xsi[7][0] + 2. * mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2], 2. * mv_xsi[7][0] * (1. - mv_xsi[7][2]) * mv_xsi[7][2], 2. * mv_xsi[7][0] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][2]) },
	  { (0.5 - 2. * mv_xsi[7][0]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2], 0., -0.5 * mv_xsi[7][0] * (2. * mv_xsi[7][0] - 1.) * (1. - 2. * mv_xsi[7][2]) },
	  { 2. * mv_xsi[7][1] * (1. - mv_xsi[7][2]) * mv_xsi[7][2], (-2. + 2. * mv_xsi[7][0] + 4. * mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2], 2. * mv_xsi[7][1] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][2]) },
	  { -2. * mv_xsi[7][1] * mv_xsi[7][2] * (1. - mv_xsi[7][2]), -2. * mv_xsi[7][0] * mv_xsi[7][2] * (1. - mv_xsi[7][2]), -2. * mv_xsi[7][0] * mv_xsi[7][1] * (1. - 2. * mv_xsi[7][2]) },
	  { 0., (0.5 - 2. * mv_xsi[7][1]) * (1. - mv_xsi[7][2]) * mv_xsi[7][2], -0.5 * mv_xsi[7][1] * (2. * mv_xsi[7][1] - 1.) * (1. - 2. * mv_xsi[7][2]) },
	  { (3. - 4. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), (3. - 4. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (2. * mv_xsi[7][2]) },
	  { (-4. + 8. * mv_xsi[7][0] + 4. * mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), 4. * mv_xsi[7][0] * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), 4. * mv_xsi[7][0] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (2. * mv_xsi[7][2]) },
	  { (1. - 4. * mv_xsi[7][0]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), 0., mv_xsi[7][0] * (2. * mv_xsi[7][0] - 1.) * (-2. * mv_xsi[7][2]) },
	  { 4. * mv_xsi[7][1] * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), (-4. + 4. * mv_xsi[7][0] + 8. * mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), 4. * mv_xsi[7][1] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (2. * mv_xsi[7][2]) },
	  { -4. * mv_xsi[7][1] * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), -4. * mv_xsi[7][0] * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), -4. * mv_xsi[7][0] * mv_xsi[7][1] * (2. * mv_xsi[7][2]) },
	  { 0., (1. - 4. * mv_xsi[7][1]) * (-1. + mv_xsi[7][2] * mv_xsi[7][2]), mv_xsi[7][1] * (2. * mv_xsi[7][1] - 1.) * (-2. * mv_xsi[7][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2], -0.5 * (3. - 4. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2], -0.5 * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. - 2. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (1. + 2. * mv_xsi[7][2]) },
	  { (2. - 4. * mv_xsi[7][0] - 2. * mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2], -2. * mv_xsi[7][0] * (1. + mv_xsi[7][2]) * mv_xsi[7][2], -2. * mv_xsi[7][0] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. + 2. * mv_xsi[7][2]) },
	  { (-0.5 + 2. * mv_xsi[7][0]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2], 0., 0.5 * mv_xsi[7][0] * (2. * mv_xsi[7][0] - 1.) * (1. + 2. * mv_xsi[7][2]) },
	  { -2. * mv_xsi[7][1] * (1. + mv_xsi[7][2]) * mv_xsi[7][2], (2. - 2. * mv_xsi[7][0] - 4. * mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2], -2. * mv_xsi[7][1] * (-1. + mv_xsi[7][0] + mv_xsi[7][1]) * (1. + 2. * mv_xsi[7][2]) },
	  { 2. * mv_xsi[7][1] * mv_xsi[7][2] * (1. + mv_xsi[7][2]), 2. * mv_xsi[7][0] * mv_xsi[7][2] * (1. + mv_xsi[7][2]), 2. * mv_xsi[7][0] * mv_xsi[7][1] * (1. + 2. * mv_xsi[7][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[7][1]) * (1. + mv_xsi[7][2]) * mv_xsi[7][2], 0.5 * mv_xsi[7][1] * (2. * mv_xsi[7][1] - 1.) * (1. + 2. * mv_xsi[7][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2], 0.5 * (3. - 4. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2], 0.5 * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][2]) },
	  { (-2. + 4. * mv_xsi[8][0] + 2. * mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2], 2. * mv_xsi[8][0] * (1. - mv_xsi[8][2]) * mv_xsi[8][2], 2. * mv_xsi[8][0] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][2]) },
	  { (0.5 - 2. * mv_xsi[8][0]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2], 0., -0.5 * mv_xsi[8][0] * (2. * mv_xsi[8][0] - 1.) * (1. - 2. * mv_xsi[8][2]) },
	  { 2. * mv_xsi[8][1] * (1. - mv_xsi[8][2]) * mv_xsi[8][2], (-2. + 2. * mv_xsi[8][0] + 4. * mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2], 2. * mv_xsi[8][1] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][2]) },
	  { -2. * mv_xsi[8][1] * mv_xsi[8][2] * (1. - mv_xsi[8][2]), -2. * mv_xsi[8][0] * mv_xsi[8][2] * (1. - mv_xsi[8][2]), -2. * mv_xsi[8][0] * mv_xsi[8][1] * (1. - 2. * mv_xsi[8][2]) },
	  { 0., (0.5 - 2. * mv_xsi[8][1]) * (1. - mv_xsi[8][2]) * mv_xsi[8][2], -0.5 * mv_xsi[8][1] * (2. * mv_xsi[8][1] - 1.) * (1. - 2. * mv_xsi[8][2]) },
	  { (3. - 4. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), (3. - 4. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (2. * mv_xsi[8][2]) },
	  { (-4. + 8. * mv_xsi[8][0] + 4. * mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), 4. * mv_xsi[8][0] * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), 4. * mv_xsi[8][0] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (2. * mv_xsi[8][2]) },
	  { (1. - 4. * mv_xsi[8][0]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), 0., mv_xsi[8][0] * (2. * mv_xsi[8][0] - 1.) * (-2. * mv_xsi[8][2]) },
	  { 4. * mv_xsi[8][1] * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), (-4. + 4. * mv_xsi[8][0] + 8. * mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), 4. * mv_xsi[8][1] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (2. * mv_xsi[8][2]) },
	  { -4. * mv_xsi[8][1] * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), -4. * mv_xsi[8][0] * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), -4. * mv_xsi[8][0] * mv_xsi[8][1] * (2. * mv_xsi[8][2]) },
	  { 0., (1. - 4. * mv_xsi[8][1]) * (-1. + mv_xsi[8][2] * mv_xsi[8][2]), mv_xsi[8][1] * (2. * mv_xsi[8][1] - 1.) * (-2. * mv_xsi[8][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2], -0.5 * (3. - 4. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2], -0.5 * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. - 2. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (1. + 2. * mv_xsi[8][2]) },
	  { (2. - 4. * mv_xsi[8][0] - 2. * mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2], -2. * mv_xsi[8][0] * (1. + mv_xsi[8][2]) * mv_xsi[8][2], -2. * mv_xsi[8][0] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. + 2. * mv_xsi[8][2]) },
	  { (-0.5 + 2. * mv_xsi[8][0]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2], 0., 0.5 * mv_xsi[8][0] * (2. * mv_xsi[8][0] - 1.) * (1. + 2. * mv_xsi[8][2]) },
	  { -2. * mv_xsi[8][1] * (1. + mv_xsi[8][2]) * mv_xsi[8][2], (2. - 2. * mv_xsi[8][0] - 4. * mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2], -2. * mv_xsi[8][1] * (-1. + mv_xsi[8][0] + mv_xsi[8][1]) * (1. + 2. * mv_xsi[8][2]) },
	  { 2. * mv_xsi[8][1] * mv_xsi[8][2] * (1. + mv_xsi[8][2]), 2. * mv_xsi[8][0] * mv_xsi[8][2] * (1. + mv_xsi[8][2]), 2. * mv_xsi[8][0] * mv_xsi[8][1] * (1. + 2. * mv_xsi[8][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[8][1]) * (1. + mv_xsi[8][2]) * mv_xsi[8][2], 0.5 * mv_xsi[8][1] * (2. * mv_xsi[8][1] - 1.) * (1. + 2. * mv_xsi[8][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2], 0.5 * (3. - 4. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2], 0.5 * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][2]) },
	  { (-2. + 4. * mv_xsi[9][0] + 2. * mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2], 2. * mv_xsi[9][0] * (1. - mv_xsi[9][2]) * mv_xsi[9][2], 2. * mv_xsi[9][0] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][2]) },
	  { (0.5 - 2. * mv_xsi[9][0]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2], 0., -0.5 * mv_xsi[9][0] * (2. * mv_xsi[9][0] - 1.) * (1. - 2. * mv_xsi[9][2]) },
	  { 2. * mv_xsi[9][1] * (1. - mv_xsi[9][2]) * mv_xsi[9][2], (-2. + 2. * mv_xsi[9][0] + 4. * mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2], 2. * mv_xsi[9][1] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][2]) },
	  { -2. * mv_xsi[9][1] * mv_xsi[9][2] * (1. - mv_xsi[9][2]), -2. * mv_xsi[9][0] * mv_xsi[9][2] * (1. - mv_xsi[9][2]), -2. * mv_xsi[9][0] * mv_xsi[9][1] * (1. - 2. * mv_xsi[9][2]) },
	  { 0., (0.5 - 2. * mv_xsi[9][1]) * (1. - mv_xsi[9][2]) * mv_xsi[9][2], -0.5 * mv_xsi[9][1] * (2. * mv_xsi[9][1] - 1.) * (1. - 2. * mv_xsi[9][2]) },
	  { (3. - 4. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), (3. - 4. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (2. * mv_xsi[9][2]) },
	  { (-4. + 8. * mv_xsi[9][0] + 4. * mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), 4. * mv_xsi[9][0] * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), 4. * mv_xsi[9][0] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (2. * mv_xsi[9][2]) },
	  { (1. - 4. * mv_xsi[9][0]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), 0., mv_xsi[9][0] * (2. * mv_xsi[9][0] - 1.) * (-2. * mv_xsi[9][2]) },
	  { 4. * mv_xsi[9][1] * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), (-4. + 4. * mv_xsi[9][0] + 8. * mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), 4. * mv_xsi[9][1] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (2. * mv_xsi[9][2]) },
	  { -4. * mv_xsi[9][1] * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), -4. * mv_xsi[9][0] * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), -4. * mv_xsi[9][0] * mv_xsi[9][1] * (2. * mv_xsi[9][2]) },
	  { 0., (1. - 4. * mv_xsi[9][1]) * (-1. + mv_xsi[9][2] * mv_xsi[9][2]), mv_xsi[9][1] * (2. * mv_xsi[9][1] - 1.) * (-2. * mv_xsi[9][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2], -0.5 * (3. - 4. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2], -0.5 * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. - 2. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (1. + 2. * mv_xsi[9][2]) },
	  { (2. - 4. * mv_xsi[9][0] - 2. * mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2], -2. * mv_xsi[9][0] * (1. + mv_xsi[9][2]) * mv_xsi[9][2], -2. * mv_xsi[9][0] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. + 2. * mv_xsi[9][2]) },
	  { (-0.5 + 2. * mv_xsi[9][0]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2], 0., 0.5 * mv_xsi[9][0] * (2. * mv_xsi[9][0] - 1.) * (1. + 2. * mv_xsi[9][2]) },
	  { -2. * mv_xsi[9][1] * (1. + mv_xsi[9][2]) * mv_xsi[9][2], (2. - 2. * mv_xsi[9][0] - 4. * mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2], -2. * mv_xsi[9][1] * (-1. + mv_xsi[9][0] + mv_xsi[9][1]) * (1. + 2. * mv_xsi[9][2]) },
	  { 2. * mv_xsi[9][1] * mv_xsi[9][2] * (1. + mv_xsi[9][2]), 2. * mv_xsi[9][0] * mv_xsi[9][2] * (1. + mv_xsi[9][2]), 2. * mv_xsi[9][0] * mv_xsi[9][1] * (1. + 2. * mv_xsi[9][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[9][1]) * (1. + mv_xsi[9][2]) * mv_xsi[9][2], 0.5 * mv_xsi[9][1] * (2. * mv_xsi[9][1] - 1.) * (1. + 2. * mv_xsi[9][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2], 0.5 * (3. - 4. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2], 0.5 * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][2]) },
	  { (-2. + 4. * mv_xsi[10][0] + 2. * mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2], 2. * mv_xsi[10][0] * (1. - mv_xsi[10][2]) * mv_xsi[10][2], 2. * mv_xsi[10][0] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][2]) },
	  { (0.5 - 2. * mv_xsi[10][0]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2], 0., -0.5 * mv_xsi[10][0] * (2. * mv_xsi[10][0] - 1.) * (1. - 2. * mv_xsi[10][2]) },
	  { 2. * mv_xsi[10][1] * (1. - mv_xsi[10][2]) * mv_xsi[10][2], (-2. + 2. * mv_xsi[10][0] + 4. * mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2], 2. * mv_xsi[10][1] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][2]) },
	  { -2. * mv_xsi[10][1] * mv_xsi[10][2] * (1. - mv_xsi[10][2]), -2. * mv_xsi[10][0] * mv_xsi[10][2] * (1. - mv_xsi[10][2]), -2. * mv_xsi[10][0] * mv_xsi[10][1] * (1. - 2. * mv_xsi[10][2]) },
	  { 0., (0.5 - 2. * mv_xsi[10][1]) * (1. - mv_xsi[10][2]) * mv_xsi[10][2], -0.5 * mv_xsi[10][1] * (2. * mv_xsi[10][1] - 1.) * (1. - 2. * mv_xsi[10][2]) },
	  { (3. - 4. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), (3. - 4. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (2. * mv_xsi[10][2]) },
	  { (-4. + 8. * mv_xsi[10][0] + 4. * mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), 4. * mv_xsi[10][0] * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), 4. * mv_xsi[10][0] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (2. * mv_xsi[10][2]) },
	  { (1. - 4. * mv_xsi[10][0]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), 0., mv_xsi[10][0] * (2. * mv_xsi[10][0] - 1.) * (-2. * mv_xsi[10][2]) },
	  { 4. * mv_xsi[10][1] * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), (-4. + 4. * mv_xsi[10][0] + 8. * mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), 4. * mv_xsi[10][1] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (2. * mv_xsi[10][2]) },
	  { -4. * mv_xsi[10][1] * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), -4. * mv_xsi[10][0] * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), -4. * mv_xsi[10][0] * mv_xsi[10][1] * (2. * mv_xsi[10][2]) },
	  { 0., (1. - 4. * mv_xsi[10][1]) * (-1. + mv_xsi[10][2] * mv_xsi[10][2]), mv_xsi[10][1] * (2. * mv_xsi[10][1] - 1.) * (-2. * mv_xsi[10][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2], -0.5 * (3. - 4. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2], -0.5 * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. - 2. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (1. + 2. * mv_xsi[10][2]) },
	  { (2. - 4. * mv_xsi[10][0] - 2. * mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2], -2. * mv_xsi[10][0] * (1. + mv_xsi[10][2]) * mv_xsi[10][2], -2. * mv_xsi[10][0] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. + 2. * mv_xsi[10][2]) },
	  { (-0.5 + 2. * mv_xsi[10][0]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2], 0., 0.5 * mv_xsi[10][0] * (2. * mv_xsi[10][0] - 1.) * (1. + 2. * mv_xsi[10][2]) },
	  { -2. * mv_xsi[10][1] * (1. + mv_xsi[10][2]) * mv_xsi[10][2], (2. - 2. * mv_xsi[10][0] - 4. * mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2], -2. * mv_xsi[10][1] * (-1. + mv_xsi[10][0] + mv_xsi[10][1]) * (1. + 2. * mv_xsi[10][2]) },
	  { 2. * mv_xsi[10][1] * mv_xsi[10][2] * (1. + mv_xsi[10][2]), 2. * mv_xsi[10][0] * mv_xsi[10][2] * (1. + mv_xsi[10][2]), 2. * mv_xsi[10][0] * mv_xsi[10][1] * (1. + 2. * mv_xsi[10][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[10][1]) * (1. + mv_xsi[10][2]) * mv_xsi[10][2], 0.5 * mv_xsi[10][1] * (2. * mv_xsi[10][1] - 1.) * (1. + 2. * mv_xsi[10][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2], 0.5 * (3. - 4. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2], 0.5 * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][2]) },
	  { (-2. + 4. * mv_xsi[11][0] + 2. * mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2], 2. * mv_xsi[11][0] * (1. - mv_xsi[11][2]) * mv_xsi[11][2], 2. * mv_xsi[11][0] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][2]) },
	  { (0.5 - 2. * mv_xsi[11][0]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2], 0., -0.5 * mv_xsi[11][0] * (2. * mv_xsi[11][0] - 1.) * (1. - 2. * mv_xsi[11][2]) },
	  { 2. * mv_xsi[11][1] * (1. - mv_xsi[11][2]) * mv_xsi[11][2], (-2. + 2. * mv_xsi[11][0] + 4. * mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2], 2. * mv_xsi[11][1] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][2]) },
	  { -2. * mv_xsi[11][1] * mv_xsi[11][2] * (1. - mv_xsi[11][2]), -2. * mv_xsi[11][0] * mv_xsi[11][2] * (1. - mv_xsi[11][2]), -2. * mv_xsi[11][0] * mv_xsi[11][1] * (1. - 2. * mv_xsi[11][2]) },
	  { 0., (0.5 - 2. * mv_xsi[11][1]) * (1. - mv_xsi[11][2]) * mv_xsi[11][2], -0.5 * mv_xsi[11][1] * (2. * mv_xsi[11][1] - 1.) * (1. - 2. * mv_xsi[11][2]) },
	  { (3. - 4. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), (3. - 4. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (2. * mv_xsi[11][2]) },
	  { (-4. + 8. * mv_xsi[11][0] + 4. * mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), 4. * mv_xsi[11][0] * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), 4. * mv_xsi[11][0] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (2. * mv_xsi[11][2]) },
	  { (1. - 4. * mv_xsi[11][0]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), 0., mv_xsi[11][0] * (2. * mv_xsi[11][0] - 1.) * (-2. * mv_xsi[11][2]) },
	  { 4. * mv_xsi[11][1] * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), (-4. + 4. * mv_xsi[11][0] + 8. * mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), 4. * mv_xsi[11][1] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (2. * mv_xsi[11][2]) },
	  { -4. * mv_xsi[11][1] * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), -4. * mv_xsi[11][0] * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), -4. * mv_xsi[11][0] * mv_xsi[11][1] * (2. * mv_xsi[11][2]) },
	  { 0., (1. - 4. * mv_xsi[11][1]) * (-1. + mv_xsi[11][2] * mv_xsi[11][2]), mv_xsi[11][1] * (2. * mv_xsi[11][1] - 1.) * (-2. * mv_xsi[11][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2], -0.5 * (3. - 4. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2], -0.5 * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. - 2. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (1. + 2. * mv_xsi[11][2]) },
	  { (2. - 4. * mv_xsi[11][0] - 2. * mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2], -2. * mv_xsi[11][0] * (1. + mv_xsi[11][2]) * mv_xsi[11][2], -2. * mv_xsi[11][0] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. + 2. * mv_xsi[11][2]) },
	  { (-0.5 + 2. * mv_xsi[11][0]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2], 0., 0.5 * mv_xsi[11][0] * (2. * mv_xsi[11][0] - 1.) * (1. + 2. * mv_xsi[11][2]) },
	  { -2. * mv_xsi[11][1] * (1. + mv_xsi[11][2]) * mv_xsi[11][2], (2. - 2. * mv_xsi[11][0] - 4. * mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2], -2. * mv_xsi[11][1] * (-1. + mv_xsi[11][0] + mv_xsi[11][1]) * (1. + 2. * mv_xsi[11][2]) },
	  { 2. * mv_xsi[11][1] * mv_xsi[11][2] * (1. + mv_xsi[11][2]), 2. * mv_xsi[11][0] * mv_xsi[11][2] * (1. + mv_xsi[11][2]), 2. * mv_xsi[11][0] * mv_xsi[11][1] * (1. + 2. * mv_xsi[11][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[11][1]) * (1. + mv_xsi[11][2]) * mv_xsi[11][2], 0.5 * mv_xsi[11][1] * (2. * mv_xsi[11][1] - 1.) * (1. + 2. * mv_xsi[11][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2], 0.5 * (3. - 4. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2], 0.5 * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][2]) },
	  { (-2. + 4. * mv_xsi[12][0] + 2. * mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2], 2. * mv_xsi[12][0] * (1. - mv_xsi[12][2]) * mv_xsi[12][2], 2. * mv_xsi[12][0] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][2]) },
	  { (0.5 - 2. * mv_xsi[12][0]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2], 0., -0.5 * mv_xsi[12][0] * (2. * mv_xsi[12][0] - 1.) * (1. - 2. * mv_xsi[12][2]) },
	  { 2. * mv_xsi[12][1] * (1. - mv_xsi[12][2]) * mv_xsi[12][2], (-2. + 2. * mv_xsi[12][0] + 4. * mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2], 2. * mv_xsi[12][1] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][2]) },
	  { -2. * mv_xsi[12][1] * mv_xsi[12][2] * (1. - mv_xsi[12][2]), -2. * mv_xsi[12][0] * mv_xsi[12][2] * (1. - mv_xsi[12][2]), -2. * mv_xsi[12][0] * mv_xsi[12][1] * (1. - 2. * mv_xsi[12][2]) },
	  { 0., (0.5 - 2. * mv_xsi[12][1]) * (1. - mv_xsi[12][2]) * mv_xsi[12][2], -0.5 * mv_xsi[12][1] * (2. * mv_xsi[12][1] - 1.) * (1. - 2. * mv_xsi[12][2]) },
	  { (3. - 4. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), (3. - 4. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (2. * mv_xsi[12][2]) },
	  { (-4. + 8. * mv_xsi[12][0] + 4. * mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), 4. * mv_xsi[12][0] * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), 4. * mv_xsi[12][0] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (2. * mv_xsi[12][2]) },
	  { (1. - 4. * mv_xsi[12][0]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), 0., mv_xsi[12][0] * (2. * mv_xsi[12][0] - 1.) * (-2. * mv_xsi[12][2]) },
	  { 4. * mv_xsi[12][1] * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), (-4. + 4. * mv_xsi[12][0] + 8. * mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), 4. * mv_xsi[12][1] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (2. * mv_xsi[12][2]) },
	  { -4. * mv_xsi[12][1] * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), -4. * mv_xsi[12][0] * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), -4. * mv_xsi[12][0] * mv_xsi[12][1] * (2. * mv_xsi[12][2]) },
	  { 0., (1. - 4. * mv_xsi[12][1]) * (-1. + mv_xsi[12][2] * mv_xsi[12][2]), mv_xsi[12][1] * (2. * mv_xsi[12][1] - 1.) * (-2. * mv_xsi[12][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2], -0.5 * (3. - 4. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2], -0.5 * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. - 2. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (1. + 2. * mv_xsi[12][2]) },
	  { (2. - 4. * mv_xsi[12][0] - 2. * mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2], -2. * mv_xsi[12][0] * (1. + mv_xsi[12][2]) * mv_xsi[12][2], -2. * mv_xsi[12][0] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. + 2. * mv_xsi[12][2]) },
	  { (-0.5 + 2. * mv_xsi[12][0]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2], 0., 0.5 * mv_xsi[12][0] * (2. * mv_xsi[12][0] - 1.) * (1. + 2. * mv_xsi[12][2]) },
	  { -2. * mv_xsi[12][1] * (1. + mv_xsi[12][2]) * mv_xsi[12][2], (2. - 2. * mv_xsi[12][0] - 4. * mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2], -2. * mv_xsi[12][1] * (-1. + mv_xsi[12][0] + mv_xsi[12][1]) * (1. + 2. * mv_xsi[12][2]) },
	  { 2. * mv_xsi[12][1] * mv_xsi[12][2] * (1. + mv_xsi[12][2]), 2. * mv_xsi[12][0] * mv_xsi[12][2] * (1. + mv_xsi[12][2]), 2. * mv_xsi[12][0] * mv_xsi[12][1] * (1. + 2. * mv_xsi[12][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[12][1]) * (1. + mv_xsi[12][2]) * mv_xsi[12][2], 0.5 * mv_xsi[12][1] * (2. * mv_xsi[12][1] - 1.) * (1. + 2. * mv_xsi[12][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2], 0.5 * (3. - 4. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2], 0.5 * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][2]) },
	  { (-2. + 4. * mv_xsi[13][0] + 2. * mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2], 2. * mv_xsi[13][0] * (1. - mv_xsi[13][2]) * mv_xsi[13][2], 2. * mv_xsi[13][0] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][2]) },
	  { (0.5 - 2. * mv_xsi[13][0]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2], 0., -0.5 * mv_xsi[13][0] * (2. * mv_xsi[13][0] - 1.) * (1. - 2. * mv_xsi[13][2]) },
	  { 2. * mv_xsi[13][1] * (1. - mv_xsi[13][2]) * mv_xsi[13][2], (-2. + 2. * mv_xsi[13][0] + 4. * mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2], 2. * mv_xsi[13][1] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][2]) },
	  { -2. * mv_xsi[13][1] * mv_xsi[13][2] * (1. - mv_xsi[13][2]), -2. * mv_xsi[13][0] * mv_xsi[13][2] * (1. - mv_xsi[13][2]), -2. * mv_xsi[13][0] * mv_xsi[13][1] * (1. - 2. * mv_xsi[13][2]) },
	  { 0., (0.5 - 2. * mv_xsi[13][1]) * (1. - mv_xsi[13][2]) * mv_xsi[13][2], -0.5 * mv_xsi[13][1] * (2. * mv_xsi[13][1] - 1.) * (1. - 2. * mv_xsi[13][2]) },
	  { (3. - 4. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), (3. - 4. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (2. * mv_xsi[13][2]) },
	  { (-4. + 8. * mv_xsi[13][0] + 4. * mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), 4. * mv_xsi[13][0] * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), 4. * mv_xsi[13][0] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (2. * mv_xsi[13][2]) },
	  { (1. - 4. * mv_xsi[13][0]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), 0., mv_xsi[13][0] * (2. * mv_xsi[13][0] - 1.) * (-2. * mv_xsi[13][2]) },
	  { 4. * mv_xsi[13][1] * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), (-4. + 4. * mv_xsi[13][0] + 8. * mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), 4. * mv_xsi[13][1] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (2. * mv_xsi[13][2]) },
	  { -4. * mv_xsi[13][1] * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), -4. * mv_xsi[13][0] * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), -4. * mv_xsi[13][0] * mv_xsi[13][1] * (2. * mv_xsi[13][2]) },
	  { 0., (1. - 4. * mv_xsi[13][1]) * (-1. + mv_xsi[13][2] * mv_xsi[13][2]), mv_xsi[13][1] * (2. * mv_xsi[13][1] - 1.) * (-2. * mv_xsi[13][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2], -0.5 * (3. - 4. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2], -0.5 * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. - 2. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (1. + 2. * mv_xsi[13][2]) },
	  { (2. - 4. * mv_xsi[13][0] - 2. * mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2], -2. * mv_xsi[13][0] * (1. + mv_xsi[13][2]) * mv_xsi[13][2], -2. * mv_xsi[13][0] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. + 2. * mv_xsi[13][2]) },
	  { (-0.5 + 2. * mv_xsi[13][0]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2], 0., 0.5 * mv_xsi[13][0] * (2. * mv_xsi[13][0] - 1.) * (1. + 2. * mv_xsi[13][2]) },
	  { -2. * mv_xsi[13][1] * (1. + mv_xsi[13][2]) * mv_xsi[13][2], (2. - 2. * mv_xsi[13][0] - 4. * mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2], -2. * mv_xsi[13][1] * (-1. + mv_xsi[13][0] + mv_xsi[13][1]) * (1. + 2. * mv_xsi[13][2]) },
	  { 2. * mv_xsi[13][1] * mv_xsi[13][2] * (1. + mv_xsi[13][2]), 2. * mv_xsi[13][0] * mv_xsi[13][2] * (1. + mv_xsi[13][2]), 2. * mv_xsi[13][0] * mv_xsi[13][1] * (1. + 2. * mv_xsi[13][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[13][1]) * (1. + mv_xsi[13][2]) * mv_xsi[13][2], 0.5 * mv_xsi[13][1] * (2. * mv_xsi[13][1] - 1.) * (1. + 2. * mv_xsi[13][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2], 0.5 * (3. - 4. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2], 0.5 * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][2]) },
	  { (-2. + 4. * mv_xsi[14][0] + 2. * mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2], 2. * mv_xsi[14][0] * (1. - mv_xsi[14][2]) * mv_xsi[14][2], 2. * mv_xsi[14][0] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][2]) },
	  { (0.5 - 2. * mv_xsi[14][0]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2], 0., -0.5 * mv_xsi[14][0] * (2. * mv_xsi[14][0] - 1.) * (1. - 2. * mv_xsi[14][2]) },
	  { 2. * mv_xsi[14][1] * (1. - mv_xsi[14][2]) * mv_xsi[14][2], (-2. + 2. * mv_xsi[14][0] + 4. * mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2], 2. * mv_xsi[14][1] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][2]) },
	  { -2. * mv_xsi[14][1] * mv_xsi[14][2] * (1. - mv_xsi[14][2]), -2. * mv_xsi[14][0] * mv_xsi[14][2] * (1. - mv_xsi[14][2]), -2. * mv_xsi[14][0] * mv_xsi[14][1] * (1. - 2. * mv_xsi[14][2]) },
	  { 0., (0.5 - 2. * mv_xsi[14][1]) * (1. - mv_xsi[14][2]) * mv_xsi[14][2], -0.5 * mv_xsi[14][1] * (2. * mv_xsi[14][1] - 1.) * (1. - 2. * mv_xsi[14][2]) },
	  { (3. - 4. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), (3. - 4. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (2. * mv_xsi[14][2]) },
	  { (-4. + 8. * mv_xsi[14][0] + 4. * mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), 4. * mv_xsi[14][0] * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), 4. * mv_xsi[14][0] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (2. * mv_xsi[14][2]) },
	  { (1. - 4. * mv_xsi[14][0]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), 0., mv_xsi[14][0] * (2. * mv_xsi[14][0] - 1.) * (-2. * mv_xsi[14][2]) },
	  { 4. * mv_xsi[14][1] * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), (-4. + 4. * mv_xsi[14][0] + 8. * mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), 4. * mv_xsi[14][1] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (2. * mv_xsi[14][2]) },
	  { -4. * mv_xsi[14][1] * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), -4. * mv_xsi[14][0] * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), -4. * mv_xsi[14][0] * mv_xsi[14][1] * (2. * mv_xsi[14][2]) },
	  { 0., (1. - 4. * mv_xsi[14][1]) * (-1. + mv_xsi[14][2] * mv_xsi[14][2]), mv_xsi[14][1] * (2. * mv_xsi[14][1] - 1.) * (-2. * mv_xsi[14][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2], -0.5 * (3. - 4. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2], -0.5 * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. - 2. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (1. + 2. * mv_xsi[14][2]) },
	  { (2. - 4. * mv_xsi[14][0] - 2. * mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2], -2. * mv_xsi[14][0] * (1. + mv_xsi[14][2]) * mv_xsi[14][2], -2. * mv_xsi[14][0] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. + 2. * mv_xsi[14][2]) },
	  { (-0.5 + 2. * mv_xsi[14][0]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2], 0., 0.5 * mv_xsi[14][0] * (2. * mv_xsi[14][0] - 1.) * (1. + 2. * mv_xsi[14][2]) },
	  { -2. * mv_xsi[14][1] * (1. + mv_xsi[14][2]) * mv_xsi[14][2], (2. - 2. * mv_xsi[14][0] - 4. * mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2], -2. * mv_xsi[14][1] * (-1. + mv_xsi[14][0] + mv_xsi[14][1]) * (1. + 2. * mv_xsi[14][2]) },
	  { 2. * mv_xsi[14][1] * mv_xsi[14][2] * (1. + mv_xsi[14][2]), 2. * mv_xsi[14][0] * mv_xsi[14][2] * (1. + mv_xsi[14][2]), 2. * mv_xsi[14][0] * mv_xsi[14][1] * (1. + 2. * mv_xsi[14][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[14][1]) * (1. + mv_xsi[14][2]) * mv_xsi[14][2], 0.5 * mv_xsi[14][1] * (2. * mv_xsi[14][1] - 1.) * (1. + 2. * mv_xsi[14][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2], 0.5 * (3. - 4. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2], 0.5 * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][2]) },
	  { (-2. + 4. * mv_xsi[15][0] + 2. * mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2], 2. * mv_xsi[15][0] * (1. - mv_xsi[15][2]) * mv_xsi[15][2], 2. * mv_xsi[15][0] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][2]) },
	  { (0.5 - 2. * mv_xsi[15][0]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2], 0., -0.5 * mv_xsi[15][0] * (2. * mv_xsi[15][0] - 1.) * (1. - 2. * mv_xsi[15][2]) },
	  { 2. * mv_xsi[15][1] * (1. - mv_xsi[15][2]) * mv_xsi[15][2], (-2. + 2. * mv_xsi[15][0] + 4. * mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2], 2. * mv_xsi[15][1] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][2]) },
	  { -2. * mv_xsi[15][1] * mv_xsi[15][2] * (1. - mv_xsi[15][2]), -2. * mv_xsi[15][0] * mv_xsi[15][2] * (1. - mv_xsi[15][2]), -2. * mv_xsi[15][0] * mv_xsi[15][1] * (1. - 2. * mv_xsi[15][2]) },
	  { 0., (0.5 - 2. * mv_xsi[15][1]) * (1. - mv_xsi[15][2]) * mv_xsi[15][2], -0.5 * mv_xsi[15][1] * (2. * mv_xsi[15][1] - 1.) * (1. - 2. * mv_xsi[15][2]) },
	  { (3. - 4. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), (3. - 4. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (2. * mv_xsi[15][2]) },
	  { (-4. + 8. * mv_xsi[15][0] + 4. * mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), 4. * mv_xsi[15][0] * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), 4. * mv_xsi[15][0] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (2. * mv_xsi[15][2]) },
	  { (1. - 4. * mv_xsi[15][0]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), 0., mv_xsi[15][0] * (2. * mv_xsi[15][0] - 1.) * (-2. * mv_xsi[15][2]) },
	  { 4. * mv_xsi[15][1] * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), (-4. + 4. * mv_xsi[15][0] + 8. * mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), 4. * mv_xsi[15][1] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (2. * mv_xsi[15][2]) },
	  { -4. * mv_xsi[15][1] * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), -4. * mv_xsi[15][0] * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), -4. * mv_xsi[15][0] * mv_xsi[15][1] * (2. * mv_xsi[15][2]) },
	  { 0., (1. - 4. * mv_xsi[15][1]) * (-1. + mv_xsi[15][2] * mv_xsi[15][2]), mv_xsi[15][1] * (2. * mv_xsi[15][1] - 1.) * (-2. * mv_xsi[15][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2], -0.5 * (3. - 4. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2], -0.5 * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. - 2. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (1. + 2. * mv_xsi[15][2]) },
	  { (2. - 4. * mv_xsi[15][0] - 2. * mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2], -2. * mv_xsi[15][0] * (1. + mv_xsi[15][2]) * mv_xsi[15][2], -2. * mv_xsi[15][0] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. + 2. * mv_xsi[15][2]) },
	  { (-0.5 + 2. * mv_xsi[15][0]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2], 0., 0.5 * mv_xsi[15][0] * (2. * mv_xsi[15][0] - 1.) * (1. + 2. * mv_xsi[15][2]) },
	  { -2. * mv_xsi[15][1] * (1. + mv_xsi[15][2]) * mv_xsi[15][2], (2. - 2. * mv_xsi[15][0] - 4. * mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2], -2. * mv_xsi[15][1] * (-1. + mv_xsi[15][0] + mv_xsi[15][1]) * (1. + 2. * mv_xsi[15][2]) },
	  { 2. * mv_xsi[15][1] * mv_xsi[15][2] * (1. + mv_xsi[15][2]), 2. * mv_xsi[15][0] * mv_xsi[15][2] * (1. + mv_xsi[15][2]), 2. * mv_xsi[15][0] * mv_xsi[15][1] * (1. + 2. * mv_xsi[15][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[15][1]) * (1. + mv_xsi[15][2]) * mv_xsi[15][2], 0.5 * mv_xsi[15][1] * (2. * mv_xsi[15][1] - 1.) * (1. + 2. * mv_xsi[15][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2], 0.5 * (3. - 4. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2], 0.5 * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][2]) },
	  { (-2. + 4. * mv_xsi[16][0] + 2. * mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2], 2. * mv_xsi[16][0] * (1. - mv_xsi[16][2]) * mv_xsi[16][2], 2. * mv_xsi[16][0] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][2]) },
	  { (0.5 - 2. * mv_xsi[16][0]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2], 0., -0.5 * mv_xsi[16][0] * (2. * mv_xsi[16][0] - 1.) * (1. - 2. * mv_xsi[16][2]) },
	  { 2. * mv_xsi[16][1] * (1. - mv_xsi[16][2]) * mv_xsi[16][2], (-2. + 2. * mv_xsi[16][0] + 4. * mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2], 2. * mv_xsi[16][1] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][2]) },
	  { -2. * mv_xsi[16][1] * mv_xsi[16][2] * (1. - mv_xsi[16][2]), -2. * mv_xsi[16][0] * mv_xsi[16][2] * (1. - mv_xsi[16][2]), -2. * mv_xsi[16][0] * mv_xsi[16][1] * (1. - 2. * mv_xsi[16][2]) },
	  { 0., (0.5 - 2. * mv_xsi[16][1]) * (1. - mv_xsi[16][2]) * mv_xsi[16][2], -0.5 * mv_xsi[16][1] * (2. * mv_xsi[16][1] - 1.) * (1. - 2. * mv_xsi[16][2]) },
	  { (3. - 4. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), (3. - 4. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (2. * mv_xsi[16][2]) },
	  { (-4. + 8. * mv_xsi[16][0] + 4. * mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), 4. * mv_xsi[16][0] * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), 4. * mv_xsi[16][0] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (2. * mv_xsi[16][2]) },
	  { (1. - 4. * mv_xsi[16][0]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), 0., mv_xsi[16][0] * (2. * mv_xsi[16][0] - 1.) * (-2. * mv_xsi[16][2]) },
	  { 4. * mv_xsi[16][1] * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), (-4. + 4. * mv_xsi[16][0] + 8. * mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), 4. * mv_xsi[16][1] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (2. * mv_xsi[16][2]) },
	  { -4. * mv_xsi[16][1] * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), -4. * mv_xsi[16][0] * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), -4. * mv_xsi[16][0] * mv_xsi[16][1] * (2. * mv_xsi[16][2]) },
	  { 0., (1. - 4. * mv_xsi[16][1]) * (-1. + mv_xsi[16][2] * mv_xsi[16][2]), mv_xsi[16][1] * (2. * mv_xsi[16][1] - 1.) * (-2. * mv_xsi[16][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2], -0.5 * (3. - 4. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2], -0.5 * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. - 2. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (1. + 2. * mv_xsi[16][2]) },
	  { (2. - 4. * mv_xsi[16][0] - 2. * mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2], -2. * mv_xsi[16][0] * (1. + mv_xsi[16][2]) * mv_xsi[16][2], -2. * mv_xsi[16][0] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. + 2. * mv_xsi[16][2]) },
	  { (-0.5 + 2. * mv_xsi[16][0]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2], 0., 0.5 * mv_xsi[16][0] * (2. * mv_xsi[16][0] - 1.) * (1. + 2. * mv_xsi[16][2]) },
	  { -2. * mv_xsi[16][1] * (1. + mv_xsi[16][2]) * mv_xsi[16][2], (2. - 2. * mv_xsi[16][0] - 4. * mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2], -2. * mv_xsi[16][1] * (-1. + mv_xsi[16][0] + mv_xsi[16][1]) * (1. + 2. * mv_xsi[16][2]) },
	  { 2. * mv_xsi[16][1] * mv_xsi[16][2] * (1. + mv_xsi[16][2]), 2. * mv_xsi[16][0] * mv_xsi[16][2] * (1. + mv_xsi[16][2]), 2. * mv_xsi[16][0] * mv_xsi[16][1] * (1. + 2. * mv_xsi[16][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[16][1]) * (1. + mv_xsi[16][2]) * mv_xsi[16][2], 0.5 * mv_xsi[16][1] * (2. * mv_xsi[16][1] - 1.) * (1. + 2. * mv_xsi[16][2]) } },

	{ { 0.5 * (3. - 4. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2], 0.5 * (3. - 4. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2], 0.5 * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][2]) },
	  { (-2. + 4. * mv_xsi[17][0] + 2. * mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2], 2. * mv_xsi[17][0] * (1. - mv_xsi[17][2]) * mv_xsi[17][2], 2. * mv_xsi[17][0] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][2]) },
	  { (0.5 - 2. * mv_xsi[17][0]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2], 0., -0.5 * mv_xsi[17][0] * (2. * mv_xsi[17][0] - 1.) * (1. - 2. * mv_xsi[17][2]) },
	  { 2. * mv_xsi[17][1] * (1. - mv_xsi[17][2]) * mv_xsi[17][2], (-2. + 2. * mv_xsi[17][0] + 4. * mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2], 2. * mv_xsi[17][1] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][2]) },
	  { -2. * mv_xsi[17][1] * mv_xsi[17][2] * (1. - mv_xsi[17][2]), -2. * mv_xsi[17][0] * mv_xsi[17][2] * (1. - mv_xsi[17][2]), -2. * mv_xsi[17][0] * mv_xsi[17][1] * (1. - 2. * mv_xsi[17][2]) },
	  { 0., (0.5 - 2. * mv_xsi[17][1]) * (1. - mv_xsi[17][2]) * mv_xsi[17][2], -0.5 * mv_xsi[17][1] * (2. * mv_xsi[17][1] - 1.) * (1. - 2. * mv_xsi[17][2]) },
	  { (3. - 4. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), (3. - 4. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (2. * mv_xsi[17][2]) },
	  { (-4. + 8. * mv_xsi[17][0] + 4. * mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), 4. * mv_xsi[17][0] * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), 4. * mv_xsi[17][0] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (2. * mv_xsi[17][2]) },
	  { (1. - 4. * mv_xsi[17][0]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), 0., mv_xsi[17][0] * (2. * mv_xsi[17][0] - 1.) * (-2. * mv_xsi[17][2]) },
	  { 4. * mv_xsi[17][1] * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), (-4. + 4. * mv_xsi[17][0] + 8. * mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), 4. * mv_xsi[17][1] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (2. * mv_xsi[17][2]) },
	  { -4. * mv_xsi[17][1] * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), -4. * mv_xsi[17][0] * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), -4. * mv_xsi[17][0] * mv_xsi[17][1] * (2. * mv_xsi[17][2]) },
	  { 0., (1. - 4. * mv_xsi[17][1]) * (-1. + mv_xsi[17][2] * mv_xsi[17][2]), mv_xsi[17][1] * (2. * mv_xsi[17][1] - 1.) * (-2. * mv_xsi[17][2]) },
	  { -0.5 * (3. - 4. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2], -0.5 * (3. - 4. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2], -0.5 * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. - 2. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (1. + 2. * mv_xsi[17][2]) },
	  { (2. - 4. * mv_xsi[17][0] - 2. * mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2], -2. * mv_xsi[17][0] * (1. + mv_xsi[17][2]) * mv_xsi[17][2], -2. * mv_xsi[17][0] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. + 2. * mv_xsi[17][2]) },
	  { (-0.5 + 2. * mv_xsi[17][0]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2], 0., 0.5 * mv_xsi[17][0] * (2. * mv_xsi[17][0] - 1.) * (1. + 2. * mv_xsi[17][2]) },
	  { -2. * mv_xsi[17][1] * (1. + mv_xsi[17][2]) * mv_xsi[17][2], (2. - 2. * mv_xsi[17][0] - 4. * mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2], -2. * mv_xsi[17][1] * (-1. + mv_xsi[17][0] + mv_xsi[17][1]) * (1. + 2. * mv_xsi[17][2]) },
	  { 2. * mv_xsi[17][1] * mv_xsi[17][2] * (1. + mv_xsi[17][2]), 2. * mv_xsi[17][0] * mv_xsi[17][2] * (1. + mv_xsi[17][2]), 2. * mv_xsi[17][0] * mv_xsi[17][1] * (1. + 2. * mv_xsi[17][2]) },
	  { 0., (-0.5 + 2. * mv_xsi[17][1]) * (1. + mv_xsi[17][2]) * mv_xsi[17][2], 0.5 * mv_xsi[17][1] * (2. * mv_xsi[17][1] - 1.) * (1. + 2. * mv_xsi[17][2]) } } };
