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
// Tetrahedral solid element, with quadratic interpolation functions.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"
#include "IntegrationPoint.h"

namespace O2P2 {
	namespace Geom {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Tet10
			  *
			  * @brief Tetrahedral quadratic element with 10 nodes.
			  * @details Solid element, with quadratic interpolation functions, tetrahedral shape.
			  * Options for integration points: 4, 10, 11, 14, 15 and 24.
			  * @image html Elem_Tet10.png height=300
			  */
			class Elem_Tet10 : public ElemSolid
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Tet10() = delete;

			protected:
				/** Constructor for tetrahedral quadratic elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Tet10(std::shared_ptr<O2P2::Geom::Material>& Material)
					: ElemSolid(Material) { }

			public:
				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 2 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " "
						<< this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " "
						<< this->mv_Conect[4]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " "
						<< this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[6]->mv_index + add << " "
						<< this->mv_Conect[7]->mv_index + add << " " << this->mv_Conect[9]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " "
						<< this->mv_Conect[5]->mv_index + add << " " << this->mv_Conect[6]->mv_index + add << " "
						<< this->mv_Conect[8]->mv_index + add << " " << this->mv_Conect[9]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[4]->mv_index + add << " "
						<< this->mv_Conect[5]->mv_index + add << " " << this->mv_Conect[7]->mv_index + add << " "
						<< this->mv_Conect[8]->mv_index + add << " " << this->mv_Conect[9]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 2 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " "
						<< (4 + add) << " " << (5 + add) << " " << (6 + add) << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " "
						<< (7 + add) << " " << (8 + add) << " " << (10 + add) << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << (1 + add) << " " << (4 + add) << " " << (6 + add) << " "
						<< (7 + add) << " " << (9 + add) << " " << (10 + add) << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "2 2 " << (3 + add) << " " << (5 + add) << " " << (6 + add) << " "
						<< (8 + add) << " " << (9 + add) << " " << (10 + add) << " "
						<< this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Evaluates shape function in the point.
				std::vector<double> getShapeFcOnPoint(const double* Point) override;

				// Evaluates the derivative of shape function in the point.
				std::vector<double> getShapeDerivOnPoint(const double* Point) override;

				// Returns the number of nodes of current element.
				const int getNumNodes() const override { return mv_numNodes; }

				// Returns the number of faces of current element.
				const int getNumFaces() const override { return mv_numFaces; }

				// Verifies dimensionless coordinates from input - if it is immersed on the element.
				inline bool evaluateXsi(const double* xsi) override {
					std::array<double, mv_ElDim + 1> new_xsi = {};

					for (int i = 0; i < mv_ElDim; ++i) {
						new_xsi.at(i) = *(xsi + i);
						new_xsi.at(mv_ElDim) -= *(xsi + i);
					}
					new_xsi.at(mv_ElDim) += 1.;

					const auto [min, max] = std::minmax_element(new_xsi.begin(), new_xsi.end());
					if (*max < 1.000001 && *min > -0.000001) return true;
					return false;
				}

			private:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			protected:
				/** @brief Number of Nodes */
				static const int mv_numNodes{ 10 };

				/** @brief Number of Faces */
				static const int mv_numFaces{ 4 };
			};


			/**
			  * @class Elem_Tet10_IP
			  *
			  * @brief Tetrahedral quadratic element with 10 nodes.
			  * @details Solid element, with quadratic interpolation functions, tetrahedral shape.
			  *
			  * @tparam nIP Number of integration points. Must be: 4, 10, 11, 14, 15 or 24.
			  */
			template<int nIP>
			class Elem_Tet10_IP : public Elem_Tet10
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Tet10_IP() = delete;

			public:
				/** Constructor for tetrahedral quadratic elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Tet10_IP(std::shared_ptr<O2P2::Geom::Material>& Material)
					: Elem_Tet10(Material) { }

				// Return a vector with values on the integration points currently known in the element' nodes.
				std::vector<double> getValueOnIPs(const double* value) override;

				// Returns a pointer to the first element of the shape functions (with size [nIP][mv_numNodes]).
				double const* getShapeFc() const override { return &mv_Psi[0][0]; }

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][mv_numNodes][mv_ElDim]).
				double const* getShapeDerivative() const override { return &mv_DPsi[0][0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return mv_weight; }

				// Returns the number of integration points of current element.
				const int getNumIP() const override { return nIP; }

			private:
				// Weights for numerical integration
				static const double* mv_weight;

				// Shape functions
				static const double mv_Psi[nIP][mv_numNodes];

				// Shape functions derivative
				static const double mv_DPsi[nIP][mv_numNodes][mv_ElDim];
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
inline std::vector<double> O2P2::Geom::Elem::Elem_Tet10::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(10);

	mi_Psi.at(0) = (1. - 2. * Point[0] - 2. * Point[1] - 2. * Point[2]) * (1. - Point[0] - Point[1] - Point[2]);
	mi_Psi.at(1) = 4. * Point[0] * (1. - Point[0] - Point[1] - Point[2]);
	mi_Psi.at(2) = (2. * Point[0] - 1.) * Point[0];
	mi_Psi.at(3) = 4. * Point[1] * (1. - Point[0] - Point[1] - Point[2]);
	mi_Psi.at(4) = 4. * Point[0] * Point[1];
	mi_Psi.at(5) = (2. * Point[1] - 1.) * Point[1];
	mi_Psi.at(6) = 4. * Point[2] * (1. - Point[0] - Point[1] - Point[2]);
	mi_Psi.at(7) = 4. * Point[0] * Point[2];
	mi_Psi.at(8) = 4. * Point[1] * Point[2];
	mi_Psi.at(9) = (2. * Point[2] - 1.) * Point[2];

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Geom::Elem::Elem_Tet10::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(10 * 3);

	mi_DPsi.at(0) = -3. + 4. * (Point[0] + Point[1] + Point[2]);
	mi_DPsi.at(1) = 4. - 4. * (2. * Point[0] + Point[1] + Point[2]);
	mi_DPsi.at(2) = 4. * Point[0] - 1.;
	mi_DPsi.at(3) = -4. * Point[1];
	mi_DPsi.at(4) = 4. * Point[1];
	mi_DPsi.at(5) = 0.;
	mi_DPsi.at(6) = -4. * Point[2];
	mi_DPsi.at(7) = 4. * Point[2];
	mi_DPsi.at(8) = 0.;
	mi_DPsi.at(9) = 0.;

	mi_DPsi.at(10) = -3. + 4. * (Point[0] + Point[1] + Point[2]);
	mi_DPsi.at(11) = -4. * Point[0];
	mi_DPsi.at(12) = 0.;
	mi_DPsi.at(13) = 4. - 4. * (Point[0] + 2. * Point[1] + Point[2]);
	mi_DPsi.at(14) = 4. * Point[0];
	mi_DPsi.at(15) = 4. * Point[1] - 1.;
	mi_DPsi.at(16) = -4. * Point[2];
	mi_DPsi.at(17) = 0.;
	mi_DPsi.at(18) = 4. * Point[2];
	mi_DPsi.at(19) = 0.;

	mi_DPsi.at(20) = -3. + 4. * (Point[0] + Point[1] + Point[2]);
	mi_DPsi.at(21) = -4. * Point[0];
	mi_DPsi.at(22) = 0.;
	mi_DPsi.at(23) = -4. * Point[1];
	mi_DPsi.at(24) = 0.;
	mi_DPsi.at(25) = 0.;
	mi_DPsi.at(26) = 4. - 4. * (Point[0] + Point[1] + 2. * Point[2]);
	mi_DPsi.at(27) = 4. * Point[0];
	mi_DPsi.at(28) = 4. * Point[1];
	mi_DPsi.at(29) = 4. * Point[2] - 1.;

	return mi_DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Geom::Elem::Elem_Tet10::setGeomProperties() {

	const int nVertices = 4;
	const int mi_Dim = mv_Conect.at(0)->getDIM();	// Dimensionality of vector space (2D or 3D)

	mv_Centroid = std::make_unique<double[]>(mi_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Geom::Node*, nVertices> vertices;
	vertices[0] = mv_Conect[0].get();
	vertices[1] = mv_Conect[2].get();
	vertices[2] = mv_Conect[5].get();
	vertices[3] = mv_Conect[9].get();

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
template<int nIP>
inline std::vector<double> O2P2::Geom::Elem::Elem_Tet10_IP<nIP>::getValueOnIPs(const double* value) {

	// return value
	std::vector<double> mi_valueOnIp(nIP, 0.);

	for (int i = 0; i < nIP; i++) {
		for (int j = 0; j < this->mv_numNodes; j++) {
			mi_valueOnIp.at(i) += value[i] * mv_Psi[i][j];
		}
	}

	return mi_valueOnIp;
};


// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
template<> const double* O2P2::Geom::Elem::Elem_Tet10_IP<4>::mv_weight = &Hammer3D::Wg_4P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tet10_IP<10>::mv_weight = &Hammer3D::Wg_10P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tet10_IP<11>::mv_weight = &Hammer3D::Wg_11P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tet10_IP<14>::mv_weight = &Hammer3D::Wg_14P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tet10_IP<15>::mv_weight = &Hammer3D::Wg_15P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tet10_IP<24>::mv_weight = &Hammer3D::Wg_24P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<4>::mv_Psi[4][mv_numNodes] = {
	{ (1. - 2. * Hammer3D::Qsi_4P[0][0] - 2. * Hammer3D::Qsi_4P[0][1] - 2. * Hammer3D::Qsi_4P[0][2]) * (1. - Hammer3D::Qsi_4P[0][0] - Hammer3D::Qsi_4P[0][1] - Hammer3D::Qsi_4P[0][2]),
	  4. * Hammer3D::Qsi_4P[0][0] * (1. - Hammer3D::Qsi_4P[0][0] - Hammer3D::Qsi_4P[0][1] - Hammer3D::Qsi_4P[0][2]),
	  (2. * Hammer3D::Qsi_4P[0][0] - 1.) * Hammer3D::Qsi_4P[0][0],
	  4. * Hammer3D::Qsi_4P[0][1] * (1. - Hammer3D::Qsi_4P[0][0] - Hammer3D::Qsi_4P[0][1] - Hammer3D::Qsi_4P[0][2]),
	  4. * Hammer3D::Qsi_4P[0][0] * Hammer3D::Qsi_4P[0][1],
	  (2. * Hammer3D::Qsi_4P[0][1] - 1.) * Hammer3D::Qsi_4P[0][1],
	  4. * Hammer3D::Qsi_4P[0][2] * (1. - Hammer3D::Qsi_4P[0][0] - Hammer3D::Qsi_4P[0][1] - Hammer3D::Qsi_4P[0][2]),
	  4. * Hammer3D::Qsi_4P[0][0] * Hammer3D::Qsi_4P[0][2],
	  4. * Hammer3D::Qsi_4P[0][1] * Hammer3D::Qsi_4P[0][2],
	  (2. * Hammer3D::Qsi_4P[0][2] - 1.) * Hammer3D::Qsi_4P[0][2] },

	{ (1. - 2. * Hammer3D::Qsi_4P[1][0] - 2. * Hammer3D::Qsi_4P[1][1] - 2. * Hammer3D::Qsi_4P[1][2]) * (1. - Hammer3D::Qsi_4P[1][0] - Hammer3D::Qsi_4P[1][1] - Hammer3D::Qsi_4P[1][2]),
	  4. * Hammer3D::Qsi_4P[1][0] * (1. - Hammer3D::Qsi_4P[1][0] - Hammer3D::Qsi_4P[1][1] - Hammer3D::Qsi_4P[1][2]),
	  (2. * Hammer3D::Qsi_4P[1][0] - 1.) * Hammer3D::Qsi_4P[1][0],
	  4. * Hammer3D::Qsi_4P[1][1] * (1. - Hammer3D::Qsi_4P[1][0] - Hammer3D::Qsi_4P[1][1] - Hammer3D::Qsi_4P[1][2]),
	  4. * Hammer3D::Qsi_4P[1][0] * Hammer3D::Qsi_4P[1][1],
	  (2. * Hammer3D::Qsi_4P[1][1] - 1.) * Hammer3D::Qsi_4P[1][1],
	  4. * Hammer3D::Qsi_4P[1][2] * (1. - Hammer3D::Qsi_4P[1][0] - Hammer3D::Qsi_4P[1][1] - Hammer3D::Qsi_4P[1][2]),
	  4. * Hammer3D::Qsi_4P[1][0] * Hammer3D::Qsi_4P[1][2],
	  4. * Hammer3D::Qsi_4P[1][1] * Hammer3D::Qsi_4P[1][2],
	  (2. * Hammer3D::Qsi_4P[1][2] - 1.) * Hammer3D::Qsi_4P[1][2] },

	{ (1. - 2. * Hammer3D::Qsi_4P[2][0] - 2. * Hammer3D::Qsi_4P[2][1] - 2. * Hammer3D::Qsi_4P[2][2]) * (1. - Hammer3D::Qsi_4P[2][0] - Hammer3D::Qsi_4P[2][1] - Hammer3D::Qsi_4P[2][2]),
	  4. * Hammer3D::Qsi_4P[2][0] * (1. - Hammer3D::Qsi_4P[2][0] - Hammer3D::Qsi_4P[2][1] - Hammer3D::Qsi_4P[2][2]),
	  (2. * Hammer3D::Qsi_4P[2][0] - 1.) * Hammer3D::Qsi_4P[2][0],
	  4. * Hammer3D::Qsi_4P[2][1] * (1. - Hammer3D::Qsi_4P[2][0] - Hammer3D::Qsi_4P[2][1] - Hammer3D::Qsi_4P[2][2]),
	  4. * Hammer3D::Qsi_4P[2][0] * Hammer3D::Qsi_4P[2][1],
	  (2. * Hammer3D::Qsi_4P[2][1] - 1.) * Hammer3D::Qsi_4P[2][1],
	  4. * Hammer3D::Qsi_4P[2][2] * (1. - Hammer3D::Qsi_4P[2][0] - Hammer3D::Qsi_4P[2][1] - Hammer3D::Qsi_4P[2][2]),
	  4. * Hammer3D::Qsi_4P[2][0] * Hammer3D::Qsi_4P[2][2],
	  4. * Hammer3D::Qsi_4P[2][1] * Hammer3D::Qsi_4P[2][2],
	  (2. * Hammer3D::Qsi_4P[2][2] - 1.) * Hammer3D::Qsi_4P[2][2] },

	{ (1. - 2. * Hammer3D::Qsi_4P[3][0] - 2. * Hammer3D::Qsi_4P[3][1] - 2. * Hammer3D::Qsi_4P[3][2]) * (1. - Hammer3D::Qsi_4P[3][0] - Hammer3D::Qsi_4P[3][1] - Hammer3D::Qsi_4P[3][2]),
	  4. * Hammer3D::Qsi_4P[3][0] * (1. - Hammer3D::Qsi_4P[3][0] - Hammer3D::Qsi_4P[3][1] - Hammer3D::Qsi_4P[3][2]),
	  (2. * Hammer3D::Qsi_4P[3][0] - 1.) * Hammer3D::Qsi_4P[3][0],
	  4. * Hammer3D::Qsi_4P[3][1] * (1. - Hammer3D::Qsi_4P[3][0] - Hammer3D::Qsi_4P[3][1] - Hammer3D::Qsi_4P[3][2]),
	  4. * Hammer3D::Qsi_4P[3][0] * Hammer3D::Qsi_4P[3][1],
	  (2. * Hammer3D::Qsi_4P[3][1] - 1.) * Hammer3D::Qsi_4P[3][1],
	  4. * Hammer3D::Qsi_4P[3][2] * (1. - Hammer3D::Qsi_4P[3][0] - Hammer3D::Qsi_4P[3][1] - Hammer3D::Qsi_4P[3][2]),
	  4. * Hammer3D::Qsi_4P[3][0] * Hammer3D::Qsi_4P[3][2],
	  4. * Hammer3D::Qsi_4P[3][1] * Hammer3D::Qsi_4P[3][2],
	  (2. * Hammer3D::Qsi_4P[3][2] - 1.) * Hammer3D::Qsi_4P[3][2] } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<10>::mv_Psi[10][mv_numNodes] = {
	{ (1. - 2. * Hammer3D::Qsi_10P[0][0] - 2. * Hammer3D::Qsi_10P[0][1] - 2. * Hammer3D::Qsi_10P[0][2]) * (1. - Hammer3D::Qsi_10P[0][0] - Hammer3D::Qsi_10P[0][1] - Hammer3D::Qsi_10P[0][2]),
	  4. * Hammer3D::Qsi_10P[0][0] * (1. - Hammer3D::Qsi_10P[0][0] - Hammer3D::Qsi_10P[0][1] - Hammer3D::Qsi_10P[0][2]),
	  (2. * Hammer3D::Qsi_10P[0][0] - 1.) * Hammer3D::Qsi_10P[0][0],
	  4. * Hammer3D::Qsi_10P[0][1] * (1. - Hammer3D::Qsi_10P[0][0] - Hammer3D::Qsi_10P[0][1] - Hammer3D::Qsi_10P[0][2]),
	  4. * Hammer3D::Qsi_10P[0][0] * Hammer3D::Qsi_10P[0][1],
	  (2. * Hammer3D::Qsi_10P[0][1] - 1.) * Hammer3D::Qsi_10P[0][1],
	  4. * Hammer3D::Qsi_10P[0][2] * (1. - Hammer3D::Qsi_10P[0][0] - Hammer3D::Qsi_10P[0][1] - Hammer3D::Qsi_10P[0][2]),
	  4. * Hammer3D::Qsi_10P[0][0] * Hammer3D::Qsi_10P[0][2],
	  4. * Hammer3D::Qsi_10P[0][1] * Hammer3D::Qsi_10P[0][2],
	  (2. * Hammer3D::Qsi_10P[0][2] - 1.) * Hammer3D::Qsi_10P[0][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[1][0] - 2. * Hammer3D::Qsi_10P[1][1] - 2. * Hammer3D::Qsi_10P[1][2]) * (1. - Hammer3D::Qsi_10P[1][0] - Hammer3D::Qsi_10P[1][1] - Hammer3D::Qsi_10P[1][2]),
	  4. * Hammer3D::Qsi_10P[1][0] * (1. - Hammer3D::Qsi_10P[1][0] - Hammer3D::Qsi_10P[1][1] - Hammer3D::Qsi_10P[1][2]),
	  (2. * Hammer3D::Qsi_10P[1][0] - 1.) * Hammer3D::Qsi_10P[1][0],
	  4. * Hammer3D::Qsi_10P[1][1] * (1. - Hammer3D::Qsi_10P[1][0] - Hammer3D::Qsi_10P[1][1] - Hammer3D::Qsi_10P[1][2]),
	  4. * Hammer3D::Qsi_10P[1][0] * Hammer3D::Qsi_10P[1][1],
	  (2. * Hammer3D::Qsi_10P[1][1] - 1.) * Hammer3D::Qsi_10P[1][1],
	  4. * Hammer3D::Qsi_10P[1][2] * (1. - Hammer3D::Qsi_10P[1][0] - Hammer3D::Qsi_10P[1][1] - Hammer3D::Qsi_10P[1][2]),
	  4. * Hammer3D::Qsi_10P[1][0] * Hammer3D::Qsi_10P[1][2],
	  4. * Hammer3D::Qsi_10P[1][1] * Hammer3D::Qsi_10P[1][2],
	  (2. * Hammer3D::Qsi_10P[1][2] - 1.) * Hammer3D::Qsi_10P[1][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[2][0] - 2. * Hammer3D::Qsi_10P[2][1] - 2. * Hammer3D::Qsi_10P[2][2]) * (1. - Hammer3D::Qsi_10P[2][0] - Hammer3D::Qsi_10P[2][1] - Hammer3D::Qsi_10P[2][2]),
	  4. * Hammer3D::Qsi_10P[2][0] * (1. - Hammer3D::Qsi_10P[2][0] - Hammer3D::Qsi_10P[2][1] - Hammer3D::Qsi_10P[2][2]),
	  (2. * Hammer3D::Qsi_10P[2][0] - 1.) * Hammer3D::Qsi_10P[2][0],
	  4. * Hammer3D::Qsi_10P[2][1] * (1. - Hammer3D::Qsi_10P[2][0] - Hammer3D::Qsi_10P[2][1] - Hammer3D::Qsi_10P[2][2]),
	  4. * Hammer3D::Qsi_10P[2][0] * Hammer3D::Qsi_10P[2][1],
	  (2. * Hammer3D::Qsi_10P[2][1] - 1.) * Hammer3D::Qsi_10P[2][1],
	  4. * Hammer3D::Qsi_10P[2][2] * (1. - Hammer3D::Qsi_10P[2][0] - Hammer3D::Qsi_10P[2][1] - Hammer3D::Qsi_10P[2][2]),
	  4. * Hammer3D::Qsi_10P[2][0] * Hammer3D::Qsi_10P[2][2],
	  4. * Hammer3D::Qsi_10P[2][1] * Hammer3D::Qsi_10P[2][2],
	  (2. * Hammer3D::Qsi_10P[2][2] - 1.) * Hammer3D::Qsi_10P[2][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[3][0] - 2. * Hammer3D::Qsi_10P[3][1] - 2. * Hammer3D::Qsi_10P[3][2]) * (1. - Hammer3D::Qsi_10P[3][0] - Hammer3D::Qsi_10P[3][1] - Hammer3D::Qsi_10P[3][2]),
	  4. * Hammer3D::Qsi_10P[3][0] * (1. - Hammer3D::Qsi_10P[3][0] - Hammer3D::Qsi_10P[3][1] - Hammer3D::Qsi_10P[3][2]),
	  (2. * Hammer3D::Qsi_10P[3][0] - 1.) * Hammer3D::Qsi_10P[3][0],
	  4. * Hammer3D::Qsi_10P[3][1] * (1. - Hammer3D::Qsi_10P[3][0] - Hammer3D::Qsi_10P[3][1] - Hammer3D::Qsi_10P[3][2]),
	  4. * Hammer3D::Qsi_10P[3][0] * Hammer3D::Qsi_10P[3][1],
	  (2. * Hammer3D::Qsi_10P[3][1] - 1.) * Hammer3D::Qsi_10P[3][1],
	  4. * Hammer3D::Qsi_10P[3][2] * (1. - Hammer3D::Qsi_10P[3][0] - Hammer3D::Qsi_10P[3][1] - Hammer3D::Qsi_10P[3][2]),
	  4. * Hammer3D::Qsi_10P[3][0] * Hammer3D::Qsi_10P[3][2],
	  4. * Hammer3D::Qsi_10P[3][1] * Hammer3D::Qsi_10P[3][2],
	  (2. * Hammer3D::Qsi_10P[3][2] - 1.) * Hammer3D::Qsi_10P[3][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[4][0] - 2. * Hammer3D::Qsi_10P[4][1] - 2. * Hammer3D::Qsi_10P[4][2]) * (1. - Hammer3D::Qsi_10P[4][0] - Hammer3D::Qsi_10P[4][1] - Hammer3D::Qsi_10P[4][2]),
	  4. * Hammer3D::Qsi_10P[4][0] * (1. - Hammer3D::Qsi_10P[4][0] - Hammer3D::Qsi_10P[4][1] - Hammer3D::Qsi_10P[4][2]),
	  (2. * Hammer3D::Qsi_10P[4][0] - 1.) * Hammer3D::Qsi_10P[4][0],
	  4. * Hammer3D::Qsi_10P[4][1] * (1. - Hammer3D::Qsi_10P[4][0] - Hammer3D::Qsi_10P[4][1] - Hammer3D::Qsi_10P[4][2]),
	  4. * Hammer3D::Qsi_10P[4][0] * Hammer3D::Qsi_10P[4][1],
	  (2. * Hammer3D::Qsi_10P[4][1] - 1.) * Hammer3D::Qsi_10P[4][1],
	  4. * Hammer3D::Qsi_10P[4][2] * (1. - Hammer3D::Qsi_10P[4][0] - Hammer3D::Qsi_10P[4][1] - Hammer3D::Qsi_10P[4][2]),
	  4. * Hammer3D::Qsi_10P[4][0] * Hammer3D::Qsi_10P[4][2],
	  4. * Hammer3D::Qsi_10P[4][1] * Hammer3D::Qsi_10P[4][2],
	  (2. * Hammer3D::Qsi_10P[4][2] - 1.) * Hammer3D::Qsi_10P[4][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[5][0] - 2. * Hammer3D::Qsi_10P[5][1] - 2. * Hammer3D::Qsi_10P[5][2]) * (1. - Hammer3D::Qsi_10P[5][0] - Hammer3D::Qsi_10P[5][1] - Hammer3D::Qsi_10P[5][2]),
	  4. * Hammer3D::Qsi_10P[5][0] * (1. - Hammer3D::Qsi_10P[5][0] - Hammer3D::Qsi_10P[5][1] - Hammer3D::Qsi_10P[5][2]),
	  (2. * Hammer3D::Qsi_10P[5][0] - 1.) * Hammer3D::Qsi_10P[5][0],
	  4. * Hammer3D::Qsi_10P[5][1] * (1. - Hammer3D::Qsi_10P[5][0] - Hammer3D::Qsi_10P[5][1] - Hammer3D::Qsi_10P[5][2]),
	  4. * Hammer3D::Qsi_10P[5][0] * Hammer3D::Qsi_10P[5][1],
	  (2. * Hammer3D::Qsi_10P[5][1] - 1.) * Hammer3D::Qsi_10P[5][1],
	  4. * Hammer3D::Qsi_10P[5][2] * (1. - Hammer3D::Qsi_10P[5][0] - Hammer3D::Qsi_10P[5][1] - Hammer3D::Qsi_10P[5][2]),
	  4. * Hammer3D::Qsi_10P[5][0] * Hammer3D::Qsi_10P[5][2],
	  4. * Hammer3D::Qsi_10P[5][1] * Hammer3D::Qsi_10P[5][2],
	  (2. * Hammer3D::Qsi_10P[5][2] - 1.) * Hammer3D::Qsi_10P[5][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[6][0] - 2. * Hammer3D::Qsi_10P[6][1] - 2. * Hammer3D::Qsi_10P[6][2]) * (1. - Hammer3D::Qsi_10P[6][0] - Hammer3D::Qsi_10P[6][1] - Hammer3D::Qsi_10P[6][2]),
	  4. * Hammer3D::Qsi_10P[6][0] * (1. - Hammer3D::Qsi_10P[6][0] - Hammer3D::Qsi_10P[6][1] - Hammer3D::Qsi_10P[6][2]),
	  (2. * Hammer3D::Qsi_10P[6][0] - 1.) * Hammer3D::Qsi_10P[6][0],
	  4. * Hammer3D::Qsi_10P[6][1] * (1. - Hammer3D::Qsi_10P[6][0] - Hammer3D::Qsi_10P[6][1] - Hammer3D::Qsi_10P[6][2]),
	  4. * Hammer3D::Qsi_10P[6][0] * Hammer3D::Qsi_10P[6][1],
	  (2. * Hammer3D::Qsi_10P[6][1] - 1.) * Hammer3D::Qsi_10P[6][1],
	  4. * Hammer3D::Qsi_10P[6][2] * (1. - Hammer3D::Qsi_10P[6][0] - Hammer3D::Qsi_10P[6][1] - Hammer3D::Qsi_10P[6][2]),
	  4. * Hammer3D::Qsi_10P[6][0] * Hammer3D::Qsi_10P[6][2],
	  4. * Hammer3D::Qsi_10P[6][1] * Hammer3D::Qsi_10P[6][2],
	  (2. * Hammer3D::Qsi_10P[6][2] - 1.) * Hammer3D::Qsi_10P[6][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[7][0] - 2. * Hammer3D::Qsi_10P[7][1] - 2. * Hammer3D::Qsi_10P[7][2]) * (1. - Hammer3D::Qsi_10P[7][0] - Hammer3D::Qsi_10P[7][1] - Hammer3D::Qsi_10P[7][2]),
	  4. * Hammer3D::Qsi_10P[7][0] * (1. - Hammer3D::Qsi_10P[7][0] - Hammer3D::Qsi_10P[7][1] - Hammer3D::Qsi_10P[7][2]),
	  (2. * Hammer3D::Qsi_10P[7][0] - 1.) * Hammer3D::Qsi_10P[7][0],
	  4. * Hammer3D::Qsi_10P[7][1] * (1. - Hammer3D::Qsi_10P[7][0] - Hammer3D::Qsi_10P[7][1] - Hammer3D::Qsi_10P[7][2]),
	  4. * Hammer3D::Qsi_10P[7][0] * Hammer3D::Qsi_10P[7][1],
	  (2. * Hammer3D::Qsi_10P[7][1] - 1.) * Hammer3D::Qsi_10P[7][1],
	  4. * Hammer3D::Qsi_10P[7][2] * (1. - Hammer3D::Qsi_10P[7][0] - Hammer3D::Qsi_10P[7][1] - Hammer3D::Qsi_10P[7][2]),
	  4. * Hammer3D::Qsi_10P[7][0] * Hammer3D::Qsi_10P[7][2],
	  4. * Hammer3D::Qsi_10P[7][1] * Hammer3D::Qsi_10P[7][2],
	  (2. * Hammer3D::Qsi_10P[7][2] - 1.) * Hammer3D::Qsi_10P[7][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[8][0] - 2. * Hammer3D::Qsi_10P[8][1] - 2. * Hammer3D::Qsi_10P[8][2]) * (1. - Hammer3D::Qsi_10P[8][0] - Hammer3D::Qsi_10P[8][1] - Hammer3D::Qsi_10P[8][2]),
	  4. * Hammer3D::Qsi_10P[8][0] * (1. - Hammer3D::Qsi_10P[8][0] - Hammer3D::Qsi_10P[8][1] - Hammer3D::Qsi_10P[8][2]),
	  (2. * Hammer3D::Qsi_10P[8][0] - 1.) * Hammer3D::Qsi_10P[8][0],
	  4. * Hammer3D::Qsi_10P[8][1] * (1. - Hammer3D::Qsi_10P[8][0] - Hammer3D::Qsi_10P[8][1] - Hammer3D::Qsi_10P[8][2]),
	  4. * Hammer3D::Qsi_10P[8][0] * Hammer3D::Qsi_10P[8][1],
	  (2. * Hammer3D::Qsi_10P[8][1] - 1.) * Hammer3D::Qsi_10P[8][1],
	  4. * Hammer3D::Qsi_10P[8][2] * (1. - Hammer3D::Qsi_10P[8][0] - Hammer3D::Qsi_10P[8][1] - Hammer3D::Qsi_10P[8][2]),
	  4. * Hammer3D::Qsi_10P[8][0] * Hammer3D::Qsi_10P[8][2],
	  4. * Hammer3D::Qsi_10P[8][1] * Hammer3D::Qsi_10P[8][2],
	  (2. * Hammer3D::Qsi_10P[8][2] - 1.) * Hammer3D::Qsi_10P[8][2] },

	{ (1. - 2. * Hammer3D::Qsi_10P[9][0] - 2. * Hammer3D::Qsi_10P[9][1] - 2. * Hammer3D::Qsi_10P[9][2]) * (1. - Hammer3D::Qsi_10P[9][0] - Hammer3D::Qsi_10P[9][1] - Hammer3D::Qsi_10P[9][2]),
	  4. * Hammer3D::Qsi_10P[9][0] * (1. - Hammer3D::Qsi_10P[9][0] - Hammer3D::Qsi_10P[9][1] - Hammer3D::Qsi_10P[9][2]),
	  (2. * Hammer3D::Qsi_10P[9][0] - 1.) * Hammer3D::Qsi_10P[9][0],
	  4. * Hammer3D::Qsi_10P[9][1] * (1. - Hammer3D::Qsi_10P[9][0] - Hammer3D::Qsi_10P[9][1] - Hammer3D::Qsi_10P[9][2]),
	  4. * Hammer3D::Qsi_10P[9][0] * Hammer3D::Qsi_10P[9][1],
	  (2. * Hammer3D::Qsi_10P[9][1] - 1.) * Hammer3D::Qsi_10P[9][1],
	  4. * Hammer3D::Qsi_10P[9][2] * (1. - Hammer3D::Qsi_10P[9][0] - Hammer3D::Qsi_10P[9][1] - Hammer3D::Qsi_10P[9][2]),
	  4. * Hammer3D::Qsi_10P[9][0] * Hammer3D::Qsi_10P[9][2],
	  4. * Hammer3D::Qsi_10P[9][1] * Hammer3D::Qsi_10P[9][2],
	  (2. * Hammer3D::Qsi_10P[9][2] - 1.) * Hammer3D::Qsi_10P[9][2] } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<11>::mv_Psi[11][mv_numNodes] = {
	{ (1. - 2. * Hammer3D::Qsi_11P[0][0] - 2. * Hammer3D::Qsi_11P[0][1] - 2. * Hammer3D::Qsi_11P[0][2]) * (1. - Hammer3D::Qsi_11P[0][0] - Hammer3D::Qsi_11P[0][1] - Hammer3D::Qsi_11P[0][2]),
	  4. * Hammer3D::Qsi_11P[0][0] * (1. - Hammer3D::Qsi_11P[0][0] - Hammer3D::Qsi_11P[0][1] - Hammer3D::Qsi_11P[0][2]),
	  (2. * Hammer3D::Qsi_11P[0][0] - 1.) * Hammer3D::Qsi_11P[0][0],
	  4. * Hammer3D::Qsi_11P[0][1] * (1. - Hammer3D::Qsi_11P[0][0] - Hammer3D::Qsi_11P[0][1] - Hammer3D::Qsi_11P[0][2]),
	  4. * Hammer3D::Qsi_11P[0][0] * Hammer3D::Qsi_11P[0][1],
	  (2. * Hammer3D::Qsi_11P[0][1] - 1.) * Hammer3D::Qsi_11P[0][1],
	  4. * Hammer3D::Qsi_11P[0][2] * (1. - Hammer3D::Qsi_11P[0][0] - Hammer3D::Qsi_11P[0][1] - Hammer3D::Qsi_11P[0][2]),
	  4. * Hammer3D::Qsi_11P[0][0] * Hammer3D::Qsi_11P[0][2],
	  4. * Hammer3D::Qsi_11P[0][1] * Hammer3D::Qsi_11P[0][2],
	  (2. * Hammer3D::Qsi_11P[0][2] - 1.) * Hammer3D::Qsi_11P[0][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[1][0] - 2. * Hammer3D::Qsi_11P[1][1] - 2. * Hammer3D::Qsi_11P[1][2]) * (1. - Hammer3D::Qsi_11P[1][0] - Hammer3D::Qsi_11P[1][1] - Hammer3D::Qsi_11P[1][2]),
	  4. * Hammer3D::Qsi_11P[1][0] * (1. - Hammer3D::Qsi_11P[1][0] - Hammer3D::Qsi_11P[1][1] - Hammer3D::Qsi_11P[1][2]),
	  (2. * Hammer3D::Qsi_11P[1][0] - 1.) * Hammer3D::Qsi_11P[1][0],
	  4. * Hammer3D::Qsi_11P[1][1] * (1. - Hammer3D::Qsi_11P[1][0] - Hammer3D::Qsi_11P[1][1] - Hammer3D::Qsi_11P[1][2]),
	  4. * Hammer3D::Qsi_11P[1][0] * Hammer3D::Qsi_11P[1][1],
	  (2. * Hammer3D::Qsi_11P[1][1] - 1.) * Hammer3D::Qsi_11P[1][1],
	  4. * Hammer3D::Qsi_11P[1][2] * (1. - Hammer3D::Qsi_11P[1][0] - Hammer3D::Qsi_11P[1][1] - Hammer3D::Qsi_11P[1][2]),
	  4. * Hammer3D::Qsi_11P[1][0] * Hammer3D::Qsi_11P[1][2],
	  4. * Hammer3D::Qsi_11P[1][1] * Hammer3D::Qsi_11P[1][2],
	  (2. * Hammer3D::Qsi_11P[1][2] - 1.) * Hammer3D::Qsi_11P[1][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[2][0] - 2. * Hammer3D::Qsi_11P[2][1] - 2. * Hammer3D::Qsi_11P[2][2]) * (1. - Hammer3D::Qsi_11P[2][0] - Hammer3D::Qsi_11P[2][1] - Hammer3D::Qsi_11P[2][2]),
	  4. * Hammer3D::Qsi_11P[2][0] * (1. - Hammer3D::Qsi_11P[2][0] - Hammer3D::Qsi_11P[2][1] - Hammer3D::Qsi_11P[2][2]),
	  (2. * Hammer3D::Qsi_11P[2][0] - 1.) * Hammer3D::Qsi_11P[2][0],
	  4. * Hammer3D::Qsi_11P[2][1] * (1. - Hammer3D::Qsi_11P[2][0] - Hammer3D::Qsi_11P[2][1] - Hammer3D::Qsi_11P[2][2]),
	  4. * Hammer3D::Qsi_11P[2][0] * Hammer3D::Qsi_11P[2][1],
	  (2. * Hammer3D::Qsi_11P[2][1] - 1.) * Hammer3D::Qsi_11P[2][1],
	  4. * Hammer3D::Qsi_11P[2][2] * (1. - Hammer3D::Qsi_11P[2][0] - Hammer3D::Qsi_11P[2][1] - Hammer3D::Qsi_11P[2][2]),
	  4. * Hammer3D::Qsi_11P[2][0] * Hammer3D::Qsi_11P[2][2],
	  4. * Hammer3D::Qsi_11P[2][1] * Hammer3D::Qsi_11P[2][2],
	  (2. * Hammer3D::Qsi_11P[2][2] - 1.) * Hammer3D::Qsi_11P[2][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[3][0] - 2. * Hammer3D::Qsi_11P[3][1] - 2. * Hammer3D::Qsi_11P[3][2]) * (1. - Hammer3D::Qsi_11P[3][0] - Hammer3D::Qsi_11P[3][1] - Hammer3D::Qsi_11P[3][2]),
	  4. * Hammer3D::Qsi_11P[3][0] * (1. - Hammer3D::Qsi_11P[3][0] - Hammer3D::Qsi_11P[3][1] - Hammer3D::Qsi_11P[3][2]),
	  (2. * Hammer3D::Qsi_11P[3][0] - 1.) * Hammer3D::Qsi_11P[3][0],
	  4. * Hammer3D::Qsi_11P[3][1] * (1. - Hammer3D::Qsi_11P[3][0] - Hammer3D::Qsi_11P[3][1] - Hammer3D::Qsi_11P[3][2]),
	  4. * Hammer3D::Qsi_11P[3][0] * Hammer3D::Qsi_11P[3][1],
	  (2. * Hammer3D::Qsi_11P[3][1] - 1.) * Hammer3D::Qsi_11P[3][1],
	  4. * Hammer3D::Qsi_11P[3][2] * (1. - Hammer3D::Qsi_11P[3][0] - Hammer3D::Qsi_11P[3][1] - Hammer3D::Qsi_11P[3][2]),
	  4. * Hammer3D::Qsi_11P[3][0] * Hammer3D::Qsi_11P[3][2],
	  4. * Hammer3D::Qsi_11P[3][1] * Hammer3D::Qsi_11P[3][2],
	  (2. * Hammer3D::Qsi_11P[3][2] - 1.) * Hammer3D::Qsi_11P[3][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[4][0] - 2. * Hammer3D::Qsi_11P[4][1] - 2. * Hammer3D::Qsi_11P[4][2]) * (1. - Hammer3D::Qsi_11P[4][0] - Hammer3D::Qsi_11P[4][1] - Hammer3D::Qsi_11P[4][2]),
	  4. * Hammer3D::Qsi_11P[4][0] * (1. - Hammer3D::Qsi_11P[4][0] - Hammer3D::Qsi_11P[4][1] - Hammer3D::Qsi_11P[4][2]),
	  (2. * Hammer3D::Qsi_11P[4][0] - 1.) * Hammer3D::Qsi_11P[4][0],
	  4. * Hammer3D::Qsi_11P[4][1] * (1. - Hammer3D::Qsi_11P[4][0] - Hammer3D::Qsi_11P[4][1] - Hammer3D::Qsi_11P[4][2]),
	  4. * Hammer3D::Qsi_11P[4][0] * Hammer3D::Qsi_11P[4][1],
	  (2. * Hammer3D::Qsi_11P[4][1] - 1.) * Hammer3D::Qsi_11P[4][1],
	  4. * Hammer3D::Qsi_11P[4][2] * (1. - Hammer3D::Qsi_11P[4][0] - Hammer3D::Qsi_11P[4][1] - Hammer3D::Qsi_11P[4][2]),
	  4. * Hammer3D::Qsi_11P[4][0] * Hammer3D::Qsi_11P[4][2],
	  4. * Hammer3D::Qsi_11P[4][1] * Hammer3D::Qsi_11P[4][2],
	  (2. * Hammer3D::Qsi_11P[4][2] - 1.) * Hammer3D::Qsi_11P[4][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[5][0] - 2. * Hammer3D::Qsi_11P[5][1] - 2. * Hammer3D::Qsi_11P[5][2]) * (1. - Hammer3D::Qsi_11P[5][0] - Hammer3D::Qsi_11P[5][1] - Hammer3D::Qsi_11P[5][2]),
	  4. * Hammer3D::Qsi_11P[5][0] * (1. - Hammer3D::Qsi_11P[5][0] - Hammer3D::Qsi_11P[5][1] - Hammer3D::Qsi_11P[5][2]),
	  (2. * Hammer3D::Qsi_11P[5][0] - 1.) * Hammer3D::Qsi_11P[5][0],
	  4. * Hammer3D::Qsi_11P[5][1] * (1. - Hammer3D::Qsi_11P[5][0] - Hammer3D::Qsi_11P[5][1] - Hammer3D::Qsi_11P[5][2]),
	  4. * Hammer3D::Qsi_11P[5][0] * Hammer3D::Qsi_11P[5][1],
	  (2. * Hammer3D::Qsi_11P[5][1] - 1.) * Hammer3D::Qsi_11P[5][1],
	  4. * Hammer3D::Qsi_11P[5][2] * (1. - Hammer3D::Qsi_11P[5][0] - Hammer3D::Qsi_11P[5][1] - Hammer3D::Qsi_11P[5][2]),
	  4. * Hammer3D::Qsi_11P[5][0] * Hammer3D::Qsi_11P[5][2],
	  4. * Hammer3D::Qsi_11P[5][1] * Hammer3D::Qsi_11P[5][2],
	  (2. * Hammer3D::Qsi_11P[5][2] - 1.) * Hammer3D::Qsi_11P[5][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[6][0] - 2. * Hammer3D::Qsi_11P[6][1] - 2. * Hammer3D::Qsi_11P[6][2]) * (1. - Hammer3D::Qsi_11P[6][0] - Hammer3D::Qsi_11P[6][1] - Hammer3D::Qsi_11P[6][2]),
	  4. * Hammer3D::Qsi_11P[6][0] * (1. - Hammer3D::Qsi_11P[6][0] - Hammer3D::Qsi_11P[6][1] - Hammer3D::Qsi_11P[6][2]),
	  (2. * Hammer3D::Qsi_11P[6][0] - 1.) * Hammer3D::Qsi_11P[6][0],
	  4. * Hammer3D::Qsi_11P[6][1] * (1. - Hammer3D::Qsi_11P[6][0] - Hammer3D::Qsi_11P[6][1] - Hammer3D::Qsi_11P[6][2]),
	  4. * Hammer3D::Qsi_11P[6][0] * Hammer3D::Qsi_11P[6][1],
	  (2. * Hammer3D::Qsi_11P[6][1] - 1.) * Hammer3D::Qsi_11P[6][1],
	  4. * Hammer3D::Qsi_11P[6][2] * (1. - Hammer3D::Qsi_11P[6][0] - Hammer3D::Qsi_11P[6][1] - Hammer3D::Qsi_11P[6][2]),
	  4. * Hammer3D::Qsi_11P[6][0] * Hammer3D::Qsi_11P[6][2],
	  4. * Hammer3D::Qsi_11P[6][1] * Hammer3D::Qsi_11P[6][2],
	  (2. * Hammer3D::Qsi_11P[6][2] - 1.) * Hammer3D::Qsi_11P[6][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[7][0] - 2. * Hammer3D::Qsi_11P[7][1] - 2. * Hammer3D::Qsi_11P[7][2]) * (1. - Hammer3D::Qsi_11P[7][0] - Hammer3D::Qsi_11P[7][1] - Hammer3D::Qsi_11P[7][2]),
	  4. * Hammer3D::Qsi_11P[7][0] * (1. - Hammer3D::Qsi_11P[7][0] - Hammer3D::Qsi_11P[7][1] - Hammer3D::Qsi_11P[7][2]),
	  (2. * Hammer3D::Qsi_11P[7][0] - 1.) * Hammer3D::Qsi_11P[7][0],
	  4. * Hammer3D::Qsi_11P[7][1] * (1. - Hammer3D::Qsi_11P[7][0] - Hammer3D::Qsi_11P[7][1] - Hammer3D::Qsi_11P[7][2]),
	  4. * Hammer3D::Qsi_11P[7][0] * Hammer3D::Qsi_11P[7][1],
	  (2. * Hammer3D::Qsi_11P[7][1] - 1.) * Hammer3D::Qsi_11P[7][1],
	  4. * Hammer3D::Qsi_11P[7][2] * (1. - Hammer3D::Qsi_11P[7][0] - Hammer3D::Qsi_11P[7][1] - Hammer3D::Qsi_11P[7][2]),
	  4. * Hammer3D::Qsi_11P[7][0] * Hammer3D::Qsi_11P[7][2],
	  4. * Hammer3D::Qsi_11P[7][1] * Hammer3D::Qsi_11P[7][2],
	  (2. * Hammer3D::Qsi_11P[7][2] - 1.) * Hammer3D::Qsi_11P[7][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[8][0] - 2. * Hammer3D::Qsi_11P[8][1] - 2. * Hammer3D::Qsi_11P[8][2]) * (1. - Hammer3D::Qsi_11P[8][0] - Hammer3D::Qsi_11P[8][1] - Hammer3D::Qsi_11P[8][2]),
	  4. * Hammer3D::Qsi_11P[8][0] * (1. - Hammer3D::Qsi_11P[8][0] - Hammer3D::Qsi_11P[8][1] - Hammer3D::Qsi_11P[8][2]),
	  (2. * Hammer3D::Qsi_11P[8][0] - 1.) * Hammer3D::Qsi_11P[8][0],
	  4. * Hammer3D::Qsi_11P[8][1] * (1. - Hammer3D::Qsi_11P[8][0] - Hammer3D::Qsi_11P[8][1] - Hammer3D::Qsi_11P[8][2]),
	  4. * Hammer3D::Qsi_11P[8][0] * Hammer3D::Qsi_11P[8][1],
	  (2. * Hammer3D::Qsi_11P[8][1] - 1.) * Hammer3D::Qsi_11P[8][1],
	  4. * Hammer3D::Qsi_11P[8][2] * (1. - Hammer3D::Qsi_11P[8][0] - Hammer3D::Qsi_11P[8][1] - Hammer3D::Qsi_11P[8][2]),
	  4. * Hammer3D::Qsi_11P[8][0] * Hammer3D::Qsi_11P[8][2],
	  4. * Hammer3D::Qsi_11P[8][1] * Hammer3D::Qsi_11P[8][2],
	  (2. * Hammer3D::Qsi_11P[8][2] - 1.) * Hammer3D::Qsi_11P[8][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[9][0] - 2. * Hammer3D::Qsi_11P[9][1] - 2. * Hammer3D::Qsi_11P[9][2]) * (1. - Hammer3D::Qsi_11P[9][0] - Hammer3D::Qsi_11P[9][1] - Hammer3D::Qsi_11P[9][2]),
	  4. * Hammer3D::Qsi_11P[9][0] * (1. - Hammer3D::Qsi_11P[9][0] - Hammer3D::Qsi_11P[9][1] - Hammer3D::Qsi_11P[9][2]),
	  (2. * Hammer3D::Qsi_11P[9][0] - 1.) * Hammer3D::Qsi_11P[9][0],
	  4. * Hammer3D::Qsi_11P[9][1] * (1. - Hammer3D::Qsi_11P[9][0] - Hammer3D::Qsi_11P[9][1] - Hammer3D::Qsi_11P[9][2]),
	  4. * Hammer3D::Qsi_11P[9][0] * Hammer3D::Qsi_11P[9][1],
	  (2. * Hammer3D::Qsi_11P[9][1] - 1.) * Hammer3D::Qsi_11P[9][1],
	  4. * Hammer3D::Qsi_11P[9][2] * (1. - Hammer3D::Qsi_11P[9][0] - Hammer3D::Qsi_11P[9][1] - Hammer3D::Qsi_11P[9][2]),
	  4. * Hammer3D::Qsi_11P[9][0] * Hammer3D::Qsi_11P[9][2],
	  4. * Hammer3D::Qsi_11P[9][1] * Hammer3D::Qsi_11P[9][2],
	  (2. * Hammer3D::Qsi_11P[9][2] - 1.) * Hammer3D::Qsi_11P[9][2] },

	{ (1. - 2. * Hammer3D::Qsi_11P[10][0] - 2. * Hammer3D::Qsi_11P[10][1] - 2. * Hammer3D::Qsi_11P[10][2])* (1. - Hammer3D::Qsi_11P[10][0] - Hammer3D::Qsi_11P[10][1] - Hammer3D::Qsi_11P[10][2]),
		4. * Hammer3D::Qsi_11P[10][0] * (1. - Hammer3D::Qsi_11P[10][0] - Hammer3D::Qsi_11P[10][1] - Hammer3D::Qsi_11P[10][2]),
		(2. * Hammer3D::Qsi_11P[10][0] - 1.)* Hammer3D::Qsi_11P[10][0],
		4. * Hammer3D::Qsi_11P[10][1] * (1. - Hammer3D::Qsi_11P[10][0] - Hammer3D::Qsi_11P[10][1] - Hammer3D::Qsi_11P[10][2]),
		4. * Hammer3D::Qsi_11P[10][0] * Hammer3D::Qsi_11P[10][1],
		(2. * Hammer3D::Qsi_11P[10][1] - 1.)* Hammer3D::Qsi_11P[10][1],
		4. * Hammer3D::Qsi_11P[10][2] * (1. - Hammer3D::Qsi_11P[10][0] - Hammer3D::Qsi_11P[10][1] - Hammer3D::Qsi_11P[10][2]),
		4. * Hammer3D::Qsi_11P[10][0] * Hammer3D::Qsi_11P[10][2],
		4. * Hammer3D::Qsi_11P[10][1] * Hammer3D::Qsi_11P[10][2],
		(2. * Hammer3D::Qsi_11P[10][2] - 1.)* Hammer3D::Qsi_11P[10][2] } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<14>::mv_Psi[14][mv_numNodes] = {
	{ (1. - 2. * Hammer3D::Qsi_14P[0][0] - 2. * Hammer3D::Qsi_14P[0][1] - 2. * Hammer3D::Qsi_14P[0][2]) * (1. - Hammer3D::Qsi_14P[0][0] - Hammer3D::Qsi_14P[0][1] - Hammer3D::Qsi_14P[0][2]),
	  4. * Hammer3D::Qsi_14P[0][0] * (1. - Hammer3D::Qsi_14P[0][0] - Hammer3D::Qsi_14P[0][1] - Hammer3D::Qsi_14P[0][2]),
	  (2. * Hammer3D::Qsi_14P[0][0] - 1.) * Hammer3D::Qsi_14P[0][0],
	  4. * Hammer3D::Qsi_14P[0][1] * (1. - Hammer3D::Qsi_14P[0][0] - Hammer3D::Qsi_14P[0][1] - Hammer3D::Qsi_14P[0][2]),
	  4. * Hammer3D::Qsi_14P[0][0] * Hammer3D::Qsi_14P[0][1],
	  (2. * Hammer3D::Qsi_14P[0][1] - 1.) * Hammer3D::Qsi_14P[0][1],
	  4. * Hammer3D::Qsi_14P[0][2] * (1. - Hammer3D::Qsi_14P[0][0] - Hammer3D::Qsi_14P[0][1] - Hammer3D::Qsi_14P[0][2]),
	  4. * Hammer3D::Qsi_14P[0][0] * Hammer3D::Qsi_14P[0][2],
	  4. * Hammer3D::Qsi_14P[0][1] * Hammer3D::Qsi_14P[0][2],
	  (2. * Hammer3D::Qsi_14P[0][2] - 1.) * Hammer3D::Qsi_14P[0][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[1][0] - 2. * Hammer3D::Qsi_14P[1][1] - 2. * Hammer3D::Qsi_14P[1][2]) * (1. - Hammer3D::Qsi_14P[1][0] - Hammer3D::Qsi_14P[1][1] - Hammer3D::Qsi_14P[1][2]),
	  4. * Hammer3D::Qsi_14P[1][0] * (1. - Hammer3D::Qsi_14P[1][0] - Hammer3D::Qsi_14P[1][1] - Hammer3D::Qsi_14P[1][2]),
	  (2. * Hammer3D::Qsi_14P[1][0] - 1.) * Hammer3D::Qsi_14P[1][0],
	  4. * Hammer3D::Qsi_14P[1][1] * (1. - Hammer3D::Qsi_14P[1][0] - Hammer3D::Qsi_14P[1][1] - Hammer3D::Qsi_14P[1][2]),
	  4. * Hammer3D::Qsi_14P[1][0] * Hammer3D::Qsi_14P[1][1],
	  (2. * Hammer3D::Qsi_14P[1][1] - 1.) * Hammer3D::Qsi_14P[1][1],
	  4. * Hammer3D::Qsi_14P[1][2] * (1. - Hammer3D::Qsi_14P[1][0] - Hammer3D::Qsi_14P[1][1] - Hammer3D::Qsi_14P[1][2]),
	  4. * Hammer3D::Qsi_14P[1][0] * Hammer3D::Qsi_14P[1][2],
	  4. * Hammer3D::Qsi_14P[1][1] * Hammer3D::Qsi_14P[1][2],
	  (2. * Hammer3D::Qsi_14P[1][2] - 1.) * Hammer3D::Qsi_14P[1][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[2][0] - 2. * Hammer3D::Qsi_14P[2][1] - 2. * Hammer3D::Qsi_14P[2][2]) * (1. - Hammer3D::Qsi_14P[2][0] - Hammer3D::Qsi_14P[2][1] - Hammer3D::Qsi_14P[2][2]),
	  4. * Hammer3D::Qsi_14P[2][0] * (1. - Hammer3D::Qsi_14P[2][0] - Hammer3D::Qsi_14P[2][1] - Hammer3D::Qsi_14P[2][2]),
	  (2. * Hammer3D::Qsi_14P[2][0] - 1.) * Hammer3D::Qsi_14P[2][0],
	  4. * Hammer3D::Qsi_14P[2][1] * (1. - Hammer3D::Qsi_14P[2][0] - Hammer3D::Qsi_14P[2][1] - Hammer3D::Qsi_14P[2][2]),
	  4. * Hammer3D::Qsi_14P[2][0] * Hammer3D::Qsi_14P[2][1],
	  (2. * Hammer3D::Qsi_14P[2][1] - 1.) * Hammer3D::Qsi_14P[2][1],
	  4. * Hammer3D::Qsi_14P[2][2] * (1. - Hammer3D::Qsi_14P[2][0] - Hammer3D::Qsi_14P[2][1] - Hammer3D::Qsi_14P[2][2]),
	  4. * Hammer3D::Qsi_14P[2][0] * Hammer3D::Qsi_14P[2][2],
	  4. * Hammer3D::Qsi_14P[2][1] * Hammer3D::Qsi_14P[2][2],
	  (2. * Hammer3D::Qsi_14P[2][2] - 1.) * Hammer3D::Qsi_14P[2][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[3][0] - 2. * Hammer3D::Qsi_14P[3][1] - 2. * Hammer3D::Qsi_14P[3][2]) * (1. - Hammer3D::Qsi_14P[3][0] - Hammer3D::Qsi_14P[3][1] - Hammer3D::Qsi_14P[3][2]),
	  4. * Hammer3D::Qsi_14P[3][0] * (1. - Hammer3D::Qsi_14P[3][0] - Hammer3D::Qsi_14P[3][1] - Hammer3D::Qsi_14P[3][2]),
	  (2. * Hammer3D::Qsi_14P[3][0] - 1.) * Hammer3D::Qsi_14P[3][0],
	  4. * Hammer3D::Qsi_14P[3][1] * (1. - Hammer3D::Qsi_14P[3][0] - Hammer3D::Qsi_14P[3][1] - Hammer3D::Qsi_14P[3][2]),
	  4. * Hammer3D::Qsi_14P[3][0] * Hammer3D::Qsi_14P[3][1],
	  (2. * Hammer3D::Qsi_14P[3][1] - 1.) * Hammer3D::Qsi_14P[3][1],
	  4. * Hammer3D::Qsi_14P[3][2] * (1. - Hammer3D::Qsi_14P[3][0] - Hammer3D::Qsi_14P[3][1] - Hammer3D::Qsi_14P[3][2]),
	  4. * Hammer3D::Qsi_14P[3][0] * Hammer3D::Qsi_14P[3][2],
	  4. * Hammer3D::Qsi_14P[3][1] * Hammer3D::Qsi_14P[3][2],
	  (2. * Hammer3D::Qsi_14P[3][2] - 1.) * Hammer3D::Qsi_14P[3][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[4][0] - 2. * Hammer3D::Qsi_14P[4][1] - 2. * Hammer3D::Qsi_14P[4][2]) * (1. - Hammer3D::Qsi_14P[4][0] - Hammer3D::Qsi_14P[4][1] - Hammer3D::Qsi_14P[4][2]),
	  4. * Hammer3D::Qsi_14P[4][0] * (1. - Hammer3D::Qsi_14P[4][0] - Hammer3D::Qsi_14P[4][1] - Hammer3D::Qsi_14P[4][2]),
	  (2. * Hammer3D::Qsi_14P[4][0] - 1.) * Hammer3D::Qsi_14P[4][0],
	  4. * Hammer3D::Qsi_14P[4][1] * (1. - Hammer3D::Qsi_14P[4][0] - Hammer3D::Qsi_14P[4][1] - Hammer3D::Qsi_14P[4][2]),
	  4. * Hammer3D::Qsi_14P[4][0] * Hammer3D::Qsi_14P[4][1],
	  (2. * Hammer3D::Qsi_14P[4][1] - 1.) * Hammer3D::Qsi_14P[4][1],
	  4. * Hammer3D::Qsi_14P[4][2] * (1. - Hammer3D::Qsi_14P[4][0] - Hammer3D::Qsi_14P[4][1] - Hammer3D::Qsi_14P[4][2]),
	  4. * Hammer3D::Qsi_14P[4][0] * Hammer3D::Qsi_14P[4][2],
	  4. * Hammer3D::Qsi_14P[4][1] * Hammer3D::Qsi_14P[4][2],
	  (2. * Hammer3D::Qsi_14P[4][2] - 1.) * Hammer3D::Qsi_14P[4][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[5][0] - 2. * Hammer3D::Qsi_14P[5][1] - 2. * Hammer3D::Qsi_14P[5][2]) * (1. - Hammer3D::Qsi_14P[5][0] - Hammer3D::Qsi_14P[5][1] - Hammer3D::Qsi_14P[5][2]),
	  4. * Hammer3D::Qsi_14P[5][0] * (1. - Hammer3D::Qsi_14P[5][0] - Hammer3D::Qsi_14P[5][1] - Hammer3D::Qsi_14P[5][2]),
	  (2. * Hammer3D::Qsi_14P[5][0] - 1.) * Hammer3D::Qsi_14P[5][0],
	  4. * Hammer3D::Qsi_14P[5][1] * (1. - Hammer3D::Qsi_14P[5][0] - Hammer3D::Qsi_14P[5][1] - Hammer3D::Qsi_14P[5][2]),
	  4. * Hammer3D::Qsi_14P[5][0] * Hammer3D::Qsi_14P[5][1],
	  (2. * Hammer3D::Qsi_14P[5][1] - 1.) * Hammer3D::Qsi_14P[5][1],
	  4. * Hammer3D::Qsi_14P[5][2] * (1. - Hammer3D::Qsi_14P[5][0] - Hammer3D::Qsi_14P[5][1] - Hammer3D::Qsi_14P[5][2]),
	  4. * Hammer3D::Qsi_14P[5][0] * Hammer3D::Qsi_14P[5][2],
	  4. * Hammer3D::Qsi_14P[5][1] * Hammer3D::Qsi_14P[5][2],
	  (2. * Hammer3D::Qsi_14P[5][2] - 1.) * Hammer3D::Qsi_14P[5][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[6][0] - 2. * Hammer3D::Qsi_14P[6][1] - 2. * Hammer3D::Qsi_14P[6][2]) * (1. - Hammer3D::Qsi_14P[6][0] - Hammer3D::Qsi_14P[6][1] - Hammer3D::Qsi_14P[6][2]),
	  4. * Hammer3D::Qsi_14P[6][0] * (1. - Hammer3D::Qsi_14P[6][0] - Hammer3D::Qsi_14P[6][1] - Hammer3D::Qsi_14P[6][2]),
	  (2. * Hammer3D::Qsi_14P[6][0] - 1.) * Hammer3D::Qsi_14P[6][0],
	  4. * Hammer3D::Qsi_14P[6][1] * (1. - Hammer3D::Qsi_14P[6][0] - Hammer3D::Qsi_14P[6][1] - Hammer3D::Qsi_14P[6][2]),
	  4. * Hammer3D::Qsi_14P[6][0] * Hammer3D::Qsi_14P[6][1],
	  (2. * Hammer3D::Qsi_14P[6][1] - 1.) * Hammer3D::Qsi_14P[6][1],
	  4. * Hammer3D::Qsi_14P[6][2] * (1. - Hammer3D::Qsi_14P[6][0] - Hammer3D::Qsi_14P[6][1] - Hammer3D::Qsi_14P[6][2]),
	  4. * Hammer3D::Qsi_14P[6][0] * Hammer3D::Qsi_14P[6][2],
	  4. * Hammer3D::Qsi_14P[6][1] * Hammer3D::Qsi_14P[6][2],
	  (2. * Hammer3D::Qsi_14P[6][2] - 1.) * Hammer3D::Qsi_14P[6][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[7][0] - 2. * Hammer3D::Qsi_14P[7][1] - 2. * Hammer3D::Qsi_14P[7][2]) * (1. - Hammer3D::Qsi_14P[7][0] - Hammer3D::Qsi_14P[7][1] - Hammer3D::Qsi_14P[7][2]),
	  4. * Hammer3D::Qsi_14P[7][0] * (1. - Hammer3D::Qsi_14P[7][0] - Hammer3D::Qsi_14P[7][1] - Hammer3D::Qsi_14P[7][2]),
	  (2. * Hammer3D::Qsi_14P[7][0] - 1.) * Hammer3D::Qsi_14P[7][0],
	  4. * Hammer3D::Qsi_14P[7][1] * (1. - Hammer3D::Qsi_14P[7][0] - Hammer3D::Qsi_14P[7][1] - Hammer3D::Qsi_14P[7][2]),
	  4. * Hammer3D::Qsi_14P[7][0] * Hammer3D::Qsi_14P[7][1],
	  (2. * Hammer3D::Qsi_14P[7][1] - 1.) * Hammer3D::Qsi_14P[7][1],
	  4. * Hammer3D::Qsi_14P[7][2] * (1. - Hammer3D::Qsi_14P[7][0] - Hammer3D::Qsi_14P[7][1] - Hammer3D::Qsi_14P[7][2]),
	  4. * Hammer3D::Qsi_14P[7][0] * Hammer3D::Qsi_14P[7][2],
	  4. * Hammer3D::Qsi_14P[7][1] * Hammer3D::Qsi_14P[7][2],
	  (2. * Hammer3D::Qsi_14P[7][2] - 1.) * Hammer3D::Qsi_14P[7][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[8][0] - 2. * Hammer3D::Qsi_14P[8][1] - 2. * Hammer3D::Qsi_14P[8][2]) * (1. - Hammer3D::Qsi_14P[8][0] - Hammer3D::Qsi_14P[8][1] - Hammer3D::Qsi_14P[8][2]),
	  4. * Hammer3D::Qsi_14P[8][0] * (1. - Hammer3D::Qsi_14P[8][0] - Hammer3D::Qsi_14P[8][1] - Hammer3D::Qsi_14P[8][2]),
	  (2. * Hammer3D::Qsi_14P[8][0] - 1.) * Hammer3D::Qsi_14P[8][0],
	  4. * Hammer3D::Qsi_14P[8][1] * (1. - Hammer3D::Qsi_14P[8][0] - Hammer3D::Qsi_14P[8][1] - Hammer3D::Qsi_14P[8][2]),
	  4. * Hammer3D::Qsi_14P[8][0] * Hammer3D::Qsi_14P[8][1],
	  (2. * Hammer3D::Qsi_14P[8][1] - 1.) * Hammer3D::Qsi_14P[8][1],
	  4. * Hammer3D::Qsi_14P[8][2] * (1. - Hammer3D::Qsi_14P[8][0] - Hammer3D::Qsi_14P[8][1] - Hammer3D::Qsi_14P[8][2]),
	  4. * Hammer3D::Qsi_14P[8][0] * Hammer3D::Qsi_14P[8][2],
	  4. * Hammer3D::Qsi_14P[8][1] * Hammer3D::Qsi_14P[8][2],
	  (2. * Hammer3D::Qsi_14P[8][2] - 1.) * Hammer3D::Qsi_14P[8][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[9][0] - 2. * Hammer3D::Qsi_14P[9][1] - 2. * Hammer3D::Qsi_14P[9][2]) * (1. - Hammer3D::Qsi_14P[9][0] - Hammer3D::Qsi_14P[9][1] - Hammer3D::Qsi_14P[9][2]),
	  4. * Hammer3D::Qsi_14P[9][0] * (1. - Hammer3D::Qsi_14P[9][0] - Hammer3D::Qsi_14P[9][1] - Hammer3D::Qsi_14P[9][2]),
	  (2. * Hammer3D::Qsi_14P[9][0] - 1.) * Hammer3D::Qsi_14P[9][0],
	  4. * Hammer3D::Qsi_14P[9][1] * (1. - Hammer3D::Qsi_14P[9][0] - Hammer3D::Qsi_14P[9][1] - Hammer3D::Qsi_14P[9][2]),
	  4. * Hammer3D::Qsi_14P[9][0] * Hammer3D::Qsi_14P[9][1],
	  (2. * Hammer3D::Qsi_14P[9][1] - 1.) * Hammer3D::Qsi_14P[9][1],
	  4. * Hammer3D::Qsi_14P[9][2] * (1. - Hammer3D::Qsi_14P[9][0] - Hammer3D::Qsi_14P[9][1] - Hammer3D::Qsi_14P[9][2]),
	  4. * Hammer3D::Qsi_14P[9][0] * Hammer3D::Qsi_14P[9][2],
	  4. * Hammer3D::Qsi_14P[9][1] * Hammer3D::Qsi_14P[9][2],
	  (2. * Hammer3D::Qsi_14P[9][2] - 1.) * Hammer3D::Qsi_14P[9][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[10][0] - 2. * Hammer3D::Qsi_14P[10][1] - 2. * Hammer3D::Qsi_14P[10][2]) * (1. - Hammer3D::Qsi_14P[10][0] - Hammer3D::Qsi_14P[10][1] - Hammer3D::Qsi_14P[10][2]),
		4. * Hammer3D::Qsi_14P[10][0] * (1. - Hammer3D::Qsi_14P[10][0] - Hammer3D::Qsi_14P[10][1] - Hammer3D::Qsi_14P[10][2]),
		(2. * Hammer3D::Qsi_14P[10][0] - 1.) * Hammer3D::Qsi_14P[10][0],
		4. * Hammer3D::Qsi_14P[10][1] * (1. - Hammer3D::Qsi_14P[10][0] - Hammer3D::Qsi_14P[10][1] - Hammer3D::Qsi_14P[10][2]),
		4. * Hammer3D::Qsi_14P[10][0] * Hammer3D::Qsi_14P[10][1],
		(2. * Hammer3D::Qsi_14P[10][1] - 1.) * Hammer3D::Qsi_14P[10][1],
		4. * Hammer3D::Qsi_14P[10][2] * (1. - Hammer3D::Qsi_14P[10][0] - Hammer3D::Qsi_14P[10][1] - Hammer3D::Qsi_14P[10][2]),
		4. * Hammer3D::Qsi_14P[10][0] * Hammer3D::Qsi_14P[10][2],
		4. * Hammer3D::Qsi_14P[10][1] * Hammer3D::Qsi_14P[10][2],
		(2. * Hammer3D::Qsi_14P[10][2] - 1.) * Hammer3D::Qsi_14P[10][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[11][0] - 2. * Hammer3D::Qsi_14P[11][1] - 2. * Hammer3D::Qsi_14P[11][2]) * (1. - Hammer3D::Qsi_14P[11][0] - Hammer3D::Qsi_14P[11][1] - Hammer3D::Qsi_14P[11][2]),
		4. * Hammer3D::Qsi_14P[11][0] * (1. - Hammer3D::Qsi_14P[11][0] - Hammer3D::Qsi_14P[11][1] - Hammer3D::Qsi_14P[11][2]),
		(2. * Hammer3D::Qsi_14P[11][0] - 1.) * Hammer3D::Qsi_14P[11][0],
		4. * Hammer3D::Qsi_14P[11][1] * (1. - Hammer3D::Qsi_14P[11][0] - Hammer3D::Qsi_14P[11][1] - Hammer3D::Qsi_14P[11][2]),
		4. * Hammer3D::Qsi_14P[11][0] * Hammer3D::Qsi_14P[11][1],
		(2. * Hammer3D::Qsi_14P[11][1] - 1.) * Hammer3D::Qsi_14P[11][1],
		4. * Hammer3D::Qsi_14P[11][2] * (1. - Hammer3D::Qsi_14P[11][0] - Hammer3D::Qsi_14P[11][1] - Hammer3D::Qsi_14P[11][2]),
		4. * Hammer3D::Qsi_14P[11][0] * Hammer3D::Qsi_14P[11][2],
		4. * Hammer3D::Qsi_14P[11][1] * Hammer3D::Qsi_14P[11][2],
		(2. * Hammer3D::Qsi_14P[11][2] - 1.) * Hammer3D::Qsi_14P[11][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[12][0] - 2. * Hammer3D::Qsi_14P[12][1] - 2. * Hammer3D::Qsi_14P[12][2]) * (1. - Hammer3D::Qsi_14P[12][0] - Hammer3D::Qsi_14P[12][1] - Hammer3D::Qsi_14P[12][2]),
		4. * Hammer3D::Qsi_14P[12][0] * (1. - Hammer3D::Qsi_14P[12][0] - Hammer3D::Qsi_14P[12][1] - Hammer3D::Qsi_14P[12][2]),
		(2. * Hammer3D::Qsi_14P[12][0] - 1.) * Hammer3D::Qsi_14P[12][0],
		4. * Hammer3D::Qsi_14P[12][1] * (1. - Hammer3D::Qsi_14P[12][0] - Hammer3D::Qsi_14P[12][1] - Hammer3D::Qsi_14P[12][2]),
		4. * Hammer3D::Qsi_14P[12][0] * Hammer3D::Qsi_14P[12][1],
		(2. * Hammer3D::Qsi_14P[12][1] - 1.) * Hammer3D::Qsi_14P[12][1],
		4. * Hammer3D::Qsi_14P[12][2] * (1. - Hammer3D::Qsi_14P[12][0] - Hammer3D::Qsi_14P[12][1] - Hammer3D::Qsi_14P[12][2]),
		4. * Hammer3D::Qsi_14P[12][0] * Hammer3D::Qsi_14P[12][2],
		4. * Hammer3D::Qsi_14P[12][1] * Hammer3D::Qsi_14P[12][2],
		(2. * Hammer3D::Qsi_14P[12][2] - 1.) * Hammer3D::Qsi_14P[12][2] },

	{ (1. - 2. * Hammer3D::Qsi_14P[13][0] - 2. * Hammer3D::Qsi_14P[13][1] - 2. * Hammer3D::Qsi_14P[13][2]) * (1. - Hammer3D::Qsi_14P[13][0] - Hammer3D::Qsi_14P[13][1] - Hammer3D::Qsi_14P[13][2]),
		4. * Hammer3D::Qsi_14P[13][0] * (1. - Hammer3D::Qsi_14P[13][0] - Hammer3D::Qsi_14P[13][1] - Hammer3D::Qsi_14P[13][2]),
		(2. * Hammer3D::Qsi_14P[13][0] - 1.) * Hammer3D::Qsi_14P[13][0],
		4. * Hammer3D::Qsi_14P[13][1] * (1. - Hammer3D::Qsi_14P[13][0] - Hammer3D::Qsi_14P[13][1] - Hammer3D::Qsi_14P[13][2]),
		4. * Hammer3D::Qsi_14P[13][0] * Hammer3D::Qsi_14P[13][1],
		(2. * Hammer3D::Qsi_14P[13][1] - 1.) * Hammer3D::Qsi_14P[13][1],
		4. * Hammer3D::Qsi_14P[13][2] * (1. - Hammer3D::Qsi_14P[13][0] - Hammer3D::Qsi_14P[13][1] - Hammer3D::Qsi_14P[13][2]),
		4. * Hammer3D::Qsi_14P[13][0] * Hammer3D::Qsi_14P[13][2],
		4. * Hammer3D::Qsi_14P[13][1] * Hammer3D::Qsi_14P[13][2],
		(2. * Hammer3D::Qsi_14P[13][2] - 1.) * Hammer3D::Qsi_14P[13][2] } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<15>::mv_Psi[15][mv_numNodes] = {
	{ (1. - 2. * Hammer3D::Qsi_15P[0][0] - 2. * Hammer3D::Qsi_15P[0][1] - 2. * Hammer3D::Qsi_15P[0][2]) * (1. - Hammer3D::Qsi_15P[0][0] - Hammer3D::Qsi_15P[0][1] - Hammer3D::Qsi_15P[0][2]),
	  4. * Hammer3D::Qsi_15P[0][0] * (1. - Hammer3D::Qsi_15P[0][0] - Hammer3D::Qsi_15P[0][1] - Hammer3D::Qsi_15P[0][2]),
	  (2. * Hammer3D::Qsi_15P[0][0] - 1.) * Hammer3D::Qsi_15P[0][0],
	  4. * Hammer3D::Qsi_15P[0][1] * (1. - Hammer3D::Qsi_15P[0][0] - Hammer3D::Qsi_15P[0][1] - Hammer3D::Qsi_15P[0][2]),
	  4. * Hammer3D::Qsi_15P[0][0] * Hammer3D::Qsi_15P[0][1],
	  (2. * Hammer3D::Qsi_15P[0][1] - 1.) * Hammer3D::Qsi_15P[0][1],
	  4. * Hammer3D::Qsi_15P[0][2] * (1. - Hammer3D::Qsi_15P[0][0] - Hammer3D::Qsi_15P[0][1] - Hammer3D::Qsi_15P[0][2]),
	  4. * Hammer3D::Qsi_15P[0][0] * Hammer3D::Qsi_15P[0][2],
	  4. * Hammer3D::Qsi_15P[0][1] * Hammer3D::Qsi_15P[0][2],
	  (2. * Hammer3D::Qsi_15P[0][2] - 1.) * Hammer3D::Qsi_15P[0][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[1][0] - 2. * Hammer3D::Qsi_15P[1][1] - 2. * Hammer3D::Qsi_15P[1][2]) * (1. - Hammer3D::Qsi_15P[1][0] - Hammer3D::Qsi_15P[1][1] - Hammer3D::Qsi_15P[1][2]),
	  4. * Hammer3D::Qsi_15P[1][0] * (1. - Hammer3D::Qsi_15P[1][0] - Hammer3D::Qsi_15P[1][1] - Hammer3D::Qsi_15P[1][2]),
	  (2. * Hammer3D::Qsi_15P[1][0] - 1.) * Hammer3D::Qsi_15P[1][0],
	  4. * Hammer3D::Qsi_15P[1][1] * (1. - Hammer3D::Qsi_15P[1][0] - Hammer3D::Qsi_15P[1][1] - Hammer3D::Qsi_15P[1][2]),
	  4. * Hammer3D::Qsi_15P[1][0] * Hammer3D::Qsi_15P[1][1],
	  (2. * Hammer3D::Qsi_15P[1][1] - 1.) * Hammer3D::Qsi_15P[1][1],
	  4. * Hammer3D::Qsi_15P[1][2] * (1. - Hammer3D::Qsi_15P[1][0] - Hammer3D::Qsi_15P[1][1] - Hammer3D::Qsi_15P[1][2]),
	  4. * Hammer3D::Qsi_15P[1][0] * Hammer3D::Qsi_15P[1][2],
	  4. * Hammer3D::Qsi_15P[1][1] * Hammer3D::Qsi_15P[1][2],
	  (2. * Hammer3D::Qsi_15P[1][2] - 1.) * Hammer3D::Qsi_15P[1][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[2][0] - 2. * Hammer3D::Qsi_15P[2][1] - 2. * Hammer3D::Qsi_15P[2][2]) * (1. - Hammer3D::Qsi_15P[2][0] - Hammer3D::Qsi_15P[2][1] - Hammer3D::Qsi_15P[2][2]),
	  4. * Hammer3D::Qsi_15P[2][0] * (1. - Hammer3D::Qsi_15P[2][0] - Hammer3D::Qsi_15P[2][1] - Hammer3D::Qsi_15P[2][2]),
	  (2. * Hammer3D::Qsi_15P[2][0] - 1.) * Hammer3D::Qsi_15P[2][0],
	  4. * Hammer3D::Qsi_15P[2][1] * (1. - Hammer3D::Qsi_15P[2][0] - Hammer3D::Qsi_15P[2][1] - Hammer3D::Qsi_15P[2][2]),
	  4. * Hammer3D::Qsi_15P[2][0] * Hammer3D::Qsi_15P[2][1],
	  (2. * Hammer3D::Qsi_15P[2][1] - 1.) * Hammer3D::Qsi_15P[2][1],
	  4. * Hammer3D::Qsi_15P[2][2] * (1. - Hammer3D::Qsi_15P[2][0] - Hammer3D::Qsi_15P[2][1] - Hammer3D::Qsi_15P[2][2]),
	  4. * Hammer3D::Qsi_15P[2][0] * Hammer3D::Qsi_15P[2][2],
	  4. * Hammer3D::Qsi_15P[2][1] * Hammer3D::Qsi_15P[2][2],
	  (2. * Hammer3D::Qsi_15P[2][2] - 1.) * Hammer3D::Qsi_15P[2][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[3][0] - 2. * Hammer3D::Qsi_15P[3][1] - 2. * Hammer3D::Qsi_15P[3][2]) * (1. - Hammer3D::Qsi_15P[3][0] - Hammer3D::Qsi_15P[3][1] - Hammer3D::Qsi_15P[3][2]),
	  4. * Hammer3D::Qsi_15P[3][0] * (1. - Hammer3D::Qsi_15P[3][0] - Hammer3D::Qsi_15P[3][1] - Hammer3D::Qsi_15P[3][2]),
	  (2. * Hammer3D::Qsi_15P[3][0] - 1.) * Hammer3D::Qsi_15P[3][0],
	  4. * Hammer3D::Qsi_15P[3][1] * (1. - Hammer3D::Qsi_15P[3][0] - Hammer3D::Qsi_15P[3][1] - Hammer3D::Qsi_15P[3][2]),
	  4. * Hammer3D::Qsi_15P[3][0] * Hammer3D::Qsi_15P[3][1],
	  (2. * Hammer3D::Qsi_15P[3][1] - 1.) * Hammer3D::Qsi_15P[3][1],
	  4. * Hammer3D::Qsi_15P[3][2] * (1. - Hammer3D::Qsi_15P[3][0] - Hammer3D::Qsi_15P[3][1] - Hammer3D::Qsi_15P[3][2]),
	  4. * Hammer3D::Qsi_15P[3][0] * Hammer3D::Qsi_15P[3][2],
	  4. * Hammer3D::Qsi_15P[3][1] * Hammer3D::Qsi_15P[3][2],
	  (2. * Hammer3D::Qsi_15P[3][2] - 1.) * Hammer3D::Qsi_15P[3][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[4][0] - 2. * Hammer3D::Qsi_15P[4][1] - 2. * Hammer3D::Qsi_15P[4][2]) * (1. - Hammer3D::Qsi_15P[4][0] - Hammer3D::Qsi_15P[4][1] - Hammer3D::Qsi_15P[4][2]),
	  4. * Hammer3D::Qsi_15P[4][0] * (1. - Hammer3D::Qsi_15P[4][0] - Hammer3D::Qsi_15P[4][1] - Hammer3D::Qsi_15P[4][2]),
	  (2. * Hammer3D::Qsi_15P[4][0] - 1.) * Hammer3D::Qsi_15P[4][0],
	  4. * Hammer3D::Qsi_15P[4][1] * (1. - Hammer3D::Qsi_15P[4][0] - Hammer3D::Qsi_15P[4][1] - Hammer3D::Qsi_15P[4][2]),
	  4. * Hammer3D::Qsi_15P[4][0] * Hammer3D::Qsi_15P[4][1],
	  (2. * Hammer3D::Qsi_15P[4][1] - 1.) * Hammer3D::Qsi_15P[4][1],
	  4. * Hammer3D::Qsi_15P[4][2] * (1. - Hammer3D::Qsi_15P[4][0] - Hammer3D::Qsi_15P[4][1] - Hammer3D::Qsi_15P[4][2]),
	  4. * Hammer3D::Qsi_15P[4][0] * Hammer3D::Qsi_15P[4][2],
	  4. * Hammer3D::Qsi_15P[4][1] * Hammer3D::Qsi_15P[4][2],
	  (2. * Hammer3D::Qsi_15P[4][2] - 1.) * Hammer3D::Qsi_15P[4][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[5][0] - 2. * Hammer3D::Qsi_15P[5][1] - 2. * Hammer3D::Qsi_15P[5][2]) * (1. - Hammer3D::Qsi_15P[5][0] - Hammer3D::Qsi_15P[5][1] - Hammer3D::Qsi_15P[5][2]),
	  4. * Hammer3D::Qsi_15P[5][0] * (1. - Hammer3D::Qsi_15P[5][0] - Hammer3D::Qsi_15P[5][1] - Hammer3D::Qsi_15P[5][2]),
	  (2. * Hammer3D::Qsi_15P[5][0] - 1.) * Hammer3D::Qsi_15P[5][0],
	  4. * Hammer3D::Qsi_15P[5][1] * (1. - Hammer3D::Qsi_15P[5][0] - Hammer3D::Qsi_15P[5][1] - Hammer3D::Qsi_15P[5][2]),
	  4. * Hammer3D::Qsi_15P[5][0] * Hammer3D::Qsi_15P[5][1],
	  (2. * Hammer3D::Qsi_15P[5][1] - 1.) * Hammer3D::Qsi_15P[5][1],
	  4. * Hammer3D::Qsi_15P[5][2] * (1. - Hammer3D::Qsi_15P[5][0] - Hammer3D::Qsi_15P[5][1] - Hammer3D::Qsi_15P[5][2]),
	  4. * Hammer3D::Qsi_15P[5][0] * Hammer3D::Qsi_15P[5][2],
	  4. * Hammer3D::Qsi_15P[5][1] * Hammer3D::Qsi_15P[5][2],
	  (2. * Hammer3D::Qsi_15P[5][2] - 1.) * Hammer3D::Qsi_15P[5][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[6][0] - 2. * Hammer3D::Qsi_15P[6][1] - 2. * Hammer3D::Qsi_15P[6][2]) * (1. - Hammer3D::Qsi_15P[6][0] - Hammer3D::Qsi_15P[6][1] - Hammer3D::Qsi_15P[6][2]),
	  4. * Hammer3D::Qsi_15P[6][0] * (1. - Hammer3D::Qsi_15P[6][0] - Hammer3D::Qsi_15P[6][1] - Hammer3D::Qsi_15P[6][2]),
	  (2. * Hammer3D::Qsi_15P[6][0] - 1.) * Hammer3D::Qsi_15P[6][0],
	  4. * Hammer3D::Qsi_15P[6][1] * (1. - Hammer3D::Qsi_15P[6][0] - Hammer3D::Qsi_15P[6][1] - Hammer3D::Qsi_15P[6][2]),
	  4. * Hammer3D::Qsi_15P[6][0] * Hammer3D::Qsi_15P[6][1],
	  (2. * Hammer3D::Qsi_15P[6][1] - 1.) * Hammer3D::Qsi_15P[6][1],
	  4. * Hammer3D::Qsi_15P[6][2] * (1. - Hammer3D::Qsi_15P[6][0] - Hammer3D::Qsi_15P[6][1] - Hammer3D::Qsi_15P[6][2]),
	  4. * Hammer3D::Qsi_15P[6][0] * Hammer3D::Qsi_15P[6][2],
	  4. * Hammer3D::Qsi_15P[6][1] * Hammer3D::Qsi_15P[6][2],
	  (2. * Hammer3D::Qsi_15P[6][2] - 1.) * Hammer3D::Qsi_15P[6][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[7][0] - 2. * Hammer3D::Qsi_15P[7][1] - 2. * Hammer3D::Qsi_15P[7][2]) * (1. - Hammer3D::Qsi_15P[7][0] - Hammer3D::Qsi_15P[7][1] - Hammer3D::Qsi_15P[7][2]),
	  4. * Hammer3D::Qsi_15P[7][0] * (1. - Hammer3D::Qsi_15P[7][0] - Hammer3D::Qsi_15P[7][1] - Hammer3D::Qsi_15P[7][2]),
	  (2. * Hammer3D::Qsi_15P[7][0] - 1.) * Hammer3D::Qsi_15P[7][0],
	  4. * Hammer3D::Qsi_15P[7][1] * (1. - Hammer3D::Qsi_15P[7][0] - Hammer3D::Qsi_15P[7][1] - Hammer3D::Qsi_15P[7][2]),
	  4. * Hammer3D::Qsi_15P[7][0] * Hammer3D::Qsi_15P[7][1],
	  (2. * Hammer3D::Qsi_15P[7][1] - 1.) * Hammer3D::Qsi_15P[7][1],
	  4. * Hammer3D::Qsi_15P[7][2] * (1. - Hammer3D::Qsi_15P[7][0] - Hammer3D::Qsi_15P[7][1] - Hammer3D::Qsi_15P[7][2]),
	  4. * Hammer3D::Qsi_15P[7][0] * Hammer3D::Qsi_15P[7][2],
	  4. * Hammer3D::Qsi_15P[7][1] * Hammer3D::Qsi_15P[7][2],
	  (2. * Hammer3D::Qsi_15P[7][2] - 1.) * Hammer3D::Qsi_15P[7][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[8][0] - 2. * Hammer3D::Qsi_15P[8][1] - 2. * Hammer3D::Qsi_15P[8][2]) * (1. - Hammer3D::Qsi_15P[8][0] - Hammer3D::Qsi_15P[8][1] - Hammer3D::Qsi_15P[8][2]),
	  4. * Hammer3D::Qsi_15P[8][0] * (1. - Hammer3D::Qsi_15P[8][0] - Hammer3D::Qsi_15P[8][1] - Hammer3D::Qsi_15P[8][2]),
	  (2. * Hammer3D::Qsi_15P[8][0] - 1.) * Hammer3D::Qsi_15P[8][0],
	  4. * Hammer3D::Qsi_15P[8][1] * (1. - Hammer3D::Qsi_15P[8][0] - Hammer3D::Qsi_15P[8][1] - Hammer3D::Qsi_15P[8][2]),
	  4. * Hammer3D::Qsi_15P[8][0] * Hammer3D::Qsi_15P[8][1],
	  (2. * Hammer3D::Qsi_15P[8][1] - 1.) * Hammer3D::Qsi_15P[8][1],
	  4. * Hammer3D::Qsi_15P[8][2] * (1. - Hammer3D::Qsi_15P[8][0] - Hammer3D::Qsi_15P[8][1] - Hammer3D::Qsi_15P[8][2]),
	  4. * Hammer3D::Qsi_15P[8][0] * Hammer3D::Qsi_15P[8][2],
	  4. * Hammer3D::Qsi_15P[8][1] * Hammer3D::Qsi_15P[8][2],
	  (2. * Hammer3D::Qsi_15P[8][2] - 1.) * Hammer3D::Qsi_15P[8][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[9][0] - 2. * Hammer3D::Qsi_15P[9][1] - 2. * Hammer3D::Qsi_15P[9][2]) * (1. - Hammer3D::Qsi_15P[9][0] - Hammer3D::Qsi_15P[9][1] - Hammer3D::Qsi_15P[9][2]),
	  4. * Hammer3D::Qsi_15P[9][0] * (1. - Hammer3D::Qsi_15P[9][0] - Hammer3D::Qsi_15P[9][1] - Hammer3D::Qsi_15P[9][2]),
	  (2. * Hammer3D::Qsi_15P[9][0] - 1.) * Hammer3D::Qsi_15P[9][0],
	  4. * Hammer3D::Qsi_15P[9][1] * (1. - Hammer3D::Qsi_15P[9][0] - Hammer3D::Qsi_15P[9][1] - Hammer3D::Qsi_15P[9][2]),
	  4. * Hammer3D::Qsi_15P[9][0] * Hammer3D::Qsi_15P[9][1],
	  (2. * Hammer3D::Qsi_15P[9][1] - 1.) * Hammer3D::Qsi_15P[9][1],
	  4. * Hammer3D::Qsi_15P[9][2] * (1. - Hammer3D::Qsi_15P[9][0] - Hammer3D::Qsi_15P[9][1] - Hammer3D::Qsi_15P[9][2]),
	  4. * Hammer3D::Qsi_15P[9][0] * Hammer3D::Qsi_15P[9][2],
	  4. * Hammer3D::Qsi_15P[9][1] * Hammer3D::Qsi_15P[9][2],
	  (2. * Hammer3D::Qsi_15P[9][2] - 1.) * Hammer3D::Qsi_15P[9][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[10][0] - 2. * Hammer3D::Qsi_15P[10][1] - 2. * Hammer3D::Qsi_15P[10][2]) * (1. - Hammer3D::Qsi_15P[10][0] - Hammer3D::Qsi_15P[10][1] - Hammer3D::Qsi_15P[10][2]),
		4. * Hammer3D::Qsi_15P[10][0] * (1. - Hammer3D::Qsi_15P[10][0] - Hammer3D::Qsi_15P[10][1] - Hammer3D::Qsi_15P[10][2]),
		(2. * Hammer3D::Qsi_15P[10][0] - 1.) * Hammer3D::Qsi_15P[10][0],
		4. * Hammer3D::Qsi_15P[10][1] * (1. - Hammer3D::Qsi_15P[10][0] - Hammer3D::Qsi_15P[10][1] - Hammer3D::Qsi_15P[10][2]),
		4. * Hammer3D::Qsi_15P[10][0] * Hammer3D::Qsi_15P[10][1],
		(2. * Hammer3D::Qsi_15P[10][1] - 1.) * Hammer3D::Qsi_15P[10][1],
		4. * Hammer3D::Qsi_15P[10][2] * (1. - Hammer3D::Qsi_15P[10][0] - Hammer3D::Qsi_15P[10][1] - Hammer3D::Qsi_15P[10][2]),
		4. * Hammer3D::Qsi_15P[10][0] * Hammer3D::Qsi_15P[10][2],
		4. * Hammer3D::Qsi_15P[10][1] * Hammer3D::Qsi_15P[10][2],
		(2. * Hammer3D::Qsi_15P[10][2] - 1.) * Hammer3D::Qsi_15P[10][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[11][0] - 2. * Hammer3D::Qsi_15P[11][1] - 2. * Hammer3D::Qsi_15P[11][2]) * (1. - Hammer3D::Qsi_15P[11][0] - Hammer3D::Qsi_15P[11][1] - Hammer3D::Qsi_15P[11][2]),
		4. * Hammer3D::Qsi_15P[11][0] * (1. - Hammer3D::Qsi_15P[11][0] - Hammer3D::Qsi_15P[11][1] - Hammer3D::Qsi_15P[11][2]),
		(2. * Hammer3D::Qsi_15P[11][0] - 1.) * Hammer3D::Qsi_15P[11][0],
		4. * Hammer3D::Qsi_15P[11][1] * (1. - Hammer3D::Qsi_15P[11][0] - Hammer3D::Qsi_15P[11][1] - Hammer3D::Qsi_15P[11][2]),
		4. * Hammer3D::Qsi_15P[11][0] * Hammer3D::Qsi_15P[11][1],
		(2. * Hammer3D::Qsi_15P[11][1] - 1.) * Hammer3D::Qsi_15P[11][1],
		4. * Hammer3D::Qsi_15P[11][2] * (1. - Hammer3D::Qsi_15P[11][0] - Hammer3D::Qsi_15P[11][1] - Hammer3D::Qsi_15P[11][2]),
		4. * Hammer3D::Qsi_15P[11][0] * Hammer3D::Qsi_15P[11][2],
		4. * Hammer3D::Qsi_15P[11][1] * Hammer3D::Qsi_15P[11][2],
		(2. * Hammer3D::Qsi_15P[11][2] - 1.) * Hammer3D::Qsi_15P[11][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[12][0] - 2. * Hammer3D::Qsi_15P[12][1] - 2. * Hammer3D::Qsi_15P[12][2]) * (1. - Hammer3D::Qsi_15P[12][0] - Hammer3D::Qsi_15P[12][1] - Hammer3D::Qsi_15P[12][2]),
		4. * Hammer3D::Qsi_15P[12][0] * (1. - Hammer3D::Qsi_15P[12][0] - Hammer3D::Qsi_15P[12][1] - Hammer3D::Qsi_15P[12][2]),
		(2. * Hammer3D::Qsi_15P[12][0] - 1.) * Hammer3D::Qsi_15P[12][0],
		4. * Hammer3D::Qsi_15P[12][1] * (1. - Hammer3D::Qsi_15P[12][0] - Hammer3D::Qsi_15P[12][1] - Hammer3D::Qsi_15P[12][2]),
		4. * Hammer3D::Qsi_15P[12][0] * Hammer3D::Qsi_15P[12][1],
		(2. * Hammer3D::Qsi_15P[12][1] - 1.) * Hammer3D::Qsi_15P[12][1],
		4. * Hammer3D::Qsi_15P[12][2] * (1. - Hammer3D::Qsi_15P[12][0] - Hammer3D::Qsi_15P[12][1] - Hammer3D::Qsi_15P[12][2]),
		4. * Hammer3D::Qsi_15P[12][0] * Hammer3D::Qsi_15P[12][2],
		4. * Hammer3D::Qsi_15P[12][1] * Hammer3D::Qsi_15P[12][2],
		(2. * Hammer3D::Qsi_15P[12][2] - 1.) * Hammer3D::Qsi_15P[12][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[13][0] - 2. * Hammer3D::Qsi_15P[13][1] - 2. * Hammer3D::Qsi_15P[13][2]) * (1. - Hammer3D::Qsi_15P[13][0] - Hammer3D::Qsi_15P[13][1] - Hammer3D::Qsi_15P[13][2]),
		4. * Hammer3D::Qsi_15P[13][0] * (1. - Hammer3D::Qsi_15P[13][0] - Hammer3D::Qsi_15P[13][1] - Hammer3D::Qsi_15P[13][2]),
		(2. * Hammer3D::Qsi_15P[13][0] - 1.) * Hammer3D::Qsi_15P[13][0],
		4. * Hammer3D::Qsi_15P[13][1] * (1. - Hammer3D::Qsi_15P[13][0] - Hammer3D::Qsi_15P[13][1] - Hammer3D::Qsi_15P[13][2]),
		4. * Hammer3D::Qsi_15P[13][0] * Hammer3D::Qsi_15P[13][1],
		(2. * Hammer3D::Qsi_15P[13][1] - 1.) * Hammer3D::Qsi_15P[13][1],
		4. * Hammer3D::Qsi_15P[13][2] * (1. - Hammer3D::Qsi_15P[13][0] - Hammer3D::Qsi_15P[13][1] - Hammer3D::Qsi_15P[13][2]),
		4. * Hammer3D::Qsi_15P[13][0] * Hammer3D::Qsi_15P[13][2],
		4. * Hammer3D::Qsi_15P[13][1] * Hammer3D::Qsi_15P[13][2],
		(2. * Hammer3D::Qsi_15P[13][2] - 1.) * Hammer3D::Qsi_15P[13][2] },

	{ (1. - 2. * Hammer3D::Qsi_15P[14][0] - 2. * Hammer3D::Qsi_15P[14][1] - 2. * Hammer3D::Qsi_15P[14][2]) * (1. - Hammer3D::Qsi_15P[14][0] - Hammer3D::Qsi_15P[14][1] - Hammer3D::Qsi_15P[14][2]),
		4. * Hammer3D::Qsi_15P[14][0] * (1. - Hammer3D::Qsi_15P[14][0] - Hammer3D::Qsi_15P[14][1] - Hammer3D::Qsi_15P[14][2]),
		(2. * Hammer3D::Qsi_15P[14][0] - 1.) * Hammer3D::Qsi_15P[14][0],
		4. * Hammer3D::Qsi_15P[14][1] * (1. - Hammer3D::Qsi_15P[14][0] - Hammer3D::Qsi_15P[14][1] - Hammer3D::Qsi_15P[14][2]),
		4. * Hammer3D::Qsi_15P[14][0] * Hammer3D::Qsi_15P[14][1],
		(2. * Hammer3D::Qsi_15P[14][1] - 1.) * Hammer3D::Qsi_15P[14][1],
		4. * Hammer3D::Qsi_15P[14][2] * (1. - Hammer3D::Qsi_15P[14][0] - Hammer3D::Qsi_15P[14][1] - Hammer3D::Qsi_15P[14][2]),
		4. * Hammer3D::Qsi_15P[14][0] * Hammer3D::Qsi_15P[14][2],
		4. * Hammer3D::Qsi_15P[14][1] * Hammer3D::Qsi_15P[14][2],
		(2. * Hammer3D::Qsi_15P[14][2] - 1.) * Hammer3D::Qsi_15P[14][2] } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<24>::mv_Psi[24][mv_numNodes] = {
	{ (1. - 2. * Hammer3D::Qsi_24P[0][0] - 2. * Hammer3D::Qsi_24P[0][1] - 2. * Hammer3D::Qsi_24P[0][2]) * (1. - Hammer3D::Qsi_24P[0][0] - Hammer3D::Qsi_24P[0][1] - Hammer3D::Qsi_24P[0][2]),
	  4. * Hammer3D::Qsi_24P[0][0] * (1. - Hammer3D::Qsi_24P[0][0] - Hammer3D::Qsi_24P[0][1] - Hammer3D::Qsi_24P[0][2]),
	  (2. * Hammer3D::Qsi_24P[0][0] - 1.) * Hammer3D::Qsi_24P[0][0],
	  4. * Hammer3D::Qsi_24P[0][1] * (1. - Hammer3D::Qsi_24P[0][0] - Hammer3D::Qsi_24P[0][1] - Hammer3D::Qsi_24P[0][2]),
	  4. * Hammer3D::Qsi_24P[0][0] * Hammer3D::Qsi_24P[0][1],
	  (2. * Hammer3D::Qsi_24P[0][1] - 1.) * Hammer3D::Qsi_24P[0][1],
	  4. * Hammer3D::Qsi_24P[0][2] * (1. - Hammer3D::Qsi_24P[0][0] - Hammer3D::Qsi_24P[0][1] - Hammer3D::Qsi_24P[0][2]),
	  4. * Hammer3D::Qsi_24P[0][0] * Hammer3D::Qsi_24P[0][2],
	  4. * Hammer3D::Qsi_24P[0][1] * Hammer3D::Qsi_24P[0][2],
	  (2. * Hammer3D::Qsi_24P[0][2] - 1.) * Hammer3D::Qsi_24P[0][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[1][0] - 2. * Hammer3D::Qsi_24P[1][1] - 2. * Hammer3D::Qsi_24P[1][2]) * (1. - Hammer3D::Qsi_24P[1][0] - Hammer3D::Qsi_24P[1][1] - Hammer3D::Qsi_24P[1][2]),
	  4. * Hammer3D::Qsi_24P[1][0] * (1. - Hammer3D::Qsi_24P[1][0] - Hammer3D::Qsi_24P[1][1] - Hammer3D::Qsi_24P[1][2]),
	  (2. * Hammer3D::Qsi_24P[1][0] - 1.) * Hammer3D::Qsi_24P[1][0],
	  4. * Hammer3D::Qsi_24P[1][1] * (1. - Hammer3D::Qsi_24P[1][0] - Hammer3D::Qsi_24P[1][1] - Hammer3D::Qsi_24P[1][2]),
	  4. * Hammer3D::Qsi_24P[1][0] * Hammer3D::Qsi_24P[1][1],
	  (2. * Hammer3D::Qsi_24P[1][1] - 1.) * Hammer3D::Qsi_24P[1][1],
	  4. * Hammer3D::Qsi_24P[1][2] * (1. - Hammer3D::Qsi_24P[1][0] - Hammer3D::Qsi_24P[1][1] - Hammer3D::Qsi_24P[1][2]),
	  4. * Hammer3D::Qsi_24P[1][0] * Hammer3D::Qsi_24P[1][2],
	  4. * Hammer3D::Qsi_24P[1][1] * Hammer3D::Qsi_24P[1][2],
	  (2. * Hammer3D::Qsi_24P[1][2] - 1.) * Hammer3D::Qsi_24P[1][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[2][0] - 2. * Hammer3D::Qsi_24P[2][1] - 2. * Hammer3D::Qsi_24P[2][2]) * (1. - Hammer3D::Qsi_24P[2][0] - Hammer3D::Qsi_24P[2][1] - Hammer3D::Qsi_24P[2][2]),
	  4. * Hammer3D::Qsi_24P[2][0] * (1. - Hammer3D::Qsi_24P[2][0] - Hammer3D::Qsi_24P[2][1] - Hammer3D::Qsi_24P[2][2]),
	  (2. * Hammer3D::Qsi_24P[2][0] - 1.) * Hammer3D::Qsi_24P[2][0],
	  4. * Hammer3D::Qsi_24P[2][1] * (1. - Hammer3D::Qsi_24P[2][0] - Hammer3D::Qsi_24P[2][1] - Hammer3D::Qsi_24P[2][2]),
	  4. * Hammer3D::Qsi_24P[2][0] * Hammer3D::Qsi_24P[2][1],
	  (2. * Hammer3D::Qsi_24P[2][1] - 1.) * Hammer3D::Qsi_24P[2][1],
	  4. * Hammer3D::Qsi_24P[2][2] * (1. - Hammer3D::Qsi_24P[2][0] - Hammer3D::Qsi_24P[2][1] - Hammer3D::Qsi_24P[2][2]),
	  4. * Hammer3D::Qsi_24P[2][0] * Hammer3D::Qsi_24P[2][2],
	  4. * Hammer3D::Qsi_24P[2][1] * Hammer3D::Qsi_24P[2][2],
	  (2. * Hammer3D::Qsi_24P[2][2] - 1.) * Hammer3D::Qsi_24P[2][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[3][0] - 2. * Hammer3D::Qsi_24P[3][1] - 2. * Hammer3D::Qsi_24P[3][2]) * (1. - Hammer3D::Qsi_24P[3][0] - Hammer3D::Qsi_24P[3][1] - Hammer3D::Qsi_24P[3][2]),
	  4. * Hammer3D::Qsi_24P[3][0] * (1. - Hammer3D::Qsi_24P[3][0] - Hammer3D::Qsi_24P[3][1] - Hammer3D::Qsi_24P[3][2]),
	  (2. * Hammer3D::Qsi_24P[3][0] - 1.) * Hammer3D::Qsi_24P[3][0],
	  4. * Hammer3D::Qsi_24P[3][1] * (1. - Hammer3D::Qsi_24P[3][0] - Hammer3D::Qsi_24P[3][1] - Hammer3D::Qsi_24P[3][2]),
	  4. * Hammer3D::Qsi_24P[3][0] * Hammer3D::Qsi_24P[3][1],
	  (2. * Hammer3D::Qsi_24P[3][1] - 1.) * Hammer3D::Qsi_24P[3][1],
	  4. * Hammer3D::Qsi_24P[3][2] * (1. - Hammer3D::Qsi_24P[3][0] - Hammer3D::Qsi_24P[3][1] - Hammer3D::Qsi_24P[3][2]),
	  4. * Hammer3D::Qsi_24P[3][0] * Hammer3D::Qsi_24P[3][2],
	  4. * Hammer3D::Qsi_24P[3][1] * Hammer3D::Qsi_24P[3][2],
	  (2. * Hammer3D::Qsi_24P[3][2] - 1.) * Hammer3D::Qsi_24P[3][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[4][0] - 2. * Hammer3D::Qsi_24P[4][1] - 2. * Hammer3D::Qsi_24P[4][2]) * (1. - Hammer3D::Qsi_24P[4][0] - Hammer3D::Qsi_24P[4][1] - Hammer3D::Qsi_24P[4][2]),
	  4. * Hammer3D::Qsi_24P[4][0] * (1. - Hammer3D::Qsi_24P[4][0] - Hammer3D::Qsi_24P[4][1] - Hammer3D::Qsi_24P[4][2]),
	  (2. * Hammer3D::Qsi_24P[4][0] - 1.) * Hammer3D::Qsi_24P[4][0],
	  4. * Hammer3D::Qsi_24P[4][1] * (1. - Hammer3D::Qsi_24P[4][0] - Hammer3D::Qsi_24P[4][1] - Hammer3D::Qsi_24P[4][2]),
	  4. * Hammer3D::Qsi_24P[4][0] * Hammer3D::Qsi_24P[4][1],
	  (2. * Hammer3D::Qsi_24P[4][1] - 1.) * Hammer3D::Qsi_24P[4][1],
	  4. * Hammer3D::Qsi_24P[4][2] * (1. - Hammer3D::Qsi_24P[4][0] - Hammer3D::Qsi_24P[4][1] - Hammer3D::Qsi_24P[4][2]),
	  4. * Hammer3D::Qsi_24P[4][0] * Hammer3D::Qsi_24P[4][2],
	  4. * Hammer3D::Qsi_24P[4][1] * Hammer3D::Qsi_24P[4][2],
	  (2. * Hammer3D::Qsi_24P[4][2] - 1.) * Hammer3D::Qsi_24P[4][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[5][0] - 2. * Hammer3D::Qsi_24P[5][1] - 2. * Hammer3D::Qsi_24P[5][2]) * (1. - Hammer3D::Qsi_24P[5][0] - Hammer3D::Qsi_24P[5][1] - Hammer3D::Qsi_24P[5][2]),
	  4. * Hammer3D::Qsi_24P[5][0] * (1. - Hammer3D::Qsi_24P[5][0] - Hammer3D::Qsi_24P[5][1] - Hammer3D::Qsi_24P[5][2]),
	  (2. * Hammer3D::Qsi_24P[5][0] - 1.) * Hammer3D::Qsi_24P[5][0],
	  4. * Hammer3D::Qsi_24P[5][1] * (1. - Hammer3D::Qsi_24P[5][0] - Hammer3D::Qsi_24P[5][1] - Hammer3D::Qsi_24P[5][2]),
	  4. * Hammer3D::Qsi_24P[5][0] * Hammer3D::Qsi_24P[5][1],
	  (2. * Hammer3D::Qsi_24P[5][1] - 1.) * Hammer3D::Qsi_24P[5][1],
	  4. * Hammer3D::Qsi_24P[5][2] * (1. - Hammer3D::Qsi_24P[5][0] - Hammer3D::Qsi_24P[5][1] - Hammer3D::Qsi_24P[5][2]),
	  4. * Hammer3D::Qsi_24P[5][0] * Hammer3D::Qsi_24P[5][2],
	  4. * Hammer3D::Qsi_24P[5][1] * Hammer3D::Qsi_24P[5][2],
	  (2. * Hammer3D::Qsi_24P[5][2] - 1.) * Hammer3D::Qsi_24P[5][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[6][0] - 2. * Hammer3D::Qsi_24P[6][1] - 2. * Hammer3D::Qsi_24P[6][2]) * (1. - Hammer3D::Qsi_24P[6][0] - Hammer3D::Qsi_24P[6][1] - Hammer3D::Qsi_24P[6][2]),
	  4. * Hammer3D::Qsi_24P[6][0] * (1. - Hammer3D::Qsi_24P[6][0] - Hammer3D::Qsi_24P[6][1] - Hammer3D::Qsi_24P[6][2]),
	  (2. * Hammer3D::Qsi_24P[6][0] - 1.) * Hammer3D::Qsi_24P[6][0],
	  4. * Hammer3D::Qsi_24P[6][1] * (1. - Hammer3D::Qsi_24P[6][0] - Hammer3D::Qsi_24P[6][1] - Hammer3D::Qsi_24P[6][2]),
	  4. * Hammer3D::Qsi_24P[6][0] * Hammer3D::Qsi_24P[6][1],
	  (2. * Hammer3D::Qsi_24P[6][1] - 1.) * Hammer3D::Qsi_24P[6][1],
	  4. * Hammer3D::Qsi_24P[6][2] * (1. - Hammer3D::Qsi_24P[6][0] - Hammer3D::Qsi_24P[6][1] - Hammer3D::Qsi_24P[6][2]),
	  4. * Hammer3D::Qsi_24P[6][0] * Hammer3D::Qsi_24P[6][2],
	  4. * Hammer3D::Qsi_24P[6][1] * Hammer3D::Qsi_24P[6][2],
	  (2. * Hammer3D::Qsi_24P[6][2] - 1.) * Hammer3D::Qsi_24P[6][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[7][0] - 2. * Hammer3D::Qsi_24P[7][1] - 2. * Hammer3D::Qsi_24P[7][2]) * (1. - Hammer3D::Qsi_24P[7][0] - Hammer3D::Qsi_24P[7][1] - Hammer3D::Qsi_24P[7][2]),
	  4. * Hammer3D::Qsi_24P[7][0] * (1. - Hammer3D::Qsi_24P[7][0] - Hammer3D::Qsi_24P[7][1] - Hammer3D::Qsi_24P[7][2]),
	  (2. * Hammer3D::Qsi_24P[7][0] - 1.) * Hammer3D::Qsi_24P[7][0],
	  4. * Hammer3D::Qsi_24P[7][1] * (1. - Hammer3D::Qsi_24P[7][0] - Hammer3D::Qsi_24P[7][1] - Hammer3D::Qsi_24P[7][2]),
	  4. * Hammer3D::Qsi_24P[7][0] * Hammer3D::Qsi_24P[7][1],
	  (2. * Hammer3D::Qsi_24P[7][1] - 1.) * Hammer3D::Qsi_24P[7][1],
	  4. * Hammer3D::Qsi_24P[7][2] * (1. - Hammer3D::Qsi_24P[7][0] - Hammer3D::Qsi_24P[7][1] - Hammer3D::Qsi_24P[7][2]),
	  4. * Hammer3D::Qsi_24P[7][0] * Hammer3D::Qsi_24P[7][2],
	  4. * Hammer3D::Qsi_24P[7][1] * Hammer3D::Qsi_24P[7][2],
	  (2. * Hammer3D::Qsi_24P[7][2] - 1.) * Hammer3D::Qsi_24P[7][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[8][0] - 2. * Hammer3D::Qsi_24P[8][1] - 2. * Hammer3D::Qsi_24P[8][2]) * (1. - Hammer3D::Qsi_24P[8][0] - Hammer3D::Qsi_24P[8][1] - Hammer3D::Qsi_24P[8][2]),
	  4. * Hammer3D::Qsi_24P[8][0] * (1. - Hammer3D::Qsi_24P[8][0] - Hammer3D::Qsi_24P[8][1] - Hammer3D::Qsi_24P[8][2]),
	  (2. * Hammer3D::Qsi_24P[8][0] - 1.) * Hammer3D::Qsi_24P[8][0],
	  4. * Hammer3D::Qsi_24P[8][1] * (1. - Hammer3D::Qsi_24P[8][0] - Hammer3D::Qsi_24P[8][1] - Hammer3D::Qsi_24P[8][2]),
	  4. * Hammer3D::Qsi_24P[8][0] * Hammer3D::Qsi_24P[8][1],
	  (2. * Hammer3D::Qsi_24P[8][1] - 1.) * Hammer3D::Qsi_24P[8][1],
	  4. * Hammer3D::Qsi_24P[8][2] * (1. - Hammer3D::Qsi_24P[8][0] - Hammer3D::Qsi_24P[8][1] - Hammer3D::Qsi_24P[8][2]),
	  4. * Hammer3D::Qsi_24P[8][0] * Hammer3D::Qsi_24P[8][2],
	  4. * Hammer3D::Qsi_24P[8][1] * Hammer3D::Qsi_24P[8][2],
	  (2. * Hammer3D::Qsi_24P[8][2] - 1.) * Hammer3D::Qsi_24P[8][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[9][0] - 2. * Hammer3D::Qsi_24P[9][1] - 2. * Hammer3D::Qsi_24P[9][2]) * (1. - Hammer3D::Qsi_24P[9][0] - Hammer3D::Qsi_24P[9][1] - Hammer3D::Qsi_24P[9][2]),
	  4. * Hammer3D::Qsi_24P[9][0] * (1. - Hammer3D::Qsi_24P[9][0] - Hammer3D::Qsi_24P[9][1] - Hammer3D::Qsi_24P[9][2]),
	  (2. * Hammer3D::Qsi_24P[9][0] - 1.) * Hammer3D::Qsi_24P[9][0],
	  4. * Hammer3D::Qsi_24P[9][1] * (1. - Hammer3D::Qsi_24P[9][0] - Hammer3D::Qsi_24P[9][1] - Hammer3D::Qsi_24P[9][2]),
	  4. * Hammer3D::Qsi_24P[9][0] * Hammer3D::Qsi_24P[9][1],
	  (2. * Hammer3D::Qsi_24P[9][1] - 1.) * Hammer3D::Qsi_24P[9][1],
	  4. * Hammer3D::Qsi_24P[9][2] * (1. - Hammer3D::Qsi_24P[9][0] - Hammer3D::Qsi_24P[9][1] - Hammer3D::Qsi_24P[9][2]),
	  4. * Hammer3D::Qsi_24P[9][0] * Hammer3D::Qsi_24P[9][2],
	  4. * Hammer3D::Qsi_24P[9][1] * Hammer3D::Qsi_24P[9][2],
	  (2. * Hammer3D::Qsi_24P[9][2] - 1.) * Hammer3D::Qsi_24P[9][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[10][0] - 2. * Hammer3D::Qsi_24P[10][1] - 2. * Hammer3D::Qsi_24P[10][2]) * (1. - Hammer3D::Qsi_24P[10][0] - Hammer3D::Qsi_24P[10][1] - Hammer3D::Qsi_24P[10][2]),
		4. * Hammer3D::Qsi_24P[10][0] * (1. - Hammer3D::Qsi_24P[10][0] - Hammer3D::Qsi_24P[10][1] - Hammer3D::Qsi_24P[10][2]),
		(2. * Hammer3D::Qsi_24P[10][0] - 1.) * Hammer3D::Qsi_24P[10][0],
		4. * Hammer3D::Qsi_24P[10][1] * (1. - Hammer3D::Qsi_24P[10][0] - Hammer3D::Qsi_24P[10][1] - Hammer3D::Qsi_24P[10][2]),
		4. * Hammer3D::Qsi_24P[10][0] * Hammer3D::Qsi_24P[10][1],
		(2. * Hammer3D::Qsi_24P[10][1] - 1.) * Hammer3D::Qsi_24P[10][1],
		4. * Hammer3D::Qsi_24P[10][2] * (1. - Hammer3D::Qsi_24P[10][0] - Hammer3D::Qsi_24P[10][1] - Hammer3D::Qsi_24P[10][2]),
		4. * Hammer3D::Qsi_24P[10][0] * Hammer3D::Qsi_24P[10][2],
		4. * Hammer3D::Qsi_24P[10][1] * Hammer3D::Qsi_24P[10][2],
		(2. * Hammer3D::Qsi_24P[10][2] - 1.) * Hammer3D::Qsi_24P[10][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[11][0] - 2. * Hammer3D::Qsi_24P[11][1] - 2. * Hammer3D::Qsi_24P[11][2]) * (1. - Hammer3D::Qsi_24P[11][0] - Hammer3D::Qsi_24P[11][1] - Hammer3D::Qsi_24P[11][2]),
		4. * Hammer3D::Qsi_24P[11][0] * (1. - Hammer3D::Qsi_24P[11][0] - Hammer3D::Qsi_24P[11][1] - Hammer3D::Qsi_24P[11][2]),
		(2. * Hammer3D::Qsi_24P[11][0] - 1.) * Hammer3D::Qsi_24P[11][0],
		4. * Hammer3D::Qsi_24P[11][1] * (1. - Hammer3D::Qsi_24P[11][0] - Hammer3D::Qsi_24P[11][1] - Hammer3D::Qsi_24P[11][2]),
		4. * Hammer3D::Qsi_24P[11][0] * Hammer3D::Qsi_24P[11][1],
		(2. * Hammer3D::Qsi_24P[11][1] - 1.) * Hammer3D::Qsi_24P[11][1],
		4. * Hammer3D::Qsi_24P[11][2] * (1. - Hammer3D::Qsi_24P[11][0] - Hammer3D::Qsi_24P[11][1] - Hammer3D::Qsi_24P[11][2]),
		4. * Hammer3D::Qsi_24P[11][0] * Hammer3D::Qsi_24P[11][2],
		4. * Hammer3D::Qsi_24P[11][1] * Hammer3D::Qsi_24P[11][2],
		(2. * Hammer3D::Qsi_24P[11][2] - 1.) * Hammer3D::Qsi_24P[11][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[12][0] - 2. * Hammer3D::Qsi_24P[12][1] - 2. * Hammer3D::Qsi_24P[12][2]) * (1. - Hammer3D::Qsi_24P[12][0] - Hammer3D::Qsi_24P[12][1] - Hammer3D::Qsi_24P[12][2]),
		4. * Hammer3D::Qsi_24P[12][0] * (1. - Hammer3D::Qsi_24P[12][0] - Hammer3D::Qsi_24P[12][1] - Hammer3D::Qsi_24P[12][2]),
		(2. * Hammer3D::Qsi_24P[12][0] - 1.) * Hammer3D::Qsi_24P[12][0],
		4. * Hammer3D::Qsi_24P[12][1] * (1. - Hammer3D::Qsi_24P[12][0] - Hammer3D::Qsi_24P[12][1] - Hammer3D::Qsi_24P[12][2]),
		4. * Hammer3D::Qsi_24P[12][0] * Hammer3D::Qsi_24P[12][1],
		(2. * Hammer3D::Qsi_24P[12][1] - 1.) * Hammer3D::Qsi_24P[12][1],
		4. * Hammer3D::Qsi_24P[12][2] * (1. - Hammer3D::Qsi_24P[12][0] - Hammer3D::Qsi_24P[12][1] - Hammer3D::Qsi_24P[12][2]),
		4. * Hammer3D::Qsi_24P[12][0] * Hammer3D::Qsi_24P[12][2],
		4. * Hammer3D::Qsi_24P[12][1] * Hammer3D::Qsi_24P[12][2],
		(2. * Hammer3D::Qsi_24P[12][2] - 1.) * Hammer3D::Qsi_24P[12][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[13][0] - 2. * Hammer3D::Qsi_24P[13][1] - 2. * Hammer3D::Qsi_24P[13][2]) * (1. - Hammer3D::Qsi_24P[13][0] - Hammer3D::Qsi_24P[13][1] - Hammer3D::Qsi_24P[13][2]),
		4. * Hammer3D::Qsi_24P[13][0] * (1. - Hammer3D::Qsi_24P[13][0] - Hammer3D::Qsi_24P[13][1] - Hammer3D::Qsi_24P[13][2]),
		(2. * Hammer3D::Qsi_24P[13][0] - 1.) * Hammer3D::Qsi_24P[13][0],
		4. * Hammer3D::Qsi_24P[13][1] * (1. - Hammer3D::Qsi_24P[13][0] - Hammer3D::Qsi_24P[13][1] - Hammer3D::Qsi_24P[13][2]),
		4. * Hammer3D::Qsi_24P[13][0] * Hammer3D::Qsi_24P[13][1],
		(2. * Hammer3D::Qsi_24P[13][1] - 1.) * Hammer3D::Qsi_24P[13][1],
		4. * Hammer3D::Qsi_24P[13][2] * (1. - Hammer3D::Qsi_24P[13][0] - Hammer3D::Qsi_24P[13][1] - Hammer3D::Qsi_24P[13][2]),
		4. * Hammer3D::Qsi_24P[13][0] * Hammer3D::Qsi_24P[13][2],
		4. * Hammer3D::Qsi_24P[13][1] * Hammer3D::Qsi_24P[13][2],
		(2. * Hammer3D::Qsi_24P[13][2] - 1.) * Hammer3D::Qsi_24P[13][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[14][0] - 2. * Hammer3D::Qsi_24P[14][1] - 2. * Hammer3D::Qsi_24P[14][2]) * (1. - Hammer3D::Qsi_24P[14][0] - Hammer3D::Qsi_24P[14][1] - Hammer3D::Qsi_24P[14][2]),
		4. * Hammer3D::Qsi_24P[14][0] * (1. - Hammer3D::Qsi_24P[14][0] - Hammer3D::Qsi_24P[14][1] - Hammer3D::Qsi_24P[14][2]),
		(2. * Hammer3D::Qsi_24P[14][0] - 1.) * Hammer3D::Qsi_24P[14][0],
		4. * Hammer3D::Qsi_24P[14][1] * (1. - Hammer3D::Qsi_24P[14][0] - Hammer3D::Qsi_24P[14][1] - Hammer3D::Qsi_24P[14][2]),
		4. * Hammer3D::Qsi_24P[14][0] * Hammer3D::Qsi_24P[14][1],
		(2. * Hammer3D::Qsi_24P[14][1] - 1.) * Hammer3D::Qsi_24P[14][1],
		4. * Hammer3D::Qsi_24P[14][2] * (1. - Hammer3D::Qsi_24P[14][0] - Hammer3D::Qsi_24P[14][1] - Hammer3D::Qsi_24P[14][2]),
		4. * Hammer3D::Qsi_24P[14][0] * Hammer3D::Qsi_24P[14][2],
		4. * Hammer3D::Qsi_24P[14][1] * Hammer3D::Qsi_24P[14][2],
		(2. * Hammer3D::Qsi_24P[14][2] - 1.) * Hammer3D::Qsi_24P[14][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[15][0] - 2. * Hammer3D::Qsi_24P[15][1] - 2. * Hammer3D::Qsi_24P[15][2]) * (1. - Hammer3D::Qsi_24P[15][0] - Hammer3D::Qsi_24P[15][1] - Hammer3D::Qsi_24P[15][2]),
		4. * Hammer3D::Qsi_24P[15][0] * (1. - Hammer3D::Qsi_24P[15][0] - Hammer3D::Qsi_24P[15][1] - Hammer3D::Qsi_24P[15][2]),
		(2. * Hammer3D::Qsi_24P[15][0] - 1.) * Hammer3D::Qsi_24P[15][0],
		4. * Hammer3D::Qsi_24P[15][1] * (1. - Hammer3D::Qsi_24P[15][0] - Hammer3D::Qsi_24P[15][1] - Hammer3D::Qsi_24P[15][2]),
		4. * Hammer3D::Qsi_24P[15][0] * Hammer3D::Qsi_24P[15][1],
		(2. * Hammer3D::Qsi_24P[15][1] - 1.) * Hammer3D::Qsi_24P[15][1],
		4. * Hammer3D::Qsi_24P[15][2] * (1. - Hammer3D::Qsi_24P[15][0] - Hammer3D::Qsi_24P[15][1] - Hammer3D::Qsi_24P[15][2]),
		4. * Hammer3D::Qsi_24P[15][0] * Hammer3D::Qsi_24P[15][2],
		4. * Hammer3D::Qsi_24P[15][1] * Hammer3D::Qsi_24P[15][2],
		(2. * Hammer3D::Qsi_24P[15][2] - 1.) * Hammer3D::Qsi_24P[15][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[16][0] - 2. * Hammer3D::Qsi_24P[16][1] - 2. * Hammer3D::Qsi_24P[16][2]) * (1. - Hammer3D::Qsi_24P[16][0] - Hammer3D::Qsi_24P[16][1] - Hammer3D::Qsi_24P[16][2]),
		4. * Hammer3D::Qsi_24P[16][0] * (1. - Hammer3D::Qsi_24P[16][0] - Hammer3D::Qsi_24P[16][1] - Hammer3D::Qsi_24P[16][2]),
		(2. * Hammer3D::Qsi_24P[16][0] - 1.) * Hammer3D::Qsi_24P[16][0],
		4. * Hammer3D::Qsi_24P[16][1] * (1. - Hammer3D::Qsi_24P[16][0] - Hammer3D::Qsi_24P[16][1] - Hammer3D::Qsi_24P[16][2]),
		4. * Hammer3D::Qsi_24P[16][0] * Hammer3D::Qsi_24P[16][1],
		(2. * Hammer3D::Qsi_24P[16][1] - 1.) * Hammer3D::Qsi_24P[16][1],
		4. * Hammer3D::Qsi_24P[16][2] * (1. - Hammer3D::Qsi_24P[16][0] - Hammer3D::Qsi_24P[16][1] - Hammer3D::Qsi_24P[16][2]),
		4. * Hammer3D::Qsi_24P[16][0] * Hammer3D::Qsi_24P[16][2],
		4. * Hammer3D::Qsi_24P[16][1] * Hammer3D::Qsi_24P[16][2],
		(2. * Hammer3D::Qsi_24P[16][2] - 1.) * Hammer3D::Qsi_24P[16][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[17][0] - 2. * Hammer3D::Qsi_24P[17][1] - 2. * Hammer3D::Qsi_24P[17][2]) * (1. - Hammer3D::Qsi_24P[17][0] - Hammer3D::Qsi_24P[17][1] - Hammer3D::Qsi_24P[17][2]),
		4. * Hammer3D::Qsi_24P[17][0] * (1. - Hammer3D::Qsi_24P[17][0] - Hammer3D::Qsi_24P[17][1] - Hammer3D::Qsi_24P[17][2]),
		(2. * Hammer3D::Qsi_24P[17][0] - 1.) * Hammer3D::Qsi_24P[17][0],
		4. * Hammer3D::Qsi_24P[17][1] * (1. - Hammer3D::Qsi_24P[17][0] - Hammer3D::Qsi_24P[17][1] - Hammer3D::Qsi_24P[17][2]),
		4. * Hammer3D::Qsi_24P[17][0] * Hammer3D::Qsi_24P[17][1],
		(2. * Hammer3D::Qsi_24P[17][1] - 1.) * Hammer3D::Qsi_24P[17][1],
		4. * Hammer3D::Qsi_24P[17][2] * (1. - Hammer3D::Qsi_24P[17][0] - Hammer3D::Qsi_24P[17][1] - Hammer3D::Qsi_24P[17][2]),
		4. * Hammer3D::Qsi_24P[17][0] * Hammer3D::Qsi_24P[17][2],
		4. * Hammer3D::Qsi_24P[17][1] * Hammer3D::Qsi_24P[17][2],
		(2. * Hammer3D::Qsi_24P[17][2] - 1.) * Hammer3D::Qsi_24P[17][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[18][0] - 2. * Hammer3D::Qsi_24P[18][1] - 2. * Hammer3D::Qsi_24P[18][2]) * (1. - Hammer3D::Qsi_24P[18][0] - Hammer3D::Qsi_24P[18][1] - Hammer3D::Qsi_24P[18][2]),
		4. * Hammer3D::Qsi_24P[18][0] * (1. - Hammer3D::Qsi_24P[18][0] - Hammer3D::Qsi_24P[18][1] - Hammer3D::Qsi_24P[18][2]),
		(2. * Hammer3D::Qsi_24P[18][0] - 1.) * Hammer3D::Qsi_24P[18][0],
		4. * Hammer3D::Qsi_24P[18][1] * (1. - Hammer3D::Qsi_24P[18][0] - Hammer3D::Qsi_24P[18][1] - Hammer3D::Qsi_24P[18][2]),
		4. * Hammer3D::Qsi_24P[18][0] * Hammer3D::Qsi_24P[18][1],
		(2. * Hammer3D::Qsi_24P[18][1] - 1.) * Hammer3D::Qsi_24P[18][1],
		4. * Hammer3D::Qsi_24P[18][2] * (1. - Hammer3D::Qsi_24P[18][0] - Hammer3D::Qsi_24P[18][1] - Hammer3D::Qsi_24P[18][2]),
		4. * Hammer3D::Qsi_24P[18][0] * Hammer3D::Qsi_24P[18][2],
		4. * Hammer3D::Qsi_24P[18][1] * Hammer3D::Qsi_24P[18][2],
		(2. * Hammer3D::Qsi_24P[18][2] - 1.) * Hammer3D::Qsi_24P[18][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[19][0] - 2. * Hammer3D::Qsi_24P[19][1] - 2. * Hammer3D::Qsi_24P[19][2]) * (1. - Hammer3D::Qsi_24P[19][0] - Hammer3D::Qsi_24P[19][1] - Hammer3D::Qsi_24P[19][2]),
		4. * Hammer3D::Qsi_24P[19][0] * (1. - Hammer3D::Qsi_24P[19][0] - Hammer3D::Qsi_24P[19][1] - Hammer3D::Qsi_24P[19][2]),
		(2. * Hammer3D::Qsi_24P[19][0] - 1.) * Hammer3D::Qsi_24P[19][0],
		4. * Hammer3D::Qsi_24P[19][1] * (1. - Hammer3D::Qsi_24P[19][0] - Hammer3D::Qsi_24P[19][1] - Hammer3D::Qsi_24P[19][2]),
		4. * Hammer3D::Qsi_24P[19][0] * Hammer3D::Qsi_24P[19][1],
		(2. * Hammer3D::Qsi_24P[19][1] - 1.) * Hammer3D::Qsi_24P[19][1],
		4. * Hammer3D::Qsi_24P[19][2] * (1. - Hammer3D::Qsi_24P[19][0] - Hammer3D::Qsi_24P[19][1] - Hammer3D::Qsi_24P[19][2]),
		4. * Hammer3D::Qsi_24P[19][0] * Hammer3D::Qsi_24P[19][2],
		4. * Hammer3D::Qsi_24P[19][1] * Hammer3D::Qsi_24P[19][2],
		(2. * Hammer3D::Qsi_24P[19][2] - 1.) * Hammer3D::Qsi_24P[19][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[20][0] - 2. * Hammer3D::Qsi_24P[20][1] - 2. * Hammer3D::Qsi_24P[20][2]) * (1. - Hammer3D::Qsi_24P[20][0] - Hammer3D::Qsi_24P[20][1] - Hammer3D::Qsi_24P[20][2]),
		4. * Hammer3D::Qsi_24P[20][0] * (1. - Hammer3D::Qsi_24P[20][0] - Hammer3D::Qsi_24P[20][1] - Hammer3D::Qsi_24P[20][2]),
		(2. * Hammer3D::Qsi_24P[20][0] - 1.) * Hammer3D::Qsi_24P[20][0],
		4. * Hammer3D::Qsi_24P[20][1] * (1. - Hammer3D::Qsi_24P[20][0] - Hammer3D::Qsi_24P[20][1] - Hammer3D::Qsi_24P[20][2]),
		4. * Hammer3D::Qsi_24P[20][0] * Hammer3D::Qsi_24P[20][1],
		(2. * Hammer3D::Qsi_24P[20][1] - 1.) * Hammer3D::Qsi_24P[20][1],
		4. * Hammer3D::Qsi_24P[20][2] * (1. - Hammer3D::Qsi_24P[20][0] - Hammer3D::Qsi_24P[20][1] - Hammer3D::Qsi_24P[20][2]),
		4. * Hammer3D::Qsi_24P[20][0] * Hammer3D::Qsi_24P[20][2],
		4. * Hammer3D::Qsi_24P[20][1] * Hammer3D::Qsi_24P[20][2],
		(2. * Hammer3D::Qsi_24P[20][2] - 1.) * Hammer3D::Qsi_24P[20][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[21][0] - 2. * Hammer3D::Qsi_24P[21][1] - 2. * Hammer3D::Qsi_24P[21][2]) * (1. - Hammer3D::Qsi_24P[21][0] - Hammer3D::Qsi_24P[21][1] - Hammer3D::Qsi_24P[21][2]),
		4. * Hammer3D::Qsi_24P[21][0] * (1. - Hammer3D::Qsi_24P[21][0] - Hammer3D::Qsi_24P[21][1] - Hammer3D::Qsi_24P[21][2]),
		(2. * Hammer3D::Qsi_24P[21][0] - 1.) * Hammer3D::Qsi_24P[21][0],
		4. * Hammer3D::Qsi_24P[21][1] * (1. - Hammer3D::Qsi_24P[21][0] - Hammer3D::Qsi_24P[21][1] - Hammer3D::Qsi_24P[21][2]),
		4. * Hammer3D::Qsi_24P[21][0] * Hammer3D::Qsi_24P[21][1],
		(2. * Hammer3D::Qsi_24P[21][1] - 1.) * Hammer3D::Qsi_24P[21][1],
		4. * Hammer3D::Qsi_24P[21][2] * (1. - Hammer3D::Qsi_24P[21][0] - Hammer3D::Qsi_24P[21][1] - Hammer3D::Qsi_24P[21][2]),
		4. * Hammer3D::Qsi_24P[21][0] * Hammer3D::Qsi_24P[21][2],
		4. * Hammer3D::Qsi_24P[21][1] * Hammer3D::Qsi_24P[21][2],
		(2. * Hammer3D::Qsi_24P[21][2] - 1.) * Hammer3D::Qsi_24P[21][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[22][0] - 2. * Hammer3D::Qsi_24P[22][1] - 2. * Hammer3D::Qsi_24P[22][2]) * (1. - Hammer3D::Qsi_24P[22][0] - Hammer3D::Qsi_24P[22][1] - Hammer3D::Qsi_24P[22][2]),
		4. * Hammer3D::Qsi_24P[22][0] * (1. - Hammer3D::Qsi_24P[22][0] - Hammer3D::Qsi_24P[22][1] - Hammer3D::Qsi_24P[22][2]),
		(2. * Hammer3D::Qsi_24P[22][0] - 1.) * Hammer3D::Qsi_24P[22][0],
		4. * Hammer3D::Qsi_24P[22][1] * (1. - Hammer3D::Qsi_24P[22][0] - Hammer3D::Qsi_24P[22][1] - Hammer3D::Qsi_24P[22][2]),
		4. * Hammer3D::Qsi_24P[22][0] * Hammer3D::Qsi_24P[22][1],
		(2. * Hammer3D::Qsi_24P[22][1] - 1.) * Hammer3D::Qsi_24P[22][1],
		4. * Hammer3D::Qsi_24P[22][2] * (1. - Hammer3D::Qsi_24P[22][0] - Hammer3D::Qsi_24P[22][1] - Hammer3D::Qsi_24P[22][2]),
		4. * Hammer3D::Qsi_24P[22][0] * Hammer3D::Qsi_24P[22][2],
		4. * Hammer3D::Qsi_24P[22][1] * Hammer3D::Qsi_24P[22][2],
		(2. * Hammer3D::Qsi_24P[22][2] - 1.) * Hammer3D::Qsi_24P[22][2] },

	{ (1. - 2. * Hammer3D::Qsi_24P[23][0] - 2. * Hammer3D::Qsi_24P[23][1] - 2. * Hammer3D::Qsi_24P[23][2]) * (1. - Hammer3D::Qsi_24P[23][0] - Hammer3D::Qsi_24P[23][1] - Hammer3D::Qsi_24P[23][2]),
		4. * Hammer3D::Qsi_24P[23][0] * (1. - Hammer3D::Qsi_24P[23][0] - Hammer3D::Qsi_24P[23][1] - Hammer3D::Qsi_24P[23][2]),
		(2. * Hammer3D::Qsi_24P[23][0] - 1.) * Hammer3D::Qsi_24P[23][0],
		4. * Hammer3D::Qsi_24P[23][1] * (1. - Hammer3D::Qsi_24P[23][0] - Hammer3D::Qsi_24P[23][1] - Hammer3D::Qsi_24P[23][2]),
		4. * Hammer3D::Qsi_24P[23][0] * Hammer3D::Qsi_24P[23][1],
		(2. * Hammer3D::Qsi_24P[23][1] - 1.) * Hammer3D::Qsi_24P[23][1],
		4. * Hammer3D::Qsi_24P[23][2] * (1. - Hammer3D::Qsi_24P[23][0] - Hammer3D::Qsi_24P[23][1] - Hammer3D::Qsi_24P[23][2]),
		4. * Hammer3D::Qsi_24P[23][0] * Hammer3D::Qsi_24P[23][2],
		4. * Hammer3D::Qsi_24P[23][1] * Hammer3D::Qsi_24P[23][2],
		(2. * Hammer3D::Qsi_24P[23][2] - 1.) * Hammer3D::Qsi_24P[23][2] } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<4>::mv_DPsi[4][mv_numNodes][mv_ElDim] = {
	{ { -3. + 4. * (Hammer3D::Qsi_4P[0][0] + Hammer3D::Qsi_4P[0][1] + Hammer3D::Qsi_4P[0][2]), -3. + 4. * (Hammer3D::Qsi_4P[0][0] + Hammer3D::Qsi_4P[0][1] + Hammer3D::Qsi_4P[0][2]), -3. + 4. * (Hammer3D::Qsi_4P[0][0] + Hammer3D::Qsi_4P[0][1] + Hammer3D::Qsi_4P[0][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_4P[0][0] + Hammer3D::Qsi_4P[0][1] + Hammer3D::Qsi_4P[0][2]), -4. * Hammer3D::Qsi_4P[0][0], -4. * Hammer3D::Qsi_4P[0][0] },
	  { 4. * Hammer3D::Qsi_4P[0][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_4P[0][1], 4. - 4. * (Hammer3D::Qsi_4P[0][0] + 2. * Hammer3D::Qsi_4P[0][1] + Hammer3D::Qsi_4P[0][2]), -4. * Hammer3D::Qsi_4P[0][1] },
	  { 4. * Hammer3D::Qsi_4P[0][1], 4. * Hammer3D::Qsi_4P[0][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_4P[0][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_4P[0][2], -4. * Hammer3D::Qsi_4P[0][2], 4. - 4. * (Hammer3D::Qsi_4P[0][0] + Hammer3D::Qsi_4P[0][1] + 2. * Hammer3D::Qsi_4P[0][2]) },
	  { 4. * Hammer3D::Qsi_4P[0][2], 0., 4. * Hammer3D::Qsi_4P[0][0] },
	  { 0., 4. * Hammer3D::Qsi_4P[0][2], 4. * Hammer3D::Qsi_4P[0][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_4P[0][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_4P[1][0] + Hammer3D::Qsi_4P[1][1] + Hammer3D::Qsi_4P[1][2]), -3. + 4. * (Hammer3D::Qsi_4P[1][0] + Hammer3D::Qsi_4P[1][1] + Hammer3D::Qsi_4P[1][2]), -3. + 4. * (Hammer3D::Qsi_4P[1][0] + Hammer3D::Qsi_4P[1][1] + Hammer3D::Qsi_4P[1][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_4P[1][0] + Hammer3D::Qsi_4P[1][1] + Hammer3D::Qsi_4P[1][2]), -4. * Hammer3D::Qsi_4P[1][0], -4. * Hammer3D::Qsi_4P[1][0] },
	  { 4. * Hammer3D::Qsi_4P[1][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_4P[1][1], 4. - 4. * (Hammer3D::Qsi_4P[1][0] + 2. * Hammer3D::Qsi_4P[1][1] + Hammer3D::Qsi_4P[1][2]), -4. * Hammer3D::Qsi_4P[1][1] },
	  { 4. * Hammer3D::Qsi_4P[1][1], 4. * Hammer3D::Qsi_4P[1][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_4P[1][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_4P[1][2], -4. * Hammer3D::Qsi_4P[1][2], 4. - 4. * (Hammer3D::Qsi_4P[1][0] + Hammer3D::Qsi_4P[1][1] + 2. * Hammer3D::Qsi_4P[1][2]) },
	  { 4. * Hammer3D::Qsi_4P[1][2], 0., 4. * Hammer3D::Qsi_4P[1][0] },
	  { 0., 4. * Hammer3D::Qsi_4P[1][2], 4. * Hammer3D::Qsi_4P[1][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_4P[1][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_4P[2][0] + Hammer3D::Qsi_4P[2][1] + Hammer3D::Qsi_4P[2][2]), -3. + 4. * (Hammer3D::Qsi_4P[2][0] + Hammer3D::Qsi_4P[2][1] + Hammer3D::Qsi_4P[2][2]), -3. + 4. * (Hammer3D::Qsi_4P[2][0] + Hammer3D::Qsi_4P[2][1] + Hammer3D::Qsi_4P[2][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_4P[2][0] + Hammer3D::Qsi_4P[2][1] + Hammer3D::Qsi_4P[2][2]), -4. * Hammer3D::Qsi_4P[2][0], -4. * Hammer3D::Qsi_4P[2][0] },
	  { 4. * Hammer3D::Qsi_4P[2][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_4P[2][1], 4. - 4. * (Hammer3D::Qsi_4P[2][0] + 2. * Hammer3D::Qsi_4P[2][1] + Hammer3D::Qsi_4P[2][2]), -4. * Hammer3D::Qsi_4P[2][1] },
	  { 4. * Hammer3D::Qsi_4P[2][1], 4. * Hammer3D::Qsi_4P[2][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_4P[2][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_4P[2][2], -4. * Hammer3D::Qsi_4P[2][2], 4. - 4. * (Hammer3D::Qsi_4P[2][0] + Hammer3D::Qsi_4P[2][1] + 2. * Hammer3D::Qsi_4P[2][2]) },
	  { 4. * Hammer3D::Qsi_4P[2][2], 0., 4. * Hammer3D::Qsi_4P[2][0] },
	  { 0., 4. * Hammer3D::Qsi_4P[2][2], 4. * Hammer3D::Qsi_4P[2][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_4P[2][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_4P[3][0] + Hammer3D::Qsi_4P[3][1] + Hammer3D::Qsi_4P[3][2]), -3. + 4. * (Hammer3D::Qsi_4P[3][0] + Hammer3D::Qsi_4P[3][1] + Hammer3D::Qsi_4P[3][2]), -3. + 4. * (Hammer3D::Qsi_4P[3][0] + Hammer3D::Qsi_4P[3][1] + Hammer3D::Qsi_4P[3][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_4P[3][0] + Hammer3D::Qsi_4P[3][1] + Hammer3D::Qsi_4P[3][2]), -4. * Hammer3D::Qsi_4P[3][0], -4. * Hammer3D::Qsi_4P[3][0] },
	  { 4. * Hammer3D::Qsi_4P[3][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_4P[3][1], 4. - 4. * (Hammer3D::Qsi_4P[3][0] + 2. * Hammer3D::Qsi_4P[3][1] + Hammer3D::Qsi_4P[3][2]), -4. * Hammer3D::Qsi_4P[3][1] },
	  { 4. * Hammer3D::Qsi_4P[3][1], 4. * Hammer3D::Qsi_4P[3][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_4P[3][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_4P[3][2], -4. * Hammer3D::Qsi_4P[3][2], 4. - 4. * (Hammer3D::Qsi_4P[3][0] + Hammer3D::Qsi_4P[3][1] + 2. * Hammer3D::Qsi_4P[3][2]) },
	  { 4. * Hammer3D::Qsi_4P[3][2], 0., 4. * Hammer3D::Qsi_4P[3][0] },
	  { 0., 4. * Hammer3D::Qsi_4P[3][2], 4. * Hammer3D::Qsi_4P[3][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_4P[3][2] - 1.} } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<10>::mv_DPsi[10][mv_numNodes][mv_ElDim] = {
	{ { -3. + 4. * (Hammer3D::Qsi_10P[0][0] + Hammer3D::Qsi_10P[0][1] + Hammer3D::Qsi_10P[0][2]), -3. + 4. * (Hammer3D::Qsi_10P[0][0] + Hammer3D::Qsi_10P[0][1] + Hammer3D::Qsi_10P[0][2]), -3. + 4. * (Hammer3D::Qsi_10P[0][0] + Hammer3D::Qsi_10P[0][1] + Hammer3D::Qsi_10P[0][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[0][0] + Hammer3D::Qsi_10P[0][1] + Hammer3D::Qsi_10P[0][2]), -4. * Hammer3D::Qsi_10P[0][0], -4. * Hammer3D::Qsi_10P[0][0] },
	  { 4. * Hammer3D::Qsi_10P[0][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[0][1], 4. - 4. * (Hammer3D::Qsi_10P[0][0] + 2. * Hammer3D::Qsi_10P[0][1] + Hammer3D::Qsi_10P[0][2]), -4. * Hammer3D::Qsi_10P[0][1] },
	  { 4. * Hammer3D::Qsi_10P[0][1], 4. * Hammer3D::Qsi_10P[0][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[0][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[0][2], -4. * Hammer3D::Qsi_10P[0][2], 4. - 4. * (Hammer3D::Qsi_10P[0][0] + Hammer3D::Qsi_10P[0][1] + 2. * Hammer3D::Qsi_10P[0][2]) },
	  { 4. * Hammer3D::Qsi_10P[0][2], 0., 4. * Hammer3D::Qsi_10P[0][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[0][2], 4. * Hammer3D::Qsi_10P[0][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[0][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[1][0] + Hammer3D::Qsi_10P[1][1] + Hammer3D::Qsi_10P[1][2]), -3. + 4. * (Hammer3D::Qsi_10P[1][0] + Hammer3D::Qsi_10P[1][1] + Hammer3D::Qsi_10P[1][2]), -3. + 4. * (Hammer3D::Qsi_10P[1][0] + Hammer3D::Qsi_10P[1][1] + Hammer3D::Qsi_10P[1][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[1][0] + Hammer3D::Qsi_10P[1][1] + Hammer3D::Qsi_10P[1][2]), -4. * Hammer3D::Qsi_10P[1][0], -4. * Hammer3D::Qsi_10P[1][0] },
	  { 4. * Hammer3D::Qsi_10P[1][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[1][1], 4. - 4. * (Hammer3D::Qsi_10P[1][0] + 2. * Hammer3D::Qsi_10P[1][1] + Hammer3D::Qsi_10P[1][2]), -4. * Hammer3D::Qsi_10P[1][1] },
	  { 4. * Hammer3D::Qsi_10P[1][1], 4. * Hammer3D::Qsi_10P[1][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[1][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[1][2], -4. * Hammer3D::Qsi_10P[1][2], 4. - 4. * (Hammer3D::Qsi_10P[1][0] + Hammer3D::Qsi_10P[1][1] + 2. * Hammer3D::Qsi_10P[1][2]) },
	  { 4. * Hammer3D::Qsi_10P[1][2], 0., 4. * Hammer3D::Qsi_10P[1][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[1][2], 4. * Hammer3D::Qsi_10P[1][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[1][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[2][0] + Hammer3D::Qsi_10P[2][1] + Hammer3D::Qsi_10P[2][2]), -3. + 4. * (Hammer3D::Qsi_10P[2][0] + Hammer3D::Qsi_10P[2][1] + Hammer3D::Qsi_10P[2][2]), -3. + 4. * (Hammer3D::Qsi_10P[2][0] + Hammer3D::Qsi_10P[2][1] + Hammer3D::Qsi_10P[2][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[2][0] + Hammer3D::Qsi_10P[2][1] + Hammer3D::Qsi_10P[2][2]), -4. * Hammer3D::Qsi_10P[2][0], -4. * Hammer3D::Qsi_10P[2][0] },
	  { 4. * Hammer3D::Qsi_10P[2][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[2][1], 4. - 4. * (Hammer3D::Qsi_10P[2][0] + 2. * Hammer3D::Qsi_10P[2][1] + Hammer3D::Qsi_10P[2][2]), -4. * Hammer3D::Qsi_10P[2][1] },
	  { 4. * Hammer3D::Qsi_10P[2][1], 4. * Hammer3D::Qsi_10P[2][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[2][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[2][2], -4. * Hammer3D::Qsi_10P[2][2], 4. - 4. * (Hammer3D::Qsi_10P[2][0] + Hammer3D::Qsi_10P[2][1] + 2. * Hammer3D::Qsi_10P[2][2]) },
	  { 4. * Hammer3D::Qsi_10P[2][2], 0., 4. * Hammer3D::Qsi_10P[2][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[2][2], 4. * Hammer3D::Qsi_10P[2][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[2][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[3][0] + Hammer3D::Qsi_10P[3][1] + Hammer3D::Qsi_10P[3][2]), -3. + 4. * (Hammer3D::Qsi_10P[3][0] + Hammer3D::Qsi_10P[3][1] + Hammer3D::Qsi_10P[3][2]), -3. + 4. * (Hammer3D::Qsi_10P[3][0] + Hammer3D::Qsi_10P[3][1] + Hammer3D::Qsi_10P[3][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[3][0] + Hammer3D::Qsi_10P[3][1] + Hammer3D::Qsi_10P[3][2]), -4. * Hammer3D::Qsi_10P[3][0], -4. * Hammer3D::Qsi_10P[3][0] },
	  { 4. * Hammer3D::Qsi_10P[3][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[3][1], 4. - 4. * (Hammer3D::Qsi_10P[3][0] + 2. * Hammer3D::Qsi_10P[3][1] + Hammer3D::Qsi_10P[3][2]), -4. * Hammer3D::Qsi_10P[3][1] },
	  { 4. * Hammer3D::Qsi_10P[3][1], 4. * Hammer3D::Qsi_10P[3][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[3][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[3][2], -4. * Hammer3D::Qsi_10P[3][2], 4. - 4. * (Hammer3D::Qsi_10P[3][0] + Hammer3D::Qsi_10P[3][1] + 2. * Hammer3D::Qsi_10P[3][2]) },
	  { 4. * Hammer3D::Qsi_10P[3][2], 0., 4. * Hammer3D::Qsi_10P[3][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[3][2], 4. * Hammer3D::Qsi_10P[3][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[3][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[4][0] + Hammer3D::Qsi_10P[4][1] + Hammer3D::Qsi_10P[4][2]), -3. + 4. * (Hammer3D::Qsi_10P[4][0] + Hammer3D::Qsi_10P[4][1] + Hammer3D::Qsi_10P[4][2]), -3. + 4. * (Hammer3D::Qsi_10P[4][0] + Hammer3D::Qsi_10P[4][1] + Hammer3D::Qsi_10P[4][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[4][0] + Hammer3D::Qsi_10P[4][1] + Hammer3D::Qsi_10P[4][2]), -4. * Hammer3D::Qsi_10P[4][0], -4. * Hammer3D::Qsi_10P[4][0] },
	  { 4. * Hammer3D::Qsi_10P[4][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[4][1], 4. - 4. * (Hammer3D::Qsi_10P[4][0] + 2. * Hammer3D::Qsi_10P[4][1] + Hammer3D::Qsi_10P[4][2]), -4. * Hammer3D::Qsi_10P[4][1] },
	  { 4. * Hammer3D::Qsi_10P[4][1], 4. * Hammer3D::Qsi_10P[4][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[4][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[4][2], -4. * Hammer3D::Qsi_10P[4][2], 4. - 4. * (Hammer3D::Qsi_10P[4][0] + Hammer3D::Qsi_10P[4][1] + 2. * Hammer3D::Qsi_10P[4][2]) },
	  { 4. * Hammer3D::Qsi_10P[4][2], 0., 4. * Hammer3D::Qsi_10P[4][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[4][2], 4. * Hammer3D::Qsi_10P[4][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[4][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[5][0] + Hammer3D::Qsi_10P[5][1] + Hammer3D::Qsi_10P[5][2]), -3. + 4. * (Hammer3D::Qsi_10P[5][0] + Hammer3D::Qsi_10P[5][1] + Hammer3D::Qsi_10P[5][2]), -3. + 4. * (Hammer3D::Qsi_10P[5][0] + Hammer3D::Qsi_10P[5][1] + Hammer3D::Qsi_10P[5][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[5][0] + Hammer3D::Qsi_10P[5][1] + Hammer3D::Qsi_10P[5][2]), -4. * Hammer3D::Qsi_10P[5][0], -4. * Hammer3D::Qsi_10P[5][0] },
	  { 4. * Hammer3D::Qsi_10P[5][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[5][1], 4. - 4. * (Hammer3D::Qsi_10P[5][0] + 2. * Hammer3D::Qsi_10P[5][1] + Hammer3D::Qsi_10P[5][2]), -4. * Hammer3D::Qsi_10P[5][1] },
	  { 4. * Hammer3D::Qsi_10P[5][1], 4. * Hammer3D::Qsi_10P[5][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[5][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[5][2], -4. * Hammer3D::Qsi_10P[5][2], 4. - 4. * (Hammer3D::Qsi_10P[5][0] + Hammer3D::Qsi_10P[5][1] + 2. * Hammer3D::Qsi_10P[5][2]) },
	  { 4. * Hammer3D::Qsi_10P[5][2], 0., 4. * Hammer3D::Qsi_10P[5][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[5][2], 4. * Hammer3D::Qsi_10P[5][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[5][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[6][0] + Hammer3D::Qsi_10P[6][1] + Hammer3D::Qsi_10P[6][2]), -3. + 4. * (Hammer3D::Qsi_10P[6][0] + Hammer3D::Qsi_10P[6][1] + Hammer3D::Qsi_10P[6][2]), -3. + 4. * (Hammer3D::Qsi_10P[6][0] + Hammer3D::Qsi_10P[6][1] + Hammer3D::Qsi_10P[6][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[6][0] + Hammer3D::Qsi_10P[6][1] + Hammer3D::Qsi_10P[6][2]), -4. * Hammer3D::Qsi_10P[6][0], -4. * Hammer3D::Qsi_10P[6][0] },
	  { 4. * Hammer3D::Qsi_10P[6][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[6][1], 4. - 4. * (Hammer3D::Qsi_10P[6][0] + 2. * Hammer3D::Qsi_10P[6][1] + Hammer3D::Qsi_10P[6][2]), -4. * Hammer3D::Qsi_10P[6][1] },
	  { 4. * Hammer3D::Qsi_10P[6][1], 4. * Hammer3D::Qsi_10P[6][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[6][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[6][2], -4. * Hammer3D::Qsi_10P[6][2], 4. - 4. * (Hammer3D::Qsi_10P[6][0] + Hammer3D::Qsi_10P[6][1] + 2. * Hammer3D::Qsi_10P[6][2]) },
	  { 4. * Hammer3D::Qsi_10P[6][2], 0., 4. * Hammer3D::Qsi_10P[6][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[6][2], 4. * Hammer3D::Qsi_10P[6][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[6][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[7][0] + Hammer3D::Qsi_10P[7][1] + Hammer3D::Qsi_10P[7][2]), -3. + 4. * (Hammer3D::Qsi_10P[7][0] + Hammer3D::Qsi_10P[7][1] + Hammer3D::Qsi_10P[7][2]), -3. + 4. * (Hammer3D::Qsi_10P[7][0] + Hammer3D::Qsi_10P[7][1] + Hammer3D::Qsi_10P[7][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[7][0] + Hammer3D::Qsi_10P[7][1] + Hammer3D::Qsi_10P[7][2]), -4. * Hammer3D::Qsi_10P[7][0], -4. * Hammer3D::Qsi_10P[7][0] },
	  { 4. * Hammer3D::Qsi_10P[7][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[7][1], 4. - 4. * (Hammer3D::Qsi_10P[7][0] + 2. * Hammer3D::Qsi_10P[7][1] + Hammer3D::Qsi_10P[7][2]), -4. * Hammer3D::Qsi_10P[7][1] },
	  { 4. * Hammer3D::Qsi_10P[7][1], 4. * Hammer3D::Qsi_10P[7][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[7][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[7][2], -4. * Hammer3D::Qsi_10P[7][2], 4. - 4. * (Hammer3D::Qsi_10P[7][0] + Hammer3D::Qsi_10P[7][1] + 2. * Hammer3D::Qsi_10P[7][2]) },
	  { 4. * Hammer3D::Qsi_10P[7][2], 0., 4. * Hammer3D::Qsi_10P[7][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[7][2], 4. * Hammer3D::Qsi_10P[7][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[7][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[8][0] + Hammer3D::Qsi_10P[8][1] + Hammer3D::Qsi_10P[8][2]), -3. + 4. * (Hammer3D::Qsi_10P[8][0] + Hammer3D::Qsi_10P[8][1] + Hammer3D::Qsi_10P[8][2]), -3. + 4. * (Hammer3D::Qsi_10P[8][0] + Hammer3D::Qsi_10P[8][1] + Hammer3D::Qsi_10P[8][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[8][0] + Hammer3D::Qsi_10P[8][1] + Hammer3D::Qsi_10P[8][2]), -4. * Hammer3D::Qsi_10P[8][0], -4. * Hammer3D::Qsi_10P[8][0] },
	  { 4. * Hammer3D::Qsi_10P[8][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[8][1], 4. - 4. * (Hammer3D::Qsi_10P[8][0] + 2. * Hammer3D::Qsi_10P[8][1] + Hammer3D::Qsi_10P[8][2]), -4. * Hammer3D::Qsi_10P[8][1] },
	  { 4. * Hammer3D::Qsi_10P[8][1], 4. * Hammer3D::Qsi_10P[8][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[8][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[8][2], -4. * Hammer3D::Qsi_10P[8][2], 4. - 4. * (Hammer3D::Qsi_10P[8][0] + Hammer3D::Qsi_10P[8][1] + 2. * Hammer3D::Qsi_10P[8][2]) },
	  { 4. * Hammer3D::Qsi_10P[8][2], 0., 4. * Hammer3D::Qsi_10P[8][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[8][2], 4. * Hammer3D::Qsi_10P[8][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[8][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_10P[9][0] + Hammer3D::Qsi_10P[9][1] + Hammer3D::Qsi_10P[9][2]), -3. + 4. * (Hammer3D::Qsi_10P[9][0] + Hammer3D::Qsi_10P[9][1] + Hammer3D::Qsi_10P[9][2]), -3. + 4. * (Hammer3D::Qsi_10P[9][0] + Hammer3D::Qsi_10P[9][1] + Hammer3D::Qsi_10P[9][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_10P[9][0] + Hammer3D::Qsi_10P[9][1] + Hammer3D::Qsi_10P[9][2]), -4. * Hammer3D::Qsi_10P[9][0], -4. * Hammer3D::Qsi_10P[9][0] },
	  { 4. * Hammer3D::Qsi_10P[9][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_10P[9][1], 4. - 4. * (Hammer3D::Qsi_10P[9][0] + 2. * Hammer3D::Qsi_10P[9][1] + Hammer3D::Qsi_10P[9][2]), -4. * Hammer3D::Qsi_10P[9][1] },
	  { 4. * Hammer3D::Qsi_10P[9][1], 4. * Hammer3D::Qsi_10P[9][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_10P[9][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_10P[9][2], -4. * Hammer3D::Qsi_10P[9][2], 4. - 4. * (Hammer3D::Qsi_10P[9][0] + Hammer3D::Qsi_10P[9][1] + 2. * Hammer3D::Qsi_10P[9][2]) },
	  { 4. * Hammer3D::Qsi_10P[9][2], 0., 4. * Hammer3D::Qsi_10P[9][0] },
	  { 0., 4. * Hammer3D::Qsi_10P[9][2], 4. * Hammer3D::Qsi_10P[9][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_10P[9][2] - 1.} } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<11>::mv_DPsi[11][mv_numNodes][mv_ElDim] = {
	{ { -3. + 4. * (Hammer3D::Qsi_11P[0][0] + Hammer3D::Qsi_11P[0][1] + Hammer3D::Qsi_11P[0][2]), -3. + 4. * (Hammer3D::Qsi_11P[0][0] + Hammer3D::Qsi_11P[0][1] + Hammer3D::Qsi_11P[0][2]), -3. + 4. * (Hammer3D::Qsi_11P[0][0] + Hammer3D::Qsi_11P[0][1] + Hammer3D::Qsi_11P[0][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[0][0] + Hammer3D::Qsi_11P[0][1] + Hammer3D::Qsi_11P[0][2]), -4. * Hammer3D::Qsi_11P[0][0], -4. * Hammer3D::Qsi_11P[0][0] },
	  { 4. * Hammer3D::Qsi_11P[0][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[0][1], 4. - 4. * (Hammer3D::Qsi_11P[0][0] + 2. * Hammer3D::Qsi_11P[0][1] + Hammer3D::Qsi_11P[0][2]), -4. * Hammer3D::Qsi_11P[0][1] },
	  { 4. * Hammer3D::Qsi_11P[0][1], 4. * Hammer3D::Qsi_11P[0][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[0][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[0][2], -4. * Hammer3D::Qsi_11P[0][2], 4. - 4. * (Hammer3D::Qsi_11P[0][0] + Hammer3D::Qsi_11P[0][1] + 2. * Hammer3D::Qsi_11P[0][2]) },
	  { 4. * Hammer3D::Qsi_11P[0][2], 0., 4. * Hammer3D::Qsi_11P[0][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[0][2], 4. * Hammer3D::Qsi_11P[0][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[0][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[1][0] + Hammer3D::Qsi_11P[1][1] + Hammer3D::Qsi_11P[1][2]), -3. + 4. * (Hammer3D::Qsi_11P[1][0] + Hammer3D::Qsi_11P[1][1] + Hammer3D::Qsi_11P[1][2]), -3. + 4. * (Hammer3D::Qsi_11P[1][0] + Hammer3D::Qsi_11P[1][1] + Hammer3D::Qsi_11P[1][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[1][0] + Hammer3D::Qsi_11P[1][1] + Hammer3D::Qsi_11P[1][2]), -4. * Hammer3D::Qsi_11P[1][0], -4. * Hammer3D::Qsi_11P[1][0] },
	  { 4. * Hammer3D::Qsi_11P[1][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[1][1], 4. - 4. * (Hammer3D::Qsi_11P[1][0] + 2. * Hammer3D::Qsi_11P[1][1] + Hammer3D::Qsi_11P[1][2]), -4. * Hammer3D::Qsi_11P[1][1] },
	  { 4. * Hammer3D::Qsi_11P[1][1], 4. * Hammer3D::Qsi_11P[1][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[1][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[1][2], -4. * Hammer3D::Qsi_11P[1][2], 4. - 4. * (Hammer3D::Qsi_11P[1][0] + Hammer3D::Qsi_11P[1][1] + 2. * Hammer3D::Qsi_11P[1][2]) },
	  { 4. * Hammer3D::Qsi_11P[1][2], 0., 4. * Hammer3D::Qsi_11P[1][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[1][2], 4. * Hammer3D::Qsi_11P[1][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[1][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[2][0] + Hammer3D::Qsi_11P[2][1] + Hammer3D::Qsi_11P[2][2]), -3. + 4. * (Hammer3D::Qsi_11P[2][0] + Hammer3D::Qsi_11P[2][1] + Hammer3D::Qsi_11P[2][2]), -3. + 4. * (Hammer3D::Qsi_11P[2][0] + Hammer3D::Qsi_11P[2][1] + Hammer3D::Qsi_11P[2][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[2][0] + Hammer3D::Qsi_11P[2][1] + Hammer3D::Qsi_11P[2][2]), -4. * Hammer3D::Qsi_11P[2][0], -4. * Hammer3D::Qsi_11P[2][0] },
	  { 4. * Hammer3D::Qsi_11P[2][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[2][1], 4. - 4. * (Hammer3D::Qsi_11P[2][0] + 2. * Hammer3D::Qsi_11P[2][1] + Hammer3D::Qsi_11P[2][2]), -4. * Hammer3D::Qsi_11P[2][1] },
	  { 4. * Hammer3D::Qsi_11P[2][1], 4. * Hammer3D::Qsi_11P[2][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[2][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[2][2], -4. * Hammer3D::Qsi_11P[2][2], 4. - 4. * (Hammer3D::Qsi_11P[2][0] + Hammer3D::Qsi_11P[2][1] + 2. * Hammer3D::Qsi_11P[2][2]) },
	  { 4. * Hammer3D::Qsi_11P[2][2], 0., 4. * Hammer3D::Qsi_11P[2][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[2][2], 4. * Hammer3D::Qsi_11P[2][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[2][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[3][0] + Hammer3D::Qsi_11P[3][1] + Hammer3D::Qsi_11P[3][2]), -3. + 4. * (Hammer3D::Qsi_11P[3][0] + Hammer3D::Qsi_11P[3][1] + Hammer3D::Qsi_11P[3][2]), -3. + 4. * (Hammer3D::Qsi_11P[3][0] + Hammer3D::Qsi_11P[3][1] + Hammer3D::Qsi_11P[3][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[3][0] + Hammer3D::Qsi_11P[3][1] + Hammer3D::Qsi_11P[3][2]), -4. * Hammer3D::Qsi_11P[3][0], -4. * Hammer3D::Qsi_11P[3][0] },
	  { 4. * Hammer3D::Qsi_11P[3][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[3][1], 4. - 4. * (Hammer3D::Qsi_11P[3][0] + 2. * Hammer3D::Qsi_11P[3][1] + Hammer3D::Qsi_11P[3][2]), -4. * Hammer3D::Qsi_11P[3][1] },
	  { 4. * Hammer3D::Qsi_11P[3][1], 4. * Hammer3D::Qsi_11P[3][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[3][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[3][2], -4. * Hammer3D::Qsi_11P[3][2], 4. - 4. * (Hammer3D::Qsi_11P[3][0] + Hammer3D::Qsi_11P[3][1] + 2. * Hammer3D::Qsi_11P[3][2]) },
	  { 4. * Hammer3D::Qsi_11P[3][2], 0., 4. * Hammer3D::Qsi_11P[3][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[3][2], 4. * Hammer3D::Qsi_11P[3][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[3][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[4][0] + Hammer3D::Qsi_11P[4][1] + Hammer3D::Qsi_11P[4][2]), -3. + 4. * (Hammer3D::Qsi_11P[4][0] + Hammer3D::Qsi_11P[4][1] + Hammer3D::Qsi_11P[4][2]), -3. + 4. * (Hammer3D::Qsi_11P[4][0] + Hammer3D::Qsi_11P[4][1] + Hammer3D::Qsi_11P[4][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[4][0] + Hammer3D::Qsi_11P[4][1] + Hammer3D::Qsi_11P[4][2]), -4. * Hammer3D::Qsi_11P[4][0], -4. * Hammer3D::Qsi_11P[4][0] },
	  { 4. * Hammer3D::Qsi_11P[4][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[4][1], 4. - 4. * (Hammer3D::Qsi_11P[4][0] + 2. * Hammer3D::Qsi_11P[4][1] + Hammer3D::Qsi_11P[4][2]), -4. * Hammer3D::Qsi_11P[4][1] },
	  { 4. * Hammer3D::Qsi_11P[4][1], 4. * Hammer3D::Qsi_11P[4][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[4][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[4][2], -4. * Hammer3D::Qsi_11P[4][2], 4. - 4. * (Hammer3D::Qsi_11P[4][0] + Hammer3D::Qsi_11P[4][1] + 2. * Hammer3D::Qsi_11P[4][2]) },
	  { 4. * Hammer3D::Qsi_11P[4][2], 0., 4. * Hammer3D::Qsi_11P[4][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[4][2], 4. * Hammer3D::Qsi_11P[4][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[4][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[5][0] + Hammer3D::Qsi_11P[5][1] + Hammer3D::Qsi_11P[5][2]), -3. + 4. * (Hammer3D::Qsi_11P[5][0] + Hammer3D::Qsi_11P[5][1] + Hammer3D::Qsi_11P[5][2]), -3. + 4. * (Hammer3D::Qsi_11P[5][0] + Hammer3D::Qsi_11P[5][1] + Hammer3D::Qsi_11P[5][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[5][0] + Hammer3D::Qsi_11P[5][1] + Hammer3D::Qsi_11P[5][2]), -4. * Hammer3D::Qsi_11P[5][0], -4. * Hammer3D::Qsi_11P[5][0] },
	  { 4. * Hammer3D::Qsi_11P[5][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[5][1], 4. - 4. * (Hammer3D::Qsi_11P[5][0] + 2. * Hammer3D::Qsi_11P[5][1] + Hammer3D::Qsi_11P[5][2]), -4. * Hammer3D::Qsi_11P[5][1] },
	  { 4. * Hammer3D::Qsi_11P[5][1], 4. * Hammer3D::Qsi_11P[5][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[5][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[5][2], -4. * Hammer3D::Qsi_11P[5][2], 4. - 4. * (Hammer3D::Qsi_11P[5][0] + Hammer3D::Qsi_11P[5][1] + 2. * Hammer3D::Qsi_11P[5][2]) },
	  { 4. * Hammer3D::Qsi_11P[5][2], 0., 4. * Hammer3D::Qsi_11P[5][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[5][2], 4. * Hammer3D::Qsi_11P[5][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[5][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[6][0] + Hammer3D::Qsi_11P[6][1] + Hammer3D::Qsi_11P[6][2]), -3. + 4. * (Hammer3D::Qsi_11P[6][0] + Hammer3D::Qsi_11P[6][1] + Hammer3D::Qsi_11P[6][2]), -3. + 4. * (Hammer3D::Qsi_11P[6][0] + Hammer3D::Qsi_11P[6][1] + Hammer3D::Qsi_11P[6][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[6][0] + Hammer3D::Qsi_11P[6][1] + Hammer3D::Qsi_11P[6][2]), -4. * Hammer3D::Qsi_11P[6][0], -4. * Hammer3D::Qsi_11P[6][0] },
	  { 4. * Hammer3D::Qsi_11P[6][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[6][1], 4. - 4. * (Hammer3D::Qsi_11P[6][0] + 2. * Hammer3D::Qsi_11P[6][1] + Hammer3D::Qsi_11P[6][2]), -4. * Hammer3D::Qsi_11P[6][1] },
	  { 4. * Hammer3D::Qsi_11P[6][1], 4. * Hammer3D::Qsi_11P[6][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[6][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[6][2], -4. * Hammer3D::Qsi_11P[6][2], 4. - 4. * (Hammer3D::Qsi_11P[6][0] + Hammer3D::Qsi_11P[6][1] + 2. * Hammer3D::Qsi_11P[6][2]) },
	  { 4. * Hammer3D::Qsi_11P[6][2], 0., 4. * Hammer3D::Qsi_11P[6][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[6][2], 4. * Hammer3D::Qsi_11P[6][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[6][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[7][0] + Hammer3D::Qsi_11P[7][1] + Hammer3D::Qsi_11P[7][2]), -3. + 4. * (Hammer3D::Qsi_11P[7][0] + Hammer3D::Qsi_11P[7][1] + Hammer3D::Qsi_11P[7][2]), -3. + 4. * (Hammer3D::Qsi_11P[7][0] + Hammer3D::Qsi_11P[7][1] + Hammer3D::Qsi_11P[7][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[7][0] + Hammer3D::Qsi_11P[7][1] + Hammer3D::Qsi_11P[7][2]), -4. * Hammer3D::Qsi_11P[7][0], -4. * Hammer3D::Qsi_11P[7][0] },
	  { 4. * Hammer3D::Qsi_11P[7][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[7][1], 4. - 4. * (Hammer3D::Qsi_11P[7][0] + 2. * Hammer3D::Qsi_11P[7][1] + Hammer3D::Qsi_11P[7][2]), -4. * Hammer3D::Qsi_11P[7][1] },
	  { 4. * Hammer3D::Qsi_11P[7][1], 4. * Hammer3D::Qsi_11P[7][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[7][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[7][2], -4. * Hammer3D::Qsi_11P[7][2], 4. - 4. * (Hammer3D::Qsi_11P[7][0] + Hammer3D::Qsi_11P[7][1] + 2. * Hammer3D::Qsi_11P[7][2]) },
	  { 4. * Hammer3D::Qsi_11P[7][2], 0., 4. * Hammer3D::Qsi_11P[7][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[7][2], 4. * Hammer3D::Qsi_11P[7][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[7][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[8][0] + Hammer3D::Qsi_11P[8][1] + Hammer3D::Qsi_11P[8][2]), -3. + 4. * (Hammer3D::Qsi_11P[8][0] + Hammer3D::Qsi_11P[8][1] + Hammer3D::Qsi_11P[8][2]), -3. + 4. * (Hammer3D::Qsi_11P[8][0] + Hammer3D::Qsi_11P[8][1] + Hammer3D::Qsi_11P[8][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[8][0] + Hammer3D::Qsi_11P[8][1] + Hammer3D::Qsi_11P[8][2]), -4. * Hammer3D::Qsi_11P[8][0], -4. * Hammer3D::Qsi_11P[8][0] },
	  { 4. * Hammer3D::Qsi_11P[8][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[8][1], 4. - 4. * (Hammer3D::Qsi_11P[8][0] + 2. * Hammer3D::Qsi_11P[8][1] + Hammer3D::Qsi_11P[8][2]), -4. * Hammer3D::Qsi_11P[8][1] },
	  { 4. * Hammer3D::Qsi_11P[8][1], 4. * Hammer3D::Qsi_11P[8][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[8][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[8][2], -4. * Hammer3D::Qsi_11P[8][2], 4. - 4. * (Hammer3D::Qsi_11P[8][0] + Hammer3D::Qsi_11P[8][1] + 2. * Hammer3D::Qsi_11P[8][2]) },
	  { 4. * Hammer3D::Qsi_11P[8][2], 0., 4. * Hammer3D::Qsi_11P[8][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[8][2], 4. * Hammer3D::Qsi_11P[8][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[8][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[9][0] + Hammer3D::Qsi_11P[9][1] + Hammer3D::Qsi_11P[9][2]), -3. + 4. * (Hammer3D::Qsi_11P[9][0] + Hammer3D::Qsi_11P[9][1] + Hammer3D::Qsi_11P[9][2]), -3. + 4. * (Hammer3D::Qsi_11P[9][0] + Hammer3D::Qsi_11P[9][1] + Hammer3D::Qsi_11P[9][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[9][0] + Hammer3D::Qsi_11P[9][1] + Hammer3D::Qsi_11P[9][2]), -4. * Hammer3D::Qsi_11P[9][0], -4. * Hammer3D::Qsi_11P[9][0] },
	  { 4. * Hammer3D::Qsi_11P[9][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[9][1], 4. - 4. * (Hammer3D::Qsi_11P[9][0] + 2. * Hammer3D::Qsi_11P[9][1] + Hammer3D::Qsi_11P[9][2]), -4. * Hammer3D::Qsi_11P[9][1] },
	  { 4. * Hammer3D::Qsi_11P[9][1], 4. * Hammer3D::Qsi_11P[9][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[9][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[9][2], -4. * Hammer3D::Qsi_11P[9][2], 4. - 4. * (Hammer3D::Qsi_11P[9][0] + Hammer3D::Qsi_11P[9][1] + 2. * Hammer3D::Qsi_11P[9][2]) },
	  { 4. * Hammer3D::Qsi_11P[9][2], 0., 4. * Hammer3D::Qsi_11P[9][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[9][2], 4. * Hammer3D::Qsi_11P[9][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[9][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_11P[10][0] + Hammer3D::Qsi_11P[10][1] + Hammer3D::Qsi_11P[10][2]), -3. + 4. * (Hammer3D::Qsi_11P[10][0] + Hammer3D::Qsi_11P[10][1] + Hammer3D::Qsi_11P[10][2]), -3. + 4. * (Hammer3D::Qsi_11P[10][0] + Hammer3D::Qsi_11P[10][1] + Hammer3D::Qsi_11P[10][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_11P[10][0] + Hammer3D::Qsi_11P[10][1] + Hammer3D::Qsi_11P[10][2]), -4. * Hammer3D::Qsi_11P[10][0], -4. * Hammer3D::Qsi_11P[10][0] },
	  { 4. * Hammer3D::Qsi_11P[10][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_11P[10][1], 4. - 4. * (Hammer3D::Qsi_11P[10][0] + 2. * Hammer3D::Qsi_11P[10][1] + Hammer3D::Qsi_11P[10][2]), -4. * Hammer3D::Qsi_11P[10][1] },
	  { 4. * Hammer3D::Qsi_11P[10][1], 4. * Hammer3D::Qsi_11P[10][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_11P[10][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_11P[10][2], -4. * Hammer3D::Qsi_11P[10][2], 4. - 4. * (Hammer3D::Qsi_11P[10][0] + Hammer3D::Qsi_11P[10][1] + 2. * Hammer3D::Qsi_11P[10][2]) },
	  { 4. * Hammer3D::Qsi_11P[10][2], 0., 4. * Hammer3D::Qsi_11P[10][0] },
	  { 0., 4. * Hammer3D::Qsi_11P[10][2], 4. * Hammer3D::Qsi_11P[10][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_11P[10][2] - 1.} } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<14>::mv_DPsi[14][mv_numNodes][mv_ElDim] = {
	{ { -3. + 4. * (Hammer3D::Qsi_14P[0][0] + Hammer3D::Qsi_14P[0][1] + Hammer3D::Qsi_14P[0][2]), -3. + 4. * (Hammer3D::Qsi_14P[0][0] + Hammer3D::Qsi_14P[0][1] + Hammer3D::Qsi_14P[0][2]), -3. + 4. * (Hammer3D::Qsi_14P[0][0] + Hammer3D::Qsi_14P[0][1] + Hammer3D::Qsi_14P[0][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[0][0] + Hammer3D::Qsi_14P[0][1] + Hammer3D::Qsi_14P[0][2]), -4. * Hammer3D::Qsi_14P[0][0], -4. * Hammer3D::Qsi_14P[0][0] },
	  { 4. * Hammer3D::Qsi_14P[0][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[0][1], 4. - 4. * (Hammer3D::Qsi_14P[0][0] + 2. * Hammer3D::Qsi_14P[0][1] + Hammer3D::Qsi_14P[0][2]), -4. * Hammer3D::Qsi_14P[0][1] },
	  { 4. * Hammer3D::Qsi_14P[0][1], 4. * Hammer3D::Qsi_14P[0][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[0][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[0][2], -4. * Hammer3D::Qsi_14P[0][2], 4. - 4. * (Hammer3D::Qsi_14P[0][0] + Hammer3D::Qsi_14P[0][1] + 2. * Hammer3D::Qsi_14P[0][2]) },
	  { 4. * Hammer3D::Qsi_14P[0][2], 0., 4. * Hammer3D::Qsi_14P[0][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[0][2], 4. * Hammer3D::Qsi_14P[0][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[0][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[1][0] + Hammer3D::Qsi_14P[1][1] + Hammer3D::Qsi_14P[1][2]), -3. + 4. * (Hammer3D::Qsi_14P[1][0] + Hammer3D::Qsi_14P[1][1] + Hammer3D::Qsi_14P[1][2]), -3. + 4. * (Hammer3D::Qsi_14P[1][0] + Hammer3D::Qsi_14P[1][1] + Hammer3D::Qsi_14P[1][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[1][0] + Hammer3D::Qsi_14P[1][1] + Hammer3D::Qsi_14P[1][2]), -4. * Hammer3D::Qsi_14P[1][0], -4. * Hammer3D::Qsi_14P[1][0] },
	  { 4. * Hammer3D::Qsi_14P[1][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[1][1], 4. - 4. * (Hammer3D::Qsi_14P[1][0] + 2. * Hammer3D::Qsi_14P[1][1] + Hammer3D::Qsi_14P[1][2]), -4. * Hammer3D::Qsi_14P[1][1] },
	  { 4. * Hammer3D::Qsi_14P[1][1], 4. * Hammer3D::Qsi_14P[1][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[1][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[1][2], -4. * Hammer3D::Qsi_14P[1][2], 4. - 4. * (Hammer3D::Qsi_14P[1][0] + Hammer3D::Qsi_14P[1][1] + 2. * Hammer3D::Qsi_14P[1][2]) },
	  { 4. * Hammer3D::Qsi_14P[1][2], 0., 4. * Hammer3D::Qsi_14P[1][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[1][2], 4. * Hammer3D::Qsi_14P[1][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[1][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[2][0] + Hammer3D::Qsi_14P[2][1] + Hammer3D::Qsi_14P[2][2]), -3. + 4. * (Hammer3D::Qsi_14P[2][0] + Hammer3D::Qsi_14P[2][1] + Hammer3D::Qsi_14P[2][2]), -3. + 4. * (Hammer3D::Qsi_14P[2][0] + Hammer3D::Qsi_14P[2][1] + Hammer3D::Qsi_14P[2][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[2][0] + Hammer3D::Qsi_14P[2][1] + Hammer3D::Qsi_14P[2][2]), -4. * Hammer3D::Qsi_14P[2][0], -4. * Hammer3D::Qsi_14P[2][0] },
	  { 4. * Hammer3D::Qsi_14P[2][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[2][1], 4. - 4. * (Hammer3D::Qsi_14P[2][0] + 2. * Hammer3D::Qsi_14P[2][1] + Hammer3D::Qsi_14P[2][2]), -4. * Hammer3D::Qsi_14P[2][1] },
	  { 4. * Hammer3D::Qsi_14P[2][1], 4. * Hammer3D::Qsi_14P[2][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[2][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[2][2], -4. * Hammer3D::Qsi_14P[2][2], 4. - 4. * (Hammer3D::Qsi_14P[2][0] + Hammer3D::Qsi_14P[2][1] + 2. * Hammer3D::Qsi_14P[2][2]) },
	  { 4. * Hammer3D::Qsi_14P[2][2], 0., 4. * Hammer3D::Qsi_14P[2][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[2][2], 4. * Hammer3D::Qsi_14P[2][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[2][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[3][0] + Hammer3D::Qsi_14P[3][1] + Hammer3D::Qsi_14P[3][2]), -3. + 4. * (Hammer3D::Qsi_14P[3][0] + Hammer3D::Qsi_14P[3][1] + Hammer3D::Qsi_14P[3][2]), -3. + 4. * (Hammer3D::Qsi_14P[3][0] + Hammer3D::Qsi_14P[3][1] + Hammer3D::Qsi_14P[3][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[3][0] + Hammer3D::Qsi_14P[3][1] + Hammer3D::Qsi_14P[3][2]), -4. * Hammer3D::Qsi_14P[3][0], -4. * Hammer3D::Qsi_14P[3][0] },
	  { 4. * Hammer3D::Qsi_14P[3][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[3][1], 4. - 4. * (Hammer3D::Qsi_14P[3][0] + 2. * Hammer3D::Qsi_14P[3][1] + Hammer3D::Qsi_14P[3][2]), -4. * Hammer3D::Qsi_14P[3][1] },
	  { 4. * Hammer3D::Qsi_14P[3][1], 4. * Hammer3D::Qsi_14P[3][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[3][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[3][2], -4. * Hammer3D::Qsi_14P[3][2], 4. - 4. * (Hammer3D::Qsi_14P[3][0] + Hammer3D::Qsi_14P[3][1] + 2. * Hammer3D::Qsi_14P[3][2]) },
	  { 4. * Hammer3D::Qsi_14P[3][2], 0., 4. * Hammer3D::Qsi_14P[3][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[3][2], 4. * Hammer3D::Qsi_14P[3][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[3][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[4][0] + Hammer3D::Qsi_14P[4][1] + Hammer3D::Qsi_14P[4][2]), -3. + 4. * (Hammer3D::Qsi_14P[4][0] + Hammer3D::Qsi_14P[4][1] + Hammer3D::Qsi_14P[4][2]), -3. + 4. * (Hammer3D::Qsi_14P[4][0] + Hammer3D::Qsi_14P[4][1] + Hammer3D::Qsi_14P[4][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[4][0] + Hammer3D::Qsi_14P[4][1] + Hammer3D::Qsi_14P[4][2]), -4. * Hammer3D::Qsi_14P[4][0], -4. * Hammer3D::Qsi_14P[4][0] },
	  { 4. * Hammer3D::Qsi_14P[4][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[4][1], 4. - 4. * (Hammer3D::Qsi_14P[4][0] + 2. * Hammer3D::Qsi_14P[4][1] + Hammer3D::Qsi_14P[4][2]), -4. * Hammer3D::Qsi_14P[4][1] },
	  { 4. * Hammer3D::Qsi_14P[4][1], 4. * Hammer3D::Qsi_14P[4][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[4][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[4][2], -4. * Hammer3D::Qsi_14P[4][2], 4. - 4. * (Hammer3D::Qsi_14P[4][0] + Hammer3D::Qsi_14P[4][1] + 2. * Hammer3D::Qsi_14P[4][2]) },
	  { 4. * Hammer3D::Qsi_14P[4][2], 0., 4. * Hammer3D::Qsi_14P[4][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[4][2], 4. * Hammer3D::Qsi_14P[4][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[4][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[5][0] + Hammer3D::Qsi_14P[5][1] + Hammer3D::Qsi_14P[5][2]), -3. + 4. * (Hammer3D::Qsi_14P[5][0] + Hammer3D::Qsi_14P[5][1] + Hammer3D::Qsi_14P[5][2]), -3. + 4. * (Hammer3D::Qsi_14P[5][0] + Hammer3D::Qsi_14P[5][1] + Hammer3D::Qsi_14P[5][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[5][0] + Hammer3D::Qsi_14P[5][1] + Hammer3D::Qsi_14P[5][2]), -4. * Hammer3D::Qsi_14P[5][0], -4. * Hammer3D::Qsi_14P[5][0] },
	  { 4. * Hammer3D::Qsi_14P[5][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[5][1], 4. - 4. * (Hammer3D::Qsi_14P[5][0] + 2. * Hammer3D::Qsi_14P[5][1] + Hammer3D::Qsi_14P[5][2]), -4. * Hammer3D::Qsi_14P[5][1] },
	  { 4. * Hammer3D::Qsi_14P[5][1], 4. * Hammer3D::Qsi_14P[5][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[5][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[5][2], -4. * Hammer3D::Qsi_14P[5][2], 4. - 4. * (Hammer3D::Qsi_14P[5][0] + Hammer3D::Qsi_14P[5][1] + 2. * Hammer3D::Qsi_14P[5][2]) },
	  { 4. * Hammer3D::Qsi_14P[5][2], 0., 4. * Hammer3D::Qsi_14P[5][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[5][2], 4. * Hammer3D::Qsi_14P[5][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[5][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[6][0] + Hammer3D::Qsi_14P[6][1] + Hammer3D::Qsi_14P[6][2]), -3. + 4. * (Hammer3D::Qsi_14P[6][0] + Hammer3D::Qsi_14P[6][1] + Hammer3D::Qsi_14P[6][2]), -3. + 4. * (Hammer3D::Qsi_14P[6][0] + Hammer3D::Qsi_14P[6][1] + Hammer3D::Qsi_14P[6][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[6][0] + Hammer3D::Qsi_14P[6][1] + Hammer3D::Qsi_14P[6][2]), -4. * Hammer3D::Qsi_14P[6][0], -4. * Hammer3D::Qsi_14P[6][0] },
	  { 4. * Hammer3D::Qsi_14P[6][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[6][1], 4. - 4. * (Hammer3D::Qsi_14P[6][0] + 2. * Hammer3D::Qsi_14P[6][1] + Hammer3D::Qsi_14P[6][2]), -4. * Hammer3D::Qsi_14P[6][1] },
	  { 4. * Hammer3D::Qsi_14P[6][1], 4. * Hammer3D::Qsi_14P[6][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[6][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[6][2], -4. * Hammer3D::Qsi_14P[6][2], 4. - 4. * (Hammer3D::Qsi_14P[6][0] + Hammer3D::Qsi_14P[6][1] + 2. * Hammer3D::Qsi_14P[6][2]) },
	  { 4. * Hammer3D::Qsi_14P[6][2], 0., 4. * Hammer3D::Qsi_14P[6][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[6][2], 4. * Hammer3D::Qsi_14P[6][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[6][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[7][0] + Hammer3D::Qsi_14P[7][1] + Hammer3D::Qsi_14P[7][2]), -3. + 4. * (Hammer3D::Qsi_14P[7][0] + Hammer3D::Qsi_14P[7][1] + Hammer3D::Qsi_14P[7][2]), -3. + 4. * (Hammer3D::Qsi_14P[7][0] + Hammer3D::Qsi_14P[7][1] + Hammer3D::Qsi_14P[7][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[7][0] + Hammer3D::Qsi_14P[7][1] + Hammer3D::Qsi_14P[7][2]), -4. * Hammer3D::Qsi_14P[7][0], -4. * Hammer3D::Qsi_14P[7][0] },
	  { 4. * Hammer3D::Qsi_14P[7][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[7][1], 4. - 4. * (Hammer3D::Qsi_14P[7][0] + 2. * Hammer3D::Qsi_14P[7][1] + Hammer3D::Qsi_14P[7][2]), -4. * Hammer3D::Qsi_14P[7][1] },
	  { 4. * Hammer3D::Qsi_14P[7][1], 4. * Hammer3D::Qsi_14P[7][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[7][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[7][2], -4. * Hammer3D::Qsi_14P[7][2], 4. - 4. * (Hammer3D::Qsi_14P[7][0] + Hammer3D::Qsi_14P[7][1] + 2. * Hammer3D::Qsi_14P[7][2]) },
	  { 4. * Hammer3D::Qsi_14P[7][2], 0., 4. * Hammer3D::Qsi_14P[7][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[7][2], 4. * Hammer3D::Qsi_14P[7][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[7][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[8][0] + Hammer3D::Qsi_14P[8][1] + Hammer3D::Qsi_14P[8][2]), -3. + 4. * (Hammer3D::Qsi_14P[8][0] + Hammer3D::Qsi_14P[8][1] + Hammer3D::Qsi_14P[8][2]), -3. + 4. * (Hammer3D::Qsi_14P[8][0] + Hammer3D::Qsi_14P[8][1] + Hammer3D::Qsi_14P[8][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[8][0] + Hammer3D::Qsi_14P[8][1] + Hammer3D::Qsi_14P[8][2]), -4. * Hammer3D::Qsi_14P[8][0], -4. * Hammer3D::Qsi_14P[8][0] },
	  { 4. * Hammer3D::Qsi_14P[8][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[8][1], 4. - 4. * (Hammer3D::Qsi_14P[8][0] + 2. * Hammer3D::Qsi_14P[8][1] + Hammer3D::Qsi_14P[8][2]), -4. * Hammer3D::Qsi_14P[8][1] },
	  { 4. * Hammer3D::Qsi_14P[8][1], 4. * Hammer3D::Qsi_14P[8][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[8][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[8][2], -4. * Hammer3D::Qsi_14P[8][2], 4. - 4. * (Hammer3D::Qsi_14P[8][0] + Hammer3D::Qsi_14P[8][1] + 2. * Hammer3D::Qsi_14P[8][2]) },
	  { 4. * Hammer3D::Qsi_14P[8][2], 0., 4. * Hammer3D::Qsi_14P[8][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[8][2], 4. * Hammer3D::Qsi_14P[8][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[8][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[9][0] + Hammer3D::Qsi_14P[9][1] + Hammer3D::Qsi_14P[9][2]), -3. + 4. * (Hammer3D::Qsi_14P[9][0] + Hammer3D::Qsi_14P[9][1] + Hammer3D::Qsi_14P[9][2]), -3. + 4. * (Hammer3D::Qsi_14P[9][0] + Hammer3D::Qsi_14P[9][1] + Hammer3D::Qsi_14P[9][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[9][0] + Hammer3D::Qsi_14P[9][1] + Hammer3D::Qsi_14P[9][2]), -4. * Hammer3D::Qsi_14P[9][0], -4. * Hammer3D::Qsi_14P[9][0] },
	  { 4. * Hammer3D::Qsi_14P[9][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[9][1], 4. - 4. * (Hammer3D::Qsi_14P[9][0] + 2. * Hammer3D::Qsi_14P[9][1] + Hammer3D::Qsi_14P[9][2]), -4. * Hammer3D::Qsi_14P[9][1] },
	  { 4. * Hammer3D::Qsi_14P[9][1], 4. * Hammer3D::Qsi_14P[9][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[9][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[9][2], -4. * Hammer3D::Qsi_14P[9][2], 4. - 4. * (Hammer3D::Qsi_14P[9][0] + Hammer3D::Qsi_14P[9][1] + 2. * Hammer3D::Qsi_14P[9][2]) },
	  { 4. * Hammer3D::Qsi_14P[9][2], 0., 4. * Hammer3D::Qsi_14P[9][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[9][2], 4. * Hammer3D::Qsi_14P[9][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[9][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[10][0] + Hammer3D::Qsi_14P[10][1] + Hammer3D::Qsi_14P[10][2]), -3. + 4. * (Hammer3D::Qsi_14P[10][0] + Hammer3D::Qsi_14P[10][1] + Hammer3D::Qsi_14P[10][2]), -3. + 4. * (Hammer3D::Qsi_14P[10][0] + Hammer3D::Qsi_14P[10][1] + Hammer3D::Qsi_14P[10][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[10][0] + Hammer3D::Qsi_14P[10][1] + Hammer3D::Qsi_14P[10][2]), -4. * Hammer3D::Qsi_14P[10][0], -4. * Hammer3D::Qsi_14P[10][0] },
	  { 4. * Hammer3D::Qsi_14P[10][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[10][1], 4. - 4. * (Hammer3D::Qsi_14P[10][0] + 2. * Hammer3D::Qsi_14P[10][1] + Hammer3D::Qsi_14P[10][2]), -4. * Hammer3D::Qsi_14P[10][1] },
	  { 4. * Hammer3D::Qsi_14P[10][1], 4. * Hammer3D::Qsi_14P[10][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[10][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[10][2], -4. * Hammer3D::Qsi_14P[10][2], 4. - 4. * (Hammer3D::Qsi_14P[10][0] + Hammer3D::Qsi_14P[10][1] + 2. * Hammer3D::Qsi_14P[10][2]) },
	  { 4. * Hammer3D::Qsi_14P[10][2], 0., 4. * Hammer3D::Qsi_14P[10][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[10][2], 4. * Hammer3D::Qsi_14P[10][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[10][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[11][0] + Hammer3D::Qsi_14P[11][1] + Hammer3D::Qsi_14P[11][2]), -3. + 4. * (Hammer3D::Qsi_14P[11][0] + Hammer3D::Qsi_14P[11][1] + Hammer3D::Qsi_14P[11][2]), -3. + 4. * (Hammer3D::Qsi_14P[11][0] + Hammer3D::Qsi_14P[11][1] + Hammer3D::Qsi_14P[11][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[11][0] + Hammer3D::Qsi_14P[11][1] + Hammer3D::Qsi_14P[11][2]), -4. * Hammer3D::Qsi_14P[11][0], -4. * Hammer3D::Qsi_14P[11][0] },
	  { 4. * Hammer3D::Qsi_14P[11][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[11][1], 4. - 4. * (Hammer3D::Qsi_14P[11][0] + 2. * Hammer3D::Qsi_14P[11][1] + Hammer3D::Qsi_14P[11][2]), -4. * Hammer3D::Qsi_14P[11][1] },
	  { 4. * Hammer3D::Qsi_14P[11][1], 4. * Hammer3D::Qsi_14P[11][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[11][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[11][2], -4. * Hammer3D::Qsi_14P[11][2], 4. - 4. * (Hammer3D::Qsi_14P[11][0] + Hammer3D::Qsi_14P[11][1] + 2. * Hammer3D::Qsi_14P[11][2]) },
	  { 4. * Hammer3D::Qsi_14P[11][2], 0., 4. * Hammer3D::Qsi_14P[11][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[11][2], 4. * Hammer3D::Qsi_14P[11][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[11][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[12][0] + Hammer3D::Qsi_14P[12][1] + Hammer3D::Qsi_14P[12][2]), -3. + 4. * (Hammer3D::Qsi_14P[12][0] + Hammer3D::Qsi_14P[12][1] + Hammer3D::Qsi_14P[12][2]), -3. + 4. * (Hammer3D::Qsi_14P[12][0] + Hammer3D::Qsi_14P[12][1] + Hammer3D::Qsi_14P[12][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[12][0] + Hammer3D::Qsi_14P[12][1] + Hammer3D::Qsi_14P[12][2]), -4. * Hammer3D::Qsi_14P[12][0], -4. * Hammer3D::Qsi_14P[12][0] },
	  { 4. * Hammer3D::Qsi_14P[12][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[12][1], 4. - 4. * (Hammer3D::Qsi_14P[12][0] + 2. * Hammer3D::Qsi_14P[12][1] + Hammer3D::Qsi_14P[12][2]), -4. * Hammer3D::Qsi_14P[12][1] },
	  { 4. * Hammer3D::Qsi_14P[12][1], 4. * Hammer3D::Qsi_14P[12][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[12][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[12][2], -4. * Hammer3D::Qsi_14P[12][2], 4. - 4. * (Hammer3D::Qsi_14P[12][0] + Hammer3D::Qsi_14P[12][1] + 2. * Hammer3D::Qsi_14P[12][2]) },
	  { 4. * Hammer3D::Qsi_14P[12][2], 0., 4. * Hammer3D::Qsi_14P[12][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[12][2], 4. * Hammer3D::Qsi_14P[12][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[12][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_14P[13][0] + Hammer3D::Qsi_14P[13][1] + Hammer3D::Qsi_14P[13][2]), -3. + 4. * (Hammer3D::Qsi_14P[13][0] + Hammer3D::Qsi_14P[13][1] + Hammer3D::Qsi_14P[13][2]), -3. + 4. * (Hammer3D::Qsi_14P[13][0] + Hammer3D::Qsi_14P[13][1] + Hammer3D::Qsi_14P[13][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_14P[13][0] + Hammer3D::Qsi_14P[13][1] + Hammer3D::Qsi_14P[13][2]), -4. * Hammer3D::Qsi_14P[13][0], -4. * Hammer3D::Qsi_14P[13][0] },
	  { 4. * Hammer3D::Qsi_14P[13][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_14P[13][1], 4. - 4. * (Hammer3D::Qsi_14P[13][0] + 2. * Hammer3D::Qsi_14P[13][1] + Hammer3D::Qsi_14P[13][2]), -4. * Hammer3D::Qsi_14P[13][1] },
	  { 4. * Hammer3D::Qsi_14P[13][1], 4. * Hammer3D::Qsi_14P[13][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_14P[13][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_14P[13][2], -4. * Hammer3D::Qsi_14P[13][2], 4. - 4. * (Hammer3D::Qsi_14P[13][0] + Hammer3D::Qsi_14P[13][1] + 2. * Hammer3D::Qsi_14P[13][2]) },
	  { 4. * Hammer3D::Qsi_14P[13][2], 0., 4. * Hammer3D::Qsi_14P[13][0] },
	  { 0., 4. * Hammer3D::Qsi_14P[13][2], 4. * Hammer3D::Qsi_14P[13][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_14P[13][2] - 1.} } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<15>::mv_DPsi[15][mv_numNodes][mv_ElDim] = {
	{ { -3. + 4. * (Hammer3D::Qsi_15P[0][0] + Hammer3D::Qsi_15P[0][1] + Hammer3D::Qsi_15P[0][2]), -3. + 4. * (Hammer3D::Qsi_15P[0][0] + Hammer3D::Qsi_15P[0][1] + Hammer3D::Qsi_15P[0][2]), -3. + 4. * (Hammer3D::Qsi_15P[0][0] + Hammer3D::Qsi_15P[0][1] + Hammer3D::Qsi_15P[0][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[0][0] + Hammer3D::Qsi_15P[0][1] + Hammer3D::Qsi_15P[0][2]), -4. * Hammer3D::Qsi_15P[0][0], -4. * Hammer3D::Qsi_15P[0][0] },
	  { 4. * Hammer3D::Qsi_15P[0][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[0][1], 4. - 4. * (Hammer3D::Qsi_15P[0][0] + 2. * Hammer3D::Qsi_15P[0][1] + Hammer3D::Qsi_15P[0][2]), -4. * Hammer3D::Qsi_15P[0][1] },
	  { 4. * Hammer3D::Qsi_15P[0][1], 4. * Hammer3D::Qsi_15P[0][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[0][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[0][2], -4. * Hammer3D::Qsi_15P[0][2], 4. - 4. * (Hammer3D::Qsi_15P[0][0] + Hammer3D::Qsi_15P[0][1] + 2. * Hammer3D::Qsi_15P[0][2]) },
	  { 4. * Hammer3D::Qsi_15P[0][2], 0., 4. * Hammer3D::Qsi_15P[0][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[0][2], 4. * Hammer3D::Qsi_15P[0][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[0][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[1][0] + Hammer3D::Qsi_15P[1][1] + Hammer3D::Qsi_15P[1][2]), -3. + 4. * (Hammer3D::Qsi_15P[1][0] + Hammer3D::Qsi_15P[1][1] + Hammer3D::Qsi_15P[1][2]), -3. + 4. * (Hammer3D::Qsi_15P[1][0] + Hammer3D::Qsi_15P[1][1] + Hammer3D::Qsi_15P[1][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[1][0] + Hammer3D::Qsi_15P[1][1] + Hammer3D::Qsi_15P[1][2]), -4. * Hammer3D::Qsi_15P[1][0], -4. * Hammer3D::Qsi_15P[1][0] },
	  { 4. * Hammer3D::Qsi_15P[1][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[1][1], 4. - 4. * (Hammer3D::Qsi_15P[1][0] + 2. * Hammer3D::Qsi_15P[1][1] + Hammer3D::Qsi_15P[1][2]), -4. * Hammer3D::Qsi_15P[1][1] },
	  { 4. * Hammer3D::Qsi_15P[1][1], 4. * Hammer3D::Qsi_15P[1][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[1][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[1][2], -4. * Hammer3D::Qsi_15P[1][2], 4. - 4. * (Hammer3D::Qsi_15P[1][0] + Hammer3D::Qsi_15P[1][1] + 2. * Hammer3D::Qsi_15P[1][2]) },
	  { 4. * Hammer3D::Qsi_15P[1][2], 0., 4. * Hammer3D::Qsi_15P[1][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[1][2], 4. * Hammer3D::Qsi_15P[1][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[1][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[2][0] + Hammer3D::Qsi_15P[2][1] + Hammer3D::Qsi_15P[2][2]), -3. + 4. * (Hammer3D::Qsi_15P[2][0] + Hammer3D::Qsi_15P[2][1] + Hammer3D::Qsi_15P[2][2]), -3. + 4. * (Hammer3D::Qsi_15P[2][0] + Hammer3D::Qsi_15P[2][1] + Hammer3D::Qsi_15P[2][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[2][0] + Hammer3D::Qsi_15P[2][1] + Hammer3D::Qsi_15P[2][2]), -4. * Hammer3D::Qsi_15P[2][0], -4. * Hammer3D::Qsi_15P[2][0] },
	  { 4. * Hammer3D::Qsi_15P[2][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[2][1], 4. - 4. * (Hammer3D::Qsi_15P[2][0] + 2. * Hammer3D::Qsi_15P[2][1] + Hammer3D::Qsi_15P[2][2]), -4. * Hammer3D::Qsi_15P[2][1] },
	  { 4. * Hammer3D::Qsi_15P[2][1], 4. * Hammer3D::Qsi_15P[2][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[2][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[2][2], -4. * Hammer3D::Qsi_15P[2][2], 4. - 4. * (Hammer3D::Qsi_15P[2][0] + Hammer3D::Qsi_15P[2][1] + 2. * Hammer3D::Qsi_15P[2][2]) },
	  { 4. * Hammer3D::Qsi_15P[2][2], 0., 4. * Hammer3D::Qsi_15P[2][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[2][2], 4. * Hammer3D::Qsi_15P[2][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[2][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[3][0] + Hammer3D::Qsi_15P[3][1] + Hammer3D::Qsi_15P[3][2]), -3. + 4. * (Hammer3D::Qsi_15P[3][0] + Hammer3D::Qsi_15P[3][1] + Hammer3D::Qsi_15P[3][2]), -3. + 4. * (Hammer3D::Qsi_15P[3][0] + Hammer3D::Qsi_15P[3][1] + Hammer3D::Qsi_15P[3][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[3][0] + Hammer3D::Qsi_15P[3][1] + Hammer3D::Qsi_15P[3][2]), -4. * Hammer3D::Qsi_15P[3][0], -4. * Hammer3D::Qsi_15P[3][0] },
	  { 4. * Hammer3D::Qsi_15P[3][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[3][1], 4. - 4. * (Hammer3D::Qsi_15P[3][0] + 2. * Hammer3D::Qsi_15P[3][1] + Hammer3D::Qsi_15P[3][2]), -4. * Hammer3D::Qsi_15P[3][1] },
	  { 4. * Hammer3D::Qsi_15P[3][1], 4. * Hammer3D::Qsi_15P[3][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[3][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[3][2], -4. * Hammer3D::Qsi_15P[3][2], 4. - 4. * (Hammer3D::Qsi_15P[3][0] + Hammer3D::Qsi_15P[3][1] + 2. * Hammer3D::Qsi_15P[3][2]) },
	  { 4. * Hammer3D::Qsi_15P[3][2], 0., 4. * Hammer3D::Qsi_15P[3][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[3][2], 4. * Hammer3D::Qsi_15P[3][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[3][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[4][0] + Hammer3D::Qsi_15P[4][1] + Hammer3D::Qsi_15P[4][2]), -3. + 4. * (Hammer3D::Qsi_15P[4][0] + Hammer3D::Qsi_15P[4][1] + Hammer3D::Qsi_15P[4][2]), -3. + 4. * (Hammer3D::Qsi_15P[4][0] + Hammer3D::Qsi_15P[4][1] + Hammer3D::Qsi_15P[4][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[4][0] + Hammer3D::Qsi_15P[4][1] + Hammer3D::Qsi_15P[4][2]), -4. * Hammer3D::Qsi_15P[4][0], -4. * Hammer3D::Qsi_15P[4][0] },
	  { 4. * Hammer3D::Qsi_15P[4][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[4][1], 4. - 4. * (Hammer3D::Qsi_15P[4][0] + 2. * Hammer3D::Qsi_15P[4][1] + Hammer3D::Qsi_15P[4][2]), -4. * Hammer3D::Qsi_15P[4][1] },
	  { 4. * Hammer3D::Qsi_15P[4][1], 4. * Hammer3D::Qsi_15P[4][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[4][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[4][2], -4. * Hammer3D::Qsi_15P[4][2], 4. - 4. * (Hammer3D::Qsi_15P[4][0] + Hammer3D::Qsi_15P[4][1] + 2. * Hammer3D::Qsi_15P[4][2]) },
	  { 4. * Hammer3D::Qsi_15P[4][2], 0., 4. * Hammer3D::Qsi_15P[4][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[4][2], 4. * Hammer3D::Qsi_15P[4][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[4][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[5][0] + Hammer3D::Qsi_15P[5][1] + Hammer3D::Qsi_15P[5][2]), -3. + 4. * (Hammer3D::Qsi_15P[5][0] + Hammer3D::Qsi_15P[5][1] + Hammer3D::Qsi_15P[5][2]), -3. + 4. * (Hammer3D::Qsi_15P[5][0] + Hammer3D::Qsi_15P[5][1] + Hammer3D::Qsi_15P[5][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[5][0] + Hammer3D::Qsi_15P[5][1] + Hammer3D::Qsi_15P[5][2]), -4. * Hammer3D::Qsi_15P[5][0], -4. * Hammer3D::Qsi_15P[5][0] },
	  { 4. * Hammer3D::Qsi_15P[5][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[5][1], 4. - 4. * (Hammer3D::Qsi_15P[5][0] + 2. * Hammer3D::Qsi_15P[5][1] + Hammer3D::Qsi_15P[5][2]), -4. * Hammer3D::Qsi_15P[5][1] },
	  { 4. * Hammer3D::Qsi_15P[5][1], 4. * Hammer3D::Qsi_15P[5][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[5][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[5][2], -4. * Hammer3D::Qsi_15P[5][2], 4. - 4. * (Hammer3D::Qsi_15P[5][0] + Hammer3D::Qsi_15P[5][1] + 2. * Hammer3D::Qsi_15P[5][2]) },
	  { 4. * Hammer3D::Qsi_15P[5][2], 0., 4. * Hammer3D::Qsi_15P[5][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[5][2], 4. * Hammer3D::Qsi_15P[5][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[5][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[6][0] + Hammer3D::Qsi_15P[6][1] + Hammer3D::Qsi_15P[6][2]), -3. + 4. * (Hammer3D::Qsi_15P[6][0] + Hammer3D::Qsi_15P[6][1] + Hammer3D::Qsi_15P[6][2]), -3. + 4. * (Hammer3D::Qsi_15P[6][0] + Hammer3D::Qsi_15P[6][1] + Hammer3D::Qsi_15P[6][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[6][0] + Hammer3D::Qsi_15P[6][1] + Hammer3D::Qsi_15P[6][2]), -4. * Hammer3D::Qsi_15P[6][0], -4. * Hammer3D::Qsi_15P[6][0] },
	  { 4. * Hammer3D::Qsi_15P[6][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[6][1], 4. - 4. * (Hammer3D::Qsi_15P[6][0] + 2. * Hammer3D::Qsi_15P[6][1] + Hammer3D::Qsi_15P[6][2]), -4. * Hammer3D::Qsi_15P[6][1] },
	  { 4. * Hammer3D::Qsi_15P[6][1], 4. * Hammer3D::Qsi_15P[6][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[6][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[6][2], -4. * Hammer3D::Qsi_15P[6][2], 4. - 4. * (Hammer3D::Qsi_15P[6][0] + Hammer3D::Qsi_15P[6][1] + 2. * Hammer3D::Qsi_15P[6][2]) },
	  { 4. * Hammer3D::Qsi_15P[6][2], 0., 4. * Hammer3D::Qsi_15P[6][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[6][2], 4. * Hammer3D::Qsi_15P[6][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[6][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[7][0] + Hammer3D::Qsi_15P[7][1] + Hammer3D::Qsi_15P[7][2]), -3. + 4. * (Hammer3D::Qsi_15P[7][0] + Hammer3D::Qsi_15P[7][1] + Hammer3D::Qsi_15P[7][2]), -3. + 4. * (Hammer3D::Qsi_15P[7][0] + Hammer3D::Qsi_15P[7][1] + Hammer3D::Qsi_15P[7][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[7][0] + Hammer3D::Qsi_15P[7][1] + Hammer3D::Qsi_15P[7][2]), -4. * Hammer3D::Qsi_15P[7][0], -4. * Hammer3D::Qsi_15P[7][0] },
	  { 4. * Hammer3D::Qsi_15P[7][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[7][1], 4. - 4. * (Hammer3D::Qsi_15P[7][0] + 2. * Hammer3D::Qsi_15P[7][1] + Hammer3D::Qsi_15P[7][2]), -4. * Hammer3D::Qsi_15P[7][1] },
	  { 4. * Hammer3D::Qsi_15P[7][1], 4. * Hammer3D::Qsi_15P[7][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[7][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[7][2], -4. * Hammer3D::Qsi_15P[7][2], 4. - 4. * (Hammer3D::Qsi_15P[7][0] + Hammer3D::Qsi_15P[7][1] + 2. * Hammer3D::Qsi_15P[7][2]) },
	  { 4. * Hammer3D::Qsi_15P[7][2], 0., 4. * Hammer3D::Qsi_15P[7][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[7][2], 4. * Hammer3D::Qsi_15P[7][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[7][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[8][0] + Hammer3D::Qsi_15P[8][1] + Hammer3D::Qsi_15P[8][2]), -3. + 4. * (Hammer3D::Qsi_15P[8][0] + Hammer3D::Qsi_15P[8][1] + Hammer3D::Qsi_15P[8][2]), -3. + 4. * (Hammer3D::Qsi_15P[8][0] + Hammer3D::Qsi_15P[8][1] + Hammer3D::Qsi_15P[8][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[8][0] + Hammer3D::Qsi_15P[8][1] + Hammer3D::Qsi_15P[8][2]), -4. * Hammer3D::Qsi_15P[8][0], -4. * Hammer3D::Qsi_15P[8][0] },
	  { 4. * Hammer3D::Qsi_15P[8][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[8][1], 4. - 4. * (Hammer3D::Qsi_15P[8][0] + 2. * Hammer3D::Qsi_15P[8][1] + Hammer3D::Qsi_15P[8][2]), -4. * Hammer3D::Qsi_15P[8][1] },
	  { 4. * Hammer3D::Qsi_15P[8][1], 4. * Hammer3D::Qsi_15P[8][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[8][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[8][2], -4. * Hammer3D::Qsi_15P[8][2], 4. - 4. * (Hammer3D::Qsi_15P[8][0] + Hammer3D::Qsi_15P[8][1] + 2. * Hammer3D::Qsi_15P[8][2]) },
	  { 4. * Hammer3D::Qsi_15P[8][2], 0., 4. * Hammer3D::Qsi_15P[8][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[8][2], 4. * Hammer3D::Qsi_15P[8][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[8][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[9][0] + Hammer3D::Qsi_15P[9][1] + Hammer3D::Qsi_15P[9][2]), -3. + 4. * (Hammer3D::Qsi_15P[9][0] + Hammer3D::Qsi_15P[9][1] + Hammer3D::Qsi_15P[9][2]), -3. + 4. * (Hammer3D::Qsi_15P[9][0] + Hammer3D::Qsi_15P[9][1] + Hammer3D::Qsi_15P[9][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[9][0] + Hammer3D::Qsi_15P[9][1] + Hammer3D::Qsi_15P[9][2]), -4. * Hammer3D::Qsi_15P[9][0], -4. * Hammer3D::Qsi_15P[9][0] },
	  { 4. * Hammer3D::Qsi_15P[9][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[9][1], 4. - 4. * (Hammer3D::Qsi_15P[9][0] + 2. * Hammer3D::Qsi_15P[9][1] + Hammer3D::Qsi_15P[9][2]), -4. * Hammer3D::Qsi_15P[9][1] },
	  { 4. * Hammer3D::Qsi_15P[9][1], 4. * Hammer3D::Qsi_15P[9][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[9][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[9][2], -4. * Hammer3D::Qsi_15P[9][2], 4. - 4. * (Hammer3D::Qsi_15P[9][0] + Hammer3D::Qsi_15P[9][1] + 2. * Hammer3D::Qsi_15P[9][2]) },
	  { 4. * Hammer3D::Qsi_15P[9][2], 0., 4. * Hammer3D::Qsi_15P[9][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[9][2], 4. * Hammer3D::Qsi_15P[9][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[9][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[10][0] + Hammer3D::Qsi_15P[10][1] + Hammer3D::Qsi_15P[10][2]), -3. + 4. * (Hammer3D::Qsi_15P[10][0] + Hammer3D::Qsi_15P[10][1] + Hammer3D::Qsi_15P[10][2]), -3. + 4. * (Hammer3D::Qsi_15P[10][0] + Hammer3D::Qsi_15P[10][1] + Hammer3D::Qsi_15P[10][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[10][0] + Hammer3D::Qsi_15P[10][1] + Hammer3D::Qsi_15P[10][2]), -4. * Hammer3D::Qsi_15P[10][0], -4. * Hammer3D::Qsi_15P[10][0] },
	  { 4. * Hammer3D::Qsi_15P[10][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[10][1], 4. - 4. * (Hammer3D::Qsi_15P[10][0] + 2. * Hammer3D::Qsi_15P[10][1] + Hammer3D::Qsi_15P[10][2]), -4. * Hammer3D::Qsi_15P[10][1] },
	  { 4. * Hammer3D::Qsi_15P[10][1], 4. * Hammer3D::Qsi_15P[10][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[10][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[10][2], -4. * Hammer3D::Qsi_15P[10][2], 4. - 4. * (Hammer3D::Qsi_15P[10][0] + Hammer3D::Qsi_15P[10][1] + 2. * Hammer3D::Qsi_15P[10][2]) },
	  { 4. * Hammer3D::Qsi_15P[10][2], 0., 4. * Hammer3D::Qsi_15P[10][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[10][2], 4. * Hammer3D::Qsi_15P[10][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[10][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[11][0] + Hammer3D::Qsi_15P[11][1] + Hammer3D::Qsi_15P[11][2]), -3. + 4. * (Hammer3D::Qsi_15P[11][0] + Hammer3D::Qsi_15P[11][1] + Hammer3D::Qsi_15P[11][2]), -3. + 4. * (Hammer3D::Qsi_15P[11][0] + Hammer3D::Qsi_15P[11][1] + Hammer3D::Qsi_15P[11][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[11][0] + Hammer3D::Qsi_15P[11][1] + Hammer3D::Qsi_15P[11][2]), -4. * Hammer3D::Qsi_15P[11][0], -4. * Hammer3D::Qsi_15P[11][0] },
	  { 4. * Hammer3D::Qsi_15P[11][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[11][1], 4. - 4. * (Hammer3D::Qsi_15P[11][0] + 2. * Hammer3D::Qsi_15P[11][1] + Hammer3D::Qsi_15P[11][2]), -4. * Hammer3D::Qsi_15P[11][1] },
	  { 4. * Hammer3D::Qsi_15P[11][1], 4. * Hammer3D::Qsi_15P[11][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[11][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[11][2], -4. * Hammer3D::Qsi_15P[11][2], 4. - 4. * (Hammer3D::Qsi_15P[11][0] + Hammer3D::Qsi_15P[11][1] + 2. * Hammer3D::Qsi_15P[11][2]) },
	  { 4. * Hammer3D::Qsi_15P[11][2], 0., 4. * Hammer3D::Qsi_15P[11][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[11][2], 4. * Hammer3D::Qsi_15P[11][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[11][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[12][0] + Hammer3D::Qsi_15P[12][1] + Hammer3D::Qsi_15P[12][2]), -3. + 4. * (Hammer3D::Qsi_15P[12][0] + Hammer3D::Qsi_15P[12][1] + Hammer3D::Qsi_15P[12][2]), -3. + 4. * (Hammer3D::Qsi_15P[12][0] + Hammer3D::Qsi_15P[12][1] + Hammer3D::Qsi_15P[12][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[12][0] + Hammer3D::Qsi_15P[12][1] + Hammer3D::Qsi_15P[12][2]), -4. * Hammer3D::Qsi_15P[12][0], -4. * Hammer3D::Qsi_15P[12][0] },
	  { 4. * Hammer3D::Qsi_15P[12][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[12][1], 4. - 4. * (Hammer3D::Qsi_15P[12][0] + 2. * Hammer3D::Qsi_15P[12][1] + Hammer3D::Qsi_15P[12][2]), -4. * Hammer3D::Qsi_15P[12][1] },
	  { 4. * Hammer3D::Qsi_15P[12][1], 4. * Hammer3D::Qsi_15P[12][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[12][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[12][2], -4. * Hammer3D::Qsi_15P[12][2], 4. - 4. * (Hammer3D::Qsi_15P[12][0] + Hammer3D::Qsi_15P[12][1] + 2. * Hammer3D::Qsi_15P[12][2]) },
	  { 4. * Hammer3D::Qsi_15P[12][2], 0., 4. * Hammer3D::Qsi_15P[12][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[12][2], 4. * Hammer3D::Qsi_15P[12][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[12][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[13][0] + Hammer3D::Qsi_15P[13][1] + Hammer3D::Qsi_15P[13][2]), -3. + 4. * (Hammer3D::Qsi_15P[13][0] + Hammer3D::Qsi_15P[13][1] + Hammer3D::Qsi_15P[13][2]), -3. + 4. * (Hammer3D::Qsi_15P[13][0] + Hammer3D::Qsi_15P[13][1] + Hammer3D::Qsi_15P[13][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[13][0] + Hammer3D::Qsi_15P[13][1] + Hammer3D::Qsi_15P[13][2]), -4. * Hammer3D::Qsi_15P[13][0], -4. * Hammer3D::Qsi_15P[13][0] },
	  { 4. * Hammer3D::Qsi_15P[13][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[13][1], 4. - 4. * (Hammer3D::Qsi_15P[13][0] + 2. * Hammer3D::Qsi_15P[13][1] + Hammer3D::Qsi_15P[13][2]), -4. * Hammer3D::Qsi_15P[13][1] },
	  { 4. * Hammer3D::Qsi_15P[13][1], 4. * Hammer3D::Qsi_15P[13][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[13][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[13][2], -4. * Hammer3D::Qsi_15P[13][2], 4. - 4. * (Hammer3D::Qsi_15P[13][0] + Hammer3D::Qsi_15P[13][1] + 2. * Hammer3D::Qsi_15P[13][2]) },
	  { 4. * Hammer3D::Qsi_15P[13][2], 0., 4. * Hammer3D::Qsi_15P[13][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[13][2], 4. * Hammer3D::Qsi_15P[13][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[13][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_15P[14][0] + Hammer3D::Qsi_15P[14][1] + Hammer3D::Qsi_15P[14][2]), -3. + 4. * (Hammer3D::Qsi_15P[14][0] + Hammer3D::Qsi_15P[14][1] + Hammer3D::Qsi_15P[14][2]), -3. + 4. * (Hammer3D::Qsi_15P[14][0] + Hammer3D::Qsi_15P[14][1] + Hammer3D::Qsi_15P[14][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_15P[14][0] + Hammer3D::Qsi_15P[14][1] + Hammer3D::Qsi_15P[14][2]), -4. * Hammer3D::Qsi_15P[14][0], -4. * Hammer3D::Qsi_15P[14][0] },
	  { 4. * Hammer3D::Qsi_15P[14][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_15P[14][1], 4. - 4. * (Hammer3D::Qsi_15P[14][0] + 2. * Hammer3D::Qsi_15P[14][1] + Hammer3D::Qsi_15P[14][2]), -4. * Hammer3D::Qsi_15P[14][1] },
	  { 4. * Hammer3D::Qsi_15P[14][1], 4. * Hammer3D::Qsi_15P[14][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_15P[14][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_15P[14][2], -4. * Hammer3D::Qsi_15P[14][2], 4. - 4. * (Hammer3D::Qsi_15P[14][0] + Hammer3D::Qsi_15P[14][1] + 2. * Hammer3D::Qsi_15P[14][2]) },
	  { 4. * Hammer3D::Qsi_15P[14][2], 0., 4. * Hammer3D::Qsi_15P[14][0] },
	  { 0., 4. * Hammer3D::Qsi_15P[14][2], 4. * Hammer3D::Qsi_15P[14][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_15P[14][2] - 1.} } };

template<> const double O2P2::Geom::Elem::Elem_Tet10_IP<24>::mv_DPsi[24][mv_numNodes][mv_ElDim] = {
	{ { -3. + 4. * (Hammer3D::Qsi_24P[0][0] + Hammer3D::Qsi_24P[0][1] + Hammer3D::Qsi_24P[0][2]), -3. + 4. * (Hammer3D::Qsi_24P[0][0] + Hammer3D::Qsi_24P[0][1] + Hammer3D::Qsi_24P[0][2]), -3. + 4. * (Hammer3D::Qsi_24P[0][0] + Hammer3D::Qsi_24P[0][1] + Hammer3D::Qsi_24P[0][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[0][0] + Hammer3D::Qsi_24P[0][1] + Hammer3D::Qsi_24P[0][2]), -4. * Hammer3D::Qsi_24P[0][0], -4. * Hammer3D::Qsi_24P[0][0] },
	  { 4. * Hammer3D::Qsi_24P[0][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[0][1], 4. - 4. * (Hammer3D::Qsi_24P[0][0] + 2. * Hammer3D::Qsi_24P[0][1] + Hammer3D::Qsi_24P[0][2]), -4. * Hammer3D::Qsi_24P[0][1] },
	  { 4. * Hammer3D::Qsi_24P[0][1], 4. * Hammer3D::Qsi_24P[0][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[0][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[0][2], -4. * Hammer3D::Qsi_24P[0][2], 4. - 4. * (Hammer3D::Qsi_24P[0][0] + Hammer3D::Qsi_24P[0][1] + 2. * Hammer3D::Qsi_24P[0][2]) },
	  { 4. * Hammer3D::Qsi_24P[0][2], 0., 4. * Hammer3D::Qsi_24P[0][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[0][2], 4. * Hammer3D::Qsi_24P[0][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[0][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[1][0] + Hammer3D::Qsi_24P[1][1] + Hammer3D::Qsi_24P[1][2]), -3. + 4. * (Hammer3D::Qsi_24P[1][0] + Hammer3D::Qsi_24P[1][1] + Hammer3D::Qsi_24P[1][2]), -3. + 4. * (Hammer3D::Qsi_24P[1][0] + Hammer3D::Qsi_24P[1][1] + Hammer3D::Qsi_24P[1][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[1][0] + Hammer3D::Qsi_24P[1][1] + Hammer3D::Qsi_24P[1][2]), -4. * Hammer3D::Qsi_24P[1][0], -4. * Hammer3D::Qsi_24P[1][0] },
	  { 4. * Hammer3D::Qsi_24P[1][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[1][1], 4. - 4. * (Hammer3D::Qsi_24P[1][0] + 2. * Hammer3D::Qsi_24P[1][1] + Hammer3D::Qsi_24P[1][2]), -4. * Hammer3D::Qsi_24P[1][1] },
	  { 4. * Hammer3D::Qsi_24P[1][1], 4. * Hammer3D::Qsi_24P[1][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[1][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[1][2], -4. * Hammer3D::Qsi_24P[1][2], 4. - 4. * (Hammer3D::Qsi_24P[1][0] + Hammer3D::Qsi_24P[1][1] + 2. * Hammer3D::Qsi_24P[1][2]) },
	  { 4. * Hammer3D::Qsi_24P[1][2], 0., 4. * Hammer3D::Qsi_24P[1][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[1][2], 4. * Hammer3D::Qsi_24P[1][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[1][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[2][0] + Hammer3D::Qsi_24P[2][1] + Hammer3D::Qsi_24P[2][2]), -3. + 4. * (Hammer3D::Qsi_24P[2][0] + Hammer3D::Qsi_24P[2][1] + Hammer3D::Qsi_24P[2][2]), -3. + 4. * (Hammer3D::Qsi_24P[2][0] + Hammer3D::Qsi_24P[2][1] + Hammer3D::Qsi_24P[2][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[2][0] + Hammer3D::Qsi_24P[2][1] + Hammer3D::Qsi_24P[2][2]), -4. * Hammer3D::Qsi_24P[2][0], -4. * Hammer3D::Qsi_24P[2][0] },
	  { 4. * Hammer3D::Qsi_24P[2][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[2][1], 4. - 4. * (Hammer3D::Qsi_24P[2][0] + 2. * Hammer3D::Qsi_24P[2][1] + Hammer3D::Qsi_24P[2][2]), -4. * Hammer3D::Qsi_24P[2][1] },
	  { 4. * Hammer3D::Qsi_24P[2][1], 4. * Hammer3D::Qsi_24P[2][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[2][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[2][2], -4. * Hammer3D::Qsi_24P[2][2], 4. - 4. * (Hammer3D::Qsi_24P[2][0] + Hammer3D::Qsi_24P[2][1] + 2. * Hammer3D::Qsi_24P[2][2]) },
	  { 4. * Hammer3D::Qsi_24P[2][2], 0., 4. * Hammer3D::Qsi_24P[2][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[2][2], 4. * Hammer3D::Qsi_24P[2][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[2][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[3][0] + Hammer3D::Qsi_24P[3][1] + Hammer3D::Qsi_24P[3][2]), -3. + 4. * (Hammer3D::Qsi_24P[3][0] + Hammer3D::Qsi_24P[3][1] + Hammer3D::Qsi_24P[3][2]), -3. + 4. * (Hammer3D::Qsi_24P[3][0] + Hammer3D::Qsi_24P[3][1] + Hammer3D::Qsi_24P[3][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[3][0] + Hammer3D::Qsi_24P[3][1] + Hammer3D::Qsi_24P[3][2]), -4. * Hammer3D::Qsi_24P[3][0], -4. * Hammer3D::Qsi_24P[3][0] },
	  { 4. * Hammer3D::Qsi_24P[3][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[3][1], 4. - 4. * (Hammer3D::Qsi_24P[3][0] + 2. * Hammer3D::Qsi_24P[3][1] + Hammer3D::Qsi_24P[3][2]), -4. * Hammer3D::Qsi_24P[3][1] },
	  { 4. * Hammer3D::Qsi_24P[3][1], 4. * Hammer3D::Qsi_24P[3][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[3][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[3][2], -4. * Hammer3D::Qsi_24P[3][2], 4. - 4. * (Hammer3D::Qsi_24P[3][0] + Hammer3D::Qsi_24P[3][1] + 2. * Hammer3D::Qsi_24P[3][2]) },
	  { 4. * Hammer3D::Qsi_24P[3][2], 0., 4. * Hammer3D::Qsi_24P[3][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[3][2], 4. * Hammer3D::Qsi_24P[3][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[3][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[4][0] + Hammer3D::Qsi_24P[4][1] + Hammer3D::Qsi_24P[4][2]), -3. + 4. * (Hammer3D::Qsi_24P[4][0] + Hammer3D::Qsi_24P[4][1] + Hammer3D::Qsi_24P[4][2]), -3. + 4. * (Hammer3D::Qsi_24P[4][0] + Hammer3D::Qsi_24P[4][1] + Hammer3D::Qsi_24P[4][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[4][0] + Hammer3D::Qsi_24P[4][1] + Hammer3D::Qsi_24P[4][2]), -4. * Hammer3D::Qsi_24P[4][0], -4. * Hammer3D::Qsi_24P[4][0] },
	  { 4. * Hammer3D::Qsi_24P[4][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[4][1], 4. - 4. * (Hammer3D::Qsi_24P[4][0] + 2. * Hammer3D::Qsi_24P[4][1] + Hammer3D::Qsi_24P[4][2]), -4. * Hammer3D::Qsi_24P[4][1] },
	  { 4. * Hammer3D::Qsi_24P[4][1], 4. * Hammer3D::Qsi_24P[4][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[4][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[4][2], -4. * Hammer3D::Qsi_24P[4][2], 4. - 4. * (Hammer3D::Qsi_24P[4][0] + Hammer3D::Qsi_24P[4][1] + 2. * Hammer3D::Qsi_24P[4][2]) },
	  { 4. * Hammer3D::Qsi_24P[4][2], 0., 4. * Hammer3D::Qsi_24P[4][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[4][2], 4. * Hammer3D::Qsi_24P[4][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[4][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[5][0] + Hammer3D::Qsi_24P[5][1] + Hammer3D::Qsi_24P[5][2]), -3. + 4. * (Hammer3D::Qsi_24P[5][0] + Hammer3D::Qsi_24P[5][1] + Hammer3D::Qsi_24P[5][2]), -3. + 4. * (Hammer3D::Qsi_24P[5][0] + Hammer3D::Qsi_24P[5][1] + Hammer3D::Qsi_24P[5][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[5][0] + Hammer3D::Qsi_24P[5][1] + Hammer3D::Qsi_24P[5][2]), -4. * Hammer3D::Qsi_24P[5][0], -4. * Hammer3D::Qsi_24P[5][0] },
	  { 4. * Hammer3D::Qsi_24P[5][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[5][1], 4. - 4. * (Hammer3D::Qsi_24P[5][0] + 2. * Hammer3D::Qsi_24P[5][1] + Hammer3D::Qsi_24P[5][2]), -4. * Hammer3D::Qsi_24P[5][1] },
	  { 4. * Hammer3D::Qsi_24P[5][1], 4. * Hammer3D::Qsi_24P[5][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[5][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[5][2], -4. * Hammer3D::Qsi_24P[5][2], 4. - 4. * (Hammer3D::Qsi_24P[5][0] + Hammer3D::Qsi_24P[5][1] + 2. * Hammer3D::Qsi_24P[5][2]) },
	  { 4. * Hammer3D::Qsi_24P[5][2], 0., 4. * Hammer3D::Qsi_24P[5][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[5][2], 4. * Hammer3D::Qsi_24P[5][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[5][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[6][0] + Hammer3D::Qsi_24P[6][1] + Hammer3D::Qsi_24P[6][2]), -3. + 4. * (Hammer3D::Qsi_24P[6][0] + Hammer3D::Qsi_24P[6][1] + Hammer3D::Qsi_24P[6][2]), -3. + 4. * (Hammer3D::Qsi_24P[6][0] + Hammer3D::Qsi_24P[6][1] + Hammer3D::Qsi_24P[6][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[6][0] + Hammer3D::Qsi_24P[6][1] + Hammer3D::Qsi_24P[6][2]), -4. * Hammer3D::Qsi_24P[6][0], -4. * Hammer3D::Qsi_24P[6][0] },
	  { 4. * Hammer3D::Qsi_24P[6][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[6][1], 4. - 4. * (Hammer3D::Qsi_24P[6][0] + 2. * Hammer3D::Qsi_24P[6][1] + Hammer3D::Qsi_24P[6][2]), -4. * Hammer3D::Qsi_24P[6][1] },
	  { 4. * Hammer3D::Qsi_24P[6][1], 4. * Hammer3D::Qsi_24P[6][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[6][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[6][2], -4. * Hammer3D::Qsi_24P[6][2], 4. - 4. * (Hammer3D::Qsi_24P[6][0] + Hammer3D::Qsi_24P[6][1] + 2. * Hammer3D::Qsi_24P[6][2]) },
	  { 4. * Hammer3D::Qsi_24P[6][2], 0., 4. * Hammer3D::Qsi_24P[6][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[6][2], 4. * Hammer3D::Qsi_24P[6][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[6][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[7][0] + Hammer3D::Qsi_24P[7][1] + Hammer3D::Qsi_24P[7][2]), -3. + 4. * (Hammer3D::Qsi_24P[7][0] + Hammer3D::Qsi_24P[7][1] + Hammer3D::Qsi_24P[7][2]), -3. + 4. * (Hammer3D::Qsi_24P[7][0] + Hammer3D::Qsi_24P[7][1] + Hammer3D::Qsi_24P[7][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[7][0] + Hammer3D::Qsi_24P[7][1] + Hammer3D::Qsi_24P[7][2]), -4. * Hammer3D::Qsi_24P[7][0], -4. * Hammer3D::Qsi_24P[7][0] },
	  { 4. * Hammer3D::Qsi_24P[7][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[7][1], 4. - 4. * (Hammer3D::Qsi_24P[7][0] + 2. * Hammer3D::Qsi_24P[7][1] + Hammer3D::Qsi_24P[7][2]), -4. * Hammer3D::Qsi_24P[7][1] },
	  { 4. * Hammer3D::Qsi_24P[7][1], 4. * Hammer3D::Qsi_24P[7][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[7][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[7][2], -4. * Hammer3D::Qsi_24P[7][2], 4. - 4. * (Hammer3D::Qsi_24P[7][0] + Hammer3D::Qsi_24P[7][1] + 2. * Hammer3D::Qsi_24P[7][2]) },
	  { 4. * Hammer3D::Qsi_24P[7][2], 0., 4. * Hammer3D::Qsi_24P[7][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[7][2], 4. * Hammer3D::Qsi_24P[7][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[7][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[8][0] + Hammer3D::Qsi_24P[8][1] + Hammer3D::Qsi_24P[8][2]), -3. + 4. * (Hammer3D::Qsi_24P[8][0] + Hammer3D::Qsi_24P[8][1] + Hammer3D::Qsi_24P[8][2]), -3. + 4. * (Hammer3D::Qsi_24P[8][0] + Hammer3D::Qsi_24P[8][1] + Hammer3D::Qsi_24P[8][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[8][0] + Hammer3D::Qsi_24P[8][1] + Hammer3D::Qsi_24P[8][2]), -4. * Hammer3D::Qsi_24P[8][0], -4. * Hammer3D::Qsi_24P[8][0] },
	  { 4. * Hammer3D::Qsi_24P[8][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[8][1], 4. - 4. * (Hammer3D::Qsi_24P[8][0] + 2. * Hammer3D::Qsi_24P[8][1] + Hammer3D::Qsi_24P[8][2]), -4. * Hammer3D::Qsi_24P[8][1] },
	  { 4. * Hammer3D::Qsi_24P[8][1], 4. * Hammer3D::Qsi_24P[8][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[8][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[8][2], -4. * Hammer3D::Qsi_24P[8][2], 4. - 4. * (Hammer3D::Qsi_24P[8][0] + Hammer3D::Qsi_24P[8][1] + 2. * Hammer3D::Qsi_24P[8][2]) },
	  { 4. * Hammer3D::Qsi_24P[8][2], 0., 4. * Hammer3D::Qsi_24P[8][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[8][2], 4. * Hammer3D::Qsi_24P[8][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[8][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[9][0] + Hammer3D::Qsi_24P[9][1] + Hammer3D::Qsi_24P[9][2]), -3. + 4. * (Hammer3D::Qsi_24P[9][0] + Hammer3D::Qsi_24P[9][1] + Hammer3D::Qsi_24P[9][2]), -3. + 4. * (Hammer3D::Qsi_24P[9][0] + Hammer3D::Qsi_24P[9][1] + Hammer3D::Qsi_24P[9][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[9][0] + Hammer3D::Qsi_24P[9][1] + Hammer3D::Qsi_24P[9][2]), -4. * Hammer3D::Qsi_24P[9][0], -4. * Hammer3D::Qsi_24P[9][0] },
	  { 4. * Hammer3D::Qsi_24P[9][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[9][1], 4. - 4. * (Hammer3D::Qsi_24P[9][0] + 2. * Hammer3D::Qsi_24P[9][1] + Hammer3D::Qsi_24P[9][2]), -4. * Hammer3D::Qsi_24P[9][1] },
	  { 4. * Hammer3D::Qsi_24P[9][1], 4. * Hammer3D::Qsi_24P[9][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[9][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[9][2], -4. * Hammer3D::Qsi_24P[9][2], 4. - 4. * (Hammer3D::Qsi_24P[9][0] + Hammer3D::Qsi_24P[9][1] + 2. * Hammer3D::Qsi_24P[9][2]) },
	  { 4. * Hammer3D::Qsi_24P[9][2], 0., 4. * Hammer3D::Qsi_24P[9][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[9][2], 4. * Hammer3D::Qsi_24P[9][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[9][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[10][0] + Hammer3D::Qsi_24P[10][1] + Hammer3D::Qsi_24P[10][2]), -3. + 4. * (Hammer3D::Qsi_24P[10][0] + Hammer3D::Qsi_24P[10][1] + Hammer3D::Qsi_24P[10][2]), -3. + 4. * (Hammer3D::Qsi_24P[10][0] + Hammer3D::Qsi_24P[10][1] + Hammer3D::Qsi_24P[10][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[10][0] + Hammer3D::Qsi_24P[10][1] + Hammer3D::Qsi_24P[10][2]), -4. * Hammer3D::Qsi_24P[10][0], -4. * Hammer3D::Qsi_24P[10][0] },
	  { 4. * Hammer3D::Qsi_24P[10][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[10][1], 4. - 4. * (Hammer3D::Qsi_24P[10][0] + 2. * Hammer3D::Qsi_24P[10][1] + Hammer3D::Qsi_24P[10][2]), -4. * Hammer3D::Qsi_24P[10][1] },
	  { 4. * Hammer3D::Qsi_24P[10][1], 4. * Hammer3D::Qsi_24P[10][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[10][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[10][2], -4. * Hammer3D::Qsi_24P[10][2], 4. - 4. * (Hammer3D::Qsi_24P[10][0] + Hammer3D::Qsi_24P[10][1] + 2. * Hammer3D::Qsi_24P[10][2]) },
	  { 4. * Hammer3D::Qsi_24P[10][2], 0., 4. * Hammer3D::Qsi_24P[10][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[10][2], 4. * Hammer3D::Qsi_24P[10][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[10][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[11][0] + Hammer3D::Qsi_24P[11][1] + Hammer3D::Qsi_24P[11][2]), -3. + 4. * (Hammer3D::Qsi_24P[11][0] + Hammer3D::Qsi_24P[11][1] + Hammer3D::Qsi_24P[11][2]), -3. + 4. * (Hammer3D::Qsi_24P[11][0] + Hammer3D::Qsi_24P[11][1] + Hammer3D::Qsi_24P[11][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[11][0] + Hammer3D::Qsi_24P[11][1] + Hammer3D::Qsi_24P[11][2]), -4. * Hammer3D::Qsi_24P[11][0], -4. * Hammer3D::Qsi_24P[11][0] },
	  { 4. * Hammer3D::Qsi_24P[11][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[11][1], 4. - 4. * (Hammer3D::Qsi_24P[11][0] + 2. * Hammer3D::Qsi_24P[11][1] + Hammer3D::Qsi_24P[11][2]), -4. * Hammer3D::Qsi_24P[11][1] },
	  { 4. * Hammer3D::Qsi_24P[11][1], 4. * Hammer3D::Qsi_24P[11][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[11][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[11][2], -4. * Hammer3D::Qsi_24P[11][2], 4. - 4. * (Hammer3D::Qsi_24P[11][0] + Hammer3D::Qsi_24P[11][1] + 2. * Hammer3D::Qsi_24P[11][2]) },
	  { 4. * Hammer3D::Qsi_24P[11][2], 0., 4. * Hammer3D::Qsi_24P[11][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[11][2], 4. * Hammer3D::Qsi_24P[11][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[11][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[12][0] + Hammer3D::Qsi_24P[12][1] + Hammer3D::Qsi_24P[12][2]), -3. + 4. * (Hammer3D::Qsi_24P[12][0] + Hammer3D::Qsi_24P[12][1] + Hammer3D::Qsi_24P[12][2]), -3. + 4. * (Hammer3D::Qsi_24P[12][0] + Hammer3D::Qsi_24P[12][1] + Hammer3D::Qsi_24P[12][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[12][0] + Hammer3D::Qsi_24P[12][1] + Hammer3D::Qsi_24P[12][2]), -4. * Hammer3D::Qsi_24P[12][0], -4. * Hammer3D::Qsi_24P[12][0] },
	  { 4. * Hammer3D::Qsi_24P[12][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[12][1], 4. - 4. * (Hammer3D::Qsi_24P[12][0] + 2. * Hammer3D::Qsi_24P[12][1] + Hammer3D::Qsi_24P[12][2]), -4. * Hammer3D::Qsi_24P[12][1] },
	  { 4. * Hammer3D::Qsi_24P[12][1], 4. * Hammer3D::Qsi_24P[12][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[12][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[12][2], -4. * Hammer3D::Qsi_24P[12][2], 4. - 4. * (Hammer3D::Qsi_24P[12][0] + Hammer3D::Qsi_24P[12][1] + 2. * Hammer3D::Qsi_24P[12][2]) },
	  { 4. * Hammer3D::Qsi_24P[12][2], 0., 4. * Hammer3D::Qsi_24P[12][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[12][2], 4. * Hammer3D::Qsi_24P[12][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[12][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[13][0] + Hammer3D::Qsi_24P[13][1] + Hammer3D::Qsi_24P[13][2]), -3. + 4. * (Hammer3D::Qsi_24P[13][0] + Hammer3D::Qsi_24P[13][1] + Hammer3D::Qsi_24P[13][2]), -3. + 4. * (Hammer3D::Qsi_24P[13][0] + Hammer3D::Qsi_24P[13][1] + Hammer3D::Qsi_24P[13][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[13][0] + Hammer3D::Qsi_24P[13][1] + Hammer3D::Qsi_24P[13][2]), -4. * Hammer3D::Qsi_24P[13][0], -4. * Hammer3D::Qsi_24P[13][0] },
	  { 4. * Hammer3D::Qsi_24P[13][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[13][1], 4. - 4. * (Hammer3D::Qsi_24P[13][0] + 2. * Hammer3D::Qsi_24P[13][1] + Hammer3D::Qsi_24P[13][2]), -4. * Hammer3D::Qsi_24P[13][1] },
	  { 4. * Hammer3D::Qsi_24P[13][1], 4. * Hammer3D::Qsi_24P[13][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[13][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[13][2], -4. * Hammer3D::Qsi_24P[13][2], 4. - 4. * (Hammer3D::Qsi_24P[13][0] + Hammer3D::Qsi_24P[13][1] + 2. * Hammer3D::Qsi_24P[13][2]) },
	  { 4. * Hammer3D::Qsi_24P[13][2], 0., 4. * Hammer3D::Qsi_24P[13][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[13][2], 4. * Hammer3D::Qsi_24P[13][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[13][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[14][0] + Hammer3D::Qsi_24P[14][1] + Hammer3D::Qsi_24P[14][2]), -3. + 4. * (Hammer3D::Qsi_24P[14][0] + Hammer3D::Qsi_24P[14][1] + Hammer3D::Qsi_24P[14][2]), -3. + 4. * (Hammer3D::Qsi_24P[14][0] + Hammer3D::Qsi_24P[14][1] + Hammer3D::Qsi_24P[14][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[14][0] + Hammer3D::Qsi_24P[14][1] + Hammer3D::Qsi_24P[14][2]), -4. * Hammer3D::Qsi_24P[14][0], -4. * Hammer3D::Qsi_24P[14][0] },
	  { 4. * Hammer3D::Qsi_24P[14][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[14][1], 4. - 4. * (Hammer3D::Qsi_24P[14][0] + 2. * Hammer3D::Qsi_24P[14][1] + Hammer3D::Qsi_24P[14][2]), -4. * Hammer3D::Qsi_24P[14][1] },
	  { 4. * Hammer3D::Qsi_24P[14][1], 4. * Hammer3D::Qsi_24P[14][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[14][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[14][2], -4. * Hammer3D::Qsi_24P[14][2], 4. - 4. * (Hammer3D::Qsi_24P[14][0] + Hammer3D::Qsi_24P[14][1] + 2. * Hammer3D::Qsi_24P[14][2]) },
	  { 4. * Hammer3D::Qsi_24P[14][2], 0., 4. * Hammer3D::Qsi_24P[14][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[14][2], 4. * Hammer3D::Qsi_24P[14][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[14][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[15][0] + Hammer3D::Qsi_24P[15][1] + Hammer3D::Qsi_24P[15][2]), -3. + 4. * (Hammer3D::Qsi_24P[15][0] + Hammer3D::Qsi_24P[15][1] + Hammer3D::Qsi_24P[15][2]), -3. + 4. * (Hammer3D::Qsi_24P[15][0] + Hammer3D::Qsi_24P[15][1] + Hammer3D::Qsi_24P[15][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[15][0] + Hammer3D::Qsi_24P[15][1] + Hammer3D::Qsi_24P[15][2]), -4. * Hammer3D::Qsi_24P[15][0], -4. * Hammer3D::Qsi_24P[15][0] },
	  { 4. * Hammer3D::Qsi_24P[15][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[15][1], 4. - 4. * (Hammer3D::Qsi_24P[15][0] + 2. * Hammer3D::Qsi_24P[15][1] + Hammer3D::Qsi_24P[15][2]), -4. * Hammer3D::Qsi_24P[15][1] },
	  { 4. * Hammer3D::Qsi_24P[15][1], 4. * Hammer3D::Qsi_24P[15][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[15][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[15][2], -4. * Hammer3D::Qsi_24P[15][2], 4. - 4. * (Hammer3D::Qsi_24P[15][0] + Hammer3D::Qsi_24P[15][1] + 2. * Hammer3D::Qsi_24P[15][2]) },
	  { 4. * Hammer3D::Qsi_24P[15][2], 0., 4. * Hammer3D::Qsi_24P[15][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[15][2], 4. * Hammer3D::Qsi_24P[15][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[15][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[16][0] + Hammer3D::Qsi_24P[16][1] + Hammer3D::Qsi_24P[16][2]), -3. + 4. * (Hammer3D::Qsi_24P[16][0] + Hammer3D::Qsi_24P[16][1] + Hammer3D::Qsi_24P[16][2]), -3. + 4. * (Hammer3D::Qsi_24P[16][0] + Hammer3D::Qsi_24P[16][1] + Hammer3D::Qsi_24P[16][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[16][0] + Hammer3D::Qsi_24P[16][1] + Hammer3D::Qsi_24P[16][2]), -4. * Hammer3D::Qsi_24P[16][0], -4. * Hammer3D::Qsi_24P[16][0] },
	  { 4. * Hammer3D::Qsi_24P[16][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[16][1], 4. - 4. * (Hammer3D::Qsi_24P[16][0] + 2. * Hammer3D::Qsi_24P[16][1] + Hammer3D::Qsi_24P[16][2]), -4. * Hammer3D::Qsi_24P[16][1] },
	  { 4. * Hammer3D::Qsi_24P[16][1], 4. * Hammer3D::Qsi_24P[16][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[16][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[16][2], -4. * Hammer3D::Qsi_24P[16][2], 4. - 4. * (Hammer3D::Qsi_24P[16][0] + Hammer3D::Qsi_24P[16][1] + 2. * Hammer3D::Qsi_24P[16][2]) },
	  { 4. * Hammer3D::Qsi_24P[16][2], 0., 4. * Hammer3D::Qsi_24P[16][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[16][2], 4. * Hammer3D::Qsi_24P[16][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[16][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[17][0] + Hammer3D::Qsi_24P[17][1] + Hammer3D::Qsi_24P[17][2]), -3. + 4. * (Hammer3D::Qsi_24P[17][0] + Hammer3D::Qsi_24P[17][1] + Hammer3D::Qsi_24P[17][2]), -3. + 4. * (Hammer3D::Qsi_24P[17][0] + Hammer3D::Qsi_24P[17][1] + Hammer3D::Qsi_24P[17][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[17][0] + Hammer3D::Qsi_24P[17][1] + Hammer3D::Qsi_24P[17][2]), -4. * Hammer3D::Qsi_24P[17][0], -4. * Hammer3D::Qsi_24P[17][0] },
	  { 4. * Hammer3D::Qsi_24P[17][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[17][1], 4. - 4. * (Hammer3D::Qsi_24P[17][0] + 2. * Hammer3D::Qsi_24P[17][1] + Hammer3D::Qsi_24P[17][2]), -4. * Hammer3D::Qsi_24P[17][1] },
	  { 4. * Hammer3D::Qsi_24P[17][1], 4. * Hammer3D::Qsi_24P[17][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[17][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[17][2], -4. * Hammer3D::Qsi_24P[17][2], 4. - 4. * (Hammer3D::Qsi_24P[17][0] + Hammer3D::Qsi_24P[17][1] + 2. * Hammer3D::Qsi_24P[17][2]) },
	  { 4. * Hammer3D::Qsi_24P[17][2], 0., 4. * Hammer3D::Qsi_24P[17][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[17][2], 4. * Hammer3D::Qsi_24P[17][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[17][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[18][0] + Hammer3D::Qsi_24P[18][1] + Hammer3D::Qsi_24P[18][2]), -3. + 4. * (Hammer3D::Qsi_24P[18][0] + Hammer3D::Qsi_24P[18][1] + Hammer3D::Qsi_24P[18][2]), -3. + 4. * (Hammer3D::Qsi_24P[18][0] + Hammer3D::Qsi_24P[18][1] + Hammer3D::Qsi_24P[18][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[18][0] + Hammer3D::Qsi_24P[18][1] + Hammer3D::Qsi_24P[18][2]), -4. * Hammer3D::Qsi_24P[18][0], -4. * Hammer3D::Qsi_24P[18][0] },
	  { 4. * Hammer3D::Qsi_24P[18][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[18][1], 4. - 4. * (Hammer3D::Qsi_24P[18][0] + 2. * Hammer3D::Qsi_24P[18][1] + Hammer3D::Qsi_24P[18][2]), -4. * Hammer3D::Qsi_24P[18][1] },
	  { 4. * Hammer3D::Qsi_24P[18][1], 4. * Hammer3D::Qsi_24P[18][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[18][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[18][2], -4. * Hammer3D::Qsi_24P[18][2], 4. - 4. * (Hammer3D::Qsi_24P[18][0] + Hammer3D::Qsi_24P[18][1] + 2. * Hammer3D::Qsi_24P[18][2]) },
	  { 4. * Hammer3D::Qsi_24P[18][2], 0., 4. * Hammer3D::Qsi_24P[18][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[18][2], 4. * Hammer3D::Qsi_24P[18][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[18][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[19][0] + Hammer3D::Qsi_24P[19][1] + Hammer3D::Qsi_24P[19][2]), -3. + 4. * (Hammer3D::Qsi_24P[19][0] + Hammer3D::Qsi_24P[19][1] + Hammer3D::Qsi_24P[19][2]), -3. + 4. * (Hammer3D::Qsi_24P[19][0] + Hammer3D::Qsi_24P[19][1] + Hammer3D::Qsi_24P[19][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[19][0] + Hammer3D::Qsi_24P[19][1] + Hammer3D::Qsi_24P[19][2]), -4. * Hammer3D::Qsi_24P[19][0], -4. * Hammer3D::Qsi_24P[19][0] },
	  { 4. * Hammer3D::Qsi_24P[19][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[19][1], 4. - 4. * (Hammer3D::Qsi_24P[19][0] + 2. * Hammer3D::Qsi_24P[19][1] + Hammer3D::Qsi_24P[19][2]), -4. * Hammer3D::Qsi_24P[19][1] },
	  { 4. * Hammer3D::Qsi_24P[19][1], 4. * Hammer3D::Qsi_24P[19][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[19][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[19][2], -4. * Hammer3D::Qsi_24P[19][2], 4. - 4. * (Hammer3D::Qsi_24P[19][0] + Hammer3D::Qsi_24P[19][1] + 2. * Hammer3D::Qsi_24P[19][2]) },
	  { 4. * Hammer3D::Qsi_24P[19][2], 0., 4. * Hammer3D::Qsi_24P[19][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[19][2], 4. * Hammer3D::Qsi_24P[19][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[19][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[20][0] + Hammer3D::Qsi_24P[20][1] + Hammer3D::Qsi_24P[20][2]), -3. + 4. * (Hammer3D::Qsi_24P[20][0] + Hammer3D::Qsi_24P[20][1] + Hammer3D::Qsi_24P[20][2]), -3. + 4. * (Hammer3D::Qsi_24P[20][0] + Hammer3D::Qsi_24P[20][1] + Hammer3D::Qsi_24P[20][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[20][0] + Hammer3D::Qsi_24P[20][1] + Hammer3D::Qsi_24P[20][2]), -4. * Hammer3D::Qsi_24P[20][0], -4. * Hammer3D::Qsi_24P[20][0] },
	  { 4. * Hammer3D::Qsi_24P[20][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[20][1], 4. - 4. * (Hammer3D::Qsi_24P[20][0] + 2. * Hammer3D::Qsi_24P[20][1] + Hammer3D::Qsi_24P[20][2]), -4. * Hammer3D::Qsi_24P[20][1] },
	  { 4. * Hammer3D::Qsi_24P[20][1], 4. * Hammer3D::Qsi_24P[20][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[20][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[20][2], -4. * Hammer3D::Qsi_24P[20][2], 4. - 4. * (Hammer3D::Qsi_24P[20][0] + Hammer3D::Qsi_24P[20][1] + 2. * Hammer3D::Qsi_24P[20][2]) },
	  { 4. * Hammer3D::Qsi_24P[20][2], 0., 4. * Hammer3D::Qsi_24P[20][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[20][2], 4. * Hammer3D::Qsi_24P[20][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[20][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[21][0] + Hammer3D::Qsi_24P[21][1] + Hammer3D::Qsi_24P[21][2]), -3. + 4. * (Hammer3D::Qsi_24P[21][0] + Hammer3D::Qsi_24P[21][1] + Hammer3D::Qsi_24P[21][2]), -3. + 4. * (Hammer3D::Qsi_24P[21][0] + Hammer3D::Qsi_24P[21][1] + Hammer3D::Qsi_24P[21][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[21][0] + Hammer3D::Qsi_24P[21][1] + Hammer3D::Qsi_24P[21][2]), -4. * Hammer3D::Qsi_24P[21][0], -4. * Hammer3D::Qsi_24P[21][0] },
	  { 4. * Hammer3D::Qsi_24P[21][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[21][1], 4. - 4. * (Hammer3D::Qsi_24P[21][0] + 2. * Hammer3D::Qsi_24P[21][1] + Hammer3D::Qsi_24P[21][2]), -4. * Hammer3D::Qsi_24P[21][1] },
	  { 4. * Hammer3D::Qsi_24P[21][1], 4. * Hammer3D::Qsi_24P[21][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[21][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[21][2], -4. * Hammer3D::Qsi_24P[21][2], 4. - 4. * (Hammer3D::Qsi_24P[21][0] + Hammer3D::Qsi_24P[21][1] + 2. * Hammer3D::Qsi_24P[21][2]) },
	  { 4. * Hammer3D::Qsi_24P[21][2], 0., 4. * Hammer3D::Qsi_24P[21][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[21][2], 4. * Hammer3D::Qsi_24P[21][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[21][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[22][0] + Hammer3D::Qsi_24P[22][1] + Hammer3D::Qsi_24P[22][2]), -3. + 4. * (Hammer3D::Qsi_24P[22][0] + Hammer3D::Qsi_24P[22][1] + Hammer3D::Qsi_24P[22][2]), -3. + 4. * (Hammer3D::Qsi_24P[22][0] + Hammer3D::Qsi_24P[22][1] + Hammer3D::Qsi_24P[22][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[22][0] + Hammer3D::Qsi_24P[22][1] + Hammer3D::Qsi_24P[22][2]), -4. * Hammer3D::Qsi_24P[22][0], -4. * Hammer3D::Qsi_24P[22][0] },
	  { 4. * Hammer3D::Qsi_24P[22][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[22][1], 4. - 4. * (Hammer3D::Qsi_24P[22][0] + 2. * Hammer3D::Qsi_24P[22][1] + Hammer3D::Qsi_24P[22][2]), -4. * Hammer3D::Qsi_24P[22][1] },
	  { 4. * Hammer3D::Qsi_24P[22][1], 4. * Hammer3D::Qsi_24P[22][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[22][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[22][2], -4. * Hammer3D::Qsi_24P[22][2], 4. - 4. * (Hammer3D::Qsi_24P[22][0] + Hammer3D::Qsi_24P[22][1] + 2. * Hammer3D::Qsi_24P[22][2]) },
	  { 4. * Hammer3D::Qsi_24P[22][2], 0., 4. * Hammer3D::Qsi_24P[22][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[22][2], 4. * Hammer3D::Qsi_24P[22][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[22][2] - 1.} },

	{ { -3. + 4. * (Hammer3D::Qsi_24P[23][0] + Hammer3D::Qsi_24P[23][1] + Hammer3D::Qsi_24P[23][2]), -3. + 4. * (Hammer3D::Qsi_24P[23][0] + Hammer3D::Qsi_24P[23][1] + Hammer3D::Qsi_24P[23][2]), -3. + 4. * (Hammer3D::Qsi_24P[23][0] + Hammer3D::Qsi_24P[23][1] + Hammer3D::Qsi_24P[23][2]) },
	  { 4. - 4. * (2. * Hammer3D::Qsi_24P[23][0] + Hammer3D::Qsi_24P[23][1] + Hammer3D::Qsi_24P[23][2]), -4. * Hammer3D::Qsi_24P[23][0], -4. * Hammer3D::Qsi_24P[23][0] },
	  { 4. * Hammer3D::Qsi_24P[23][0] - 1., 0., 0. },
	  { -4. * Hammer3D::Qsi_24P[23][1], 4. - 4. * (Hammer3D::Qsi_24P[23][0] + 2. * Hammer3D::Qsi_24P[23][1] + Hammer3D::Qsi_24P[23][2]), -4. * Hammer3D::Qsi_24P[23][1] },
	  { 4. * Hammer3D::Qsi_24P[23][1], 4. * Hammer3D::Qsi_24P[23][0], 0. },
	  { 0., 4. * Hammer3D::Qsi_24P[23][1] - 1., 0. },
	  { -4. * Hammer3D::Qsi_24P[23][2], -4. * Hammer3D::Qsi_24P[23][2], 4. - 4. * (Hammer3D::Qsi_24P[23][0] + Hammer3D::Qsi_24P[23][1] + 2. * Hammer3D::Qsi_24P[23][2]) },
	  { 4. * Hammer3D::Qsi_24P[23][2], 0., 4. * Hammer3D::Qsi_24P[23][0] },
	  { 0., 4. * Hammer3D::Qsi_24P[23][2], 4. * Hammer3D::Qsi_24P[23][1] },
	  { 0., 0., 4. * Hammer3D::Qsi_24P[23][2] - 1.} } };
