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
// Plane element, with linear interpolation functions, rectangular shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"
#include "IntegrationPoint.h"

namespace O2P2 {
	namespace Prep {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Rect4
			  *
			  * @brief Quadrangular linear element with 4 nodes.
			  * @details Plane element, with linear interpolation functions, rectangular shaped.
			  * Options for integration points: 4, 9 and 16.
			  * @image html Elem_Quad4.png height=300
			  */
			class Elem_Rect4 : public ElementPlane
			{
			private:
				Elem_Rect4() = delete;

			protected:
				/** Constructor for rectangular linear elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				Elem_Rect4(std::shared_ptr<O2P2::Prep::Material>& Material, std::shared_ptr<O2P2::Prep::Section>& Section)
					: ElementPlane(Material, Section) { }

			public:
				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "3 1 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " "
						<< this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "3 1 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (4 + add)
						<< " " << this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Evaluates shape function in the point.
				std::vector<double> getShapeFcOnPoint(const double* Point) override;

				// Evaluates the derivative of shape function in the point.
				std::vector<double> getShapeDerivOnPoint(const double* Point) override;

				// Returns the number of nodes of current element.
				int getNumNodes() override { return mv_numNodes; }

				// Returns the number of faces of current element.
				int getNumFaces() override { return mv_numFaces; }

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates.
				  */
				inline bool evaluateXsi(const std::array<double, mv_Dim> xsi) override {
					const auto [min, max] = std::minmax_element(xsi.begin(), xsi.end());
					if (*max < 1.000001 && *min > -1.000001) return true;
					return false;
				}

			protected:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			protected:
				/** @brief Number of Nodes */
				static const int mv_numNodes{ 4 };

				/** @brief Number of Faces */
				static const int mv_numFaces{ 1 };
			};

			/**
			  * @class Elem_Rect4_IP
			  *
			  * @brief Quadrangular linear element with 4 nodes.
			  * @details Plane element, with linear interpolation functions, rectangular shaped.
			  *
			  * @tparam nIP Number of integration points. Must be: 4, 9 or 16.
			  */
			template<int nIP>
			class Elem_Rect4_IP : public Elem_Rect4
			{
			private:
				Elem_Rect4_IP() = delete;

			public:
				/** Constructor for rectangular linear elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				explicit Elem_Rect4_IP(std::shared_ptr<O2P2::Prep::Material>& Material, std::shared_ptr<O2P2::Prep::Section>& Section)
					: Elem_Rect4(Material, Section) { }

				// Return a vector with values on the integration points currently known in the element' nodes.
				std::vector<double> getValueOnIPs(const double* value) override;

				// Returns a pointer to the first element of the shape functions (with size [nIP][mv_numNodes]).
				double const* getShapeFc() const override { return &mv_Psi[0][0]; }

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][mv_numNodes][mv_Dim]).
				double const* getShapeDerivative() const override { return &mv_DPsi[0][0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return mv_weight; }

				// Returns the number of integration points of current element.
				int getNumIP() override { return nIP; }

			private:
				// Weights for numerical integration
				static const double* mv_weight;

				// Shape functions
				static const double mv_Psi[nIP][mv_numNodes];

				// Shape functions derivative
				static const double mv_DPsi[nIP][mv_numNodes][mv_Dim];
			};
		} // End of Elem Namespace
	} // End of Prep Namespace
} // End of O2P2 Namespace


// ================================================================================================
//
// Implementation of Member Function: getShapeFcOnPoint
// Shape functions evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Prep::Elem::Elem_Rect4::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = 0.25 * (1. - Point[0]) * (1. - Point[1]);
	mi_Psi.at(1) = 0.25 * (1. + Point[0]) * (1. - Point[1]);
	mi_Psi.at(2) = 0.25 * (1. - Point[0]) * (1. + Point[1]);
	mi_Psi.at(3) = 0.25 * (1. + Point[0]) * (1. + Point[1]);

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Prep::Elem::Elem_Rect4::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 2);

	mi_DPsi.at(0) = -0.25 * (1. - Point[1]);
	mi_DPsi.at(1) =  0.25 * (1. - Point[1]);
	mi_DPsi.at(2) = -0.25 * (1. + Point[1]);
	mi_DPsi.at(3) =  0.25 * (1. + Point[1]);

	mi_DPsi.at(4) = -0.25 * (1. - Point[0]);
	mi_DPsi.at(5) = -0.25 * (1. + Point[0]);
	mi_DPsi.at(6) =  0.25 * (1. - Point[0]);
	mi_DPsi.at(7) =  0.25 * (1. + Point[0]);

	return mi_DPsi;
};

// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Prep::Elem::Elem_Rect4::setGeomProperties() {

	const int nVertices = 4;

	// Allocate an array with size mv_Dim to which mv_Centroid points to.
	mv_Centroid = std::make_unique<double[]>(mv_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Prep::Node<mv_Dim>*, nVertices> vertices;
	vertices[0] = mv_Conect[0].get();
	vertices[1] = mv_Conect[1].get();
	vertices[2] = mv_Conect[2].get();
	vertices[3] = mv_Conect[3].get();

	// Memory requested by make_unique is not empty
	for (int i = 0; i < mv_Dim; i++) mv_Centroid[i] = 0.;

	for (auto& node : vertices) {
		std::array<double, mv_Dim> x = node->getInitPos();

		for (int i = 0; i < mv_Dim; i++) mv_Centroid[i] += x[i];
	}

	// Finishing up
	for (int i = 0; i < mv_Dim; i++) mv_Centroid[i] /= nVertices;

	// Distance from centroid to vertices
	double dist[nVertices] = {};
	int i = 0;

	for (auto& node : vertices) {
		std::array<double, mv_Dim> x = node->getInitPos();

		for (int j = 0; j < mv_Dim; j++) {
			dist[i] += (mv_Centroid[j] - x[j]) * (mv_Centroid[j] - x[j]);
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
inline std::vector<double> O2P2::Prep::Elem::Elem_Rect4_IP<nIP>::getValueOnIPs(const double* value) {

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
template<> const double* O2P2::Prep::Elem::Elem_Rect4_IP<4>::mv_weight = &Gauss2D::Wg_4P[0];
template<> const double* O2P2::Prep::Elem::Elem_Rect4_IP<9>::mv_weight = &Gauss2D::Wg_9P[0];
template<> const double* O2P2::Prep::Elem::Elem_Rect4_IP<16>::mv_weight = &Gauss2D::Wg_16P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double O2P2::Prep::Elem::Elem_Rect4_IP<4>::mv_Psi[4][mv_numNodes] = {
	{ 0.25 * (1. - Gauss2D::Qsi_4P[0][0]) * (1. - Gauss2D::Qsi_4P[0][1]), 0.25 * (1. + Gauss2D::Qsi_4P[0][0]) * (1. - Gauss2D::Qsi_4P[0][1]),
	  0.25 * (1. - Gauss2D::Qsi_4P[0][0]) * (1. + Gauss2D::Qsi_4P[0][1]), 0.25 * (1. + Gauss2D::Qsi_4P[0][0]) * (1. + Gauss2D::Qsi_4P[0][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_4P[1][0]) * (1. - Gauss2D::Qsi_4P[1][1]), 0.25 * (1. + Gauss2D::Qsi_4P[1][0]) * (1. - Gauss2D::Qsi_4P[1][1]),
	  0.25 * (1. - Gauss2D::Qsi_4P[1][0]) * (1. + Gauss2D::Qsi_4P[1][1]), 0.25 * (1. + Gauss2D::Qsi_4P[1][0]) * (1. + Gauss2D::Qsi_4P[1][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_4P[2][0]) * (1. - Gauss2D::Qsi_4P[2][1]), 0.25 * (1. + Gauss2D::Qsi_4P[2][0]) * (1. - Gauss2D::Qsi_4P[2][1]),
	  0.25 * (1. - Gauss2D::Qsi_4P[2][0]) * (1. + Gauss2D::Qsi_4P[2][1]), 0.25 * (1. + Gauss2D::Qsi_4P[2][0]) * (1. + Gauss2D::Qsi_4P[2][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_4P[3][0]) * (1. - Gauss2D::Qsi_4P[3][1]), 0.25 * (1. + Gauss2D::Qsi_4P[3][0]) * (1. - Gauss2D::Qsi_4P[3][1]),
	  0.25 * (1. - Gauss2D::Qsi_4P[3][0]) * (1. + Gauss2D::Qsi_4P[3][1]), 0.25 * (1. + Gauss2D::Qsi_4P[3][0]) * (1. + Gauss2D::Qsi_4P[3][1]) } };

template<> const double O2P2::Prep::Elem::Elem_Rect4_IP<9>::mv_Psi[9][mv_numNodes] = {
	{ 0.25 * (1. - Gauss2D::Qsi_9P[0][0]) * (1. - Gauss2D::Qsi_9P[0][1]), 0.25 * (1. + Gauss2D::Qsi_9P[0][0]) * (1. - Gauss2D::Qsi_9P[0][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[0][0]) * (1. + Gauss2D::Qsi_9P[0][1]), 0.25 * (1. + Gauss2D::Qsi_9P[0][0]) * (1. + Gauss2D::Qsi_9P[0][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[1][0]) * (1. - Gauss2D::Qsi_9P[1][1]), 0.25 * (1. + Gauss2D::Qsi_9P[1][0]) * (1. - Gauss2D::Qsi_9P[1][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[1][0]) * (1. + Gauss2D::Qsi_9P[1][1]), 0.25 * (1. + Gauss2D::Qsi_9P[1][0]) * (1. + Gauss2D::Qsi_9P[1][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[2][0]) * (1. - Gauss2D::Qsi_9P[2][1]), 0.25 * (1. + Gauss2D::Qsi_9P[2][0]) * (1. - Gauss2D::Qsi_9P[2][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[2][0]) * (1. + Gauss2D::Qsi_9P[2][1]), 0.25 * (1. + Gauss2D::Qsi_9P[2][0]) * (1. + Gauss2D::Qsi_9P[2][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[3][0]) * (1. - Gauss2D::Qsi_9P[3][1]), 0.25 * (1. + Gauss2D::Qsi_9P[3][0]) * (1. - Gauss2D::Qsi_9P[3][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[3][0]) * (1. + Gauss2D::Qsi_9P[3][1]), 0.25 * (1. + Gauss2D::Qsi_9P[3][0]) * (1. + Gauss2D::Qsi_9P[3][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[4][0]) * (1. - Gauss2D::Qsi_9P[4][1]), 0.25 * (1. + Gauss2D::Qsi_9P[4][0]) * (1. - Gauss2D::Qsi_9P[4][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[4][0]) * (1. + Gauss2D::Qsi_9P[4][1]), 0.25 * (1. + Gauss2D::Qsi_9P[4][0]) * (1. + Gauss2D::Qsi_9P[4][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[5][0]) * (1. - Gauss2D::Qsi_9P[5][1]), 0.25 * (1. + Gauss2D::Qsi_9P[5][0]) * (1. - Gauss2D::Qsi_9P[5][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[5][0]) * (1. + Gauss2D::Qsi_9P[5][1]), 0.25 * (1. + Gauss2D::Qsi_9P[5][0]) * (1. + Gauss2D::Qsi_9P[5][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[6][0]) * (1. - Gauss2D::Qsi_9P[6][1]), 0.25 * (1. + Gauss2D::Qsi_9P[6][0]) * (1. - Gauss2D::Qsi_9P[6][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[6][0]) * (1. + Gauss2D::Qsi_9P[6][1]), 0.25 * (1. + Gauss2D::Qsi_9P[6][0]) * (1. + Gauss2D::Qsi_9P[6][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[7][0]) * (1. - Gauss2D::Qsi_9P[7][1]), 0.25 * (1. + Gauss2D::Qsi_9P[7][0]) * (1. - Gauss2D::Qsi_9P[7][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[7][0]) * (1. + Gauss2D::Qsi_9P[7][1]), 0.25 * (1. + Gauss2D::Qsi_9P[7][0]) * (1. + Gauss2D::Qsi_9P[7][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_9P[8][0]) * (1. - Gauss2D::Qsi_9P[8][1]), 0.25 * (1. + Gauss2D::Qsi_9P[8][0]) * (1. - Gauss2D::Qsi_9P[8][1]),
	  0.25 * (1. - Gauss2D::Qsi_9P[8][0]) * (1. + Gauss2D::Qsi_9P[8][1]), 0.25 * (1. + Gauss2D::Qsi_9P[8][0]) * (1. + Gauss2D::Qsi_9P[8][1]) } };

template<> const double O2P2::Prep::Elem::Elem_Rect4_IP<16>::mv_Psi[16][mv_numNodes] = {
	{ 0.25 * (1. - Gauss2D::Qsi_16P[0][0]) * (1. - Gauss2D::Qsi_16P[0][1]), 0.25 * (1. + Gauss2D::Qsi_16P[0][0]) * (1. - Gauss2D::Qsi_16P[0][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[0][0]) * (1. + Gauss2D::Qsi_16P[0][1]), 0.25 * (1. + Gauss2D::Qsi_16P[0][0]) * (1. + Gauss2D::Qsi_16P[0][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[1][0]) * (1. - Gauss2D::Qsi_16P[1][1]), 0.25 * (1. + Gauss2D::Qsi_16P[1][0]) * (1. - Gauss2D::Qsi_16P[1][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[1][0]) * (1. + Gauss2D::Qsi_16P[1][1]), 0.25 * (1. + Gauss2D::Qsi_16P[1][0]) * (1. + Gauss2D::Qsi_16P[1][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[2][0]) * (1. - Gauss2D::Qsi_16P[2][1]), 0.25 * (1. + Gauss2D::Qsi_16P[2][0]) * (1. - Gauss2D::Qsi_16P[2][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[2][0]) * (1. + Gauss2D::Qsi_16P[2][1]), 0.25 * (1. + Gauss2D::Qsi_16P[2][0]) * (1. + Gauss2D::Qsi_16P[2][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[3][0]) * (1. - Gauss2D::Qsi_16P[3][1]), 0.25 * (1. + Gauss2D::Qsi_16P[3][0]) * (1. - Gauss2D::Qsi_16P[3][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[3][0]) * (1. + Gauss2D::Qsi_16P[3][1]), 0.25 * (1. + Gauss2D::Qsi_16P[3][0]) * (1. + Gauss2D::Qsi_16P[3][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[4][0]) * (1. - Gauss2D::Qsi_16P[4][1]), 0.25 * (1. + Gauss2D::Qsi_16P[4][0]) * (1. - Gauss2D::Qsi_16P[4][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[4][0]) * (1. + Gauss2D::Qsi_16P[4][1]), 0.25 * (1. + Gauss2D::Qsi_16P[4][0]) * (1. + Gauss2D::Qsi_16P[4][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[5][0]) * (1. - Gauss2D::Qsi_16P[5][1]), 0.25 * (1. + Gauss2D::Qsi_16P[5][0]) * (1. - Gauss2D::Qsi_16P[5][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[5][0]) * (1. + Gauss2D::Qsi_16P[5][1]), 0.25 * (1. + Gauss2D::Qsi_16P[5][0]) * (1. + Gauss2D::Qsi_16P[5][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[6][0]) * (1. - Gauss2D::Qsi_16P[6][1]), 0.25 * (1. + Gauss2D::Qsi_16P[6][0]) * (1. - Gauss2D::Qsi_16P[6][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[6][0]) * (1. + Gauss2D::Qsi_16P[6][1]), 0.25 * (1. + Gauss2D::Qsi_16P[6][0]) * (1. + Gauss2D::Qsi_16P[6][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[7][0]) * (1. - Gauss2D::Qsi_16P[7][1]), 0.25 * (1. + Gauss2D::Qsi_16P[7][0]) * (1. - Gauss2D::Qsi_16P[7][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[7][0]) * (1. + Gauss2D::Qsi_16P[7][1]), 0.25 * (1. + Gauss2D::Qsi_16P[7][0]) * (1. + Gauss2D::Qsi_16P[7][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[8][0]) * (1. - Gauss2D::Qsi_16P[8][1]), 0.25 * (1. + Gauss2D::Qsi_16P[8][0]) * (1. - Gauss2D::Qsi_16P[8][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[8][0]) * (1. + Gauss2D::Qsi_16P[8][1]), 0.25 * (1. + Gauss2D::Qsi_16P[8][0]) * (1. + Gauss2D::Qsi_16P[8][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[9][0]) * (1. - Gauss2D::Qsi_16P[9][1]), 0.25 * (1. + Gauss2D::Qsi_16P[9][0]) * (1. - Gauss2D::Qsi_16P[9][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[9][0]) * (1. + Gauss2D::Qsi_16P[9][1]), 0.25 * (1. + Gauss2D::Qsi_16P[9][0]) * (1. + Gauss2D::Qsi_16P[9][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[10][0]) * (1. - Gauss2D::Qsi_16P[10][1]), 0.25 * (1. + Gauss2D::Qsi_16P[10][0]) * (1. - Gauss2D::Qsi_16P[10][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[10][0]) * (1. + Gauss2D::Qsi_16P[10][1]), 0.25 * (1. + Gauss2D::Qsi_16P[10][0]) * (1. + Gauss2D::Qsi_16P[10][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[11][0]) * (1. - Gauss2D::Qsi_16P[11][1]), 0.25 * (1. + Gauss2D::Qsi_16P[11][0]) * (1. - Gauss2D::Qsi_16P[11][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[11][0]) * (1. + Gauss2D::Qsi_16P[11][1]), 0.25 * (1. + Gauss2D::Qsi_16P[11][0]) * (1. + Gauss2D::Qsi_16P[11][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[12][0]) * (1. - Gauss2D::Qsi_16P[12][1]), 0.25 * (1. + Gauss2D::Qsi_16P[12][0]) * (1. - Gauss2D::Qsi_16P[12][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[12][0]) * (1. + Gauss2D::Qsi_16P[12][1]), 0.25 * (1. + Gauss2D::Qsi_16P[12][0]) * (1. + Gauss2D::Qsi_16P[12][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[13][0]) * (1. - Gauss2D::Qsi_16P[13][1]), 0.25 * (1. + Gauss2D::Qsi_16P[13][0]) * (1. - Gauss2D::Qsi_16P[13][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[13][0]) * (1. + Gauss2D::Qsi_16P[13][1]), 0.25 * (1. + Gauss2D::Qsi_16P[13][0]) * (1. + Gauss2D::Qsi_16P[13][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[14][0]) * (1. - Gauss2D::Qsi_16P[14][1]), 0.25 * (1. + Gauss2D::Qsi_16P[14][0]) * (1. - Gauss2D::Qsi_16P[14][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[14][0]) * (1. + Gauss2D::Qsi_16P[14][1]), 0.25 * (1. + Gauss2D::Qsi_16P[14][0]) * (1. + Gauss2D::Qsi_16P[14][1]) },
	{ 0.25 * (1. - Gauss2D::Qsi_16P[15][0]) * (1. - Gauss2D::Qsi_16P[15][1]), 0.25 * (1. + Gauss2D::Qsi_16P[15][0]) * (1. - Gauss2D::Qsi_16P[15][1]),
	  0.25 * (1. - Gauss2D::Qsi_16P[15][0]) * (1. + Gauss2D::Qsi_16P[15][1]), 0.25 * (1. + Gauss2D::Qsi_16P[15][0]) * (1. + Gauss2D::Qsi_16P[15][1]) } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
template<> const double O2P2::Prep::Elem::Elem_Rect4_IP<4>::mv_DPsi[4][mv_numNodes][mv_Dim] = {
	{ { -0.25 * (1. - Gauss2D::Qsi_4P[0][1]), -0.25 * (1. - Gauss2D::Qsi_4P[0][0]) }, {  0.25 * (1. - Gauss2D::Qsi_4P[0][1]), -0.25 * (1. + Gauss2D::Qsi_4P[0][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_4P[0][1]),  0.25 * (1. - Gauss2D::Qsi_4P[0][0]) }, {  0.25 * (1. + Gauss2D::Qsi_4P[0][1]),  0.25 * (1. + Gauss2D::Qsi_4P[0][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_4P[1][1]), -0.25 * (1. - Gauss2D::Qsi_4P[1][0]) }, {  0.25 * (1. - Gauss2D::Qsi_4P[1][1]), -0.25 * (1. + Gauss2D::Qsi_4P[1][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_4P[1][1]),  0.25 * (1. - Gauss2D::Qsi_4P[1][0]) }, {  0.25 * (1. + Gauss2D::Qsi_4P[1][1]),  0.25 * (1. + Gauss2D::Qsi_4P[1][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_4P[2][1]), -0.25 * (1. - Gauss2D::Qsi_4P[2][0]) }, {  0.25 * (1. - Gauss2D::Qsi_4P[2][1]), -0.25 * (1. + Gauss2D::Qsi_4P[2][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_4P[2][1]),  0.25 * (1. - Gauss2D::Qsi_4P[2][0]) }, {  0.25 * (1. + Gauss2D::Qsi_4P[2][1]),  0.25 * (1. + Gauss2D::Qsi_4P[2][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_4P[3][1]), -0.25 * (1. - Gauss2D::Qsi_4P[3][0]) }, {  0.25 * (1. - Gauss2D::Qsi_4P[3][1]), -0.25 * (1. + Gauss2D::Qsi_4P[3][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_4P[3][1]),  0.25 * (1. - Gauss2D::Qsi_4P[3][0]) }, {  0.25 * (1. + Gauss2D::Qsi_4P[3][1]),  0.25 * (1. + Gauss2D::Qsi_4P[3][0]) } } };

template<> const double O2P2::Prep::Elem::Elem_Rect4_IP<9>::mv_DPsi[9][mv_numNodes][mv_Dim] = {
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[0][1]), -0.25 * (1. - Gauss2D::Qsi_9P[0][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[0][1]), -0.25 * (1. + Gauss2D::Qsi_9P[0][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[0][1]),  0.25 * (1. - Gauss2D::Qsi_9P[0][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[0][1]),  0.25 * (1. + Gauss2D::Qsi_9P[0][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[1][1]), -0.25 * (1. - Gauss2D::Qsi_9P[1][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[1][1]), -0.25 * (1. + Gauss2D::Qsi_9P[1][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[1][1]),  0.25 * (1. - Gauss2D::Qsi_9P[1][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[1][1]),  0.25 * (1. + Gauss2D::Qsi_9P[1][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[2][1]), -0.25 * (1. - Gauss2D::Qsi_9P[2][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[2][1]), -0.25 * (1. + Gauss2D::Qsi_9P[2][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[2][1]),  0.25 * (1. - Gauss2D::Qsi_9P[2][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[2][1]),  0.25 * (1. + Gauss2D::Qsi_9P[2][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[3][1]), -0.25 * (1. - Gauss2D::Qsi_9P[3][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[3][1]), -0.25 * (1. + Gauss2D::Qsi_9P[3][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[3][1]),  0.25 * (1. - Gauss2D::Qsi_9P[3][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[3][1]),  0.25 * (1. + Gauss2D::Qsi_9P[3][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[4][1]), -0.25 * (1. - Gauss2D::Qsi_9P[4][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[4][1]), -0.25 * (1. + Gauss2D::Qsi_9P[4][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[4][1]),  0.25 * (1. - Gauss2D::Qsi_9P[4][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[4][1]),  0.25 * (1. + Gauss2D::Qsi_9P[4][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[5][1]), -0.25 * (1. - Gauss2D::Qsi_9P[5][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[5][1]), -0.25 * (1. + Gauss2D::Qsi_9P[5][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[5][1]),  0.25 * (1. - Gauss2D::Qsi_9P[5][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[5][1]),  0.25 * (1. + Gauss2D::Qsi_9P[5][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[6][1]), -0.25 * (1. - Gauss2D::Qsi_9P[6][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[6][1]), -0.25 * (1. + Gauss2D::Qsi_9P[6][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[6][1]),  0.25 * (1. - Gauss2D::Qsi_9P[6][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[6][1]),  0.25 * (1. + Gauss2D::Qsi_9P[6][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[7][1]), -0.25 * (1. - Gauss2D::Qsi_9P[7][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[7][1]), -0.25 * (1. + Gauss2D::Qsi_9P[7][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[7][1]),  0.25 * (1. - Gauss2D::Qsi_9P[7][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[7][1]),  0.25 * (1. + Gauss2D::Qsi_9P[7][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_9P[8][1]), -0.25 * (1. - Gauss2D::Qsi_9P[8][0]) }, {  0.25 * (1. - Gauss2D::Qsi_9P[8][1]), -0.25 * (1. + Gauss2D::Qsi_9P[8][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_9P[8][1]),  0.25 * (1. - Gauss2D::Qsi_9P[8][0]) }, {  0.25 * (1. + Gauss2D::Qsi_9P[8][1]),  0.25 * (1. + Gauss2D::Qsi_9P[8][0]) } } };

template<> const double O2P2::Prep::Elem::Elem_Rect4_IP<16>::mv_DPsi[16][mv_numNodes][mv_Dim] = {
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[0][1]), -0.25 * (1. - Gauss2D::Qsi_16P[0][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[0][1]), -0.25 * (1. + Gauss2D::Qsi_16P[0][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[0][1]),  0.25 * (1. - Gauss2D::Qsi_16P[0][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[0][1]),  0.25 * (1. + Gauss2D::Qsi_16P[0][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[1][1]), -0.25 * (1. - Gauss2D::Qsi_16P[1][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[1][1]), -0.25 * (1. + Gauss2D::Qsi_16P[1][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[1][1]),  0.25 * (1. - Gauss2D::Qsi_16P[1][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[1][1]),  0.25 * (1. + Gauss2D::Qsi_16P[1][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[2][1]), -0.25 * (1. - Gauss2D::Qsi_16P[2][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[2][1]), -0.25 * (1. + Gauss2D::Qsi_16P[2][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[2][1]),  0.25 * (1. - Gauss2D::Qsi_16P[2][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[2][1]),  0.25 * (1. + Gauss2D::Qsi_16P[2][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[3][1]), -0.25 * (1. - Gauss2D::Qsi_16P[3][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[3][1]), -0.25 * (1. + Gauss2D::Qsi_16P[3][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[3][1]),  0.25 * (1. - Gauss2D::Qsi_16P[3][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[3][1]),  0.25 * (1. + Gauss2D::Qsi_16P[3][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[4][1]), -0.25 * (1. - Gauss2D::Qsi_16P[4][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[4][1]), -0.25 * (1. + Gauss2D::Qsi_16P[4][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[4][1]),  0.25 * (1. - Gauss2D::Qsi_16P[4][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[4][1]),  0.25 * (1. + Gauss2D::Qsi_16P[4][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[5][1]), -0.25 * (1. - Gauss2D::Qsi_16P[5][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[5][1]), -0.25 * (1. + Gauss2D::Qsi_16P[5][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[5][1]),  0.25 * (1. - Gauss2D::Qsi_16P[5][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[5][1]),  0.25 * (1. + Gauss2D::Qsi_16P[5][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[6][1]), -0.25 * (1. - Gauss2D::Qsi_16P[6][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[6][1]), -0.25 * (1. + Gauss2D::Qsi_16P[6][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[6][1]),  0.25 * (1. - Gauss2D::Qsi_16P[6][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[6][1]),  0.25 * (1. + Gauss2D::Qsi_16P[6][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[7][1]), -0.25 * (1. - Gauss2D::Qsi_16P[7][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[7][1]), -0.25 * (1. + Gauss2D::Qsi_16P[7][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[7][1]),  0.25 * (1. - Gauss2D::Qsi_16P[7][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[7][1]),  0.25 * (1. + Gauss2D::Qsi_16P[7][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[8][1]), -0.25 * (1. - Gauss2D::Qsi_16P[8][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[8][1]), -0.25 * (1. + Gauss2D::Qsi_16P[8][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[8][1]),  0.25 * (1. - Gauss2D::Qsi_16P[8][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[8][1]),  0.25 * (1. + Gauss2D::Qsi_16P[8][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[9][1]), -0.25 * (1. - Gauss2D::Qsi_16P[9][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[9][1]), -0.25 * (1. + Gauss2D::Qsi_16P[9][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[9][1]),  0.25 * (1. - Gauss2D::Qsi_16P[9][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[9][1]),  0.25 * (1. + Gauss2D::Qsi_16P[9][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[10][1]), -0.25 * (1. - Gauss2D::Qsi_16P[10][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[10][1]), -0.25 * (1. + Gauss2D::Qsi_16P[10][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[10][1]),  0.25 * (1. - Gauss2D::Qsi_16P[10][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[10][1]),  0.25 * (1. + Gauss2D::Qsi_16P[10][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[11][1]), -0.25 * (1. - Gauss2D::Qsi_16P[11][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[11][1]), -0.25 * (1. + Gauss2D::Qsi_16P[11][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[11][1]),  0.25 * (1. - Gauss2D::Qsi_16P[11][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[11][1]),  0.25 * (1. + Gauss2D::Qsi_16P[11][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[12][1]), -0.25 * (1. - Gauss2D::Qsi_16P[12][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[12][1]), -0.25 * (1. + Gauss2D::Qsi_16P[12][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[12][1]),  0.25 * (1. - Gauss2D::Qsi_16P[12][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[12][1]),  0.25 * (1. + Gauss2D::Qsi_16P[12][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[13][1]), -0.25 * (1. - Gauss2D::Qsi_16P[13][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[13][1]), -0.25 * (1. + Gauss2D::Qsi_16P[13][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[13][1]),  0.25 * (1. - Gauss2D::Qsi_16P[13][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[13][1]),  0.25 * (1. + Gauss2D::Qsi_16P[13][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[14][1]), -0.25 * (1. - Gauss2D::Qsi_16P[14][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[14][1]), -0.25 * (1. + Gauss2D::Qsi_16P[14][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[14][1]),  0.25 * (1. - Gauss2D::Qsi_16P[14][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[14][1]),  0.25 * (1. + Gauss2D::Qsi_16P[14][0]) } },
	{ { -0.25 * (1. - Gauss2D::Qsi_16P[15][1]), -0.25 * (1. - Gauss2D::Qsi_16P[15][0]) }, {  0.25 * (1. - Gauss2D::Qsi_16P[15][1]), -0.25 * (1. + Gauss2D::Qsi_16P[15][0]) },
	  { -0.25 * (1. + Gauss2D::Qsi_16P[15][1]),  0.25 * (1. - Gauss2D::Qsi_16P[15][0]) }, {  0.25 * (1. + Gauss2D::Qsi_16P[15][1]),  0.25 * (1. + Gauss2D::Qsi_16P[15][0]) } } };
