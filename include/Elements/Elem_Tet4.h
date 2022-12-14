// ================================================================================================
// 
// This file is part of O2P2, an object oriented environment for the positional FEM
//
// Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved.
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
//
// Tetrahedral solid element, with linear interpolation functions.
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
			  * @class Elem_Tet4
			  *
			  * @brief Tetrahedral linear element with 4 nodes.
			  * @details Solid element, with linear interpolation functions, tetrahedral shape.
			  * Options for integration points: 1, 4, 10, 11, 14, 15 and 24.
			  * @image html Elem_Tet4.png height=300
			  */
			class Elem_Tet4 : public ElementSolid
			{
			private:
				Elem_Tet4() = delete;

			protected:
				/** Constructor for tetrahedral linear elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Tet4(std::shared_ptr<O2P2::Prep::Material>& Material)
					: ElementSolid(Material) { }

			public:
				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->m_Mat->m_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 1 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << (1 + add) << " " << (2 + add) << " " << (4 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << (1 + add) << " " << (3 + add) << " " << (4 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << (2 + add) << " " << (3 + add) << " " << (4 + add) << " " << this->m_Mat->m_index << "\n";
					return msg.str();
				}

				// Evaluates shape function in the point.
				Eigen::VectorXd getShapeFcOnPoint(const double* Point) override;

				// Evaluates the derivative of shape function in the point.
				Eigen::MatrixXd getShapeDerivOnPoint(const double* Point) override;

				// Returns the number of nodes of current element.
				int getNumNodes() override { return m_NumNodes; }

				// Returns the number of faces of current element.
				int getNumFaces() override { return m_NumFaces; }

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates.
				  */
				bool evaluateXsi(const std::array<double, m_Dim> xsi) override {

					std::array<double, m_Dim + 1> new_xsi = {};

					for (int i = 0; i < m_Dim; ++i) {
						new_xsi.at(i) = xsi.at(i);
						new_xsi.at(m_Dim) -= xsi.at(i);
					}
					new_xsi.at(m_Dim) += 1.;

					const auto [min, max] = std::minmax_element(new_xsi.begin(), new_xsi.end());
					if (*max < 1.000001 && *min > -0.000001) return true;
					return false;
				}

			protected:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			protected:
				/** @brief Number of Nodes */
				static const int m_NumNodes{ 4 };

				/** @brief Number of Faces */
				static const int m_NumFaces{ 4 };
			};


			/**
			  * @class Elem_Tet4_IP
			  *
			  * @brief Tetrahedral linear element with 4 nodes.
			  * @details Solid element, with linear interpolation functions, tetrahedral shape.
			  *
			  * @tparam nIP Number of integration points. Must be: 1, 4, 10, 11, 14, 15 or 24.
			  */
			template<int nIP>
			class Elem_Tet4_IP : public Elem_Tet4
			{
			private:
				Elem_Tet4_IP() = delete;

			public:
				/** Constructor for tetrahedral linear elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Tet4_IP(std::shared_ptr<O2P2::Prep::Material>& Material)
					: Elem_Tet4(Material) { }

				// Return a vector with values on the integration points currently known in the element' nodes.
				Eigen::VectorXd getValueOnIPs(const double* value) override;

				// Returns a pointer to the first element of the shape functions (with size [nIP][m_NumNodes]).
				double const* getShapeFc() const override { return &m_Psi[0][0]; }

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][m_NumNodes][m_Dim]).
				double const* getShapeDerivative() const override { return &m_DPsi[0][0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return m_weight; }

				// Returns the number of integration points of current element.
				int getNumIP() override { return nIP; }

			private:
				// Weights for numerical integration
				static const double* m_weight;

				// Shape functions
				static const double m_Psi[nIP][m_NumNodes];

				// Shape functions derivative
				static const double m_DPsi[nIP][m_NumNodes][m_Dim];
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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Tet4::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(4);

	Psi(0) = 1. - Point[0] - Point[1] - Point[2];
	Psi(1) = Point[0];
	Psi(2) = Point[1];
	Psi(3) = Point[2];

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd O2P2::Prep::Elem::Elem_Tet4::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(4, 3);

	DPsi(0, 0) = -1.;
	DPsi(1, 0) =  1.;
	DPsi(2, 0) =  0.;
	DPsi(3, 0) =  0.;

	DPsi(0, 1) = -1.;
	DPsi(1, 1) =  0.;
	DPsi(2, 1) =  1.;
	DPsi(3, 1) =  0.;

	DPsi(0, 2) = -1.;
	DPsi(1, 2) =  0.;
	DPsi(2, 2) =  0.;
	DPsi(3, 2) =  1.;

	return DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Prep::Elem::Elem_Tet4::setGeomProperties() {

	const int nVertices = 4;

	// Allocate an array with size m_Dim to which m_Centroid points to.
	m_Centroid = std::make_unique<double[]>(m_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Prep::Node<m_Dim>*, nVertices> vertices;
	vertices[0] = v_Conect[0].get();
	vertices[1] = v_Conect[1].get();
	vertices[2] = v_Conect[2].get();
	vertices[3] = v_Conect[3].get();

	// Memory requested by make_unique is not empty
	for (int i = 0; i < m_Dim; i++) m_Centroid[i] = 0.;

	for (auto& node : vertices) {
		std::array<double, m_Dim> x = node->getInitPos();

		for (int i = 0; i < m_Dim; i++) m_Centroid[i] += x[i];
	}

	// Finishing up
	for (int i = 0; i < m_Dim; i++) m_Centroid[i] /= nVertices;

	// Distance from centroid to vertices
	double dist[nVertices] = {};
	int i = 0;

	for (auto& node : vertices) {
		std::array<double, m_Dim> x = node->getInitPos();

		for (int j = 0; j < m_Dim; j++) {
			dist[i] += (m_Centroid[j] - x[j]) * (m_Centroid[j] - x[j]);
		}
		dist[i] = std::sqrt(dist[i]);

		i++;
	}

	// Since centroid is not the circumcenter, the radius is related to the minimum bounding circle
	m_Radius = *std::max_element(dist, dist + nVertices);
};


// ================================================================================================
//
// Implementation of Member Function: getValueOnIPs
// Return the values on the integration points currently known in the element' nodes
// 
// ================================================================================================
template<int nIP>
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Tet4_IP<nIP>::getValueOnIPs(const double* value) {

	// return value
	Eigen::VectorXd valueOnIp = Eigen::VectorXd::Zero(this->m_NumNodes);

	for (int i = 0; i < nIP; i++) {
		for (int j = 0; j < this->m_NumNodes; j++) {
			valueOnIp(i) += value[i] * m_Psi[i][j];
		}
	}

	return valueOnIp;
};


// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<1>::m_weight = &Hammer3D::Wg_1P[0];
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<4>::m_weight = &Hammer3D::Wg_4P[0];
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<10>::m_weight = &Hammer3D::Wg_10P[0];
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<11>::m_weight = &Hammer3D::Wg_11P[0];
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<14>::m_weight = &Hammer3D::Wg_14P[0];
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<15>::m_weight = &Hammer3D::Wg_15P[0];
template<> const double* O2P2::Prep::Elem::Elem_Tet4_IP<24>::m_weight = &Hammer3D::Wg_24P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<1>::m_Psi[1][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_1P[0][0] - Hammer3D::Qsi_1P[0][1] - Hammer3D::Qsi_1P[0][2], Hammer3D::Qsi_1P[0][0], Hammer3D::Qsi_1P[0][1], Hammer3D::Qsi_1P[0][2] } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<4>::m_Psi[4][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_4P[0][0] - Hammer3D::Qsi_4P[0][1] - Hammer3D::Qsi_4P[0][2], Hammer3D::Qsi_4P[0][0], Hammer3D::Qsi_4P[0][1], Hammer3D::Qsi_4P[0][2] },
	{ 1. - Hammer3D::Qsi_4P[1][0] - Hammer3D::Qsi_4P[1][1] - Hammer3D::Qsi_4P[1][2], Hammer3D::Qsi_4P[1][0], Hammer3D::Qsi_4P[1][1], Hammer3D::Qsi_4P[1][2] },
	{ 1. - Hammer3D::Qsi_4P[2][0] - Hammer3D::Qsi_4P[2][1] - Hammer3D::Qsi_4P[2][2], Hammer3D::Qsi_4P[2][0], Hammer3D::Qsi_4P[2][1], Hammer3D::Qsi_4P[2][2] },
	{ 1. - Hammer3D::Qsi_4P[3][0] - Hammer3D::Qsi_4P[3][1] - Hammer3D::Qsi_4P[3][2], Hammer3D::Qsi_4P[3][0], Hammer3D::Qsi_4P[3][1], Hammer3D::Qsi_4P[3][2] } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<10>::m_Psi[10][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_10P[0][0] - Hammer3D::Qsi_10P[0][1] - Hammer3D::Qsi_10P[0][2], Hammer3D::Qsi_10P[0][0], Hammer3D::Qsi_10P[0][1], Hammer3D::Qsi_10P[0][2] },
	{ 1. - Hammer3D::Qsi_10P[1][0] - Hammer3D::Qsi_10P[1][1] - Hammer3D::Qsi_10P[1][2], Hammer3D::Qsi_10P[1][0], Hammer3D::Qsi_10P[1][1], Hammer3D::Qsi_10P[1][2] },
	{ 1. - Hammer3D::Qsi_10P[2][0] - Hammer3D::Qsi_10P[2][1] - Hammer3D::Qsi_10P[2][2], Hammer3D::Qsi_10P[2][0], Hammer3D::Qsi_10P[2][1], Hammer3D::Qsi_10P[2][2] },
	{ 1. - Hammer3D::Qsi_10P[3][0] - Hammer3D::Qsi_10P[3][1] - Hammer3D::Qsi_10P[3][2], Hammer3D::Qsi_10P[3][0], Hammer3D::Qsi_10P[3][1], Hammer3D::Qsi_10P[3][2] },
	{ 1. - Hammer3D::Qsi_10P[4][0] - Hammer3D::Qsi_10P[4][1] - Hammer3D::Qsi_10P[4][2], Hammer3D::Qsi_10P[4][0], Hammer3D::Qsi_10P[4][1], Hammer3D::Qsi_10P[4][2] },
	{ 1. - Hammer3D::Qsi_10P[5][0] - Hammer3D::Qsi_10P[5][1] - Hammer3D::Qsi_10P[5][2], Hammer3D::Qsi_10P[5][0], Hammer3D::Qsi_10P[5][1], Hammer3D::Qsi_10P[5][2] },
	{ 1. - Hammer3D::Qsi_10P[6][0] - Hammer3D::Qsi_10P[6][1] - Hammer3D::Qsi_10P[6][2], Hammer3D::Qsi_10P[6][0], Hammer3D::Qsi_10P[6][1], Hammer3D::Qsi_10P[6][2] },
	{ 1. - Hammer3D::Qsi_10P[7][0] - Hammer3D::Qsi_10P[7][1] - Hammer3D::Qsi_10P[7][2], Hammer3D::Qsi_10P[7][0], Hammer3D::Qsi_10P[7][1], Hammer3D::Qsi_10P[7][2] },
	{ 1. - Hammer3D::Qsi_10P[8][0] - Hammer3D::Qsi_10P[8][1] - Hammer3D::Qsi_10P[8][2], Hammer3D::Qsi_10P[8][0], Hammer3D::Qsi_10P[8][1], Hammer3D::Qsi_10P[8][2] },
	{ 1. - Hammer3D::Qsi_10P[9][0] - Hammer3D::Qsi_10P[9][1] - Hammer3D::Qsi_10P[9][2], Hammer3D::Qsi_10P[9][0], Hammer3D::Qsi_10P[9][1], Hammer3D::Qsi_10P[9][2] } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<11>::m_Psi[11][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_11P[0][0] - Hammer3D::Qsi_11P[0][1] - Hammer3D::Qsi_11P[0][2], Hammer3D::Qsi_11P[0][0], Hammer3D::Qsi_11P[0][1], Hammer3D::Qsi_11P[0][2] },
	{ 1. - Hammer3D::Qsi_11P[1][0] - Hammer3D::Qsi_11P[1][1] - Hammer3D::Qsi_11P[1][2], Hammer3D::Qsi_11P[1][0], Hammer3D::Qsi_11P[1][1], Hammer3D::Qsi_11P[1][2] },
	{ 1. - Hammer3D::Qsi_11P[2][0] - Hammer3D::Qsi_11P[2][1] - Hammer3D::Qsi_11P[2][2], Hammer3D::Qsi_11P[2][0], Hammer3D::Qsi_11P[2][1], Hammer3D::Qsi_11P[2][2] },
	{ 1. - Hammer3D::Qsi_11P[3][0] - Hammer3D::Qsi_11P[3][1] - Hammer3D::Qsi_11P[3][2], Hammer3D::Qsi_11P[3][0], Hammer3D::Qsi_11P[3][1], Hammer3D::Qsi_11P[3][2] },
	{ 1. - Hammer3D::Qsi_11P[4][0] - Hammer3D::Qsi_11P[4][1] - Hammer3D::Qsi_11P[4][2], Hammer3D::Qsi_11P[4][0], Hammer3D::Qsi_11P[4][1], Hammer3D::Qsi_11P[4][2] },
	{ 1. - Hammer3D::Qsi_11P[5][0] - Hammer3D::Qsi_11P[5][1] - Hammer3D::Qsi_11P[5][2], Hammer3D::Qsi_11P[5][0], Hammer3D::Qsi_11P[5][1], Hammer3D::Qsi_11P[5][2] },
	{ 1. - Hammer3D::Qsi_11P[6][0] - Hammer3D::Qsi_11P[6][1] - Hammer3D::Qsi_11P[6][2], Hammer3D::Qsi_11P[6][0], Hammer3D::Qsi_11P[6][1], Hammer3D::Qsi_11P[6][2] },
	{ 1. - Hammer3D::Qsi_11P[7][0] - Hammer3D::Qsi_11P[7][1] - Hammer3D::Qsi_11P[7][2], Hammer3D::Qsi_11P[7][0], Hammer3D::Qsi_11P[7][1], Hammer3D::Qsi_11P[7][2] },
	{ 1. - Hammer3D::Qsi_11P[8][0] - Hammer3D::Qsi_11P[8][1] - Hammer3D::Qsi_11P[8][2], Hammer3D::Qsi_11P[8][0], Hammer3D::Qsi_11P[8][1], Hammer3D::Qsi_11P[8][2] },
	{ 1. - Hammer3D::Qsi_11P[9][0] - Hammer3D::Qsi_11P[9][1] - Hammer3D::Qsi_11P[9][2], Hammer3D::Qsi_11P[9][0], Hammer3D::Qsi_11P[9][1], Hammer3D::Qsi_11P[9][2] },
	{ 1. - Hammer3D::Qsi_11P[10][0] - Hammer3D::Qsi_11P[10][1] - Hammer3D::Qsi_11P[10][2], Hammer3D::Qsi_11P[10][0], Hammer3D::Qsi_11P[10][1], Hammer3D::Qsi_11P[10][2] } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<14>::m_Psi[14][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_14P[0][0] - Hammer3D::Qsi_14P[0][1] - Hammer3D::Qsi_14P[0][2], Hammer3D::Qsi_14P[0][0], Hammer3D::Qsi_14P[0][1], Hammer3D::Qsi_14P[0][2] },
	{ 1. - Hammer3D::Qsi_14P[1][0] - Hammer3D::Qsi_14P[1][1] - Hammer3D::Qsi_14P[1][2], Hammer3D::Qsi_14P[1][0], Hammer3D::Qsi_14P[1][1], Hammer3D::Qsi_14P[1][2] },
	{ 1. - Hammer3D::Qsi_14P[2][0] - Hammer3D::Qsi_14P[2][1] - Hammer3D::Qsi_14P[2][2], Hammer3D::Qsi_14P[2][0], Hammer3D::Qsi_14P[2][1], Hammer3D::Qsi_14P[2][2] },
	{ 1. - Hammer3D::Qsi_14P[3][0] - Hammer3D::Qsi_14P[3][1] - Hammer3D::Qsi_14P[3][2], Hammer3D::Qsi_14P[3][0], Hammer3D::Qsi_14P[3][1], Hammer3D::Qsi_14P[3][2] },
	{ 1. - Hammer3D::Qsi_14P[4][0] - Hammer3D::Qsi_14P[4][1] - Hammer3D::Qsi_14P[4][2], Hammer3D::Qsi_14P[4][0], Hammer3D::Qsi_14P[4][1], Hammer3D::Qsi_14P[4][2] },
	{ 1. - Hammer3D::Qsi_14P[5][0] - Hammer3D::Qsi_14P[5][1] - Hammer3D::Qsi_14P[5][2], Hammer3D::Qsi_14P[5][0], Hammer3D::Qsi_14P[5][1], Hammer3D::Qsi_14P[5][2] },
	{ 1. - Hammer3D::Qsi_14P[6][0] - Hammer3D::Qsi_14P[6][1] - Hammer3D::Qsi_14P[6][2], Hammer3D::Qsi_14P[6][0], Hammer3D::Qsi_14P[6][1], Hammer3D::Qsi_14P[6][2] },
	{ 1. - Hammer3D::Qsi_14P[7][0] - Hammer3D::Qsi_14P[7][1] - Hammer3D::Qsi_14P[7][2], Hammer3D::Qsi_14P[7][0], Hammer3D::Qsi_14P[7][1], Hammer3D::Qsi_14P[7][2] },
	{ 1. - Hammer3D::Qsi_14P[8][0] - Hammer3D::Qsi_14P[8][1] - Hammer3D::Qsi_14P[8][2], Hammer3D::Qsi_14P[8][0], Hammer3D::Qsi_14P[8][1], Hammer3D::Qsi_14P[8][2] },
	{ 1. - Hammer3D::Qsi_14P[9][0] - Hammer3D::Qsi_14P[9][1] - Hammer3D::Qsi_14P[9][2], Hammer3D::Qsi_14P[9][0], Hammer3D::Qsi_14P[9][1], Hammer3D::Qsi_14P[9][2] },
	{ 1. - Hammer3D::Qsi_14P[10][0] - Hammer3D::Qsi_14P[10][1] - Hammer3D::Qsi_14P[10][2], Hammer3D::Qsi_14P[10][0], Hammer3D::Qsi_14P[10][1], Hammer3D::Qsi_14P[10][2] },
	{ 1. - Hammer3D::Qsi_14P[11][0] - Hammer3D::Qsi_14P[11][1] - Hammer3D::Qsi_14P[11][2], Hammer3D::Qsi_14P[11][0], Hammer3D::Qsi_14P[11][1], Hammer3D::Qsi_14P[11][2] },
	{ 1. - Hammer3D::Qsi_14P[12][0] - Hammer3D::Qsi_14P[12][1] - Hammer3D::Qsi_14P[12][2], Hammer3D::Qsi_14P[12][0], Hammer3D::Qsi_14P[12][1], Hammer3D::Qsi_14P[12][2] },
	{ 1. - Hammer3D::Qsi_14P[13][0] - Hammer3D::Qsi_14P[13][1] - Hammer3D::Qsi_14P[13][2], Hammer3D::Qsi_14P[13][0], Hammer3D::Qsi_14P[13][1], Hammer3D::Qsi_14P[13][2] } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<15>::m_Psi[15][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_15P[0][0] - Hammer3D::Qsi_15P[0][1] - Hammer3D::Qsi_15P[0][2], Hammer3D::Qsi_15P[0][0], Hammer3D::Qsi_15P[0][1], Hammer3D::Qsi_15P[0][2] },
	{ 1. - Hammer3D::Qsi_15P[1][0] - Hammer3D::Qsi_15P[1][1] - Hammer3D::Qsi_15P[1][2], Hammer3D::Qsi_15P[1][0], Hammer3D::Qsi_15P[1][1], Hammer3D::Qsi_15P[1][2] },
	{ 1. - Hammer3D::Qsi_15P[2][0] - Hammer3D::Qsi_15P[2][1] - Hammer3D::Qsi_15P[2][2], Hammer3D::Qsi_15P[2][0], Hammer3D::Qsi_15P[2][1], Hammer3D::Qsi_15P[2][2] },
	{ 1. - Hammer3D::Qsi_15P[3][0] - Hammer3D::Qsi_15P[3][1] - Hammer3D::Qsi_15P[3][2], Hammer3D::Qsi_15P[3][0], Hammer3D::Qsi_15P[3][1], Hammer3D::Qsi_15P[3][2] },
	{ 1. - Hammer3D::Qsi_15P[4][0] - Hammer3D::Qsi_15P[4][1] - Hammer3D::Qsi_15P[4][2], Hammer3D::Qsi_15P[4][0], Hammer3D::Qsi_15P[4][1], Hammer3D::Qsi_15P[4][2] },
	{ 1. - Hammer3D::Qsi_15P[5][0] - Hammer3D::Qsi_15P[5][1] - Hammer3D::Qsi_15P[5][2], Hammer3D::Qsi_15P[5][0], Hammer3D::Qsi_15P[5][1], Hammer3D::Qsi_15P[5][2] },
	{ 1. - Hammer3D::Qsi_15P[6][0] - Hammer3D::Qsi_15P[6][1] - Hammer3D::Qsi_15P[6][2], Hammer3D::Qsi_15P[6][0], Hammer3D::Qsi_15P[6][1], Hammer3D::Qsi_15P[6][2] },
	{ 1. - Hammer3D::Qsi_15P[7][0] - Hammer3D::Qsi_15P[7][1] - Hammer3D::Qsi_15P[7][2], Hammer3D::Qsi_15P[7][0], Hammer3D::Qsi_15P[7][1], Hammer3D::Qsi_15P[7][2] },
	{ 1. - Hammer3D::Qsi_15P[8][0] - Hammer3D::Qsi_15P[8][1] - Hammer3D::Qsi_15P[8][2], Hammer3D::Qsi_15P[8][0], Hammer3D::Qsi_15P[8][1], Hammer3D::Qsi_15P[8][2] },
	{ 1. - Hammer3D::Qsi_15P[9][0] - Hammer3D::Qsi_15P[9][1] - Hammer3D::Qsi_15P[9][2], Hammer3D::Qsi_15P[9][0], Hammer3D::Qsi_15P[9][1], Hammer3D::Qsi_15P[9][2] },
	{ 1. - Hammer3D::Qsi_15P[10][0] - Hammer3D::Qsi_15P[10][1] - Hammer3D::Qsi_15P[10][2], Hammer3D::Qsi_15P[10][0], Hammer3D::Qsi_15P[10][1], Hammer3D::Qsi_15P[10][2] },
	{ 1. - Hammer3D::Qsi_15P[11][0] - Hammer3D::Qsi_15P[11][1] - Hammer3D::Qsi_15P[11][2], Hammer3D::Qsi_15P[11][0], Hammer3D::Qsi_15P[11][1], Hammer3D::Qsi_15P[11][2] },
	{ 1. - Hammer3D::Qsi_15P[12][0] - Hammer3D::Qsi_15P[12][1] - Hammer3D::Qsi_15P[12][2], Hammer3D::Qsi_15P[12][0], Hammer3D::Qsi_15P[12][1], Hammer3D::Qsi_15P[12][2] },
	{ 1. - Hammer3D::Qsi_15P[13][0] - Hammer3D::Qsi_15P[13][1] - Hammer3D::Qsi_15P[13][2], Hammer3D::Qsi_15P[13][0], Hammer3D::Qsi_15P[13][1], Hammer3D::Qsi_15P[13][2] },
	{ 1. - Hammer3D::Qsi_15P[14][0] - Hammer3D::Qsi_15P[14][1] - Hammer3D::Qsi_15P[14][2], Hammer3D::Qsi_15P[14][0], Hammer3D::Qsi_15P[14][1], Hammer3D::Qsi_15P[14][2] } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<24>::m_Psi[24][m_NumNodes] = {
	{ 1. - Hammer3D::Qsi_24P[0][0] - Hammer3D::Qsi_24P[0][1] - Hammer3D::Qsi_24P[0][2], Hammer3D::Qsi_24P[0][0], Hammer3D::Qsi_24P[0][1], Hammer3D::Qsi_24P[0][2] },
	{ 1. - Hammer3D::Qsi_24P[1][0] - Hammer3D::Qsi_24P[1][1] - Hammer3D::Qsi_24P[1][2], Hammer3D::Qsi_24P[1][0], Hammer3D::Qsi_24P[1][1], Hammer3D::Qsi_24P[1][2] },
	{ 1. - Hammer3D::Qsi_24P[2][0] - Hammer3D::Qsi_24P[2][1] - Hammer3D::Qsi_24P[2][2], Hammer3D::Qsi_24P[2][0], Hammer3D::Qsi_24P[2][1], Hammer3D::Qsi_24P[2][2] },
	{ 1. - Hammer3D::Qsi_24P[3][0] - Hammer3D::Qsi_24P[3][1] - Hammer3D::Qsi_24P[3][2], Hammer3D::Qsi_24P[3][0], Hammer3D::Qsi_24P[3][1], Hammer3D::Qsi_24P[3][2] },
	{ 1. - Hammer3D::Qsi_24P[4][0] - Hammer3D::Qsi_24P[4][1] - Hammer3D::Qsi_24P[4][2], Hammer3D::Qsi_24P[4][0], Hammer3D::Qsi_24P[4][1], Hammer3D::Qsi_24P[4][2] },
	{ 1. - Hammer3D::Qsi_24P[5][0] - Hammer3D::Qsi_24P[5][1] - Hammer3D::Qsi_24P[5][2], Hammer3D::Qsi_24P[5][0], Hammer3D::Qsi_24P[5][1], Hammer3D::Qsi_24P[5][2] },
	{ 1. - Hammer3D::Qsi_24P[6][0] - Hammer3D::Qsi_24P[6][1] - Hammer3D::Qsi_24P[6][2], Hammer3D::Qsi_24P[6][0], Hammer3D::Qsi_24P[6][1], Hammer3D::Qsi_24P[6][2] },
	{ 1. - Hammer3D::Qsi_24P[7][0] - Hammer3D::Qsi_24P[7][1] - Hammer3D::Qsi_24P[7][2], Hammer3D::Qsi_24P[7][0], Hammer3D::Qsi_24P[7][1], Hammer3D::Qsi_24P[7][2] },
	{ 1. - Hammer3D::Qsi_24P[8][0] - Hammer3D::Qsi_24P[8][1] - Hammer3D::Qsi_24P[8][2], Hammer3D::Qsi_24P[8][0], Hammer3D::Qsi_24P[8][1], Hammer3D::Qsi_24P[8][2] },
	{ 1. - Hammer3D::Qsi_24P[9][0] - Hammer3D::Qsi_24P[9][1] - Hammer3D::Qsi_24P[9][2], Hammer3D::Qsi_24P[9][0], Hammer3D::Qsi_24P[9][1], Hammer3D::Qsi_24P[9][2] },
	{ 1. - Hammer3D::Qsi_24P[10][0] - Hammer3D::Qsi_24P[10][1] - Hammer3D::Qsi_24P[10][2], Hammer3D::Qsi_24P[10][0], Hammer3D::Qsi_24P[10][1], Hammer3D::Qsi_24P[10][2] },
	{ 1. - Hammer3D::Qsi_24P[11][0] - Hammer3D::Qsi_24P[11][1] - Hammer3D::Qsi_24P[11][2], Hammer3D::Qsi_24P[11][0], Hammer3D::Qsi_24P[11][1], Hammer3D::Qsi_24P[11][2] },
	{ 1. - Hammer3D::Qsi_24P[12][0] - Hammer3D::Qsi_24P[12][1] - Hammer3D::Qsi_24P[12][2], Hammer3D::Qsi_24P[12][0], Hammer3D::Qsi_24P[12][1], Hammer3D::Qsi_24P[12][2] },
	{ 1. - Hammer3D::Qsi_24P[13][0] - Hammer3D::Qsi_24P[13][1] - Hammer3D::Qsi_24P[13][2], Hammer3D::Qsi_24P[13][0], Hammer3D::Qsi_24P[13][1], Hammer3D::Qsi_24P[13][2] },
	{ 1. - Hammer3D::Qsi_24P[14][0] - Hammer3D::Qsi_24P[14][1] - Hammer3D::Qsi_24P[14][2], Hammer3D::Qsi_24P[14][0], Hammer3D::Qsi_24P[14][1], Hammer3D::Qsi_24P[14][2] },
	{ 1. - Hammer3D::Qsi_24P[15][0] - Hammer3D::Qsi_24P[15][1] - Hammer3D::Qsi_24P[15][2], Hammer3D::Qsi_24P[15][0], Hammer3D::Qsi_24P[15][1], Hammer3D::Qsi_24P[15][2] },
	{ 1. - Hammer3D::Qsi_24P[16][0] - Hammer3D::Qsi_24P[16][1] - Hammer3D::Qsi_24P[16][2], Hammer3D::Qsi_24P[16][0], Hammer3D::Qsi_24P[16][1], Hammer3D::Qsi_24P[16][2] },
	{ 1. - Hammer3D::Qsi_24P[17][0] - Hammer3D::Qsi_24P[17][1] - Hammer3D::Qsi_24P[17][2], Hammer3D::Qsi_24P[17][0], Hammer3D::Qsi_24P[17][1], Hammer3D::Qsi_24P[17][2] },
	{ 1. - Hammer3D::Qsi_24P[18][0] - Hammer3D::Qsi_24P[18][1] - Hammer3D::Qsi_24P[18][2], Hammer3D::Qsi_24P[18][0], Hammer3D::Qsi_24P[18][1], Hammer3D::Qsi_24P[18][2] },
	{ 1. - Hammer3D::Qsi_24P[19][0] - Hammer3D::Qsi_24P[19][1] - Hammer3D::Qsi_24P[19][2], Hammer3D::Qsi_24P[19][0], Hammer3D::Qsi_24P[19][1], Hammer3D::Qsi_24P[19][2] },
	{ 1. - Hammer3D::Qsi_24P[20][0] - Hammer3D::Qsi_24P[20][1] - Hammer3D::Qsi_24P[20][2], Hammer3D::Qsi_24P[20][0], Hammer3D::Qsi_24P[20][1], Hammer3D::Qsi_24P[20][2] },
	{ 1. - Hammer3D::Qsi_24P[21][0] - Hammer3D::Qsi_24P[21][1] - Hammer3D::Qsi_24P[21][2], Hammer3D::Qsi_24P[21][0], Hammer3D::Qsi_24P[21][1], Hammer3D::Qsi_24P[21][2] },
	{ 1. - Hammer3D::Qsi_24P[22][0] - Hammer3D::Qsi_24P[22][1] - Hammer3D::Qsi_24P[22][2], Hammer3D::Qsi_24P[22][0], Hammer3D::Qsi_24P[22][1], Hammer3D::Qsi_24P[22][2] },
	{ 1. - Hammer3D::Qsi_24P[23][0] - Hammer3D::Qsi_24P[23][1] - Hammer3D::Qsi_24P[23][2], Hammer3D::Qsi_24P[23][0], Hammer3D::Qsi_24P[23][1], Hammer3D::Qsi_24P[23][2] } };


// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<1>::m_DPsi[1][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<4>::m_DPsi[4][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<10>::m_DPsi[10][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<11>::m_DPsi[11][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<14>::m_DPsi[14][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<15>::m_DPsi[15][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };

template<> const double O2P2::Prep::Elem::Elem_Tet4_IP<24>::m_DPsi[24][m_NumNodes][m_Dim] = {
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } },
	{ { -1., -1., -1. } , { 1., 0., 0. }, { 0., 1., 0. } , { 0., 0., 1., } } };
