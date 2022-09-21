// ================================================================================================
//
// This file is part of PosFEM++, an object oriented environment for the positional FEM
//
// Copyright(C) 2021 Rogerio Carrazedo - All Rights Reserved.
// 
// This source code form is subject to the terms of the GNU General Public License.
// If a copy of GPL was not distributed with this file, you can obtain one at
// https://www.gnu.org/licenses/.
// 
// ================================================================================================
//
// Solid element, with linear interpolation functions, prism shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

namespace O2P2 {
	namespace Prep {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Pri6
			  *
			  * @brief Prismatic linear element with 6 nodes.
			  * @details Solid element, with linear interpolation functions, prism shape.
			  * @image html Elem_Pri6.png height=300
			  */
			class Elem_Pri6 : public ElementSolid
			{
			private:
				Elem_Pri6() = delete;

			public:
				/** Constructor for prism linear elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Pri6(std::shared_ptr<O2P2::Prep::Material>& Material)
					: ElementSolid(Material) { }

				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << this->v_Conect[3]->m_index + add << " " << this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->v_Conect[4]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->v_Conect[5]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " " << this->m_Mat->m_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;

					msg << "2 1 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "2 1 " << (4 + add) << " " << (5 + add) << " " << (6 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << (1 + add) << " " << (2 + add) << " " << (4 + add) << " " << (5 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << (1 + add) << " " << (3 + add) << " " << (4 + add) << " " << (6 + add) << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << (2 + add) << " " << (3 + add) << " " << (5 + add) << " " << (6 + add) << " " << this->m_Mat->m_index << "\n";

					return msg.str();
				}

				// Evaluates shape function in the point.
				Eigen::VectorXd getShapeFcOnPoint(const double* Point) override;

				// Evaluates the derivative of shape function in the point.
				Eigen::MatrixXd getShapeDerivOnPoint(const double* Point) override;

				// Return a vector with values on the integration points currently known in the element' nodes.
				Eigen::VectorXd getValueOnIPs(const double* value) override;

				// Returns a pointer to the first element of the shape functions (with size [nIP][m_NumNodes]).
				double const* getShapeFc() const override { return &m_Psi[0][0]; }

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][m_NumNodes][m_Dim]).
				double const* getShapeDerivative() const override { return &m_DPsi[0][0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return &m_weight[0]; }

				// Returns the number of nodes of current element.
				int getNumNodes() override { return m_NumNodes; }

				// Returns the number of faces of current element.
				int getNumFaces() override { return m_NumFaces; }

				// Returns the number of integration points of current element.
				int getNumIP() override { return m_NumIP; }

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates.
				  */
				bool evaluateXsi(const std::array<double, m_Dim> xsi) override {

					std::array<double, m_Dim + 1> new_xsi = {};

					for (int i = 0; i < m_Dim - 1; ++i) {
						new_xsi.at(i) = xsi.at(i);
						new_xsi.at(m_Dim - 1) -= xsi.at(i);
					}
					new_xsi.at(m_Dim - 1) += 1.;
					new_xsi.at(m_Dim) = xsi.at(m_Dim - 1);

					if (*std::max_element(new_xsi.begin(), new_xsi.end()) < 1.000001 && *std::min_element(new_xsi.begin(), new_xsi.end() - 1) > -0.000001 && new_xsi.at(m_Dim) > -1.000001) return true;
					return false;
				}

			private:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			private:
				/** @brief Number of Nodes */
				static const int m_NumNodes{ 6 };

				/** @brief Number of Integration Points */
				static const int m_NumIP{ 8 };

				/** @brief Number of Faces */
				static const int m_NumFaces{ 5 };

				/** @brief Weights for numerical integration */
				static const double m_weight[m_NumIP];

				/** @brief Shape functions */
				static const double m_Psi[m_NumIP][m_NumNodes];

				/** @brief Shape functions derivative */
				static const double m_DPsi[m_NumIP][m_NumNodes][m_Dim];

				/** @brief Integration points */
				static const double m_xsi[m_NumIP][m_Dim];
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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Pri6::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(6);

	Psi(0) = 0.5 * (1. - Point[0] - Point[1]) * (1. - Point[2]);
	Psi(1) = 0.5 * Point[0] * (1. - Point[2]);
	Psi(2) = 0.5 * Point[1] * (1. - Point[2]);
	Psi(3) = 0.5 * (1. - Point[0] - Point[1]) * (1. + Point[2]);
	Psi(4) = 0.5 * Point[0] * (1. + Point[2]);
	Psi(5) = 0.5 * Point[1] * (1. + Point[2]);

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd O2P2::Prep::Elem::Elem_Pri6::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(6, 3);

	DPsi(0, 0) = -0.5 * (1. - Point[2]);
	DPsi(1, 0) = 0.5 * (1. - Point[2]);
	DPsi(2, 0) = 0.;
	DPsi(3, 0) = -0.5 * (1. + Point[2]);
	DPsi(4, 0) = 0.5 * (1. + Point[2]);
	DPsi(5, 0) = 0.;

	DPsi(0, 1) = -0.5 * (1. - Point[2]);
	DPsi(1, 1) = 0.;
	DPsi(2, 1) = 0.5 * (1. - Point[2]);
	DPsi(3, 1) = -0.5 * (1. + Point[2]);
	DPsi(4, 1) = 0.;
	DPsi(5, 1) = 0.5 * (1. + Point[2]);

	DPsi(0, 2) = -0.5 * (1. - Point[0] - Point[1]);
	DPsi(1, 2) = -0.5 * Point[0];
	DPsi(2, 2) = -0.5 * Point[1];
	DPsi(3, 2) = 0.5 * (1. - Point[0] -Point[1]);
	DPsi(4, 2) = 0.5 * Point[0];
	DPsi(5, 2) = 0.5 * Point[1];

	return DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Prep::Elem::Elem_Pri6::setGeomProperties() {

	const int nVertices = 6;

	// Allocate an array with size m_Dim to which m_Centroid points to.
	m_Centroid = std::make_unique<double[]>(m_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Prep::Node<m_Dim>*, nVertices> vertices;
	vertices[0] = v_Conect[0].get();
	vertices[1] = v_Conect[1].get();
	vertices[2] = v_Conect[2].get();
	vertices[3] = v_Conect[3].get();
	vertices[4] = v_Conect[4].get();
	vertices[5] = v_Conect[5].get();

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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Pri6::getValueOnIPs(const double* value) {

	// return value
	Eigen::VectorXd valueOnIp = Eigen::VectorXd::Zero(m_NumNodes);

	for (int i = 0; i < m_NumIP; i++) {
		for (int j = 0; j < this->m_NumNodes; j++) {
			valueOnIp(i) += value[i] * m_Psi[i][j];
		}
	}

	return valueOnIp;
};


// ================================================================================================
//
// Integration Points
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri6::m_xsi[m_NumIP][m_Dim] =
	{ { 1. / 3., 1. / 3. , -0.5773502691896257645091488 },
	  { 3. / 5., 1. / 5. , -0.5773502691896257645091488 },
	  { 1. / 5., 3. / 5. , -0.5773502691896257645091488 },
	  { 1. / 5., 1. / 5. , -0.5773502691896257645091488 },
	  { 1. / 3., 1. / 3. ,  0.5773502691896257645091488 },
	  { 3. / 5., 1. / 5. ,  0.5773502691896257645091488 },
	  { 1. / 5., 3. / 5. ,  0.5773502691896257645091488 },
	  { 1. / 5., 1. / 5. ,  0.5773502691896257645091488 } };

// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri6::m_weight[m_NumIP] = { -9.0 / 32.0, 25.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0, -9.0 / 32.0, 25.0 / 96.0, 25.0 / 96.0, 25.0 / 96.0 };

// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri6::m_Psi[m_NumIP][m_NumNodes] = {
	{ 0.5 * (1. - m_xsi[0][0] - m_xsi[0][1]) * (1. - m_xsi[0][2]), 0.5 * m_xsi[0][0] * (1. - m_xsi[0][2]), 0.5 * m_xsi[0][1] * (1. - m_xsi[0][2]),
	  0.5 * (1. - m_xsi[0][0] - m_xsi[0][1]) * (1. + m_xsi[0][2]), 0.5 * m_xsi[0][0] * (1. + m_xsi[0][2]), 0.5 * m_xsi[0][1] * (1. + m_xsi[0][2]) },
	{ 0.5 * (1. - m_xsi[1][0] - m_xsi[1][1]) * (1. - m_xsi[1][2]), 0.5 * m_xsi[1][0] * (1. - m_xsi[1][2]), 0.5 * m_xsi[1][1] * (1. - m_xsi[1][2]),
	  0.5 * (1. - m_xsi[1][0] - m_xsi[1][1]) * (1. + m_xsi[1][2]), 0.5 * m_xsi[1][0] * (1. + m_xsi[1][2]), 0.5 * m_xsi[1][1] * (1. + m_xsi[1][2]) },
	{ 0.5 * (1. - m_xsi[2][0] - m_xsi[2][1]) * (1. - m_xsi[2][2]), 0.5 * m_xsi[2][0] * (1. - m_xsi[2][2]), 0.5 * m_xsi[2][1] * (1. - m_xsi[2][2]),
	  0.5 * (1. - m_xsi[2][0] - m_xsi[2][1]) * (1. + m_xsi[2][2]), 0.5 * m_xsi[2][0] * (1. + m_xsi[2][2]), 0.5 * m_xsi[2][1] * (1. + m_xsi[2][2]) },
	{ 0.5 * (1. - m_xsi[3][0] - m_xsi[3][1]) * (1. - m_xsi[3][2]), 0.5 * m_xsi[3][0] * (1. - m_xsi[3][2]), 0.5 * m_xsi[3][1] * (1. - m_xsi[3][2]),
	  0.5 * (1. - m_xsi[3][0] - m_xsi[3][1]) * (1. + m_xsi[3][2]), 0.5 * m_xsi[3][0] * (1. + m_xsi[3][2]), 0.5 * m_xsi[3][1] * (1. + m_xsi[3][2]) },
	{ 0.5 * (1. - m_xsi[4][0] - m_xsi[4][1]) * (1. - m_xsi[4][2]), 0.5 * m_xsi[4][0] * (1. - m_xsi[4][2]), 0.5 * m_xsi[4][1] * (1. - m_xsi[4][2]),
	  0.5 * (1. - m_xsi[4][0] - m_xsi[4][1]) * (1. + m_xsi[4][2]), 0.5 * m_xsi[4][0] * (1. + m_xsi[4][2]), 0.5 * m_xsi[4][1] * (1. + m_xsi[4][2]) },
	{ 0.5 * (1. - m_xsi[5][0] - m_xsi[5][1]) * (1. - m_xsi[5][2]), 0.5 * m_xsi[5][0] * (1. - m_xsi[5][2]), 0.5 * m_xsi[5][1] * (1. - m_xsi[5][2]),
	  0.5 * (1. - m_xsi[5][0] - m_xsi[5][1]) * (1. + m_xsi[5][2]), 0.5 * m_xsi[5][0] * (1. + m_xsi[5][2]), 0.5 * m_xsi[5][1] * (1. + m_xsi[5][2]) },
	{ 0.5 * (1. - m_xsi[6][0] - m_xsi[6][1]) * (1. - m_xsi[6][2]), 0.5 * m_xsi[6][0] * (1. - m_xsi[6][2]), 0.5 * m_xsi[6][1] * (1. - m_xsi[6][2]),
	  0.5 * (1. - m_xsi[6][0] - m_xsi[6][1]) * (1. + m_xsi[6][2]), 0.5 * m_xsi[6][0] * (1. + m_xsi[6][2]), 0.5 * m_xsi[6][1] * (1. + m_xsi[6][2]) },
	{ 0.5 * (1. - m_xsi[7][0] - m_xsi[7][1]) * (1. - m_xsi[7][2]), 0.5 * m_xsi[7][0] * (1. - m_xsi[7][2]), 0.5 * m_xsi[7][1] * (1. - m_xsi[7][2]),
	  0.5 * (1. - m_xsi[7][0] - m_xsi[7][1]) * (1. + m_xsi[7][2]), 0.5 * m_xsi[7][0] * (1. + m_xsi[7][2]), 0.5 * m_xsi[7][1] * (1. + m_xsi[7][2]) } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri6::m_DPsi[m_NumIP][m_NumNodes][m_Dim] =
	{ { { -0.5 * (1. - m_xsi[0][2]), -0.5 * (1. - m_xsi[0][2]), -0.5 * (1. - m_xsi[0][0] - m_xsi[0][1]) },
		{  0.5 * (1. - m_xsi[0][2]), 0., -0.5 * m_xsi[0][0] },
		{  0., 0.5 * (1. - m_xsi[0][2]), -0.5 * m_xsi[0][1] },
		{ -0.5 * (1. + m_xsi[0][2]), -0.5 * (1. + m_xsi[0][2]), 0.5 * (1. - m_xsi[0][0] - m_xsi[0][1]) },
		{  0.5 * (1. + m_xsi[0][2]), 0., 0.5 * m_xsi[0][0] },
		{  0., 0.5 * (1. + m_xsi[0][2]), 0.5 * m_xsi[0][1] } },

	  { { -0.5 * (1. - m_xsi[1][2]), -0.5 * (1. - m_xsi[1][2]), -0.5 * (1. - m_xsi[1][0] - m_xsi[1][1]) },
		{  0.5 * (1. - m_xsi[1][2]), 0., -0.5 * m_xsi[1][0] },
		{  0., 0.5 * (1. - m_xsi[1][2]), -0.5 * m_xsi[1][1] },
		{ -0.5 * (1. + m_xsi[1][2]), -0.5 * (1. + m_xsi[1][2]), 0.5 * (1. - m_xsi[1][0] - m_xsi[1][1]) },
		{  0.5 * (1. + m_xsi[1][2]), 0., 0.5 * m_xsi[1][0] },
		{  0., 0.5 * (1. + m_xsi[1][2]), 0.5 * m_xsi[1][1] } },

	  { { -0.5 * (1. - m_xsi[2][2]), -0.5 * (1. - m_xsi[2][2]), -0.5 * (1. - m_xsi[2][0] - m_xsi[2][1]) },
		{  0.5 * (1. - m_xsi[2][2]), 0., -0.5 * m_xsi[2][0] },
		{  0., 0.5 * (1. - m_xsi[2][2]), -0.5 * m_xsi[2][1] },
		{ -0.5 * (1. + m_xsi[2][2]), -0.5 * (1. + m_xsi[2][2]), 0.5 * (1. - m_xsi[2][0] - m_xsi[2][1]) },
		{  0.5 * (1. + m_xsi[2][2]), 0., 0.5 * m_xsi[2][0] },
		{  0., 0.5 * (1. + m_xsi[2][2]), 0.5 * m_xsi[2][1] } },

	  { { -0.5 * (1. - m_xsi[3][2]), -0.5 * (1. - m_xsi[3][2]), -0.5 * (1. - m_xsi[3][0] - m_xsi[3][1]) },
		{  0.5 * (1. - m_xsi[3][2]), 0., -0.5 * m_xsi[3][0] },
		{  0., 0.5 * (1. - m_xsi[3][2]), -0.5 * m_xsi[3][1] },
		{ -0.5 * (1. + m_xsi[3][2]), -0.5 * (1. + m_xsi[3][2]), 0.5 * (1. - m_xsi[3][0] - m_xsi[3][1]) },
		{  0.5 * (1. + m_xsi[3][2]), 0., 0.5 * m_xsi[3][0] },
		{  0., 0.5 * (1. + m_xsi[3][2]), 0.5 * m_xsi[3][1] } },

	  { { -0.5 * (1. - m_xsi[4][2]), -0.5 * (1. - m_xsi[4][2]), -0.5 * (1. - m_xsi[4][0] - m_xsi[4][1]) },
		{  0.5 * (1. - m_xsi[4][2]), 0., -0.5 * m_xsi[4][0] },
		{  0., 0.5 * (1. - m_xsi[4][2]), -0.5 * m_xsi[4][1] },
		{ -0.5 * (1. + m_xsi[4][2]), -0.5 * (1. + m_xsi[4][2]), 0.5 * (1. - m_xsi[4][0] - m_xsi[4][1]) },
		{  0.5 * (1. + m_xsi[4][2]), 0., 0.5 * m_xsi[4][0] },
		{  0., 0.5 * (1. + m_xsi[4][2]), 0.5 * m_xsi[4][1] } },

	  { { -0.5 * (1. - m_xsi[5][2]), -0.5 * (1. - m_xsi[5][2]), -0.5 * (1. - m_xsi[5][0] - m_xsi[5][1]) },
		{  0.5 * (1. - m_xsi[5][2]), 0., -0.5 * m_xsi[5][0] },
		{  0., 0.5 * (1. - m_xsi[5][2]), -0.5 * m_xsi[5][1] },
		{ -0.5 * (1. + m_xsi[5][2]), -0.5 * (1. + m_xsi[5][2]), 0.5 * (1. - m_xsi[5][0] - m_xsi[5][1]) },
		{  0.5 * (1. + m_xsi[5][2]), 0., 0.5 * m_xsi[5][0] },
		{  0., 0.5 * (1. + m_xsi[5][2]), 0.5 * m_xsi[5][1] } },

	  { { -0.5 * (1. - m_xsi[6][2]), -0.5 * (1. - m_xsi[6][2]), -0.5 * (1. - m_xsi[6][0] - m_xsi[6][1]) },
		{  0.5 * (1. - m_xsi[6][2]), 0., -0.5 * m_xsi[6][0] },
		{  0., 0.5 * (1. - m_xsi[6][2]), -0.5 * m_xsi[6][1] },
		{ -0.5 * (1. + m_xsi[6][2]), -0.5 * (1. + m_xsi[6][2]), 0.5 * (1. - m_xsi[6][0] - m_xsi[6][1]) },
		{  0.5 * (1. + m_xsi[6][2]), 0., 0.5 * m_xsi[6][0] },
		{  0., 0.5 * (1. + m_xsi[6][2]), 0.5 * m_xsi[6][1] } },

	  { { -0.5 * (1. - m_xsi[7][2]), -0.5 * (1. - m_xsi[7][2]), -0.5 * (1. - m_xsi[7][0] - m_xsi[7][1]) },
		{  0.5 * (1. - m_xsi[7][2]), 0., -0.5 * m_xsi[7][0] },
		{  0., 0.5 * (1. - m_xsi[7][2]), -0.5 * m_xsi[7][1] },
		{ -0.5 * (1. + m_xsi[7][2]), -0.5 * (1. + m_xsi[7][2]), 0.5 * (1. - m_xsi[7][0] - m_xsi[7][1]) },
		{  0.5 * (1. + m_xsi[7][2]), 0., 0.5 * m_xsi[7][0] },
		{  0., 0.5 * (1. + m_xsi[7][2]), 0.5 * m_xsi[7][1] } } };
