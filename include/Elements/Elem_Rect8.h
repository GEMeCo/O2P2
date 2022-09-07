// ================================================================================================
//
// This file is part of O2P2, an object oriented environment for the positional FEM
//
// Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved.
// 
// This source code form is subject to the terms of 
// Creative Commons Attribution-NonCommerical 4.0 International license
// 
// ================================================================================================
//
// Plane element, with cubic / linear interpolation functions, rectangular shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

namespace O2P2 {
	namespace Prep {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Rect8
			  *
			  * @brief Quadrangular cubic / linear element with 8 nodes.
			  * @details Plane element, with cubic interpolation functions in one direction and linear function in the other, rectangular shaped.
			  * @image html Elem_Quad8.png height=300
			  */
			class Elem_Rect8 : public ElementPlane
			{
			private:
				Elem_Rect8() = delete;

			public:
				/** Constructor for rectangular cubic / linear elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				explicit Elem_Rect8(std::shared_ptr<O2P2::Prep::Material>& Material, std::shared_ptr<O2P2::Prep::Section>& Section)
					: ElementPlane(Material, Section) { }

				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;

					// It is divided into 3 linear elements
					msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " "
						<< this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " "
						<< this->v_Conect[5]->m_index + add << " " << this->v_Conect[6]->m_index + add << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " "
						<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[7]->m_index + add << " "
						<< this->m_Mat->m_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "3 1 " << (1 + add) << " " << (2 + add) << " " << (5 + add) << " " << (6 + add)
						<< " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << (2 + add) << " " << (3 + add) << " " << (6 + add) << " " << (7 + add)
						<< " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << (3 + add) << " " << (4 + add) << " " << (7 + add) << " " << (8 + add)
						<< " " << this->m_Mat->m_index << "\n";
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
					const auto [min, max] = std::minmax_element(xsi.begin(), xsi.end());
					if (*max < 1.000001 && *min > -1.000001) return true;
					return false;
				}

			private:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			private:
				/** @brief Number of Nodes */
				static const int m_NumNodes{ 8 };

				/** @brief Number of Integration Points */
				static const int m_NumIP{ 8 };

				/** @brief Number of Faces (for output purposes) */
				static const int m_NumFaces{ 3 };

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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Rect8::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(8);

	Psi(0) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] - Point[0] * Point[1] + Point[0] + Point[1] - 1.);
	Psi(1) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] + 27. * Point[0] * Point[1] - 27. * Point[0] - 9. * Point[1] + 9.);
	Psi(2) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] - 27. * Point[0] * Point[1] + 27. * Point[0] - 9. * Point[1] + 9.);
	Psi(3) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] + Point[0] * Point[1] - Point[0] + Point[1] - 1.);
	Psi(4) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] + Point[0] * Point[1] + Point[0] - Point[1] - 1.);
	Psi(5) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 27. * Point[0] + 9. * Point[1] + 9.);
	Psi(6) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] + 27. * Point[0] * Point[1] + 27. * Point[0] + 9. * Point[1] + 9.);
	Psi(7) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] - Point[0] * Point[1] - Point[0] - Point[1] - 1.);

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd O2P2::Prep::Elem::Elem_Rect8::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(8, 2);

	DPsi(0, 0) = 0.03125 * (27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] - 18. * Point[0] * Point[1] + 18. * Point[0] - Point[1] + 1.);
	DPsi(1, 0) = 0.03125 * (-81. * Point[0] * Point[0] * Point[1] + 81. * Point[0] * Point[0] + 18. * Point[0] * Point[1] - 18. * Point[0] + 27. * Point[1] - 27.);
	DPsi(2, 0) = 0.03125 * (81. * Point[0] * Point[0] * Point[1] - 81. * Point[0] * Point[0] + 18. * Point[0] * Point[1] - 18. * Point[0] - 27. * Point[1] + 27.);
	DPsi(3, 0) = 0.03125 * (-27. * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] - 18. * Point[0] * Point[1] + 18. * Point[0] + Point[1] - 1.);
	DPsi(4, 0) = 0.03125 * (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 18. * Point[0] + Point[1] + 1.);
	DPsi(5, 0) = 0.03125 * (81. * Point[0] * Point[0] * Point[1] + 81. * Point[0] * Point[0] - 18. * Point[0] * Point[1] - 18. * Point[0] - 27. * Point[1] - 27.);
	DPsi(6, 0) = 0.03125 * (-81. * Point[0] * Point[0] * Point[1] - 81. * Point[0] * Point[0] - 18. * Point[0] * Point[1] - 18. * Point[0] + 27. * Point[1] + 27.);
	DPsi(7, 0) = 0.03125 * (27. * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 18. * Point[0] - Point[1] - 1.);

	DPsi(0, 1) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] - Point[0] + 1.);
	DPsi(1, 1) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] + 27. * Point[0] - 9.);
	DPsi(2, 1) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] - 27. * Point[0] - 9.);
	DPsi(3, 1) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] + Point[0] + 1.);
	DPsi(4, 1) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] + Point[0] - 1.);
	DPsi(5, 1) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] - 27. * Point[0] + 9.);
	DPsi(6, 1) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] + 27. * Point[0] + 9.);
	DPsi(7, 1) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] - Point[0] - 1.);

	return DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Prep::Elem::Elem_Rect8::setGeomProperties() {

	const int nVertices = 4;

	// Allocate an array with size m_Dim to which m_Centroid points to.
	m_Centroid = std::make_unique<double[]>(m_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Prep::Node<m_Dim>*, nVertices> vertices;
	vertices[0] = v_Conect[0].get();
	vertices[1] = v_Conect[3].get();
	vertices[2] = v_Conect[4].get();
	vertices[3] = v_Conect[7].get();

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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Rect8::getValueOnIPs(const double* value) {

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
inline const double O2P2::Prep::Elem::Elem_Rect8::m_xsi[m_NumIP][m_Dim] = {
	{ -0.861136311594053, -0.577350269189626 } ,
	{ -0.861136311594053,  0.577350269189626 } ,
	{  0.861136311594053, -0.577350269189626 } ,
	{  0.861136311594053,  0.577350269189626 } ,
	{ -0.339981043584856, -0.577350269189626 } ,
	{ -0.339981043584856,  0.577350269189626 } ,
	{  0.339981043584856, -0.577350269189626 } ,
	{  0.339981043584856,  0.577350269189626 } };

// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Rect8::m_weight[m_NumIP] = { 0.347854845137454, 0.347854845137454, 0.347854845137454, 0.347854845137454, 0.652145154862546, 0.652145154862546, 0.652145154862546, 0.652145154862546 };

// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Rect8::m_Psi[m_NumIP][m_NumNodes] = {
	{ 0.03125 * (9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] - m_xsi[0][0] * m_xsi[0][1] + m_xsi[0][0] + m_xsi[0][1] - 1.),
	  0.03125 * (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 9. * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] - 9. * m_xsi[0][1] + 9.),
	  0.03125 * (27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 9. * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] - 9. * m_xsi[0][1] + 9.),
	  0.03125 * (-9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0] * m_xsi[0][1] - m_xsi[0][0] + m_xsi[0][1] - 1.),
	  0.03125 * (-9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0] * m_xsi[0][1] + m_xsi[0][0] - m_xsi[0][1] - 1.),
	  0.03125 * (27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 9. * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] + 9. * m_xsi[0][1] + 9.),
	  0.03125 * (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 9. * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] + 9. * m_xsi[0][1] + 9.),
	  0.03125 * (9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] - m_xsi[0][0] * m_xsi[0][1] - m_xsi[0][0] - m_xsi[0][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] - m_xsi[1][0] * m_xsi[1][1] + m_xsi[1][0] + m_xsi[1][1] - 1.),
	  0.03125 * (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 9. * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] - 9. * m_xsi[1][1] + 9.),
	  0.03125 * (27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 9. * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] - 9. * m_xsi[1][1] + 9.),
	  0.03125 * (-9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0] * m_xsi[1][1] - m_xsi[1][0] + m_xsi[1][1] - 1.),
	  0.03125 * (-9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0] * m_xsi[1][1] + m_xsi[1][0] - m_xsi[1][1] - 1.),
	  0.03125 * (27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 9. * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] + 9. * m_xsi[1][1] + 9.),
	  0.03125 * (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 9. * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] + 9. * m_xsi[1][1] + 9.),
	  0.03125 * (9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] - m_xsi[1][0] * m_xsi[1][1] - m_xsi[1][0] - m_xsi[1][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] - m_xsi[2][0] * m_xsi[2][1] + m_xsi[2][0] + m_xsi[2][1] - 1.),
	  0.03125 * (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 9. * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] - 9. * m_xsi[2][1] + 9.),
	  0.03125 * (27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 9. * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] - 9. * m_xsi[2][1] + 9.),
	  0.03125 * (-9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0] * m_xsi[2][1] - m_xsi[2][0] + m_xsi[2][1] - 1.),
	  0.03125 * (-9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0] * m_xsi[2][1] + m_xsi[2][0] - m_xsi[2][1] - 1.),
	  0.03125 * (27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 9. * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] + 9. * m_xsi[2][1] + 9.),
	  0.03125 * (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 9. * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] + 9. * m_xsi[2][1] + 9.),
	  0.03125 * (9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] - m_xsi[2][0] * m_xsi[2][1] - m_xsi[2][0] - m_xsi[2][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] - m_xsi[3][0] * m_xsi[3][1] + m_xsi[3][0] + m_xsi[3][1] - 1.),
	  0.03125 * (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 9. * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] - 9. * m_xsi[3][1] + 9.),
	  0.03125 * (27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 9. * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] - 9. * m_xsi[3][1] + 9.),
	  0.03125 * (-9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0] * m_xsi[3][1] - m_xsi[3][0] + m_xsi[3][1] - 1.),
	  0.03125 * (-9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0] * m_xsi[3][1] + m_xsi[3][0] - m_xsi[3][1] - 1.),
	  0.03125 * (27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 9. * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] + 9. * m_xsi[3][1] + 9.),
	  0.03125 * (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 9. * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] + 9. * m_xsi[3][1] + 9.),
	  0.03125 * (9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] - m_xsi[3][0] * m_xsi[3][1] - m_xsi[3][0] - m_xsi[3][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] - m_xsi[4][0] * m_xsi[4][1] + m_xsi[4][0] + m_xsi[4][1] - 1.),
	  0.03125 * (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 9. * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] - 9. * m_xsi[4][1] + 9.),
	  0.03125 * (27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 9. * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] - 9. * m_xsi[4][1] + 9.),
	  0.03125 * (-9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0] * m_xsi[4][1] - m_xsi[4][0] + m_xsi[4][1] - 1.),
	  0.03125 * (-9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0] * m_xsi[4][1] + m_xsi[4][0] - m_xsi[4][1] - 1.),
	  0.03125 * (27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 9. * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] + 9. * m_xsi[4][1] + 9.),
	  0.03125 * (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 9. * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] + 9. * m_xsi[4][1] + 9.),
	  0.03125 * (9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] - m_xsi[4][0] * m_xsi[4][1] - m_xsi[4][0] - m_xsi[4][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] - m_xsi[5][0] * m_xsi[5][1] + m_xsi[5][0] + m_xsi[5][1] - 1.),
	  0.03125 * (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 9. * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] - 9. * m_xsi[5][1] + 9.),
	  0.03125 * (27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 9. * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] - 9. * m_xsi[5][1] + 9.),
	  0.03125 * (-9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0] * m_xsi[5][1] - m_xsi[5][0] + m_xsi[5][1] - 1.),
	  0.03125 * (-9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0] * m_xsi[5][1] + m_xsi[5][0] - m_xsi[5][1] - 1.),
	  0.03125 * (27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 9. * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] + 9. * m_xsi[5][1] + 9.),
	  0.03125 * (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 9. * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] + 9. * m_xsi[5][1] + 9.),
	  0.03125 * (9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] - m_xsi[5][0] * m_xsi[5][1] - m_xsi[5][0] - m_xsi[5][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] - m_xsi[6][0] * m_xsi[6][1] + m_xsi[6][0] + m_xsi[6][1] - 1.),
	  0.03125 * (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 9. * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] - 9. * m_xsi[6][1] + 9.),
	  0.03125 * (27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 9. * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] - 9. * m_xsi[6][1] + 9.),
	  0.03125 * (-9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0] * m_xsi[6][1] - m_xsi[6][0] + m_xsi[6][1] - 1.),
	  0.03125 * (-9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0] * m_xsi[6][1] + m_xsi[6][0] - m_xsi[6][1] - 1.),
	  0.03125 * (27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 9. * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] + 9. * m_xsi[6][1] + 9.),
	  0.03125 * (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 9. * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] + 9. * m_xsi[6][1] + 9.),
	  0.03125 * (9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] - m_xsi[6][0] * m_xsi[6][1] - m_xsi[6][0] - m_xsi[6][1] - 1.) },

	{ 0.03125 * (9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] - m_xsi[7][0] * m_xsi[7][1] + m_xsi[7][0] + m_xsi[7][1] - 1.),
	  0.03125 * (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 9. * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] - 9. * m_xsi[7][1] + 9.),
	  0.03125 * (27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 9. * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] - 9. * m_xsi[7][1] + 9.),
	  0.03125 * (-9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0] * m_xsi[7][1] - m_xsi[7][0] + m_xsi[7][1] - 1.),
	  0.03125 * (-9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0] * m_xsi[7][1] + m_xsi[7][0] - m_xsi[7][1] - 1.),
	  0.03125 * (27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 9. * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] + 9. * m_xsi[7][1] + 9.),
	  0.03125 * (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 9. * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] + 9. * m_xsi[7][1] + 9.),
	  0.03125 * (9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] - m_xsi[7][0] * m_xsi[7][1] - m_xsi[7][0] - m_xsi[7][1] - 1.) } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Rect8::m_DPsi[m_NumIP][m_NumNodes][m_Dim] = {
	{ { 0.03125 * (27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][0] - 18. * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] - m_xsi[0][1] + 1.),
		0.03125 * (9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] - m_xsi[0][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 81. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] - 18. * m_xsi[0][0] + 27. * m_xsi[0][1] - 27.),
		0.03125 * (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 81. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] - 18. * m_xsi[0][0] - 27. * m_xsi[0][1] + 27.),
		0.03125 * (27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][0] - 18. * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] + m_xsi[0][1] - 1.),
		0.03125 * (-9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] + m_xsi[0][1] + 1.),
		0.03125 * (-9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 81. * m_xsi[0][0] * m_xsi[0][0] - 18. * m_xsi[0][0] * m_xsi[0][1] - 18. * m_xsi[0][0] - 27. * m_xsi[0][1] - 27.),
		0.03125 * (27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 81. * m_xsi[0][0] * m_xsi[0][0] - 18. * m_xsi[0][0] * m_xsi[0][1] - 18. * m_xsi[0][0] + 27. * m_xsi[0][1] + 27.),
		0.03125 * (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] - m_xsi[0][1] - 1.),
		0.03125 * (9. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 9. * m_xsi[0][0] * m_xsi[0][0] - m_xsi[0][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][0] - 18. * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] - m_xsi[1][1] + 1.),
		0.03125 * (9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] - m_xsi[1][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 81. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] - 18. * m_xsi[1][0] + 27. * m_xsi[1][1] - 27.),
		0.03125 * (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 81. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] - 18. * m_xsi[1][0] - 27. * m_xsi[1][1] + 27.),
		0.03125 * (27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][0] - 18. * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] + m_xsi[1][1] - 1.),
		0.03125 * (-9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] + m_xsi[1][1] + 1.),
		0.03125 * (-9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 81. * m_xsi[1][0] * m_xsi[1][0] - 18. * m_xsi[1][0] * m_xsi[1][1] - 18. * m_xsi[1][0] - 27. * m_xsi[1][1] - 27.),
		0.03125 * (27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 81. * m_xsi[1][0] * m_xsi[1][0] - 18. * m_xsi[1][0] * m_xsi[1][1] - 18. * m_xsi[1][0] + 27. * m_xsi[1][1] + 27.),
		0.03125 * (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] - m_xsi[1][1] - 1.),
		0.03125 * (9. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 9. * m_xsi[1][0] * m_xsi[1][0] - m_xsi[1][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][0] - 18. * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] - m_xsi[2][1] + 1.),
		0.03125 * (9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] - m_xsi[2][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 81. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] - 18. * m_xsi[2][0] + 27. * m_xsi[2][1] - 27.),
		0.03125 * (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 81. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] - 18. * m_xsi[2][0] - 27. * m_xsi[2][1] + 27.),
		0.03125 * (27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][0] - 18. * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] + m_xsi[2][1] - 1.),
		0.03125 * (-9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] + m_xsi[2][1] + 1.),
		0.03125 * (-9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 81. * m_xsi[2][0] * m_xsi[2][0] - 18. * m_xsi[2][0] * m_xsi[2][1] - 18. * m_xsi[2][0] - 27. * m_xsi[2][1] - 27.),
		0.03125 * (27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 81. * m_xsi[2][0] * m_xsi[2][0] - 18. * m_xsi[2][0] * m_xsi[2][1] - 18. * m_xsi[2][0] + 27. * m_xsi[2][1] + 27.),
		0.03125 * (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] - m_xsi[2][1] - 1.),
		0.03125 * (9. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 9. * m_xsi[2][0] * m_xsi[2][0] - m_xsi[2][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][0] - 18. * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] - m_xsi[3][1] + 1.),
		0.03125 * (9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] - m_xsi[3][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 81. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] - 18. * m_xsi[3][0] + 27. * m_xsi[3][1] - 27.),
		0.03125 * (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 81. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] - 18. * m_xsi[3][0] - 27. * m_xsi[3][1] + 27.),
		0.03125 * (27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][0] - 18. * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] + m_xsi[3][1] - 1.),
		0.03125 * (-9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] + m_xsi[3][1] + 1.),
		0.03125 * (-9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 81. * m_xsi[3][0] * m_xsi[3][0] - 18. * m_xsi[3][0] * m_xsi[3][1] - 18. * m_xsi[3][0] - 27. * m_xsi[3][1] - 27.),
		0.03125 * (27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 81. * m_xsi[3][0] * m_xsi[3][0] - 18. * m_xsi[3][0] * m_xsi[3][1] - 18. * m_xsi[3][0] + 27. * m_xsi[3][1] + 27.),
		0.03125 * (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] - m_xsi[3][1] - 1.),
		0.03125 * (9. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 9. * m_xsi[3][0] * m_xsi[3][0] - m_xsi[3][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][0] - 18. * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] - m_xsi[4][1] + 1.),
		0.03125 * (9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] - m_xsi[4][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 81. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] - 18. * m_xsi[4][0] + 27. * m_xsi[4][1] - 27.),
		0.03125 * (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 81. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] - 18. * m_xsi[4][0] - 27. * m_xsi[4][1] + 27.),
		0.03125 * (27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][0] - 18. * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] + m_xsi[4][1] - 1.),
		0.03125 * (-9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] + m_xsi[4][1] + 1.),
		0.03125 * (-9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 81. * m_xsi[4][0] * m_xsi[4][0] - 18. * m_xsi[4][0] * m_xsi[4][1] - 18. * m_xsi[4][0] - 27. * m_xsi[4][1] - 27.),
		0.03125 * (27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 81. * m_xsi[4][0] * m_xsi[4][0] - 18. * m_xsi[4][0] * m_xsi[4][1] - 18. * m_xsi[4][0] + 27. * m_xsi[4][1] + 27.),
		0.03125 * (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] - m_xsi[4][1] - 1.),
		0.03125 * (9. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 9. * m_xsi[4][0] * m_xsi[4][0] - m_xsi[4][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][0] - 18. * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] - m_xsi[5][1] + 1.),
		0.03125 * (9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] - m_xsi[5][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 81. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] - 18. * m_xsi[5][0] + 27. * m_xsi[5][1] - 27.),
		0.03125 * (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 81. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] - 18. * m_xsi[5][0] - 27. * m_xsi[5][1] + 27.),
		0.03125 * (27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][0] - 18. * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] + m_xsi[5][1] - 1.),
		0.03125 * (-9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] + m_xsi[5][1] + 1.),
		0.03125 * (-9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 81. * m_xsi[5][0] * m_xsi[5][0] - 18. * m_xsi[5][0] * m_xsi[5][1] - 18. * m_xsi[5][0] - 27. * m_xsi[5][1] - 27.),
		0.03125 * (27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 81. * m_xsi[5][0] * m_xsi[5][0] - 18. * m_xsi[5][0] * m_xsi[5][1] - 18. * m_xsi[5][0] + 27. * m_xsi[5][1] + 27.),
		0.03125 * (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] - m_xsi[5][1] - 1.),
		0.03125 * (9. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 9. * m_xsi[5][0] * m_xsi[5][0] - m_xsi[5][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][0] - 18. * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] - m_xsi[6][1] + 1.),
		0.03125 * (9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] - m_xsi[6][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 81. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] - 18. * m_xsi[6][0] + 27. * m_xsi[6][1] - 27.),
		0.03125 * (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 81. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] - 18. * m_xsi[6][0] - 27. * m_xsi[6][1] + 27.),
		0.03125 * (27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][0] - 18. * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] + m_xsi[6][1] - 1.),
		0.03125 * (-9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] + m_xsi[6][1] + 1.),
		0.03125 * (-9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 81. * m_xsi[6][0] * m_xsi[6][0] - 18. * m_xsi[6][0] * m_xsi[6][1] - 18. * m_xsi[6][0] - 27. * m_xsi[6][1] - 27.),
		0.03125 * (27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 81. * m_xsi[6][0] * m_xsi[6][0] - 18. * m_xsi[6][0] * m_xsi[6][1] - 18. * m_xsi[6][0] + 27. * m_xsi[6][1] + 27.),
		0.03125 * (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] - m_xsi[6][1] - 1.),
		0.03125 * (9. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 9. * m_xsi[6][0] * m_xsi[6][0] - m_xsi[6][0] - 1.) } },

	{ { 0.03125 * (27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][0] - 18. * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] - m_xsi[7][1] + 1.),
		0.03125 * (9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] - m_xsi[7][0] + 1.) },
	  { 0.03125 * (-81. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 81. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] - 18. * m_xsi[7][0] + 27. * m_xsi[7][1] - 27.),
		0.03125 * (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] - 9.) },
	  { 0.03125 * (81. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 81. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] - 18. * m_xsi[7][0] - 27. * m_xsi[7][1] + 27.),
		0.03125 * (27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] - 9.) },
	  { 0.03125 * (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][0] - 18. * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] + m_xsi[7][1] - 1.),
		0.03125 * (-9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0] + 1.) },
	  { 0.03125 * (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] + m_xsi[7][1] + 1.),
		0.03125 * (-9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0] - 1.) },
	  { 0.03125 * (81. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 81. * m_xsi[7][0] * m_xsi[7][0] - 18. * m_xsi[7][0] * m_xsi[7][1] - 18. * m_xsi[7][0] - 27. * m_xsi[7][1] - 27.),
		0.03125 * (27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] + 9.) },
	  { 0.03125 * (-81. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 81. * m_xsi[7][0] * m_xsi[7][0] - 18. * m_xsi[7][0] * m_xsi[7][1] - 18. * m_xsi[7][0] + 27. * m_xsi[7][1] + 27.),
		0.03125 * (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] + 9.) },
	  { 0.03125 * (27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] - m_xsi[7][1] - 1.),
		0.03125 * (9. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 9. * m_xsi[7][0] * m_xsi[7][0] - m_xsi[7][0] - 1.) } } };
