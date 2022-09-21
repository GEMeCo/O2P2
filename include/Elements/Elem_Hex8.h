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
// Hexahedral solid element, with linear interpolation functions.
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
			  * @class Elem_Hex8
			  *
			  * @brief Hexahedral linear element with 8 nodes.
			  * @details Solid element, with linear interpolation functions, hexahedral shape.
			  * Options for integration points: 8, 27 and 64.
			  * @image html Elem_Hex8.png height=300
			  */
			class Elem_Hex8 : public ElementSolid
			{
			private:
				Elem_Hex8() = delete;

			protected:
				/** Constructor for hexahedral linear elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Hex8(std::shared_ptr<O2P2::Prep::Material>& Material)
					: ElementSolid(Material) { }

			public:
				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " "
						<< this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " "
						<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[7]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " "
						<< this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " "
						<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[7]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[2]->m_index + add << " "
						<< this->v_Conect[4]->m_index + add << " " << this->v_Conect[6]->m_index + add << " " << this->m_Mat->m_index << "\n";
					msg << "3 1 " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[3]->m_index + add << " "
						<< this->v_Conect[5]->m_index + add << " " << this->v_Conect[7]->m_index + add << " " << this->m_Mat->m_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "3 1 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (4 + add) << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << (5 + add) << " " << (6 + add) << " " << (7 + add) << " " << (8 + add) << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << (1 + add) << " " << (2 + add) << " " << (5 + add) << " " << (6 + add) << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << (3 + add) << " " << (4 + add) << " " << (7 + add) << " " << (8 + add) << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << (1 + add) << " " << (3 + add) << " " << (5 + add) << " " << (7 + add) << " "
						<< this->m_Mat->m_index << "\n";
					msg << "3 1 " << (2 + add) << " " << (4 + add) << " " << (6 + add) << " " << (8 + add) << " "
						<< this->m_Mat->m_index << "\n";
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
					const auto [min, max] = std::minmax_element(xsi.begin(), xsi.end());
					if (*max < 1.000001 && *min > -1.000001) return true;
					return false;
				}

			protected:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			protected:
				/** @brief Number of Nodes */
				static const int m_NumNodes{ 8 };

				/** @brief Number of Faces */
				static const int m_NumFaces{ 6 };
			};


			/**
			  * @class Elem_Hex8_IP
			  *
			  * @brief Hexahedral linear element with 8 nodes.
			  * @details Solid element, with linear interpolation functions, hexahedral shape.
			  *
			  * @tparam nIP Number of integration points. Must be: 8, 27 or 64.
			  */
			template<int nIP>
			class Elem_Hex8_IP : public Elem_Hex8
			{
			private:
				Elem_Hex8_IP() = delete;

			public:
				/** Constructor for hexahedral linear elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Hex8_IP(std::shared_ptr<O2P2::Prep::Material>& Material)
					: Elem_Hex8(Material) { }

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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Hex8::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(8);

	Psi(0) = 0.125 * (1. - Point[0]) * (1. - Point[1]) * (1. - Point[2]);
	Psi(1) = 0.125 * (1. + Point[0]) * (1. - Point[1]) * (1. - Point[2]);
	Psi(2) = 0.125 * (1. - Point[0]) * (1. + Point[1]) * (1. - Point[2]);
	Psi(3) = 0.125 * (1. + Point[0]) * (1. + Point[1]) * (1. - Point[2]);
	Psi(4) = 0.125 * (1. - Point[0]) * (1. - Point[1]) * (1. + Point[2]);
	Psi(5) = 0.125 * (1. + Point[0]) * (1. - Point[1]) * (1. + Point[2]);
	Psi(6) = 0.125 * (1. - Point[0]) * (1. + Point[1]) * (1. + Point[2]);
	Psi(7) = 0.125 * (1. + Point[0]) * (1. + Point[1]) * (1. + Point[2]);

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd O2P2::Prep::Elem::Elem_Hex8::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(8, 3);

	DPsi(0, 0) = -0.125*(1.-Point[1])*(1.-Point[2]);
	DPsi(1, 0) =  0.125*(1.-Point[1])*(1.-Point[2]);
	DPsi(2, 0) = -0.125*(1.+Point[1])*(1.-Point[2]);
	DPsi(3, 0) =  0.125*(1.+Point[1])*(1.-Point[2]);
	DPsi(4, 0) = -0.125*(1.-Point[1])*(1.+Point[2]);
	DPsi(5, 0) =  0.125*(1.-Point[1])*(1.+Point[2]);
	DPsi(6, 0) = -0.125*(1.+Point[1])*(1.+Point[2]);
	DPsi(7, 0) =  0.125*(1.+Point[1])*(1.+Point[2]);

	DPsi(0, 1) = -0.125*(1.-Point[0])*(1.-Point[2]);
	DPsi(1, 1) = -0.125*(1.+Point[0])*(1.-Point[2]);
	DPsi(2, 1) =  0.125*(1.-Point[0])*(1.-Point[2]);
	DPsi(3, 1) =  0.125*(1.+Point[0])*(1.-Point[2]);
	DPsi(4, 1) = -0.125*(1.-Point[0])*(1.+Point[2]);
	DPsi(5, 1) = -0.125*(1.+Point[0])*(1.+Point[2]);
	DPsi(6, 1) =  0.125*(1.-Point[0])*(1.+Point[2]);
	DPsi(7, 1) =  0.125*(1.+Point[0])*(1.+Point[2]);

	DPsi(0, 2) = -0.125*(1.-Point[0])*(1.-Point[1]);
	DPsi(1, 2) = -0.125*(1.+Point[0])*(1.-Point[1]);
	DPsi(2, 2) = -0.125*(1.-Point[0])*(1.+Point[1]);
	DPsi(3, 2) = -0.125*(1.+Point[0])*(1.+Point[1]);
	DPsi(4, 2) =  0.125*(1.-Point[0])*(1.-Point[1]);
	DPsi(5, 2) =  0.125*(1.+Point[0])*(1.-Point[1]);
	DPsi(6, 2) =  0.125*(1.-Point[0])*(1.+Point[1]);
	DPsi(7, 2) =  0.125*(1.+Point[0])*(1.+Point[1]);

	return DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Prep::Elem::Elem_Hex8::setGeomProperties() {

	const int nVertices = 8;

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
	vertices[6] = v_Conect[6].get();
	vertices[7] = v_Conect[7].get();

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
inline Eigen::VectorXd O2P2::Prep::Elem::Elem_Hex8_IP<nIP>::getValueOnIPs(const double* value) {

	// return value
	Eigen::VectorXd valueOnIp = Eigen::VectorXd::Zero(nIP);

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
template<> const double* O2P2::Prep::Elem::Elem_Hex8_IP<8>::m_weight =  &Gauss3D::Wg_8P[0];
template<> const double* O2P2::Prep::Elem::Elem_Hex8_IP<27>::m_weight = &Gauss3D::Wg_27P[0];
template<> const double* O2P2::Prep::Elem::Elem_Hex8_IP<64>::m_weight = &Gauss3D::Wg_64P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double O2P2::Prep::Elem::Elem_Hex8_IP<8>::m_Psi[8][m_NumNodes] = {
	{ 0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]) } };

template<> const double O2P2::Prep::Elem::Elem_Hex8_IP<27>::m_Psi[27][m_NumNodes] = {
	{ 0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]) } };

template<> const double O2P2::Prep::Elem::Elem_Hex8_IP<64>::m_Psi[64][m_NumNodes] = {
	{ 0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]) },

	{ 0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]),
	  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]) } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
template<> const double O2P2::Prep::Elem::Elem_Hex8_IP<8>::m_DPsi[8][m_NumNodes][m_Dim] = {
	{ { -0.125 * (1. - Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]), -0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][2]), -0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]), -0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][2]), -0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]),  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][2]), -0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[0][1]) * (1. - Gauss3D::Qsi_8P[0][2]),  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][2]), -0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]), -0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][2]),  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]), -0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][2]),  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. - Gauss3D::Qsi_8P[0][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]),  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][2]),  0.125 * (1. - Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[0][1]) * (1. + Gauss3D::Qsi_8P[0][2]),  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][2]),  0.125 * (1. + Gauss3D::Qsi_8P[0][0]) * (1. + Gauss3D::Qsi_8P[0][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]), -0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][2]), -0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]), -0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][2]), -0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]),  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][2]), -0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[1][1]) * (1. - Gauss3D::Qsi_8P[1][2]),  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][2]), -0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]), -0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][2]),  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]), -0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][2]),  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. - Gauss3D::Qsi_8P[1][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]),  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][2]),  0.125 * (1. - Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[1][1]) * (1. + Gauss3D::Qsi_8P[1][2]),  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][2]),  0.125 * (1. + Gauss3D::Qsi_8P[1][0]) * (1. + Gauss3D::Qsi_8P[1][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]), -0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][2]), -0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]), -0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][2]), -0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]),  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][2]), -0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[2][1]) * (1. - Gauss3D::Qsi_8P[2][2]),  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][2]), -0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]), -0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][2]),  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]), -0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][2]),  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. - Gauss3D::Qsi_8P[2][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]),  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][2]),  0.125 * (1. - Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[2][1]) * (1. + Gauss3D::Qsi_8P[2][2]),  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][2]),  0.125 * (1. + Gauss3D::Qsi_8P[2][0]) * (1. + Gauss3D::Qsi_8P[2][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]), -0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][2]), -0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]), -0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][2]), -0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]),  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][2]), -0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[3][1]) * (1. - Gauss3D::Qsi_8P[3][2]),  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][2]), -0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]), -0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][2]),  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]), -0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][2]),  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. - Gauss3D::Qsi_8P[3][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]),  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][2]),  0.125 * (1. - Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[3][1]) * (1. + Gauss3D::Qsi_8P[3][2]),  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][2]),  0.125 * (1. + Gauss3D::Qsi_8P[3][0]) * (1. + Gauss3D::Qsi_8P[3][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]), -0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][2]), -0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]), -0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][2]), -0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]),  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][2]), -0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[4][1]) * (1. - Gauss3D::Qsi_8P[4][2]),  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][2]), -0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]), -0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][2]),  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]), -0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][2]),  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. - Gauss3D::Qsi_8P[4][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]),  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][2]),  0.125 * (1. - Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[4][1]) * (1. + Gauss3D::Qsi_8P[4][2]),  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][2]),  0.125 * (1. + Gauss3D::Qsi_8P[4][0]) * (1. + Gauss3D::Qsi_8P[4][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]), -0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][2]), -0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]), -0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][2]), -0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]),  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][2]), -0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[5][1]) * (1. - Gauss3D::Qsi_8P[5][2]),  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][2]), -0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]), -0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][2]),  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]), -0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][2]),  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. - Gauss3D::Qsi_8P[5][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]),  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][2]),  0.125 * (1. - Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[5][1]) * (1. + Gauss3D::Qsi_8P[5][2]),  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][2]),  0.125 * (1. + Gauss3D::Qsi_8P[5][0]) * (1. + Gauss3D::Qsi_8P[5][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]), -0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][2]), -0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]), -0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][2]), -0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]),  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][2]), -0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[6][1]) * (1. - Gauss3D::Qsi_8P[6][2]),  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][2]), -0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]), -0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][2]),  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]), -0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][2]),  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. - Gauss3D::Qsi_8P[6][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]),  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][2]),  0.125 * (1. - Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[6][1]) * (1. + Gauss3D::Qsi_8P[6][2]),  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][2]),  0.125 * (1. + Gauss3D::Qsi_8P[6][0]) * (1. + Gauss3D::Qsi_8P[6][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]), -0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][2]), -0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]), -0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][2]), -0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]),  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][2]), -0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[7][1]) * (1. - Gauss3D::Qsi_8P[7][2]),  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][2]), -0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]), -0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][2]),  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]), -0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][2]),  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. - Gauss3D::Qsi_8P[7][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]),  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][2]),  0.125 * (1. - Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_8P[7][1]) * (1. + Gauss3D::Qsi_8P[7][2]),  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][2]),  0.125 * (1. + Gauss3D::Qsi_8P[7][0]) * (1. + Gauss3D::Qsi_8P[7][1]) } } };

template<> const double O2P2::Prep::Elem::Elem_Hex8_IP<27>::m_DPsi[27][m_NumNodes][m_Dim] = {
	{ { -0.125 * (1. - Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]), -0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][2]), -0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]), -0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][2]), -0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]),  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][2]), -0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[0][1]) * (1. - Gauss3D::Qsi_27P[0][2]),  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][2]), -0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]), -0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][2]),  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]), -0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][2]),  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. - Gauss3D::Qsi_27P[0][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]),  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][2]),  0.125 * (1. - Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[0][1]) * (1. + Gauss3D::Qsi_27P[0][2]),  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][2]),  0.125 * (1. + Gauss3D::Qsi_27P[0][0]) * (1. + Gauss3D::Qsi_27P[0][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]), -0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][2]), -0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]), -0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][2]), -0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]),  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][2]), -0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[1][1]) * (1. - Gauss3D::Qsi_27P[1][2]),  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][2]), -0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]), -0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][2]),  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]), -0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][2]),  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. - Gauss3D::Qsi_27P[1][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]),  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][2]),  0.125 * (1. - Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[1][1]) * (1. + Gauss3D::Qsi_27P[1][2]),  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][2]),  0.125 * (1. + Gauss3D::Qsi_27P[1][0]) * (1. + Gauss3D::Qsi_27P[1][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]), -0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][2]), -0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]), -0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][2]), -0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]),  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][2]), -0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[2][1]) * (1. - Gauss3D::Qsi_27P[2][2]),  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][2]), -0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]), -0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][2]),  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]), -0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][2]),  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. - Gauss3D::Qsi_27P[2][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]),  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][2]),  0.125 * (1. - Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[2][1]) * (1. + Gauss3D::Qsi_27P[2][2]),  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][2]),  0.125 * (1. + Gauss3D::Qsi_27P[2][0]) * (1. + Gauss3D::Qsi_27P[2][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]), -0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][2]), -0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]), -0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][2]), -0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]),  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][2]), -0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[3][1]) * (1. - Gauss3D::Qsi_27P[3][2]),  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][2]), -0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]), -0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][2]),  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]), -0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][2]),  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. - Gauss3D::Qsi_27P[3][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]),  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][2]),  0.125 * (1. - Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[3][1]) * (1. + Gauss3D::Qsi_27P[3][2]),  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][2]),  0.125 * (1. + Gauss3D::Qsi_27P[3][0]) * (1. + Gauss3D::Qsi_27P[3][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]), -0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][2]), -0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]), -0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][2]), -0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]),  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][2]), -0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[4][1]) * (1. - Gauss3D::Qsi_27P[4][2]),  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][2]), -0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]), -0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][2]),  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]), -0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][2]),  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. - Gauss3D::Qsi_27P[4][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]),  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][2]),  0.125 * (1. - Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[4][1]) * (1. + Gauss3D::Qsi_27P[4][2]),  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][2]),  0.125 * (1. + Gauss3D::Qsi_27P[4][0]) * (1. + Gauss3D::Qsi_27P[4][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]), -0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][2]), -0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]), -0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][2]), -0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]),  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][2]), -0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[5][1]) * (1. - Gauss3D::Qsi_27P[5][2]),  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][2]), -0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]), -0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][2]),  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]), -0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][2]),  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. - Gauss3D::Qsi_27P[5][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]),  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][2]),  0.125 * (1. - Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[5][1]) * (1. + Gauss3D::Qsi_27P[5][2]),  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][2]),  0.125 * (1. + Gauss3D::Qsi_27P[5][0]) * (1. + Gauss3D::Qsi_27P[5][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]), -0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][2]), -0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]), -0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][2]), -0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]),  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][2]), -0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[6][1]) * (1. - Gauss3D::Qsi_27P[6][2]),  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][2]), -0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]), -0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][2]),  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]), -0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][2]),  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. - Gauss3D::Qsi_27P[6][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]),  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][2]),  0.125 * (1. - Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[6][1]) * (1. + Gauss3D::Qsi_27P[6][2]),  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][2]),  0.125 * (1. + Gauss3D::Qsi_27P[6][0]) * (1. + Gauss3D::Qsi_27P[6][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]), -0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][2]), -0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]), -0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][2]), -0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]),  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][2]), -0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[7][1]) * (1. - Gauss3D::Qsi_27P[7][2]),  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][2]), -0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]), -0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][2]),  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]), -0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][2]),  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. - Gauss3D::Qsi_27P[7][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]),  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][2]),  0.125 * (1. - Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[7][1]) * (1. + Gauss3D::Qsi_27P[7][2]),  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][2]),  0.125 * (1. + Gauss3D::Qsi_27P[7][0]) * (1. + Gauss3D::Qsi_27P[7][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]), -0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][2]), -0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]), -0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][2]), -0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]),  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][2]), -0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[8][1]) * (1. - Gauss3D::Qsi_27P[8][2]),  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][2]), -0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]), -0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][2]),  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]), -0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][2]),  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. - Gauss3D::Qsi_27P[8][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]),  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][2]),  0.125 * (1. - Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[8][1]) * (1. + Gauss3D::Qsi_27P[8][2]),  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][2]),  0.125 * (1. + Gauss3D::Qsi_27P[8][0]) * (1. + Gauss3D::Qsi_27P[8][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]), -0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][2]), -0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]), -0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][2]), -0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]),  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][2]), -0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[9][1]) * (1. - Gauss3D::Qsi_27P[9][2]),  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][2]), -0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]), -0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][2]),  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]), -0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][2]),  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. - Gauss3D::Qsi_27P[9][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]),  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][2]),  0.125 * (1. - Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[9][1]) * (1. + Gauss3D::Qsi_27P[9][2]),  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][2]),  0.125 * (1. + Gauss3D::Qsi_27P[9][0]) * (1. + Gauss3D::Qsi_27P[9][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]), -0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][2]), -0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]), -0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][2]), -0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]),  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][2]), -0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[10][1]) * (1. - Gauss3D::Qsi_27P[10][2]),  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][2]), -0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]), -0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][2]),  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]), -0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][2]),  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. - Gauss3D::Qsi_27P[10][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]),  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][2]),  0.125 * (1. - Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[10][1]) * (1. + Gauss3D::Qsi_27P[10][2]),  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][2]),  0.125 * (1. + Gauss3D::Qsi_27P[10][0]) * (1. + Gauss3D::Qsi_27P[10][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]), -0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][2]), -0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]), -0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][2]), -0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]),  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][2]), -0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[11][1]) * (1. - Gauss3D::Qsi_27P[11][2]),  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][2]), -0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]), -0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][2]),  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]), -0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][2]),  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. - Gauss3D::Qsi_27P[11][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]),  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][2]),  0.125 * (1. - Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[11][1]) * (1. + Gauss3D::Qsi_27P[11][2]),  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][2]),  0.125 * (1. + Gauss3D::Qsi_27P[11][0]) * (1. + Gauss3D::Qsi_27P[11][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]), -0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][2]), -0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]), -0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][2]), -0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]),  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][2]), -0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[12][1]) * (1. - Gauss3D::Qsi_27P[12][2]),  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][2]), -0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]), -0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][2]),  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]), -0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][2]),  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. - Gauss3D::Qsi_27P[12][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]),  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][2]),  0.125 * (1. - Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[12][1]) * (1. + Gauss3D::Qsi_27P[12][2]),  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][2]),  0.125 * (1. + Gauss3D::Qsi_27P[12][0]) * (1. + Gauss3D::Qsi_27P[12][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]), -0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][2]), -0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]), -0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][2]), -0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]),  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][2]), -0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[13][1]) * (1. - Gauss3D::Qsi_27P[13][2]),  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][2]), -0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]), -0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][2]),  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]), -0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][2]),  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. - Gauss3D::Qsi_27P[13][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]),  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][2]),  0.125 * (1. - Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[13][1]) * (1. + Gauss3D::Qsi_27P[13][2]),  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][2]),  0.125 * (1. + Gauss3D::Qsi_27P[13][0]) * (1. + Gauss3D::Qsi_27P[13][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]), -0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][2]), -0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]), -0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][2]), -0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]),  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][2]), -0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[14][1]) * (1. - Gauss3D::Qsi_27P[14][2]),  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][2]), -0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]), -0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][2]),  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]), -0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][2]),  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. - Gauss3D::Qsi_27P[14][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]),  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][2]),  0.125 * (1. - Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[14][1]) * (1. + Gauss3D::Qsi_27P[14][2]),  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][2]),  0.125 * (1. + Gauss3D::Qsi_27P[14][0]) * (1. + Gauss3D::Qsi_27P[14][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]), -0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][2]), -0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]), -0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][2]), -0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]),  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][2]), -0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[15][1]) * (1. - Gauss3D::Qsi_27P[15][2]),  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][2]), -0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]), -0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][2]),  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]), -0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][2]),  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. - Gauss3D::Qsi_27P[15][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]),  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][2]),  0.125 * (1. - Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[15][1]) * (1. + Gauss3D::Qsi_27P[15][2]),  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][2]),  0.125 * (1. + Gauss3D::Qsi_27P[15][0]) * (1. + Gauss3D::Qsi_27P[15][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]), -0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][2]), -0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]), -0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][2]), -0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]),  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][2]), -0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[16][1]) * (1. - Gauss3D::Qsi_27P[16][2]),  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][2]), -0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]), -0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][2]),  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]), -0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][2]),  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. - Gauss3D::Qsi_27P[16][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]),  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][2]),  0.125 * (1. - Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[16][1]) * (1. + Gauss3D::Qsi_27P[16][2]),  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][2]),  0.125 * (1. + Gauss3D::Qsi_27P[16][0]) * (1. + Gauss3D::Qsi_27P[16][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]), -0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][2]), -0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]), -0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][2]), -0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]),  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][2]), -0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[17][1]) * (1. - Gauss3D::Qsi_27P[17][2]),  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][2]), -0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]), -0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][2]),  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]), -0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][2]),  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. - Gauss3D::Qsi_27P[17][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]),  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][2]),  0.125 * (1. - Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[17][1]) * (1. + Gauss3D::Qsi_27P[17][2]),  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][2]),  0.125 * (1. + Gauss3D::Qsi_27P[17][0]) * (1. + Gauss3D::Qsi_27P[17][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]), -0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][2]), -0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]), -0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][2]), -0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]),  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][2]), -0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[18][1]) * (1. - Gauss3D::Qsi_27P[18][2]),  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][2]), -0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]), -0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][2]),  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]), -0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][2]),  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. - Gauss3D::Qsi_27P[18][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]),  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][2]),  0.125 * (1. - Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[18][1]) * (1. + Gauss3D::Qsi_27P[18][2]),  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][2]),  0.125 * (1. + Gauss3D::Qsi_27P[18][0]) * (1. + Gauss3D::Qsi_27P[18][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]), -0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][2]), -0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]), -0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][2]), -0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]),  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][2]), -0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[19][1]) * (1. - Gauss3D::Qsi_27P[19][2]),  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][2]), -0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]), -0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][2]),  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]), -0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][2]),  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. - Gauss3D::Qsi_27P[19][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]),  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][2]),  0.125 * (1. - Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[19][1]) * (1. + Gauss3D::Qsi_27P[19][2]),  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][2]),  0.125 * (1. + Gauss3D::Qsi_27P[19][0]) * (1. + Gauss3D::Qsi_27P[19][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]), -0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][2]), -0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]), -0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][2]), -0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]),  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][2]), -0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[20][1]) * (1. - Gauss3D::Qsi_27P[20][2]),  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][2]), -0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]), -0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][2]),  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]), -0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][2]),  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. - Gauss3D::Qsi_27P[20][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]),  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][2]),  0.125 * (1. - Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[20][1]) * (1. + Gauss3D::Qsi_27P[20][2]),  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][2]),  0.125 * (1. + Gauss3D::Qsi_27P[20][0]) * (1. + Gauss3D::Qsi_27P[20][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]), -0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][2]), -0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]), -0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][2]), -0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]),  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][2]), -0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[21][1]) * (1. - Gauss3D::Qsi_27P[21][2]),  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][2]), -0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]), -0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][2]),  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]), -0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][2]),  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. - Gauss3D::Qsi_27P[21][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]),  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][2]),  0.125 * (1. - Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[21][1]) * (1. + Gauss3D::Qsi_27P[21][2]),  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][2]),  0.125 * (1. + Gauss3D::Qsi_27P[21][0]) * (1. + Gauss3D::Qsi_27P[21][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]), -0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][2]), -0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]), -0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][2]), -0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]),  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][2]), -0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[22][1]) * (1. - Gauss3D::Qsi_27P[22][2]),  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][2]), -0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]), -0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][2]),  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]), -0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][2]),  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. - Gauss3D::Qsi_27P[22][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]),  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][2]),  0.125 * (1. - Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[22][1]) * (1. + Gauss3D::Qsi_27P[22][2]),  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][2]),  0.125 * (1. + Gauss3D::Qsi_27P[22][0]) * (1. + Gauss3D::Qsi_27P[22][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]), -0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][2]), -0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]), -0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][2]), -0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]),  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][2]), -0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[23][1]) * (1. - Gauss3D::Qsi_27P[23][2]),  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][2]), -0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]), -0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][2]),  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]), -0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][2]),  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. - Gauss3D::Qsi_27P[23][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]),  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][2]),  0.125 * (1. - Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[23][1]) * (1. + Gauss3D::Qsi_27P[23][2]),  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][2]),  0.125 * (1. + Gauss3D::Qsi_27P[23][0]) * (1. + Gauss3D::Qsi_27P[23][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]), -0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][2]), -0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]), -0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][2]), -0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]),  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][2]), -0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[24][1]) * (1. - Gauss3D::Qsi_27P[24][2]),  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][2]), -0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]), -0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][2]),  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]), -0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][2]),  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. - Gauss3D::Qsi_27P[24][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]),  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][2]),  0.125 * (1. - Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[24][1]) * (1. + Gauss3D::Qsi_27P[24][2]),  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][2]),  0.125 * (1. + Gauss3D::Qsi_27P[24][0]) * (1. + Gauss3D::Qsi_27P[24][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]), -0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][2]), -0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]), -0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][2]), -0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]),  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][2]), -0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[25][1]) * (1. - Gauss3D::Qsi_27P[25][2]),  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][2]), -0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]), -0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][2]),  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]), -0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][2]),  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. - Gauss3D::Qsi_27P[25][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]),  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][2]),  0.125 * (1. - Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[25][1]) * (1. + Gauss3D::Qsi_27P[25][2]),  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][2]),  0.125 * (1. + Gauss3D::Qsi_27P[25][0]) * (1. + Gauss3D::Qsi_27P[25][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]), -0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][2]), -0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]), -0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][2]), -0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]),  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][2]), -0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[26][1]) * (1. - Gauss3D::Qsi_27P[26][2]),  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][2]), -0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]), -0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][2]),  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]), -0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][2]),  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. - Gauss3D::Qsi_27P[26][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]),  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][2]),  0.125 * (1. - Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_27P[26][1]) * (1. + Gauss3D::Qsi_27P[26][2]),  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][2]),  0.125 * (1. + Gauss3D::Qsi_27P[26][0]) * (1. + Gauss3D::Qsi_27P[26][1]) } } };

template<> const double O2P2::Prep::Elem::Elem_Hex8_IP<64>::m_DPsi[64][m_NumNodes][m_Dim] = {
	{ { -0.125 * (1. - Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]), -0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][2]), -0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]), -0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][2]), -0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]),  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][2]), -0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[0][1]) * (1. - Gauss3D::Qsi_64P[0][2]),  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][2]), -0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]), -0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][2]),  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]), -0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][2]),  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. - Gauss3D::Qsi_64P[0][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]),  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][2]),  0.125 * (1. - Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[0][1]) * (1. + Gauss3D::Qsi_64P[0][2]),  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][2]),  0.125 * (1. + Gauss3D::Qsi_64P[0][0]) * (1. + Gauss3D::Qsi_64P[0][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]), -0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][2]), -0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]), -0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][2]), -0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]),  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][2]), -0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[1][1]) * (1. - Gauss3D::Qsi_64P[1][2]),  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][2]), -0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]), -0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][2]),  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]), -0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][2]),  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. - Gauss3D::Qsi_64P[1][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]),  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][2]),  0.125 * (1. - Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[1][1]) * (1. + Gauss3D::Qsi_64P[1][2]),  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][2]),  0.125 * (1. + Gauss3D::Qsi_64P[1][0]) * (1. + Gauss3D::Qsi_64P[1][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]), -0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][2]), -0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]), -0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][2]), -0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]),  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][2]), -0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[2][1]) * (1. - Gauss3D::Qsi_64P[2][2]),  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][2]), -0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]), -0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][2]),  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]), -0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][2]),  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. - Gauss3D::Qsi_64P[2][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]),  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][2]),  0.125 * (1. - Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[2][1]) * (1. + Gauss3D::Qsi_64P[2][2]),  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][2]),  0.125 * (1. + Gauss3D::Qsi_64P[2][0]) * (1. + Gauss3D::Qsi_64P[2][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]), -0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][2]), -0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]), -0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][2]), -0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]),  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][2]), -0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[3][1]) * (1. - Gauss3D::Qsi_64P[3][2]),  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][2]), -0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]), -0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][2]),  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]), -0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][2]),  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. - Gauss3D::Qsi_64P[3][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]),  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][2]),  0.125 * (1. - Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[3][1]) * (1. + Gauss3D::Qsi_64P[3][2]),  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][2]),  0.125 * (1. + Gauss3D::Qsi_64P[3][0]) * (1. + Gauss3D::Qsi_64P[3][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]), -0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][2]), -0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]), -0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][2]), -0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]),  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][2]), -0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[4][1]) * (1. - Gauss3D::Qsi_64P[4][2]),  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][2]), -0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]), -0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][2]),  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]), -0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][2]),  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. - Gauss3D::Qsi_64P[4][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]),  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][2]),  0.125 * (1. - Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[4][1]) * (1. + Gauss3D::Qsi_64P[4][2]),  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][2]),  0.125 * (1. + Gauss3D::Qsi_64P[4][0]) * (1. + Gauss3D::Qsi_64P[4][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]), -0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][2]), -0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]), -0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][2]), -0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]),  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][2]), -0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[5][1]) * (1. - Gauss3D::Qsi_64P[5][2]),  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][2]), -0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]), -0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][2]),  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]), -0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][2]),  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. - Gauss3D::Qsi_64P[5][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]),  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][2]),  0.125 * (1. - Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[5][1]) * (1. + Gauss3D::Qsi_64P[5][2]),  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][2]),  0.125 * (1. + Gauss3D::Qsi_64P[5][0]) * (1. + Gauss3D::Qsi_64P[5][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]), -0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][2]), -0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]), -0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][2]), -0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]),  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][2]), -0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[6][1]) * (1. - Gauss3D::Qsi_64P[6][2]),  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][2]), -0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]), -0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][2]),  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]), -0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][2]),  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. - Gauss3D::Qsi_64P[6][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]),  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][2]),  0.125 * (1. - Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[6][1]) * (1. + Gauss3D::Qsi_64P[6][2]),  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][2]),  0.125 * (1. + Gauss3D::Qsi_64P[6][0]) * (1. + Gauss3D::Qsi_64P[6][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]), -0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][2]), -0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]), -0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][2]), -0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]),  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][2]), -0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[7][1]) * (1. - Gauss3D::Qsi_64P[7][2]),  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][2]), -0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]), -0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][2]),  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]), -0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][2]),  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. - Gauss3D::Qsi_64P[7][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]),  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][2]),  0.125 * (1. - Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[7][1]) * (1. + Gauss3D::Qsi_64P[7][2]),  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][2]),  0.125 * (1. + Gauss3D::Qsi_64P[7][0]) * (1. + Gauss3D::Qsi_64P[7][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]), -0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][2]), -0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]), -0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][2]), -0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]),  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][2]), -0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[8][1]) * (1. - Gauss3D::Qsi_64P[8][2]),  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][2]), -0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]), -0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][2]),  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]), -0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][2]),  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. - Gauss3D::Qsi_64P[8][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]),  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][2]),  0.125 * (1. - Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[8][1]) * (1. + Gauss3D::Qsi_64P[8][2]),  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][2]),  0.125 * (1. + Gauss3D::Qsi_64P[8][0]) * (1. + Gauss3D::Qsi_64P[8][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]), -0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][2]), -0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]), -0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][2]), -0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]),  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][2]), -0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[9][1]) * (1. - Gauss3D::Qsi_64P[9][2]),  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][2]), -0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]), -0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][2]),  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]), -0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][2]),  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. - Gauss3D::Qsi_64P[9][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]),  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][2]),  0.125 * (1. - Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[9][1]) * (1. + Gauss3D::Qsi_64P[9][2]),  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][2]),  0.125 * (1. + Gauss3D::Qsi_64P[9][0]) * (1. + Gauss3D::Qsi_64P[9][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]), -0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][2]), -0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]), -0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][2]), -0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]),  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][2]), -0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[10][1]) * (1. - Gauss3D::Qsi_64P[10][2]),  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][2]), -0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]), -0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][2]),  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]), -0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][2]),  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. - Gauss3D::Qsi_64P[10][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]),  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][2]),  0.125 * (1. - Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[10][1]) * (1. + Gauss3D::Qsi_64P[10][2]),  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][2]),  0.125 * (1. + Gauss3D::Qsi_64P[10][0]) * (1. + Gauss3D::Qsi_64P[10][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]), -0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][2]), -0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]), -0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][2]), -0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]),  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][2]), -0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[11][1]) * (1. - Gauss3D::Qsi_64P[11][2]),  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][2]), -0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]), -0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][2]),  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]), -0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][2]),  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. - Gauss3D::Qsi_64P[11][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]),  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][2]),  0.125 * (1. - Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[11][1]) * (1. + Gauss3D::Qsi_64P[11][2]),  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][2]),  0.125 * (1. + Gauss3D::Qsi_64P[11][0]) * (1. + Gauss3D::Qsi_64P[11][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]), -0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][2]), -0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]), -0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][2]), -0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]),  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][2]), -0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[12][1]) * (1. - Gauss3D::Qsi_64P[12][2]),  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][2]), -0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]), -0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][2]),  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]), -0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][2]),  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. - Gauss3D::Qsi_64P[12][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]),  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][2]),  0.125 * (1. - Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[12][1]) * (1. + Gauss3D::Qsi_64P[12][2]),  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][2]),  0.125 * (1. + Gauss3D::Qsi_64P[12][0]) * (1. + Gauss3D::Qsi_64P[12][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[13][1]) * (1. - Gauss3D::Qsi_64P[13][2]), -0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][2]), -0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[13][1]) * (1. - Gauss3D::Qsi_64P[13][2]), -0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][2]), -0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[13][1]) * (1. - Gauss3D::Qsi_64P[13][2]),  0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][2]), -0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[13][1]) * (1. - Gauss3D::Qsi_64P[13][2]),  0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][2]), -0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[13][1]) * (1. + Gauss3D::Qsi_64P[13][2]), -0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][2]),  0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[13][1]) * (1. + Gauss3D::Qsi_64P[13][2]), -0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][2]),  0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. - Gauss3D::Qsi_64P[13][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[13][1]) * (1. + Gauss3D::Qsi_64P[13][2]),  0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][2]),  0.125 * (1. - Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[13][1]) * (1. + Gauss3D::Qsi_64P[13][2]),  0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][2]),  0.125 * (1. + Gauss3D::Qsi_64P[13][0]) * (1. + Gauss3D::Qsi_64P[13][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]), -0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][2]), -0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]), -0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][2]), -0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]),  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][2]), -0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[14][1]) * (1. - Gauss3D::Qsi_64P[14][2]),  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][2]), -0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]), -0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][2]),  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]), -0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][2]),  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. - Gauss3D::Qsi_64P[14][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]),  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][2]),  0.125 * (1. - Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[14][1]) * (1. + Gauss3D::Qsi_64P[14][2]),  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][2]),  0.125 * (1. + Gauss3D::Qsi_64P[14][0]) * (1. + Gauss3D::Qsi_64P[14][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]), -0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][2]), -0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]), -0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][2]), -0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]),  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][2]), -0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[15][1]) * (1. - Gauss3D::Qsi_64P[15][2]),  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][2]), -0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]), -0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][2]),  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]), -0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][2]),  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. - Gauss3D::Qsi_64P[15][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]),  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][2]),  0.125 * (1. - Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[15][1]) * (1. + Gauss3D::Qsi_64P[15][2]),  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][2]),  0.125 * (1. + Gauss3D::Qsi_64P[15][0]) * (1. + Gauss3D::Qsi_64P[15][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]), -0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][2]), -0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]), -0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][2]), -0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]),  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][2]), -0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[16][1]) * (1. - Gauss3D::Qsi_64P[16][2]),  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][2]), -0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]), -0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][2]),  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]), -0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][2]),  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. - Gauss3D::Qsi_64P[16][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]),  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][2]),  0.125 * (1. - Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[16][1]) * (1. + Gauss3D::Qsi_64P[16][2]),  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][2]),  0.125 * (1. + Gauss3D::Qsi_64P[16][0]) * (1. + Gauss3D::Qsi_64P[16][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]), -0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][2]), -0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]), -0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][2]), -0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]),  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][2]), -0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[17][1]) * (1. - Gauss3D::Qsi_64P[17][2]),  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][2]), -0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]), -0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][2]),  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]), -0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][2]),  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. - Gauss3D::Qsi_64P[17][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]),  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][2]),  0.125 * (1. - Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[17][1]) * (1. + Gauss3D::Qsi_64P[17][2]),  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][2]),  0.125 * (1. + Gauss3D::Qsi_64P[17][0]) * (1. + Gauss3D::Qsi_64P[17][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]), -0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][2]), -0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]), -0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][2]), -0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]),  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][2]), -0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[18][1]) * (1. - Gauss3D::Qsi_64P[18][2]),  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][2]), -0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]), -0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][2]),  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]), -0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][2]),  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. - Gauss3D::Qsi_64P[18][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]),  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][2]),  0.125 * (1. - Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[18][1]) * (1. + Gauss3D::Qsi_64P[18][2]),  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][2]),  0.125 * (1. + Gauss3D::Qsi_64P[18][0]) * (1. + Gauss3D::Qsi_64P[18][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]), -0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][2]), -0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]), -0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][2]), -0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]),  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][2]), -0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[19][1]) * (1. - Gauss3D::Qsi_64P[19][2]),  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][2]), -0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]), -0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][2]),  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]), -0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][2]),  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. - Gauss3D::Qsi_64P[19][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]),  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][2]),  0.125 * (1. - Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[19][1]) * (1. + Gauss3D::Qsi_64P[19][2]),  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][2]),  0.125 * (1. + Gauss3D::Qsi_64P[19][0]) * (1. + Gauss3D::Qsi_64P[19][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]), -0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][2]), -0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]), -0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][2]), -0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]),  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][2]), -0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[20][1]) * (1. - Gauss3D::Qsi_64P[20][2]),  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][2]), -0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]), -0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][2]),  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]), -0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][2]),  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. - Gauss3D::Qsi_64P[20][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]),  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][2]),  0.125 * (1. - Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[20][1]) * (1. + Gauss3D::Qsi_64P[20][2]),  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][2]),  0.125 * (1. + Gauss3D::Qsi_64P[20][0]) * (1. + Gauss3D::Qsi_64P[20][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]), -0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][2]), -0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]), -0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][2]), -0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]),  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][2]), -0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[21][1]) * (1. - Gauss3D::Qsi_64P[21][2]),  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][2]), -0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]), -0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][2]),  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]), -0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][2]),  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. - Gauss3D::Qsi_64P[21][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]),  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][2]),  0.125 * (1. - Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[21][1]) * (1. + Gauss3D::Qsi_64P[21][2]),  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][2]),  0.125 * (1. + Gauss3D::Qsi_64P[21][0]) * (1. + Gauss3D::Qsi_64P[21][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]), -0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][2]), -0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]), -0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][2]), -0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]),  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][2]), -0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[22][1]) * (1. - Gauss3D::Qsi_64P[22][2]),  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][2]), -0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]), -0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][2]),  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]), -0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][2]),  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. - Gauss3D::Qsi_64P[22][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]),  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][2]),  0.125 * (1. - Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[22][1]) * (1. + Gauss3D::Qsi_64P[22][2]),  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][2]),  0.125 * (1. + Gauss3D::Qsi_64P[22][0]) * (1. + Gauss3D::Qsi_64P[22][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]), -0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][2]), -0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]), -0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][2]), -0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]),  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][2]), -0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[23][1]) * (1. - Gauss3D::Qsi_64P[23][2]),  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][2]), -0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]), -0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][2]),  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]), -0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][2]),  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. - Gauss3D::Qsi_64P[23][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]),  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][2]),  0.125 * (1. - Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[23][1]) * (1. + Gauss3D::Qsi_64P[23][2]),  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][2]),  0.125 * (1. + Gauss3D::Qsi_64P[23][0]) * (1. + Gauss3D::Qsi_64P[23][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]), -0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][2]), -0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]), -0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][2]), -0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]),  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][2]), -0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[24][1]) * (1. - Gauss3D::Qsi_64P[24][2]),  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][2]), -0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]), -0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][2]),  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]), -0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][2]),  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. - Gauss3D::Qsi_64P[24][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]),  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][2]),  0.125 * (1. - Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[24][1]) * (1. + Gauss3D::Qsi_64P[24][2]),  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][2]),  0.125 * (1. + Gauss3D::Qsi_64P[24][0]) * (1. + Gauss3D::Qsi_64P[24][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]), -0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][2]), -0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]), -0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][2]), -0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]),  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][2]), -0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[25][1]) * (1. - Gauss3D::Qsi_64P[25][2]),  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][2]), -0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]), -0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][2]),  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]), -0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][2]),  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. - Gauss3D::Qsi_64P[25][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]),  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][2]),  0.125 * (1. - Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[25][1]) * (1. + Gauss3D::Qsi_64P[25][2]),  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][2]),  0.125 * (1. + Gauss3D::Qsi_64P[25][0]) * (1. + Gauss3D::Qsi_64P[25][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]), -0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][2]), -0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]), -0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][2]), -0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]),  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][2]), -0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[26][1]) * (1. - Gauss3D::Qsi_64P[26][2]),  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][2]), -0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]), -0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][2]),  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]), -0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][2]),  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. - Gauss3D::Qsi_64P[26][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]),  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][2]),  0.125 * (1. - Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[26][1]) * (1. + Gauss3D::Qsi_64P[26][2]),  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][2]),  0.125 * (1. + Gauss3D::Qsi_64P[26][0]) * (1. + Gauss3D::Qsi_64P[26][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]), -0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][2]), -0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]), -0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][2]), -0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]),  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][2]), -0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[27][1]) * (1. - Gauss3D::Qsi_64P[27][2]),  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][2]), -0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]), -0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][2]),  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]), -0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][2]),  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. - Gauss3D::Qsi_64P[27][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]),  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][2]),  0.125 * (1. - Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[27][1]) * (1. + Gauss3D::Qsi_64P[27][2]),  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][2]),  0.125 * (1. + Gauss3D::Qsi_64P[27][0]) * (1. + Gauss3D::Qsi_64P[27][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]), -0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][2]), -0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]), -0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][2]), -0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]),  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][2]), -0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[28][1]) * (1. - Gauss3D::Qsi_64P[28][2]),  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][2]), -0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]), -0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][2]),  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]), -0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][2]),  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. - Gauss3D::Qsi_64P[28][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]),  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][2]),  0.125 * (1. - Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[28][1]) * (1. + Gauss3D::Qsi_64P[28][2]),  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][2]),  0.125 * (1. + Gauss3D::Qsi_64P[28][0]) * (1. + Gauss3D::Qsi_64P[28][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]), -0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][2]), -0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]), -0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][2]), -0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]),  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][2]), -0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[29][1]) * (1. - Gauss3D::Qsi_64P[29][2]),  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][2]), -0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]), -0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][2]),  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]), -0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][2]),  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. - Gauss3D::Qsi_64P[29][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]),  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][2]),  0.125 * (1. - Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[29][1]) * (1. + Gauss3D::Qsi_64P[29][2]),  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][2]),  0.125 * (1. + Gauss3D::Qsi_64P[29][0]) * (1. + Gauss3D::Qsi_64P[29][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]), -0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][2]), -0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]), -0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][2]), -0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]),  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][2]), -0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[30][1]) * (1. - Gauss3D::Qsi_64P[30][2]),  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][2]), -0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]), -0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][2]),  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]), -0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][2]),  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. - Gauss3D::Qsi_64P[30][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]),  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][2]),  0.125 * (1. - Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[30][1]) * (1. + Gauss3D::Qsi_64P[30][2]),  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][2]),  0.125 * (1. + Gauss3D::Qsi_64P[30][0]) * (1. + Gauss3D::Qsi_64P[30][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]), -0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][2]), -0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]), -0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][2]), -0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]),  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][2]), -0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[31][1]) * (1. - Gauss3D::Qsi_64P[31][2]),  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][2]), -0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]), -0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][2]),  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]), -0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][2]),  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. - Gauss3D::Qsi_64P[31][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]),  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][2]),  0.125 * (1. - Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[31][1]) * (1. + Gauss3D::Qsi_64P[31][2]),  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][2]),  0.125 * (1. + Gauss3D::Qsi_64P[31][0]) * (1. + Gauss3D::Qsi_64P[31][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]), -0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][2]), -0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]), -0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][2]), -0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]),  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][2]), -0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[32][1]) * (1. - Gauss3D::Qsi_64P[32][2]),  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][2]), -0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]), -0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][2]),  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]), -0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][2]),  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. - Gauss3D::Qsi_64P[32][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]),  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][2]),  0.125 * (1. - Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[32][1]) * (1. + Gauss3D::Qsi_64P[32][2]),  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][2]),  0.125 * (1. + Gauss3D::Qsi_64P[32][0]) * (1. + Gauss3D::Qsi_64P[32][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]), -0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][2]), -0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]), -0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][2]), -0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]),  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][2]), -0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[33][1]) * (1. - Gauss3D::Qsi_64P[33][2]),  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][2]), -0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]), -0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][2]),  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]), -0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][2]),  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. - Gauss3D::Qsi_64P[33][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]),  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][2]),  0.125 * (1. - Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[33][1]) * (1. + Gauss3D::Qsi_64P[33][2]),  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][2]),  0.125 * (1. + Gauss3D::Qsi_64P[33][0]) * (1. + Gauss3D::Qsi_64P[33][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]), -0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][2]), -0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]), -0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][2]), -0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]),  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][2]), -0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[34][1]) * (1. - Gauss3D::Qsi_64P[34][2]),  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][2]), -0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]), -0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][2]),  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]), -0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][2]),  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. - Gauss3D::Qsi_64P[34][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]),  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][2]),  0.125 * (1. - Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[34][1]) * (1. + Gauss3D::Qsi_64P[34][2]),  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][2]),  0.125 * (1. + Gauss3D::Qsi_64P[34][0]) * (1. + Gauss3D::Qsi_64P[34][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]), -0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][2]), -0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]), -0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][2]), -0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]),  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][2]), -0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[35][1]) * (1. - Gauss3D::Qsi_64P[35][2]),  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][2]), -0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]), -0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][2]),  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]), -0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][2]),  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. - Gauss3D::Qsi_64P[35][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]),  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][2]),  0.125 * (1. - Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[35][1]) * (1. + Gauss3D::Qsi_64P[35][2]),  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][2]),  0.125 * (1. + Gauss3D::Qsi_64P[35][0]) * (1. + Gauss3D::Qsi_64P[35][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]), -0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][2]), -0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]), -0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][2]), -0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]),  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][2]), -0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[36][1]) * (1. - Gauss3D::Qsi_64P[36][2]),  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][2]), -0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]), -0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][2]),  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]), -0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][2]),  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. - Gauss3D::Qsi_64P[36][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]),  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][2]),  0.125 * (1. - Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[36][1]) * (1. + Gauss3D::Qsi_64P[36][2]),  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][2]),  0.125 * (1. + Gauss3D::Qsi_64P[36][0]) * (1. + Gauss3D::Qsi_64P[36][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]), -0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][2]), -0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]), -0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][2]), -0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]),  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][2]), -0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[37][1]) * (1. - Gauss3D::Qsi_64P[37][2]),  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][2]), -0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]), -0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][2]),  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]), -0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][2]),  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. - Gauss3D::Qsi_64P[37][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]),  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][2]),  0.125 * (1. - Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[37][1]) * (1. + Gauss3D::Qsi_64P[37][2]),  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][2]),  0.125 * (1. + Gauss3D::Qsi_64P[37][0]) * (1. + Gauss3D::Qsi_64P[37][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]), -0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][2]), -0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]), -0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][2]), -0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]),  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][2]), -0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[38][1]) * (1. - Gauss3D::Qsi_64P[38][2]),  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][2]), -0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]), -0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][2]),  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]), -0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][2]),  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. - Gauss3D::Qsi_64P[38][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]),  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][2]),  0.125 * (1. - Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[38][1]) * (1. + Gauss3D::Qsi_64P[38][2]),  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][2]),  0.125 * (1. + Gauss3D::Qsi_64P[38][0]) * (1. + Gauss3D::Qsi_64P[38][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]), -0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][2]), -0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]), -0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][2]), -0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]),  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][2]), -0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[39][1]) * (1. - Gauss3D::Qsi_64P[39][2]),  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][2]), -0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]), -0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][2]),  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]), -0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][2]),  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. - Gauss3D::Qsi_64P[39][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]),  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][2]),  0.125 * (1. - Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[39][1]) * (1. + Gauss3D::Qsi_64P[39][2]),  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][2]),  0.125 * (1. + Gauss3D::Qsi_64P[39][0]) * (1. + Gauss3D::Qsi_64P[39][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]), -0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][2]), -0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]), -0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][2]), -0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][2]), -0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[40][1]) * (1. - Gauss3D::Qsi_64P[40][2]),  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][2]), -0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]), -0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][2]),  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]), -0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][2]),  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. - Gauss3D::Qsi_64P[40][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][2]),  0.125 * (1. - Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[40][1]) * (1. + Gauss3D::Qsi_64P[40][2]),  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][2]),  0.125 * (1. + Gauss3D::Qsi_64P[40][0]) * (1. + Gauss3D::Qsi_64P[40][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]), -0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][2]), -0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]), -0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][2]), -0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]),  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][2]), -0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[41][1]) * (1. - Gauss3D::Qsi_64P[41][2]),  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][2]), -0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]), -0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][2]),  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]), -0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][2]),  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. - Gauss3D::Qsi_64P[41][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]),  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][2]),  0.125 * (1. - Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[41][1]) * (1. + Gauss3D::Qsi_64P[41][2]),  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][2]),  0.125 * (1. + Gauss3D::Qsi_64P[41][0]) * (1. + Gauss3D::Qsi_64P[41][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]), -0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][2]), -0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]), -0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][2]), -0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]),  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][2]), -0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[42][1]) * (1. - Gauss3D::Qsi_64P[42][2]),  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][2]), -0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]), -0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][2]),  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]), -0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][2]),  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. - Gauss3D::Qsi_64P[42][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]),  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][2]),  0.125 * (1. - Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[42][1]) * (1. + Gauss3D::Qsi_64P[42][2]),  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][2]),  0.125 * (1. + Gauss3D::Qsi_64P[42][0]) * (1. + Gauss3D::Qsi_64P[42][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]), -0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][2]), -0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]), -0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][2]), -0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]),  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][2]), -0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[43][1]) * (1. - Gauss3D::Qsi_64P[43][2]),  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][2]), -0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]), -0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][2]),  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]), -0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][2]),  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. - Gauss3D::Qsi_64P[43][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]),  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][2]),  0.125 * (1. - Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[43][1]) * (1. + Gauss3D::Qsi_64P[43][2]),  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][2]),  0.125 * (1. + Gauss3D::Qsi_64P[43][0]) * (1. + Gauss3D::Qsi_64P[43][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]), -0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][2]), -0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]), -0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][2]), -0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]),  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][2]), -0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[44][1]) * (1. - Gauss3D::Qsi_64P[44][2]),  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][2]), -0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]), -0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][2]),  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]), -0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][2]),  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. - Gauss3D::Qsi_64P[44][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]),  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][2]),  0.125 * (1. - Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[44][1]) * (1. + Gauss3D::Qsi_64P[44][2]),  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][2]),  0.125 * (1. + Gauss3D::Qsi_64P[44][0]) * (1. + Gauss3D::Qsi_64P[44][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]), -0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][2]), -0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]), -0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][2]), -0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]),  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][2]), -0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[45][1]) * (1. - Gauss3D::Qsi_64P[45][2]),  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][2]), -0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]), -0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][2]),  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]), -0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][2]),  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. - Gauss3D::Qsi_64P[45][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]),  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][2]),  0.125 * (1. - Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[45][1]) * (1. + Gauss3D::Qsi_64P[45][2]),  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][2]),  0.125 * (1. + Gauss3D::Qsi_64P[45][0]) * (1. + Gauss3D::Qsi_64P[45][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]), -0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][2]), -0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]), -0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][2]), -0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]),  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][2]), -0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[46][1]) * (1. - Gauss3D::Qsi_64P[46][2]),  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][2]), -0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]), -0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][2]),  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]), -0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][2]),  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. - Gauss3D::Qsi_64P[46][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]),  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][2]),  0.125 * (1. - Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[46][1]) * (1. + Gauss3D::Qsi_64P[46][2]),  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][2]),  0.125 * (1. + Gauss3D::Qsi_64P[46][0]) * (1. + Gauss3D::Qsi_64P[46][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]), -0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][2]), -0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]), -0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][2]), -0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]),  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][2]), -0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[47][1]) * (1. - Gauss3D::Qsi_64P[47][2]),  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][2]), -0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]), -0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][2]),  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]), -0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][2]),  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. - Gauss3D::Qsi_64P[47][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]),  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][2]),  0.125 * (1. - Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[47][1]) * (1. + Gauss3D::Qsi_64P[47][2]),  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][2]),  0.125 * (1. + Gauss3D::Qsi_64P[47][0]) * (1. + Gauss3D::Qsi_64P[47][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]), -0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][2]), -0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]), -0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][2]), -0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]),  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][2]), -0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[48][1]) * (1. - Gauss3D::Qsi_64P[48][2]),  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][2]), -0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]), -0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][2]),  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]), -0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][2]),  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. - Gauss3D::Qsi_64P[48][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]),  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][2]),  0.125 * (1. - Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[48][1]) * (1. + Gauss3D::Qsi_64P[48][2]),  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][2]),  0.125 * (1. + Gauss3D::Qsi_64P[48][0]) * (1. + Gauss3D::Qsi_64P[48][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]), -0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][2]), -0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]), -0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][2]), -0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]),  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][2]), -0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[49][1]) * (1. - Gauss3D::Qsi_64P[49][2]),  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][2]), -0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]), -0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][2]),  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]), -0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][2]),  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. - Gauss3D::Qsi_64P[49][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]),  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][2]),  0.125 * (1. - Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[49][1]) * (1. + Gauss3D::Qsi_64P[49][2]),  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][2]),  0.125 * (1. + Gauss3D::Qsi_64P[49][0]) * (1. + Gauss3D::Qsi_64P[49][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]), -0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][2]), -0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]), -0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][2]), -0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]),  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][2]), -0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[50][1]) * (1. - Gauss3D::Qsi_64P[50][2]),  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][2]), -0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]), -0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][2]),  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]), -0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][2]),  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. - Gauss3D::Qsi_64P[50][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]),  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][2]),  0.125 * (1. - Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[50][1]) * (1. + Gauss3D::Qsi_64P[50][2]),  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][2]),  0.125 * (1. + Gauss3D::Qsi_64P[50][0]) * (1. + Gauss3D::Qsi_64P[50][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]), -0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][2]), -0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]), -0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][2]), -0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]),  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][2]), -0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[51][1]) * (1. - Gauss3D::Qsi_64P[51][2]),  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][2]), -0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]), -0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][2]),  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]), -0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][2]),  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. - Gauss3D::Qsi_64P[51][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]),  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][2]),  0.125 * (1. - Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[51][1]) * (1. + Gauss3D::Qsi_64P[51][2]),  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][2]),  0.125 * (1. + Gauss3D::Qsi_64P[51][0]) * (1. + Gauss3D::Qsi_64P[51][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]), -0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][2]), -0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]), -0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][2]), -0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]),  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][2]), -0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[52][1]) * (1. - Gauss3D::Qsi_64P[52][2]),  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][2]), -0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]), -0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][2]),  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]), -0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][2]),  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. - Gauss3D::Qsi_64P[52][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]),  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][2]),  0.125 * (1. - Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[52][1]) * (1. + Gauss3D::Qsi_64P[52][2]),  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][2]),  0.125 * (1. + Gauss3D::Qsi_64P[52][0]) * (1. + Gauss3D::Qsi_64P[52][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]), -0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][2]), -0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]), -0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][2]), -0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]),  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][2]), -0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[53][1]) * (1. - Gauss3D::Qsi_64P[53][2]),  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][2]), -0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]), -0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][2]),  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]), -0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][2]),  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. - Gauss3D::Qsi_64P[53][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]),  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][2]),  0.125 * (1. - Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[53][1]) * (1. + Gauss3D::Qsi_64P[53][2]),  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][2]),  0.125 * (1. + Gauss3D::Qsi_64P[53][0]) * (1. + Gauss3D::Qsi_64P[53][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]), -0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][2]), -0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]), -0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][2]), -0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]),  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][2]), -0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[54][1]) * (1. - Gauss3D::Qsi_64P[54][2]),  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][2]), -0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]), -0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][2]),  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]), -0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][2]),  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. - Gauss3D::Qsi_64P[54][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]),  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][2]),  0.125 * (1. - Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[54][1]) * (1. + Gauss3D::Qsi_64P[54][2]),  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][2]),  0.125 * (1. + Gauss3D::Qsi_64P[54][0]) * (1. + Gauss3D::Qsi_64P[54][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]), -0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][2]), -0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]), -0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][2]), -0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]),  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][2]), -0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[55][1]) * (1. - Gauss3D::Qsi_64P[55][2]),  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][2]), -0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]), -0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][2]),  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]), -0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][2]),  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. - Gauss3D::Qsi_64P[55][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]),  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][2]),  0.125 * (1. - Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[55][1]) * (1. + Gauss3D::Qsi_64P[55][2]),  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][2]),  0.125 * (1. + Gauss3D::Qsi_64P[55][0]) * (1. + Gauss3D::Qsi_64P[55][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]), -0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][2]), -0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]), -0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][2]), -0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]),  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][2]), -0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[56][1]) * (1. - Gauss3D::Qsi_64P[56][2]),  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][2]), -0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]), -0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][2]),  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]), -0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][2]),  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. - Gauss3D::Qsi_64P[56][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]),  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][2]),  0.125 * (1. - Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[56][1]) * (1. + Gauss3D::Qsi_64P[56][2]),  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][2]),  0.125 * (1. + Gauss3D::Qsi_64P[56][0]) * (1. + Gauss3D::Qsi_64P[56][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]), -0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][2]), -0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]), -0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][2]), -0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]),  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][2]), -0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[57][1]) * (1. - Gauss3D::Qsi_64P[57][2]),  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][2]), -0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]), -0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][2]),  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]), -0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][2]),  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. - Gauss3D::Qsi_64P[57][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]),  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][2]),  0.125 * (1. - Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[57][1]) * (1. + Gauss3D::Qsi_64P[57][2]),  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][2]),  0.125 * (1. + Gauss3D::Qsi_64P[57][0]) * (1. + Gauss3D::Qsi_64P[57][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]), -0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][2]), -0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]), -0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][2]), -0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]),  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][2]), -0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[58][1]) * (1. - Gauss3D::Qsi_64P[58][2]),  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][2]), -0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]), -0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][2]),  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]), -0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][2]),  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. - Gauss3D::Qsi_64P[58][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]),  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][2]),  0.125 * (1. - Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[58][1]) * (1. + Gauss3D::Qsi_64P[58][2]),  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][2]),  0.125 * (1. + Gauss3D::Qsi_64P[58][0]) * (1. + Gauss3D::Qsi_64P[58][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]), -0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][2]), -0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]), -0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][2]), -0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]),  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][2]), -0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[59][1]) * (1. - Gauss3D::Qsi_64P[59][2]),  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][2]), -0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]), -0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][2]),  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]), -0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][2]),  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. - Gauss3D::Qsi_64P[59][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]),  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][2]),  0.125 * (1. - Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[59][1]) * (1. + Gauss3D::Qsi_64P[59][2]),  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][2]),  0.125 * (1. + Gauss3D::Qsi_64P[59][0]) * (1. + Gauss3D::Qsi_64P[59][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]), -0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][2]), -0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]), -0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][2]), -0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]),  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][2]), -0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[60][1]) * (1. - Gauss3D::Qsi_64P[60][2]),  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][2]), -0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]), -0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][2]),  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]), -0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][2]),  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. - Gauss3D::Qsi_64P[60][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]),  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][2]),  0.125 * (1. - Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[60][1]) * (1. + Gauss3D::Qsi_64P[60][2]),  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][2]),  0.125 * (1. + Gauss3D::Qsi_64P[60][0]) * (1. + Gauss3D::Qsi_64P[60][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]), -0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][2]), -0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]), -0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][2]), -0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]),  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][2]), -0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[61][1]) * (1. - Gauss3D::Qsi_64P[61][2]),  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][2]), -0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]), -0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][2]),  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]), -0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][2]),  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. - Gauss3D::Qsi_64P[61][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]),  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][2]),  0.125 * (1. - Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[61][1]) * (1. + Gauss3D::Qsi_64P[61][2]),  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][2]),  0.125 * (1. + Gauss3D::Qsi_64P[61][0]) * (1. + Gauss3D::Qsi_64P[61][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]), -0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][2]), -0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]), -0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][2]), -0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]),  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][2]), -0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[62][1]) * (1. - Gauss3D::Qsi_64P[62][2]),  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][2]), -0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]), -0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][2]),  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]), -0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][2]),  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. - Gauss3D::Qsi_64P[62][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]),  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][2]),  0.125 * (1. - Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[62][1]) * (1. + Gauss3D::Qsi_64P[62][2]),  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][2]),  0.125 * (1. + Gauss3D::Qsi_64P[62][0]) * (1. + Gauss3D::Qsi_64P[62][1]) } },

	{ { -0.125 * (1. - Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]), -0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][2]), -0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]), -0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][2]), -0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]),  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][2]), -0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[63][1]) * (1. - Gauss3D::Qsi_64P[63][2]),  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][2]), -0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) },
	  { -0.125 * (1. - Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]), -0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][2]),  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) },
	  {  0.125 * (1. - Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]), -0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][2]),  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. - Gauss3D::Qsi_64P[63][1]) },
	  { -0.125 * (1. + Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]),  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][2]),  0.125 * (1. - Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) },
	  {  0.125 * (1. + Gauss3D::Qsi_64P[63][1]) * (1. + Gauss3D::Qsi_64P[63][2]),  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][2]),  0.125 * (1. + Gauss3D::Qsi_64P[63][0]) * (1. + Gauss3D::Qsi_64P[63][1]) } } };
