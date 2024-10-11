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
// Plane element, with cubic interpolation functions, triangular shaped.
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
			  * @class Elem_Tri10
			  *
			  * @brief Triangular cubic element with 10 nodes.
			  * @details Plane element, with cubic interpolation functions, triangular shaped.
			  * Options for integration points: 6, 7, 12 and 13.
			  * @image html Elem_Tri10.png height=300
			  * @note Minimum number of integration points: 6.
			  */
			class Elem_Tri10 : public ElemPlane
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Tri10() = delete;

			protected:
				/** Constructor for triangular cubic elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to PlaneSection class.
				  */
				explicit Elem_Tri10(std::shared_ptr<O2P2::Geom::Material>& Material, std::shared_ptr<O2P2::Geom::PlaneSection>& Section)
					: ElemPlane(Material, Section) { }

			public:
				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 3 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " "
						<< this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " "
						<< this->mv_Conect[4]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " "
						<< this->mv_Conect[6]->mv_index + add << " " << this->mv_Conect[7]->mv_index + add << " "
						<< this->mv_Conect[8]->mv_index + add << " " << this->mv_Conect[9]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "2 3 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " "
						<< (4 + add) << " " << (5 + add) << " " << (6 + add) << " "
						<< (7 + add) << " " << (8 + add) << " " << (9 + add) << " "
						<< (10 + add) << " " << this->mv_Mat->mv_index << "\n";
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

			protected:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			protected:
				/** @brief Number of Nodes */
				static const int mv_numNodes{ 10 };

				/** @brief Number of Faces */
				static const int mv_numFaces{ 1 };
			};


			/**
			  * @class Elem_Tri10_IP
			  *
			  * @brief Triangular cubic element with 10 nodes.
			  * @details Plane element, with cubic interpolation functions, triangular shaped.
			  *
			  * @tparam nIP Number of integration points. Must be: 6, 7, 12 or 13.
			  */
			template<int nIP>
			class Elem_Tri10_IP : public Elem_Tri10
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Tri10_IP() = delete;

			public:
				/** Constructor for triangular cubic elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to PlaneSection class.
				  */
				explicit Elem_Tri10_IP(std::shared_ptr<O2P2::Geom::Material>& Material, std::shared_ptr<O2P2::Geom::PlaneSection>& Section)
					: Elem_Tri10(Material, Section) { }

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
inline std::vector<double> O2P2::Geom::Elem::Elem_Tri10::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(10);

	mi_Psi.at(0) = 1. - 5.5 * Point[0] - 5.5 * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1];
	mi_Psi.at(1) = 9. * Point[0] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1];
	mi_Psi.at(2) = -4.5 * Point[0] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1];
	mi_Psi.at(3) = Point[0] - 4.5 * Point[0] * Point[0] + 4.5 * Point[0] * Point[0] * Point[0];
	mi_Psi.at(4) = 9. * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1];
	mi_Psi.at(5) = 27. * Point[0] * Point[1] - 27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1];
	mi_Psi.at(6) = -4.5 * Point[0] * Point[1] + 13.5 * Point[0] * Point[0] * Point[1];
	mi_Psi.at(7) = -4.5 * Point[1] + 18. * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1];
	mi_Psi.at(8) = -4.5 * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1];
	mi_Psi.at(9) = Point[1] - 4.5 * Point[1] * Point[1] + 4.5 * Point[1] * Point[1] * Point[1];

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Geom::Elem::Elem_Tri10::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(10 * 2);

	mi_DPsi.at(0) = -5.5 + 18. * Point[0] + 18. * Point[1] - 13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1];
	mi_DPsi.at(1) = 9. - 45. * Point[0] - 22.5 * Point[1] + 40.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 13.5 * Point[1] * Point[1];
	mi_DPsi.at(2) = -4.5 + 36. * Point[0] + 4.5 * Point[1] - 40.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1];
	mi_DPsi.at(3) = 1. - 9. * Point[0] + 13.5 * Point[0] * Point[0];
	mi_DPsi.at(4) = -22.5 * Point[1] + 27. * Point[0] * Point[1] + 27. * Point[1] * Point[1];
	mi_DPsi.at(5) = 27. * Point[1] - 54. * Point[0] * Point[1] - 27. * Point[1] * Point[1];
	mi_DPsi.at(6) = -4.5 * Point[1] + 27. * Point[0] * Point[1];
	mi_DPsi.at(7) = 4.5 * Point[1] - 13.5 * Point[1] * Point[1];
	mi_DPsi.at(8) = -4.5 * Point[1] + 13.5 * Point[1] * Point[1];
	mi_DPsi.at(9) = 0.;

	mi_DPsi.at(10) = -5.5 + 18. * Point[0] + 18. * Point[1] - 13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1];
	mi_DPsi.at(11) = -22.5 * Point[0] + 27. * Point[0] * Point[0] + 27. * Point[0] * Point[1];
	mi_DPsi.at(12) = 4.5 * Point[0] - 13.5 * Point[0] * Point[0];
	mi_DPsi.at(13) = 0.;
	mi_DPsi.at(14) = 9. - 22.5 * Point[0] - 45. * Point[1] + 13.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 40.5 * Point[1] * Point[1];
	mi_DPsi.at(15) = 27. * Point[0] - 27. * Point[0] * Point[0] - 54. * Point[0] * Point[1];
	mi_DPsi.at(16) = -4.5 * Point[0] + 13.5 * Point[0] * Point[0];
	mi_DPsi.at(17) = -4.5 + 36. * Point[1] + 4.5 * Point[0] - 27. * Point[0] * Point[1] - 40.5 * Point[1] * Point[1];
	mi_DPsi.at(18) = -4.5 * Point[0] + 27. * Point[0] * Point[1];
	mi_DPsi.at(19) = 1. - 9. * Point[1] + 13.5 * Point[1] * Point[1];

	return mi_DPsi;
};

// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Geom::Elem::Elem_Tri10::setGeomProperties() {

	const int nVertices = 3;
	const int mi_Dim = mv_Conect.at(0)->getDIM();	// Dimensionality of vector space (2D or 3D)

	mv_Centroid = std::make_unique<double[]>(mi_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Geom::Node*, nVertices> vertices;
	vertices[0] = mv_Conect[0].get();
	vertices[1] = mv_Conect[3].get();
	vertices[2] = mv_Conect[9].get();

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
inline std::vector<double> O2P2::Geom::Elem::Elem_Tri10_IP<nIP>::getValueOnIPs(const double* value) {

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
template<> const double* O2P2::Geom::Elem::Elem_Tri10_IP<6>::mv_weight = &Hammer2D::Wg_6P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tri10_IP<7>::mv_weight = &Hammer2D::Wg_7P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tri10_IP<12>::mv_weight = &Hammer2D::Wg_12P[0];
template<> const double* O2P2::Geom::Elem::Elem_Tri10_IP<13>::mv_weight = &Hammer2D::Wg_13P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<6>::mv_Psi[6][mv_numNodes] = {
	{ 1. - 5.5 * Hammer2D::Qsi_6P[0][0] - 5.5 * Hammer2D::Qsi_6P[0][1] + 9. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 18. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 9. * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] - 4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] - 4.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] ,
	  9. * Hammer2D::Qsi_6P[0][0] - 22.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 22.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] ,
	  -4.5 * Hammer2D::Qsi_6P[0][0] + 18. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] ,
	  Hammer2D::Qsi_6P[0][0] - 4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] ,
	  9. * Hammer2D::Qsi_6P[0][1] - 22.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 22.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] ,
	  27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] ,
	  -4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] ,
	  -4.5 * Hammer2D::Qsi_6P[0][1] + 18. * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] + 4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] ,
	  -4.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] ,
	  Hammer2D::Qsi_6P[0][1] - 4.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] + 4.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_6P[1][0] - 5.5 * Hammer2D::Qsi_6P[1][1] + 9. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 18. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 9. * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] - 4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] - 4.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] ,
	  9. * Hammer2D::Qsi_6P[1][0] - 22.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 22.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] ,
	  -4.5 * Hammer2D::Qsi_6P[1][0] + 18. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] ,
	  Hammer2D::Qsi_6P[1][0] - 4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] ,
	  9. * Hammer2D::Qsi_6P[1][1] - 22.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 22.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] ,
	  27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] ,
	  -4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] ,
	  -4.5 * Hammer2D::Qsi_6P[1][1] + 18. * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] + 4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] ,
	  -4.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] ,
	  Hammer2D::Qsi_6P[1][1] - 4.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] + 4.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_6P[2][0] - 5.5 * Hammer2D::Qsi_6P[2][1] + 9. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 18. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 9. * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] - 4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] - 4.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] ,
	  9. * Hammer2D::Qsi_6P[2][0] - 22.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 22.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] ,
	  -4.5 * Hammer2D::Qsi_6P[2][0] + 18. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] ,
	  Hammer2D::Qsi_6P[2][0] - 4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] ,
	  9. * Hammer2D::Qsi_6P[2][1] - 22.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 22.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] ,
	  27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] ,
	  -4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] ,
	  -4.5 * Hammer2D::Qsi_6P[2][1] + 18. * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] + 4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] ,
	  -4.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] ,
	  Hammer2D::Qsi_6P[2][1] - 4.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] + 4.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_6P[3][0] - 5.5 * Hammer2D::Qsi_6P[3][1] + 9. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 18. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 9. * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] - 4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] - 4.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] ,
	  9. * Hammer2D::Qsi_6P[3][0] - 22.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 22.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] ,
	  -4.5 * Hammer2D::Qsi_6P[3][0] + 18. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] ,
	  Hammer2D::Qsi_6P[3][0] - 4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] ,
	  9. * Hammer2D::Qsi_6P[3][1] - 22.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 22.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] ,
	  27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] ,
	  -4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] ,
	  -4.5 * Hammer2D::Qsi_6P[3][1] + 18. * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] + 4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] ,
	  -4.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] ,
	  Hammer2D::Qsi_6P[3][1] - 4.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] + 4.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_6P[4][0] - 5.5 * Hammer2D::Qsi_6P[4][1] + 9. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 18. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 9. * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] - 4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] - 4.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] ,
	  9. * Hammer2D::Qsi_6P[4][0] - 22.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 22.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] ,
	  -4.5 * Hammer2D::Qsi_6P[4][0] + 18. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] ,
	  Hammer2D::Qsi_6P[4][0] - 4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] ,
	  9. * Hammer2D::Qsi_6P[4][1] - 22.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 22.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] ,
	  27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] ,
	  -4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] ,
	  -4.5 * Hammer2D::Qsi_6P[4][1] + 18. * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] + 4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] ,
	  -4.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] ,
	  Hammer2D::Qsi_6P[4][1] - 4.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] + 4.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_6P[5][0] - 5.5 * Hammer2D::Qsi_6P[5][1] + 9. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 18. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 9. * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] - 4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] - 4.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] ,
	  9. * Hammer2D::Qsi_6P[5][0] - 22.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 22.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] ,
	  -4.5 * Hammer2D::Qsi_6P[5][0] + 18. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] ,
	  Hammer2D::Qsi_6P[5][0] - 4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] ,
	  9. * Hammer2D::Qsi_6P[5][1] - 22.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 22.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] ,
	  27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] ,
	  -4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] ,
	  -4.5 * Hammer2D::Qsi_6P[5][1] + 18. * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] + 4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] ,
	  -4.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] ,
	  Hammer2D::Qsi_6P[5][1] - 4.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] + 4.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] } };

template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<7>::mv_Psi[7][mv_numNodes] = {
	{ 1. - 5.5 * Hammer2D::Qsi_7P[0][0] - 5.5 * Hammer2D::Qsi_7P[0][1] + 9. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 18. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 9. * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] - 4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] - 4.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] ,
	  9. * Hammer2D::Qsi_7P[0][0] - 22.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 22.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] ,
	  -4.5 * Hammer2D::Qsi_7P[0][0] + 18. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] ,
	  Hammer2D::Qsi_7P[0][0] - 4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] ,
	  9. * Hammer2D::Qsi_7P[0][1] - 22.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 22.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] ,
	  27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] ,
	  -4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] ,
	  -4.5 * Hammer2D::Qsi_7P[0][1] + 18. * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] + 4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] ,
	  -4.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] ,
	  Hammer2D::Qsi_7P[0][1] - 4.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] + 4.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_7P[1][0] - 5.5 * Hammer2D::Qsi_7P[1][1] + 9. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 18. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 9. * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] - 4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] - 4.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] ,
	  9. * Hammer2D::Qsi_7P[1][0] - 22.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 22.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] ,
	  -4.5 * Hammer2D::Qsi_7P[1][0] + 18. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] ,
	  Hammer2D::Qsi_7P[1][0] - 4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] ,
	  9. * Hammer2D::Qsi_7P[1][1] - 22.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 22.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] ,
	  27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] ,
	  -4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] ,
	  -4.5 * Hammer2D::Qsi_7P[1][1] + 18. * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] + 4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] ,
	  -4.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] ,
	  Hammer2D::Qsi_7P[1][1] - 4.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] + 4.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_7P[2][0] - 5.5 * Hammer2D::Qsi_7P[2][1] + 9. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 18. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 9. * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] - 4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] - 4.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] ,
	  9. * Hammer2D::Qsi_7P[2][0] - 22.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 22.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] ,
	  -4.5 * Hammer2D::Qsi_7P[2][0] + 18. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] ,
	  Hammer2D::Qsi_7P[2][0] - 4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] ,
	  9. * Hammer2D::Qsi_7P[2][1] - 22.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 22.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] ,
	  27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] ,
	  -4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] ,
	  -4.5 * Hammer2D::Qsi_7P[2][1] + 18. * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] + 4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] ,
	  -4.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] ,
	  Hammer2D::Qsi_7P[2][1] - 4.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] + 4.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_7P[3][0] - 5.5 * Hammer2D::Qsi_7P[3][1] + 9. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 18. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 9. * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] - 4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] - 4.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] ,
	  9. * Hammer2D::Qsi_7P[3][0] - 22.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 22.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] ,
	  -4.5 * Hammer2D::Qsi_7P[3][0] + 18. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] ,
	  Hammer2D::Qsi_7P[3][0] - 4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] ,
	  9. * Hammer2D::Qsi_7P[3][1] - 22.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 22.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] ,
	  27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] ,
	  -4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] ,
	  -4.5 * Hammer2D::Qsi_7P[3][1] + 18. * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] + 4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] ,
	  -4.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] ,
	  Hammer2D::Qsi_7P[3][1] - 4.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] + 4.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_7P[4][0] - 5.5 * Hammer2D::Qsi_7P[4][1] + 9. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 18. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 9. * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] - 4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] - 4.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] ,
	  9. * Hammer2D::Qsi_7P[4][0] - 22.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 22.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] ,
	  -4.5 * Hammer2D::Qsi_7P[4][0] + 18. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] ,
	  Hammer2D::Qsi_7P[4][0] - 4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] ,
	  9. * Hammer2D::Qsi_7P[4][1] - 22.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 22.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] ,
	  27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] ,
	  -4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] ,
	  -4.5 * Hammer2D::Qsi_7P[4][1] + 18. * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] + 4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] ,
	  -4.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] ,
	  Hammer2D::Qsi_7P[4][1] - 4.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] + 4.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_7P[5][0] - 5.5 * Hammer2D::Qsi_7P[5][1] + 9. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 18. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 9. * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] - 4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] - 4.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] ,
	  9. * Hammer2D::Qsi_7P[5][0] - 22.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 22.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] ,
	  -4.5 * Hammer2D::Qsi_7P[5][0] + 18. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] ,
	  Hammer2D::Qsi_7P[5][0] - 4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] ,
	  9. * Hammer2D::Qsi_7P[5][1] - 22.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 22.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] ,
	  27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] ,
	  -4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] ,
	  -4.5 * Hammer2D::Qsi_7P[5][1] + 18. * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] + 4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] ,
	  -4.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] ,
	  Hammer2D::Qsi_7P[5][1] - 4.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] + 4.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_7P[6][0] - 5.5 * Hammer2D::Qsi_7P[6][1] + 9. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 18. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 9. * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] - 4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] - 4.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] ,
	  9. * Hammer2D::Qsi_7P[6][0] - 22.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 22.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] ,
	  -4.5 * Hammer2D::Qsi_7P[6][0] + 18. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] ,
	  Hammer2D::Qsi_7P[6][0] - 4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] ,
	  9. * Hammer2D::Qsi_7P[6][1] - 22.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 22.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] ,
	  27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] ,
	  -4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] ,
	  -4.5 * Hammer2D::Qsi_7P[6][1] + 18. * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] + 4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] ,
	  -4.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] ,
	  Hammer2D::Qsi_7P[6][1] - 4.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] + 4.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] } };

template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<12>::mv_Psi[12][mv_numNodes] = {
	{ 1. - 5.5 * Hammer2D::Qsi_12P[0][0] - 5.5 * Hammer2D::Qsi_12P[0][1] + 9. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 18. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 9. * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] - 4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] - 4.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] ,
	  9. * Hammer2D::Qsi_12P[0][0] - 22.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 22.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] ,
	  -4.5 * Hammer2D::Qsi_12P[0][0] + 18. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] ,
	  Hammer2D::Qsi_12P[0][0] - 4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] ,
	  9. * Hammer2D::Qsi_12P[0][1] - 22.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 22.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] ,
	  27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] ,
	  -4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] ,
	  -4.5 * Hammer2D::Qsi_12P[0][1] + 18. * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] + 4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] ,
	  -4.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] ,
	  Hammer2D::Qsi_12P[0][1] - 4.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] + 4.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[1][0] - 5.5 * Hammer2D::Qsi_12P[1][1] + 9. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 18. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 9. * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] - 4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] - 4.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] ,
	  9. * Hammer2D::Qsi_12P[1][0] - 22.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 22.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] ,
	  -4.5 * Hammer2D::Qsi_12P[1][0] + 18. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] ,
	  Hammer2D::Qsi_12P[1][0] - 4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] ,
	  9. * Hammer2D::Qsi_12P[1][1] - 22.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 22.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] ,
	  27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] ,
	  -4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] ,
	  -4.5 * Hammer2D::Qsi_12P[1][1] + 18. * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] + 4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] ,
	  -4.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] ,
	  Hammer2D::Qsi_12P[1][1] - 4.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] + 4.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[2][0] - 5.5 * Hammer2D::Qsi_12P[2][1] + 9. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 18. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 9. * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] - 4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] - 4.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] ,
	  9. * Hammer2D::Qsi_12P[2][0] - 22.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 22.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] ,
	  -4.5 * Hammer2D::Qsi_12P[2][0] + 18. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] ,
	  Hammer2D::Qsi_12P[2][0] - 4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] ,
	  9. * Hammer2D::Qsi_12P[2][1] - 22.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 22.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] ,
	  27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] ,
	  -4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] ,
	  -4.5 * Hammer2D::Qsi_12P[2][1] + 18. * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] + 4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] ,
	  -4.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] ,
	  Hammer2D::Qsi_12P[2][1] - 4.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] + 4.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[3][0] - 5.5 * Hammer2D::Qsi_12P[3][1] + 9. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 18. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 9. * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] - 4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] - 4.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] ,
	  9. * Hammer2D::Qsi_12P[3][0] - 22.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 22.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] ,
	  -4.5 * Hammer2D::Qsi_12P[3][0] + 18. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] ,
	  Hammer2D::Qsi_12P[3][0] - 4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] ,
	  9. * Hammer2D::Qsi_12P[3][1] - 22.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 22.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] ,
	  27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] ,
	  -4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] ,
	  -4.5 * Hammer2D::Qsi_12P[3][1] + 18. * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] + 4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] ,
	  -4.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] ,
	  Hammer2D::Qsi_12P[3][1] - 4.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] + 4.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[4][0] - 5.5 * Hammer2D::Qsi_12P[4][1] + 9. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 18. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 9. * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] - 4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] - 4.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] ,
	  9. * Hammer2D::Qsi_12P[4][0] - 22.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 22.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] ,
	  -4.5 * Hammer2D::Qsi_12P[4][0] + 18. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] ,
	  Hammer2D::Qsi_12P[4][0] - 4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] ,
	  9. * Hammer2D::Qsi_12P[4][1] - 22.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 22.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] ,
	  27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] ,
	  -4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] ,
	  -4.5 * Hammer2D::Qsi_12P[4][1] + 18. * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] + 4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] ,
	  -4.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] ,
	  Hammer2D::Qsi_12P[4][1] - 4.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] + 4.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[5][0] - 5.5 * Hammer2D::Qsi_12P[5][1] + 9. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 18. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 9. * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] - 4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] - 4.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] ,
	  9. * Hammer2D::Qsi_12P[5][0] - 22.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 22.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] ,
	  -4.5 * Hammer2D::Qsi_12P[5][0] + 18. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] ,
	  Hammer2D::Qsi_12P[5][0] - 4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] ,
	  9. * Hammer2D::Qsi_12P[5][1] - 22.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 22.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] ,
	  27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] ,
	  -4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] ,
	  -4.5 * Hammer2D::Qsi_12P[5][1] + 18. * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] + 4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] ,
	  -4.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] ,
	  Hammer2D::Qsi_12P[5][1] - 4.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] + 4.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[6][0] - 5.5 * Hammer2D::Qsi_12P[6][1] + 9. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 18. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 9. * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] - 4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] - 4.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] ,
	  9. * Hammer2D::Qsi_12P[6][0] - 22.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 22.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] ,
	  -4.5 * Hammer2D::Qsi_12P[6][0] + 18. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] ,
	  Hammer2D::Qsi_12P[6][0] - 4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] ,
	  9. * Hammer2D::Qsi_12P[6][1] - 22.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 22.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] ,
	  27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] ,
	  -4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] ,
	  -4.5 * Hammer2D::Qsi_12P[6][1] + 18. * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] + 4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] ,
	  -4.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] ,
	  Hammer2D::Qsi_12P[6][1] - 4.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] + 4.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[7][0] - 5.5 * Hammer2D::Qsi_12P[7][1] + 9. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 18. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 9. * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] - 4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] - 4.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] ,
	  9. * Hammer2D::Qsi_12P[7][0] - 22.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 22.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] ,
	  -4.5 * Hammer2D::Qsi_12P[7][0] + 18. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] ,
	  Hammer2D::Qsi_12P[7][0] - 4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] ,
	  9. * Hammer2D::Qsi_12P[7][1] - 22.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 22.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] ,
	  27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] ,
	  -4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] ,
	  -4.5 * Hammer2D::Qsi_12P[7][1] + 18. * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] + 4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] ,
	  -4.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] ,
	  Hammer2D::Qsi_12P[7][1] - 4.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] + 4.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[8][0] - 5.5 * Hammer2D::Qsi_12P[8][1] + 9. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 18. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 9. * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] - 4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] - 4.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] ,
	  9. * Hammer2D::Qsi_12P[8][0] - 22.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 22.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] ,
	  -4.5 * Hammer2D::Qsi_12P[8][0] + 18. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] ,
	  Hammer2D::Qsi_12P[8][0] - 4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] ,
	  9. * Hammer2D::Qsi_12P[8][1] - 22.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 22.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] ,
	  27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] ,
	  -4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] ,
	  -4.5 * Hammer2D::Qsi_12P[8][1] + 18. * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] + 4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] ,
	  -4.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] ,
	  Hammer2D::Qsi_12P[8][1] - 4.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] + 4.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[9][0] - 5.5 * Hammer2D::Qsi_12P[9][1] + 9. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 18. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 9. * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] - 4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] - 4.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] ,
	  9. * Hammer2D::Qsi_12P[9][0] - 22.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 22.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] ,
	  -4.5 * Hammer2D::Qsi_12P[9][0] + 18. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] ,
	  Hammer2D::Qsi_12P[9][0] - 4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] ,
	  9. * Hammer2D::Qsi_12P[9][1] - 22.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 22.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] ,
	  27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] ,
	  -4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] ,
	  -4.5 * Hammer2D::Qsi_12P[9][1] + 18. * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] + 4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] ,
	  -4.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] ,
	  Hammer2D::Qsi_12P[9][1] - 4.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] + 4.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[10][0] - 5.5 * Hammer2D::Qsi_12P[10][1] + 9. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 18. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 9. * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] - 4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] - 4.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] ,
	  9. * Hammer2D::Qsi_12P[10][0] - 22.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 22.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] ,
	  -4.5 * Hammer2D::Qsi_12P[10][0] + 18. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] ,
	  Hammer2D::Qsi_12P[10][0] - 4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] ,
	  9. * Hammer2D::Qsi_12P[10][1] - 22.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 22.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] ,
	  27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] ,
	  -4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] ,
	  -4.5 * Hammer2D::Qsi_12P[10][1] + 18. * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] + 4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] ,
	  -4.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] ,
	  Hammer2D::Qsi_12P[10][1] - 4.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] + 4.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_12P[11][0] - 5.5 * Hammer2D::Qsi_12P[11][1] + 9. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 18. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 9. * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] - 4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] - 4.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] ,
	  9. * Hammer2D::Qsi_12P[11][0] - 22.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 22.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] ,
	  -4.5 * Hammer2D::Qsi_12P[11][0] + 18. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] ,
	  Hammer2D::Qsi_12P[11][0] - 4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] ,
	  9. * Hammer2D::Qsi_12P[11][1] - 22.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 22.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] ,
	  27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] ,
	  -4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] ,
	  -4.5 * Hammer2D::Qsi_12P[11][1] + 18. * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] + 4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] ,
	  -4.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] ,
	  Hammer2D::Qsi_12P[11][1] - 4.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] + 4.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] } };

template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<13>::mv_Psi[13][mv_numNodes] = {
	{ 1. - 5.5 * Hammer2D::Qsi_13P[0][0] - 5.5 * Hammer2D::Qsi_13P[0][1] + 9. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 18. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 9. * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] - 4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] - 4.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] ,
	  9. * Hammer2D::Qsi_13P[0][0] - 22.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 22.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] ,
	  -4.5 * Hammer2D::Qsi_13P[0][0] + 18. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] ,
	  Hammer2D::Qsi_13P[0][0] - 4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] ,
	  9. * Hammer2D::Qsi_13P[0][1] - 22.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 22.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] ,
	  27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] ,
	  -4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] ,
	  -4.5 * Hammer2D::Qsi_13P[0][1] + 18. * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] + 4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] ,
	  -4.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] ,
	  Hammer2D::Qsi_13P[0][1] - 4.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] + 4.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[1][0] - 5.5 * Hammer2D::Qsi_13P[1][1] + 9. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 18. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 9. * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] - 4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] - 4.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] ,
	  9. * Hammer2D::Qsi_13P[1][0] - 22.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 22.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] ,
	  -4.5 * Hammer2D::Qsi_13P[1][0] + 18. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] ,
	  Hammer2D::Qsi_13P[1][0] - 4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] ,
	  9. * Hammer2D::Qsi_13P[1][1] - 22.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 22.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] ,
	  27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] ,
	  -4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] ,
	  -4.5 * Hammer2D::Qsi_13P[1][1] + 18. * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] + 4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] ,
	  -4.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] ,
	  Hammer2D::Qsi_13P[1][1] - 4.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] + 4.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[2][0] - 5.5 * Hammer2D::Qsi_13P[2][1] + 9. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 18. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 9. * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] - 4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] - 4.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] ,
	  9. * Hammer2D::Qsi_13P[2][0] - 22.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 22.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] ,
	  -4.5 * Hammer2D::Qsi_13P[2][0] + 18. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] ,
	  Hammer2D::Qsi_13P[2][0] - 4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] ,
	  9. * Hammer2D::Qsi_13P[2][1] - 22.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 22.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] ,
	  27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] ,
	  -4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] ,
	  -4.5 * Hammer2D::Qsi_13P[2][1] + 18. * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] + 4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] ,
	  -4.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] ,
	  Hammer2D::Qsi_13P[2][1] - 4.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] + 4.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[3][0] - 5.5 * Hammer2D::Qsi_13P[3][1] + 9. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 18. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 9. * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] - 4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] - 4.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] ,
	  9. * Hammer2D::Qsi_13P[3][0] - 22.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 22.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] ,
	  -4.5 * Hammer2D::Qsi_13P[3][0] + 18. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] ,
	  Hammer2D::Qsi_13P[3][0] - 4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] ,
	  9. * Hammer2D::Qsi_13P[3][1] - 22.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 22.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] ,
	  27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] ,
	  -4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] ,
	  -4.5 * Hammer2D::Qsi_13P[3][1] + 18. * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] + 4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] ,
	  -4.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] ,
	  Hammer2D::Qsi_13P[3][1] - 4.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] + 4.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[4][0] - 5.5 * Hammer2D::Qsi_13P[4][1] + 9. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 18. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 9. * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] - 4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] - 4.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] ,
	  9. * Hammer2D::Qsi_13P[4][0] - 22.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 22.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] ,
	  -4.5 * Hammer2D::Qsi_13P[4][0] + 18. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] ,
	  Hammer2D::Qsi_13P[4][0] - 4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] ,
	  9. * Hammer2D::Qsi_13P[4][1] - 22.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 22.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] ,
	  27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] ,
	  -4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] ,
	  -4.5 * Hammer2D::Qsi_13P[4][1] + 18. * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] + 4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] ,
	  -4.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] ,
	  Hammer2D::Qsi_13P[4][1] - 4.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] + 4.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[5][0] - 5.5 * Hammer2D::Qsi_13P[5][1] + 9. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 18. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 9. * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] - 4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] - 4.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] ,
	  9. * Hammer2D::Qsi_13P[5][0] - 22.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 22.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] ,
	  -4.5 * Hammer2D::Qsi_13P[5][0] + 18. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] ,
	  Hammer2D::Qsi_13P[5][0] - 4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] ,
	  9. * Hammer2D::Qsi_13P[5][1] - 22.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 22.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] ,
	  27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] ,
	  -4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] ,
	  -4.5 * Hammer2D::Qsi_13P[5][1] + 18. * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] + 4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] ,
	  -4.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] ,
	  Hammer2D::Qsi_13P[5][1] - 4.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] + 4.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[6][0] - 5.5 * Hammer2D::Qsi_13P[6][1] + 9. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 18. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 9. * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] - 4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] - 4.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] ,
	  9. * Hammer2D::Qsi_13P[6][0] - 22.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 22.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] ,
	  -4.5 * Hammer2D::Qsi_13P[6][0] + 18. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] ,
	  Hammer2D::Qsi_13P[6][0] - 4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] ,
	  9. * Hammer2D::Qsi_13P[6][1] - 22.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 22.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] ,
	  27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] ,
	  -4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] ,
	  -4.5 * Hammer2D::Qsi_13P[6][1] + 18. * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] + 4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] ,
	  -4.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] ,
	  Hammer2D::Qsi_13P[6][1] - 4.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] + 4.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[7][0] - 5.5 * Hammer2D::Qsi_13P[7][1] + 9. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 18. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 9. * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] - 4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] - 4.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] ,
	  9. * Hammer2D::Qsi_13P[7][0] - 22.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 22.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] ,
	  -4.5 * Hammer2D::Qsi_13P[7][0] + 18. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] ,
	  Hammer2D::Qsi_13P[7][0] - 4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] ,
	  9. * Hammer2D::Qsi_13P[7][1] - 22.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 22.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] ,
	  27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] ,
	  -4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] ,
	  -4.5 * Hammer2D::Qsi_13P[7][1] + 18. * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] + 4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] ,
	  -4.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] ,
	  Hammer2D::Qsi_13P[7][1] - 4.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] + 4.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[8][0] - 5.5 * Hammer2D::Qsi_13P[8][1] + 9. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 18. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 9. * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] - 4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] - 4.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] ,
	  9. * Hammer2D::Qsi_13P[8][0] - 22.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 22.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] ,
	  -4.5 * Hammer2D::Qsi_13P[8][0] + 18. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] ,
	  Hammer2D::Qsi_13P[8][0] - 4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] ,
	  9. * Hammer2D::Qsi_13P[8][1] - 22.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 22.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] ,
	  27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] ,
	  -4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] ,
	  -4.5 * Hammer2D::Qsi_13P[8][1] + 18. * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] + 4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] ,
	  -4.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] ,
	  Hammer2D::Qsi_13P[8][1] - 4.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] + 4.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[9][0] - 5.5 * Hammer2D::Qsi_13P[9][1] + 9. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 18. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 9. * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] - 4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] - 4.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] ,
	  9. * Hammer2D::Qsi_13P[9][0] - 22.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 22.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] ,
	  -4.5 * Hammer2D::Qsi_13P[9][0] + 18. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] ,
	  Hammer2D::Qsi_13P[9][0] - 4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] ,
	  9. * Hammer2D::Qsi_13P[9][1] - 22.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 22.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] ,
	  27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] ,
	  -4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] ,
	  -4.5 * Hammer2D::Qsi_13P[9][1] + 18. * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] + 4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] ,
	  -4.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] ,
	  Hammer2D::Qsi_13P[9][1] - 4.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] + 4.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[10][0] - 5.5 * Hammer2D::Qsi_13P[10][1] + 9. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 18. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 9. * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] - 4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] - 4.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] ,
	  9. * Hammer2D::Qsi_13P[10][0] - 22.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 22.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] ,
	  -4.5 * Hammer2D::Qsi_13P[10][0] + 18. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] ,
	  Hammer2D::Qsi_13P[10][0] - 4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] ,
	  9. * Hammer2D::Qsi_13P[10][1] - 22.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 22.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] ,
	  27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] ,
	  -4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] ,
	  -4.5 * Hammer2D::Qsi_13P[10][1] + 18. * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] + 4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] ,
	  -4.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] ,
	  Hammer2D::Qsi_13P[10][1] - 4.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] + 4.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[11][0] - 5.5 * Hammer2D::Qsi_13P[11][1] + 9. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 18. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 9. * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] - 4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] - 4.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] ,
	  9. * Hammer2D::Qsi_13P[11][0] - 22.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 22.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] ,
	  -4.5 * Hammer2D::Qsi_13P[11][0] + 18. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] ,
	  Hammer2D::Qsi_13P[11][0] - 4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] ,
	  9. * Hammer2D::Qsi_13P[11][1] - 22.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 22.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] ,
	  27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] ,
	  -4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] ,
	  -4.5 * Hammer2D::Qsi_13P[11][1] + 18. * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] + 4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] ,
	  -4.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] ,
	  Hammer2D::Qsi_13P[11][1] - 4.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] + 4.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] } ,

	{ 1. - 5.5 * Hammer2D::Qsi_13P[12][0] - 5.5 * Hammer2D::Qsi_13P[12][1] + 9. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 18. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 9. * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] - 4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] - 4.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] ,
	  9. * Hammer2D::Qsi_13P[12][0] - 22.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 22.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] ,
	  -4.5 * Hammer2D::Qsi_13P[12][0] + 18. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] ,
	  Hammer2D::Qsi_13P[12][0] - 4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] ,
	  9. * Hammer2D::Qsi_13P[12][1] - 22.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 22.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] ,
	  27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] ,
	  -4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] ,
	  -4.5 * Hammer2D::Qsi_13P[12][1] + 18. * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] + 4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] ,
	  -4.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] ,
	  Hammer2D::Qsi_13P[12][1] - 4.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] + 4.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<6>::mv_DPsi[6][mv_numNodes][mv_ElDim] = {
	{ { -5.5 + 18. * Hammer2D::Qsi_6P[0][0] + 18. * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] , -5.5 + 18. * Hammer2D::Qsi_6P[0][0] + 18. * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_6P[0][0] - 22.5 * Hammer2D::Qsi_6P[0][1] + 40.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 54. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] , -22.5 * Hammer2D::Qsi_6P[0][0] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_6P[0][0] + 4.5 * Hammer2D::Qsi_6P[0][1] - 40.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] , 4.5 * Hammer2D::Qsi_6P[0][0] - 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_6P[0][0] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_6P[0][1] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 27. * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] , 9. - 22.5 * Hammer2D::Qsi_6P[0][0] - 45. * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] + 54. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] + 40.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] } ,
	  { 27. * Hammer2D::Qsi_6P[0][1] - 54. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 27. * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] , 27. * Hammer2D::Qsi_6P[0][0] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] - 54. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[0][1] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] , -4.5 * Hammer2D::Qsi_6P[0][0] + 13.5 * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][0] } ,
	  {  4.5 * Hammer2D::Qsi_6P[0][1] - 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] , -4.5 + 36. * Hammer2D::Qsi_6P[0][1] + 4.5 * Hammer2D::Qsi_6P[0][0] - 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] - 40.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] , -4.5 * Hammer2D::Qsi_6P[0][0] + 27. * Hammer2D::Qsi_6P[0][0] * Hammer2D::Qsi_6P[0][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_6P[0][1] + 13.5 * Hammer2D::Qsi_6P[0][1] * Hammer2D::Qsi_6P[0][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_6P[1][0] + 18. * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] , -5.5 + 18. * Hammer2D::Qsi_6P[1][0] + 18. * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_6P[1][0] - 22.5 * Hammer2D::Qsi_6P[1][1] + 40.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 54. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] , -22.5 * Hammer2D::Qsi_6P[1][0] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_6P[1][0] + 4.5 * Hammer2D::Qsi_6P[1][1] - 40.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] , 4.5 * Hammer2D::Qsi_6P[1][0] - 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_6P[1][0] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_6P[1][1] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 27. * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] , 9. - 22.5 * Hammer2D::Qsi_6P[1][0] - 45. * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] + 54. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] + 40.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] } ,
	  { 27. * Hammer2D::Qsi_6P[1][1] - 54. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 27. * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] , 27. * Hammer2D::Qsi_6P[1][0] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] - 54. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[1][1] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] , -4.5 * Hammer2D::Qsi_6P[1][0] + 13.5 * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][0] } ,
	  {  4.5 * Hammer2D::Qsi_6P[1][1] - 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] , -4.5 + 36. * Hammer2D::Qsi_6P[1][1] + 4.5 * Hammer2D::Qsi_6P[1][0] - 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] - 40.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] , -4.5 * Hammer2D::Qsi_6P[1][0] + 27. * Hammer2D::Qsi_6P[1][0] * Hammer2D::Qsi_6P[1][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_6P[1][1] + 13.5 * Hammer2D::Qsi_6P[1][1] * Hammer2D::Qsi_6P[1][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_6P[2][0] + 18. * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] , -5.5 + 18. * Hammer2D::Qsi_6P[2][0] + 18. * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_6P[2][0] - 22.5 * Hammer2D::Qsi_6P[2][1] + 40.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 54. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] , -22.5 * Hammer2D::Qsi_6P[2][0] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_6P[2][0] + 4.5 * Hammer2D::Qsi_6P[2][1] - 40.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] , 4.5 * Hammer2D::Qsi_6P[2][0] - 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_6P[2][0] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_6P[2][1] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 27. * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] , 9. - 22.5 * Hammer2D::Qsi_6P[2][0] - 45. * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] + 54. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] + 40.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] } ,
	  { 27. * Hammer2D::Qsi_6P[2][1] - 54. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 27. * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] , 27. * Hammer2D::Qsi_6P[2][0] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] - 54. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[2][1] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] , -4.5 * Hammer2D::Qsi_6P[2][0] + 13.5 * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][0] } ,
	  {  4.5 * Hammer2D::Qsi_6P[2][1] - 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] , -4.5 + 36. * Hammer2D::Qsi_6P[2][1] + 4.5 * Hammer2D::Qsi_6P[2][0] - 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] - 40.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] , -4.5 * Hammer2D::Qsi_6P[2][0] + 27. * Hammer2D::Qsi_6P[2][0] * Hammer2D::Qsi_6P[2][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_6P[2][1] + 13.5 * Hammer2D::Qsi_6P[2][1] * Hammer2D::Qsi_6P[2][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_6P[3][0] + 18. * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] , -5.5 + 18. * Hammer2D::Qsi_6P[3][0] + 18. * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_6P[3][0] - 22.5 * Hammer2D::Qsi_6P[3][1] + 40.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 54. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] , -22.5 * Hammer2D::Qsi_6P[3][0] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_6P[3][0] + 4.5 * Hammer2D::Qsi_6P[3][1] - 40.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] , 4.5 * Hammer2D::Qsi_6P[3][0] - 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_6P[3][0] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_6P[3][1] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 27. * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] , 9. - 22.5 * Hammer2D::Qsi_6P[3][0] - 45. * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] + 54. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] + 40.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] } ,
	  { 27. * Hammer2D::Qsi_6P[3][1] - 54. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 27. * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] , 27. * Hammer2D::Qsi_6P[3][0] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] - 54. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[3][1] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] , -4.5 * Hammer2D::Qsi_6P[3][0] + 13.5 * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][0] } ,
	  {  4.5 * Hammer2D::Qsi_6P[3][1] - 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] , -4.5 + 36. * Hammer2D::Qsi_6P[3][1] + 4.5 * Hammer2D::Qsi_6P[3][0] - 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] - 40.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] , -4.5 * Hammer2D::Qsi_6P[3][0] + 27. * Hammer2D::Qsi_6P[3][0] * Hammer2D::Qsi_6P[3][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_6P[3][1] + 13.5 * Hammer2D::Qsi_6P[3][1] * Hammer2D::Qsi_6P[3][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_6P[4][0] + 18. * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] , -5.5 + 18. * Hammer2D::Qsi_6P[4][0] + 18. * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_6P[4][0] - 22.5 * Hammer2D::Qsi_6P[4][1] + 40.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 54. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] , -22.5 * Hammer2D::Qsi_6P[4][0] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_6P[4][0] + 4.5 * Hammer2D::Qsi_6P[4][1] - 40.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] , 4.5 * Hammer2D::Qsi_6P[4][0] - 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_6P[4][0] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_6P[4][1] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 27. * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] , 9. - 22.5 * Hammer2D::Qsi_6P[4][0] - 45. * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] + 54. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] + 40.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] } ,
	  { 27. * Hammer2D::Qsi_6P[4][1] - 54. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 27. * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] , 27. * Hammer2D::Qsi_6P[4][0] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] - 54. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[4][1] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] , -4.5 * Hammer2D::Qsi_6P[4][0] + 13.5 * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][0] } ,
	  {  4.5 * Hammer2D::Qsi_6P[4][1] - 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] , -4.5 + 36. * Hammer2D::Qsi_6P[4][1] + 4.5 * Hammer2D::Qsi_6P[4][0] - 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] - 40.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] , -4.5 * Hammer2D::Qsi_6P[4][0] + 27. * Hammer2D::Qsi_6P[4][0] * Hammer2D::Qsi_6P[4][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_6P[4][1] + 13.5 * Hammer2D::Qsi_6P[4][1] * Hammer2D::Qsi_6P[4][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_6P[5][0] + 18. * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] , -5.5 + 18. * Hammer2D::Qsi_6P[5][0] + 18. * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_6P[5][0] - 22.5 * Hammer2D::Qsi_6P[5][1] + 40.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 54. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] , -22.5 * Hammer2D::Qsi_6P[5][0] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_6P[5][0] + 4.5 * Hammer2D::Qsi_6P[5][1] - 40.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] , 4.5 * Hammer2D::Qsi_6P[5][0] - 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_6P[5][0] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_6P[5][1] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 27. * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] , 9. - 22.5 * Hammer2D::Qsi_6P[5][0] - 45. * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] + 54. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] + 40.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] } ,
	  { 27. * Hammer2D::Qsi_6P[5][1] - 54. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 27. * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] , 27. * Hammer2D::Qsi_6P[5][0] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] - 54. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[5][1] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] , -4.5 * Hammer2D::Qsi_6P[5][0] + 13.5 * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][0] } ,
	  {  4.5 * Hammer2D::Qsi_6P[5][1] - 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] , -4.5 + 36. * Hammer2D::Qsi_6P[5][1] + 4.5 * Hammer2D::Qsi_6P[5][0] - 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] - 40.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] , -4.5 * Hammer2D::Qsi_6P[5][0] + 27. * Hammer2D::Qsi_6P[5][0] * Hammer2D::Qsi_6P[5][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_6P[5][1] + 13.5 * Hammer2D::Qsi_6P[5][1] * Hammer2D::Qsi_6P[5][1] } } };

template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<7>::mv_DPsi[7][mv_numNodes][mv_ElDim] = {
	{ { -5.5 + 18. * Hammer2D::Qsi_7P[0][0] + 18. * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] , -5.5 + 18. * Hammer2D::Qsi_7P[0][0] + 18. * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[0][0] - 22.5 * Hammer2D::Qsi_7P[0][1] + 40.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 54. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] , -22.5 * Hammer2D::Qsi_7P[0][0] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[0][0] + 4.5 * Hammer2D::Qsi_7P[0][1] - 40.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] , 4.5 * Hammer2D::Qsi_7P[0][0] - 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[0][0] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[0][1] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 27. * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] , 9. - 22.5 * Hammer2D::Qsi_7P[0][0] - 45. * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] + 54. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] + 40.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] } ,
	  { 27. * Hammer2D::Qsi_7P[0][1] - 54. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 27. * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] , 27. * Hammer2D::Qsi_7P[0][0] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] - 54. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[0][1] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] , -4.5 * Hammer2D::Qsi_7P[0][0] + 13.5 * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[0][1] - 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] , -4.5 + 36. * Hammer2D::Qsi_7P[0][1] + 4.5 * Hammer2D::Qsi_7P[0][0] - 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] - 40.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] , -4.5 * Hammer2D::Qsi_7P[0][0] + 27. * Hammer2D::Qsi_7P[0][0] * Hammer2D::Qsi_7P[0][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[0][1] + 13.5 * Hammer2D::Qsi_7P[0][1] * Hammer2D::Qsi_7P[0][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_7P[1][0] + 18. * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] , -5.5 + 18. * Hammer2D::Qsi_7P[1][0] + 18. * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[1][0] - 22.5 * Hammer2D::Qsi_7P[1][1] + 40.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 54. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] , -22.5 * Hammer2D::Qsi_7P[1][0] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[1][0] + 4.5 * Hammer2D::Qsi_7P[1][1] - 40.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] , 4.5 * Hammer2D::Qsi_7P[1][0] - 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[1][0] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[1][1] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 27. * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] , 9. - 22.5 * Hammer2D::Qsi_7P[1][0] - 45. * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] + 54. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] + 40.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] } ,
	  { 27. * Hammer2D::Qsi_7P[1][1] - 54. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 27. * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] , 27. * Hammer2D::Qsi_7P[1][0] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] - 54. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[1][1] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] , -4.5 * Hammer2D::Qsi_7P[1][0] + 13.5 * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[1][1] - 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] , -4.5 + 36. * Hammer2D::Qsi_7P[1][1] + 4.5 * Hammer2D::Qsi_7P[1][0] - 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] - 40.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] , -4.5 * Hammer2D::Qsi_7P[1][0] + 27. * Hammer2D::Qsi_7P[1][0] * Hammer2D::Qsi_7P[1][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[1][1] + 13.5 * Hammer2D::Qsi_7P[1][1] * Hammer2D::Qsi_7P[1][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_7P[2][0] + 18. * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] , -5.5 + 18. * Hammer2D::Qsi_7P[2][0] + 18. * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[2][0] - 22.5 * Hammer2D::Qsi_7P[2][1] + 40.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 54. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] , -22.5 * Hammer2D::Qsi_7P[2][0] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[2][0] + 4.5 * Hammer2D::Qsi_7P[2][1] - 40.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] , 4.5 * Hammer2D::Qsi_7P[2][0] - 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[2][0] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[2][1] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 27. * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] , 9. - 22.5 * Hammer2D::Qsi_7P[2][0] - 45. * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] + 54. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] + 40.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] } ,
	  { 27. * Hammer2D::Qsi_7P[2][1] - 54. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 27. * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] , 27. * Hammer2D::Qsi_7P[2][0] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] - 54. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[2][1] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] , -4.5 * Hammer2D::Qsi_7P[2][0] + 13.5 * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[2][1] - 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] , -4.5 + 36. * Hammer2D::Qsi_7P[2][1] + 4.5 * Hammer2D::Qsi_7P[2][0] - 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] - 40.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] , -4.5 * Hammer2D::Qsi_7P[2][0] + 27. * Hammer2D::Qsi_7P[2][0] * Hammer2D::Qsi_7P[2][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[2][1] + 13.5 * Hammer2D::Qsi_7P[2][1] * Hammer2D::Qsi_7P[2][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_7P[3][0] + 18. * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] , -5.5 + 18. * Hammer2D::Qsi_7P[3][0] + 18. * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[3][0] - 22.5 * Hammer2D::Qsi_7P[3][1] + 40.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 54. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] , -22.5 * Hammer2D::Qsi_7P[3][0] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[3][0] + 4.5 * Hammer2D::Qsi_7P[3][1] - 40.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] , 4.5 * Hammer2D::Qsi_7P[3][0] - 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[3][0] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[3][1] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 27. * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] , 9. - 22.5 * Hammer2D::Qsi_7P[3][0] - 45. * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] + 54. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] + 40.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] } ,
	  { 27. * Hammer2D::Qsi_7P[3][1] - 54. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 27. * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] , 27. * Hammer2D::Qsi_7P[3][0] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] - 54. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[3][1] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] , -4.5 * Hammer2D::Qsi_7P[3][0] + 13.5 * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[3][1] - 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] , -4.5 + 36. * Hammer2D::Qsi_7P[3][1] + 4.5 * Hammer2D::Qsi_7P[3][0] - 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] - 40.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] , -4.5 * Hammer2D::Qsi_7P[3][0] + 27. * Hammer2D::Qsi_7P[3][0] * Hammer2D::Qsi_7P[3][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[3][1] + 13.5 * Hammer2D::Qsi_7P[3][1] * Hammer2D::Qsi_7P[3][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_7P[4][0] + 18. * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] , -5.5 + 18. * Hammer2D::Qsi_7P[4][0] + 18. * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[4][0] - 22.5 * Hammer2D::Qsi_7P[4][1] + 40.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 54. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] , -22.5 * Hammer2D::Qsi_7P[4][0] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[4][0] + 4.5 * Hammer2D::Qsi_7P[4][1] - 40.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] , 4.5 * Hammer2D::Qsi_7P[4][0] - 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[4][0] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[4][1] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 27. * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] , 9. - 22.5 * Hammer2D::Qsi_7P[4][0] - 45. * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] + 54. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] + 40.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] } ,
	  { 27. * Hammer2D::Qsi_7P[4][1] - 54. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 27. * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] , 27. * Hammer2D::Qsi_7P[4][0] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] - 54. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[4][1] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] , -4.5 * Hammer2D::Qsi_7P[4][0] + 13.5 * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[4][1] - 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] , -4.5 + 36. * Hammer2D::Qsi_7P[4][1] + 4.5 * Hammer2D::Qsi_7P[4][0] - 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] - 40.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] , -4.5 * Hammer2D::Qsi_7P[4][0] + 27. * Hammer2D::Qsi_7P[4][0] * Hammer2D::Qsi_7P[4][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[4][1] + 13.5 * Hammer2D::Qsi_7P[4][1] * Hammer2D::Qsi_7P[4][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_7P[5][0] + 18. * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] , -5.5 + 18. * Hammer2D::Qsi_7P[5][0] + 18. * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[5][0] - 22.5 * Hammer2D::Qsi_7P[5][1] + 40.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 54. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] , -22.5 * Hammer2D::Qsi_7P[5][0] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[5][0] + 4.5 * Hammer2D::Qsi_7P[5][1] - 40.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] , 4.5 * Hammer2D::Qsi_7P[5][0] - 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[5][0] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[5][1] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 27. * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] , 9. - 22.5 * Hammer2D::Qsi_7P[5][0] - 45. * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] + 54. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] + 40.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] } ,
	  { 27. * Hammer2D::Qsi_7P[5][1] - 54. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 27. * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] , 27. * Hammer2D::Qsi_7P[5][0] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] - 54. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[5][1] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] , -4.5 * Hammer2D::Qsi_7P[5][0] + 13.5 * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[5][1] - 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] , -4.5 + 36. * Hammer2D::Qsi_7P[5][1] + 4.5 * Hammer2D::Qsi_7P[5][0] - 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] - 40.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] , -4.5 * Hammer2D::Qsi_7P[5][0] + 27. * Hammer2D::Qsi_7P[5][0] * Hammer2D::Qsi_7P[5][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[5][1] + 13.5 * Hammer2D::Qsi_7P[5][1] * Hammer2D::Qsi_7P[5][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_7P[6][0] + 18. * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] , -5.5 + 18. * Hammer2D::Qsi_7P[6][0] + 18. * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_7P[6][0] - 22.5 * Hammer2D::Qsi_7P[6][1] + 40.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 54. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] , -22.5 * Hammer2D::Qsi_7P[6][0] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_7P[6][0] + 4.5 * Hammer2D::Qsi_7P[6][1] - 40.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] , 4.5 * Hammer2D::Qsi_7P[6][0] - 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_7P[6][0] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_7P[6][1] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 27. * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] , 9. - 22.5 * Hammer2D::Qsi_7P[6][0] - 45. * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] + 54. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] + 40.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] } ,
	  { 27. * Hammer2D::Qsi_7P[6][1] - 54. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 27. * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] , 27. * Hammer2D::Qsi_7P[6][0] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] - 54. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[6][1] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] , -4.5 * Hammer2D::Qsi_7P[6][0] + 13.5 * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][0] } ,
	  {  4.5 * Hammer2D::Qsi_7P[6][1] - 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] , -4.5 + 36. * Hammer2D::Qsi_7P[6][1] + 4.5 * Hammer2D::Qsi_7P[6][0] - 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] - 40.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] } ,
	  { -4.5 * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] , -4.5 * Hammer2D::Qsi_7P[6][0] + 27. * Hammer2D::Qsi_7P[6][0] * Hammer2D::Qsi_7P[6][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_7P[6][1] + 13.5 * Hammer2D::Qsi_7P[6][1] * Hammer2D::Qsi_7P[6][1] } } };

template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<12>::mv_DPsi[12][mv_numNodes][mv_ElDim] = {
	{ { -5.5 + 18. * Hammer2D::Qsi_12P[0][0] + 18. * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] , -5.5 + 18. * Hammer2D::Qsi_12P[0][0] + 18. * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[0][0] - 22.5 * Hammer2D::Qsi_12P[0][1] + 40.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 54. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] , -22.5 * Hammer2D::Qsi_12P[0][0] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[0][0] + 4.5 * Hammer2D::Qsi_12P[0][1] - 40.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] , 4.5 * Hammer2D::Qsi_12P[0][0] - 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[0][0] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[0][1] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 27. * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] , 9. - 22.5 * Hammer2D::Qsi_12P[0][0] - 45. * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] + 54. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] + 40.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] } ,
	  { 27. * Hammer2D::Qsi_12P[0][1] - 54. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 27. * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] , 27. * Hammer2D::Qsi_12P[0][0] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] - 54. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[0][1] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] , -4.5 * Hammer2D::Qsi_12P[0][0] + 13.5 * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[0][1] - 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] , -4.5 + 36. * Hammer2D::Qsi_12P[0][1] + 4.5 * Hammer2D::Qsi_12P[0][0] - 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] - 40.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] , -4.5 * Hammer2D::Qsi_12P[0][0] + 27. * Hammer2D::Qsi_12P[0][0] * Hammer2D::Qsi_12P[0][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[0][1] + 13.5 * Hammer2D::Qsi_12P[0][1] * Hammer2D::Qsi_12P[0][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[1][0] + 18. * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] , -5.5 + 18. * Hammer2D::Qsi_12P[1][0] + 18. * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[1][0] - 22.5 * Hammer2D::Qsi_12P[1][1] + 40.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 54. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] , -22.5 * Hammer2D::Qsi_12P[1][0] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[1][0] + 4.5 * Hammer2D::Qsi_12P[1][1] - 40.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] , 4.5 * Hammer2D::Qsi_12P[1][0] - 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[1][0] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[1][1] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 27. * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] , 9. - 22.5 * Hammer2D::Qsi_12P[1][0] - 45. * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] + 54. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] + 40.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] } ,
	  { 27. * Hammer2D::Qsi_12P[1][1] - 54. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 27. * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] , 27. * Hammer2D::Qsi_12P[1][0] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] - 54. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[1][1] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] , -4.5 * Hammer2D::Qsi_12P[1][0] + 13.5 * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[1][1] - 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] , -4.5 + 36. * Hammer2D::Qsi_12P[1][1] + 4.5 * Hammer2D::Qsi_12P[1][0] - 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] - 40.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] , -4.5 * Hammer2D::Qsi_12P[1][0] + 27. * Hammer2D::Qsi_12P[1][0] * Hammer2D::Qsi_12P[1][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[1][1] + 13.5 * Hammer2D::Qsi_12P[1][1] * Hammer2D::Qsi_12P[1][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[2][0] + 18. * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] , -5.5 + 18. * Hammer2D::Qsi_12P[2][0] + 18. * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[2][0] - 22.5 * Hammer2D::Qsi_12P[2][1] + 40.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 54. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] , -22.5 * Hammer2D::Qsi_12P[2][0] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[2][0] + 4.5 * Hammer2D::Qsi_12P[2][1] - 40.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] , 4.5 * Hammer2D::Qsi_12P[2][0] - 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[2][0] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[2][1] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 27. * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] , 9. - 22.5 * Hammer2D::Qsi_12P[2][0] - 45. * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] + 54. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] + 40.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] } ,
	  { 27. * Hammer2D::Qsi_12P[2][1] - 54. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 27. * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] , 27. * Hammer2D::Qsi_12P[2][0] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] - 54. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[2][1] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] , -4.5 * Hammer2D::Qsi_12P[2][0] + 13.5 * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[2][1] - 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] , -4.5 + 36. * Hammer2D::Qsi_12P[2][1] + 4.5 * Hammer2D::Qsi_12P[2][0] - 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] - 40.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] , -4.5 * Hammer2D::Qsi_12P[2][0] + 27. * Hammer2D::Qsi_12P[2][0] * Hammer2D::Qsi_12P[2][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[2][1] + 13.5 * Hammer2D::Qsi_12P[2][1] * Hammer2D::Qsi_12P[2][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[3][0] + 18. * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] , -5.5 + 18. * Hammer2D::Qsi_12P[3][0] + 18. * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[3][0] - 22.5 * Hammer2D::Qsi_12P[3][1] + 40.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 54. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] , -22.5 * Hammer2D::Qsi_12P[3][0] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[3][0] + 4.5 * Hammer2D::Qsi_12P[3][1] - 40.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] , 4.5 * Hammer2D::Qsi_12P[3][0] - 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[3][0] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[3][1] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 27. * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] , 9. - 22.5 * Hammer2D::Qsi_12P[3][0] - 45. * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] + 54. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] + 40.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] } ,
	  { 27. * Hammer2D::Qsi_12P[3][1] - 54. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 27. * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] , 27. * Hammer2D::Qsi_12P[3][0] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] - 54. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[3][1] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] , -4.5 * Hammer2D::Qsi_12P[3][0] + 13.5 * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[3][1] - 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] , -4.5 + 36. * Hammer2D::Qsi_12P[3][1] + 4.5 * Hammer2D::Qsi_12P[3][0] - 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] - 40.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] , -4.5 * Hammer2D::Qsi_12P[3][0] + 27. * Hammer2D::Qsi_12P[3][0] * Hammer2D::Qsi_12P[3][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[3][1] + 13.5 * Hammer2D::Qsi_12P[3][1] * Hammer2D::Qsi_12P[3][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[4][0] + 18. * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] , -5.5 + 18. * Hammer2D::Qsi_12P[4][0] + 18. * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[4][0] - 22.5 * Hammer2D::Qsi_12P[4][1] + 40.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 54. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] , -22.5 * Hammer2D::Qsi_12P[4][0] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[4][0] + 4.5 * Hammer2D::Qsi_12P[4][1] - 40.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] , 4.5 * Hammer2D::Qsi_12P[4][0] - 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[4][0] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[4][1] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 27. * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] , 9. - 22.5 * Hammer2D::Qsi_12P[4][0] - 45. * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] + 54. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] + 40.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] } ,
	  { 27. * Hammer2D::Qsi_12P[4][1] - 54. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 27. * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] , 27. * Hammer2D::Qsi_12P[4][0] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] - 54. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[4][1] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] , -4.5 * Hammer2D::Qsi_12P[4][0] + 13.5 * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[4][1] - 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] , -4.5 + 36. * Hammer2D::Qsi_12P[4][1] + 4.5 * Hammer2D::Qsi_12P[4][0] - 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] - 40.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] , -4.5 * Hammer2D::Qsi_12P[4][0] + 27. * Hammer2D::Qsi_12P[4][0] * Hammer2D::Qsi_12P[4][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[4][1] + 13.5 * Hammer2D::Qsi_12P[4][1] * Hammer2D::Qsi_12P[4][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[5][0] + 18. * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] , -5.5 + 18. * Hammer2D::Qsi_12P[5][0] + 18. * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[5][0] - 22.5 * Hammer2D::Qsi_12P[5][1] + 40.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 54. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] , -22.5 * Hammer2D::Qsi_12P[5][0] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[5][0] + 4.5 * Hammer2D::Qsi_12P[5][1] - 40.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] , 4.5 * Hammer2D::Qsi_12P[5][0] - 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[5][0] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[5][1] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 27. * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] , 9. - 22.5 * Hammer2D::Qsi_12P[5][0] - 45. * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] + 54. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] + 40.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] } ,
	  { 27. * Hammer2D::Qsi_12P[5][1] - 54. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 27. * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] , 27. * Hammer2D::Qsi_12P[5][0] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] - 54. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[5][1] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] , -4.5 * Hammer2D::Qsi_12P[5][0] + 13.5 * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[5][1] - 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] , -4.5 + 36. * Hammer2D::Qsi_12P[5][1] + 4.5 * Hammer2D::Qsi_12P[5][0] - 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] - 40.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] , -4.5 * Hammer2D::Qsi_12P[5][0] + 27. * Hammer2D::Qsi_12P[5][0] * Hammer2D::Qsi_12P[5][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[5][1] + 13.5 * Hammer2D::Qsi_12P[5][1] * Hammer2D::Qsi_12P[5][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[6][0] + 18. * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] , -5.5 + 18. * Hammer2D::Qsi_12P[6][0] + 18. * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[6][0] - 22.5 * Hammer2D::Qsi_12P[6][1] + 40.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 54. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] , -22.5 * Hammer2D::Qsi_12P[6][0] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[6][0] + 4.5 * Hammer2D::Qsi_12P[6][1] - 40.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] , 4.5 * Hammer2D::Qsi_12P[6][0] - 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[6][0] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[6][1] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 27. * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] , 9. - 22.5 * Hammer2D::Qsi_12P[6][0] - 45. * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] + 54. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] + 40.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] } ,
	  { 27. * Hammer2D::Qsi_12P[6][1] - 54. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 27. * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] , 27. * Hammer2D::Qsi_12P[6][0] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] - 54. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[6][1] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] , -4.5 * Hammer2D::Qsi_12P[6][0] + 13.5 * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[6][1] - 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] , -4.5 + 36. * Hammer2D::Qsi_12P[6][1] + 4.5 * Hammer2D::Qsi_12P[6][0] - 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] - 40.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] , -4.5 * Hammer2D::Qsi_12P[6][0] + 27. * Hammer2D::Qsi_12P[6][0] * Hammer2D::Qsi_12P[6][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[6][1] + 13.5 * Hammer2D::Qsi_12P[6][1] * Hammer2D::Qsi_12P[6][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[7][0] + 18. * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] , -5.5 + 18. * Hammer2D::Qsi_12P[7][0] + 18. * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[7][0] - 22.5 * Hammer2D::Qsi_12P[7][1] + 40.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 54. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] , -22.5 * Hammer2D::Qsi_12P[7][0] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[7][0] + 4.5 * Hammer2D::Qsi_12P[7][1] - 40.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] , 4.5 * Hammer2D::Qsi_12P[7][0] - 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[7][0] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[7][1] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 27. * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] , 9. - 22.5 * Hammer2D::Qsi_12P[7][0] - 45. * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] + 54. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] + 40.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] } ,
	  { 27. * Hammer2D::Qsi_12P[7][1] - 54. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 27. * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] , 27. * Hammer2D::Qsi_12P[7][0] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] - 54. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[7][1] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] , -4.5 * Hammer2D::Qsi_12P[7][0] + 13.5 * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[7][1] - 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] , -4.5 + 36. * Hammer2D::Qsi_12P[7][1] + 4.5 * Hammer2D::Qsi_12P[7][0] - 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] - 40.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] , -4.5 * Hammer2D::Qsi_12P[7][0] + 27. * Hammer2D::Qsi_12P[7][0] * Hammer2D::Qsi_12P[7][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[7][1] + 13.5 * Hammer2D::Qsi_12P[7][1] * Hammer2D::Qsi_12P[7][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[8][0] + 18. * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] , -5.5 + 18. * Hammer2D::Qsi_12P[8][0] + 18. * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[8][0] - 22.5 * Hammer2D::Qsi_12P[8][1] + 40.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 54. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] , -22.5 * Hammer2D::Qsi_12P[8][0] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[8][0] + 4.5 * Hammer2D::Qsi_12P[8][1] - 40.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] , 4.5 * Hammer2D::Qsi_12P[8][0] - 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[8][0] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[8][1] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 27. * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] , 9. - 22.5 * Hammer2D::Qsi_12P[8][0] - 45. * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] + 54. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] + 40.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] } ,
	  { 27. * Hammer2D::Qsi_12P[8][1] - 54. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 27. * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] , 27. * Hammer2D::Qsi_12P[8][0] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] - 54. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[8][1] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] , -4.5 * Hammer2D::Qsi_12P[8][0] + 13.5 * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[8][1] - 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] , -4.5 + 36. * Hammer2D::Qsi_12P[8][1] + 4.5 * Hammer2D::Qsi_12P[8][0] - 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] - 40.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] , -4.5 * Hammer2D::Qsi_12P[8][0] + 27. * Hammer2D::Qsi_12P[8][0] * Hammer2D::Qsi_12P[8][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[8][1] + 13.5 * Hammer2D::Qsi_12P[8][1] * Hammer2D::Qsi_12P[8][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[9][0] + 18. * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] , -5.5 + 18. * Hammer2D::Qsi_12P[9][0] + 18. * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[9][0] - 22.5 * Hammer2D::Qsi_12P[9][1] + 40.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 54. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] , -22.5 * Hammer2D::Qsi_12P[9][0] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[9][0] + 4.5 * Hammer2D::Qsi_12P[9][1] - 40.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] , 4.5 * Hammer2D::Qsi_12P[9][0] - 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[9][0] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[9][1] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 27. * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] , 9. - 22.5 * Hammer2D::Qsi_12P[9][0] - 45. * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] + 54. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] + 40.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] } ,
	  { 27. * Hammer2D::Qsi_12P[9][1] - 54. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 27. * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] , 27. * Hammer2D::Qsi_12P[9][0] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] - 54. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[9][1] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] , -4.5 * Hammer2D::Qsi_12P[9][0] + 13.5 * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[9][1] - 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] , -4.5 + 36. * Hammer2D::Qsi_12P[9][1] + 4.5 * Hammer2D::Qsi_12P[9][0] - 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] - 40.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] , -4.5 * Hammer2D::Qsi_12P[9][0] + 27. * Hammer2D::Qsi_12P[9][0] * Hammer2D::Qsi_12P[9][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[9][1] + 13.5 * Hammer2D::Qsi_12P[9][1] * Hammer2D::Qsi_12P[9][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[10][0] + 18. * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] , -5.5 + 18. * Hammer2D::Qsi_12P[10][0] + 18. * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[10][0] - 22.5 * Hammer2D::Qsi_12P[10][1] + 40.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 54. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] , -22.5 * Hammer2D::Qsi_12P[10][0] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[10][0] + 4.5 * Hammer2D::Qsi_12P[10][1] - 40.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] , 4.5 * Hammer2D::Qsi_12P[10][0] - 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[10][0] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[10][1] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 27. * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] , 9. - 22.5 * Hammer2D::Qsi_12P[10][0] - 45. * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] + 54. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] + 40.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] } ,
	  { 27. * Hammer2D::Qsi_12P[10][1] - 54. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 27. * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] , 27. * Hammer2D::Qsi_12P[10][0] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] - 54. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[10][1] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] , -4.5 * Hammer2D::Qsi_12P[10][0] + 13.5 * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[10][1] - 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] , -4.5 + 36. * Hammer2D::Qsi_12P[10][1] + 4.5 * Hammer2D::Qsi_12P[10][0] - 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] - 40.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] , -4.5 * Hammer2D::Qsi_12P[10][0] + 27. * Hammer2D::Qsi_12P[10][0] * Hammer2D::Qsi_12P[10][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[10][1] + 13.5 * Hammer2D::Qsi_12P[10][1] * Hammer2D::Qsi_12P[10][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_12P[11][0] + 18. * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] , -5.5 + 18. * Hammer2D::Qsi_12P[11][0] + 18. * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_12P[11][0] - 22.5 * Hammer2D::Qsi_12P[11][1] + 40.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 54. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] , -22.5 * Hammer2D::Qsi_12P[11][0] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_12P[11][0] + 4.5 * Hammer2D::Qsi_12P[11][1] - 40.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] , 4.5 * Hammer2D::Qsi_12P[11][0] - 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_12P[11][0] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_12P[11][1] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 27. * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] , 9. - 22.5 * Hammer2D::Qsi_12P[11][0] - 45. * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] + 54. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] + 40.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] } ,
	  { 27. * Hammer2D::Qsi_12P[11][1] - 54. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 27. * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] , 27. * Hammer2D::Qsi_12P[11][0] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] - 54. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[11][1] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] , -4.5 * Hammer2D::Qsi_12P[11][0] + 13.5 * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][0] } ,
	  {  4.5 * Hammer2D::Qsi_12P[11][1] - 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] , -4.5 + 36. * Hammer2D::Qsi_12P[11][1] + 4.5 * Hammer2D::Qsi_12P[11][0] - 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] - 40.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] } ,
	  { -4.5 * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] , -4.5 * Hammer2D::Qsi_12P[11][0] + 27. * Hammer2D::Qsi_12P[11][0] * Hammer2D::Qsi_12P[11][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_12P[11][1] + 13.5 * Hammer2D::Qsi_12P[11][1] * Hammer2D::Qsi_12P[11][1] } } };

template<> const double O2P2::Geom::Elem::Elem_Tri10_IP<13>::mv_DPsi[13][mv_numNodes][mv_ElDim] = {
	{ { -5.5 + 18. * Hammer2D::Qsi_13P[0][0] + 18. * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] , -5.5 + 18. * Hammer2D::Qsi_13P[0][0] + 18. * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[0][0] - 22.5 * Hammer2D::Qsi_13P[0][1] + 40.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 54. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] , -22.5 * Hammer2D::Qsi_13P[0][0] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[0][0] + 4.5 * Hammer2D::Qsi_13P[0][1] - 40.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] , 4.5 * Hammer2D::Qsi_13P[0][0] - 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[0][0] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[0][1] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 27. * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] , 9. - 22.5 * Hammer2D::Qsi_13P[0][0] - 45. * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] + 54. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] + 40.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] } ,
	  { 27. * Hammer2D::Qsi_13P[0][1] - 54. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 27. * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] , 27. * Hammer2D::Qsi_13P[0][0] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] - 54. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[0][1] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] , -4.5 * Hammer2D::Qsi_13P[0][0] + 13.5 * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[0][1] - 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] , -4.5 + 36. * Hammer2D::Qsi_13P[0][1] + 4.5 * Hammer2D::Qsi_13P[0][0] - 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] - 40.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] , -4.5 * Hammer2D::Qsi_13P[0][0] + 27. * Hammer2D::Qsi_13P[0][0] * Hammer2D::Qsi_13P[0][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[0][1] + 13.5 * Hammer2D::Qsi_13P[0][1] * Hammer2D::Qsi_13P[0][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[1][0] + 18. * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] , -5.5 + 18. * Hammer2D::Qsi_13P[1][0] + 18. * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[1][0] - 22.5 * Hammer2D::Qsi_13P[1][1] + 40.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 54. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] , -22.5 * Hammer2D::Qsi_13P[1][0] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[1][0] + 4.5 * Hammer2D::Qsi_13P[1][1] - 40.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] , 4.5 * Hammer2D::Qsi_13P[1][0] - 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[1][0] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[1][1] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 27. * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] , 9. - 22.5 * Hammer2D::Qsi_13P[1][0] - 45. * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] + 54. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] + 40.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] } ,
	  { 27. * Hammer2D::Qsi_13P[1][1] - 54. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 27. * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] , 27. * Hammer2D::Qsi_13P[1][0] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] - 54. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[1][1] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] , -4.5 * Hammer2D::Qsi_13P[1][0] + 13.5 * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[1][1] - 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] , -4.5 + 36. * Hammer2D::Qsi_13P[1][1] + 4.5 * Hammer2D::Qsi_13P[1][0] - 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] - 40.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] , -4.5 * Hammer2D::Qsi_13P[1][0] + 27. * Hammer2D::Qsi_13P[1][0] * Hammer2D::Qsi_13P[1][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[1][1] + 13.5 * Hammer2D::Qsi_13P[1][1] * Hammer2D::Qsi_13P[1][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[2][0] + 18. * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] , -5.5 + 18. * Hammer2D::Qsi_13P[2][0] + 18. * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[2][0] - 22.5 * Hammer2D::Qsi_13P[2][1] + 40.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 54. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] , -22.5 * Hammer2D::Qsi_13P[2][0] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[2][0] + 4.5 * Hammer2D::Qsi_13P[2][1] - 40.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] , 4.5 * Hammer2D::Qsi_13P[2][0] - 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[2][0] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[2][1] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 27. * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] , 9. - 22.5 * Hammer2D::Qsi_13P[2][0] - 45. * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] + 54. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] + 40.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] } ,
	  { 27. * Hammer2D::Qsi_13P[2][1] - 54. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 27. * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] , 27. * Hammer2D::Qsi_13P[2][0] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] - 54. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[2][1] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] , -4.5 * Hammer2D::Qsi_13P[2][0] + 13.5 * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[2][1] - 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] , -4.5 + 36. * Hammer2D::Qsi_13P[2][1] + 4.5 * Hammer2D::Qsi_13P[2][0] - 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] - 40.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] , -4.5 * Hammer2D::Qsi_13P[2][0] + 27. * Hammer2D::Qsi_13P[2][0] * Hammer2D::Qsi_13P[2][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[2][1] + 13.5 * Hammer2D::Qsi_13P[2][1] * Hammer2D::Qsi_13P[2][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[3][0] + 18. * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] , -5.5 + 18. * Hammer2D::Qsi_13P[3][0] + 18. * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[3][0] - 22.5 * Hammer2D::Qsi_13P[3][1] + 40.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 54. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] , -22.5 * Hammer2D::Qsi_13P[3][0] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[3][0] + 4.5 * Hammer2D::Qsi_13P[3][1] - 40.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] , 4.5 * Hammer2D::Qsi_13P[3][0] - 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[3][0] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[3][1] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 27. * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] , 9. - 22.5 * Hammer2D::Qsi_13P[3][0] - 45. * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] + 54. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] + 40.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] } ,
	  { 27. * Hammer2D::Qsi_13P[3][1] - 54. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 27. * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] , 27. * Hammer2D::Qsi_13P[3][0] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] - 54. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[3][1] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] , -4.5 * Hammer2D::Qsi_13P[3][0] + 13.5 * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[3][1] - 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] , -4.5 + 36. * Hammer2D::Qsi_13P[3][1] + 4.5 * Hammer2D::Qsi_13P[3][0] - 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] - 40.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] , -4.5 * Hammer2D::Qsi_13P[3][0] + 27. * Hammer2D::Qsi_13P[3][0] * Hammer2D::Qsi_13P[3][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[3][1] + 13.5 * Hammer2D::Qsi_13P[3][1] * Hammer2D::Qsi_13P[3][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[4][0] + 18. * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] , -5.5 + 18. * Hammer2D::Qsi_13P[4][0] + 18. * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[4][0] - 22.5 * Hammer2D::Qsi_13P[4][1] + 40.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 54. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] , -22.5 * Hammer2D::Qsi_13P[4][0] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[4][0] + 4.5 * Hammer2D::Qsi_13P[4][1] - 40.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] , 4.5 * Hammer2D::Qsi_13P[4][0] - 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[4][0] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[4][1] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 27. * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] , 9. - 22.5 * Hammer2D::Qsi_13P[4][0] - 45. * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] + 54. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] + 40.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] } ,
	  { 27. * Hammer2D::Qsi_13P[4][1] - 54. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 27. * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] , 27. * Hammer2D::Qsi_13P[4][0] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] - 54. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[4][1] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] , -4.5 * Hammer2D::Qsi_13P[4][0] + 13.5 * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[4][1] - 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] , -4.5 + 36. * Hammer2D::Qsi_13P[4][1] + 4.5 * Hammer2D::Qsi_13P[4][0] - 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] - 40.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] , -4.5 * Hammer2D::Qsi_13P[4][0] + 27. * Hammer2D::Qsi_13P[4][0] * Hammer2D::Qsi_13P[4][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[4][1] + 13.5 * Hammer2D::Qsi_13P[4][1] * Hammer2D::Qsi_13P[4][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[5][0] + 18. * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] , -5.5 + 18. * Hammer2D::Qsi_13P[5][0] + 18. * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[5][0] - 22.5 * Hammer2D::Qsi_13P[5][1] + 40.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 54. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] , -22.5 * Hammer2D::Qsi_13P[5][0] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[5][0] + 4.5 * Hammer2D::Qsi_13P[5][1] - 40.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] , 4.5 * Hammer2D::Qsi_13P[5][0] - 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[5][0] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[5][1] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 27. * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] , 9. - 22.5 * Hammer2D::Qsi_13P[5][0] - 45. * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] + 54. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] + 40.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] } ,
	  { 27. * Hammer2D::Qsi_13P[5][1] - 54. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 27. * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] , 27. * Hammer2D::Qsi_13P[5][0] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] - 54. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[5][1] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] , -4.5 * Hammer2D::Qsi_13P[5][0] + 13.5 * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[5][1] - 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] , -4.5 + 36. * Hammer2D::Qsi_13P[5][1] + 4.5 * Hammer2D::Qsi_13P[5][0] - 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] - 40.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] , -4.5 * Hammer2D::Qsi_13P[5][0] + 27. * Hammer2D::Qsi_13P[5][0] * Hammer2D::Qsi_13P[5][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[5][1] + 13.5 * Hammer2D::Qsi_13P[5][1] * Hammer2D::Qsi_13P[5][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[6][0] + 18. * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] , -5.5 + 18. * Hammer2D::Qsi_13P[6][0] + 18. * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[6][0] - 22.5 * Hammer2D::Qsi_13P[6][1] + 40.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 54. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] , -22.5 * Hammer2D::Qsi_13P[6][0] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[6][0] + 4.5 * Hammer2D::Qsi_13P[6][1] - 40.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] , 4.5 * Hammer2D::Qsi_13P[6][0] - 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[6][0] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[6][1] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 27. * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] , 9. - 22.5 * Hammer2D::Qsi_13P[6][0] - 45. * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] + 54. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] + 40.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] } ,
	  { 27. * Hammer2D::Qsi_13P[6][1] - 54. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 27. * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] , 27. * Hammer2D::Qsi_13P[6][0] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] - 54. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[6][1] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] , -4.5 * Hammer2D::Qsi_13P[6][0] + 13.5 * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[6][1] - 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] , -4.5 + 36. * Hammer2D::Qsi_13P[6][1] + 4.5 * Hammer2D::Qsi_13P[6][0] - 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] - 40.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] , -4.5 * Hammer2D::Qsi_13P[6][0] + 27. * Hammer2D::Qsi_13P[6][0] * Hammer2D::Qsi_13P[6][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[6][1] + 13.5 * Hammer2D::Qsi_13P[6][1] * Hammer2D::Qsi_13P[6][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[7][0] + 18. * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] , -5.5 + 18. * Hammer2D::Qsi_13P[7][0] + 18. * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[7][0] - 22.5 * Hammer2D::Qsi_13P[7][1] + 40.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 54. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] , -22.5 * Hammer2D::Qsi_13P[7][0] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[7][0] + 4.5 * Hammer2D::Qsi_13P[7][1] - 40.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] , 4.5 * Hammer2D::Qsi_13P[7][0] - 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[7][0] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[7][1] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 27. * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] , 9. - 22.5 * Hammer2D::Qsi_13P[7][0] - 45. * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] + 54. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] + 40.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] } ,
	  { 27. * Hammer2D::Qsi_13P[7][1] - 54. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 27. * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] , 27. * Hammer2D::Qsi_13P[7][0] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] - 54. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[7][1] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] , -4.5 * Hammer2D::Qsi_13P[7][0] + 13.5 * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[7][1] - 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] , -4.5 + 36. * Hammer2D::Qsi_13P[7][1] + 4.5 * Hammer2D::Qsi_13P[7][0] - 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] - 40.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] , -4.5 * Hammer2D::Qsi_13P[7][0] + 27. * Hammer2D::Qsi_13P[7][0] * Hammer2D::Qsi_13P[7][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[7][1] + 13.5 * Hammer2D::Qsi_13P[7][1] * Hammer2D::Qsi_13P[7][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[8][0] + 18. * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] , -5.5 + 18. * Hammer2D::Qsi_13P[8][0] + 18. * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[8][0] - 22.5 * Hammer2D::Qsi_13P[8][1] + 40.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 54. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] , -22.5 * Hammer2D::Qsi_13P[8][0] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[8][0] + 4.5 * Hammer2D::Qsi_13P[8][1] - 40.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] , 4.5 * Hammer2D::Qsi_13P[8][0] - 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[8][0] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[8][1] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 27. * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] , 9. - 22.5 * Hammer2D::Qsi_13P[8][0] - 45. * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] + 54. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] + 40.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] } ,
	  { 27. * Hammer2D::Qsi_13P[8][1] - 54. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 27. * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] , 27. * Hammer2D::Qsi_13P[8][0] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] - 54. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[8][1] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] , -4.5 * Hammer2D::Qsi_13P[8][0] + 13.5 * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[8][1] - 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] , -4.5 + 36. * Hammer2D::Qsi_13P[8][1] + 4.5 * Hammer2D::Qsi_13P[8][0] - 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] - 40.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] , -4.5 * Hammer2D::Qsi_13P[8][0] + 27. * Hammer2D::Qsi_13P[8][0] * Hammer2D::Qsi_13P[8][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[8][1] + 13.5 * Hammer2D::Qsi_13P[8][1] * Hammer2D::Qsi_13P[8][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[9][0] + 18. * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] , -5.5 + 18. * Hammer2D::Qsi_13P[9][0] + 18. * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[9][0] - 22.5 * Hammer2D::Qsi_13P[9][1] + 40.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 54. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] , -22.5 * Hammer2D::Qsi_13P[9][0] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[9][0] + 4.5 * Hammer2D::Qsi_13P[9][1] - 40.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] , 4.5 * Hammer2D::Qsi_13P[9][0] - 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[9][0] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[9][1] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 27. * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] , 9. - 22.5 * Hammer2D::Qsi_13P[9][0] - 45. * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] + 54. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] + 40.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] } ,
	  { 27. * Hammer2D::Qsi_13P[9][1] - 54. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 27. * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] , 27. * Hammer2D::Qsi_13P[9][0] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] - 54. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[9][1] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] , -4.5 * Hammer2D::Qsi_13P[9][0] + 13.5 * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[9][1] - 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] , -4.5 + 36. * Hammer2D::Qsi_13P[9][1] + 4.5 * Hammer2D::Qsi_13P[9][0] - 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] - 40.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] , -4.5 * Hammer2D::Qsi_13P[9][0] + 27. * Hammer2D::Qsi_13P[9][0] * Hammer2D::Qsi_13P[9][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[9][1] + 13.5 * Hammer2D::Qsi_13P[9][1] * Hammer2D::Qsi_13P[9][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[10][0] + 18. * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] , -5.5 + 18. * Hammer2D::Qsi_13P[10][0] + 18. * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[10][0] - 22.5 * Hammer2D::Qsi_13P[10][1] + 40.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 54. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] , -22.5 * Hammer2D::Qsi_13P[10][0] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[10][0] + 4.5 * Hammer2D::Qsi_13P[10][1] - 40.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] , 4.5 * Hammer2D::Qsi_13P[10][0] - 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[10][0] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[10][1] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 27. * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] , 9. - 22.5 * Hammer2D::Qsi_13P[10][0] - 45. * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] + 54. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] + 40.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] } ,
	  { 27. * Hammer2D::Qsi_13P[10][1] - 54. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 27. * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] , 27. * Hammer2D::Qsi_13P[10][0] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] - 54. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[10][1] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] , -4.5 * Hammer2D::Qsi_13P[10][0] + 13.5 * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[10][1] - 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] , -4.5 + 36. * Hammer2D::Qsi_13P[10][1] + 4.5 * Hammer2D::Qsi_13P[10][0] - 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] - 40.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] , -4.5 * Hammer2D::Qsi_13P[10][0] + 27. * Hammer2D::Qsi_13P[10][0] * Hammer2D::Qsi_13P[10][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[10][1] + 13.5 * Hammer2D::Qsi_13P[10][1] * Hammer2D::Qsi_13P[10][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[11][0] + 18. * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] , -5.5 + 18. * Hammer2D::Qsi_13P[11][0] + 18. * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[11][0] - 22.5 * Hammer2D::Qsi_13P[11][1] + 40.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 54. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] , -22.5 * Hammer2D::Qsi_13P[11][0] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[11][0] + 4.5 * Hammer2D::Qsi_13P[11][1] - 40.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] , 4.5 * Hammer2D::Qsi_13P[11][0] - 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[11][0] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[11][1] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 27. * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] , 9. - 22.5 * Hammer2D::Qsi_13P[11][0] - 45. * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] + 54. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] + 40.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] } ,
	  { 27. * Hammer2D::Qsi_13P[11][1] - 54. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 27. * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] , 27. * Hammer2D::Qsi_13P[11][0] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] - 54. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[11][1] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] , -4.5 * Hammer2D::Qsi_13P[11][0] + 13.5 * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[11][1] - 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] , -4.5 + 36. * Hammer2D::Qsi_13P[11][1] + 4.5 * Hammer2D::Qsi_13P[11][0] - 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] - 40.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] , -4.5 * Hammer2D::Qsi_13P[11][0] + 27. * Hammer2D::Qsi_13P[11][0] * Hammer2D::Qsi_13P[11][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[11][1] + 13.5 * Hammer2D::Qsi_13P[11][1] * Hammer2D::Qsi_13P[11][1] } } ,

	{ { -5.5 + 18. * Hammer2D::Qsi_13P[12][0] + 18. * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] , -5.5 + 18. * Hammer2D::Qsi_13P[12][0] + 18. * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] } ,
	  { 9. - 45. * Hammer2D::Qsi_13P[12][0] - 22.5 * Hammer2D::Qsi_13P[12][1] + 40.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 54. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] , -22.5 * Hammer2D::Qsi_13P[12][0] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] } ,
	  { -4.5 + 36. * Hammer2D::Qsi_13P[12][0] + 4.5 * Hammer2D::Qsi_13P[12][1] - 40.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] , 4.5 * Hammer2D::Qsi_13P[12][0] - 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] } ,
	  { 1. - 9. * Hammer2D::Qsi_13P[12][0] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] , 0. } ,
	  { -22.5 * Hammer2D::Qsi_13P[12][1] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 27. * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] , 9. - 22.5 * Hammer2D::Qsi_13P[12][0] - 45. * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] + 54. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] + 40.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] } ,
	  { 27. * Hammer2D::Qsi_13P[12][1] - 54. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 27. * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] , 27. * Hammer2D::Qsi_13P[12][0] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] - 54. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[12][1] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] , -4.5 * Hammer2D::Qsi_13P[12][0] + 13.5 * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][0] } ,
	  {  4.5 * Hammer2D::Qsi_13P[12][1] - 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] , -4.5 + 36. * Hammer2D::Qsi_13P[12][1] + 4.5 * Hammer2D::Qsi_13P[12][0] - 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] - 40.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] } ,
	  { -4.5 * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] , -4.5 * Hammer2D::Qsi_13P[12][0] + 27. * Hammer2D::Qsi_13P[12][0] * Hammer2D::Qsi_13P[12][1] } ,
	  { 0. , 1. - 9. * Hammer2D::Qsi_13P[12][1] + 13.5 * Hammer2D::Qsi_13P[12][1] * Hammer2D::Qsi_13P[12][1] } } };
