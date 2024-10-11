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
// Plane element, with cubic / linear interpolation functions, rectangular shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

namespace O2P2 {
	namespace Geom {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Rect8
			  *
			  * @brief Quadrangular cubic / linear element with 8 nodes.
			  * @details Plane element, with cubic interpolation functions in one direction and linear function in the other, rectangular shaped.
			  * @image html Elem_Quad8.png height=300
			  */
			class Elem_Rect8 : public ElemPlane
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Rect8() = delete;

			public:
				/** Constructor for rectangular cubic / linear elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to PlaneSection class.
				  */
				explicit Elem_Rect8(std::shared_ptr<O2P2::Geom::Material>& Material, std::shared_ptr<O2P2::Geom::PlaneSection>& Section)
					: ElemPlane(Material, Section) { }

				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;

					// It is divided into 3 linear elements
					msg << "3 1 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " "
						<< this->mv_Conect[4]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << this->mv_Conect[1]->mv_index + add << " " << this->mv_Conect[2]->mv_index + add << " "
						<< this->mv_Conect[5]->mv_index + add << " " << this->mv_Conect[6]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << this->mv_Conect[2]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " "
						<< this->mv_Conect[6]->mv_index + add << " " << this->mv_Conect[7]->mv_index + add << " "
						<< this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "3 1 " << (1 + add) << " " << (2 + add) << " " << (5 + add) << " " << (6 + add)
						<< " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << (2 + add) << " " << (3 + add) << " " << (6 + add) << " " << (7 + add)
						<< " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << (3 + add) << " " << (4 + add) << " " << (7 + add) << " " << (8 + add)
						<< " " << this->mv_Mat->mv_index << "\n";
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
					std::array<double, mv_ElDim> new_xsi = {};
					for (int i = 0; i < mv_ElDim; ++i) {
						new_xsi.at(i) = *(xsi + i);
					}

					const auto [min, max] = std::minmax_element(new_xsi.begin(), new_xsi.end());
					if (*max < 1.000001 && *min > -1.000001) return true;
					return false;
				}

			private:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			private:
				/** @brief Number of Nodes */
				static const int mv_numNodes{ 8 };

				/** @brief Number of Integration Points */
				static const int mv_numIP{ 8 };

				/** @brief Number of Faces (for output purposes) */
				static const int mv_numFaces{ 3 };

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
inline std::vector<double> O2P2::Geom::Elem::Elem_Rect8::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(8);

	mi_Psi.at(0) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] - Point[0] * Point[1] + Point[0] + Point[1] - 1.);
	mi_Psi.at(1) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] + 27. * Point[0] * Point[1] - 27. * Point[0] - 9. * Point[1] + 9.);
	mi_Psi.at(2) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] - 27. * Point[0] * Point[1] + 27. * Point[0] - 9. * Point[1] + 9.);
	mi_Psi.at(3) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] + Point[0] * Point[1] - Point[0] + Point[1] - 1.);
	mi_Psi.at(4) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] + Point[0] * Point[1] + Point[0] - Point[1] - 1.);
	mi_Psi.at(5) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 27. * Point[0] + 9. * Point[1] + 9.);
	mi_Psi.at(6) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] * Point[1] - 9. * Point[0] * Point[0] + 27. * Point[0] * Point[1] + 27. * Point[0] + 9. * Point[1] + 9.);
	mi_Psi.at(7) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] * Point[1] + 9. * Point[0] * Point[0] - Point[0] * Point[1] - Point[0] - Point[1] - 1.);

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Geom::Elem::Elem_Rect8::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(8 * 2);

	mi_DPsi.at(0) = 0.03125 * (27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] - 18. * Point[0] * Point[1] + 18. * Point[0] - Point[1] + 1.);
	mi_DPsi.at(1) = 0.03125 * (-81. * Point[0] * Point[0] * Point[1] + 81. * Point[0] * Point[0] + 18. * Point[0] * Point[1] - 18. * Point[0] + 27. * Point[1] - 27.);
	mi_DPsi.at(2) = 0.03125 * (81. * Point[0] * Point[0] * Point[1] - 81. * Point[0] * Point[0] + 18. * Point[0] * Point[1] - 18. * Point[0] - 27. * Point[1] + 27.);
	mi_DPsi.at(3) = 0.03125 * (-27. * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] - 18. * Point[0] * Point[1] + 18. * Point[0] + Point[1] - 1.);
	mi_DPsi.at(4) = 0.03125 * (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 18. * Point[0] + Point[1] + 1.);
	mi_DPsi.at(5) = 0.03125 * (81. * Point[0] * Point[0] * Point[1] + 81. * Point[0] * Point[0] - 18. * Point[0] * Point[1] - 18. * Point[0] - 27. * Point[1] - 27.);
	mi_DPsi.at(6) = 0.03125 * (-81. * Point[0] * Point[0] * Point[1] - 81. * Point[0] * Point[0] - 18. * Point[0] * Point[1] - 18. * Point[0] + 27. * Point[1] + 27.);
	mi_DPsi.at(7) = 0.03125 * (27. * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 18. * Point[0] - Point[1] - 1.);

	mi_DPsi.at(8) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] - Point[0] + 1.);
	mi_DPsi.at(9) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] + 27. * Point[0] - 9.);
	mi_DPsi.at(10) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] - 27. * Point[0] - 9.);
	mi_DPsi.at(11) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] + Point[0] + 1.);
	mi_DPsi.at(12) = 0.03125 * (-9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] + Point[0] - 1.);
	mi_DPsi.at(13) = 0.03125 * (27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] - 27. * Point[0] + 9.);
	mi_DPsi.at(14) = 0.03125 * (-27. * Point[0] * Point[0] * Point[0] - 9. * Point[0] * Point[0] + 27. * Point[0] + 9.);
	mi_DPsi.at(15) = 0.03125 * (9. * Point[0] * Point[0] * Point[0] + 9. * Point[0] * Point[0] - Point[0] - 1.);

	return mi_DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Geom::Elem::Elem_Rect8::setGeomProperties() {

	const int nVertices = 4;
	const int mi_Dim = mv_Conect.at(0)->getDIM();	// Dimensionality of vector space (2D or 3D)

	mv_Centroid = std::make_unique<double[]>(mi_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Geom::Node*, nVertices> vertices;
	vertices[0] = mv_Conect[0].get();
	vertices[1] = mv_Conect[3].get();
	vertices[2] = mv_Conect[4].get();
	vertices[3] = mv_Conect[7].get();

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
inline std::vector<double> O2P2::Geom::Elem::Elem_Rect8::getValueOnIPs(const double* value) {

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
inline const double O2P2::Geom::Elem::Elem_Rect8::mv_xsi[mv_numIP][mv_ElDim] = {
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
inline const double O2P2::Geom::Elem::Elem_Rect8::mv_weight[mv_numIP] = { 0.347854845137454, 0.347854845137454, 0.347854845137454, 0.347854845137454, 0.652145154862546, 0.652145154862546, 0.652145154862546, 0.652145154862546 };

// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double O2P2::Geom::Elem::Elem_Rect8::mv_Psi[mv_numIP][mv_numNodes] = {
	{ 0.03125 * (9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 9. * mv_xsi[0][0] * mv_xsi[0][0] - mv_xsi[0][0] * mv_xsi[0][1] + mv_xsi[0][0] + mv_xsi[0][1] - 1.),
	  0.03125 * (-27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 9. * mv_xsi[0][0] * mv_xsi[0][0] + 27. * mv_xsi[0][0] * mv_xsi[0][1] - 27. * mv_xsi[0][0] - 9. * mv_xsi[0][1] + 9.),
	  0.03125 * (27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 9. * mv_xsi[0][0] * mv_xsi[0][0] - 27. * mv_xsi[0][0] * mv_xsi[0][1] + 27. * mv_xsi[0][0] - 9. * mv_xsi[0][1] + 9.),
	  0.03125 * (-9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 9. * mv_xsi[0][0] * mv_xsi[0][0] + mv_xsi[0][0] * mv_xsi[0][1] - mv_xsi[0][0] + mv_xsi[0][1] - 1.),
	  0.03125 * (-9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 9. * mv_xsi[0][0] * mv_xsi[0][0] + mv_xsi[0][0] * mv_xsi[0][1] + mv_xsi[0][0] - mv_xsi[0][1] - 1.),
	  0.03125 * (27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 9. * mv_xsi[0][0] * mv_xsi[0][0] - 27. * mv_xsi[0][0] * mv_xsi[0][1] - 27. * mv_xsi[0][0] + 9. * mv_xsi[0][1] + 9.),
	  0.03125 * (-27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 9. * mv_xsi[0][0] * mv_xsi[0][0] + 27. * mv_xsi[0][0] * mv_xsi[0][1] + 27. * mv_xsi[0][0] + 9. * mv_xsi[0][1] + 9.),
	  0.03125 * (9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 9. * mv_xsi[0][0] * mv_xsi[0][0] - mv_xsi[0][0] * mv_xsi[0][1] - mv_xsi[0][0] - mv_xsi[0][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 9. * mv_xsi[1][0] * mv_xsi[1][0] - mv_xsi[1][0] * mv_xsi[1][1] + mv_xsi[1][0] + mv_xsi[1][1] - 1.),
	  0.03125 * (-27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 9. * mv_xsi[1][0] * mv_xsi[1][0] + 27. * mv_xsi[1][0] * mv_xsi[1][1] - 27. * mv_xsi[1][0] - 9. * mv_xsi[1][1] + 9.),
	  0.03125 * (27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 9. * mv_xsi[1][0] * mv_xsi[1][0] - 27. * mv_xsi[1][0] * mv_xsi[1][1] + 27. * mv_xsi[1][0] - 9. * mv_xsi[1][1] + 9.),
	  0.03125 * (-9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 9. * mv_xsi[1][0] * mv_xsi[1][0] + mv_xsi[1][0] * mv_xsi[1][1] - mv_xsi[1][0] + mv_xsi[1][1] - 1.),
	  0.03125 * (-9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 9. * mv_xsi[1][0] * mv_xsi[1][0] + mv_xsi[1][0] * mv_xsi[1][1] + mv_xsi[1][0] - mv_xsi[1][1] - 1.),
	  0.03125 * (27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 9. * mv_xsi[1][0] * mv_xsi[1][0] - 27. * mv_xsi[1][0] * mv_xsi[1][1] - 27. * mv_xsi[1][0] + 9. * mv_xsi[1][1] + 9.),
	  0.03125 * (-27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 9. * mv_xsi[1][0] * mv_xsi[1][0] + 27. * mv_xsi[1][0] * mv_xsi[1][1] + 27. * mv_xsi[1][0] + 9. * mv_xsi[1][1] + 9.),
	  0.03125 * (9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 9. * mv_xsi[1][0] * mv_xsi[1][0] - mv_xsi[1][0] * mv_xsi[1][1] - mv_xsi[1][0] - mv_xsi[1][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 9. * mv_xsi[2][0] * mv_xsi[2][0] - mv_xsi[2][0] * mv_xsi[2][1] + mv_xsi[2][0] + mv_xsi[2][1] - 1.),
	  0.03125 * (-27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 9. * mv_xsi[2][0] * mv_xsi[2][0] + 27. * mv_xsi[2][0] * mv_xsi[2][1] - 27. * mv_xsi[2][0] - 9. * mv_xsi[2][1] + 9.),
	  0.03125 * (27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 9. * mv_xsi[2][0] * mv_xsi[2][0] - 27. * mv_xsi[2][0] * mv_xsi[2][1] + 27. * mv_xsi[2][0] - 9. * mv_xsi[2][1] + 9.),
	  0.03125 * (-9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 9. * mv_xsi[2][0] * mv_xsi[2][0] + mv_xsi[2][0] * mv_xsi[2][1] - mv_xsi[2][0] + mv_xsi[2][1] - 1.),
	  0.03125 * (-9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 9. * mv_xsi[2][0] * mv_xsi[2][0] + mv_xsi[2][0] * mv_xsi[2][1] + mv_xsi[2][0] - mv_xsi[2][1] - 1.),
	  0.03125 * (27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 9. * mv_xsi[2][0] * mv_xsi[2][0] - 27. * mv_xsi[2][0] * mv_xsi[2][1] - 27. * mv_xsi[2][0] + 9. * mv_xsi[2][1] + 9.),
	  0.03125 * (-27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 9. * mv_xsi[2][0] * mv_xsi[2][0] + 27. * mv_xsi[2][0] * mv_xsi[2][1] + 27. * mv_xsi[2][0] + 9. * mv_xsi[2][1] + 9.),
	  0.03125 * (9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 9. * mv_xsi[2][0] * mv_xsi[2][0] - mv_xsi[2][0] * mv_xsi[2][1] - mv_xsi[2][0] - mv_xsi[2][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 9. * mv_xsi[3][0] * mv_xsi[3][0] - mv_xsi[3][0] * mv_xsi[3][1] + mv_xsi[3][0] + mv_xsi[3][1] - 1.),
	  0.03125 * (-27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 9. * mv_xsi[3][0] * mv_xsi[3][0] + 27. * mv_xsi[3][0] * mv_xsi[3][1] - 27. * mv_xsi[3][0] - 9. * mv_xsi[3][1] + 9.),
	  0.03125 * (27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 9. * mv_xsi[3][0] * mv_xsi[3][0] - 27. * mv_xsi[3][0] * mv_xsi[3][1] + 27. * mv_xsi[3][0] - 9. * mv_xsi[3][1] + 9.),
	  0.03125 * (-9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 9. * mv_xsi[3][0] * mv_xsi[3][0] + mv_xsi[3][0] * mv_xsi[3][1] - mv_xsi[3][0] + mv_xsi[3][1] - 1.),
	  0.03125 * (-9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 9. * mv_xsi[3][0] * mv_xsi[3][0] + mv_xsi[3][0] * mv_xsi[3][1] + mv_xsi[3][0] - mv_xsi[3][1] - 1.),
	  0.03125 * (27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 9. * mv_xsi[3][0] * mv_xsi[3][0] - 27. * mv_xsi[3][0] * mv_xsi[3][1] - 27. * mv_xsi[3][0] + 9. * mv_xsi[3][1] + 9.),
	  0.03125 * (-27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 9. * mv_xsi[3][0] * mv_xsi[3][0] + 27. * mv_xsi[3][0] * mv_xsi[3][1] + 27. * mv_xsi[3][0] + 9. * mv_xsi[3][1] + 9.),
	  0.03125 * (9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 9. * mv_xsi[3][0] * mv_xsi[3][0] - mv_xsi[3][0] * mv_xsi[3][1] - mv_xsi[3][0] - mv_xsi[3][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 9. * mv_xsi[4][0] * mv_xsi[4][0] - mv_xsi[4][0] * mv_xsi[4][1] + mv_xsi[4][0] + mv_xsi[4][1] - 1.),
	  0.03125 * (-27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 9. * mv_xsi[4][0] * mv_xsi[4][0] + 27. * mv_xsi[4][0] * mv_xsi[4][1] - 27. * mv_xsi[4][0] - 9. * mv_xsi[4][1] + 9.),
	  0.03125 * (27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 9. * mv_xsi[4][0] * mv_xsi[4][0] - 27. * mv_xsi[4][0] * mv_xsi[4][1] + 27. * mv_xsi[4][0] - 9. * mv_xsi[4][1] + 9.),
	  0.03125 * (-9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 9. * mv_xsi[4][0] * mv_xsi[4][0] + mv_xsi[4][0] * mv_xsi[4][1] - mv_xsi[4][0] + mv_xsi[4][1] - 1.),
	  0.03125 * (-9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 9. * mv_xsi[4][0] * mv_xsi[4][0] + mv_xsi[4][0] * mv_xsi[4][1] + mv_xsi[4][0] - mv_xsi[4][1] - 1.),
	  0.03125 * (27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 9. * mv_xsi[4][0] * mv_xsi[4][0] - 27. * mv_xsi[4][0] * mv_xsi[4][1] - 27. * mv_xsi[4][0] + 9. * mv_xsi[4][1] + 9.),
	  0.03125 * (-27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 9. * mv_xsi[4][0] * mv_xsi[4][0] + 27. * mv_xsi[4][0] * mv_xsi[4][1] + 27. * mv_xsi[4][0] + 9. * mv_xsi[4][1] + 9.),
	  0.03125 * (9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 9. * mv_xsi[4][0] * mv_xsi[4][0] - mv_xsi[4][0] * mv_xsi[4][1] - mv_xsi[4][0] - mv_xsi[4][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 9. * mv_xsi[5][0] * mv_xsi[5][0] - mv_xsi[5][0] * mv_xsi[5][1] + mv_xsi[5][0] + mv_xsi[5][1] - 1.),
	  0.03125 * (-27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 9. * mv_xsi[5][0] * mv_xsi[5][0] + 27. * mv_xsi[5][0] * mv_xsi[5][1] - 27. * mv_xsi[5][0] - 9. * mv_xsi[5][1] + 9.),
	  0.03125 * (27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 9. * mv_xsi[5][0] * mv_xsi[5][0] - 27. * mv_xsi[5][0] * mv_xsi[5][1] + 27. * mv_xsi[5][0] - 9. * mv_xsi[5][1] + 9.),
	  0.03125 * (-9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 9. * mv_xsi[5][0] * mv_xsi[5][0] + mv_xsi[5][0] * mv_xsi[5][1] - mv_xsi[5][0] + mv_xsi[5][1] - 1.),
	  0.03125 * (-9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 9. * mv_xsi[5][0] * mv_xsi[5][0] + mv_xsi[5][0] * mv_xsi[5][1] + mv_xsi[5][0] - mv_xsi[5][1] - 1.),
	  0.03125 * (27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 9. * mv_xsi[5][0] * mv_xsi[5][0] - 27. * mv_xsi[5][0] * mv_xsi[5][1] - 27. * mv_xsi[5][0] + 9. * mv_xsi[5][1] + 9.),
	  0.03125 * (-27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 9. * mv_xsi[5][0] * mv_xsi[5][0] + 27. * mv_xsi[5][0] * mv_xsi[5][1] + 27. * mv_xsi[5][0] + 9. * mv_xsi[5][1] + 9.),
	  0.03125 * (9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 9. * mv_xsi[5][0] * mv_xsi[5][0] - mv_xsi[5][0] * mv_xsi[5][1] - mv_xsi[5][0] - mv_xsi[5][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 9. * mv_xsi[6][0] * mv_xsi[6][0] - mv_xsi[6][0] * mv_xsi[6][1] + mv_xsi[6][0] + mv_xsi[6][1] - 1.),
	  0.03125 * (-27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 9. * mv_xsi[6][0] * mv_xsi[6][0] + 27. * mv_xsi[6][0] * mv_xsi[6][1] - 27. * mv_xsi[6][0] - 9. * mv_xsi[6][1] + 9.),
	  0.03125 * (27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 9. * mv_xsi[6][0] * mv_xsi[6][0] - 27. * mv_xsi[6][0] * mv_xsi[6][1] + 27. * mv_xsi[6][0] - 9. * mv_xsi[6][1] + 9.),
	  0.03125 * (-9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 9. * mv_xsi[6][0] * mv_xsi[6][0] + mv_xsi[6][0] * mv_xsi[6][1] - mv_xsi[6][0] + mv_xsi[6][1] - 1.),
	  0.03125 * (-9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 9. * mv_xsi[6][0] * mv_xsi[6][0] + mv_xsi[6][0] * mv_xsi[6][1] + mv_xsi[6][0] - mv_xsi[6][1] - 1.),
	  0.03125 * (27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 9. * mv_xsi[6][0] * mv_xsi[6][0] - 27. * mv_xsi[6][0] * mv_xsi[6][1] - 27. * mv_xsi[6][0] + 9. * mv_xsi[6][1] + 9.),
	  0.03125 * (-27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 9. * mv_xsi[6][0] * mv_xsi[6][0] + 27. * mv_xsi[6][0] * mv_xsi[6][1] + 27. * mv_xsi[6][0] + 9. * mv_xsi[6][1] + 9.),
	  0.03125 * (9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 9. * mv_xsi[6][0] * mv_xsi[6][0] - mv_xsi[6][0] * mv_xsi[6][1] - mv_xsi[6][0] - mv_xsi[6][1] - 1.) },

	{ 0.03125 * (9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 9. * mv_xsi[7][0] * mv_xsi[7][0] - mv_xsi[7][0] * mv_xsi[7][1] + mv_xsi[7][0] + mv_xsi[7][1] - 1.),
	  0.03125 * (-27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 9. * mv_xsi[7][0] * mv_xsi[7][0] + 27. * mv_xsi[7][0] * mv_xsi[7][1] - 27. * mv_xsi[7][0] - 9. * mv_xsi[7][1] + 9.),
	  0.03125 * (27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 9. * mv_xsi[7][0] * mv_xsi[7][0] - 27. * mv_xsi[7][0] * mv_xsi[7][1] + 27. * mv_xsi[7][0] - 9. * mv_xsi[7][1] + 9.),
	  0.03125 * (-9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 9. * mv_xsi[7][0] * mv_xsi[7][0] + mv_xsi[7][0] * mv_xsi[7][1] - mv_xsi[7][0] + mv_xsi[7][1] - 1.),
	  0.03125 * (-9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 9. * mv_xsi[7][0] * mv_xsi[7][0] + mv_xsi[7][0] * mv_xsi[7][1] + mv_xsi[7][0] - mv_xsi[7][1] - 1.),
	  0.03125 * (27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 9. * mv_xsi[7][0] * mv_xsi[7][0] - 27. * mv_xsi[7][0] * mv_xsi[7][1] - 27. * mv_xsi[7][0] + 9. * mv_xsi[7][1] + 9.),
	  0.03125 * (-27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 9. * mv_xsi[7][0] * mv_xsi[7][0] + 27. * mv_xsi[7][0] * mv_xsi[7][1] + 27. * mv_xsi[7][0] + 9. * mv_xsi[7][1] + 9.),
	  0.03125 * (9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 9. * mv_xsi[7][0] * mv_xsi[7][0] - mv_xsi[7][0] * mv_xsi[7][1] - mv_xsi[7][0] - mv_xsi[7][1] - 1.) } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
inline const double O2P2::Geom::Elem::Elem_Rect8::mv_DPsi[mv_numIP][mv_numNodes][mv_ElDim] = {
	{ { 0.03125 * (27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 27. * mv_xsi[0][0] * mv_xsi[0][0] - 18. * mv_xsi[0][0] * mv_xsi[0][1] + 18. * mv_xsi[0][0] - mv_xsi[0][1] + 1.),
		0.03125 * (9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] - mv_xsi[0][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 81. * mv_xsi[0][0] * mv_xsi[0][0] + 18. * mv_xsi[0][0] * mv_xsi[0][1] - 18. * mv_xsi[0][0] + 27. * mv_xsi[0][1] - 27.),
		0.03125 * (-27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] + 27. * mv_xsi[0][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 81. * mv_xsi[0][0] * mv_xsi[0][0] + 18. * mv_xsi[0][0] * mv_xsi[0][1] - 18. * mv_xsi[0][0] - 27. * mv_xsi[0][1] + 27.),
		0.03125 * (27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] - 27. * mv_xsi[0][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 27. * mv_xsi[0][0] * mv_xsi[0][0] - 18. * mv_xsi[0][0] * mv_xsi[0][1] + 18. * mv_xsi[0][0] + mv_xsi[0][1] - 1.),
		0.03125 * (-9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] + mv_xsi[0][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 27. * mv_xsi[0][0] * mv_xsi[0][0] + 18. * mv_xsi[0][0] * mv_xsi[0][1] + 18. * mv_xsi[0][0] + mv_xsi[0][1] + 1.),
		0.03125 * (-9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] + mv_xsi[0][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 81. * mv_xsi[0][0] * mv_xsi[0][0] - 18. * mv_xsi[0][0] * mv_xsi[0][1] - 18. * mv_xsi[0][0] - 27. * mv_xsi[0][1] - 27.),
		0.03125 * (27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] - 27. * mv_xsi[0][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] - 81. * mv_xsi[0][0] * mv_xsi[0][0] - 18. * mv_xsi[0][0] * mv_xsi[0][1] - 18. * mv_xsi[0][0] + 27. * mv_xsi[0][1] + 27.),
		0.03125 * (-27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] - 9. * mv_xsi[0][0] * mv_xsi[0][0] + 27. * mv_xsi[0][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][1] + 27. * mv_xsi[0][0] * mv_xsi[0][0] + 18. * mv_xsi[0][0] * mv_xsi[0][1] + 18. * mv_xsi[0][0] - mv_xsi[0][1] - 1.),
		0.03125 * (9. * mv_xsi[0][0] * mv_xsi[0][0] * mv_xsi[0][0] + 9. * mv_xsi[0][0] * mv_xsi[0][0] - mv_xsi[0][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 27. * mv_xsi[1][0] * mv_xsi[1][0] - 18. * mv_xsi[1][0] * mv_xsi[1][1] + 18. * mv_xsi[1][0] - mv_xsi[1][1] + 1.),
		0.03125 * (9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] - mv_xsi[1][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 81. * mv_xsi[1][0] * mv_xsi[1][0] + 18. * mv_xsi[1][0] * mv_xsi[1][1] - 18. * mv_xsi[1][0] + 27. * mv_xsi[1][1] - 27.),
		0.03125 * (-27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] + 27. * mv_xsi[1][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 81. * mv_xsi[1][0] * mv_xsi[1][0] + 18. * mv_xsi[1][0] * mv_xsi[1][1] - 18. * mv_xsi[1][0] - 27. * mv_xsi[1][1] + 27.),
		0.03125 * (27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] - 27. * mv_xsi[1][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 27. * mv_xsi[1][0] * mv_xsi[1][0] - 18. * mv_xsi[1][0] * mv_xsi[1][1] + 18. * mv_xsi[1][0] + mv_xsi[1][1] - 1.),
		0.03125 * (-9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] + mv_xsi[1][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 27. * mv_xsi[1][0] * mv_xsi[1][0] + 18. * mv_xsi[1][0] * mv_xsi[1][1] + 18. * mv_xsi[1][0] + mv_xsi[1][1] + 1.),
		0.03125 * (-9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] + mv_xsi[1][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 81. * mv_xsi[1][0] * mv_xsi[1][0] - 18. * mv_xsi[1][0] * mv_xsi[1][1] - 18. * mv_xsi[1][0] - 27. * mv_xsi[1][1] - 27.),
		0.03125 * (27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] - 27. * mv_xsi[1][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] - 81. * mv_xsi[1][0] * mv_xsi[1][0] - 18. * mv_xsi[1][0] * mv_xsi[1][1] - 18. * mv_xsi[1][0] + 27. * mv_xsi[1][1] + 27.),
		0.03125 * (-27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] - 9. * mv_xsi[1][0] * mv_xsi[1][0] + 27. * mv_xsi[1][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][1] + 27. * mv_xsi[1][0] * mv_xsi[1][0] + 18. * mv_xsi[1][0] * mv_xsi[1][1] + 18. * mv_xsi[1][0] - mv_xsi[1][1] - 1.),
		0.03125 * (9. * mv_xsi[1][0] * mv_xsi[1][0] * mv_xsi[1][0] + 9. * mv_xsi[1][0] * mv_xsi[1][0] - mv_xsi[1][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 27. * mv_xsi[2][0] * mv_xsi[2][0] - 18. * mv_xsi[2][0] * mv_xsi[2][1] + 18. * mv_xsi[2][0] - mv_xsi[2][1] + 1.),
		0.03125 * (9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] - mv_xsi[2][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 81. * mv_xsi[2][0] * mv_xsi[2][0] + 18. * mv_xsi[2][0] * mv_xsi[2][1] - 18. * mv_xsi[2][0] + 27. * mv_xsi[2][1] - 27.),
		0.03125 * (-27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] + 27. * mv_xsi[2][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 81. * mv_xsi[2][0] * mv_xsi[2][0] + 18. * mv_xsi[2][0] * mv_xsi[2][1] - 18. * mv_xsi[2][0] - 27. * mv_xsi[2][1] + 27.),
		0.03125 * (27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] - 27. * mv_xsi[2][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 27. * mv_xsi[2][0] * mv_xsi[2][0] - 18. * mv_xsi[2][0] * mv_xsi[2][1] + 18. * mv_xsi[2][0] + mv_xsi[2][1] - 1.),
		0.03125 * (-9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] + mv_xsi[2][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 27. * mv_xsi[2][0] * mv_xsi[2][0] + 18. * mv_xsi[2][0] * mv_xsi[2][1] + 18. * mv_xsi[2][0] + mv_xsi[2][1] + 1.),
		0.03125 * (-9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] + mv_xsi[2][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 81. * mv_xsi[2][0] * mv_xsi[2][0] - 18. * mv_xsi[2][0] * mv_xsi[2][1] - 18. * mv_xsi[2][0] - 27. * mv_xsi[2][1] - 27.),
		0.03125 * (27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] - 27. * mv_xsi[2][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] - 81. * mv_xsi[2][0] * mv_xsi[2][0] - 18. * mv_xsi[2][0] * mv_xsi[2][1] - 18. * mv_xsi[2][0] + 27. * mv_xsi[2][1] + 27.),
		0.03125 * (-27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] - 9. * mv_xsi[2][0] * mv_xsi[2][0] + 27. * mv_xsi[2][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][1] + 27. * mv_xsi[2][0] * mv_xsi[2][0] + 18. * mv_xsi[2][0] * mv_xsi[2][1] + 18. * mv_xsi[2][0] - mv_xsi[2][1] - 1.),
		0.03125 * (9. * mv_xsi[2][0] * mv_xsi[2][0] * mv_xsi[2][0] + 9. * mv_xsi[2][0] * mv_xsi[2][0] - mv_xsi[2][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 27. * mv_xsi[3][0] * mv_xsi[3][0] - 18. * mv_xsi[3][0] * mv_xsi[3][1] + 18. * mv_xsi[3][0] - mv_xsi[3][1] + 1.),
		0.03125 * (9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] - mv_xsi[3][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 81. * mv_xsi[3][0] * mv_xsi[3][0] + 18. * mv_xsi[3][0] * mv_xsi[3][1] - 18. * mv_xsi[3][0] + 27. * mv_xsi[3][1] - 27.),
		0.03125 * (-27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] + 27. * mv_xsi[3][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 81. * mv_xsi[3][0] * mv_xsi[3][0] + 18. * mv_xsi[3][0] * mv_xsi[3][1] - 18. * mv_xsi[3][0] - 27. * mv_xsi[3][1] + 27.),
		0.03125 * (27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] - 27. * mv_xsi[3][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 27. * mv_xsi[3][0] * mv_xsi[3][0] - 18. * mv_xsi[3][0] * mv_xsi[3][1] + 18. * mv_xsi[3][0] + mv_xsi[3][1] - 1.),
		0.03125 * (-9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] + mv_xsi[3][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 27. * mv_xsi[3][0] * mv_xsi[3][0] + 18. * mv_xsi[3][0] * mv_xsi[3][1] + 18. * mv_xsi[3][0] + mv_xsi[3][1] + 1.),
		0.03125 * (-9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] + mv_xsi[3][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 81. * mv_xsi[3][0] * mv_xsi[3][0] - 18. * mv_xsi[3][0] * mv_xsi[3][1] - 18. * mv_xsi[3][0] - 27. * mv_xsi[3][1] - 27.),
		0.03125 * (27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] - 27. * mv_xsi[3][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] - 81. * mv_xsi[3][0] * mv_xsi[3][0] - 18. * mv_xsi[3][0] * mv_xsi[3][1] - 18. * mv_xsi[3][0] + 27. * mv_xsi[3][1] + 27.),
		0.03125 * (-27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] - 9. * mv_xsi[3][0] * mv_xsi[3][0] + 27. * mv_xsi[3][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][1] + 27. * mv_xsi[3][0] * mv_xsi[3][0] + 18. * mv_xsi[3][0] * mv_xsi[3][1] + 18. * mv_xsi[3][0] - mv_xsi[3][1] - 1.),
		0.03125 * (9. * mv_xsi[3][0] * mv_xsi[3][0] * mv_xsi[3][0] + 9. * mv_xsi[3][0] * mv_xsi[3][0] - mv_xsi[3][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 27. * mv_xsi[4][0] * mv_xsi[4][0] - 18. * mv_xsi[4][0] * mv_xsi[4][1] + 18. * mv_xsi[4][0] - mv_xsi[4][1] + 1.),
		0.03125 * (9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] - mv_xsi[4][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 81. * mv_xsi[4][0] * mv_xsi[4][0] + 18. * mv_xsi[4][0] * mv_xsi[4][1] - 18. * mv_xsi[4][0] + 27. * mv_xsi[4][1] - 27.),
		0.03125 * (-27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] + 27. * mv_xsi[4][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 81. * mv_xsi[4][0] * mv_xsi[4][0] + 18. * mv_xsi[4][0] * mv_xsi[4][1] - 18. * mv_xsi[4][0] - 27. * mv_xsi[4][1] + 27.),
		0.03125 * (27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] - 27. * mv_xsi[4][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 27. * mv_xsi[4][0] * mv_xsi[4][0] - 18. * mv_xsi[4][0] * mv_xsi[4][1] + 18. * mv_xsi[4][0] + mv_xsi[4][1] - 1.),
		0.03125 * (-9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] + mv_xsi[4][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 27. * mv_xsi[4][0] * mv_xsi[4][0] + 18. * mv_xsi[4][0] * mv_xsi[4][1] + 18. * mv_xsi[4][0] + mv_xsi[4][1] + 1.),
		0.03125 * (-9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] + mv_xsi[4][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 81. * mv_xsi[4][0] * mv_xsi[4][0] - 18. * mv_xsi[4][0] * mv_xsi[4][1] - 18. * mv_xsi[4][0] - 27. * mv_xsi[4][1] - 27.),
		0.03125 * (27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] - 27. * mv_xsi[4][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] - 81. * mv_xsi[4][0] * mv_xsi[4][0] - 18. * mv_xsi[4][0] * mv_xsi[4][1] - 18. * mv_xsi[4][0] + 27. * mv_xsi[4][1] + 27.),
		0.03125 * (-27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] - 9. * mv_xsi[4][0] * mv_xsi[4][0] + 27. * mv_xsi[4][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][1] + 27. * mv_xsi[4][0] * mv_xsi[4][0] + 18. * mv_xsi[4][0] * mv_xsi[4][1] + 18. * mv_xsi[4][0] - mv_xsi[4][1] - 1.),
		0.03125 * (9. * mv_xsi[4][0] * mv_xsi[4][0] * mv_xsi[4][0] + 9. * mv_xsi[4][0] * mv_xsi[4][0] - mv_xsi[4][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 27. * mv_xsi[5][0] * mv_xsi[5][0] - 18. * mv_xsi[5][0] * mv_xsi[5][1] + 18. * mv_xsi[5][0] - mv_xsi[5][1] + 1.),
		0.03125 * (9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] - mv_xsi[5][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 81. * mv_xsi[5][0] * mv_xsi[5][0] + 18. * mv_xsi[5][0] * mv_xsi[5][1] - 18. * mv_xsi[5][0] + 27. * mv_xsi[5][1] - 27.),
		0.03125 * (-27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] + 27. * mv_xsi[5][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 81. * mv_xsi[5][0] * mv_xsi[5][0] + 18. * mv_xsi[5][0] * mv_xsi[5][1] - 18. * mv_xsi[5][0] - 27. * mv_xsi[5][1] + 27.),
		0.03125 * (27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] - 27. * mv_xsi[5][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 27. * mv_xsi[5][0] * mv_xsi[5][0] - 18. * mv_xsi[5][0] * mv_xsi[5][1] + 18. * mv_xsi[5][0] + mv_xsi[5][1] - 1.),
		0.03125 * (-9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] + mv_xsi[5][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 27. * mv_xsi[5][0] * mv_xsi[5][0] + 18. * mv_xsi[5][0] * mv_xsi[5][1] + 18. * mv_xsi[5][0] + mv_xsi[5][1] + 1.),
		0.03125 * (-9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] + mv_xsi[5][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 81. * mv_xsi[5][0] * mv_xsi[5][0] - 18. * mv_xsi[5][0] * mv_xsi[5][1] - 18. * mv_xsi[5][0] - 27. * mv_xsi[5][1] - 27.),
		0.03125 * (27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] - 27. * mv_xsi[5][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] - 81. * mv_xsi[5][0] * mv_xsi[5][0] - 18. * mv_xsi[5][0] * mv_xsi[5][1] - 18. * mv_xsi[5][0] + 27. * mv_xsi[5][1] + 27.),
		0.03125 * (-27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] - 9. * mv_xsi[5][0] * mv_xsi[5][0] + 27. * mv_xsi[5][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][1] + 27. * mv_xsi[5][0] * mv_xsi[5][0] + 18. * mv_xsi[5][0] * mv_xsi[5][1] + 18. * mv_xsi[5][0] - mv_xsi[5][1] - 1.),
		0.03125 * (9. * mv_xsi[5][0] * mv_xsi[5][0] * mv_xsi[5][0] + 9. * mv_xsi[5][0] * mv_xsi[5][0] - mv_xsi[5][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 27. * mv_xsi[6][0] * mv_xsi[6][0] - 18. * mv_xsi[6][0] * mv_xsi[6][1] + 18. * mv_xsi[6][0] - mv_xsi[6][1] + 1.),
		0.03125 * (9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] - mv_xsi[6][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 81. * mv_xsi[6][0] * mv_xsi[6][0] + 18. * mv_xsi[6][0] * mv_xsi[6][1] - 18. * mv_xsi[6][0] + 27. * mv_xsi[6][1] - 27.),
		0.03125 * (-27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] + 27. * mv_xsi[6][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 81. * mv_xsi[6][0] * mv_xsi[6][0] + 18. * mv_xsi[6][0] * mv_xsi[6][1] - 18. * mv_xsi[6][0] - 27. * mv_xsi[6][1] + 27.),
		0.03125 * (27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] - 27. * mv_xsi[6][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 27. * mv_xsi[6][0] * mv_xsi[6][0] - 18. * mv_xsi[6][0] * mv_xsi[6][1] + 18. * mv_xsi[6][0] + mv_xsi[6][1] - 1.),
		0.03125 * (-9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] + mv_xsi[6][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 27. * mv_xsi[6][0] * mv_xsi[6][0] + 18. * mv_xsi[6][0] * mv_xsi[6][1] + 18. * mv_xsi[6][0] + mv_xsi[6][1] + 1.),
		0.03125 * (-9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] + mv_xsi[6][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 81. * mv_xsi[6][0] * mv_xsi[6][0] - 18. * mv_xsi[6][0] * mv_xsi[6][1] - 18. * mv_xsi[6][0] - 27. * mv_xsi[6][1] - 27.),
		0.03125 * (27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] - 27. * mv_xsi[6][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] - 81. * mv_xsi[6][0] * mv_xsi[6][0] - 18. * mv_xsi[6][0] * mv_xsi[6][1] - 18. * mv_xsi[6][0] + 27. * mv_xsi[6][1] + 27.),
		0.03125 * (-27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] - 9. * mv_xsi[6][0] * mv_xsi[6][0] + 27. * mv_xsi[6][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][1] + 27. * mv_xsi[6][0] * mv_xsi[6][0] + 18. * mv_xsi[6][0] * mv_xsi[6][1] + 18. * mv_xsi[6][0] - mv_xsi[6][1] - 1.),
		0.03125 * (9. * mv_xsi[6][0] * mv_xsi[6][0] * mv_xsi[6][0] + 9. * mv_xsi[6][0] * mv_xsi[6][0] - mv_xsi[6][0] - 1.) } },

	{ { 0.03125 * (27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 27. * mv_xsi[7][0] * mv_xsi[7][0] - 18. * mv_xsi[7][0] * mv_xsi[7][1] + 18. * mv_xsi[7][0] - mv_xsi[7][1] + 1.),
		0.03125 * (9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] - mv_xsi[7][0] + 1.) },
	  { 0.03125 * (-81. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 81. * mv_xsi[7][0] * mv_xsi[7][0] + 18. * mv_xsi[7][0] * mv_xsi[7][1] - 18. * mv_xsi[7][0] + 27. * mv_xsi[7][1] - 27.),
		0.03125 * (-27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] + 27. * mv_xsi[7][0] - 9.) },
	  { 0.03125 * (81. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 81. * mv_xsi[7][0] * mv_xsi[7][0] + 18. * mv_xsi[7][0] * mv_xsi[7][1] - 18. * mv_xsi[7][0] - 27. * mv_xsi[7][1] + 27.),
		0.03125 * (27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] - 27. * mv_xsi[7][0] - 9.) },
	  { 0.03125 * (-27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 27. * mv_xsi[7][0] * mv_xsi[7][0] - 18. * mv_xsi[7][0] * mv_xsi[7][1] + 18. * mv_xsi[7][0] + mv_xsi[7][1] - 1.),
		0.03125 * (-9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] + mv_xsi[7][0] + 1.) },
	  { 0.03125 * (-27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 27. * mv_xsi[7][0] * mv_xsi[7][0] + 18. * mv_xsi[7][0] * mv_xsi[7][1] + 18. * mv_xsi[7][0] + mv_xsi[7][1] + 1.),
		0.03125 * (-9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] + mv_xsi[7][0] - 1.) },
	  { 0.03125 * (81. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 81. * mv_xsi[7][0] * mv_xsi[7][0] - 18. * mv_xsi[7][0] * mv_xsi[7][1] - 18. * mv_xsi[7][0] - 27. * mv_xsi[7][1] - 27.),
		0.03125 * (27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] - 27. * mv_xsi[7][0] + 9.) },
	  { 0.03125 * (-81. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] - 81. * mv_xsi[7][0] * mv_xsi[7][0] - 18. * mv_xsi[7][0] * mv_xsi[7][1] - 18. * mv_xsi[7][0] + 27. * mv_xsi[7][1] + 27.),
		0.03125 * (-27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] - 9. * mv_xsi[7][0] * mv_xsi[7][0] + 27. * mv_xsi[7][0] + 9.) },
	  { 0.03125 * (27. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][1] + 27. * mv_xsi[7][0] * mv_xsi[7][0] + 18. * mv_xsi[7][0] * mv_xsi[7][1] + 18. * mv_xsi[7][0] - mv_xsi[7][1] - 1.),
		0.03125 * (9. * mv_xsi[7][0] * mv_xsi[7][0] * mv_xsi[7][0] + 9. * mv_xsi[7][0] * mv_xsi[7][0] - mv_xsi[7][0] - 1.) } } };
