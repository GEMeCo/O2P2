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
// Linear fiber element, with 2 to 4 nodes.
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
			  * @class Elem_Lin
			  *
			  * @brief Linear element with 2 to 4 nodes.
			  * @details Linear bar element, with Lagrangian interpolation functions
			  * Options for integration points: 1, 2 and 3.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  * @tparam nNodes Number of nodes. Must be: 2, 3 or 4.
			  * @tparam nIP Number of integration points. Must be: 2, 3 or 4.
			  */
			template<int nDim, int nNodes, int nIP>
			class Elem_Lin : public ElemLinear<nDim>
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Elem_Lin() = delete;

			public:
				/** Constructor for linear bar elements.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				explicit Elem_Lin(std::shared_ptr<O2P2::Geom::Material>& Material, std::shared_ptr<O2P2::Geom::CrossSection>& Section)
					: ElemLinear<nDim>(Material, Section) { }

				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "1 " << nNodes - 1 << " ";
					for (int i = 0; i < nNodes; ++i) {
						msg << this->mv_Conect[i]->mv_index + add << " ";
					}
					msg << this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;
					msg << "1 " << nNodes - 1 << " ";
					for (int i = 0; i < nNodes; ++i) {
						msg << i + 1 + add << " ";
					}
					msg << this->mv_Mat->mv_index << "\n";
					return msg.str();
				}

				// Evaluates shape function in the point.
				std::vector<double> getShapeFcOnPoint(const double* Point) override;

				// Evaluates the derivative of shape function in the point.
				std::vector<double> getShapeDerivOnPoint(const double* Point) override;

				// Return a vector with values on the integration points currently known in the element' nodes.
				std::vector<double> getValueOnIPs(const double* value) override;

				// Returns a pointer to the first element of the shape functions (with size [nIP][nNodes]).
				double const* getShapeFc() const override { return &mv_Psi[0][0]; }

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][nNodes][nDim]).
				double const* getShapeDerivative() const override { return &mv_DPsi[0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return mv_weight; }

				// Returns the number of nodes of current element.
				const int getNumNodes() const override { return nNodes; }

				// Returns the number of faces of current element.
				const int getNumFaces() const override { return 1; }

				// Returns the number of integration points of current element.
				const int getNumIP() const override { return nIP; }

			private:
				// Evaluate approximated initial length. It is going to be stored in mv_Raidus. Must be called after setting the conectivity.
				void setGeomProperties() override;

			private:
				/** @brief Weights for numerical integration */
				static const double* mv_weight;

				/** @brief Shape functions */
				static const double mv_Psi[nIP][nNodes];

				/** @brief Shape functions derivative */
				static const double mv_DPsi[nIP][nNodes];
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
template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 2, 2>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(2);

	mi_Psi.at(0) = 0.5 - 0.5 * Point[0];
	mi_Psi.at(1) = 0.5 + 0.5 * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 2, 3>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(2);

	mi_Psi.at(0) = 0.5 - 0.5 * Point[0];
	mi_Psi.at(1) = 0.5 + 0.5 * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 2, 4>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(2);

	mi_Psi.at(0) = 0.5 - 0.5 * Point[0];
	mi_Psi.at(1) = 0.5 + 0.5 * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 2, 2>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(2);

	mi_Psi.at(0) = 0.5 - 0.5 * Point[0];
	mi_Psi.at(1) = 0.5 + 0.5 * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 2, 3>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(2);

	mi_Psi.at(0) = 0.5 - 0.5 * Point[0];
	mi_Psi.at(1) = 0.5 + 0.5 * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 2, 4>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(2);

	mi_Psi.at(0) = 0.5 - 0.5 * Point[0];
	mi_Psi.at(1) = 0.5 + 0.5 * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 3, 2>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(3);

	mi_Psi.at(0) = 0.5 * (Point[0] - 1.) * Point[0];
	mi_Psi.at(1) = (1. + Point[0]) * (1. - Point[0]);
	mi_Psi.at(2) = 0.5 * (Point[0] + 1.) * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 3, 3>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(3);

	mi_Psi.at(0) = 0.5 * (Point[0] - 1.) * Point[0];
	mi_Psi.at(1) = (1. + Point[0]) * (1. - Point[0]);
	mi_Psi.at(2) = 0.5 * (Point[0] + 1.) * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 3, 4>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(3);

	mi_Psi.at(0) = 0.5 * (Point[0] - 1.) * Point[0];
	mi_Psi.at(1) = (1. + Point[0]) * (1. - Point[0]);
	mi_Psi.at(2) = 0.5 * (Point[0] + 1.) * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 3, 2>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(3);

	mi_Psi.at(0) = 0.5 * (Point[0] - 1.) * Point[0];
	mi_Psi.at(1) = (1. + Point[0]) * (1. - Point[0]);
	mi_Psi.at(2) = 0.5 * (Point[0] + 1.) * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 3, 3>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(3);

	mi_Psi.at(0) = 0.5 * (Point[0] - 1.) * Point[0];
	mi_Psi.at(1) = (1. + Point[0]) * (1. - Point[0]);
	mi_Psi.at(2) = 0.5 * (Point[0] + 1.) * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 3, 4>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(3);

	mi_Psi.at(0) = 0.5 * (Point[0] - 1.) * Point[0];
	mi_Psi.at(1) = (1. + Point[0]) * (1. - Point[0]);
	mi_Psi.at(2) = 0.5 * (Point[0] + 1.) * Point[0];

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 4, 2>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = -0.0625 + Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(1) = 0.5625 - (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 + (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(2) = 0.5625 + (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 - (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(3) = -0.0625 - Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 + (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 4, 3>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = -0.0625 + Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(1) = 0.5625 - (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 + (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(2) = 0.5625 + (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 - (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(3) = -0.0625 - Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 + (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 4, 4>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = -0.0625 + Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(1) = 0.5625 - (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 + (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(2) = 0.5625 + (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 - (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(3) = -0.0625 - Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 + (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 4, 2>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = -0.0625 + Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(1) = 0.5625 - (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 + (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(2) = 0.5625 + (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 - (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(3) = -0.0625 - Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 + (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 4, 3>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = -0.0625 + Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(1) = 0.5625 - (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 + (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(2) = 0.5625 + (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 - (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(3) = -0.0625 - Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 + (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;

	return mi_Psi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 4, 4>::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(4);

	mi_Psi.at(0) = -0.0625 + Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(1) = 0.5625 - (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 + (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(2) = 0.5625 + (27.0 * Point[0]) / 16.0 - (9.0 * Point[0] * Point[0]) / 16.0 - (27.0 * Point[0] * Point[0] * Point[0]) / 16.0;
	mi_Psi.at(3) = -0.0625 - Point[0] / 16.0 + (9.0 * Point[0] * Point[0]) / 16.0 + (9.0 * Point[0] * Point[0] * Point[0]) / 16.0;

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 2, 2>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(2 * 1);

	mi_DPsi.at(0) = -0.5;
	mi_DPsi.at(1) = 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 2, 3>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(2 * 1);

	mi_DPsi.at(0) = -0.5;
	mi_DPsi.at(1) = 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 2, 4>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(2 * 1);

	mi_DPsi.at(0) = -0.5;
	mi_DPsi.at(1) = 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 2, 2>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(2 * 1);

	mi_DPsi.at(0) = -0.5;
	mi_DPsi.at(1) = 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 2, 3>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(2 * 1);

	mi_DPsi.at(0) = -0.5;
	mi_DPsi.at(1) = 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 2, 4>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(2 * 1);

	mi_DPsi.at(0) = -0.5;
	mi_DPsi.at(1) = 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 3, 2>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(3 * 1);

	mi_DPsi.at(0) = Point[0] - 0.5;
	mi_DPsi.at(1) = -2. * Point[0];
	mi_DPsi.at(2) = Point[0] + 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 3, 3>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(3 * 1);

	mi_DPsi.at(0) = Point[0] - 0.5;
	mi_DPsi.at(1) = -2. * Point[0];
	mi_DPsi.at(2) = Point[0] + 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 3, 4>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(3 * 1);

	mi_DPsi.at(0) = Point[0] - 0.5;
	mi_DPsi.at(1) = -2. * Point[0];
	mi_DPsi.at(2) = Point[0] + 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 3, 2>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(3 * 1);

	mi_DPsi.at(0) = Point[0] - 0.5;
	mi_DPsi.at(1) = -2. * Point[0];
	mi_DPsi.at(2) = Point[0] + 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 3, 3>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(3 * 1);

	mi_DPsi.at(0) = Point[0] - 0.5;
	mi_DPsi.at(1) = -2. * Point[0];
	mi_DPsi.at(2) = Point[0] + 0.5;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 3, 4>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(3 * 1);

	mi_DPsi.at(0) = Point[0] - 0.5;
	mi_DPsi.at(1) = -2. * Point[0];
	mi_DPsi.at(2) = Point[0] + 0.5;

	return mi_DPsi;
};


template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 4, 2>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 1);

	mi_DPsi.at(0) = 0.0625 + (9.0 * Point[0]) / 8.0 - (27.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(1) = -1.6875 - (9.0 * Point[0]) / 8.0 + (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(2) = 1.6875 - (9.0 * Point[0]) / 8.0 - (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(3) = -0.0625 + (9.0 * Point[0]) / 8.0 + (27.0 * Point[0] * Point[0]) / 16.0;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 4, 3>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 1);

	mi_DPsi.at(0) = 0.0625 + (9.0 * Point[0]) / 8.0 - (27.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(1) = -1.6875 - (9.0 * Point[0]) / 8.0 + (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(2) = 1.6875 - (9.0 * Point[0]) / 8.0 - (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(3) = -0.0625 + (9.0 * Point[0]) / 8.0 + (27.0 * Point[0] * Point[0]) / 16.0;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<2, 4, 4>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 1);

	mi_DPsi.at(0) = 0.0625 + (9.0 * Point[0]) / 8.0 - (27.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(1) = -1.6875 - (9.0 * Point[0]) / 8.0 + (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(2) = 1.6875 - (9.0 * Point[0]) / 8.0 - (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(3) = -0.0625 + (9.0 * Point[0]) / 8.0 + (27.0 * Point[0] * Point[0]) / 16.0;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 4, 2>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 1);

	mi_DPsi.at(0) = 0.0625 + (9.0 * Point[0]) / 8.0 - (27.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(1) = -1.6875 - (9.0 * Point[0]) / 8.0 + (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(2) = 1.6875 - (9.0 * Point[0]) / 8.0 - (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(3) = -0.0625 + (9.0 * Point[0]) / 8.0 + (27.0 * Point[0] * Point[0]) / 16.0;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 4, 3>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 1);

	mi_DPsi.at(0) = 0.0625 + (9.0 * Point[0]) / 8.0 - (27.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(1) = -1.6875 - (9.0 * Point[0]) / 8.0 + (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(2) = 1.6875 - (9.0 * Point[0]) / 8.0 - (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(3) = -0.0625 + (9.0 * Point[0]) / 8.0 + (27.0 * Point[0] * Point[0]) / 16.0;

	return mi_DPsi;
};

template<> inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<3, 4, 4>::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(4 * 1);

	mi_DPsi.at(0) = 0.0625 + (9.0 * Point[0]) / 8.0 - (27.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(1) = -1.6875 - (9.0 * Point[0]) / 8.0 + (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(2) = 1.6875 - (9.0 * Point[0]) / 8.0 - (81.0 * Point[0] * Point[0]) / 16.0;
	mi_DPsi.at(3) = -0.0625 + (9.0 * Point[0]) / 8.0 + (27.0 * Point[0] * Point[0]) / 16.0;

	return mi_DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
template<int nDim, int nNodes, int nIP>
inline void O2P2::Geom::Elem::Elem_Lin<nDim, nNodes, nIP>::setGeomProperties() {

	// Allocate an array with size mv_ElDim to which mv_Centroid points to.
	this->mv_Centroid = std::make_unique<double[]>(nDim);

	// Memory requested by make_unique is not empty
	for (int i = 0; i < nDim; ++i) this->mv_Centroid[i] = 0.;

	for (int i = 0; i < nNodes; ++i) {
		for (int j = 0; j < nDim; ++j) {
			this->mv_Centroid[j] += this->mv_Conect[i]->getInitPos()[j];
		}
	}

	// Finishing up
	for (int i = 0; i < nDim; i++) this->mv_Centroid[i] /= nNodes;

	// Approximated length (sum of the distance between nodes, thus sum of linear segments)
	this->mv_Radius = 0.;
	for (int i = 1; i < nNodes; ++i) {
		for (int j = 0; j < nDim; ++j) {
			this->mv_Radius += (this->mv_Conect[i]->getInitPos()[j] - this->mv_Conect[i - 1]->getInitPos()[j]) * (this->mv_Conect[i]->getInitPos()[j] - this->mv_Conect[i - 1]->getInitPos()[j]);
		}
	}
	this->mv_Radius = std::sqrt(this->mv_Radius);
};


// ================================================================================================
//
// Implementation of Member Function: getValueOnIPs
// Return the values on the integration points currently known in the element' nodes
// 
// ================================================================================================
template<int nDim, int nNodes, int nIP>
inline std::vector<double> O2P2::Geom::Elem::Elem_Lin<nDim, nNodes, nIP>::getValueOnIPs(const double* value) {

	// return value
	std::vector<double> mi_valueOnIp(nIP, 0.);

	for (int i = 0; i < nIP; i++) {
		for (int j = 0; j < nNodes; j++) {
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
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 2, 2>::mv_weight = &Gauss1D::Wg_2P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 3, 2>::mv_weight = &Gauss1D::Wg_2P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 4, 2>::mv_weight = &Gauss1D::Wg_2P[0];

template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 2, 3>::mv_weight = &Gauss1D::Wg_3P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 3, 3>::mv_weight = &Gauss1D::Wg_3P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 4, 3>::mv_weight = &Gauss1D::Wg_3P[0];

template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 2, 4>::mv_weight = &Gauss1D::Wg_4P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 3, 4>::mv_weight = &Gauss1D::Wg_4P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<2, 4, 4>::mv_weight = &Gauss1D::Wg_4P[0];

template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 2, 2>::mv_weight = &Gauss1D::Wg_2P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 3, 2>::mv_weight = &Gauss1D::Wg_2P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 4, 2>::mv_weight = &Gauss1D::Wg_2P[0];

template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 2, 3>::mv_weight = &Gauss1D::Wg_3P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 3, 3>::mv_weight = &Gauss1D::Wg_3P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 4, 3>::mv_weight = &Gauss1D::Wg_3P[0];

template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 2, 4>::mv_weight = &Gauss1D::Wg_4P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 3, 4>::mv_weight = &Gauss1D::Wg_4P[0];
template<> const double* O2P2::Geom::Elem::Elem_Lin<3, 4, 4>::mv_weight = &Gauss1D::Wg_4P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double O2P2::Geom::Elem::Elem_Lin<2, 2, 2>::mv_Psi[2][2] = {
	{ 0.5 - 0.5 * Gauss1D::Qsi_2P[0], 0.5 + 0.5 * Gauss1D::Qsi_2P[0] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_2P[1], 0.5 + 0.5 * Gauss1D::Qsi_2P[1] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 2, 3>::mv_Psi[3][2] = {
	{ 0.5 - 0.5 * Gauss1D::Qsi_3P[0], 0.5 + 0.5 * Gauss1D::Qsi_3P[0] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_3P[1], 0.5 + 0.5 * Gauss1D::Qsi_3P[1] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_3P[2], 0.5 + 0.5 * Gauss1D::Qsi_3P[2] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 2, 4>::mv_Psi[4][2] = {
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[0], 0.5 + 0.5 * Gauss1D::Qsi_4P[0] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[1], 0.5 + 0.5 * Gauss1D::Qsi_4P[1] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[2], 0.5 + 0.5 * Gauss1D::Qsi_4P[2] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[3], 0.5 + 0.5 * Gauss1D::Qsi_4P[3] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 2, 2>::mv_Psi[2][2] = {
	{ 0.5 - 0.5 * Gauss1D::Qsi_2P[0], 0.5 + 0.5 * Gauss1D::Qsi_2P[0] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_2P[1], 0.5 + 0.5 * Gauss1D::Qsi_2P[1] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 2, 3>::mv_Psi[3][2] = {
	{ 0.5 - 0.5 * Gauss1D::Qsi_3P[0], 0.5 + 0.5 * Gauss1D::Qsi_3P[0] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_3P[1], 0.5 + 0.5 * Gauss1D::Qsi_3P[1] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_3P[2], 0.5 + 0.5 * Gauss1D::Qsi_3P[2] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 2, 4>::mv_Psi[4][2] = {
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[0], 0.5 + 0.5 * Gauss1D::Qsi_4P[0] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[1], 0.5 + 0.5 * Gauss1D::Qsi_4P[1] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[2], 0.5 + 0.5 * Gauss1D::Qsi_4P[2] } ,
	{ 0.5 - 0.5 * Gauss1D::Qsi_4P[3], 0.5 + 0.5 * Gauss1D::Qsi_4P[3] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 3, 2>::mv_Psi[2][3] = {
	{ 0.5 * (Gauss1D::Qsi_2P[0] - 1.) * Gauss1D::Qsi_2P[0],  (1. + Gauss1D::Qsi_2P[0]) * (1. - Gauss1D::Qsi_2P[0]), 0.5 * (Gauss1D::Qsi_2P[0] + 1.) * Gauss1D::Qsi_2P[0] },
	{ 0.5 * (Gauss1D::Qsi_2P[1] - 1.) * Gauss1D::Qsi_2P[1],  (1. + Gauss1D::Qsi_2P[1]) * (1. - Gauss1D::Qsi_2P[1]), 0.5 * (Gauss1D::Qsi_2P[1] + 1.) * Gauss1D::Qsi_2P[1] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 3, 3>::mv_Psi[3][3] = {
	{ 0.5 * (Gauss1D::Qsi_3P[0] - 1.) * Gauss1D::Qsi_3P[0],  (1. + Gauss1D::Qsi_3P[0]) * (1. - Gauss1D::Qsi_3P[0]), 0.5 * (Gauss1D::Qsi_3P[0] + 1.) * Gauss1D::Qsi_3P[0] },
	{ 0.5 * (Gauss1D::Qsi_3P[1] - 1.) * Gauss1D::Qsi_3P[1],  (1. + Gauss1D::Qsi_3P[1]) * (1. - Gauss1D::Qsi_3P[1]), 0.5 * (Gauss1D::Qsi_3P[1] + 1.) * Gauss1D::Qsi_3P[1] },
	{ 0.5 * (Gauss1D::Qsi_3P[2] - 1.) * Gauss1D::Qsi_3P[2],  (1. + Gauss1D::Qsi_3P[2]) * (1. - Gauss1D::Qsi_3P[2]), 0.5 * (Gauss1D::Qsi_3P[2] + 1.) * Gauss1D::Qsi_3P[2] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 3, 4>::mv_Psi[4][3] = {
	{ 0.5 * (Gauss1D::Qsi_4P[0] - 1.) * Gauss1D::Qsi_4P[0],  (1. + Gauss1D::Qsi_4P[0]) * (1. - Gauss1D::Qsi_4P[0]), 0.5 * (Gauss1D::Qsi_4P[0] + 1.) * Gauss1D::Qsi_4P[0] },
	{ 0.5 * (Gauss1D::Qsi_4P[1] - 1.) * Gauss1D::Qsi_4P[1],  (1. + Gauss1D::Qsi_4P[1]) * (1. - Gauss1D::Qsi_4P[1]), 0.5 * (Gauss1D::Qsi_4P[1] + 1.) * Gauss1D::Qsi_4P[1] },
	{ 0.5 * (Gauss1D::Qsi_4P[2] - 1.) * Gauss1D::Qsi_4P[2],  (1. + Gauss1D::Qsi_4P[2]) * (1. - Gauss1D::Qsi_4P[2]), 0.5 * (Gauss1D::Qsi_4P[2] + 1.) * Gauss1D::Qsi_4P[2] },
	{ 0.5 * (Gauss1D::Qsi_4P[3] - 1.) * Gauss1D::Qsi_4P[3],  (1. + Gauss1D::Qsi_4P[3]) * (1. - Gauss1D::Qsi_4P[3]), 0.5 * (Gauss1D::Qsi_4P[3] + 1.) * Gauss1D::Qsi_4P[3] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 3, 2>::mv_Psi[2][3] = {
	{ 0.5 * (Gauss1D::Qsi_2P[0] - 1.) * Gauss1D::Qsi_2P[0],  (1. + Gauss1D::Qsi_2P[0]) * (1. - Gauss1D::Qsi_2P[0]), 0.5 * (Gauss1D::Qsi_2P[0] + 1.) * Gauss1D::Qsi_2P[0] },
	{ 0.5 * (Gauss1D::Qsi_2P[1] - 1.) * Gauss1D::Qsi_2P[1],  (1. + Gauss1D::Qsi_2P[1]) * (1. - Gauss1D::Qsi_2P[1]), 0.5 * (Gauss1D::Qsi_2P[1] + 1.) * Gauss1D::Qsi_2P[1] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 3, 3>::mv_Psi[3][3] = {
	{ 0.5 * (Gauss1D::Qsi_3P[0] - 1.) * Gauss1D::Qsi_3P[0],  (1. + Gauss1D::Qsi_3P[0]) * (1. - Gauss1D::Qsi_3P[0]), 0.5 * (Gauss1D::Qsi_3P[0] + 1.) * Gauss1D::Qsi_3P[0] },
	{ 0.5 * (Gauss1D::Qsi_3P[1] - 1.) * Gauss1D::Qsi_3P[1],  (1. + Gauss1D::Qsi_3P[1]) * (1. - Gauss1D::Qsi_3P[1]), 0.5 * (Gauss1D::Qsi_3P[1] + 1.) * Gauss1D::Qsi_3P[1] },
	{ 0.5 * (Gauss1D::Qsi_3P[2] - 1.) * Gauss1D::Qsi_3P[2],  (1. + Gauss1D::Qsi_3P[2]) * (1. - Gauss1D::Qsi_3P[2]), 0.5 * (Gauss1D::Qsi_3P[2] + 1.) * Gauss1D::Qsi_3P[2] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 3, 4>::mv_Psi[4][3] = {
	{ 0.5 * (Gauss1D::Qsi_4P[0] - 1.) * Gauss1D::Qsi_4P[0],  (1. + Gauss1D::Qsi_4P[0]) * (1. - Gauss1D::Qsi_4P[0]), 0.5 * (Gauss1D::Qsi_4P[0] + 1.) * Gauss1D::Qsi_4P[0] },
	{ 0.5 * (Gauss1D::Qsi_4P[1] - 1.) * Gauss1D::Qsi_4P[1],  (1. + Gauss1D::Qsi_4P[1]) * (1. - Gauss1D::Qsi_4P[1]), 0.5 * (Gauss1D::Qsi_4P[1] + 1.) * Gauss1D::Qsi_4P[1] },
	{ 0.5 * (Gauss1D::Qsi_4P[2] - 1.) * Gauss1D::Qsi_4P[2],  (1. + Gauss1D::Qsi_4P[2]) * (1. - Gauss1D::Qsi_4P[2]), 0.5 * (Gauss1D::Qsi_4P[2] + 1.) * Gauss1D::Qsi_4P[2] },
	{ 0.5 * (Gauss1D::Qsi_4P[3] - 1.) * Gauss1D::Qsi_4P[3],  (1. + Gauss1D::Qsi_4P[3]) * (1. - Gauss1D::Qsi_4P[3]), 0.5 * (Gauss1D::Qsi_4P[3] + 1.) * Gauss1D::Qsi_4P[3] } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 4, 2>::mv_Psi[2][4] = {
	{ -0.0625 + Gauss1D::Qsi_2P[0] / 16.0 + (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_2P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 + (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_2P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 - (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_2P[0] / 16.0 + (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 + (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_2P[1] / 16.0 + (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_2P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 + (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_2P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 - (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_2P[1] / 16.0 + (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 + (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 4, 3>::mv_Psi[3][4] = {
	{ -0.0625 + Gauss1D::Qsi_3P[0] / 16.0 + (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_3P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 + (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_3P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 - (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_3P[0] / 16.0 + (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 + (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_3P[1] / 16.0 + (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_3P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 + (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_3P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 - (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_3P[1] / 16.0 + (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 + (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_3P[2] / 16.0 + (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_3P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 + (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_3P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 - (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_3P[2] / 16.0 + (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 + (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 4, 4>::mv_Psi[4][4] = {
	{ -0.0625 + Gauss1D::Qsi_4P[0] / 16.0 + (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[0] / 16.0 + (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_4P[1] / 16.0 + (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[1] / 16.0 + (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_4P[2] / 16.0 + (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[2] / 16.0 + (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_4P[3] / 16.0 + (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[3]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[3]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[3] / 16.0 + (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 4, 2>::mv_Psi[2][4] = {
	{ -0.0625 + Gauss1D::Qsi_2P[0] / 16.0 + (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_2P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 + (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_2P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 - (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_2P[0] / 16.0 + (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 + (9.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_2P[1] / 16.0 + (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_2P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 + (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_2P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 - (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_2P[1] / 16.0 + (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 + (9.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 4, 3>::mv_Psi[3][4] = {
	{ -0.0625 + Gauss1D::Qsi_3P[0] / 16.0 + (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_3P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 + (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_3P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 - (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_3P[0] / 16.0 + (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 + (9.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_3P[1] / 16.0 + (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_3P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 + (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_3P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 - (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_3P[1] / 16.0 + (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 + (9.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_3P[2] / 16.0 + (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_3P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 + (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_3P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 - (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_3P[2] / 16.0 + (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 + (9.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 4, 4>::mv_Psi[4][4] = {
	{ -0.0625 + Gauss1D::Qsi_4P[0] / 16.0 + (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[0]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[0] / 16.0 + (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_4P[1] / 16.0 + (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[1]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[1] / 16.0 + (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_4P[2] / 16.0 + (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[2]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[2] / 16.0 + (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 },

	{ -0.0625 + Gauss1D::Qsi_4P[3] / 16.0 + (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	   0.5625 - (27.0 * Gauss1D::Qsi_4P[3]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 + (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	   0.5625 + (27.0 * Gauss1D::Qsi_4P[3]) / 16.0 - (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 - (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	  -0.0625 - Gauss1D::Qsi_4P[3] / 16.0 + (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 + (9.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 } };


// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
template<> const double O2P2::Geom::Elem::Elem_Lin<2, 2, 2>::mv_DPsi[2][2] = {
	{ -0.5, 0.5 }, { -0.5, 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 2, 3>::mv_DPsi[3][2] = {
	{ -0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 2, 4>::mv_DPsi[4][2] = {
	{ -0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 2, 2>::mv_DPsi[2][2] = {
	{ -0.5, 0.5 }, { -0.5, 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 2, 3>::mv_DPsi[3][2] = {
	{ -0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 2, 4>::mv_DPsi[4][2] = {
	{ -0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, 0.5 }, { -0.5, 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 3, 2>::mv_DPsi[2][3] = {
	{ Gauss1D::Qsi_2P[0] - 0.5, -2. * Gauss1D::Qsi_2P[0], Gauss1D::Qsi_2P[0] + 0.5 },
	{ Gauss1D::Qsi_2P[1] - 0.5, -2. * Gauss1D::Qsi_2P[1], Gauss1D::Qsi_2P[1] + 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 3, 3>::mv_DPsi[3][3] = {
	{ Gauss1D::Qsi_3P[0] - 0.5, -2. * Gauss1D::Qsi_3P[0], Gauss1D::Qsi_3P[0] + 0.5 },
	{ Gauss1D::Qsi_3P[1] - 0.5, -2. * Gauss1D::Qsi_3P[1], Gauss1D::Qsi_3P[1] + 0.5 },
	{ Gauss1D::Qsi_3P[2] - 0.5, -2. * Gauss1D::Qsi_3P[2], Gauss1D::Qsi_3P[2] + 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 3, 4>::mv_DPsi[4][3] = {
	{ Gauss1D::Qsi_4P[0] - 0.5, -2. * Gauss1D::Qsi_4P[0], Gauss1D::Qsi_4P[0] + 0.5 },
	{ Gauss1D::Qsi_4P[1] - 0.5, -2. * Gauss1D::Qsi_4P[1], Gauss1D::Qsi_4P[1] + 0.5 },
	{ Gauss1D::Qsi_4P[2] - 0.5, -2. * Gauss1D::Qsi_4P[2], Gauss1D::Qsi_4P[2] + 0.5 },
	{ Gauss1D::Qsi_4P[3] - 0.5, -2. * Gauss1D::Qsi_4P[3], Gauss1D::Qsi_4P[3] + 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 3, 2>::mv_DPsi[2][3] = {
	{ Gauss1D::Qsi_2P[0] - 0.5, -2. * Gauss1D::Qsi_2P[0], Gauss1D::Qsi_2P[0] + 0.5 },
	{ Gauss1D::Qsi_2P[1] - 0.5, -2. * Gauss1D::Qsi_2P[1], Gauss1D::Qsi_2P[1] + 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 3, 3>::mv_DPsi[3][3] = {
	{ Gauss1D::Qsi_3P[0] - 0.5, -2. * Gauss1D::Qsi_3P[0], Gauss1D::Qsi_3P[0] + 0.5 },
	{ Gauss1D::Qsi_3P[1] - 0.5, -2. * Gauss1D::Qsi_3P[1], Gauss1D::Qsi_3P[1] + 0.5 },
	{ Gauss1D::Qsi_3P[2] - 0.5, -2. * Gauss1D::Qsi_3P[2], Gauss1D::Qsi_3P[2] + 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 3, 4>::mv_DPsi[4][3] = {
	{ Gauss1D::Qsi_4P[0] - 0.5, -2. * Gauss1D::Qsi_4P[0], Gauss1D::Qsi_4P[0] + 0.5 },
	{ Gauss1D::Qsi_4P[1] - 0.5, -2. * Gauss1D::Qsi_4P[1], Gauss1D::Qsi_4P[1] + 0.5 },
	{ Gauss1D::Qsi_4P[2] - 0.5, -2. * Gauss1D::Qsi_4P[2], Gauss1D::Qsi_4P[2] + 0.5 },
	{ Gauss1D::Qsi_4P[3] - 0.5, -2. * Gauss1D::Qsi_4P[3], Gauss1D::Qsi_4P[3] + 0.5 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 4, 2>::mv_DPsi[2][4] = {
	{ 0.0625 + (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 - (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 + (81.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 - (81.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 + (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 - (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 + (81.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 - (81.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 + (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 4, 3>::mv_DPsi[3][4] = {
	{ 0.0625 + (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 - (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 + (81.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 - (81.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 + (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 - (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 + (81.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 - (81.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 + (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 - (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 + (81.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 - (81.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 + (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<2, 4, 4>::mv_DPsi[4][4] = {
	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 4, 2>::mv_DPsi[2][4] = {
	{ 0.0625 + (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 - (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 + (81.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 - (81.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_2P[0]) / 8.0 + (27.0 * Gauss1D::Qsi_2P[0] * Gauss1D::Qsi_2P[0]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 - (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 + (81.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 - (81.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_2P[1]) / 8.0 + (27.0 * Gauss1D::Qsi_2P[1] * Gauss1D::Qsi_2P[1]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 4, 3>::mv_DPsi[3][4] = {
	{ 0.0625 + (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 - (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 + (81.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 - (81.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_3P[0]) / 8.0 + (27.0 * Gauss1D::Qsi_3P[0] * Gauss1D::Qsi_3P[0]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 - (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 + (81.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 - (81.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_3P[1]) / 8.0 + (27.0 * Gauss1D::Qsi_3P[1] * Gauss1D::Qsi_3P[1]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 - (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 + (81.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 - (81.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_3P[2]) / 8.0 + (27.0 * Gauss1D::Qsi_3P[2] * Gauss1D::Qsi_3P[2]) / 16.0 } };

template<> const double O2P2::Geom::Elem::Elem_Lin<3, 4, 4>::mv_DPsi[4][4] = {
	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[0]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[0] * Gauss1D::Qsi_4P[0]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[1]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[1] * Gauss1D::Qsi_4P[1]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[2]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[2] * Gauss1D::Qsi_4P[2]) / 16.0 },

	{ 0.0625 + (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 - (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	 -1.6875 - (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 + (81.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	  1.6875 - (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 - (81.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0,
	 -0.0625 + (9.0 * Gauss1D::Qsi_4P[3]) / 8.0 + (27.0 * Gauss1D::Qsi_4P[3] * Gauss1D::Qsi_4P[3]) / 16.0 } };
