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
// Solid element, with cubic / linear interpolation functions, prism shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

namespace O2P2 {
	namespace Prep {
		namespace Elem {
			/** @ingroup Elements
			  * @class Elem_Pri20
			  *
			  * @brief Prismatic cubic / linear element with 20 nodes.
			  * @details Solid element, with cubic interpolation functions in plane direction and linear function in the extrusion direction, prism shape.
			  * @image html Elem_Pri20.png height=300
			  */
			class Elem_Pri20 : public ElementSolid
			{
			private:
				Elem_Pri20() = delete;

			public:
				/** Constructor for prism cubic / linear elements.
				  * @param Material Pointer to Material class.
				  */
				explicit Elem_Pri20(std::shared_ptr<O2P2::Prep::Material>& Material)
					: ElementSolid(Material) { }

				// Output function for AcadView, based on element index.
				const std::string printByIndex_AV(const size_t add) const override {
					std::stringstream msg;

					msg << "2 3 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[1]->mv_index + add << " " << this->mv_Conect[2]->mv_index + add << " "
						<< this->mv_Conect[3]->mv_index + add << " " << this->mv_Conect[4]->mv_index + add << " " << this->mv_Conect[5]->mv_index + add << " "
						<< this->mv_Conect[6]->mv_index + add << " " << this->mv_Conect[7]->mv_index + add << " " << this->mv_Conect[8]->mv_index + add << " "
						<< this->mv_Conect[9]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "2 3 " << this->mv_Conect[10]->mv_index + add << " " << this->mv_Conect[11]->mv_index + add << " " << this->mv_Conect[12]->mv_index + add << " "
						<< this->mv_Conect[13]->mv_index + add << " " << this->mv_Conect[14]->mv_index + add << " " << this->mv_Conect[15]->mv_index + add << " "
						<< this->mv_Conect[16]->mv_index + add << " " << this->mv_Conect[17]->mv_index + add << " " << this->mv_Conect[18]->mv_index + add << " "
						<< this->mv_Conect[19]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[3]->mv_index + add << " " << this->mv_Conect[10]->mv_index + add << " " << this->mv_Conect[13]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << this->mv_Conect[3]->mv_index + add << " " << this->mv_Conect[9]->mv_index + add << " " << this->mv_Conect[13]->mv_index + add << " " << this->mv_Conect[19]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << this->mv_Conect[9]->mv_index + add << " " << this->mv_Conect[0]->mv_index + add << " " << this->mv_Conect[19]->mv_index + add << " " << this->mv_Conect[10]->mv_index + add << " " << this->mv_Mat->mv_index << "\n";

					return msg.str();
				}

				// Output function for AcadView, based on element node number.
				const std::string printByAdder_AV(const size_t add) const override {
					std::stringstream msg;

					msg << "2 3 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (4 + add) << " " << (5 + add) << " "
						<< (6 + add) << " " << (7 + add) << " " << (8 + add) << " " << (9 + add) << " " << (10 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "2 3 " << (11 + add) << " " << (12 + add) << " " << (13 + add) << " " << (14 + add) << " " << (15 + add) << " "
						<< (16 + add) << " " << (17 + add) << " " << (18 + add) << " " << (19 + add) << " " << (20 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << (1 + add) << " " << (4 + add) << " " << (11 + add) << " " << (14 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << (4 + add) << " " << (10 + add) << " " << (14 + add) << " " << (20 + add) << " " << this->mv_Mat->mv_index << "\n";
					msg << "3 1 " << (10 + add) << " " << (1 + add) << " " << (20 + add) << " " << (11 + add) << " " << this->mv_Mat->mv_index << "\n";

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

				// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][mv_numNodes][mv_Dim]).
				double const* getShapeDerivative() const override { return &mv_DPsi[0][0][0]; }

				// Returns a pointer to the weight of the integation points (with size [nIP]).
				double const* getWeight() const override { return &mv_weight[0]; }

				// Returns the number of nodes of current element.
				int getNumNodes() override { return mv_numNodes; }

				// Returns the number of faces of current element.
				int getNumFaces() override { return mv_numFaces; }

				// Returns the number of integration points of current element.
				int getNumIP() override { return mv_numIP; }

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates.
				  */
				inline bool evaluateXsi(const std::array<double, mv_Dim> xsi) override {

					std::array<double, mv_Dim + 1> new_xsi = {};

					for (int i = 0; i < mv_Dim - 1; ++i) {
						new_xsi.at(i) = xsi.at(i);
						new_xsi.at(mv_Dim - 1) -= xsi.at(i);
					}
					new_xsi.at(mv_Dim - 1) += 1.;
					new_xsi.at(mv_Dim) = xsi.at(mv_Dim - 1);

					if (*std::max_element(new_xsi.begin(), new_xsi.end()) < 1.000001 && *std::min_element(new_xsi.begin(), new_xsi.end() - 1) > -0.000001 && new_xsi.at(mv_Dim) > -1.000001) return true;
					return false;
				}

			private:
				// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
				void setGeomProperties() override;

			private:
				/** @brief Number of Nodes */
				static const int mv_numNodes{ 20 };

				/** @brief Number of Integration Points */
				static const int mv_numIP{ 14 };

				/** @brief Number of Faces */
				static const int mv_numFaces{ 5 };

				/** @brief Weights for numerical integration */
				static const double mv_weight[mv_numIP];

				/** @brief Shape functions */
				static const double mv_Psi[mv_numIP][mv_numNodes];

				/** @brief Shape functions derivative */
				static const double mv_DPsi[mv_numIP][mv_numNodes][mv_Dim];

				/** @brief Integration points */
				static const double m_xsi[mv_numIP][mv_Dim];
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
inline std::vector<double> O2P2::Prep::Elem::Elem_Pri20::getShapeFcOnPoint(const double* Point) {
	std::vector<double> mi_Psi(20);

	mi_Psi.at(0) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(1) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(2) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(3) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(4) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(5) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(6) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(7) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(8) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(9) = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_Psi.at(10) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(11) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(12) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(13) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(14) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(15) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(16) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(17) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(18) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_Psi.at(19) = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * (0.5 * Point[2] + 0.5);

	return mi_Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline std::vector<double> O2P2::Prep::Elem::Elem_Pri20::getShapeDerivOnPoint(const double* Point) {
	std::vector<double> mi_DPsi(20 * 3);

	mi_DPsi.at( 0) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 1) = (40.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 13.5 * Point[1] * Point[1] - 45. * Point[0] - 22.5 * Point[1] + 9.) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 2) = (-40.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] + 36. * Point[0] + 4.5 * Point[1] - 4.5) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 3) = (13.5 * Point[0] * Point[0] - 9. * Point[0] + 1.) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 4) = (27. * Point[0] * Point[1] + 27. * Point[1] * Point[1] - 22.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 5) = (-54. * Point[0] * Point[1] - 27. * Point[1] * Point[1] + 27. * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 6) = (27. * Point[0] * Point[1] - 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 7) = (-13.5 * Point[1] * Point[1] + 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 8) = (13.5 * Point[1] * Point[1] - 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at( 9) = 0.;
	mi_DPsi.at(10) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(11) = (40.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 13.5 * Point[1] * Point[1] - 45. * Point[0] - 22.5 * Point[1] + 9.) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(12) = (-40.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] + 36. * Point[0] + 4.5 * Point[1] - 4.5) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(13) = (13.5 * Point[0] * Point[0] - 9. * Point[0] + 1.) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(14) = (27. * Point[0] * Point[1] + 27. * Point[1] * Point[1] - 22.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(15) = (-54. * Point[0] * Point[1] - 27. * Point[1] * Point[1] + 27. * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(16) = (27. * Point[0] * Point[1] - 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(17) = (-13.5 * Point[1] * Point[1] + 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(18) = (13.5 * Point[1] * Point[1] - 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(19) = 0.;

	mi_DPsi.at(20) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(21) = (27. * Point[0] * Point[0] + 27. * Point[0] * Point[1] - 22.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(22) = (-13.5 * Point[0] * Point[0] + 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(23) = 0.;
	mi_DPsi.at(24) = (13.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 40.5 * Point[1] * Point[1] - 22.5 * Point[0] - 45. * Point[1] + 9.) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(25) = (-27. * Point[0] * Point[0] - 54. * Point[0] * Point[1] + 27. * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(26) = (13.5 * Point[0] * Point[0] - 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(27) = (-27. * Point[0] * Point[1] - 40.5 * Point[1] * Point[1] + 4.5 * Point[0] + 36. * Point[1] - 4.5) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(28) = (27. * Point[0] * Point[1] - 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(29)  = (13.5 * Point[1] * Point[1] - 9. * Point[1] + 1.) * (-0.5 * Point[2] + 0.5);
	mi_DPsi.at(30) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(31) = (27. * Point[0] * Point[0] + 27. * Point[0] * Point[1] - 22.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(32) = (-13.5 * Point[0] * Point[0] + 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(33) = 0.;
	mi_DPsi.at(34) = (13.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 40.5 * Point[1] * Point[1] - 22.5 * Point[0] - 45. * Point[1] + 9.) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(35) = (-27. * Point[0] * Point[0] - 54. * Point[0] * Point[1] + 27. * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(36) = (13.5 * Point[0] * Point[0] - 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(37) = (-27. * Point[0] * Point[1] - 40.5 * Point[1] * Point[1] + 4.5 * Point[0] + 36. * Point[1] - 4.5) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(38) = (27. * Point[0] * Point[1] - 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	mi_DPsi.at(39) = (13.5 * Point[1] * Point[1] - 9. * Point[1] + 1.) * (0.5 * Point[2] + 0.5);

	mi_DPsi.at(40) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * -0.5;
	mi_DPsi.at(41) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * -0.5;
	mi_DPsi.at(42) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * -0.5;
	mi_DPsi.at(43) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * -0.5;
	mi_DPsi.at(44) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * -0.5;
	mi_DPsi.at(45) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * -0.5;
	mi_DPsi.at(46) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * -0.5;
	mi_DPsi.at(47) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * -0.5;
	mi_DPsi.at(48) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * -0.5;
	mi_DPsi.at(49)  = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * -0.5;
	mi_DPsi.at(50) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * 0.5;
	mi_DPsi.at(51) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * 0.5;
	mi_DPsi.at(52) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * 0.5;
	mi_DPsi.at(53) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * 0.5;
	mi_DPsi.at(54) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * 0.5;
	mi_DPsi.at(55) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * 0.5;
	mi_DPsi.at(56) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * 0.5;
	mi_DPsi.at(57) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * 0.5;
	mi_DPsi.at(58) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * 0.5;
	mi_DPsi.at(59) = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * 0.5;

	return mi_DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void O2P2::Prep::Elem::Elem_Pri20::setGeomProperties() {

	const int nVertices = 6;

	// Allocate an array with size mv_Dim to which mv_Centroid points to.
	mv_Centroid = std::make_unique<double[]>(mv_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<O2P2::Prep::Node<mv_Dim>*, nVertices> vertices;
	vertices[0] = mv_Conect[0].get();
	vertices[1] = mv_Conect[3].get();
	vertices[2] = mv_Conect[9].get();
	vertices[3] = mv_Conect[10].get();
	vertices[4] = mv_Conect[13].get();
	vertices[5] = mv_Conect[19].get();

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
inline std::vector<double> O2P2::Prep::Elem::Elem_Pri20::getValueOnIPs(const double* value) {

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
inline const double O2P2::Prep::Elem::Elem_Pri20::m_xsi[mv_numIP][mv_Dim] =
	{ { 0.333333333333333, 0.333333333333333, -0.5773502691896257645091488 },
	  { 0.797426985353087, 0.101286507323456, -0.5773502691896257645091488 },
	  { 0.101286507323456, 0.797426985353087, -0.5773502691896257645091488 },
	  { 0.101286507323456, 0.101286507323456, -0.5773502691896257645091488 },
	  { 0.470142064105115, 0.470142064105115, -0.5773502691896257645091488 },
	  { 0.059715871789770, 0.470142064105115, -0.5773502691896257645091488 },
	  { 0.470142064105115, 0.059715871789770, -0.5773502691896257645091488 },
	  { 0.333333333333333, 0.333333333333333,  0.5773502691896257645091488 },
	  { 0.797426985353087, 0.101286507323456,  0.5773502691896257645091488 },
	  { 0.101286507323456, 0.797426985353087,  0.5773502691896257645091488 },
	  { 0.101286507323456, 0.101286507323456,  0.5773502691896257645091488 },
	  { 0.470142064105115, 0.470142064105115,  0.5773502691896257645091488 },
	  { 0.059715871789770, 0.470142064105115,  0.5773502691896257645091488 },
	  { 0.470142064105115, 0.059715871789770,  0.5773502691896257645091488 } };

// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri20::mv_weight[mv_numIP] = { 0.1125, 0.0629695902724135, 0.0629695902724135, 0.0629695902724135, 0.066197076394253, 0.066197076394253, 0.066197076394253, 0.1125, 0.0629695902724135, 0.0629695902724135, 0.0629695902724135, 0.066197076394253, 0.066197076394253, 0.066197076394253 };

// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri20::mv_Psi[mv_numIP][mv_numNodes] = {
	{ (-4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][1] * m_xsi[0][1] - 5.5 * m_xsi[0][0] - 5.5 * m_xsi[0][1] + 1.) * (-0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][0] - 22.5 * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (-13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] * m_xsi[0][0] + 4.5 * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 4.5 * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][1] - 22.5 * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (-13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] + m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
	  (-4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][1] * m_xsi[0][1] - 5.5 * m_xsi[0][0] - 5.5 * m_xsi[0][1] + 1.) * (0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][0] - 22.5 * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
	  (-13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] * m_xsi[0][0] + 4.5 * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
	  (4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 4.5 * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][1] - 22.5 * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
	  (-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
	  (-13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
	  (13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
	  (4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] + m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5) },

	{ (-4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][1] * m_xsi[1][1] - 5.5 * m_xsi[1][0] - 5.5 * m_xsi[1][1] + 1.) * (-0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][0] - 22.5 * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (-13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] * m_xsi[1][0] + 4.5 * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 4.5 * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][1] - 22.5 * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (-13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] + m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
	  (-4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][1] * m_xsi[1][1] - 5.5 * m_xsi[1][0] - 5.5 * m_xsi[1][1] + 1.) * (0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][0] - 22.5 * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
	  (-13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] * m_xsi[1][0] + 4.5 * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
	  (4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 4.5 * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][1] - 22.5 * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
	  (-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
	  (-13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
	  (13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
	  (4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] + m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5) },

	{ (-4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][1] * m_xsi[2][1] - 5.5 * m_xsi[2][0] - 5.5 * m_xsi[2][1] + 1.) * (-0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][0] - 22.5 * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (-13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] * m_xsi[2][0] + 4.5 * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 4.5 * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][1] - 22.5 * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (-13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] + m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
	  (-4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][1] * m_xsi[2][1] - 5.5 * m_xsi[2][0] - 5.5 * m_xsi[2][1] + 1.) * (0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][0] - 22.5 * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
	  (-13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] * m_xsi[2][0] + 4.5 * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
	  (4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 4.5 * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][1] - 22.5 * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
	  (-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
	  (-13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
	  (13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
	  (4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] + m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5) },

	{ (-4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][1] * m_xsi[3][1] - 5.5 * m_xsi[3][0] - 5.5 * m_xsi[3][1] + 1.) * (-0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][0] - 22.5 * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (-13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] * m_xsi[3][0] + 4.5 * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 4.5 * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][1] - 22.5 * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (-13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] + m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
	  (-4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][1] * m_xsi[3][1] - 5.5 * m_xsi[3][0] - 5.5 * m_xsi[3][1] + 1.) * (0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][0] - 22.5 * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
	  (-13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] * m_xsi[3][0] + 4.5 * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
	  (4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 4.5 * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][1] - 22.5 * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
	  (-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
	  (-13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
	  (13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
	  (4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] + m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5) },

	{ (-4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][1] * m_xsi[4][1] - 5.5 * m_xsi[4][0] - 5.5 * m_xsi[4][1] + 1.) * (-0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][0] - 22.5 * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (-13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] * m_xsi[4][0] + 4.5 * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 4.5 * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][1] - 22.5 * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (-13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] + m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
	  (-4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][1] * m_xsi[4][1] - 5.5 * m_xsi[4][0] - 5.5 * m_xsi[4][1] + 1.) * (0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][0] - 22.5 * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
	  (-13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] * m_xsi[4][0] + 4.5 * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
	  (4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 4.5 * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][1] - 22.5 * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
	  (-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
	  (-13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
	  (13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
	  (4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] + m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5) },

	{ (-4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][1] * m_xsi[5][1] - 5.5 * m_xsi[5][0] - 5.5 * m_xsi[5][1] + 1.) * (-0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][0] - 22.5 * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (-13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] * m_xsi[5][0] + 4.5 * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 4.5 * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][1] - 22.5 * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (-13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] + m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
	  (-4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][1] * m_xsi[5][1] - 5.5 * m_xsi[5][0] - 5.5 * m_xsi[5][1] + 1.) * (0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][0] - 22.5 * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
	  (-13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] * m_xsi[5][0] + 4.5 * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
	  (4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 4.5 * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][1] - 22.5 * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
	  (-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
	  (-13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
	  (13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
	  (4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] + m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5) },

	{ (-4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][1] * m_xsi[6][1] - 5.5 * m_xsi[6][0] - 5.5 * m_xsi[6][1] + 1.) * (-0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][0] - 22.5 * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (-13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] * m_xsi[6][0] + 4.5 * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 4.5 * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][1] - 22.5 * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (-13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] + m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
	  (-4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][1] * m_xsi[6][1] - 5.5 * m_xsi[6][0] - 5.5 * m_xsi[6][1] + 1.) * (0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][0] - 22.5 * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
	  (-13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] * m_xsi[6][0] + 4.5 * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
	  (4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 4.5 * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][1] - 22.5 * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
	  (-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
	  (-13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
	  (13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
	  (4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] + m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5) },

	{ (-4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][1] * m_xsi[7][1] - 5.5 * m_xsi[7][0] - 5.5 * m_xsi[7][1] + 1.) * (-0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][0] - 22.5 * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (-13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] * m_xsi[7][0] + 4.5 * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 4.5 * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][1] - 22.5 * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (-13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] + m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
	  (-4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][1] * m_xsi[7][1] - 5.5 * m_xsi[7][0] - 5.5 * m_xsi[7][1] + 1.) * (0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][0] - 22.5 * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
	  (-13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] * m_xsi[7][0] + 4.5 * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
	  (4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 4.5 * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][1] - 22.5 * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
	  (-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
	  (-13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
	  (13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
	  (4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] + m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5) },

	{ (-4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][0] * m_xsi[8][0] + 18. * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][1] * m_xsi[8][1] - 5.5 * m_xsi[8][0] - 5.5 * m_xsi[8][1] + 1.) * (-0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] + 27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][0] - 22.5 * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (-13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][0] * m_xsi[8][0] + 4.5 * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 4.5 * m_xsi[8][0] * m_xsi[8][0] + m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][1] - 22.5 * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (-27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (-13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] + m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
	  (-4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][0] * m_xsi[8][0] + 18. * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][1] * m_xsi[8][1] - 5.5 * m_xsi[8][0] - 5.5 * m_xsi[8][1] + 1.) * (0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] + 27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][0] - 22.5 * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
	  (-13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][0] * m_xsi[8][0] + 4.5 * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
	  (4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 4.5 * m_xsi[8][0] * m_xsi[8][0] + m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][1] - 22.5 * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
	  (-27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
	  (-13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
	  (13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
	  (4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] + m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5) },

	{ (-4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][0] * m_xsi[9][0] + 18. * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][1] * m_xsi[9][1] - 5.5 * m_xsi[9][0] - 5.5 * m_xsi[9][1] + 1.) * (-0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] + 27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][0] - 22.5 * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (-13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][0] * m_xsi[9][0] + 4.5 * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 4.5 * m_xsi[9][0] * m_xsi[9][0] + m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][1] - 22.5 * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (-27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (-13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] + m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
	  (-4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][0] * m_xsi[9][0] + 18. * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][1] * m_xsi[9][1] - 5.5 * m_xsi[9][0] - 5.5 * m_xsi[9][1] + 1.) * (0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] + 27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][0] - 22.5 * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
	  (-13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][0] * m_xsi[9][0] + 4.5 * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
	  (4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 4.5 * m_xsi[9][0] * m_xsi[9][0] + m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][1] - 22.5 * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
	  (-27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
	  (-13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
	  (13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
	  (4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] + m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5) },

	{ (-4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][0] * m_xsi[10][0] + 18. * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][1] * m_xsi[10][1] - 5.5 * m_xsi[10][0] - 5.5 * m_xsi[10][1] + 1.) * (-0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] + 27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][0] - 22.5 * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (-13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][0] * m_xsi[10][0] + 4.5 * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 4.5 * m_xsi[10][0] * m_xsi[10][0] + m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][1] - 22.5 * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (-27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (-13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] + m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
	  (-4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][0] * m_xsi[10][0] + 18. * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][1] * m_xsi[10][1] - 5.5 * m_xsi[10][0] - 5.5 * m_xsi[10][1] + 1.) * (0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] + 27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][0] - 22.5 * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
	  (-13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][0] * m_xsi[10][0] + 4.5 * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
	  (4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 4.5 * m_xsi[10][0] * m_xsi[10][0] + m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][1] - 22.5 * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
	  (-27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
	  (-13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
	  (13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
	  (4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] + m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5) },

	{ (-4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][0] * m_xsi[11][0] + 18. * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][1] * m_xsi[11][1] - 5.5 * m_xsi[11][0] - 5.5 * m_xsi[11][1] + 1.) * (-0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] + 27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][0] - 22.5 * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (-13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][0] * m_xsi[11][0] + 4.5 * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 4.5 * m_xsi[11][0] * m_xsi[11][0] + m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][1] - 22.5 * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (-27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (-13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] + m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
	  (-4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][0] * m_xsi[11][0] + 18. * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][1] * m_xsi[11][1] - 5.5 * m_xsi[11][0] - 5.5 * m_xsi[11][1] + 1.) * (0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] + 27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][0] - 22.5 * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
	  (-13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][0] * m_xsi[11][0] + 4.5 * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
	  (4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 4.5 * m_xsi[11][0] * m_xsi[11][0] + m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][1] - 22.5 * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
	  (-27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
	  (-13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
	  (13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
	  (4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] + m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5) },

	{ (-4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][0] * m_xsi[12][0] + 18. * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][1] * m_xsi[12][1] - 5.5 * m_xsi[12][0] - 5.5 * m_xsi[12][1] + 1.) * (-0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] + 27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][0] - 22.5 * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (-13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][0] * m_xsi[12][0] + 4.5 * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 4.5 * m_xsi[12][0] * m_xsi[12][0] + m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][1] - 22.5 * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (-27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (-13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] + m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
	  (-4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][0] * m_xsi[12][0] + 18. * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][1] * m_xsi[12][1] - 5.5 * m_xsi[12][0] - 5.5 * m_xsi[12][1] + 1.) * (0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] + 27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][0] - 22.5 * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
	  (-13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][0] * m_xsi[12][0] + 4.5 * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
	  (4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 4.5 * m_xsi[12][0] * m_xsi[12][0] + m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][1] - 22.5 * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
	  (-27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
	  (-13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
	  (13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
	  (4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] + m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5) },

	{ (-4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][0] * m_xsi[13][0] + 18. * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][1] * m_xsi[13][1] - 5.5 * m_xsi[13][0] - 5.5 * m_xsi[13][1] + 1.) * (-0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] + 27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][0] - 22.5 * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (-13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][0] * m_xsi[13][0] + 4.5 * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 4.5 * m_xsi[13][0] * m_xsi[13][0] + m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][1] - 22.5 * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (-27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (-13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] + m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
	  (-4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][0] * m_xsi[13][0] + 18. * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][1] * m_xsi[13][1] - 5.5 * m_xsi[13][0] - 5.5 * m_xsi[13][1] + 1.) * (0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] + 27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][0] - 22.5 * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
	  (-13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][0] * m_xsi[13][0] + 4.5 * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
	  (4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 4.5 * m_xsi[13][0] * m_xsi[13][0] + m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][1] - 22.5 * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
	  (-27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
	  (-13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
	  (13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
	  (4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] + m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5) } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
inline const double O2P2::Prep::Elem::Elem_Pri20::mv_DPsi[mv_numIP][mv_numNodes][mv_Dim] = {
	{ { (-13.5 * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] + 18. * m_xsi[0][0] + 18. * m_xsi[0][1] - 5.5) * (-0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] + 18. * m_xsi[0][0] + 18. * m_xsi[0][1] - 5.5) * (-0.5 * m_xsi[0][2] + 0.5),
		(-4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][1] * m_xsi[0][1] - 5.5 * m_xsi[0][0] - 5.5 * m_xsi[0][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[0][0] * m_xsi[0][0] + 54. * m_xsi[0][0] * m_xsi[0][1] + 13.5 * m_xsi[0][1] * m_xsi[0][1] - 45. * m_xsi[0][0] - 22.5 * m_xsi[0][1] + 9.) * (-0.5 * m_xsi[0][2] + 0.5),
		(27. * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][1] - 22.5 * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][0] - 22.5 * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0]) * -0.5 },
	  { (-40.5 * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] + 36. * m_xsi[0][0] + 4.5 * m_xsi[0][1] - 4.5) * (-0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][0] + 4.5 * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] * m_xsi[0][0] + 4.5 * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0]) * -0.5 },
	  { (13.5 * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] + 1.) * (-0.5 * m_xsi[0][2] + 0.5),
		0.,
		(4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 4.5 * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0]) * -0.5 },
	  { (27. * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] + 54. * m_xsi[0][0] * m_xsi[0][1] + 40.5 * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] - 45. * m_xsi[0][1] + 9.) * (-0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][1] - 22.5 * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][1]) * -0.5 },
	  { (-54. * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][1] * m_xsi[0][1] + 27. * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
		(-27. * m_xsi[0][0] * m_xsi[0][0] - 54. * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
		(-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1]) * -0.5 },
	  { (27. * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] - 4.5 * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * -0.5 },
	  { (-13.5 * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
		(-27. * m_xsi[0][0] * m_xsi[0][1] - 40.5 * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][0] + 36. * m_xsi[0][1] - 4.5) * (-0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * -0.5 },
	  { (13.5 * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * (-0.5 * m_xsi[0][2] + 0.5),
		(27. * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0]) * (-0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[0][1] * m_xsi[0][1] - 9. * m_xsi[0][1] + 1.) * (-0.5 * m_xsi[0][2] + 0.5),
		(4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] + m_xsi[0][1]) * -0.5 },
	  { (-13.5 * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] + 18. * m_xsi[0][0] + 18. * m_xsi[0][1] - 5.5) * (0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] + 18. * m_xsi[0][0] + 18. * m_xsi[0][1] - 5.5) * (0.5 * m_xsi[0][2] + 0.5),
		(-4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][0] * m_xsi[0][0] + 18. * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][1] * m_xsi[0][1] - 5.5 * m_xsi[0][0] - 5.5 * m_xsi[0][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[0][0] * m_xsi[0][0] + 54. * m_xsi[0][0] * m_xsi[0][1] + 13.5 * m_xsi[0][1] * m_xsi[0][1] - 45. * m_xsi[0][0] - 22.5 * m_xsi[0][1] + 9.) * (0.5 * m_xsi[0][2] + 0.5),
		(27. * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][1] - 22.5 * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] + 27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][0] - 22.5 * m_xsi[0][0] * m_xsi[0][1] + 9. * m_xsi[0][0]) * 0.5 },
	  { (-40.5 * m_xsi[0][0] * m_xsi[0][0] - 27. * m_xsi[0][0] * m_xsi[0][1] + 36. * m_xsi[0][0] + 4.5 * m_xsi[0][1] - 4.5) * (0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][0] + 4.5 * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][0] * m_xsi[0][0] + 4.5 * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0]) * 0.5 },
	  { (13.5 * m_xsi[0][0] * m_xsi[0][0] - 9. * m_xsi[0][0] + 1.) * (0.5 * m_xsi[0][2] + 0.5),
		0.,
		(4.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][0] - 4.5 * m_xsi[0][0] * m_xsi[0][0] + m_xsi[0][0]) * 0.5 },
	  { (27. * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] + 54. * m_xsi[0][0] * m_xsi[0][1] + 40.5 * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] - 45. * m_xsi[0][1] + 9.) * (0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 22.5 * m_xsi[0][0] * m_xsi[0][1] - 22.5 * m_xsi[0][1] * m_xsi[0][1] + 9. * m_xsi[0][1]) * 0.5 },
	  { (-54. * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][1] * m_xsi[0][1] + 27. * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
		(-27. * m_xsi[0][0] * m_xsi[0][0] - 54. * m_xsi[0][0] * m_xsi[0][1] + 27. * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
		(-27. * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 27. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] + 27. * m_xsi[0][0] * m_xsi[0][1]) * 0.5 },
	  { (27. * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] - 4.5 * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * 0.5 },
	  { (-13.5 * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
		(-27. * m_xsi[0][0] * m_xsi[0][1] - 40.5 * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][0] + 36. * m_xsi[0][1] - 4.5) * (0.5 * m_xsi[0][2] + 0.5),
		(-13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 13.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] + 4.5 * m_xsi[0][0] * m_xsi[0][1] + 18. * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * 0.5 },
	  { (13.5 * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1]) * (0.5 * m_xsi[0][2] + 0.5),
		(27. * m_xsi[0][0] * m_xsi[0][1] - 4.5 * m_xsi[0][0]) * (0.5 * m_xsi[0][2] + 0.5),
		(13.5 * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][0] * m_xsi[0][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[0][1] * m_xsi[0][1] - 9. * m_xsi[0][1] + 1.) * (0.5 * m_xsi[0][2] + 0.5),
		(4.5 * m_xsi[0][1] * m_xsi[0][1] * m_xsi[0][1] - 4.5 * m_xsi[0][1] * m_xsi[0][1] + m_xsi[0][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] + 18. * m_xsi[1][0] + 18. * m_xsi[1][1] - 5.5) * (-0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] + 18. * m_xsi[1][0] + 18. * m_xsi[1][1] - 5.5) * (-0.5 * m_xsi[1][2] + 0.5),
		(-4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][1] * m_xsi[1][1] - 5.5 * m_xsi[1][0] - 5.5 * m_xsi[1][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[1][0] * m_xsi[1][0] + 54. * m_xsi[1][0] * m_xsi[1][1] + 13.5 * m_xsi[1][1] * m_xsi[1][1] - 45. * m_xsi[1][0] - 22.5 * m_xsi[1][1] + 9.) * (-0.5 * m_xsi[1][2] + 0.5),
		(27. * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][1] - 22.5 * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][0] - 22.5 * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0]) * -0.5 },
	  { (-40.5 * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] + 36. * m_xsi[1][0] + 4.5 * m_xsi[1][1] - 4.5) * (-0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][0] + 4.5 * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] * m_xsi[1][0] + 4.5 * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0]) * -0.5 },
	  { (13.5 * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] + 1.) * (-0.5 * m_xsi[1][2] + 0.5),
		0.,
		(4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 4.5 * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0]) * -0.5 },
	  { (27. * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] + 54. * m_xsi[1][0] * m_xsi[1][1] + 40.5 * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] - 45. * m_xsi[1][1] + 9.) * (-0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][1] - 22.5 * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][1]) * -0.5 },
	  { (-54. * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][1] * m_xsi[1][1] + 27. * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
		(-27. * m_xsi[1][0] * m_xsi[1][0] - 54. * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
		(-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1]) * -0.5 },
	  { (27. * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] - 4.5 * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * -0.5 },
	  { (-13.5 * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
		(-27. * m_xsi[1][0] * m_xsi[1][1] - 40.5 * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][0] + 36. * m_xsi[1][1] - 4.5) * (-0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * -0.5 },
	  { (13.5 * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * (-0.5 * m_xsi[1][2] + 0.5),
		(27. * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0]) * (-0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[1][1] * m_xsi[1][1] - 9. * m_xsi[1][1] + 1.) * (-0.5 * m_xsi[1][2] + 0.5),
		(4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] + m_xsi[1][1]) * -0.5 },
	  { (-13.5 * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] + 18. * m_xsi[1][0] + 18. * m_xsi[1][1] - 5.5) * (0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] + 18. * m_xsi[1][0] + 18. * m_xsi[1][1] - 5.5) * (0.5 * m_xsi[1][2] + 0.5),
		(-4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][0] * m_xsi[1][0] + 18. * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][1] * m_xsi[1][1] - 5.5 * m_xsi[1][0] - 5.5 * m_xsi[1][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[1][0] * m_xsi[1][0] + 54. * m_xsi[1][0] * m_xsi[1][1] + 13.5 * m_xsi[1][1] * m_xsi[1][1] - 45. * m_xsi[1][0] - 22.5 * m_xsi[1][1] + 9.) * (0.5 * m_xsi[1][2] + 0.5),
		(27. * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][1] - 22.5 * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] + 27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][0] - 22.5 * m_xsi[1][0] * m_xsi[1][1] + 9. * m_xsi[1][0]) * 0.5 },
	  { (-40.5 * m_xsi[1][0] * m_xsi[1][0] - 27. * m_xsi[1][0] * m_xsi[1][1] + 36. * m_xsi[1][0] + 4.5 * m_xsi[1][1] - 4.5) * (0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][0] + 4.5 * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][0] * m_xsi[1][0] + 4.5 * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0]) * 0.5 },
	  { (13.5 * m_xsi[1][0] * m_xsi[1][0] - 9. * m_xsi[1][0] + 1.) * (0.5 * m_xsi[1][2] + 0.5),
		0.,
		(4.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][0] - 4.5 * m_xsi[1][0] * m_xsi[1][0] + m_xsi[1][0]) * 0.5 },
	  { (27. * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] + 54. * m_xsi[1][0] * m_xsi[1][1] + 40.5 * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] - 45. * m_xsi[1][1] + 9.) * (0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 22.5 * m_xsi[1][0] * m_xsi[1][1] - 22.5 * m_xsi[1][1] * m_xsi[1][1] + 9. * m_xsi[1][1]) * 0.5 },
	  { (-54. * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][1] * m_xsi[1][1] + 27. * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
		(-27. * m_xsi[1][0] * m_xsi[1][0] - 54. * m_xsi[1][0] * m_xsi[1][1] + 27. * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
		(-27. * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 27. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] + 27. * m_xsi[1][0] * m_xsi[1][1]) * 0.5 },
	  { (27. * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] - 4.5 * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * 0.5 },
	  { (-13.5 * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
		(-27. * m_xsi[1][0] * m_xsi[1][1] - 40.5 * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][0] + 36. * m_xsi[1][1] - 4.5) * (0.5 * m_xsi[1][2] + 0.5),
		(-13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 13.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] + 4.5 * m_xsi[1][0] * m_xsi[1][1] + 18. * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * 0.5 },
	  { (13.5 * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1]) * (0.5 * m_xsi[1][2] + 0.5),
		(27. * m_xsi[1][0] * m_xsi[1][1] - 4.5 * m_xsi[1][0]) * (0.5 * m_xsi[1][2] + 0.5),
		(13.5 * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][0] * m_xsi[1][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[1][1] * m_xsi[1][1] - 9. * m_xsi[1][1] + 1.) * (0.5 * m_xsi[1][2] + 0.5),
		(4.5 * m_xsi[1][1] * m_xsi[1][1] * m_xsi[1][1] - 4.5 * m_xsi[1][1] * m_xsi[1][1] + m_xsi[1][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] + 18. * m_xsi[2][0] + 18. * m_xsi[2][1] - 5.5) * (-0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] + 18. * m_xsi[2][0] + 18. * m_xsi[2][1] - 5.5) * (-0.5 * m_xsi[2][2] + 0.5),
		(-4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][1] * m_xsi[2][1] - 5.5 * m_xsi[2][0] - 5.5 * m_xsi[2][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[2][0] * m_xsi[2][0] + 54. * m_xsi[2][0] * m_xsi[2][1] + 13.5 * m_xsi[2][1] * m_xsi[2][1] - 45. * m_xsi[2][0] - 22.5 * m_xsi[2][1] + 9.) * (-0.5 * m_xsi[2][2] + 0.5),
		(27. * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][1] - 22.5 * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][0] - 22.5 * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0]) * -0.5 },
	  { (-40.5 * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] + 36. * m_xsi[2][0] + 4.5 * m_xsi[2][1] - 4.5) * (-0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][0] + 4.5 * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] * m_xsi[2][0] + 4.5 * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0]) * -0.5 },
	  { (13.5 * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] + 1.) * (-0.5 * m_xsi[2][2] + 0.5),
		0.,
		(4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 4.5 * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0]) * -0.5 },
	  { (27. * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] + 54. * m_xsi[2][0] * m_xsi[2][1] + 40.5 * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] - 45. * m_xsi[2][1] + 9.) * (-0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][1] - 22.5 * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][1]) * -0.5 },
	  { (-54. * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][1] * m_xsi[2][1] + 27. * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
		(-27. * m_xsi[2][0] * m_xsi[2][0] - 54. * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
		(-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1]) * -0.5 },
	  { (27. * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] - 4.5 * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * -0.5 },
	  { (-13.5 * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
		(-27. * m_xsi[2][0] * m_xsi[2][1] - 40.5 * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][0] + 36. * m_xsi[2][1] - 4.5) * (-0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * -0.5 },
	  { (13.5 * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * (-0.5 * m_xsi[2][2] + 0.5),
		(27. * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0]) * (-0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[2][1] * m_xsi[2][1] - 9. * m_xsi[2][1] + 1.) * (-0.5 * m_xsi[2][2] + 0.5),
		(4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] + m_xsi[2][1]) * -0.5 },
	  { (-13.5 * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] + 18. * m_xsi[2][0] + 18. * m_xsi[2][1] - 5.5) * (0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] + 18. * m_xsi[2][0] + 18. * m_xsi[2][1] - 5.5) * (0.5 * m_xsi[2][2] + 0.5),
		(-4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][0] * m_xsi[2][0] + 18. * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][1] * m_xsi[2][1] - 5.5 * m_xsi[2][0] - 5.5 * m_xsi[2][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[2][0] * m_xsi[2][0] + 54. * m_xsi[2][0] * m_xsi[2][1] + 13.5 * m_xsi[2][1] * m_xsi[2][1] - 45. * m_xsi[2][0] - 22.5 * m_xsi[2][1] + 9.) * (0.5 * m_xsi[2][2] + 0.5),
		(27. * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][1] - 22.5 * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] + 27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][0] - 22.5 * m_xsi[2][0] * m_xsi[2][1] + 9. * m_xsi[2][0]) * 0.5 },
	  { (-40.5 * m_xsi[2][0] * m_xsi[2][0] - 27. * m_xsi[2][0] * m_xsi[2][1] + 36. * m_xsi[2][0] + 4.5 * m_xsi[2][1] - 4.5) * (0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][0] + 4.5 * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][0] * m_xsi[2][0] + 4.5 * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0]) * 0.5 },
	  { (13.5 * m_xsi[2][0] * m_xsi[2][0] - 9. * m_xsi[2][0] + 1.) * (0.5 * m_xsi[2][2] + 0.5),
		0.,
		(4.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][0] - 4.5 * m_xsi[2][0] * m_xsi[2][0] + m_xsi[2][0]) * 0.5 },
	  { (27. * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] + 54. * m_xsi[2][0] * m_xsi[2][1] + 40.5 * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] - 45. * m_xsi[2][1] + 9.) * (0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 22.5 * m_xsi[2][0] * m_xsi[2][1] - 22.5 * m_xsi[2][1] * m_xsi[2][1] + 9. * m_xsi[2][1]) * 0.5 },
	  { (-54. * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][1] * m_xsi[2][1] + 27. * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
		(-27. * m_xsi[2][0] * m_xsi[2][0] - 54. * m_xsi[2][0] * m_xsi[2][1] + 27. * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
		(-27. * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 27. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] + 27. * m_xsi[2][0] * m_xsi[2][1]) * 0.5 },
	  { (27. * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] - 4.5 * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * 0.5 },
	  { (-13.5 * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
		(-27. * m_xsi[2][0] * m_xsi[2][1] - 40.5 * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][0] + 36. * m_xsi[2][1] - 4.5) * (0.5 * m_xsi[2][2] + 0.5),
		(-13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 13.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] + 4.5 * m_xsi[2][0] * m_xsi[2][1] + 18. * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * 0.5 },
	  { (13.5 * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1]) * (0.5 * m_xsi[2][2] + 0.5),
		(27. * m_xsi[2][0] * m_xsi[2][1] - 4.5 * m_xsi[2][0]) * (0.5 * m_xsi[2][2] + 0.5),
		(13.5 * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][0] * m_xsi[2][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[2][1] * m_xsi[2][1] - 9. * m_xsi[2][1] + 1.) * (0.5 * m_xsi[2][2] + 0.5),
		(4.5 * m_xsi[2][1] * m_xsi[2][1] * m_xsi[2][1] - 4.5 * m_xsi[2][1] * m_xsi[2][1] + m_xsi[2][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] + 18. * m_xsi[3][0] + 18. * m_xsi[3][1] - 5.5) * (-0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] + 18. * m_xsi[3][0] + 18. * m_xsi[3][1] - 5.5) * (-0.5 * m_xsi[3][2] + 0.5),
		(-4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][1] * m_xsi[3][1] - 5.5 * m_xsi[3][0] - 5.5 * m_xsi[3][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[3][0] * m_xsi[3][0] + 54. * m_xsi[3][0] * m_xsi[3][1] + 13.5 * m_xsi[3][1] * m_xsi[3][1] - 45. * m_xsi[3][0] - 22.5 * m_xsi[3][1] + 9.) * (-0.5 * m_xsi[3][2] + 0.5),
		(27. * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][1] - 22.5 * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][0] - 22.5 * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0]) * -0.5 },
	  { (-40.5 * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] + 36. * m_xsi[3][0] + 4.5 * m_xsi[3][1] - 4.5) * (-0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][0] + 4.5 * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] * m_xsi[3][0] + 4.5 * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0]) * -0.5 },
	  { (13.5 * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] + 1.) * (-0.5 * m_xsi[3][2] + 0.5),
		0.,
		(4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 4.5 * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0]) * -0.5 },
	  { (27. * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] + 54. * m_xsi[3][0] * m_xsi[3][1] + 40.5 * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] - 45. * m_xsi[3][1] + 9.) * (-0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][1] - 22.5 * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][1]) * -0.5 },
	  { (-54. * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][1] * m_xsi[3][1] + 27. * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
		(-27. * m_xsi[3][0] * m_xsi[3][0] - 54. * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
		(-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1]) * -0.5 },
	  { (27. * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] - 4.5 * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * -0.5 },
	  { (-13.5 * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
		(-27. * m_xsi[3][0] * m_xsi[3][1] - 40.5 * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][0] + 36. * m_xsi[3][1] - 4.5) * (-0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * -0.5 },
	  { (13.5 * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * (-0.5 * m_xsi[3][2] + 0.5),
		(27. * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0]) * (-0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[3][1] * m_xsi[3][1] - 9. * m_xsi[3][1] + 1.) * (-0.5 * m_xsi[3][2] + 0.5),
		(4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] + m_xsi[3][1]) * -0.5 },
	  { (-13.5 * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] + 18. * m_xsi[3][0] + 18. * m_xsi[3][1] - 5.5) * (0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] + 18. * m_xsi[3][0] + 18. * m_xsi[3][1] - 5.5) * (0.5 * m_xsi[3][2] + 0.5),
		(-4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][0] * m_xsi[3][0] + 18. * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][1] * m_xsi[3][1] - 5.5 * m_xsi[3][0] - 5.5 * m_xsi[3][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[3][0] * m_xsi[3][0] + 54. * m_xsi[3][0] * m_xsi[3][1] + 13.5 * m_xsi[3][1] * m_xsi[3][1] - 45. * m_xsi[3][0] - 22.5 * m_xsi[3][1] + 9.) * (0.5 * m_xsi[3][2] + 0.5),
		(27. * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][1] - 22.5 * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] + 27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][0] - 22.5 * m_xsi[3][0] * m_xsi[3][1] + 9. * m_xsi[3][0]) * 0.5 },
	  { (-40.5 * m_xsi[3][0] * m_xsi[3][0] - 27. * m_xsi[3][0] * m_xsi[3][1] + 36. * m_xsi[3][0] + 4.5 * m_xsi[3][1] - 4.5) * (0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][0] + 4.5 * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][0] * m_xsi[3][0] + 4.5 * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0]) * 0.5 },
	  { (13.5 * m_xsi[3][0] * m_xsi[3][0] - 9. * m_xsi[3][0] + 1.) * (0.5 * m_xsi[3][2] + 0.5),
		0.,
		(4.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][0] - 4.5 * m_xsi[3][0] * m_xsi[3][0] + m_xsi[3][0]) * 0.5 },
	  { (27. * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] + 54. * m_xsi[3][0] * m_xsi[3][1] + 40.5 * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] - 45. * m_xsi[3][1] + 9.) * (0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 22.5 * m_xsi[3][0] * m_xsi[3][1] - 22.5 * m_xsi[3][1] * m_xsi[3][1] + 9. * m_xsi[3][1]) * 0.5 },
	  { (-54. * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][1] * m_xsi[3][1] + 27. * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
		(-27. * m_xsi[3][0] * m_xsi[3][0] - 54. * m_xsi[3][0] * m_xsi[3][1] + 27. * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
		(-27. * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 27. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] + 27. * m_xsi[3][0] * m_xsi[3][1]) * 0.5 },
	  { (27. * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] - 4.5 * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * 0.5 },
	  { (-13.5 * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
		(-27. * m_xsi[3][0] * m_xsi[3][1] - 40.5 * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][0] + 36. * m_xsi[3][1] - 4.5) * (0.5 * m_xsi[3][2] + 0.5),
		(-13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 13.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] + 4.5 * m_xsi[3][0] * m_xsi[3][1] + 18. * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * 0.5 },
	  { (13.5 * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1]) * (0.5 * m_xsi[3][2] + 0.5),
		(27. * m_xsi[3][0] * m_xsi[3][1] - 4.5 * m_xsi[3][0]) * (0.5 * m_xsi[3][2] + 0.5),
		(13.5 * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][0] * m_xsi[3][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[3][1] * m_xsi[3][1] - 9. * m_xsi[3][1] + 1.) * (0.5 * m_xsi[3][2] + 0.5),
		(4.5 * m_xsi[3][1] * m_xsi[3][1] * m_xsi[3][1] - 4.5 * m_xsi[3][1] * m_xsi[3][1] + m_xsi[3][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] + 18. * m_xsi[4][0] + 18. * m_xsi[4][1] - 5.5) * (-0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] + 18. * m_xsi[4][0] + 18. * m_xsi[4][1] - 5.5) * (-0.5 * m_xsi[4][2] + 0.5),
		(-4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][1] * m_xsi[4][1] - 5.5 * m_xsi[4][0] - 5.5 * m_xsi[4][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[4][0] * m_xsi[4][0] + 54. * m_xsi[4][0] * m_xsi[4][1] + 13.5 * m_xsi[4][1] * m_xsi[4][1] - 45. * m_xsi[4][0] - 22.5 * m_xsi[4][1] + 9.) * (-0.5 * m_xsi[4][2] + 0.5),
		(27. * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][1] - 22.5 * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][0] - 22.5 * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0]) * -0.5 },
	  { (-40.5 * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] + 36. * m_xsi[4][0] + 4.5 * m_xsi[4][1] - 4.5) * (-0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][0] + 4.5 * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] * m_xsi[4][0] + 4.5 * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0]) * -0.5 },
	  { (13.5 * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] + 1.) * (-0.5 * m_xsi[4][2] + 0.5),
		0.,
		(4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 4.5 * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0]) * -0.5 },
	  { (27. * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] + 54. * m_xsi[4][0] * m_xsi[4][1] + 40.5 * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] - 45. * m_xsi[4][1] + 9.) * (-0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][1] - 22.5 * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][1]) * -0.5 },
	  { (-54. * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][1] * m_xsi[4][1] + 27. * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
		(-27. * m_xsi[4][0] * m_xsi[4][0] - 54. * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
		(-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1]) * -0.5 },
	  { (27. * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] - 4.5 * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * -0.5 },
	  { (-13.5 * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
		(-27. * m_xsi[4][0] * m_xsi[4][1] - 40.5 * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][0] + 36. * m_xsi[4][1] - 4.5) * (-0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * -0.5 },
	  { (13.5 * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * (-0.5 * m_xsi[4][2] + 0.5),
		(27. * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0]) * (-0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[4][1] * m_xsi[4][1] - 9. * m_xsi[4][1] + 1.) * (-0.5 * m_xsi[4][2] + 0.5),
		(4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] + m_xsi[4][1]) * -0.5 },
	  { (-13.5 * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] + 18. * m_xsi[4][0] + 18. * m_xsi[4][1] - 5.5) * (0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] + 18. * m_xsi[4][0] + 18. * m_xsi[4][1] - 5.5) * (0.5 * m_xsi[4][2] + 0.5),
		(-4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][0] * m_xsi[4][0] + 18. * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][1] * m_xsi[4][1] - 5.5 * m_xsi[4][0] - 5.5 * m_xsi[4][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[4][0] * m_xsi[4][0] + 54. * m_xsi[4][0] * m_xsi[4][1] + 13.5 * m_xsi[4][1] * m_xsi[4][1] - 45. * m_xsi[4][0] - 22.5 * m_xsi[4][1] + 9.) * (0.5 * m_xsi[4][2] + 0.5),
		(27. * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][1] - 22.5 * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] + 27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][0] - 22.5 * m_xsi[4][0] * m_xsi[4][1] + 9. * m_xsi[4][0]) * 0.5 },
	  { (-40.5 * m_xsi[4][0] * m_xsi[4][0] - 27. * m_xsi[4][0] * m_xsi[4][1] + 36. * m_xsi[4][0] + 4.5 * m_xsi[4][1] - 4.5) * (0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][0] + 4.5 * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][0] * m_xsi[4][0] + 4.5 * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0]) * 0.5 },
	  { (13.5 * m_xsi[4][0] * m_xsi[4][0] - 9. * m_xsi[4][0] + 1.) * (0.5 * m_xsi[4][2] + 0.5),
		0.,
		(4.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][0] - 4.5 * m_xsi[4][0] * m_xsi[4][0] + m_xsi[4][0]) * 0.5 },
	  { (27. * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] + 54. * m_xsi[4][0] * m_xsi[4][1] + 40.5 * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] - 45. * m_xsi[4][1] + 9.) * (0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 22.5 * m_xsi[4][0] * m_xsi[4][1] - 22.5 * m_xsi[4][1] * m_xsi[4][1] + 9. * m_xsi[4][1]) * 0.5 },
	  { (-54. * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][1] * m_xsi[4][1] + 27. * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
		(-27. * m_xsi[4][0] * m_xsi[4][0] - 54. * m_xsi[4][0] * m_xsi[4][1] + 27. * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
		(-27. * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 27. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] + 27. * m_xsi[4][0] * m_xsi[4][1]) * 0.5 },
	  { (27. * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] - 4.5 * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * 0.5 },
	  { (-13.5 * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
		(-27. * m_xsi[4][0] * m_xsi[4][1] - 40.5 * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][0] + 36. * m_xsi[4][1] - 4.5) * (0.5 * m_xsi[4][2] + 0.5),
		(-13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 13.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] + 4.5 * m_xsi[4][0] * m_xsi[4][1] + 18. * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * 0.5 },
	  { (13.5 * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1]) * (0.5 * m_xsi[4][2] + 0.5),
		(27. * m_xsi[4][0] * m_xsi[4][1] - 4.5 * m_xsi[4][0]) * (0.5 * m_xsi[4][2] + 0.5),
		(13.5 * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][0] * m_xsi[4][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[4][1] * m_xsi[4][1] - 9. * m_xsi[4][1] + 1.) * (0.5 * m_xsi[4][2] + 0.5),
		(4.5 * m_xsi[4][1] * m_xsi[4][1] * m_xsi[4][1] - 4.5 * m_xsi[4][1] * m_xsi[4][1] + m_xsi[4][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] + 18. * m_xsi[5][0] + 18. * m_xsi[5][1] - 5.5) * (-0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] + 18. * m_xsi[5][0] + 18. * m_xsi[5][1] - 5.5) * (-0.5 * m_xsi[5][2] + 0.5),
		(-4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][1] * m_xsi[5][1] - 5.5 * m_xsi[5][0] - 5.5 * m_xsi[5][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[5][0] * m_xsi[5][0] + 54. * m_xsi[5][0] * m_xsi[5][1] + 13.5 * m_xsi[5][1] * m_xsi[5][1] - 45. * m_xsi[5][0] - 22.5 * m_xsi[5][1] + 9.) * (-0.5 * m_xsi[5][2] + 0.5),
		(27. * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][1] - 22.5 * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][0] - 22.5 * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0]) * -0.5 },
	  { (-40.5 * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] + 36. * m_xsi[5][0] + 4.5 * m_xsi[5][1] - 4.5) * (-0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][0] + 4.5 * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] * m_xsi[5][0] + 4.5 * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0]) * -0.5 },
	  { (13.5 * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] + 1.) * (-0.5 * m_xsi[5][2] + 0.5),
		0.,
		(4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 4.5 * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0]) * -0.5 },
	  { (27. * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] + 54. * m_xsi[5][0] * m_xsi[5][1] + 40.5 * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] - 45. * m_xsi[5][1] + 9.) * (-0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][1] - 22.5 * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][1]) * -0.5 },
	  { (-54. * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][1] * m_xsi[5][1] + 27. * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
		(-27. * m_xsi[5][0] * m_xsi[5][0] - 54. * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
		(-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1]) * -0.5 },
	  { (27. * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] - 4.5 * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * -0.5 },
	  { (-13.5 * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
		(-27. * m_xsi[5][0] * m_xsi[5][1] - 40.5 * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][0] + 36. * m_xsi[5][1] - 4.5) * (-0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * -0.5 },
	  { (13.5 * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * (-0.5 * m_xsi[5][2] + 0.5),
		(27. * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0]) * (-0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[5][1] * m_xsi[5][1] - 9. * m_xsi[5][1] + 1.) * (-0.5 * m_xsi[5][2] + 0.5),
		(4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] + m_xsi[5][1]) * -0.5 },
	  { (-13.5 * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] + 18. * m_xsi[5][0] + 18. * m_xsi[5][1] - 5.5) * (0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] + 18. * m_xsi[5][0] + 18. * m_xsi[5][1] - 5.5) * (0.5 * m_xsi[5][2] + 0.5),
		(-4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][0] * m_xsi[5][0] + 18. * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][1] * m_xsi[5][1] - 5.5 * m_xsi[5][0] - 5.5 * m_xsi[5][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[5][0] * m_xsi[5][0] + 54. * m_xsi[5][0] * m_xsi[5][1] + 13.5 * m_xsi[5][1] * m_xsi[5][1] - 45. * m_xsi[5][0] - 22.5 * m_xsi[5][1] + 9.) * (0.5 * m_xsi[5][2] + 0.5),
		(27. * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][1] - 22.5 * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] + 27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][0] - 22.5 * m_xsi[5][0] * m_xsi[5][1] + 9. * m_xsi[5][0]) * 0.5 },
	  { (-40.5 * m_xsi[5][0] * m_xsi[5][0] - 27. * m_xsi[5][0] * m_xsi[5][1] + 36. * m_xsi[5][0] + 4.5 * m_xsi[5][1] - 4.5) * (0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][0] + 4.5 * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][0] * m_xsi[5][0] + 4.5 * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0]) * 0.5 },
	  { (13.5 * m_xsi[5][0] * m_xsi[5][0] - 9. * m_xsi[5][0] + 1.) * (0.5 * m_xsi[5][2] + 0.5),
		0.,
		(4.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][0] - 4.5 * m_xsi[5][0] * m_xsi[5][0] + m_xsi[5][0]) * 0.5 },
	  { (27. * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] + 54. * m_xsi[5][0] * m_xsi[5][1] + 40.5 * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] - 45. * m_xsi[5][1] + 9.) * (0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 22.5 * m_xsi[5][0] * m_xsi[5][1] - 22.5 * m_xsi[5][1] * m_xsi[5][1] + 9. * m_xsi[5][1]) * 0.5 },
	  { (-54. * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][1] * m_xsi[5][1] + 27. * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
		(-27. * m_xsi[5][0] * m_xsi[5][0] - 54. * m_xsi[5][0] * m_xsi[5][1] + 27. * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
		(-27. * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 27. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] + 27. * m_xsi[5][0] * m_xsi[5][1]) * 0.5 },
	  { (27. * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] - 4.5 * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * 0.5 },
	  { (-13.5 * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
		(-27. * m_xsi[5][0] * m_xsi[5][1] - 40.5 * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][0] + 36. * m_xsi[5][1] - 4.5) * (0.5 * m_xsi[5][2] + 0.5),
		(-13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 13.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] + 4.5 * m_xsi[5][0] * m_xsi[5][1] + 18. * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * 0.5 },
	  { (13.5 * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1]) * (0.5 * m_xsi[5][2] + 0.5),
		(27. * m_xsi[5][0] * m_xsi[5][1] - 4.5 * m_xsi[5][0]) * (0.5 * m_xsi[5][2] + 0.5),
		(13.5 * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][0] * m_xsi[5][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[5][1] * m_xsi[5][1] - 9. * m_xsi[5][1] + 1.) * (0.5 * m_xsi[5][2] + 0.5),
		(4.5 * m_xsi[5][1] * m_xsi[5][1] * m_xsi[5][1] - 4.5 * m_xsi[5][1] * m_xsi[5][1] + m_xsi[5][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] + 18. * m_xsi[6][0] + 18. * m_xsi[6][1] - 5.5) * (-0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] + 18. * m_xsi[6][0] + 18. * m_xsi[6][1] - 5.5) * (-0.5 * m_xsi[6][2] + 0.5),
		(-4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][1] * m_xsi[6][1] - 5.5 * m_xsi[6][0] - 5.5 * m_xsi[6][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[6][0] * m_xsi[6][0] + 54. * m_xsi[6][0] * m_xsi[6][1] + 13.5 * m_xsi[6][1] * m_xsi[6][1] - 45. * m_xsi[6][0] - 22.5 * m_xsi[6][1] + 9.) * (-0.5 * m_xsi[6][2] + 0.5),
		(27. * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][1] - 22.5 * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][0] - 22.5 * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0]) * -0.5 },
	  { (-40.5 * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] + 36. * m_xsi[6][0] + 4.5 * m_xsi[6][1] - 4.5) * (-0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][0] + 4.5 * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] * m_xsi[6][0] + 4.5 * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0]) * -0.5 },
	  { (13.5 * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] + 1.) * (-0.5 * m_xsi[6][2] + 0.5),
		0.,
		(4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 4.5 * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0]) * -0.5 },
	  { (27. * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] + 54. * m_xsi[6][0] * m_xsi[6][1] + 40.5 * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] - 45. * m_xsi[6][1] + 9.) * (-0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][1] - 22.5 * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][1]) * -0.5 },
	  { (-54. * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][1] * m_xsi[6][1] + 27. * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
		(-27. * m_xsi[6][0] * m_xsi[6][0] - 54. * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
		(-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1]) * -0.5 },
	  { (27. * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] - 4.5 * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * -0.5 },
	  { (-13.5 * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
		(-27. * m_xsi[6][0] * m_xsi[6][1] - 40.5 * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][0] + 36. * m_xsi[6][1] - 4.5) * (-0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * -0.5 },
	  { (13.5 * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * (-0.5 * m_xsi[6][2] + 0.5),
		(27. * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0]) * (-0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[6][1] * m_xsi[6][1] - 9. * m_xsi[6][1] + 1.) * (-0.5 * m_xsi[6][2] + 0.5),
		(4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] + m_xsi[6][1]) * -0.5 },
	  { (-13.5 * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] + 18. * m_xsi[6][0] + 18. * m_xsi[6][1] - 5.5) * (0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] + 18. * m_xsi[6][0] + 18. * m_xsi[6][1] - 5.5) * (0.5 * m_xsi[6][2] + 0.5),
		(-4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][0] * m_xsi[6][0] + 18. * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][1] * m_xsi[6][1] - 5.5 * m_xsi[6][0] - 5.5 * m_xsi[6][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[6][0] * m_xsi[6][0] + 54. * m_xsi[6][0] * m_xsi[6][1] + 13.5 * m_xsi[6][1] * m_xsi[6][1] - 45. * m_xsi[6][0] - 22.5 * m_xsi[6][1] + 9.) * (0.5 * m_xsi[6][2] + 0.5),
		(27. * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][1] - 22.5 * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] + 27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][0] - 22.5 * m_xsi[6][0] * m_xsi[6][1] + 9. * m_xsi[6][0]) * 0.5 },
	  { (-40.5 * m_xsi[6][0] * m_xsi[6][0] - 27. * m_xsi[6][0] * m_xsi[6][1] + 36. * m_xsi[6][0] + 4.5 * m_xsi[6][1] - 4.5) * (0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][0] + 4.5 * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][0] * m_xsi[6][0] + 4.5 * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0]) * 0.5 },
	  { (13.5 * m_xsi[6][0] * m_xsi[6][0] - 9. * m_xsi[6][0] + 1.) * (0.5 * m_xsi[6][2] + 0.5),
		0.,
		(4.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][0] - 4.5 * m_xsi[6][0] * m_xsi[6][0] + m_xsi[6][0]) * 0.5 },
	  { (27. * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] + 54. * m_xsi[6][0] * m_xsi[6][1] + 40.5 * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] - 45. * m_xsi[6][1] + 9.) * (0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 22.5 * m_xsi[6][0] * m_xsi[6][1] - 22.5 * m_xsi[6][1] * m_xsi[6][1] + 9. * m_xsi[6][1]) * 0.5 },
	  { (-54. * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][1] * m_xsi[6][1] + 27. * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
		(-27. * m_xsi[6][0] * m_xsi[6][0] - 54. * m_xsi[6][0] * m_xsi[6][1] + 27. * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
		(-27. * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 27. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] + 27. * m_xsi[6][0] * m_xsi[6][1]) * 0.5 },
	  { (27. * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] - 4.5 * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * 0.5 },
	  { (-13.5 * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
		(-27. * m_xsi[6][0] * m_xsi[6][1] - 40.5 * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][0] + 36. * m_xsi[6][1] - 4.5) * (0.5 * m_xsi[6][2] + 0.5),
		(-13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 13.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] + 4.5 * m_xsi[6][0] * m_xsi[6][1] + 18. * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * 0.5 },
	  { (13.5 * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1]) * (0.5 * m_xsi[6][2] + 0.5),
		(27. * m_xsi[6][0] * m_xsi[6][1] - 4.5 * m_xsi[6][0]) * (0.5 * m_xsi[6][2] + 0.5),
		(13.5 * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][0] * m_xsi[6][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[6][1] * m_xsi[6][1] - 9. * m_xsi[6][1] + 1.) * (0.5 * m_xsi[6][2] + 0.5),
		(4.5 * m_xsi[6][1] * m_xsi[6][1] * m_xsi[6][1] - 4.5 * m_xsi[6][1] * m_xsi[6][1] + m_xsi[6][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] + 18. * m_xsi[7][0] + 18. * m_xsi[7][1] - 5.5) * (-0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] + 18. * m_xsi[7][0] + 18. * m_xsi[7][1] - 5.5) * (-0.5 * m_xsi[7][2] + 0.5),
		(-4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][1] * m_xsi[7][1] - 5.5 * m_xsi[7][0] - 5.5 * m_xsi[7][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[7][0] * m_xsi[7][0] + 54. * m_xsi[7][0] * m_xsi[7][1] + 13.5 * m_xsi[7][1] * m_xsi[7][1] - 45. * m_xsi[7][0] - 22.5 * m_xsi[7][1] + 9.) * (-0.5 * m_xsi[7][2] + 0.5),
		(27. * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][1] - 22.5 * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][0] - 22.5 * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0]) * -0.5 },
	  { (-40.5 * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] + 36. * m_xsi[7][0] + 4.5 * m_xsi[7][1] - 4.5) * (-0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][0] + 4.5 * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] * m_xsi[7][0] + 4.5 * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0]) * -0.5 },
	  { (13.5 * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] + 1.) * (-0.5 * m_xsi[7][2] + 0.5),
		0.,
		(4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 4.5 * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0]) * -0.5 },
	  { (27. * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] + 54. * m_xsi[7][0] * m_xsi[7][1] + 40.5 * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] - 45. * m_xsi[7][1] + 9.) * (-0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][1] - 22.5 * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][1]) * -0.5 },
	  { (-54. * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][1] * m_xsi[7][1] + 27. * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
		(-27. * m_xsi[7][0] * m_xsi[7][0] - 54. * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
		(-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1]) * -0.5 },
	  { (27. * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] - 4.5 * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * -0.5 },
	  { (-13.5 * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
		(-27. * m_xsi[7][0] * m_xsi[7][1] - 40.5 * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][0] + 36. * m_xsi[7][1] - 4.5) * (-0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * -0.5 },
	  { (13.5 * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * (-0.5 * m_xsi[7][2] + 0.5),
		(27. * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0]) * (-0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[7][1] * m_xsi[7][1] - 9. * m_xsi[7][1] + 1.) * (-0.5 * m_xsi[7][2] + 0.5),
		(4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] + m_xsi[7][1]) * -0.5 },
	  { (-13.5 * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] + 18. * m_xsi[7][0] + 18. * m_xsi[7][1] - 5.5) * (0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] + 18. * m_xsi[7][0] + 18. * m_xsi[7][1] - 5.5) * (0.5 * m_xsi[7][2] + 0.5),
		(-4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][0] * m_xsi[7][0] + 18. * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][1] * m_xsi[7][1] - 5.5 * m_xsi[7][0] - 5.5 * m_xsi[7][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[7][0] * m_xsi[7][0] + 54. * m_xsi[7][0] * m_xsi[7][1] + 13.5 * m_xsi[7][1] * m_xsi[7][1] - 45. * m_xsi[7][0] - 22.5 * m_xsi[7][1] + 9.) * (0.5 * m_xsi[7][2] + 0.5),
		(27. * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][1] - 22.5 * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] + 27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][0] - 22.5 * m_xsi[7][0] * m_xsi[7][1] + 9. * m_xsi[7][0]) * 0.5 },
	  { (-40.5 * m_xsi[7][0] * m_xsi[7][0] - 27. * m_xsi[7][0] * m_xsi[7][1] + 36. * m_xsi[7][0] + 4.5 * m_xsi[7][1] - 4.5) * (0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][0] + 4.5 * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][0] * m_xsi[7][0] + 4.5 * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0]) * 0.5 },
	  { (13.5 * m_xsi[7][0] * m_xsi[7][0] - 9. * m_xsi[7][0] + 1.) * (0.5 * m_xsi[7][2] + 0.5),
		0.,
		(4.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][0] - 4.5 * m_xsi[7][0] * m_xsi[7][0] + m_xsi[7][0]) * 0.5 },
	  { (27. * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] + 54. * m_xsi[7][0] * m_xsi[7][1] + 40.5 * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] - 45. * m_xsi[7][1] + 9.) * (0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 22.5 * m_xsi[7][0] * m_xsi[7][1] - 22.5 * m_xsi[7][1] * m_xsi[7][1] + 9. * m_xsi[7][1]) * 0.5 },
	  { (-54. * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][1] * m_xsi[7][1] + 27. * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
		(-27. * m_xsi[7][0] * m_xsi[7][0] - 54. * m_xsi[7][0] * m_xsi[7][1] + 27. * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
		(-27. * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 27. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] + 27. * m_xsi[7][0] * m_xsi[7][1]) * 0.5 },
	  { (27. * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] - 4.5 * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * 0.5 },
	  { (-13.5 * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
		(-27. * m_xsi[7][0] * m_xsi[7][1] - 40.5 * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][0] + 36. * m_xsi[7][1] - 4.5) * (0.5 * m_xsi[7][2] + 0.5),
		(-13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 13.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] + 4.5 * m_xsi[7][0] * m_xsi[7][1] + 18. * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * 0.5 },
	  { (13.5 * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1]) * (0.5 * m_xsi[7][2] + 0.5),
		(27. * m_xsi[7][0] * m_xsi[7][1] - 4.5 * m_xsi[7][0]) * (0.5 * m_xsi[7][2] + 0.5),
		(13.5 * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][0] * m_xsi[7][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[7][1] * m_xsi[7][1] - 9. * m_xsi[7][1] + 1.) * (0.5 * m_xsi[7][2] + 0.5),
		(4.5 * m_xsi[7][1] * m_xsi[7][1] * m_xsi[7][1] - 4.5 * m_xsi[7][1] * m_xsi[7][1] + m_xsi[7][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[8][0] * m_xsi[8][0] - 27. * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] + 18. * m_xsi[8][0] + 18. * m_xsi[8][1] - 5.5) * (-0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][0] - 27. * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] + 18. * m_xsi[8][0] + 18. * m_xsi[8][1] - 5.5) * (-0.5 * m_xsi[8][2] + 0.5),
		(-4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][0] * m_xsi[8][0] + 18. * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][1] * m_xsi[8][1] - 5.5 * m_xsi[8][0] - 5.5 * m_xsi[8][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[8][0] * m_xsi[8][0] + 54. * m_xsi[8][0] * m_xsi[8][1] + 13.5 * m_xsi[8][1] * m_xsi[8][1] - 45. * m_xsi[8][0] - 22.5 * m_xsi[8][1] + 9.) * (-0.5 * m_xsi[8][2] + 0.5),
		(27. * m_xsi[8][0] * m_xsi[8][0] + 27. * m_xsi[8][0] * m_xsi[8][1] - 22.5 * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] + 27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][0] - 22.5 * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][0]) * -0.5 },
	  { (-40.5 * m_xsi[8][0] * m_xsi[8][0] - 27. * m_xsi[8][0] * m_xsi[8][1] + 36. * m_xsi[8][0] + 4.5 * m_xsi[8][1] - 4.5) * (-0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][0] + 4.5 * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][0] * m_xsi[8][0] + 4.5 * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0]) * -0.5 },
	  { (13.5 * m_xsi[8][0] * m_xsi[8][0] - 9. * m_xsi[8][0] + 1.) * (-0.5 * m_xsi[8][2] + 0.5),
		0.,
		(4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 4.5 * m_xsi[8][0] * m_xsi[8][0] + m_xsi[8][0]) * -0.5 },
	  { (27. * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] + 54. * m_xsi[8][0] * m_xsi[8][1] + 40.5 * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] - 45. * m_xsi[8][1] + 9.) * (-0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][1] - 22.5 * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][1]) * -0.5 },
	  { (-54. * m_xsi[8][0] * m_xsi[8][1] - 27. * m_xsi[8][1] * m_xsi[8][1] + 27. * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
		(-27. * m_xsi[8][0] * m_xsi[8][0] - 54. * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
		(-27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1]) * -0.5 },
	  { (27. * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] - 4.5 * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * -0.5 },
	  { (-13.5 * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
		(-27. * m_xsi[8][0] * m_xsi[8][1] - 40.5 * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][0] + 36. * m_xsi[8][1] - 4.5) * (-0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * -0.5 },
	  { (13.5 * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * (-0.5 * m_xsi[8][2] + 0.5),
		(27. * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0]) * (-0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[8][1] * m_xsi[8][1] - 9. * m_xsi[8][1] + 1.) * (-0.5 * m_xsi[8][2] + 0.5),
		(4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] + m_xsi[8][1]) * -0.5 },
	  { (-13.5 * m_xsi[8][0] * m_xsi[8][0] - 27. * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] + 18. * m_xsi[8][0] + 18. * m_xsi[8][1] - 5.5) * (0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][0] - 27. * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] + 18. * m_xsi[8][0] + 18. * m_xsi[8][1] - 5.5) * (0.5 * m_xsi[8][2] + 0.5),
		(-4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][0] * m_xsi[8][0] + 18. * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][1] * m_xsi[8][1] - 5.5 * m_xsi[8][0] - 5.5 * m_xsi[8][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[8][0] * m_xsi[8][0] + 54. * m_xsi[8][0] * m_xsi[8][1] + 13.5 * m_xsi[8][1] * m_xsi[8][1] - 45. * m_xsi[8][0] - 22.5 * m_xsi[8][1] + 9.) * (0.5 * m_xsi[8][2] + 0.5),
		(27. * m_xsi[8][0] * m_xsi[8][0] + 27. * m_xsi[8][0] * m_xsi[8][1] - 22.5 * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] + 27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][0] - 22.5 * m_xsi[8][0] * m_xsi[8][1] + 9. * m_xsi[8][0]) * 0.5 },
	  { (-40.5 * m_xsi[8][0] * m_xsi[8][0] - 27. * m_xsi[8][0] * m_xsi[8][1] + 36. * m_xsi[8][0] + 4.5 * m_xsi[8][1] - 4.5) * (0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][0] + 4.5 * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][0] * m_xsi[8][0] + 4.5 * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0]) * 0.5 },
	  { (13.5 * m_xsi[8][0] * m_xsi[8][0] - 9. * m_xsi[8][0] + 1.) * (0.5 * m_xsi[8][2] + 0.5),
		0.,
		(4.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][0] - 4.5 * m_xsi[8][0] * m_xsi[8][0] + m_xsi[8][0]) * 0.5 },
	  { (27. * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] + 54. * m_xsi[8][0] * m_xsi[8][1] + 40.5 * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] - 45. * m_xsi[8][1] + 9.) * (0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 22.5 * m_xsi[8][0] * m_xsi[8][1] - 22.5 * m_xsi[8][1] * m_xsi[8][1] + 9. * m_xsi[8][1]) * 0.5 },
	  { (-54. * m_xsi[8][0] * m_xsi[8][1] - 27. * m_xsi[8][1] * m_xsi[8][1] + 27. * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
		(-27. * m_xsi[8][0] * m_xsi[8][0] - 54. * m_xsi[8][0] * m_xsi[8][1] + 27. * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
		(-27. * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 27. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] + 27. * m_xsi[8][0] * m_xsi[8][1]) * 0.5 },
	  { (27. * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] - 4.5 * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * 0.5 },
	  { (-13.5 * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
		(-27. * m_xsi[8][0] * m_xsi[8][1] - 40.5 * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][0] + 36. * m_xsi[8][1] - 4.5) * (0.5 * m_xsi[8][2] + 0.5),
		(-13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 13.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] + 4.5 * m_xsi[8][0] * m_xsi[8][1] + 18. * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * 0.5 },
	  { (13.5 * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1]) * (0.5 * m_xsi[8][2] + 0.5),
		(27. * m_xsi[8][0] * m_xsi[8][1] - 4.5 * m_xsi[8][0]) * (0.5 * m_xsi[8][2] + 0.5),
		(13.5 * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][0] * m_xsi[8][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[8][1] * m_xsi[8][1] - 9. * m_xsi[8][1] + 1.) * (0.5 * m_xsi[8][2] + 0.5),
		(4.5 * m_xsi[8][1] * m_xsi[8][1] * m_xsi[8][1] - 4.5 * m_xsi[8][1] * m_xsi[8][1] + m_xsi[8][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[9][0] * m_xsi[9][0] - 27. * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] + 18. * m_xsi[9][0] + 18. * m_xsi[9][1] - 5.5) * (-0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][0] - 27. * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] + 18. * m_xsi[9][0] + 18. * m_xsi[9][1] - 5.5) * (-0.5 * m_xsi[9][2] + 0.5),
		(-4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][0] * m_xsi[9][0] + 18. * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][1] * m_xsi[9][1] - 5.5 * m_xsi[9][0] - 5.5 * m_xsi[9][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[9][0] * m_xsi[9][0] + 54. * m_xsi[9][0] * m_xsi[9][1] + 13.5 * m_xsi[9][1] * m_xsi[9][1] - 45. * m_xsi[9][0] - 22.5 * m_xsi[9][1] + 9.) * (-0.5 * m_xsi[9][2] + 0.5),
		(27. * m_xsi[9][0] * m_xsi[9][0] + 27. * m_xsi[9][0] * m_xsi[9][1] - 22.5 * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] + 27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][0] - 22.5 * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][0]) * -0.5 },
	  { (-40.5 * m_xsi[9][0] * m_xsi[9][0] - 27. * m_xsi[9][0] * m_xsi[9][1] + 36. * m_xsi[9][0] + 4.5 * m_xsi[9][1] - 4.5) * (-0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][0] + 4.5 * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][0] * m_xsi[9][0] + 4.5 * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0]) * -0.5 },
	  { (13.5 * m_xsi[9][0] * m_xsi[9][0] - 9. * m_xsi[9][0] + 1.) * (-0.5 * m_xsi[9][2] + 0.5),
		0.,
		(4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 4.5 * m_xsi[9][0] * m_xsi[9][0] + m_xsi[9][0]) * -0.5 },
	  { (27. * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] + 54. * m_xsi[9][0] * m_xsi[9][1] + 40.5 * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] - 45. * m_xsi[9][1] + 9.) * (-0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][1] - 22.5 * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][1]) * -0.5 },
	  { (-54. * m_xsi[9][0] * m_xsi[9][1] - 27. * m_xsi[9][1] * m_xsi[9][1] + 27. * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
		(-27. * m_xsi[9][0] * m_xsi[9][0] - 54. * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
		(-27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1]) * -0.5 },
	  { (27. * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] - 4.5 * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * -0.5 },
	  { (-13.5 * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
		(-27. * m_xsi[9][0] * m_xsi[9][1] - 40.5 * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][0] + 36. * m_xsi[9][1] - 4.5) * (-0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * -0.5 },
	  { (13.5 * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * (-0.5 * m_xsi[9][2] + 0.5),
		(27. * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0]) * (-0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[9][1] * m_xsi[9][1] - 9. * m_xsi[9][1] + 1.) * (-0.5 * m_xsi[9][2] + 0.5),
		(4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] + m_xsi[9][1]) * -0.5 },
	  { (-13.5 * m_xsi[9][0] * m_xsi[9][0] - 27. * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] + 18. * m_xsi[9][0] + 18. * m_xsi[9][1] - 5.5) * (0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][0] - 27. * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] + 18. * m_xsi[9][0] + 18. * m_xsi[9][1] - 5.5) * (0.5 * m_xsi[9][2] + 0.5),
		(-4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][0] * m_xsi[9][0] + 18. * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][1] * m_xsi[9][1] - 5.5 * m_xsi[9][0] - 5.5 * m_xsi[9][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[9][0] * m_xsi[9][0] + 54. * m_xsi[9][0] * m_xsi[9][1] + 13.5 * m_xsi[9][1] * m_xsi[9][1] - 45. * m_xsi[9][0] - 22.5 * m_xsi[9][1] + 9.) * (0.5 * m_xsi[9][2] + 0.5),
		(27. * m_xsi[9][0] * m_xsi[9][0] + 27. * m_xsi[9][0] * m_xsi[9][1] - 22.5 * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] + 27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][0] - 22.5 * m_xsi[9][0] * m_xsi[9][1] + 9. * m_xsi[9][0]) * 0.5 },
	  { (-40.5 * m_xsi[9][0] * m_xsi[9][0] - 27. * m_xsi[9][0] * m_xsi[9][1] + 36. * m_xsi[9][0] + 4.5 * m_xsi[9][1] - 4.5) * (0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][0] + 4.5 * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][0] * m_xsi[9][0] + 4.5 * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0]) * 0.5 },
	  { (13.5 * m_xsi[9][0] * m_xsi[9][0] - 9. * m_xsi[9][0] + 1.) * (0.5 * m_xsi[9][2] + 0.5),
		0.,
		(4.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][0] - 4.5 * m_xsi[9][0] * m_xsi[9][0] + m_xsi[9][0]) * 0.5 },
	  { (27. * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] + 54. * m_xsi[9][0] * m_xsi[9][1] + 40.5 * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] - 45. * m_xsi[9][1] + 9.) * (0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 22.5 * m_xsi[9][0] * m_xsi[9][1] - 22.5 * m_xsi[9][1] * m_xsi[9][1] + 9. * m_xsi[9][1]) * 0.5 },
	  { (-54. * m_xsi[9][0] * m_xsi[9][1] - 27. * m_xsi[9][1] * m_xsi[9][1] + 27. * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
		(-27. * m_xsi[9][0] * m_xsi[9][0] - 54. * m_xsi[9][0] * m_xsi[9][1] + 27. * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
		(-27. * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 27. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] + 27. * m_xsi[9][0] * m_xsi[9][1]) * 0.5 },
	  { (27. * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] - 4.5 * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * 0.5 },
	  { (-13.5 * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
		(-27. * m_xsi[9][0] * m_xsi[9][1] - 40.5 * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][0] + 36. * m_xsi[9][1] - 4.5) * (0.5 * m_xsi[9][2] + 0.5),
		(-13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 13.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] + 4.5 * m_xsi[9][0] * m_xsi[9][1] + 18. * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * 0.5 },
	  { (13.5 * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1]) * (0.5 * m_xsi[9][2] + 0.5),
		(27. * m_xsi[9][0] * m_xsi[9][1] - 4.5 * m_xsi[9][0]) * (0.5 * m_xsi[9][2] + 0.5),
		(13.5 * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][0] * m_xsi[9][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[9][1] * m_xsi[9][1] - 9. * m_xsi[9][1] + 1.) * (0.5 * m_xsi[9][2] + 0.5),
		(4.5 * m_xsi[9][1] * m_xsi[9][1] * m_xsi[9][1] - 4.5 * m_xsi[9][1] * m_xsi[9][1] + m_xsi[9][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[10][0] * m_xsi[10][0] - 27. * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] + 18. * m_xsi[10][0] + 18. * m_xsi[10][1] - 5.5) * (-0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][0] - 27. * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] + 18. * m_xsi[10][0] + 18. * m_xsi[10][1] - 5.5) * (-0.5 * m_xsi[10][2] + 0.5),
		(-4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][0] * m_xsi[10][0] + 18. * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][1] * m_xsi[10][1] - 5.5 * m_xsi[10][0] - 5.5 * m_xsi[10][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[10][0] * m_xsi[10][0] + 54. * m_xsi[10][0] * m_xsi[10][1] + 13.5 * m_xsi[10][1] * m_xsi[10][1] - 45. * m_xsi[10][0] - 22.5 * m_xsi[10][1] + 9.) * (-0.5 * m_xsi[10][2] + 0.5),
		(27. * m_xsi[10][0] * m_xsi[10][0] + 27. * m_xsi[10][0] * m_xsi[10][1] - 22.5 * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] + 27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][0] - 22.5 * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][0]) * -0.5 },
	  { (-40.5 * m_xsi[10][0] * m_xsi[10][0] - 27. * m_xsi[10][0] * m_xsi[10][1] + 36. * m_xsi[10][0] + 4.5 * m_xsi[10][1] - 4.5) * (-0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][0] + 4.5 * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][0] * m_xsi[10][0] + 4.5 * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0]) * -0.5 },
	  { (13.5 * m_xsi[10][0] * m_xsi[10][0] - 9. * m_xsi[10][0] + 1.) * (-0.5 * m_xsi[10][2] + 0.5),
		0.,
		(4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 4.5 * m_xsi[10][0] * m_xsi[10][0] + m_xsi[10][0]) * -0.5 },
	  { (27. * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] + 54. * m_xsi[10][0] * m_xsi[10][1] + 40.5 * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] - 45. * m_xsi[10][1] + 9.) * (-0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][1] - 22.5 * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][1]) * -0.5 },
	  { (-54. * m_xsi[10][0] * m_xsi[10][1] - 27. * m_xsi[10][1] * m_xsi[10][1] + 27. * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
		(-27. * m_xsi[10][0] * m_xsi[10][0] - 54. * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
		(-27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1]) * -0.5 },
	  { (27. * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] - 4.5 * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * -0.5 },
	  { (-13.5 * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
		(-27. * m_xsi[10][0] * m_xsi[10][1] - 40.5 * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][0] + 36. * m_xsi[10][1] - 4.5) * (-0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * -0.5 },
	  { (13.5 * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * (-0.5 * m_xsi[10][2] + 0.5),
		(27. * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0]) * (-0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[10][1] * m_xsi[10][1] - 9. * m_xsi[10][1] + 1.) * (-0.5 * m_xsi[10][2] + 0.5),
		(4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] + m_xsi[10][1]) * -0.5 },
	  { (-13.5 * m_xsi[10][0] * m_xsi[10][0] - 27. * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] + 18. * m_xsi[10][0] + 18. * m_xsi[10][1] - 5.5) * (0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][0] - 27. * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] + 18. * m_xsi[10][0] + 18. * m_xsi[10][1] - 5.5) * (0.5 * m_xsi[10][2] + 0.5),
		(-4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][0] * m_xsi[10][0] + 18. * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][1] * m_xsi[10][1] - 5.5 * m_xsi[10][0] - 5.5 * m_xsi[10][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[10][0] * m_xsi[10][0] + 54. * m_xsi[10][0] * m_xsi[10][1] + 13.5 * m_xsi[10][1] * m_xsi[10][1] - 45. * m_xsi[10][0] - 22.5 * m_xsi[10][1] + 9.) * (0.5 * m_xsi[10][2] + 0.5),
		(27. * m_xsi[10][0] * m_xsi[10][0] + 27. * m_xsi[10][0] * m_xsi[10][1] - 22.5 * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] + 27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][0] - 22.5 * m_xsi[10][0] * m_xsi[10][1] + 9. * m_xsi[10][0]) * 0.5 },
	  { (-40.5 * m_xsi[10][0] * m_xsi[10][0] - 27. * m_xsi[10][0] * m_xsi[10][1] + 36. * m_xsi[10][0] + 4.5 * m_xsi[10][1] - 4.5) * (0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][0] + 4.5 * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][0] * m_xsi[10][0] + 4.5 * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0]) * 0.5 },
	  { (13.5 * m_xsi[10][0] * m_xsi[10][0] - 9. * m_xsi[10][0] + 1.) * (0.5 * m_xsi[10][2] + 0.5),
		0.,
		(4.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][0] - 4.5 * m_xsi[10][0] * m_xsi[10][0] + m_xsi[10][0]) * 0.5 },
	  { (27. * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] + 54. * m_xsi[10][0] * m_xsi[10][1] + 40.5 * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] - 45. * m_xsi[10][1] + 9.) * (0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 22.5 * m_xsi[10][0] * m_xsi[10][1] - 22.5 * m_xsi[10][1] * m_xsi[10][1] + 9. * m_xsi[10][1]) * 0.5 },
	  { (-54. * m_xsi[10][0] * m_xsi[10][1] - 27. * m_xsi[10][1] * m_xsi[10][1] + 27. * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
		(-27. * m_xsi[10][0] * m_xsi[10][0] - 54. * m_xsi[10][0] * m_xsi[10][1] + 27. * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
		(-27. * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 27. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] + 27. * m_xsi[10][0] * m_xsi[10][1]) * 0.5 },
	  { (27. * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] - 4.5 * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * 0.5 },
	  { (-13.5 * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
		(-27. * m_xsi[10][0] * m_xsi[10][1] - 40.5 * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][0] + 36. * m_xsi[10][1] - 4.5) * (0.5 * m_xsi[10][2] + 0.5),
		(-13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 13.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] + 4.5 * m_xsi[10][0] * m_xsi[10][1] + 18. * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * 0.5 },
	  { (13.5 * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1]) * (0.5 * m_xsi[10][2] + 0.5),
		(27. * m_xsi[10][0] * m_xsi[10][1] - 4.5 * m_xsi[10][0]) * (0.5 * m_xsi[10][2] + 0.5),
		(13.5 * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][0] * m_xsi[10][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[10][1] * m_xsi[10][1] - 9. * m_xsi[10][1] + 1.) * (0.5 * m_xsi[10][2] + 0.5),
		(4.5 * m_xsi[10][1] * m_xsi[10][1] * m_xsi[10][1] - 4.5 * m_xsi[10][1] * m_xsi[10][1] + m_xsi[10][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[11][0] * m_xsi[11][0] - 27. * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] + 18. * m_xsi[11][0] + 18. * m_xsi[11][1] - 5.5) * (-0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][0] - 27. * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] + 18. * m_xsi[11][0] + 18. * m_xsi[11][1] - 5.5) * (-0.5 * m_xsi[11][2] + 0.5),
		(-4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][0] * m_xsi[11][0] + 18. * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][1] * m_xsi[11][1] - 5.5 * m_xsi[11][0] - 5.5 * m_xsi[11][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[11][0] * m_xsi[11][0] + 54. * m_xsi[11][0] * m_xsi[11][1] + 13.5 * m_xsi[11][1] * m_xsi[11][1] - 45. * m_xsi[11][0] - 22.5 * m_xsi[11][1] + 9.) * (-0.5 * m_xsi[11][2] + 0.5),
		(27. * m_xsi[11][0] * m_xsi[11][0] + 27. * m_xsi[11][0] * m_xsi[11][1] - 22.5 * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] + 27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][0] - 22.5 * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][0]) * -0.5 },
	  { (-40.5 * m_xsi[11][0] * m_xsi[11][0] - 27. * m_xsi[11][0] * m_xsi[11][1] + 36. * m_xsi[11][0] + 4.5 * m_xsi[11][1] - 4.5) * (-0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][0] + 4.5 * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][0] * m_xsi[11][0] + 4.5 * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0]) * -0.5 },
	  { (13.5 * m_xsi[11][0] * m_xsi[11][0] - 9. * m_xsi[11][0] + 1.) * (-0.5 * m_xsi[11][2] + 0.5),
		0.,
		(4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 4.5 * m_xsi[11][0] * m_xsi[11][0] + m_xsi[11][0]) * -0.5 },
	  { (27. * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] + 54. * m_xsi[11][0] * m_xsi[11][1] + 40.5 * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] - 45. * m_xsi[11][1] + 9.) * (-0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][1] - 22.5 * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][1]) * -0.5 },
	  { (-54. * m_xsi[11][0] * m_xsi[11][1] - 27. * m_xsi[11][1] * m_xsi[11][1] + 27. * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
		(-27. * m_xsi[11][0] * m_xsi[11][0] - 54. * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
		(-27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1]) * -0.5 },
	  { (27. * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] - 4.5 * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * -0.5 },
	  { (-13.5 * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
		(-27. * m_xsi[11][0] * m_xsi[11][1] - 40.5 * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][0] + 36. * m_xsi[11][1] - 4.5) * (-0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * -0.5 },
	  { (13.5 * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * (-0.5 * m_xsi[11][2] + 0.5),
		(27. * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0]) * (-0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[11][1] * m_xsi[11][1] - 9. * m_xsi[11][1] + 1.) * (-0.5 * m_xsi[11][2] + 0.5),
		(4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] + m_xsi[11][1]) * -0.5 },
	  { (-13.5 * m_xsi[11][0] * m_xsi[11][0] - 27. * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] + 18. * m_xsi[11][0] + 18. * m_xsi[11][1] - 5.5) * (0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][0] - 27. * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] + 18. * m_xsi[11][0] + 18. * m_xsi[11][1] - 5.5) * (0.5 * m_xsi[11][2] + 0.5),
		(-4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][0] * m_xsi[11][0] + 18. * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][1] * m_xsi[11][1] - 5.5 * m_xsi[11][0] - 5.5 * m_xsi[11][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[11][0] * m_xsi[11][0] + 54. * m_xsi[11][0] * m_xsi[11][1] + 13.5 * m_xsi[11][1] * m_xsi[11][1] - 45. * m_xsi[11][0] - 22.5 * m_xsi[11][1] + 9.) * (0.5 * m_xsi[11][2] + 0.5),
		(27. * m_xsi[11][0] * m_xsi[11][0] + 27. * m_xsi[11][0] * m_xsi[11][1] - 22.5 * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] + 27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][0] - 22.5 * m_xsi[11][0] * m_xsi[11][1] + 9. * m_xsi[11][0]) * 0.5 },
	  { (-40.5 * m_xsi[11][0] * m_xsi[11][0] - 27. * m_xsi[11][0] * m_xsi[11][1] + 36. * m_xsi[11][0] + 4.5 * m_xsi[11][1] - 4.5) * (0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][0] + 4.5 * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][0] * m_xsi[11][0] + 4.5 * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0]) * 0.5 },
	  { (13.5 * m_xsi[11][0] * m_xsi[11][0] - 9. * m_xsi[11][0] + 1.) * (0.5 * m_xsi[11][2] + 0.5),
		0.,
		(4.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][0] - 4.5 * m_xsi[11][0] * m_xsi[11][0] + m_xsi[11][0]) * 0.5 },
	  { (27. * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] + 54. * m_xsi[11][0] * m_xsi[11][1] + 40.5 * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] - 45. * m_xsi[11][1] + 9.) * (0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 22.5 * m_xsi[11][0] * m_xsi[11][1] - 22.5 * m_xsi[11][1] * m_xsi[11][1] + 9. * m_xsi[11][1]) * 0.5 },
	  { (-54. * m_xsi[11][0] * m_xsi[11][1] - 27. * m_xsi[11][1] * m_xsi[11][1] + 27. * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
		(-27. * m_xsi[11][0] * m_xsi[11][0] - 54. * m_xsi[11][0] * m_xsi[11][1] + 27. * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
		(-27. * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 27. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] + 27. * m_xsi[11][0] * m_xsi[11][1]) * 0.5 },
	  { (27. * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] - 4.5 * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * 0.5 },
	  { (-13.5 * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
		(-27. * m_xsi[11][0] * m_xsi[11][1] - 40.5 * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][0] + 36. * m_xsi[11][1] - 4.5) * (0.5 * m_xsi[11][2] + 0.5),
		(-13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 13.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] + 4.5 * m_xsi[11][0] * m_xsi[11][1] + 18. * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * 0.5 },
	  { (13.5 * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1]) * (0.5 * m_xsi[11][2] + 0.5),
		(27. * m_xsi[11][0] * m_xsi[11][1] - 4.5 * m_xsi[11][0]) * (0.5 * m_xsi[11][2] + 0.5),
		(13.5 * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][0] * m_xsi[11][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[11][1] * m_xsi[11][1] - 9. * m_xsi[11][1] + 1.) * (0.5 * m_xsi[11][2] + 0.5),
		(4.5 * m_xsi[11][1] * m_xsi[11][1] * m_xsi[11][1] - 4.5 * m_xsi[11][1] * m_xsi[11][1] + m_xsi[11][1]) * 0.5 } },

	{ { (-13.5 * m_xsi[12][0] * m_xsi[12][0] - 27. * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] + 18. * m_xsi[12][0] + 18. * m_xsi[12][1] - 5.5) * (-0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][0] - 27. * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] + 18. * m_xsi[12][0] + 18. * m_xsi[12][1] - 5.5) * (-0.5 * m_xsi[12][2] + 0.5),
		(-4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][0] * m_xsi[12][0] + 18. * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][1] * m_xsi[12][1] - 5.5 * m_xsi[12][0] - 5.5 * m_xsi[12][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[12][0] * m_xsi[12][0] + 54. * m_xsi[12][0] * m_xsi[12][1] + 13.5 * m_xsi[12][1] * m_xsi[12][1] - 45. * m_xsi[12][0] - 22.5 * m_xsi[12][1] + 9.) * (-0.5 * m_xsi[12][2] + 0.5),
		(27. * m_xsi[12][0] * m_xsi[12][0] + 27. * m_xsi[12][0] * m_xsi[12][1] - 22.5 * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] + 27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][0] - 22.5 * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][0]) * -0.5 },
	  { (-40.5 * m_xsi[12][0] * m_xsi[12][0] - 27. * m_xsi[12][0] * m_xsi[12][1] + 36. * m_xsi[12][0] + 4.5 * m_xsi[12][1] - 4.5) * (-0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][0] + 4.5 * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][0] * m_xsi[12][0] + 4.5 * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0]) * -0.5 },
	  { (13.5 * m_xsi[12][0] * m_xsi[12][0] - 9. * m_xsi[12][0] + 1.) * (-0.5 * m_xsi[12][2] + 0.5),
		0.,
		(4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 4.5 * m_xsi[12][0] * m_xsi[12][0] + m_xsi[12][0]) * -0.5 },
	  { (27. * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] + 54. * m_xsi[12][0] * m_xsi[12][1] + 40.5 * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] - 45. * m_xsi[12][1] + 9.) * (-0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][1] - 22.5 * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][1]) * -0.5 },
	  { (-54. * m_xsi[12][0] * m_xsi[12][1] - 27. * m_xsi[12][1] * m_xsi[12][1] + 27. * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
		(-27. * m_xsi[12][0] * m_xsi[12][0] - 54. * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
		(-27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1]) * -0.5 },
	  { (27. * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] - 4.5 * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * -0.5 },
	  { (-13.5 * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
		(-27. * m_xsi[12][0] * m_xsi[12][1] - 40.5 * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][0] + 36. * m_xsi[12][1] - 4.5) * (-0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * -0.5 },
	  { (13.5 * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * (-0.5 * m_xsi[12][2] + 0.5),
		(27. * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0]) * (-0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[12][1] * m_xsi[12][1] - 9. * m_xsi[12][1] + 1.) * (-0.5 * m_xsi[12][2] + 0.5),
		(4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] + m_xsi[12][1]) * -0.5 },
	  { (-13.5 * m_xsi[12][0] * m_xsi[12][0] - 27. * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] + 18. * m_xsi[12][0] + 18. * m_xsi[12][1] - 5.5) * (0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][0] - 27. * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] + 18. * m_xsi[12][0] + 18. * m_xsi[12][1] - 5.5) * (0.5 * m_xsi[12][2] + 0.5),
		(-4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][0] * m_xsi[12][0] + 18. * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][1] * m_xsi[12][1] - 5.5 * m_xsi[12][0] - 5.5 * m_xsi[12][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[12][0] * m_xsi[12][0] + 54. * m_xsi[12][0] * m_xsi[12][1] + 13.5 * m_xsi[12][1] * m_xsi[12][1] - 45. * m_xsi[12][0] - 22.5 * m_xsi[12][1] + 9.) * (0.5 * m_xsi[12][2] + 0.5),
		(27. * m_xsi[12][0] * m_xsi[12][0] + 27. * m_xsi[12][0] * m_xsi[12][1] - 22.5 * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] + 27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][0] - 22.5 * m_xsi[12][0] * m_xsi[12][1] + 9. * m_xsi[12][0]) * 0.5 },
	  { (-40.5 * m_xsi[12][0] * m_xsi[12][0] - 27. * m_xsi[12][0] * m_xsi[12][1] + 36. * m_xsi[12][0] + 4.5 * m_xsi[12][1] - 4.5) * (0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][0] + 4.5 * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][0] * m_xsi[12][0] + 4.5 * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0]) * 0.5 },
	  { (13.5 * m_xsi[12][0] * m_xsi[12][0] - 9. * m_xsi[12][0] + 1.) * (0.5 * m_xsi[12][2] + 0.5),
		0.,
		(4.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][0] - 4.5 * m_xsi[12][0] * m_xsi[12][0] + m_xsi[12][0]) * 0.5 },
	  { (27. * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] + 54. * m_xsi[12][0] * m_xsi[12][1] + 40.5 * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] - 45. * m_xsi[12][1] + 9.) * (0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 22.5 * m_xsi[12][0] * m_xsi[12][1] - 22.5 * m_xsi[12][1] * m_xsi[12][1] + 9. * m_xsi[12][1]) * 0.5 },
	  { (-54. * m_xsi[12][0] * m_xsi[12][1] - 27. * m_xsi[12][1] * m_xsi[12][1] + 27. * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
		(-27. * m_xsi[12][0] * m_xsi[12][0] - 54. * m_xsi[12][0] * m_xsi[12][1] + 27. * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
		(-27. * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 27. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] + 27. * m_xsi[12][0] * m_xsi[12][1]) * 0.5 },
	  { (27. * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] - 4.5 * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * 0.5 },
	  { (-13.5 * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
		(-27. * m_xsi[12][0] * m_xsi[12][1] - 40.5 * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][0] + 36. * m_xsi[12][1] - 4.5) * (0.5 * m_xsi[12][2] + 0.5),
		(-13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 13.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] + 4.5 * m_xsi[12][0] * m_xsi[12][1] + 18. * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * 0.5 },
	  { (13.5 * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1]) * (0.5 * m_xsi[12][2] + 0.5),
		(27. * m_xsi[12][0] * m_xsi[12][1] - 4.5 * m_xsi[12][0]) * (0.5 * m_xsi[12][2] + 0.5),
		(13.5 * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][0] * m_xsi[12][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[12][1] * m_xsi[12][1] - 9. * m_xsi[12][1] + 1.) * (0.5 * m_xsi[12][2] + 0.5),
		(4.5 * m_xsi[12][1] * m_xsi[12][1] * m_xsi[12][1] - 4.5 * m_xsi[12][1] * m_xsi[12][1] + m_xsi[12][1]) * 0.5 } },


	{ { (-13.5 * m_xsi[13][0] * m_xsi[13][0] - 27. * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] + 18. * m_xsi[13][0] + 18. * m_xsi[13][1] - 5.5) * (-0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][0] - 27. * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] + 18. * m_xsi[13][0] + 18. * m_xsi[13][1] - 5.5) * (-0.5 * m_xsi[13][2] + 0.5),
		(-4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][0] * m_xsi[13][0] + 18. * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][1] * m_xsi[13][1] - 5.5 * m_xsi[13][0] - 5.5 * m_xsi[13][1] + 1.) * -0.5 },
	  { (40.5 * m_xsi[13][0] * m_xsi[13][0] + 54. * m_xsi[13][0] * m_xsi[13][1] + 13.5 * m_xsi[13][1] * m_xsi[13][1] - 45. * m_xsi[13][0] - 22.5 * m_xsi[13][1] + 9.) * (-0.5 * m_xsi[13][2] + 0.5),
		(27. * m_xsi[13][0] * m_xsi[13][0] + 27. * m_xsi[13][0] * m_xsi[13][1] - 22.5 * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] + 27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][0] - 22.5 * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][0]) * -0.5 },
	  { (-40.5 * m_xsi[13][0] * m_xsi[13][0] - 27. * m_xsi[13][0] * m_xsi[13][1] + 36. * m_xsi[13][0] + 4.5 * m_xsi[13][1] - 4.5) * (-0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][0] + 4.5 * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][0] * m_xsi[13][0] + 4.5 * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0]) * -0.5 },
	  { (13.5 * m_xsi[13][0] * m_xsi[13][0] - 9. * m_xsi[13][0] + 1.) * (-0.5 * m_xsi[13][2] + 0.5),
		0.,
		(4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 4.5 * m_xsi[13][0] * m_xsi[13][0] + m_xsi[13][0]) * -0.5 },
	  { (27. * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] + 54. * m_xsi[13][0] * m_xsi[13][1] + 40.5 * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] - 45. * m_xsi[13][1] + 9.) * (-0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][1] - 22.5 * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][1]) * -0.5 },
	  { (-54. * m_xsi[13][0] * m_xsi[13][1] - 27. * m_xsi[13][1] * m_xsi[13][1] + 27. * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
		(-27. * m_xsi[13][0] * m_xsi[13][0] - 54. * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
		(-27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1]) * -0.5 },
	  { (27. * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] - 4.5 * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * -0.5 },
	  { (-13.5 * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
		(-27. * m_xsi[13][0] * m_xsi[13][1] - 40.5 * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][0] + 36. * m_xsi[13][1] - 4.5) * (-0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * -0.5 },
	  { (13.5 * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * (-0.5 * m_xsi[13][2] + 0.5),
		(27. * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0]) * (-0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * -0.5 },
	  { 0.,
		(13.5 * m_xsi[13][1] * m_xsi[13][1] - 9. * m_xsi[13][1] + 1.) * (-0.5 * m_xsi[13][2] + 0.5),
		(4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] + m_xsi[13][1]) * -0.5 },
	  { (-13.5 * m_xsi[13][0] * m_xsi[13][0] - 27. * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] + 18. * m_xsi[13][0] + 18. * m_xsi[13][1] - 5.5) * (0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][0] - 27. * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] + 18. * m_xsi[13][0] + 18. * m_xsi[13][1] - 5.5) * (0.5 * m_xsi[13][2] + 0.5),
		(-4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][0] * m_xsi[13][0] + 18. * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][1] * m_xsi[13][1] - 5.5 * m_xsi[13][0] - 5.5 * m_xsi[13][1] + 1.) * 0.5 },
	  { (40.5 * m_xsi[13][0] * m_xsi[13][0] + 54. * m_xsi[13][0] * m_xsi[13][1] + 13.5 * m_xsi[13][1] * m_xsi[13][1] - 45. * m_xsi[13][0] - 22.5 * m_xsi[13][1] + 9.) * (0.5 * m_xsi[13][2] + 0.5),
		(27. * m_xsi[13][0] * m_xsi[13][0] + 27. * m_xsi[13][0] * m_xsi[13][1] - 22.5 * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] + 27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][0] - 22.5 * m_xsi[13][0] * m_xsi[13][1] + 9. * m_xsi[13][0]) * 0.5 },
	  { (-40.5 * m_xsi[13][0] * m_xsi[13][0] - 27. * m_xsi[13][0] * m_xsi[13][1] + 36. * m_xsi[13][0] + 4.5 * m_xsi[13][1] - 4.5) * (0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][0] + 4.5 * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][0] * m_xsi[13][0] + 4.5 * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0]) * 0.5 },
	  { (13.5 * m_xsi[13][0] * m_xsi[13][0] - 9. * m_xsi[13][0] + 1.) * (0.5 * m_xsi[13][2] + 0.5),
		0.,
		(4.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][0] - 4.5 * m_xsi[13][0] * m_xsi[13][0] + m_xsi[13][0]) * 0.5 },
	  { (27. * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] + 54. * m_xsi[13][0] * m_xsi[13][1] + 40.5 * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] - 45. * m_xsi[13][1] + 9.) * (0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 22.5 * m_xsi[13][0] * m_xsi[13][1] - 22.5 * m_xsi[13][1] * m_xsi[13][1] + 9. * m_xsi[13][1]) * 0.5 },
	  { (-54. * m_xsi[13][0] * m_xsi[13][1] - 27. * m_xsi[13][1] * m_xsi[13][1] + 27. * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
		(-27. * m_xsi[13][0] * m_xsi[13][0] - 54. * m_xsi[13][0] * m_xsi[13][1] + 27. * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
		(-27. * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 27. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] + 27. * m_xsi[13][0] * m_xsi[13][1]) * 0.5 },
	  { (27. * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] - 4.5 * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * 0.5 },
	  { (-13.5 * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
		(-27. * m_xsi[13][0] * m_xsi[13][1] - 40.5 * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][0] + 36. * m_xsi[13][1] - 4.5) * (0.5 * m_xsi[13][2] + 0.5),
		(-13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 13.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] + 4.5 * m_xsi[13][0] * m_xsi[13][1] + 18. * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * 0.5 },
	  { (13.5 * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1]) * (0.5 * m_xsi[13][2] + 0.5),
		(27. * m_xsi[13][0] * m_xsi[13][1] - 4.5 * m_xsi[13][0]) * (0.5 * m_xsi[13][2] + 0.5),
		(13.5 * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][0] * m_xsi[13][1]) * 0.5 },
	  { 0.,
		(13.5 * m_xsi[13][1] * m_xsi[13][1] - 9. * m_xsi[13][1] + 1.) * (0.5 * m_xsi[13][2] + 0.5),
		(4.5 * m_xsi[13][1] * m_xsi[13][1] * m_xsi[13][1] - 4.5 * m_xsi[13][1] * m_xsi[13][1] + m_xsi[13][1]) * 0.5 } } };

