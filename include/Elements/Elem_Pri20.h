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
// Solid element, with cubic / linear interpolation functions, prism shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

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
	explicit Elem_Pri20(std::shared_ptr<Material>& Material)
		: ElementSolid(Material) { };

	// Output function for AcadView, based on element index.
	const std::string printByIndex_AV(const size_t add) const override {
		std::stringstream msg;

		msg << "2 3 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " "
			<< this->v_Conect[3]->m_index + add << " " << this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " "
			<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[7]->m_index + add << " " << this->v_Conect[8]->m_index + add << " "
			<< this->v_Conect[9]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "2 3 " << this->v_Conect[10]->m_index + add << " " << this->v_Conect[11]->m_index + add << " " << this->v_Conect[12]->m_index + add << " "
			<< this->v_Conect[13]->m_index + add << " " << this->v_Conect[14]->m_index + add << " " << this->v_Conect[15]->m_index + add << " "
			<< this->v_Conect[16]->m_index + add << " " << this->v_Conect[17]->m_index + add << " " << this->v_Conect[18]->m_index + add << " "
			<< this->v_Conect[19]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "3 1 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->v_Conect[10]->m_index + add << " " << this->v_Conect[13]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "3 1 " << this->v_Conect[3]->m_index + add << " " << this->v_Conect[9]->m_index + add << " " << this->v_Conect[13]->m_index + add << " " << this->v_Conect[19]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "3 1 " << this->v_Conect[9]->m_index + add << " " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[19]->m_index + add << " " << this->v_Conect[10]->m_index + add << " " << this->m_Mat->m_index << "\n";

		return msg.str();
	};

	// Output function for AcadView, based on element node number.
	const std::string printByAdder_AV(const size_t add) const override {
		std::stringstream msg;

		msg << "2 3 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (4 + add) << " " << (5 + add) << " "
			<< (6 + add) << " " << (7 + add) << " " << (8 + add) << " " << (9 + add) << " " << (10 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "2 3 " << (11 + add) << " " << (12 + add) << " " << (13 + add) << " " << (14 + add) << " " << (15 + add) << " "
			<< (16 + add) << " " << (17 + add) << " " << (18 + add) << " " << (19 + add) << " " << (20 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "3 1 " << (1 + add) << " " << (4 + add) << " " << (11 + add) << " " << (14 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "3 1 " << (4 + add) << " " << (10 + add) << " " << (14 + add) << " " << (20 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "3 1 " << (10 + add) << " " << (1 + add) << " " << (20 + add) << " " << (11 + add) << " " << this->m_Mat->m_index << "\n";

		return msg.str();
	};

	// Evaluates shape function in the point.
	Eigen::VectorXd getShapeFcOnPoint(const double* Point) override;

	// Evaluates the derivative of shape function in the point.
	Eigen::MatrixXd getShapeDerivOnPoint(const double* Point) override;

	// Returns a pointer to the first element of the shape functions (with size [nIP][m_NumNodes]).
	double const* getShapeFc() const override { return &m_Psi[0][0]; };

	// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][m_NumNodes][m_Dim]).
	double const* getShapeDerivative() const override { return &m_DPsi[0][0][0]; };

	// Returns a pointer to the weight of the integation points (with size [nIP]).
	double const* getWeight() const override { return &m_weight[0]; };

	// Returns the number of nodes of current element.
	int getNumNodes() override { return m_NumNodes; };

	// Returns the number of faces of current element.
	int getNumFaces() override { return m_NumFaces; };

	// Returns the number of integration points of current element.
	int getNumIP() override { return m_NumIP; };

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

		if (*std::max_element(new_xsi.begin(), new_xsi.end()) < 1.000001 && *std::min_element(new_xsi.begin(), new_xsi.end()-1) > -0.000001 && new_xsi.at(m_Dim) > -1.000001) return true;
		return false;
	};

private:
	// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
	void setGeomProperties() override;

private:
	/** @brief Number of Nodes */
	static const int m_NumNodes{ 20 };

	/** @brief Number of Integration Points */
	static const int m_NumIP{ 14 };

	/** @brief Number of Faces */
	static const int m_NumFaces{ 5 };

	/** @brief Weights for numerical integration */
	static const double m_weight[m_NumIP];

	/** @brief Shape functions */
	static const double m_Psi[m_NumIP][m_NumNodes];

	/** @brief Shape functions derivative */
	static const double m_DPsi[m_NumIP][m_NumNodes][m_Dim];

	// Since the number of Integration points is fixed, the shape functions are also constant (static)
	static const double m_xsi[m_NumIP][m_Dim];
};


// ================================================================================================
//
// Implementation of Member Function: getShapeFcOnPoint
// Shape functions evaluated on Point
// 
// ================================================================================================
inline Eigen::VectorXd Elem_Pri20::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(20);

	Psi(0) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * (-0.5 * Point[2] + 0.5);
	Psi(1) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * (-0.5 * Point[2] + 0.5);
	Psi(2) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	Psi(3) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * (-0.5 * Point[2] + 0.5);
	Psi(4) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * (-0.5 * Point[2] + 0.5);
	Psi(5) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * (-0.5 * Point[2] + 0.5);
	Psi(6) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * (-0.5 * Point[2] + 0.5);
	Psi(7) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	Psi(8) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * (-0.5 * Point[2] + 0.5);
	Psi(9) = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * (-0.5 * Point[2] + 0.5);
	Psi(10) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * (0.5 * Point[2] + 0.5);
	Psi(11) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * (0.5 * Point[2] + 0.5);
	Psi(12) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	Psi(13) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * (0.5 * Point[2] + 0.5);
	Psi(14) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * (0.5 * Point[2] + 0.5);
	Psi(15) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * (0.5 * Point[2] + 0.5);
	Psi(16) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * (0.5 * Point[2] + 0.5);
	Psi(17) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	Psi(18) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * (0.5 * Point[2] + 0.5);
	Psi(19) = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * (0.5 * Point[2] + 0.5);

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd Elem_Pri20::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(20, 3);

	DPsi(0, 0) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (-0.5 * Point[2] + 0.5);
	DPsi(1, 0) = (40.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 13.5 * Point[1] * Point[1] - 45. * Point[0] - 22.5 * Point[1] + 9.) * (-0.5 * Point[2] + 0.5);
	DPsi(2, 0) = (-40.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] + 36. * Point[0] + 4.5 * Point[1] - 4.5) * (-0.5 * Point[2] + 0.5);
	DPsi(3, 0) = (13.5 * Point[0] * Point[0] - 9. * Point[0] + 1.) * (-0.5 * Point[2] + 0.5);
	DPsi(4, 0) = (27. * Point[0] * Point[1] + 27. * Point[1] * Point[1] - 22.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	DPsi(5, 0) = (-54. * Point[0] * Point[1] - 27. * Point[1] * Point[1] + 27. * Point[1]) * (-0.5 * Point[2] + 0.5);
	DPsi(6, 0) = (27. * Point[0] * Point[1] - 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	DPsi(7, 0) = (-13.5 * Point[1] * Point[1] + 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	DPsi(8, 0) = (13.5 * Point[1] * Point[1] - 4.5 * Point[1]) * (-0.5 * Point[2] + 0.5);
	DPsi(9, 0) = 0.;
	DPsi(10, 0) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (0.5 * Point[2] + 0.5);
	DPsi(11, 0) = (40.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 13.5 * Point[1] * Point[1] - 45. * Point[0] - 22.5 * Point[1] + 9.) * (0.5 * Point[2] + 0.5);
	DPsi(12, 0) = (-40.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] + 36. * Point[0] + 4.5 * Point[1] - 4.5) * (0.5 * Point[2] + 0.5);
	DPsi(13, 0) = (13.5 * Point[0] * Point[0] - 9. * Point[0] + 1.) * (0.5 * Point[2] + 0.5);
	DPsi(14, 0) = (27. * Point[0] * Point[1] + 27. * Point[1] * Point[1] - 22.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	DPsi(15, 0) = (-54. * Point[0] * Point[1] - 27. * Point[1] * Point[1] + 27. * Point[1]) * (0.5 * Point[2] + 0.5);
	DPsi(16, 0) = (27. * Point[0] * Point[1] - 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	DPsi(17, 0) = (-13.5 * Point[1] * Point[1] + 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	DPsi(18, 0) = (13.5 * Point[1] * Point[1] - 4.5 * Point[1]) * (0.5 * Point[2] + 0.5);
	DPsi(19, 0) = 0.;

	DPsi(0,  1) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (-0.5 * Point[2] + 0.5);
	DPsi(1,  1) = (27. * Point[0] * Point[0] + 27. * Point[0] * Point[1] - 22.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	DPsi(2,  1) = (-13.5 * Point[0] * Point[0] + 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	DPsi(3,  1) = 0.;
	DPsi(4,  1) = (13.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 40.5 * Point[1] * Point[1] - 22.5 * Point[0] - 45. * Point[1] + 9.) * (-0.5 * Point[2] + 0.5);
	DPsi(5,  1) = (-27. * Point[0] * Point[0] - 54. * Point[0] * Point[1] + 27. * Point[0]) * (-0.5 * Point[2] + 0.5);
	DPsi(6,  1) = (13.5 * Point[0] * Point[0] - 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	DPsi(7,  1) = (-27. * Point[0] * Point[1] - 40.5 * Point[1] * Point[1] + 4.5 * Point[0] + 36. * Point[1] - 4.5) * (-0.5 * Point[2] + 0.5);
	DPsi(8,  1) = (27. * Point[0] * Point[1] - 4.5 * Point[0]) * (-0.5 * Point[2] + 0.5);
	DPsi(9,  1)  = (13.5 * Point[1] * Point[1] - 9. * Point[1] + 1.) * (-0.5 * Point[2] + 0.5);
	DPsi(10, 1) = (-13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1] + 18. * Point[0] + 18. * Point[1] - 5.5) * (0.5 * Point[2] + 0.5);
	DPsi(11, 1) = (27. * Point[0] * Point[0] + 27. * Point[0] * Point[1] - 22.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	DPsi(12, 1) = (-13.5 * Point[0] * Point[0] + 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	DPsi(13, 1) = 0.;
	DPsi(14, 1) = (13.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 40.5 * Point[1] * Point[1] - 22.5 * Point[0] - 45. * Point[1] + 9.) * (0.5 * Point[2] + 0.5);
	DPsi(15, 1) = (-27. * Point[0] * Point[0] - 54. * Point[0] * Point[1] + 27. * Point[0]) * (0.5 * Point[2] + 0.5);
	DPsi(16, 1) = (13.5 * Point[0] * Point[0] - 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	DPsi(17, 1) = (-27. * Point[0] * Point[1] - 40.5 * Point[1] * Point[1] + 4.5 * Point[0] + 36. * Point[1] - 4.5) * (0.5 * Point[2] + 0.5);
	DPsi(18, 1) = (27. * Point[0] * Point[1] - 4.5 * Point[0]) * (0.5 * Point[2] + 0.5);
	DPsi(19, 1) = (13.5 * Point[1] * Point[1] - 9. * Point[1] + 1.) * (0.5 * Point[2] + 0.5);

	DPsi(0,  2) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * -0.5;
	DPsi(1,  2) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * -0.5;
	DPsi(2,  2) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * -0.5;
	DPsi(3,  2) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * -0.5;
	DPsi(4,  2) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * -0.5;
	DPsi(5,  2) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * -0.5;
	DPsi(6,  2) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * -0.5;
	DPsi(7,  2) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * -0.5;
	DPsi(8,  2) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * -0.5;
	DPsi(9,  2)  = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * -0.5;
	DPsi(10, 2) = (-4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 5.5 * Point[0] - 5.5 * Point[1] + 1.) * 0.5;
	DPsi(11, 2) = (13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 9. * Point[0]) * 0.5;
	DPsi(12, 2) = (-13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 4.5 * Point[0]) * 0.5;
	DPsi(13, 2) = (4.5 * Point[0] * Point[0] * Point[0] - 4.5 * Point[0] * Point[0] + Point[0]) * 0.5;
	DPsi(14, 2) = (13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 9. * Point[1]) * 0.5;
	DPsi(15, 2) = (-27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1] + 27. * Point[0] * Point[1]) * 0.5;
	DPsi(16, 2) = (13.5 * Point[0] * Point[0] * Point[1] - 4.5 * Point[0] * Point[1]) * 0.5;
	DPsi(17, 2) = (-13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] + 18. * Point[1] * Point[1] - 4.5 * Point[1]) * 0.5;
	DPsi(18, 2) = (13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[0] * Point[1]) * 0.5;
	DPsi(19, 2) = (4.5 * Point[1] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] + Point[1]) * 0.5;

	return DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void Elem_Pri20::setGeomProperties() {

	const int nVertices = 6;

	// Allocate an array with size m_Dim to which m_Centroid points to.
	m_Centroid = std::make_unique<double[]>(m_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<Node<m_Dim>*, nVertices> vertices;
	vertices[0] = v_Conect[0].get();
	vertices[1] = v_Conect[3].get();
	vertices[2] = v_Conect[9].get();
	vertices[3] = v_Conect[10].get();
	vertices[4] = v_Conect[13].get();
	vertices[5] = v_Conect[19].get();

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
// Integration Points
//
// ================================================================================================
inline const double Elem_Pri20::m_xsi[m_NumIP][m_Dim] =
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
inline const double Elem_Pri20::m_weight[m_NumIP] = { 0.1125, 0.0629695902724135, 0.0629695902724135, 0.0629695902724135, 0.066197076394253, 0.066197076394253, 0.066197076394253, 0.1125, 0.0629695902724135, 0.0629695902724135, 0.0629695902724135, 0.066197076394253, 0.066197076394253, 0.066197076394253 };

// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double Elem_Pri20::m_Psi[m_NumIP][m_NumNodes] = {
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
inline const double Elem_Pri20::m_DPsi[m_NumIP][m_NumNodes][m_Dim] = {
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

