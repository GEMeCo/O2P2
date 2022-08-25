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
// Solid element, with quadratic interpolation functions, prism shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"

/** @ingroup Elements
  * @class Elem_Pri18
  *
  * @brief Prismatic quadratic element with 18 nodes.
  * @details Solid element, with quadratic interpolation functions, prism shape.
  * @image html Elem_Pri18.png height=300
  */
class Elem_Pri18 : public ElementSolid
{
private:
	Elem_Pri18() = delete;

public:
	/** Constructor for prism quadratic elements.
	  * @param Material Pointer to Material class.
	  */
	explicit Elem_Pri18(std::shared_ptr<Material>& Material)
		: ElementSolid(Material) { };

	// Output function for AcadView, based on element index.
	const std::string printByIndex_AV(const size_t add) const override {
		std::stringstream msg;

		msg << "2 2 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " "	<< this->v_Conect[2]->m_index + add << " "
			<< this->v_Conect[3]->m_index + add << " " << this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "2 2 " << this->v_Conect[12]->m_index + add << " " << this->v_Conect[13]->m_index + add << " " << this->v_Conect[14]->m_index + add << " "
			<< this->v_Conect[15]->m_index + add << " " << this->v_Conect[16]->m_index + add << " " << this->v_Conect[17]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "3 2 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " " << this->v_Conect[2]->m_index + add << " "
			<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[7]->m_index + add << " " << this->v_Conect[8]->m_index + add << " "
			<< this->v_Conect[12]->m_index + add << " " << this->v_Conect[13]->m_index + add << " " << this->v_Conect[14]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "3 2 " << this->v_Conect[2]->m_index + add << " " << this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " "
			<< this->v_Conect[8]->m_index + add << " " << this->v_Conect[10]->m_index + add << " " << this->v_Conect[11]->m_index + add << " "
			<< this->v_Conect[14]->m_index + add << " " << this->v_Conect[16]->m_index + add << " " << this->v_Conect[17]->m_index + add << " " << this->m_Mat->m_index << "\n";
		msg << "3 2 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[3]->m_index + add << " " << this->v_Conect[5]->m_index + add << " "
			<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[9]->m_index + add << " " << this->v_Conect[11]->m_index + add << " "
			<< this->v_Conect[12]->m_index + add << " " << this->v_Conect[15]->m_index + add << " " << this->v_Conect[17]->m_index + add << " " << this->m_Mat->m_index << "\n";

		return msg.str();
	};

	// Output function for AcadView, based on element node number.
	const std::string printByAdder_AV(const size_t add) const override {
		std::stringstream msg;

		msg << "2 2 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (4 + add) << " " 
			<< (5 + add) << " " << (6 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "2 2 " << (13 + add) << " " << (14 + add) << " " << (15 + add) << " " << (16 + add) << " " 
			<< (17 + add) << " " << (18 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "3 2 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " " << (7 + add) << " " << (8 + add) << " "
			<< (9 + add) << " " << (13 + add) << " " << (14 + add) << " " << (15 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "3 2 " << (3 + add) << " " << (5 + add) << " " << (6 + add) << " " << (9 + add) << " " << (11 + add) << " "
			<< (12 + add) << " " << (15 + add) << " " << (17 + add) << " " << (18 + add) << " " << this->m_Mat->m_index << "\n";
		msg << "3 2 " << (1 + add) << " " << (4 + add) << " " << (6 + add) << " " << (7 + add) << " " << (10 + add) << " "
			<< (12 + add) << " " << (13 + add) << " " << (16 + add) << " " << (18 + add) << " " << this->m_Mat->m_index << "\n";

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
	static const int m_NumNodes{ 18 };

	/** @brief Number of Integration Points */
	static const int m_NumIP{ 18 };

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
inline Eigen::VectorXd Elem_Pri18::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(18);

	Psi(0) = 0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. - Point[2]) * Point[2];
	Psi(1) = 2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. - Point[2]) * Point[2];
	Psi(2) = -0.5 * Point[0] * (2. * Point[0] - 1.) * (1. - Point[2]) * Point[2];
	Psi(3) = 2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. - Point[2]) * Point[2];
	Psi(4) = -2. * Point[0] * Point[1] * Point[2] * (1. - Point[2]);
	Psi(5) = -0.5 * Point[1] * (2. * Point[1] - 1.) * (1. - Point[2]) * Point[2];
	Psi(6) = (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (-1. + Point[2] * Point[2]);
	Psi(7) = 4. * Point[0] * (-1. + Point[0] + Point[1]) * (-1. + Point[2] * Point[2]);
	Psi(8) = Point[0] * (2. * Point[0] - 1.) * (1. - Point[2] * Point[2]);
	Psi(9)  = 4. * Point[1] * (-1. + Point[0] + Point[1]) * (-1. + Point[2] * Point[2]);
	Psi(10) = -4. * Point[0] * Point[1] * (-1. + Point[2] * Point[2]);
	Psi(11) = Point[1] * (2. * Point[1] - 1.) * (1. - Point[2] * Point[2]);
	Psi(12) = -0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. + Point[2]) * Point[2];
	Psi(13) = -2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. + Point[2]) * Point[2];
	Psi(14) = 0.5 * Point[0] * (2. * Point[0] - 1.) * (1. + Point[2]) * Point[2];
	Psi(15) = -2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. + Point[2]) * Point[2];
	Psi(16) = 2. * Point[0] * Point[1] * (1. + Point[2]) * Point[2];
	Psi(17) = 0.5 * Point[1] * (2. * Point[1] - 1.) * (1. + Point[2]) * Point[2];

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd Elem_Pri18::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(18, 3);

	DPsi(0, 0) = 0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. - Point[2]) * Point[2];
	DPsi(1, 0) = (-2. + 4. * Point[0] + 2. * Point[1]) * (1. - Point[2]) * Point[2];
	DPsi(2, 0) = (0.5 - 2. * Point[0]) * (1. - Point[2]) * Point[2];
	DPsi(3, 0) = 2. * Point[1] * (1. - Point[2]) * Point[2];
	DPsi(4, 0) = -2. * Point[1] * Point[2] * (1. - Point[2]);
	DPsi(5, 0) = 0.;
	DPsi(6, 0) = (3. - 4. * Point[0] - 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	DPsi(7, 0) = (-4. + 8. * Point[0] + 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	DPsi(8, 0) = (1. - 4. * Point[0]) * (-1. + Point[2] * Point[2]);
	DPsi(9, 0) =  4. * Point[1] * (-1. + Point[2] * Point[2]);
	DPsi(10, 0) = -4. * Point[1] * (-1. + Point[2] * Point[2]);
	DPsi(11, 0) = 0.;
	DPsi(12, 0) = -0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. + Point[2]) * Point[2];
	DPsi(13, 0) = (2. - 4. * Point[0] - 2. * Point[1]) * (1. + Point[2]) * Point[2];
	DPsi(14, 0) = (-0.5 + 2. * Point[0]) * (1. + Point[2]) * Point[2];
	DPsi(15, 0) = -2. * Point[1] * (1. + Point[2]) * Point[2];
	DPsi(16, 0) = 2. * Point[1] * Point[2] * (1. + Point[2]);
	DPsi(17, 0) = 0.;

	DPsi(0, 1) = 0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. - Point[2]) * Point[2];
	DPsi(1, 1) = 2. * Point[0] * (1. - Point[2]) * Point[2];
	DPsi(2, 1) = 0.;
	DPsi(3, 1) = (-2. + 2. * Point[0] + 4. * Point[1]) * (1. - Point[2]) * Point[2];
	DPsi(4, 1) = -2. * Point[0] * Point[2] * (1. - Point[2]);
	DPsi(5, 1) = (0.5 - 2. * Point[1]) * (1. - Point[2]) * Point[2];
	DPsi(6, 1) = (3. - 4. * Point[0] - 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	DPsi(7, 1) = 4. * Point[0] * (-1. + Point[2] * Point[2]);
	DPsi(8, 1) = 0.;
	DPsi(9, 1) =  (-4. + 4. * Point[0] + 8. * Point[1]) * (-1. + Point[2] * Point[2]);
	DPsi(10, 1) = -4. * Point[0] * (-1. + Point[2] * Point[2]);
	DPsi(11, 1) = (1. - 4. * Point[1]) * (-1. + Point[2] * Point[2]);
	DPsi(12, 1) = -0.5 * (3. - 4. * Point[0] - 4. * Point[1]) * (1. + Point[2]) * Point[2];
	DPsi(13, 1) = -2. * Point[0] * (1. + Point[2]) * Point[2];
	DPsi(14, 1) = 0.;
	DPsi(15, 1) = (2. - 2. * Point[0] - 4. * Point[1]) * (1. + Point[2]) * Point[2];
	DPsi(16, 1) = 2. * Point[0] * Point[2] * (1. + Point[2]);
	DPsi(17, 1) = (-0.5 + 2. * Point[1]) * (1. + Point[2]) * Point[2];

	DPsi(0, 2) = 0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. - 2. * Point[2]);
	DPsi(1, 2) = 2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[2]);
	DPsi(2, 2) = -0.5 * Point[0] * (2. * Point[0] - 1.) * (1. - 2. * Point[2]);
	DPsi(3, 2) = 2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[2]);
	DPsi(4, 2) = -2. * Point[0] * Point[1] * (1. - 2. * Point[2]);
	DPsi(5, 2) = -0.5 * Point[1] * (2. * Point[1] - 1.) * (1. - 2. * Point[2]);
	DPsi(6, 2) = (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (2. * Point[2]);
	DPsi(7, 2) = 4. * Point[0] * (-1. + Point[0] + Point[1]) * (2. * Point[2]);
	DPsi(8, 2) = Point[0] * (2. * Point[0] - 1.) * (-2. * Point[2]);
	DPsi(9, 2) =  4. * Point[1] * (-1. + Point[0] + Point[1]) * (2. * Point[2]);
	DPsi(10, 2) = -4. * Point[0] * Point[1] * (2. * Point[2]);
	DPsi(11, 2) = Point[1] * (2. * Point[1] - 1.) * (-2. * Point[2]);
	DPsi(12, 2) = -0.5 * (-1. + Point[0] + Point[1]) * (1. - 2. * Point[0] - 2. * Point[1]) * (1. + 2. * Point[2]);
	DPsi(13, 2) = -2. * Point[0] * (-1. + Point[0] + Point[1]) * (1. + 2. * Point[2]);
	DPsi(14, 2) = 0.5 * Point[0] * (2. * Point[0] - 1.) * (1. + 2. * Point[2]);
	DPsi(15, 2) = -2. * Point[1] * (-1. + Point[0] + Point[1]) * (1. + 2. * Point[2]);
	DPsi(16, 2) = 2. * Point[0] * Point[1] * (1. + 2. * Point[2]);
	DPsi(17, 2) = 0.5 * Point[1] * (2. * Point[1] - 1.) * (1. + 2. * Point[2]);

	return DPsi;
};


// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void Elem_Pri18::setGeomProperties() {

	const int nVertices = 6;

	// Allocate an array with size m_Dim to which m_Centroid points to.
	m_Centroid = std::make_unique<double[]>(m_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<Node<m_Dim>*, nVertices> vertices;
	vertices[0] = v_Conect[0].get();
	vertices[1] = v_Conect[2].get();
	vertices[2] = v_Conect[5].get();
	vertices[3] = v_Conect[12].get();
	vertices[4] = v_Conect[14].get();
	vertices[5] = v_Conect[17].get();

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
inline const double Elem_Pri18::m_xsi[m_NumIP][m_Dim] =
	{ { 0.816847572980459, 0.091576213509771, -0.774596669241483 },
	  { 0.091576213509771, 0.816847572980459, -0.774596669241483 },
	  { 0.091576213509771, 0.091576213509771, -0.774596669241483 },
	  { 0.816847572980459, 0.091576213509771,  0.774596669241483 },
	  { 0.091576213509771, 0.816847572980459,  0.774596669241483 },
	  { 0.091576213509771, 0.091576213509771,  0.774596669241483 },
	  { 0.108103018168070, 0.445948490915965, -0.774596669241483 },
	  { 0.445948490915965, 0.108103018168070, -0.774596669241483 },
	  { 0.445948490915965, 0.445948490915965, -0.774596669241483 },
	  { 0.816847572980459, 0.091576213509771,  0. },
	  { 0.091576213509771, 0.816847572980459,  0. },
	  { 0.091576213509771, 0.091576213509771,  0. },
	  { 0.108103018168070, 0.445948490915965,  0.774596669241483 },
	  { 0.445948490915965, 0.108103018168070,  0.774596669241483 },
	  { 0.445948490915965, 0.445948490915965,  0.774596669241483 },
	  { 0.108103018168070, 0.445948490915965,  0. },
	  { 0.445948490915965, 0.108103018168070,  0. },
	  { 0.445948490915965, 0.445948490915965,  0. } };

// ================================================================================================
//
// Weights for numerical integration
//
// ================================================================================================
inline const double Elem_Pri18::m_weight[m_NumIP] = { 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.03054215101536722222222, 0.06205044157722527777778, 0.06205044157722527777778, 0.06205044157722527777778, 0.04886744162458755555556, 0.04886744162458755555556, 0.04886744162458755555556, 0.06205044157722527777778, 0.06205044157722527777778, 0.06205044157722527777778, 0.09928070652356044444444, 0.09928070652356044444444, 0.09928070652356044444444 };


// ================================================================================================
//
// Shape function
//
// ================================================================================================
inline const double Elem_Pri18::m_Psi[m_NumIP][m_NumNodes] = {
	{ 0.5 * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2],
	  2. * m_xsi[0][0] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2],
	  -0.5 * m_xsi[0][0] * (2. * m_xsi[0][0] - 1.) * (1. - m_xsi[0][2]) * m_xsi[0][2],
	  2. * m_xsi[0][1] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2],
	  -2. * m_xsi[0][0] * m_xsi[0][1] * m_xsi[0][2] * (1. - m_xsi[0][2]),
	  -0.5 * m_xsi[0][1] * (2. * m_xsi[0][1] - 1.) * (1. - m_xsi[0][2]) * m_xsi[0][2],
	  (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]),
	  4. * m_xsi[0][0] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]),
	  m_xsi[0][0] * (2. * m_xsi[0][0] - 1.) * (1. - m_xsi[0][2] * m_xsi[0][2]),
	  4. * m_xsi[0][1] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]),
	  -4. * m_xsi[0][0] * m_xsi[0][1] * (-1. + m_xsi[0][2] * m_xsi[0][2]),
	  m_xsi[0][1] * (2. * m_xsi[0][1] - 1.) * (1. - m_xsi[0][2] * m_xsi[0][2]),
	  -0.5 * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2],
	  -2. * m_xsi[0][0] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2],
	  0.5 * m_xsi[0][0] * (2. * m_xsi[0][0] - 1.) * (1. + m_xsi[0][2]) * m_xsi[0][2],
	  -2. * m_xsi[0][1] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2],
	  2. * m_xsi[0][0] * m_xsi[0][1] * (1. + m_xsi[0][2]) * m_xsi[0][2],
	  0.5 * m_xsi[0][1] * (2. * m_xsi[0][1] - 1.) * (1. + m_xsi[0][2]) * m_xsi[0][2] },

	{ 0.5 * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2],
	  2. * m_xsi[1][0] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2],
	  -0.5 * m_xsi[1][0] * (2. * m_xsi[1][0] - 1.) * (1. - m_xsi[1][2]) * m_xsi[1][2],
	  2. * m_xsi[1][1] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2],
	  -2. * m_xsi[1][0] * m_xsi[1][1] * m_xsi[1][2] * (1. - m_xsi[1][2]),
	  -0.5 * m_xsi[1][1] * (2. * m_xsi[1][1] - 1.) * (1. - m_xsi[1][2]) * m_xsi[1][2],
	  (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]),
	  4. * m_xsi[1][0] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]),
	  m_xsi[1][0] * (2. * m_xsi[1][0] - 1.) * (1. - m_xsi[1][2] * m_xsi[1][2]),
	  4. * m_xsi[1][1] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]),
	  -4. * m_xsi[1][0] * m_xsi[1][1] * (-1. + m_xsi[1][2] * m_xsi[1][2]),
	  m_xsi[1][1] * (2. * m_xsi[1][1] - 1.) * (1. - m_xsi[1][2] * m_xsi[1][2]),
	  -0.5 * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2],
	  -2. * m_xsi[1][0] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2],
	  0.5 * m_xsi[1][0] * (2. * m_xsi[1][0] - 1.) * (1. + m_xsi[1][2]) * m_xsi[1][2],
	  -2. * m_xsi[1][1] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2],
	  2. * m_xsi[1][0] * m_xsi[1][1] * (1. + m_xsi[1][2]) * m_xsi[1][2],
	  0.5 * m_xsi[1][1] * (2. * m_xsi[1][1] - 1.) * (1. + m_xsi[1][2]) * m_xsi[1][2] },

	{ 0.5 * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2],
	  2. * m_xsi[2][0] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2],
	  -0.5 * m_xsi[2][0] * (2. * m_xsi[2][0] - 1.) * (1. - m_xsi[2][2]) * m_xsi[2][2],
	  2. * m_xsi[2][1] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2],
	  -2. * m_xsi[2][0] * m_xsi[2][1] * m_xsi[2][2] * (1. - m_xsi[2][2]),
	  -0.5 * m_xsi[2][1] * (2. * m_xsi[2][1] - 1.) * (1. - m_xsi[2][2]) * m_xsi[2][2],
	  (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]),
	  4. * m_xsi[2][0] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]),
	  m_xsi[2][0] * (2. * m_xsi[2][0] - 1.) * (1. - m_xsi[2][2] * m_xsi[2][2]),
	  4. * m_xsi[2][1] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]),
	  -4. * m_xsi[2][0] * m_xsi[2][1] * (-1. + m_xsi[2][2] * m_xsi[2][2]),
	  m_xsi[2][1] * (2. * m_xsi[2][1] - 1.) * (1. - m_xsi[2][2] * m_xsi[2][2]),
	  -0.5 * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2],
	  -2. * m_xsi[2][0] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2],
	  0.5 * m_xsi[2][0] * (2. * m_xsi[2][0] - 1.) * (1. + m_xsi[2][2]) * m_xsi[2][2],
	  -2. * m_xsi[2][1] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2],
	  2. * m_xsi[2][0] * m_xsi[2][1] * (1. + m_xsi[2][2]) * m_xsi[2][2],
	  0.5 * m_xsi[2][1] * (2. * m_xsi[2][1] - 1.) * (1. + m_xsi[2][2]) * m_xsi[2][2] },

	{ 0.5 * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2],
	  2. * m_xsi[3][0] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2],
	  -0.5 * m_xsi[3][0] * (2. * m_xsi[3][0] - 1.) * (1. - m_xsi[3][2]) * m_xsi[3][2],
	  2. * m_xsi[3][1] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2],
	  -2. * m_xsi[3][0] * m_xsi[3][1] * m_xsi[3][2] * (1. - m_xsi[3][2]),
	  -0.5 * m_xsi[3][1] * (2. * m_xsi[3][1] - 1.) * (1. - m_xsi[3][2]) * m_xsi[3][2],
	  (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]),
	  4. * m_xsi[3][0] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]),
	  m_xsi[3][0] * (2. * m_xsi[3][0] - 1.) * (1. - m_xsi[3][2] * m_xsi[3][2]),
	  4. * m_xsi[3][1] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]),
	  -4. * m_xsi[3][0] * m_xsi[3][1] * (-1. + m_xsi[3][2] * m_xsi[3][2]),
	  m_xsi[3][1] * (2. * m_xsi[3][1] - 1.) * (1. - m_xsi[3][2] * m_xsi[3][2]),
	  -0.5 * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2],
	  -2. * m_xsi[3][0] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2],
	  0.5 * m_xsi[3][0] * (2. * m_xsi[3][0] - 1.) * (1. + m_xsi[3][2]) * m_xsi[3][2],
	  -2. * m_xsi[3][1] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2],
	  2. * m_xsi[3][0] * m_xsi[3][1] * (1. + m_xsi[3][2]) * m_xsi[3][2],
	  0.5 * m_xsi[3][1] * (2. * m_xsi[3][1] - 1.) * (1. + m_xsi[3][2]) * m_xsi[3][2] },

	{ 0.5 * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2],
	  2. * m_xsi[4][0] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2],
	  -0.5 * m_xsi[4][0] * (2. * m_xsi[4][0] - 1.) * (1. - m_xsi[4][2]) * m_xsi[4][2],
	  2. * m_xsi[4][1] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2],
	  -2. * m_xsi[4][0] * m_xsi[4][1] * m_xsi[4][2] * (1. - m_xsi[4][2]),
	  -0.5 * m_xsi[4][1] * (2. * m_xsi[4][1] - 1.) * (1. - m_xsi[4][2]) * m_xsi[4][2],
	  (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]),
	  4. * m_xsi[4][0] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]),
	  m_xsi[4][0] * (2. * m_xsi[4][0] - 1.) * (1. - m_xsi[4][2] * m_xsi[4][2]),
	  4. * m_xsi[4][1] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]),
	  -4. * m_xsi[4][0] * m_xsi[4][1] * (-1. + m_xsi[4][2] * m_xsi[4][2]),
	  m_xsi[4][1] * (2. * m_xsi[4][1] - 1.) * (1. - m_xsi[4][2] * m_xsi[4][2]),
	  -0.5 * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2],
	  -2. * m_xsi[4][0] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2],
	  0.5 * m_xsi[4][0] * (2. * m_xsi[4][0] - 1.) * (1. + m_xsi[4][2]) * m_xsi[4][2],
	  -2. * m_xsi[4][1] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2],
	  2. * m_xsi[4][0] * m_xsi[4][1] * (1. + m_xsi[4][2]) * m_xsi[4][2],
	  0.5 * m_xsi[4][1] * (2. * m_xsi[4][1] - 1.) * (1. + m_xsi[4][2]) * m_xsi[4][2] },

	{ 0.5 * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2],
	  2. * m_xsi[5][0] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2],
	  -0.5 * m_xsi[5][0] * (2. * m_xsi[5][0] - 1.) * (1. - m_xsi[5][2]) * m_xsi[5][2],
	  2. * m_xsi[5][1] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2],
	  -2. * m_xsi[5][0] * m_xsi[5][1] * m_xsi[5][2] * (1. - m_xsi[5][2]),
	  -0.5 * m_xsi[5][1] * (2. * m_xsi[5][1] - 1.) * (1. - m_xsi[5][2]) * m_xsi[5][2],
	  (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]),
	  4. * m_xsi[5][0] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]),
	  m_xsi[5][0] * (2. * m_xsi[5][0] - 1.) * (1. - m_xsi[5][2] * m_xsi[5][2]),
	  4. * m_xsi[5][1] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]),
	  -4. * m_xsi[5][0] * m_xsi[5][1] * (-1. + m_xsi[5][2] * m_xsi[5][2]),
	  m_xsi[5][1] * (2. * m_xsi[5][1] - 1.) * (1. - m_xsi[5][2] * m_xsi[5][2]),
	  -0.5 * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2],
	  -2. * m_xsi[5][0] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2],
	  0.5 * m_xsi[5][0] * (2. * m_xsi[5][0] - 1.) * (1. + m_xsi[5][2]) * m_xsi[5][2],
	  -2. * m_xsi[5][1] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2],
	  2. * m_xsi[5][0] * m_xsi[5][1] * (1. + m_xsi[5][2]) * m_xsi[5][2],
	  0.5 * m_xsi[5][1] * (2. * m_xsi[5][1] - 1.) * (1. + m_xsi[5][2]) * m_xsi[5][2] },

	{ 0.5 * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2],
	  2. * m_xsi[6][0] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2],
	  -0.5 * m_xsi[6][0] * (2. * m_xsi[6][0] - 1.) * (1. - m_xsi[6][2]) * m_xsi[6][2],
	  2. * m_xsi[6][1] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2],
	  -2. * m_xsi[6][0] * m_xsi[6][1] * m_xsi[6][2] * (1. - m_xsi[6][2]),
	  -0.5 * m_xsi[6][1] * (2. * m_xsi[6][1] - 1.) * (1. - m_xsi[6][2]) * m_xsi[6][2],
	  (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]),
	  4. * m_xsi[6][0] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]),
	  m_xsi[6][0] * (2. * m_xsi[6][0] - 1.) * (1. - m_xsi[6][2] * m_xsi[6][2]),
	  4. * m_xsi[6][1] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]),
	  -4. * m_xsi[6][0] * m_xsi[6][1] * (-1. + m_xsi[6][2] * m_xsi[6][2]),
	  m_xsi[6][1] * (2. * m_xsi[6][1] - 1.) * (1. - m_xsi[6][2] * m_xsi[6][2]),
	  -0.5 * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2],
	  -2. * m_xsi[6][0] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2],
	  0.5 * m_xsi[6][0] * (2. * m_xsi[6][0] - 1.) * (1. + m_xsi[6][2]) * m_xsi[6][2],
	  -2. * m_xsi[6][1] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2],
	  2. * m_xsi[6][0] * m_xsi[6][1] * (1. + m_xsi[6][2]) * m_xsi[6][2],
	  0.5 * m_xsi[6][1] * (2. * m_xsi[6][1] - 1.) * (1. + m_xsi[6][2]) * m_xsi[6][2] },

	{ 0.5 * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2],
	  2. * m_xsi[7][0] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2],
	  -0.5 * m_xsi[7][0] * (2. * m_xsi[7][0] - 1.) * (1. - m_xsi[7][2]) * m_xsi[7][2],
	  2. * m_xsi[7][1] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2],
	  -2. * m_xsi[7][0] * m_xsi[7][1] * m_xsi[7][2] * (1. - m_xsi[7][2]),
	  -0.5 * m_xsi[7][1] * (2. * m_xsi[7][1] - 1.) * (1. - m_xsi[7][2]) * m_xsi[7][2],
	  (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]),
	  4. * m_xsi[7][0] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]),
	  m_xsi[7][0] * (2. * m_xsi[7][0] - 1.) * (1. - m_xsi[7][2] * m_xsi[7][2]),
	  4. * m_xsi[7][1] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]),
	  -4. * m_xsi[7][0] * m_xsi[7][1] * (-1. + m_xsi[7][2] * m_xsi[7][2]),
	  m_xsi[7][1] * (2. * m_xsi[7][1] - 1.) * (1. - m_xsi[7][2] * m_xsi[7][2]),
	  -0.5 * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2],
	  -2. * m_xsi[7][0] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2],
	  0.5 * m_xsi[7][0] * (2. * m_xsi[7][0] - 1.) * (1. + m_xsi[7][2]) * m_xsi[7][2],
	  -2. * m_xsi[7][1] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2],
	  2. * m_xsi[7][0] * m_xsi[7][1] * (1. + m_xsi[7][2]) * m_xsi[7][2],
	  0.5 * m_xsi[7][1] * (2. * m_xsi[7][1] - 1.) * (1. + m_xsi[7][2]) * m_xsi[7][2] },

	{ 0.5 * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2],
	  2. * m_xsi[8][0] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2],
	  -0.5 * m_xsi[8][0] * (2. * m_xsi[8][0] - 1.) * (1. - m_xsi[8][2]) * m_xsi[8][2],
	  2. * m_xsi[8][1] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2],
	  -2. * m_xsi[8][0] * m_xsi[8][1] * m_xsi[8][2] * (1. - m_xsi[8][2]),
	  -0.5 * m_xsi[8][1] * (2. * m_xsi[8][1] - 1.) * (1. - m_xsi[8][2]) * m_xsi[8][2],
	  (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]),
	  4. * m_xsi[8][0] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]),
	  m_xsi[8][0] * (2. * m_xsi[8][0] - 1.) * (1. - m_xsi[8][2] * m_xsi[8][2]),
	  4. * m_xsi[8][1] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]),
	  -4. * m_xsi[8][0] * m_xsi[8][1] * (-1. + m_xsi[8][2] * m_xsi[8][2]),
	  m_xsi[8][1] * (2. * m_xsi[8][1] - 1.) * (1. - m_xsi[8][2] * m_xsi[8][2]),
	  -0.5 * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2],
	  -2. * m_xsi[8][0] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2],
	  0.5 * m_xsi[8][0] * (2. * m_xsi[8][0] - 1.) * (1. + m_xsi[8][2]) * m_xsi[8][2],
	  -2. * m_xsi[8][1] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2],
	  2. * m_xsi[8][0] * m_xsi[8][1] * (1. + m_xsi[8][2]) * m_xsi[8][2],
	  0.5 * m_xsi[8][1] * (2. * m_xsi[8][1] - 1.) * (1. + m_xsi[8][2]) * m_xsi[8][2] },

	{ 0.5 * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2],
	  2. * m_xsi[9][0] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2],
	  -0.5 * m_xsi[9][0] * (2. * m_xsi[9][0] - 1.) * (1. - m_xsi[9][2]) * m_xsi[9][2],
	  2. * m_xsi[9][1] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2],
	  -2. * m_xsi[9][0] * m_xsi[9][1] * m_xsi[9][2] * (1. - m_xsi[9][2]),
	  -0.5 * m_xsi[9][1] * (2. * m_xsi[9][1] - 1.) * (1. - m_xsi[9][2]) * m_xsi[9][2],
	  (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]),
	  4. * m_xsi[9][0] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]),
	  m_xsi[9][0] * (2. * m_xsi[9][0] - 1.) * (1. - m_xsi[9][2] * m_xsi[9][2]),
	  4. * m_xsi[9][1] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]),
	  -4. * m_xsi[9][0] * m_xsi[9][1] * (-1. + m_xsi[9][2] * m_xsi[9][2]),
	  m_xsi[9][1] * (2. * m_xsi[9][1] - 1.) * (1. - m_xsi[9][2] * m_xsi[9][2]),
	  -0.5 * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2],
	  -2. * m_xsi[9][0] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2],
	  0.5 * m_xsi[9][0] * (2. * m_xsi[9][0] - 1.) * (1. + m_xsi[9][2]) * m_xsi[9][2],
	  -2. * m_xsi[9][1] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2],
	  2. * m_xsi[9][0] * m_xsi[9][1] * (1. + m_xsi[9][2]) * m_xsi[9][2],
	  0.5 * m_xsi[9][1] * (2. * m_xsi[9][1] - 1.) * (1. + m_xsi[9][2]) * m_xsi[9][2] },

	{ 0.5 * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2],
	  2. * m_xsi[10][0] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2],
	  -0.5 * m_xsi[10][0] * (2. * m_xsi[10][0] - 1.) * (1. - m_xsi[10][2]) * m_xsi[10][2],
	  2. * m_xsi[10][1] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2],
	  -2. * m_xsi[10][0] * m_xsi[10][1] * m_xsi[10][2] * (1. - m_xsi[10][2]),
	  -0.5 * m_xsi[10][1] * (2. * m_xsi[10][1] - 1.) * (1. - m_xsi[10][2]) * m_xsi[10][2],
	  (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]),
	  4. * m_xsi[10][0] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]),
	  m_xsi[10][0] * (2. * m_xsi[10][0] - 1.) * (1. - m_xsi[10][2] * m_xsi[10][2]),
	  4. * m_xsi[10][1] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]),
	  -4. * m_xsi[10][0] * m_xsi[10][1] * (-1. + m_xsi[10][2] * m_xsi[10][2]),
	  m_xsi[10][1] * (2. * m_xsi[10][1] - 1.) * (1. - m_xsi[10][2] * m_xsi[10][2]),
	  -0.5 * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2],
	  -2. * m_xsi[10][0] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2],
	  0.5 * m_xsi[10][0] * (2. * m_xsi[10][0] - 1.) * (1. + m_xsi[10][2]) * m_xsi[10][2],
	  -2. * m_xsi[10][1] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2],
	  2. * m_xsi[10][0] * m_xsi[10][1] * (1. + m_xsi[10][2]) * m_xsi[10][2],
	  0.5 * m_xsi[10][1] * (2. * m_xsi[10][1] - 1.) * (1. + m_xsi[10][2]) * m_xsi[10][2] },

	{ 0.5 * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2],
	  2. * m_xsi[11][0] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2],
	  -0.5 * m_xsi[11][0] * (2. * m_xsi[11][0] - 1.) * (1. - m_xsi[11][2]) * m_xsi[11][2],
	  2. * m_xsi[11][1] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2],
	  -2. * m_xsi[11][0] * m_xsi[11][1] * m_xsi[11][2] * (1. - m_xsi[11][2]),
	  -0.5 * m_xsi[11][1] * (2. * m_xsi[11][1] - 1.) * (1. - m_xsi[11][2]) * m_xsi[11][2],
	  (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]),
	  4. * m_xsi[11][0] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]),
	  m_xsi[11][0] * (2. * m_xsi[11][0] - 1.) * (1. - m_xsi[11][2] * m_xsi[11][2]),
	  4. * m_xsi[11][1] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]),
	  -4. * m_xsi[11][0] * m_xsi[11][1] * (-1. + m_xsi[11][2] * m_xsi[11][2]),
	  m_xsi[11][1] * (2. * m_xsi[11][1] - 1.) * (1. - m_xsi[11][2] * m_xsi[11][2]),
	  -0.5 * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2],
	  -2. * m_xsi[11][0] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2],
	  0.5 * m_xsi[11][0] * (2. * m_xsi[11][0] - 1.) * (1. + m_xsi[11][2]) * m_xsi[11][2],
	  -2. * m_xsi[11][1] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2],
	  2. * m_xsi[11][0] * m_xsi[11][1] * (1. + m_xsi[11][2]) * m_xsi[11][2],
	  0.5 * m_xsi[11][1] * (2. * m_xsi[11][1] - 1.) * (1. + m_xsi[11][2]) * m_xsi[11][2] },

	{ 0.5 * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2],
	  2. * m_xsi[12][0] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2],
	  -0.5 * m_xsi[12][0] * (2. * m_xsi[12][0] - 1.) * (1. - m_xsi[12][2]) * m_xsi[12][2],
	  2. * m_xsi[12][1] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2],
	  -2. * m_xsi[12][0] * m_xsi[12][1] * m_xsi[12][2] * (1. - m_xsi[12][2]),
	  -0.5 * m_xsi[12][1] * (2. * m_xsi[12][1] - 1.) * (1. - m_xsi[12][2]) * m_xsi[12][2],
	  (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]),
	  4. * m_xsi[12][0] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]),
	  m_xsi[12][0] * (2. * m_xsi[12][0] - 1.) * (1. - m_xsi[12][2] * m_xsi[12][2]),
	  4. * m_xsi[12][1] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]),
	  -4. * m_xsi[12][0] * m_xsi[12][1] * (-1. + m_xsi[12][2] * m_xsi[12][2]),
	  m_xsi[12][1] * (2. * m_xsi[12][1] - 1.) * (1. - m_xsi[12][2] * m_xsi[12][2]),
	  -0.5 * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2],
	  -2. * m_xsi[12][0] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2],
	  0.5 * m_xsi[12][0] * (2. * m_xsi[12][0] - 1.) * (1. + m_xsi[12][2]) * m_xsi[12][2],
	  -2. * m_xsi[12][1] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2],
	  2. * m_xsi[12][0] * m_xsi[12][1] * (1. + m_xsi[12][2]) * m_xsi[12][2],
	  0.5 * m_xsi[12][1] * (2. * m_xsi[12][1] - 1.) * (1. + m_xsi[12][2]) * m_xsi[12][2] },

	{ 0.5 * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2],
	  2. * m_xsi[13][0] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2],
	  -0.5 * m_xsi[13][0] * (2. * m_xsi[13][0] - 1.) * (1. - m_xsi[13][2]) * m_xsi[13][2],
	  2. * m_xsi[13][1] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2],
	  -2. * m_xsi[13][0] * m_xsi[13][1] * m_xsi[13][2] * (1. - m_xsi[13][2]),
	  -0.5 * m_xsi[13][1] * (2. * m_xsi[13][1] - 1.) * (1. - m_xsi[13][2]) * m_xsi[13][2],
	  (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]),
	  4. * m_xsi[13][0] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]),
	  m_xsi[13][0] * (2. * m_xsi[13][0] - 1.) * (1. - m_xsi[13][2] * m_xsi[13][2]),
	  4. * m_xsi[13][1] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]),
	  -4. * m_xsi[13][0] * m_xsi[13][1] * (-1. + m_xsi[13][2] * m_xsi[13][2]),
	  m_xsi[13][1] * (2. * m_xsi[13][1] - 1.) * (1. - m_xsi[13][2] * m_xsi[13][2]),
	  -0.5 * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2],
	  -2. * m_xsi[13][0] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2],
	  0.5 * m_xsi[13][0] * (2. * m_xsi[13][0] - 1.) * (1. + m_xsi[13][2]) * m_xsi[13][2],
	  -2. * m_xsi[13][1] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2],
	  2. * m_xsi[13][0] * m_xsi[13][1] * (1. + m_xsi[13][2]) * m_xsi[13][2],
	  0.5 * m_xsi[13][1] * (2. * m_xsi[13][1] - 1.) * (1. + m_xsi[13][2]) * m_xsi[13][2] },

	{ 0.5 * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2],
	  2. * m_xsi[14][0] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2],
	  -0.5 * m_xsi[14][0] * (2. * m_xsi[14][0] - 1.) * (1. - m_xsi[14][2]) * m_xsi[14][2],
	  2. * m_xsi[14][1] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2],
	  -2. * m_xsi[14][0] * m_xsi[14][1] * m_xsi[14][2] * (1. - m_xsi[14][2]),
	  -0.5 * m_xsi[14][1] * (2. * m_xsi[14][1] - 1.) * (1. - m_xsi[14][2]) * m_xsi[14][2],
	  (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]),
	  4. * m_xsi[14][0] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]),
	  m_xsi[14][0] * (2. * m_xsi[14][0] - 1.) * (1. - m_xsi[14][2] * m_xsi[14][2]),
	  4. * m_xsi[14][1] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]),
	  -4. * m_xsi[14][0] * m_xsi[14][1] * (-1. + m_xsi[14][2] * m_xsi[14][2]),
	  m_xsi[14][1] * (2. * m_xsi[14][1] - 1.) * (1. - m_xsi[14][2] * m_xsi[14][2]),
	  -0.5 * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2],
	  -2. * m_xsi[14][0] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2],
	  0.5 * m_xsi[14][0] * (2. * m_xsi[14][0] - 1.) * (1. + m_xsi[14][2]) * m_xsi[14][2],
	  -2. * m_xsi[14][1] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2],
	  2. * m_xsi[14][0] * m_xsi[14][1] * (1. + m_xsi[14][2]) * m_xsi[14][2],
	  0.5 * m_xsi[14][1] * (2. * m_xsi[14][1] - 1.) * (1. + m_xsi[14][2]) * m_xsi[14][2] },

	{ 0.5 * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2],
	  2. * m_xsi[15][0] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2],
	  -0.5 * m_xsi[15][0] * (2. * m_xsi[15][0] - 1.) * (1. - m_xsi[15][2]) * m_xsi[15][2],
	  2. * m_xsi[15][1] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2],
	  -2. * m_xsi[15][0] * m_xsi[15][1] * m_xsi[15][2] * (1. - m_xsi[15][2]),
	  -0.5 * m_xsi[15][1] * (2. * m_xsi[15][1] - 1.) * (1. - m_xsi[15][2]) * m_xsi[15][2],
	  (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]),
	  4. * m_xsi[15][0] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]),
	  m_xsi[15][0] * (2. * m_xsi[15][0] - 1.) * (1. - m_xsi[15][2] * m_xsi[15][2]),
	  4. * m_xsi[15][1] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]),
	  -4. * m_xsi[15][0] * m_xsi[15][1] * (-1. + m_xsi[15][2] * m_xsi[15][2]),
	  m_xsi[15][1] * (2. * m_xsi[15][1] - 1.) * (1. - m_xsi[15][2] * m_xsi[15][2]),
	  -0.5 * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2],
	  -2. * m_xsi[15][0] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2],
	  0.5 * m_xsi[15][0] * (2. * m_xsi[15][0] - 1.) * (1. + m_xsi[15][2]) * m_xsi[15][2],
	  -2. * m_xsi[15][1] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2],
	  2. * m_xsi[15][0] * m_xsi[15][1] * (1. + m_xsi[15][2]) * m_xsi[15][2],
	  0.5 * m_xsi[15][1] * (2. * m_xsi[15][1] - 1.) * (1. + m_xsi[15][2]) * m_xsi[15][2] },

	{ 0.5 * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2],
	  2. * m_xsi[16][0] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2],
	  -0.5 * m_xsi[16][0] * (2. * m_xsi[16][0] - 1.) * (1. - m_xsi[16][2]) * m_xsi[16][2],
	  2. * m_xsi[16][1] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2],
	  -2. * m_xsi[16][0] * m_xsi[16][1] * m_xsi[16][2] * (1. - m_xsi[16][2]),
	  -0.5 * m_xsi[16][1] * (2. * m_xsi[16][1] - 1.) * (1. - m_xsi[16][2]) * m_xsi[16][2],
	  (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]),
	  4. * m_xsi[16][0] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]),
	  m_xsi[16][0] * (2. * m_xsi[16][0] - 1.) * (1. - m_xsi[16][2] * m_xsi[16][2]),
	  4. * m_xsi[16][1] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]),
	  -4. * m_xsi[16][0] * m_xsi[16][1] * (-1. + m_xsi[16][2] * m_xsi[16][2]),
	  m_xsi[16][1] * (2. * m_xsi[16][1] - 1.) * (1. - m_xsi[16][2] * m_xsi[16][2]),
	  -0.5 * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2],
	  -2. * m_xsi[16][0] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2],
	  0.5 * m_xsi[16][0] * (2. * m_xsi[16][0] - 1.) * (1. + m_xsi[16][2]) * m_xsi[16][2],
	  -2. * m_xsi[16][1] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2],
	  2. * m_xsi[16][0] * m_xsi[16][1] * (1. + m_xsi[16][2]) * m_xsi[16][2],
	  0.5 * m_xsi[16][1] * (2. * m_xsi[16][1] - 1.) * (1. + m_xsi[16][2]) * m_xsi[16][2] },

	{ 0.5 * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2],
	  2. * m_xsi[17][0] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2],
	  -0.5 * m_xsi[17][0] * (2. * m_xsi[17][0] - 1.) * (1. - m_xsi[17][2]) * m_xsi[17][2],
	  2. * m_xsi[17][1] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2],
	  -2. * m_xsi[17][0] * m_xsi[17][1] * m_xsi[17][2] * (1. - m_xsi[17][2]),
	  -0.5 * m_xsi[17][1] * (2. * m_xsi[17][1] - 1.) * (1. - m_xsi[17][2]) * m_xsi[17][2],
	  (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]),
	  4. * m_xsi[17][0] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]),
	  m_xsi[17][0] * (2. * m_xsi[17][0] - 1.) * (1. - m_xsi[17][2] * m_xsi[17][2]),
	  4. * m_xsi[17][1] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]),
	  -4. * m_xsi[17][0] * m_xsi[17][1] * (-1. + m_xsi[17][2] * m_xsi[17][2]),
	  m_xsi[17][1] * (2. * m_xsi[17][1] - 1.) * (1. - m_xsi[17][2] * m_xsi[17][2]),
	  -0.5 * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2],
	  -2. * m_xsi[17][0] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2],
	  0.5 * m_xsi[17][0] * (2. * m_xsi[17][0] - 1.) * (1. + m_xsi[17][2]) * m_xsi[17][2],
	  -2. * m_xsi[17][1] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2],
	  2. * m_xsi[17][0] * m_xsi[17][1] * (1. + m_xsi[17][2]) * m_xsi[17][2],
	  0.5 * m_xsi[17][1] * (2. * m_xsi[17][1] - 1.) * (1. + m_xsi[17][2]) * m_xsi[17][2] } };

// ================================================================================================
//
// Shape functions derivative
//
// ================================================================================================
inline const double Elem_Pri18::m_DPsi[m_NumIP][m_NumNodes][m_Dim] = {
	{ { 0.5 * (3. - 4. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2], 0.5 * (3. - 4. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2], 0.5 * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (1. - 2. * m_xsi[0][2]) },
	  { (-2. + 4. * m_xsi[0][0] + 2. * m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2], 2. * m_xsi[0][0] * (1. - m_xsi[0][2]) * m_xsi[0][2], 2. * m_xsi[0][0] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][2]) },
	  { (0.5 - 2. * m_xsi[0][0]) * (1. - m_xsi[0][2]) * m_xsi[0][2], 0., -0.5 * m_xsi[0][0] * (2. * m_xsi[0][0] - 1.) * (1. - 2. * m_xsi[0][2]) },
	  { 2. * m_xsi[0][1] * (1. - m_xsi[0][2]) * m_xsi[0][2], (-2. + 2. * m_xsi[0][0] + 4. * m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2], 2. * m_xsi[0][1] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][2]) },
	  { -2. * m_xsi[0][1] * m_xsi[0][2] * (1. - m_xsi[0][2]), -2. * m_xsi[0][0] * m_xsi[0][2] * (1. - m_xsi[0][2]), -2. * m_xsi[0][0] * m_xsi[0][1] * (1. - 2. * m_xsi[0][2]) },
	  { 0., (0.5 - 2. * m_xsi[0][1]) * (1. - m_xsi[0][2]) * m_xsi[0][2], -0.5 * m_xsi[0][1] * (2. * m_xsi[0][1] - 1.) * (1. - 2. * m_xsi[0][2]) },
	  { (3. - 4. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]), (3. - 4. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]), (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (2. * m_xsi[0][2]) },
	  { (-4. + 8. * m_xsi[0][0] + 4. * m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]), 4. * m_xsi[0][0] * (-1. + m_xsi[0][2] * m_xsi[0][2]), 4. * m_xsi[0][0] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (2. * m_xsi[0][2]) },
	  { (1. - 4. * m_xsi[0][0]) * (-1. + m_xsi[0][2] * m_xsi[0][2]), 0., m_xsi[0][0] * (2. * m_xsi[0][0] - 1.) * (-2. * m_xsi[0][2]) },
	  { 4. * m_xsi[0][1] * (-1. + m_xsi[0][2] * m_xsi[0][2]), (-4. + 4. * m_xsi[0][0] + 8. * m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]), 4. * m_xsi[0][1] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (2. * m_xsi[0][2]) },
	  { -4. * m_xsi[0][1] * (-1. + m_xsi[0][2] * m_xsi[0][2]), -4. * m_xsi[0][0] * (-1. + m_xsi[0][2] * m_xsi[0][2]), -4. * m_xsi[0][0] * m_xsi[0][1] * (2. * m_xsi[0][2]) },
	  { 0., (1. - 4. * m_xsi[0][1]) * (-1. + m_xsi[0][2] * m_xsi[0][2]), m_xsi[0][1] * (2. * m_xsi[0][1] - 1.) * (-2. * m_xsi[0][2]) },
	  { -0.5 * (3. - 4. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2], -0.5 * (3. - 4. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2], -0.5 * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. - 2. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (1. + 2. * m_xsi[0][2]) },
	  { (2. - 4. * m_xsi[0][0] - 2. * m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2], -2. * m_xsi[0][0] * (1. + m_xsi[0][2]) * m_xsi[0][2], -2. * m_xsi[0][0] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. + 2. * m_xsi[0][2]) },
	  { (-0.5 + 2. * m_xsi[0][0]) * (1. + m_xsi[0][2]) * m_xsi[0][2], 0., 0.5 * m_xsi[0][0] * (2. * m_xsi[0][0] - 1.) * (1. + 2. * m_xsi[0][2]) },
	  { -2. * m_xsi[0][1] * (1. + m_xsi[0][2]) * m_xsi[0][2], (2. - 2. * m_xsi[0][0] - 4. * m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2], -2. * m_xsi[0][1] * (-1. + m_xsi[0][0] + m_xsi[0][1]) * (1. + 2. * m_xsi[0][2]) },
	  { 2. * m_xsi[0][1] * m_xsi[0][2] * (1. + m_xsi[0][2]), 2. * m_xsi[0][0] * m_xsi[0][2] * (1. + m_xsi[0][2]), 2. * m_xsi[0][0] * m_xsi[0][1] * (1. + 2. * m_xsi[0][2]) },
	  { 0., (-0.5 + 2. * m_xsi[0][1]) * (1. + m_xsi[0][2]) * m_xsi[0][2], 0.5 * m_xsi[0][1] * (2. * m_xsi[0][1] - 1.) * (1. + 2. * m_xsi[0][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2], 0.5 * (3. - 4. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2], 0.5 * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (1. - 2. * m_xsi[1][2]) },
	  { (-2. + 4. * m_xsi[1][0] + 2. * m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2], 2. * m_xsi[1][0] * (1. - m_xsi[1][2]) * m_xsi[1][2], 2. * m_xsi[1][0] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][2]) },
	  { (0.5 - 2. * m_xsi[1][0]) * (1. - m_xsi[1][2]) * m_xsi[1][2], 0., -0.5 * m_xsi[1][0] * (2. * m_xsi[1][0] - 1.) * (1. - 2. * m_xsi[1][2]) },
	  { 2. * m_xsi[1][1] * (1. - m_xsi[1][2]) * m_xsi[1][2], (-2. + 2. * m_xsi[1][0] + 4. * m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2], 2. * m_xsi[1][1] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][2]) },
	  { -2. * m_xsi[1][1] * m_xsi[1][2] * (1. - m_xsi[1][2]), -2. * m_xsi[1][0] * m_xsi[1][2] * (1. - m_xsi[1][2]), -2. * m_xsi[1][0] * m_xsi[1][1] * (1. - 2. * m_xsi[1][2]) },
	  { 0., (0.5 - 2. * m_xsi[1][1]) * (1. - m_xsi[1][2]) * m_xsi[1][2], -0.5 * m_xsi[1][1] * (2. * m_xsi[1][1] - 1.) * (1. - 2. * m_xsi[1][2]) },
	  { (3. - 4. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]), (3. - 4. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]), (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (2. * m_xsi[1][2]) },
	  { (-4. + 8. * m_xsi[1][0] + 4. * m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]), 4. * m_xsi[1][0] * (-1. + m_xsi[1][2] * m_xsi[1][2]), 4. * m_xsi[1][0] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (2. * m_xsi[1][2]) },
	  { (1. - 4. * m_xsi[1][0]) * (-1. + m_xsi[1][2] * m_xsi[1][2]), 0., m_xsi[1][0] * (2. * m_xsi[1][0] - 1.) * (-2. * m_xsi[1][2]) },
	  { 4. * m_xsi[1][1] * (-1. + m_xsi[1][2] * m_xsi[1][2]), (-4. + 4. * m_xsi[1][0] + 8. * m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]), 4. * m_xsi[1][1] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (2. * m_xsi[1][2]) },
	  { -4. * m_xsi[1][1] * (-1. + m_xsi[1][2] * m_xsi[1][2]), -4. * m_xsi[1][0] * (-1. + m_xsi[1][2] * m_xsi[1][2]), -4. * m_xsi[1][0] * m_xsi[1][1] * (2. * m_xsi[1][2]) },
	  { 0., (1. - 4. * m_xsi[1][1]) * (-1. + m_xsi[1][2] * m_xsi[1][2]), m_xsi[1][1] * (2. * m_xsi[1][1] - 1.) * (-2. * m_xsi[1][2]) },
	  { -0.5 * (3. - 4. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2], -0.5 * (3. - 4. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2], -0.5 * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. - 2. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (1. + 2. * m_xsi[1][2]) },
	  { (2. - 4. * m_xsi[1][0] - 2. * m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2], -2. * m_xsi[1][0] * (1. + m_xsi[1][2]) * m_xsi[1][2], -2. * m_xsi[1][0] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. + 2. * m_xsi[1][2]) },
	  { (-0.5 + 2. * m_xsi[1][0]) * (1. + m_xsi[1][2]) * m_xsi[1][2], 0., 0.5 * m_xsi[1][0] * (2. * m_xsi[1][0] - 1.) * (1. + 2. * m_xsi[1][2]) },
	  { -2. * m_xsi[1][1] * (1. + m_xsi[1][2]) * m_xsi[1][2], (2. - 2. * m_xsi[1][0] - 4. * m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2], -2. * m_xsi[1][1] * (-1. + m_xsi[1][0] + m_xsi[1][1]) * (1. + 2. * m_xsi[1][2]) },
	  { 2. * m_xsi[1][1] * m_xsi[1][2] * (1. + m_xsi[1][2]), 2. * m_xsi[1][0] * m_xsi[1][2] * (1. + m_xsi[1][2]), 2. * m_xsi[1][0] * m_xsi[1][1] * (1. + 2. * m_xsi[1][2]) },
	  { 0., (-0.5 + 2. * m_xsi[1][1]) * (1. + m_xsi[1][2]) * m_xsi[1][2], 0.5 * m_xsi[1][1] * (2. * m_xsi[1][1] - 1.) * (1. + 2. * m_xsi[1][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2], 0.5 * (3. - 4. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2], 0.5 * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (1. - 2. * m_xsi[2][2]) },
	  { (-2. + 4. * m_xsi[2][0] + 2. * m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2], 2. * m_xsi[2][0] * (1. - m_xsi[2][2]) * m_xsi[2][2], 2. * m_xsi[2][0] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][2]) },
	  { (0.5 - 2. * m_xsi[2][0]) * (1. - m_xsi[2][2]) * m_xsi[2][2], 0., -0.5 * m_xsi[2][0] * (2. * m_xsi[2][0] - 1.) * (1. - 2. * m_xsi[2][2]) },
	  { 2. * m_xsi[2][1] * (1. - m_xsi[2][2]) * m_xsi[2][2], (-2. + 2. * m_xsi[2][0] + 4. * m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2], 2. * m_xsi[2][1] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][2]) },
	  { -2. * m_xsi[2][1] * m_xsi[2][2] * (1. - m_xsi[2][2]), -2. * m_xsi[2][0] * m_xsi[2][2] * (1. - m_xsi[2][2]), -2. * m_xsi[2][0] * m_xsi[2][1] * (1. - 2. * m_xsi[2][2]) },
	  { 0., (0.5 - 2. * m_xsi[2][1]) * (1. - m_xsi[2][2]) * m_xsi[2][2], -0.5 * m_xsi[2][1] * (2. * m_xsi[2][1] - 1.) * (1. - 2. * m_xsi[2][2]) },
	  { (3. - 4. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]), (3. - 4. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]), (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (2. * m_xsi[2][2]) },
	  { (-4. + 8. * m_xsi[2][0] + 4. * m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]), 4. * m_xsi[2][0] * (-1. + m_xsi[2][2] * m_xsi[2][2]), 4. * m_xsi[2][0] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (2. * m_xsi[2][2]) },
	  { (1. - 4. * m_xsi[2][0]) * (-1. + m_xsi[2][2] * m_xsi[2][2]), 0., m_xsi[2][0] * (2. * m_xsi[2][0] - 1.) * (-2. * m_xsi[2][2]) },
	  { 4. * m_xsi[2][1] * (-1. + m_xsi[2][2] * m_xsi[2][2]), (-4. + 4. * m_xsi[2][0] + 8. * m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]), 4. * m_xsi[2][1] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (2. * m_xsi[2][2]) },
	  { -4. * m_xsi[2][1] * (-1. + m_xsi[2][2] * m_xsi[2][2]), -4. * m_xsi[2][0] * (-1. + m_xsi[2][2] * m_xsi[2][2]), -4. * m_xsi[2][0] * m_xsi[2][1] * (2. * m_xsi[2][2]) },
	  { 0., (1. - 4. * m_xsi[2][1]) * (-1. + m_xsi[2][2] * m_xsi[2][2]), m_xsi[2][1] * (2. * m_xsi[2][1] - 1.) * (-2. * m_xsi[2][2]) },
	  { -0.5 * (3. - 4. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2], -0.5 * (3. - 4. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2], -0.5 * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. - 2. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (1. + 2. * m_xsi[2][2]) },
	  { (2. - 4. * m_xsi[2][0] - 2. * m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2], -2. * m_xsi[2][0] * (1. + m_xsi[2][2]) * m_xsi[2][2], -2. * m_xsi[2][0] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. + 2. * m_xsi[2][2]) },
	  { (-0.5 + 2. * m_xsi[2][0]) * (1. + m_xsi[2][2]) * m_xsi[2][2], 0., 0.5 * m_xsi[2][0] * (2. * m_xsi[2][0] - 1.) * (1. + 2. * m_xsi[2][2]) },
	  { -2. * m_xsi[2][1] * (1. + m_xsi[2][2]) * m_xsi[2][2], (2. - 2. * m_xsi[2][0] - 4. * m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2], -2. * m_xsi[2][1] * (-1. + m_xsi[2][0] + m_xsi[2][1]) * (1. + 2. * m_xsi[2][2]) },
	  { 2. * m_xsi[2][1] * m_xsi[2][2] * (1. + m_xsi[2][2]), 2. * m_xsi[2][0] * m_xsi[2][2] * (1. + m_xsi[2][2]), 2. * m_xsi[2][0] * m_xsi[2][1] * (1. + 2. * m_xsi[2][2]) },
	  { 0., (-0.5 + 2. * m_xsi[2][1]) * (1. + m_xsi[2][2]) * m_xsi[2][2], 0.5 * m_xsi[2][1] * (2. * m_xsi[2][1] - 1.) * (1. + 2. * m_xsi[2][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2], 0.5 * (3. - 4. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2], 0.5 * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (1. - 2. * m_xsi[3][2]) },
	  { (-2. + 4. * m_xsi[3][0] + 2. * m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2], 2. * m_xsi[3][0] * (1. - m_xsi[3][2]) * m_xsi[3][2], 2. * m_xsi[3][0] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][2]) },
	  { (0.5 - 2. * m_xsi[3][0]) * (1. - m_xsi[3][2]) * m_xsi[3][2], 0., -0.5 * m_xsi[3][0] * (2. * m_xsi[3][0] - 1.) * (1. - 2. * m_xsi[3][2]) },
	  { 2. * m_xsi[3][1] * (1. - m_xsi[3][2]) * m_xsi[3][2], (-2. + 2. * m_xsi[3][0] + 4. * m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2], 2. * m_xsi[3][1] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][2]) },
	  { -2. * m_xsi[3][1] * m_xsi[3][2] * (1. - m_xsi[3][2]), -2. * m_xsi[3][0] * m_xsi[3][2] * (1. - m_xsi[3][2]), -2. * m_xsi[3][0] * m_xsi[3][1] * (1. - 2. * m_xsi[3][2]) },
	  { 0., (0.5 - 2. * m_xsi[3][1]) * (1. - m_xsi[3][2]) * m_xsi[3][2], -0.5 * m_xsi[3][1] * (2. * m_xsi[3][1] - 1.) * (1. - 2. * m_xsi[3][2]) },
	  { (3. - 4. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]), (3. - 4. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]), (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (2. * m_xsi[3][2]) },
	  { (-4. + 8. * m_xsi[3][0] + 4. * m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]), 4. * m_xsi[3][0] * (-1. + m_xsi[3][2] * m_xsi[3][2]), 4. * m_xsi[3][0] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (2. * m_xsi[3][2]) },
	  { (1. - 4. * m_xsi[3][0]) * (-1. + m_xsi[3][2] * m_xsi[3][2]), 0., m_xsi[3][0] * (2. * m_xsi[3][0] - 1.) * (-2. * m_xsi[3][2]) },
	  { 4. * m_xsi[3][1] * (-1. + m_xsi[3][2] * m_xsi[3][2]), (-4. + 4. * m_xsi[3][0] + 8. * m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]), 4. * m_xsi[3][1] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (2. * m_xsi[3][2]) },
	  { -4. * m_xsi[3][1] * (-1. + m_xsi[3][2] * m_xsi[3][2]), -4. * m_xsi[3][0] * (-1. + m_xsi[3][2] * m_xsi[3][2]), -4. * m_xsi[3][0] * m_xsi[3][1] * (2. * m_xsi[3][2]) },
	  { 0., (1. - 4. * m_xsi[3][1]) * (-1. + m_xsi[3][2] * m_xsi[3][2]), m_xsi[3][1] * (2. * m_xsi[3][1] - 1.) * (-2. * m_xsi[3][2]) },
	  { -0.5 * (3. - 4. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2], -0.5 * (3. - 4. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2], -0.5 * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. - 2. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (1. + 2. * m_xsi[3][2]) },
	  { (2. - 4. * m_xsi[3][0] - 2. * m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2], -2. * m_xsi[3][0] * (1. + m_xsi[3][2]) * m_xsi[3][2], -2. * m_xsi[3][0] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. + 2. * m_xsi[3][2]) },
	  { (-0.5 + 2. * m_xsi[3][0]) * (1. + m_xsi[3][2]) * m_xsi[3][2], 0., 0.5 * m_xsi[3][0] * (2. * m_xsi[3][0] - 1.) * (1. + 2. * m_xsi[3][2]) },
	  { -2. * m_xsi[3][1] * (1. + m_xsi[3][2]) * m_xsi[3][2], (2. - 2. * m_xsi[3][0] - 4. * m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2], -2. * m_xsi[3][1] * (-1. + m_xsi[3][0] + m_xsi[3][1]) * (1. + 2. * m_xsi[3][2]) },
	  { 2. * m_xsi[3][1] * m_xsi[3][2] * (1. + m_xsi[3][2]), 2. * m_xsi[3][0] * m_xsi[3][2] * (1. + m_xsi[3][2]), 2. * m_xsi[3][0] * m_xsi[3][1] * (1. + 2. * m_xsi[3][2]) },
	  { 0., (-0.5 + 2. * m_xsi[3][1]) * (1. + m_xsi[3][2]) * m_xsi[3][2], 0.5 * m_xsi[3][1] * (2. * m_xsi[3][1] - 1.) * (1. + 2. * m_xsi[3][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2], 0.5 * (3. - 4. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2], 0.5 * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (1. - 2. * m_xsi[4][2]) },
	  { (-2. + 4. * m_xsi[4][0] + 2. * m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2], 2. * m_xsi[4][0] * (1. - m_xsi[4][2]) * m_xsi[4][2], 2. * m_xsi[4][0] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][2]) },
	  { (0.5 - 2. * m_xsi[4][0]) * (1. - m_xsi[4][2]) * m_xsi[4][2], 0., -0.5 * m_xsi[4][0] * (2. * m_xsi[4][0] - 1.) * (1. - 2. * m_xsi[4][2]) },
	  { 2. * m_xsi[4][1] * (1. - m_xsi[4][2]) * m_xsi[4][2], (-2. + 2. * m_xsi[4][0] + 4. * m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2], 2. * m_xsi[4][1] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][2]) },
	  { -2. * m_xsi[4][1] * m_xsi[4][2] * (1. - m_xsi[4][2]), -2. * m_xsi[4][0] * m_xsi[4][2] * (1. - m_xsi[4][2]), -2. * m_xsi[4][0] * m_xsi[4][1] * (1. - 2. * m_xsi[4][2]) },
	  { 0., (0.5 - 2. * m_xsi[4][1]) * (1. - m_xsi[4][2]) * m_xsi[4][2], -0.5 * m_xsi[4][1] * (2. * m_xsi[4][1] - 1.) * (1. - 2. * m_xsi[4][2]) },
	  { (3. - 4. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]), (3. - 4. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]), (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (2. * m_xsi[4][2]) },
	  { (-4. + 8. * m_xsi[4][0] + 4. * m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]), 4. * m_xsi[4][0] * (-1. + m_xsi[4][2] * m_xsi[4][2]), 4. * m_xsi[4][0] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (2. * m_xsi[4][2]) },
	  { (1. - 4. * m_xsi[4][0]) * (-1. + m_xsi[4][2] * m_xsi[4][2]), 0., m_xsi[4][0] * (2. * m_xsi[4][0] - 1.) * (-2. * m_xsi[4][2]) },
	  { 4. * m_xsi[4][1] * (-1. + m_xsi[4][2] * m_xsi[4][2]), (-4. + 4. * m_xsi[4][0] + 8. * m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]), 4. * m_xsi[4][1] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (2. * m_xsi[4][2]) },
	  { -4. * m_xsi[4][1] * (-1. + m_xsi[4][2] * m_xsi[4][2]), -4. * m_xsi[4][0] * (-1. + m_xsi[4][2] * m_xsi[4][2]), -4. * m_xsi[4][0] * m_xsi[4][1] * (2. * m_xsi[4][2]) },
	  { 0., (1. - 4. * m_xsi[4][1]) * (-1. + m_xsi[4][2] * m_xsi[4][2]), m_xsi[4][1] * (2. * m_xsi[4][1] - 1.) * (-2. * m_xsi[4][2]) },
	  { -0.5 * (3. - 4. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2], -0.5 * (3. - 4. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2], -0.5 * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. - 2. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (1. + 2. * m_xsi[4][2]) },
	  { (2. - 4. * m_xsi[4][0] - 2. * m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2], -2. * m_xsi[4][0] * (1. + m_xsi[4][2]) * m_xsi[4][2], -2. * m_xsi[4][0] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. + 2. * m_xsi[4][2]) },
	  { (-0.5 + 2. * m_xsi[4][0]) * (1. + m_xsi[4][2]) * m_xsi[4][2], 0., 0.5 * m_xsi[4][0] * (2. * m_xsi[4][0] - 1.) * (1. + 2. * m_xsi[4][2]) },
	  { -2. * m_xsi[4][1] * (1. + m_xsi[4][2]) * m_xsi[4][2], (2. - 2. * m_xsi[4][0] - 4. * m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2], -2. * m_xsi[4][1] * (-1. + m_xsi[4][0] + m_xsi[4][1]) * (1. + 2. * m_xsi[4][2]) },
	  { 2. * m_xsi[4][1] * m_xsi[4][2] * (1. + m_xsi[4][2]), 2. * m_xsi[4][0] * m_xsi[4][2] * (1. + m_xsi[4][2]), 2. * m_xsi[4][0] * m_xsi[4][1] * (1. + 2. * m_xsi[4][2]) },
	  { 0., (-0.5 + 2. * m_xsi[4][1]) * (1. + m_xsi[4][2]) * m_xsi[4][2], 0.5 * m_xsi[4][1] * (2. * m_xsi[4][1] - 1.) * (1. + 2. * m_xsi[4][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2], 0.5 * (3. - 4. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2], 0.5 * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (1. - 2. * m_xsi[5][2]) },
	  { (-2. + 4. * m_xsi[5][0] + 2. * m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2], 2. * m_xsi[5][0] * (1. - m_xsi[5][2]) * m_xsi[5][2], 2. * m_xsi[5][0] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][2]) },
	  { (0.5 - 2. * m_xsi[5][0]) * (1. - m_xsi[5][2]) * m_xsi[5][2], 0., -0.5 * m_xsi[5][0] * (2. * m_xsi[5][0] - 1.) * (1. - 2. * m_xsi[5][2]) },
	  { 2. * m_xsi[5][1] * (1. - m_xsi[5][2]) * m_xsi[5][2], (-2. + 2. * m_xsi[5][0] + 4. * m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2], 2. * m_xsi[5][1] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][2]) },
	  { -2. * m_xsi[5][1] * m_xsi[5][2] * (1. - m_xsi[5][2]), -2. * m_xsi[5][0] * m_xsi[5][2] * (1. - m_xsi[5][2]), -2. * m_xsi[5][0] * m_xsi[5][1] * (1. - 2. * m_xsi[5][2]) },
	  { 0., (0.5 - 2. * m_xsi[5][1]) * (1. - m_xsi[5][2]) * m_xsi[5][2], -0.5 * m_xsi[5][1] * (2. * m_xsi[5][1] - 1.) * (1. - 2. * m_xsi[5][2]) },
	  { (3. - 4. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]), (3. - 4. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]), (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (2. * m_xsi[5][2]) },
	  { (-4. + 8. * m_xsi[5][0] + 4. * m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]), 4. * m_xsi[5][0] * (-1. + m_xsi[5][2] * m_xsi[5][2]), 4. * m_xsi[5][0] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (2. * m_xsi[5][2]) },
	  { (1. - 4. * m_xsi[5][0]) * (-1. + m_xsi[5][2] * m_xsi[5][2]), 0., m_xsi[5][0] * (2. * m_xsi[5][0] - 1.) * (-2. * m_xsi[5][2]) },
	  { 4. * m_xsi[5][1] * (-1. + m_xsi[5][2] * m_xsi[5][2]), (-4. + 4. * m_xsi[5][0] + 8. * m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]), 4. * m_xsi[5][1] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (2. * m_xsi[5][2]) },
	  { -4. * m_xsi[5][1] * (-1. + m_xsi[5][2] * m_xsi[5][2]), -4. * m_xsi[5][0] * (-1. + m_xsi[5][2] * m_xsi[5][2]), -4. * m_xsi[5][0] * m_xsi[5][1] * (2. * m_xsi[5][2]) },
	  { 0., (1. - 4. * m_xsi[5][1]) * (-1. + m_xsi[5][2] * m_xsi[5][2]), m_xsi[5][1] * (2. * m_xsi[5][1] - 1.) * (-2. * m_xsi[5][2]) },
	  { -0.5 * (3. - 4. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2], -0.5 * (3. - 4. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2], -0.5 * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. - 2. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (1. + 2. * m_xsi[5][2]) },
	  { (2. - 4. * m_xsi[5][0] - 2. * m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2], -2. * m_xsi[5][0] * (1. + m_xsi[5][2]) * m_xsi[5][2], -2. * m_xsi[5][0] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. + 2. * m_xsi[5][2]) },
	  { (-0.5 + 2. * m_xsi[5][0]) * (1. + m_xsi[5][2]) * m_xsi[5][2], 0., 0.5 * m_xsi[5][0] * (2. * m_xsi[5][0] - 1.) * (1. + 2. * m_xsi[5][2]) },
	  { -2. * m_xsi[5][1] * (1. + m_xsi[5][2]) * m_xsi[5][2], (2. - 2. * m_xsi[5][0] - 4. * m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2], -2. * m_xsi[5][1] * (-1. + m_xsi[5][0] + m_xsi[5][1]) * (1. + 2. * m_xsi[5][2]) },
	  { 2. * m_xsi[5][1] * m_xsi[5][2] * (1. + m_xsi[5][2]), 2. * m_xsi[5][0] * m_xsi[5][2] * (1. + m_xsi[5][2]), 2. * m_xsi[5][0] * m_xsi[5][1] * (1. + 2. * m_xsi[5][2]) },
	  { 0., (-0.5 + 2. * m_xsi[5][1]) * (1. + m_xsi[5][2]) * m_xsi[5][2], 0.5 * m_xsi[5][1] * (2. * m_xsi[5][1] - 1.) * (1. + 2. * m_xsi[5][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2], 0.5 * (3. - 4. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2], 0.5 * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (1. - 2. * m_xsi[6][2]) },
	  { (-2. + 4. * m_xsi[6][0] + 2. * m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2], 2. * m_xsi[6][0] * (1. - m_xsi[6][2]) * m_xsi[6][2], 2. * m_xsi[6][0] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][2]) },
	  { (0.5 - 2. * m_xsi[6][0]) * (1. - m_xsi[6][2]) * m_xsi[6][2], 0., -0.5 * m_xsi[6][0] * (2. * m_xsi[6][0] - 1.) * (1. - 2. * m_xsi[6][2]) },
	  { 2. * m_xsi[6][1] * (1. - m_xsi[6][2]) * m_xsi[6][2], (-2. + 2. * m_xsi[6][0] + 4. * m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2], 2. * m_xsi[6][1] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][2]) },
	  { -2. * m_xsi[6][1] * m_xsi[6][2] * (1. - m_xsi[6][2]), -2. * m_xsi[6][0] * m_xsi[6][2] * (1. - m_xsi[6][2]), -2. * m_xsi[6][0] * m_xsi[6][1] * (1. - 2. * m_xsi[6][2]) },
	  { 0., (0.5 - 2. * m_xsi[6][1]) * (1. - m_xsi[6][2]) * m_xsi[6][2], -0.5 * m_xsi[6][1] * (2. * m_xsi[6][1] - 1.) * (1. - 2. * m_xsi[6][2]) },
	  { (3. - 4. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]), (3. - 4. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]), (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (2. * m_xsi[6][2]) },
	  { (-4. + 8. * m_xsi[6][0] + 4. * m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]), 4. * m_xsi[6][0] * (-1. + m_xsi[6][2] * m_xsi[6][2]), 4. * m_xsi[6][0] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (2. * m_xsi[6][2]) },
	  { (1. - 4. * m_xsi[6][0]) * (-1. + m_xsi[6][2] * m_xsi[6][2]), 0., m_xsi[6][0] * (2. * m_xsi[6][0] - 1.) * (-2. * m_xsi[6][2]) },
	  { 4. * m_xsi[6][1] * (-1. + m_xsi[6][2] * m_xsi[6][2]), (-4. + 4. * m_xsi[6][0] + 8. * m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]), 4. * m_xsi[6][1] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (2. * m_xsi[6][2]) },
	  { -4. * m_xsi[6][1] * (-1. + m_xsi[6][2] * m_xsi[6][2]), -4. * m_xsi[6][0] * (-1. + m_xsi[6][2] * m_xsi[6][2]), -4. * m_xsi[6][0] * m_xsi[6][1] * (2. * m_xsi[6][2]) },
	  { 0., (1. - 4. * m_xsi[6][1]) * (-1. + m_xsi[6][2] * m_xsi[6][2]), m_xsi[6][1] * (2. * m_xsi[6][1] - 1.) * (-2. * m_xsi[6][2]) },
	  { -0.5 * (3. - 4. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2], -0.5 * (3. - 4. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2], -0.5 * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. - 2. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (1. + 2. * m_xsi[6][2]) },
	  { (2. - 4. * m_xsi[6][0] - 2. * m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2], -2. * m_xsi[6][0] * (1. + m_xsi[6][2]) * m_xsi[6][2], -2. * m_xsi[6][0] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. + 2. * m_xsi[6][2]) },
	  { (-0.5 + 2. * m_xsi[6][0]) * (1. + m_xsi[6][2]) * m_xsi[6][2], 0., 0.5 * m_xsi[6][0] * (2. * m_xsi[6][0] - 1.) * (1. + 2. * m_xsi[6][2]) },
	  { -2. * m_xsi[6][1] * (1. + m_xsi[6][2]) * m_xsi[6][2], (2. - 2. * m_xsi[6][0] - 4. * m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2], -2. * m_xsi[6][1] * (-1. + m_xsi[6][0] + m_xsi[6][1]) * (1. + 2. * m_xsi[6][2]) },
	  { 2. * m_xsi[6][1] * m_xsi[6][2] * (1. + m_xsi[6][2]), 2. * m_xsi[6][0] * m_xsi[6][2] * (1. + m_xsi[6][2]), 2. * m_xsi[6][0] * m_xsi[6][1] * (1. + 2. * m_xsi[6][2]) },
	  { 0., (-0.5 + 2. * m_xsi[6][1]) * (1. + m_xsi[6][2]) * m_xsi[6][2], 0.5 * m_xsi[6][1] * (2. * m_xsi[6][1] - 1.) * (1. + 2. * m_xsi[6][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2], 0.5 * (3. - 4. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2], 0.5 * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (1. - 2. * m_xsi[7][2]) },
	  { (-2. + 4. * m_xsi[7][0] + 2. * m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2], 2. * m_xsi[7][0] * (1. - m_xsi[7][2]) * m_xsi[7][2], 2. * m_xsi[7][0] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][2]) },
	  { (0.5 - 2. * m_xsi[7][0]) * (1. - m_xsi[7][2]) * m_xsi[7][2], 0., -0.5 * m_xsi[7][0] * (2. * m_xsi[7][0] - 1.) * (1. - 2. * m_xsi[7][2]) },
	  { 2. * m_xsi[7][1] * (1. - m_xsi[7][2]) * m_xsi[7][2], (-2. + 2. * m_xsi[7][0] + 4. * m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2], 2. * m_xsi[7][1] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][2]) },
	  { -2. * m_xsi[7][1] * m_xsi[7][2] * (1. - m_xsi[7][2]), -2. * m_xsi[7][0] * m_xsi[7][2] * (1. - m_xsi[7][2]), -2. * m_xsi[7][0] * m_xsi[7][1] * (1. - 2. * m_xsi[7][2]) },
	  { 0., (0.5 - 2. * m_xsi[7][1]) * (1. - m_xsi[7][2]) * m_xsi[7][2], -0.5 * m_xsi[7][1] * (2. * m_xsi[7][1] - 1.) * (1. - 2. * m_xsi[7][2]) },
	  { (3. - 4. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]), (3. - 4. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]), (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (2. * m_xsi[7][2]) },
	  { (-4. + 8. * m_xsi[7][0] + 4. * m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]), 4. * m_xsi[7][0] * (-1. + m_xsi[7][2] * m_xsi[7][2]), 4. * m_xsi[7][0] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (2. * m_xsi[7][2]) },
	  { (1. - 4. * m_xsi[7][0]) * (-1. + m_xsi[7][2] * m_xsi[7][2]), 0., m_xsi[7][0] * (2. * m_xsi[7][0] - 1.) * (-2. * m_xsi[7][2]) },
	  { 4. * m_xsi[7][1] * (-1. + m_xsi[7][2] * m_xsi[7][2]), (-4. + 4. * m_xsi[7][0] + 8. * m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]), 4. * m_xsi[7][1] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (2. * m_xsi[7][2]) },
	  { -4. * m_xsi[7][1] * (-1. + m_xsi[7][2] * m_xsi[7][2]), -4. * m_xsi[7][0] * (-1. + m_xsi[7][2] * m_xsi[7][2]), -4. * m_xsi[7][0] * m_xsi[7][1] * (2. * m_xsi[7][2]) },
	  { 0., (1. - 4. * m_xsi[7][1]) * (-1. + m_xsi[7][2] * m_xsi[7][2]), m_xsi[7][1] * (2. * m_xsi[7][1] - 1.) * (-2. * m_xsi[7][2]) },
	  { -0.5 * (3. - 4. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2], -0.5 * (3. - 4. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2], -0.5 * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. - 2. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (1. + 2. * m_xsi[7][2]) },
	  { (2. - 4. * m_xsi[7][0] - 2. * m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2], -2. * m_xsi[7][0] * (1. + m_xsi[7][2]) * m_xsi[7][2], -2. * m_xsi[7][0] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. + 2. * m_xsi[7][2]) },
	  { (-0.5 + 2. * m_xsi[7][0]) * (1. + m_xsi[7][2]) * m_xsi[7][2], 0., 0.5 * m_xsi[7][0] * (2. * m_xsi[7][0] - 1.) * (1. + 2. * m_xsi[7][2]) },
	  { -2. * m_xsi[7][1] * (1. + m_xsi[7][2]) * m_xsi[7][2], (2. - 2. * m_xsi[7][0] - 4. * m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2], -2. * m_xsi[7][1] * (-1. + m_xsi[7][0] + m_xsi[7][1]) * (1. + 2. * m_xsi[7][2]) },
	  { 2. * m_xsi[7][1] * m_xsi[7][2] * (1. + m_xsi[7][2]), 2. * m_xsi[7][0] * m_xsi[7][2] * (1. + m_xsi[7][2]), 2. * m_xsi[7][0] * m_xsi[7][1] * (1. + 2. * m_xsi[7][2]) },
	  { 0., (-0.5 + 2. * m_xsi[7][1]) * (1. + m_xsi[7][2]) * m_xsi[7][2], 0.5 * m_xsi[7][1] * (2. * m_xsi[7][1] - 1.) * (1. + 2. * m_xsi[7][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2], 0.5 * (3. - 4. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2], 0.5 * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (1. - 2. * m_xsi[8][2]) },
	  { (-2. + 4. * m_xsi[8][0] + 2. * m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2], 2. * m_xsi[8][0] * (1. - m_xsi[8][2]) * m_xsi[8][2], 2. * m_xsi[8][0] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][2]) },
	  { (0.5 - 2. * m_xsi[8][0]) * (1. - m_xsi[8][2]) * m_xsi[8][2], 0., -0.5 * m_xsi[8][0] * (2. * m_xsi[8][0] - 1.) * (1. - 2. * m_xsi[8][2]) },
	  { 2. * m_xsi[8][1] * (1. - m_xsi[8][2]) * m_xsi[8][2], (-2. + 2. * m_xsi[8][0] + 4. * m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2], 2. * m_xsi[8][1] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][2]) },
	  { -2. * m_xsi[8][1] * m_xsi[8][2] * (1. - m_xsi[8][2]), -2. * m_xsi[8][0] * m_xsi[8][2] * (1. - m_xsi[8][2]), -2. * m_xsi[8][0] * m_xsi[8][1] * (1. - 2. * m_xsi[8][2]) },
	  { 0., (0.5 - 2. * m_xsi[8][1]) * (1. - m_xsi[8][2]) * m_xsi[8][2], -0.5 * m_xsi[8][1] * (2. * m_xsi[8][1] - 1.) * (1. - 2. * m_xsi[8][2]) },
	  { (3. - 4. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]), (3. - 4. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]), (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (2. * m_xsi[8][2]) },
	  { (-4. + 8. * m_xsi[8][0] + 4. * m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]), 4. * m_xsi[8][0] * (-1. + m_xsi[8][2] * m_xsi[8][2]), 4. * m_xsi[8][0] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (2. * m_xsi[8][2]) },
	  { (1. - 4. * m_xsi[8][0]) * (-1. + m_xsi[8][2] * m_xsi[8][2]), 0., m_xsi[8][0] * (2. * m_xsi[8][0] - 1.) * (-2. * m_xsi[8][2]) },
	  { 4. * m_xsi[8][1] * (-1. + m_xsi[8][2] * m_xsi[8][2]), (-4. + 4. * m_xsi[8][0] + 8. * m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]), 4. * m_xsi[8][1] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (2. * m_xsi[8][2]) },
	  { -4. * m_xsi[8][1] * (-1. + m_xsi[8][2] * m_xsi[8][2]), -4. * m_xsi[8][0] * (-1. + m_xsi[8][2] * m_xsi[8][2]), -4. * m_xsi[8][0] * m_xsi[8][1] * (2. * m_xsi[8][2]) },
	  { 0., (1. - 4. * m_xsi[8][1]) * (-1. + m_xsi[8][2] * m_xsi[8][2]), m_xsi[8][1] * (2. * m_xsi[8][1] - 1.) * (-2. * m_xsi[8][2]) },
	  { -0.5 * (3. - 4. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2], -0.5 * (3. - 4. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2], -0.5 * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. - 2. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (1. + 2. * m_xsi[8][2]) },
	  { (2. - 4. * m_xsi[8][0] - 2. * m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2], -2. * m_xsi[8][0] * (1. + m_xsi[8][2]) * m_xsi[8][2], -2. * m_xsi[8][0] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. + 2. * m_xsi[8][2]) },
	  { (-0.5 + 2. * m_xsi[8][0]) * (1. + m_xsi[8][2]) * m_xsi[8][2], 0., 0.5 * m_xsi[8][0] * (2. * m_xsi[8][0] - 1.) * (1. + 2. * m_xsi[8][2]) },
	  { -2. * m_xsi[8][1] * (1. + m_xsi[8][2]) * m_xsi[8][2], (2. - 2. * m_xsi[8][0] - 4. * m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2], -2. * m_xsi[8][1] * (-1. + m_xsi[8][0] + m_xsi[8][1]) * (1. + 2. * m_xsi[8][2]) },
	  { 2. * m_xsi[8][1] * m_xsi[8][2] * (1. + m_xsi[8][2]), 2. * m_xsi[8][0] * m_xsi[8][2] * (1. + m_xsi[8][2]), 2. * m_xsi[8][0] * m_xsi[8][1] * (1. + 2. * m_xsi[8][2]) },
	  { 0., (-0.5 + 2. * m_xsi[8][1]) * (1. + m_xsi[8][2]) * m_xsi[8][2], 0.5 * m_xsi[8][1] * (2. * m_xsi[8][1] - 1.) * (1. + 2. * m_xsi[8][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2], 0.5 * (3. - 4. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2], 0.5 * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (1. - 2. * m_xsi[9][2]) },
	  { (-2. + 4. * m_xsi[9][0] + 2. * m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2], 2. * m_xsi[9][0] * (1. - m_xsi[9][2]) * m_xsi[9][2], 2. * m_xsi[9][0] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][2]) },
	  { (0.5 - 2. * m_xsi[9][0]) * (1. - m_xsi[9][2]) * m_xsi[9][2], 0., -0.5 * m_xsi[9][0] * (2. * m_xsi[9][0] - 1.) * (1. - 2. * m_xsi[9][2]) },
	  { 2. * m_xsi[9][1] * (1. - m_xsi[9][2]) * m_xsi[9][2], (-2. + 2. * m_xsi[9][0] + 4. * m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2], 2. * m_xsi[9][1] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][2]) },
	  { -2. * m_xsi[9][1] * m_xsi[9][2] * (1. - m_xsi[9][2]), -2. * m_xsi[9][0] * m_xsi[9][2] * (1. - m_xsi[9][2]), -2. * m_xsi[9][0] * m_xsi[9][1] * (1. - 2. * m_xsi[9][2]) },
	  { 0., (0.5 - 2. * m_xsi[9][1]) * (1. - m_xsi[9][2]) * m_xsi[9][2], -0.5 * m_xsi[9][1] * (2. * m_xsi[9][1] - 1.) * (1. - 2. * m_xsi[9][2]) },
	  { (3. - 4. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]), (3. - 4. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]), (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (2. * m_xsi[9][2]) },
	  { (-4. + 8. * m_xsi[9][0] + 4. * m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]), 4. * m_xsi[9][0] * (-1. + m_xsi[9][2] * m_xsi[9][2]), 4. * m_xsi[9][0] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (2. * m_xsi[9][2]) },
	  { (1. - 4. * m_xsi[9][0]) * (-1. + m_xsi[9][2] * m_xsi[9][2]), 0., m_xsi[9][0] * (2. * m_xsi[9][0] - 1.) * (-2. * m_xsi[9][2]) },
	  { 4. * m_xsi[9][1] * (-1. + m_xsi[9][2] * m_xsi[9][2]), (-4. + 4. * m_xsi[9][0] + 8. * m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]), 4. * m_xsi[9][1] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (2. * m_xsi[9][2]) },
	  { -4. * m_xsi[9][1] * (-1. + m_xsi[9][2] * m_xsi[9][2]), -4. * m_xsi[9][0] * (-1. + m_xsi[9][2] * m_xsi[9][2]), -4. * m_xsi[9][0] * m_xsi[9][1] * (2. * m_xsi[9][2]) },
	  { 0., (1. - 4. * m_xsi[9][1]) * (-1. + m_xsi[9][2] * m_xsi[9][2]), m_xsi[9][1] * (2. * m_xsi[9][1] - 1.) * (-2. * m_xsi[9][2]) },
	  { -0.5 * (3. - 4. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2], -0.5 * (3. - 4. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2], -0.5 * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. - 2. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (1. + 2. * m_xsi[9][2]) },
	  { (2. - 4. * m_xsi[9][0] - 2. * m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2], -2. * m_xsi[9][0] * (1. + m_xsi[9][2]) * m_xsi[9][2], -2. * m_xsi[9][0] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. + 2. * m_xsi[9][2]) },
	  { (-0.5 + 2. * m_xsi[9][0]) * (1. + m_xsi[9][2]) * m_xsi[9][2], 0., 0.5 * m_xsi[9][0] * (2. * m_xsi[9][0] - 1.) * (1. + 2. * m_xsi[9][2]) },
	  { -2. * m_xsi[9][1] * (1. + m_xsi[9][2]) * m_xsi[9][2], (2. - 2. * m_xsi[9][0] - 4. * m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2], -2. * m_xsi[9][1] * (-1. + m_xsi[9][0] + m_xsi[9][1]) * (1. + 2. * m_xsi[9][2]) },
	  { 2. * m_xsi[9][1] * m_xsi[9][2] * (1. + m_xsi[9][2]), 2. * m_xsi[9][0] * m_xsi[9][2] * (1. + m_xsi[9][2]), 2. * m_xsi[9][0] * m_xsi[9][1] * (1. + 2. * m_xsi[9][2]) },
	  { 0., (-0.5 + 2. * m_xsi[9][1]) * (1. + m_xsi[9][2]) * m_xsi[9][2], 0.5 * m_xsi[9][1] * (2. * m_xsi[9][1] - 1.) * (1. + 2. * m_xsi[9][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2], 0.5 * (3. - 4. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2], 0.5 * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (1. - 2. * m_xsi[10][2]) },
	  { (-2. + 4. * m_xsi[10][0] + 2. * m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2], 2. * m_xsi[10][0] * (1. - m_xsi[10][2]) * m_xsi[10][2], 2. * m_xsi[10][0] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][2]) },
	  { (0.5 - 2. * m_xsi[10][0]) * (1. - m_xsi[10][2]) * m_xsi[10][2], 0., -0.5 * m_xsi[10][0] * (2. * m_xsi[10][0] - 1.) * (1. - 2. * m_xsi[10][2]) },
	  { 2. * m_xsi[10][1] * (1. - m_xsi[10][2]) * m_xsi[10][2], (-2. + 2. * m_xsi[10][0] + 4. * m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2], 2. * m_xsi[10][1] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][2]) },
	  { -2. * m_xsi[10][1] * m_xsi[10][2] * (1. - m_xsi[10][2]), -2. * m_xsi[10][0] * m_xsi[10][2] * (1. - m_xsi[10][2]), -2. * m_xsi[10][0] * m_xsi[10][1] * (1. - 2. * m_xsi[10][2]) },
	  { 0., (0.5 - 2. * m_xsi[10][1]) * (1. - m_xsi[10][2]) * m_xsi[10][2], -0.5 * m_xsi[10][1] * (2. * m_xsi[10][1] - 1.) * (1. - 2. * m_xsi[10][2]) },
	  { (3. - 4. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]), (3. - 4. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]), (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (2. * m_xsi[10][2]) },
	  { (-4. + 8. * m_xsi[10][0] + 4. * m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]), 4. * m_xsi[10][0] * (-1. + m_xsi[10][2] * m_xsi[10][2]), 4. * m_xsi[10][0] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (2. * m_xsi[10][2]) },
	  { (1. - 4. * m_xsi[10][0]) * (-1. + m_xsi[10][2] * m_xsi[10][2]), 0., m_xsi[10][0] * (2. * m_xsi[10][0] - 1.) * (-2. * m_xsi[10][2]) },
	  { 4. * m_xsi[10][1] * (-1. + m_xsi[10][2] * m_xsi[10][2]), (-4. + 4. * m_xsi[10][0] + 8. * m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]), 4. * m_xsi[10][1] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (2. * m_xsi[10][2]) },
	  { -4. * m_xsi[10][1] * (-1. + m_xsi[10][2] * m_xsi[10][2]), -4. * m_xsi[10][0] * (-1. + m_xsi[10][2] * m_xsi[10][2]), -4. * m_xsi[10][0] * m_xsi[10][1] * (2. * m_xsi[10][2]) },
	  { 0., (1. - 4. * m_xsi[10][1]) * (-1. + m_xsi[10][2] * m_xsi[10][2]), m_xsi[10][1] * (2. * m_xsi[10][1] - 1.) * (-2. * m_xsi[10][2]) },
	  { -0.5 * (3. - 4. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2], -0.5 * (3. - 4. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2], -0.5 * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. - 2. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (1. + 2. * m_xsi[10][2]) },
	  { (2. - 4. * m_xsi[10][0] - 2. * m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2], -2. * m_xsi[10][0] * (1. + m_xsi[10][2]) * m_xsi[10][2], -2. * m_xsi[10][0] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. + 2. * m_xsi[10][2]) },
	  { (-0.5 + 2. * m_xsi[10][0]) * (1. + m_xsi[10][2]) * m_xsi[10][2], 0., 0.5 * m_xsi[10][0] * (2. * m_xsi[10][0] - 1.) * (1. + 2. * m_xsi[10][2]) },
	  { -2. * m_xsi[10][1] * (1. + m_xsi[10][2]) * m_xsi[10][2], (2. - 2. * m_xsi[10][0] - 4. * m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2], -2. * m_xsi[10][1] * (-1. + m_xsi[10][0] + m_xsi[10][1]) * (1. + 2. * m_xsi[10][2]) },
	  { 2. * m_xsi[10][1] * m_xsi[10][2] * (1. + m_xsi[10][2]), 2. * m_xsi[10][0] * m_xsi[10][2] * (1. + m_xsi[10][2]), 2. * m_xsi[10][0] * m_xsi[10][1] * (1. + 2. * m_xsi[10][2]) },
	  { 0., (-0.5 + 2. * m_xsi[10][1]) * (1. + m_xsi[10][2]) * m_xsi[10][2], 0.5 * m_xsi[10][1] * (2. * m_xsi[10][1] - 1.) * (1. + 2. * m_xsi[10][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2], 0.5 * (3. - 4. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2], 0.5 * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (1. - 2. * m_xsi[11][2]) },
	  { (-2. + 4. * m_xsi[11][0] + 2. * m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2], 2. * m_xsi[11][0] * (1. - m_xsi[11][2]) * m_xsi[11][2], 2. * m_xsi[11][0] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][2]) },
	  { (0.5 - 2. * m_xsi[11][0]) * (1. - m_xsi[11][2]) * m_xsi[11][2], 0., -0.5 * m_xsi[11][0] * (2. * m_xsi[11][0] - 1.) * (1. - 2. * m_xsi[11][2]) },
	  { 2. * m_xsi[11][1] * (1. - m_xsi[11][2]) * m_xsi[11][2], (-2. + 2. * m_xsi[11][0] + 4. * m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2], 2. * m_xsi[11][1] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][2]) },
	  { -2. * m_xsi[11][1] * m_xsi[11][2] * (1. - m_xsi[11][2]), -2. * m_xsi[11][0] * m_xsi[11][2] * (1. - m_xsi[11][2]), -2. * m_xsi[11][0] * m_xsi[11][1] * (1. - 2. * m_xsi[11][2]) },
	  { 0., (0.5 - 2. * m_xsi[11][1]) * (1. - m_xsi[11][2]) * m_xsi[11][2], -0.5 * m_xsi[11][1] * (2. * m_xsi[11][1] - 1.) * (1. - 2. * m_xsi[11][2]) },
	  { (3. - 4. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]), (3. - 4. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]), (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (2. * m_xsi[11][2]) },
	  { (-4. + 8. * m_xsi[11][0] + 4. * m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]), 4. * m_xsi[11][0] * (-1. + m_xsi[11][2] * m_xsi[11][2]), 4. * m_xsi[11][0] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (2. * m_xsi[11][2]) },
	  { (1. - 4. * m_xsi[11][0]) * (-1. + m_xsi[11][2] * m_xsi[11][2]), 0., m_xsi[11][0] * (2. * m_xsi[11][0] - 1.) * (-2. * m_xsi[11][2]) },
	  { 4. * m_xsi[11][1] * (-1. + m_xsi[11][2] * m_xsi[11][2]), (-4. + 4. * m_xsi[11][0] + 8. * m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]), 4. * m_xsi[11][1] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (2. * m_xsi[11][2]) },
	  { -4. * m_xsi[11][1] * (-1. + m_xsi[11][2] * m_xsi[11][2]), -4. * m_xsi[11][0] * (-1. + m_xsi[11][2] * m_xsi[11][2]), -4. * m_xsi[11][0] * m_xsi[11][1] * (2. * m_xsi[11][2]) },
	  { 0., (1. - 4. * m_xsi[11][1]) * (-1. + m_xsi[11][2] * m_xsi[11][2]), m_xsi[11][1] * (2. * m_xsi[11][1] - 1.) * (-2. * m_xsi[11][2]) },
	  { -0.5 * (3. - 4. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2], -0.5 * (3. - 4. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2], -0.5 * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. - 2. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (1. + 2. * m_xsi[11][2]) },
	  { (2. - 4. * m_xsi[11][0] - 2. * m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2], -2. * m_xsi[11][0] * (1. + m_xsi[11][2]) * m_xsi[11][2], -2. * m_xsi[11][0] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. + 2. * m_xsi[11][2]) },
	  { (-0.5 + 2. * m_xsi[11][0]) * (1. + m_xsi[11][2]) * m_xsi[11][2], 0., 0.5 * m_xsi[11][0] * (2. * m_xsi[11][0] - 1.) * (1. + 2. * m_xsi[11][2]) },
	  { -2. * m_xsi[11][1] * (1. + m_xsi[11][2]) * m_xsi[11][2], (2. - 2. * m_xsi[11][0] - 4. * m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2], -2. * m_xsi[11][1] * (-1. + m_xsi[11][0] + m_xsi[11][1]) * (1. + 2. * m_xsi[11][2]) },
	  { 2. * m_xsi[11][1] * m_xsi[11][2] * (1. + m_xsi[11][2]), 2. * m_xsi[11][0] * m_xsi[11][2] * (1. + m_xsi[11][2]), 2. * m_xsi[11][0] * m_xsi[11][1] * (1. + 2. * m_xsi[11][2]) },
	  { 0., (-0.5 + 2. * m_xsi[11][1]) * (1. + m_xsi[11][2]) * m_xsi[11][2], 0.5 * m_xsi[11][1] * (2. * m_xsi[11][1] - 1.) * (1. + 2. * m_xsi[11][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2], 0.5 * (3. - 4. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2], 0.5 * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (1. - 2. * m_xsi[12][2]) },
	  { (-2. + 4. * m_xsi[12][0] + 2. * m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2], 2. * m_xsi[12][0] * (1. - m_xsi[12][2]) * m_xsi[12][2], 2. * m_xsi[12][0] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][2]) },
	  { (0.5 - 2. * m_xsi[12][0]) * (1. - m_xsi[12][2]) * m_xsi[12][2], 0., -0.5 * m_xsi[12][0] * (2. * m_xsi[12][0] - 1.) * (1. - 2. * m_xsi[12][2]) },
	  { 2. * m_xsi[12][1] * (1. - m_xsi[12][2]) * m_xsi[12][2], (-2. + 2. * m_xsi[12][0] + 4. * m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2], 2. * m_xsi[12][1] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][2]) },
	  { -2. * m_xsi[12][1] * m_xsi[12][2] * (1. - m_xsi[12][2]), -2. * m_xsi[12][0] * m_xsi[12][2] * (1. - m_xsi[12][2]), -2. * m_xsi[12][0] * m_xsi[12][1] * (1. - 2. * m_xsi[12][2]) },
	  { 0., (0.5 - 2. * m_xsi[12][1]) * (1. - m_xsi[12][2]) * m_xsi[12][2], -0.5 * m_xsi[12][1] * (2. * m_xsi[12][1] - 1.) * (1. - 2. * m_xsi[12][2]) },
	  { (3. - 4. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]), (3. - 4. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]), (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (2. * m_xsi[12][2]) },
	  { (-4. + 8. * m_xsi[12][0] + 4. * m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]), 4. * m_xsi[12][0] * (-1. + m_xsi[12][2] * m_xsi[12][2]), 4. * m_xsi[12][0] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (2. * m_xsi[12][2]) },
	  { (1. - 4. * m_xsi[12][0]) * (-1. + m_xsi[12][2] * m_xsi[12][2]), 0., m_xsi[12][0] * (2. * m_xsi[12][0] - 1.) * (-2. * m_xsi[12][2]) },
	  { 4. * m_xsi[12][1] * (-1. + m_xsi[12][2] * m_xsi[12][2]), (-4. + 4. * m_xsi[12][0] + 8. * m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]), 4. * m_xsi[12][1] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (2. * m_xsi[12][2]) },
	  { -4. * m_xsi[12][1] * (-1. + m_xsi[12][2] * m_xsi[12][2]), -4. * m_xsi[12][0] * (-1. + m_xsi[12][2] * m_xsi[12][2]), -4. * m_xsi[12][0] * m_xsi[12][1] * (2. * m_xsi[12][2]) },
	  { 0., (1. - 4. * m_xsi[12][1]) * (-1. + m_xsi[12][2] * m_xsi[12][2]), m_xsi[12][1] * (2. * m_xsi[12][1] - 1.) * (-2. * m_xsi[12][2]) },
	  { -0.5 * (3. - 4. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2], -0.5 * (3. - 4. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2], -0.5 * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. - 2. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (1. + 2. * m_xsi[12][2]) },
	  { (2. - 4. * m_xsi[12][0] - 2. * m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2], -2. * m_xsi[12][0] * (1. + m_xsi[12][2]) * m_xsi[12][2], -2. * m_xsi[12][0] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. + 2. * m_xsi[12][2]) },
	  { (-0.5 + 2. * m_xsi[12][0]) * (1. + m_xsi[12][2]) * m_xsi[12][2], 0., 0.5 * m_xsi[12][0] * (2. * m_xsi[12][0] - 1.) * (1. + 2. * m_xsi[12][2]) },
	  { -2. * m_xsi[12][1] * (1. + m_xsi[12][2]) * m_xsi[12][2], (2. - 2. * m_xsi[12][0] - 4. * m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2], -2. * m_xsi[12][1] * (-1. + m_xsi[12][0] + m_xsi[12][1]) * (1. + 2. * m_xsi[12][2]) },
	  { 2. * m_xsi[12][1] * m_xsi[12][2] * (1. + m_xsi[12][2]), 2. * m_xsi[12][0] * m_xsi[12][2] * (1. + m_xsi[12][2]), 2. * m_xsi[12][0] * m_xsi[12][1] * (1. + 2. * m_xsi[12][2]) },
	  { 0., (-0.5 + 2. * m_xsi[12][1]) * (1. + m_xsi[12][2]) * m_xsi[12][2], 0.5 * m_xsi[12][1] * (2. * m_xsi[12][1] - 1.) * (1. + 2. * m_xsi[12][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2], 0.5 * (3. - 4. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2], 0.5 * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (1. - 2. * m_xsi[13][2]) },
	  { (-2. + 4. * m_xsi[13][0] + 2. * m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2], 2. * m_xsi[13][0] * (1. - m_xsi[13][2]) * m_xsi[13][2], 2. * m_xsi[13][0] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][2]) },
	  { (0.5 - 2. * m_xsi[13][0]) * (1. - m_xsi[13][2]) * m_xsi[13][2], 0., -0.5 * m_xsi[13][0] * (2. * m_xsi[13][0] - 1.) * (1. - 2. * m_xsi[13][2]) },
	  { 2. * m_xsi[13][1] * (1. - m_xsi[13][2]) * m_xsi[13][2], (-2. + 2. * m_xsi[13][0] + 4. * m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2], 2. * m_xsi[13][1] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][2]) },
	  { -2. * m_xsi[13][1] * m_xsi[13][2] * (1. - m_xsi[13][2]), -2. * m_xsi[13][0] * m_xsi[13][2] * (1. - m_xsi[13][2]), -2. * m_xsi[13][0] * m_xsi[13][1] * (1. - 2. * m_xsi[13][2]) },
	  { 0., (0.5 - 2. * m_xsi[13][1]) * (1. - m_xsi[13][2]) * m_xsi[13][2], -0.5 * m_xsi[13][1] * (2. * m_xsi[13][1] - 1.) * (1. - 2. * m_xsi[13][2]) },
	  { (3. - 4. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]), (3. - 4. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]), (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (2. * m_xsi[13][2]) },
	  { (-4. + 8. * m_xsi[13][0] + 4. * m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]), 4. * m_xsi[13][0] * (-1. + m_xsi[13][2] * m_xsi[13][2]), 4. * m_xsi[13][0] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (2. * m_xsi[13][2]) },
	  { (1. - 4. * m_xsi[13][0]) * (-1. + m_xsi[13][2] * m_xsi[13][2]), 0., m_xsi[13][0] * (2. * m_xsi[13][0] - 1.) * (-2. * m_xsi[13][2]) },
	  { 4. * m_xsi[13][1] * (-1. + m_xsi[13][2] * m_xsi[13][2]), (-4. + 4. * m_xsi[13][0] + 8. * m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]), 4. * m_xsi[13][1] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (2. * m_xsi[13][2]) },
	  { -4. * m_xsi[13][1] * (-1. + m_xsi[13][2] * m_xsi[13][2]), -4. * m_xsi[13][0] * (-1. + m_xsi[13][2] * m_xsi[13][2]), -4. * m_xsi[13][0] * m_xsi[13][1] * (2. * m_xsi[13][2]) },
	  { 0., (1. - 4. * m_xsi[13][1]) * (-1. + m_xsi[13][2] * m_xsi[13][2]), m_xsi[13][1] * (2. * m_xsi[13][1] - 1.) * (-2. * m_xsi[13][2]) },
	  { -0.5 * (3. - 4. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2], -0.5 * (3. - 4. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2], -0.5 * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. - 2. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (1. + 2. * m_xsi[13][2]) },
	  { (2. - 4. * m_xsi[13][0] - 2. * m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2], -2. * m_xsi[13][0] * (1. + m_xsi[13][2]) * m_xsi[13][2], -2. * m_xsi[13][0] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. + 2. * m_xsi[13][2]) },
	  { (-0.5 + 2. * m_xsi[13][0]) * (1. + m_xsi[13][2]) * m_xsi[13][2], 0., 0.5 * m_xsi[13][0] * (2. * m_xsi[13][0] - 1.) * (1. + 2. * m_xsi[13][2]) },
	  { -2. * m_xsi[13][1] * (1. + m_xsi[13][2]) * m_xsi[13][2], (2. - 2. * m_xsi[13][0] - 4. * m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2], -2. * m_xsi[13][1] * (-1. + m_xsi[13][0] + m_xsi[13][1]) * (1. + 2. * m_xsi[13][2]) },
	  { 2. * m_xsi[13][1] * m_xsi[13][2] * (1. + m_xsi[13][2]), 2. * m_xsi[13][0] * m_xsi[13][2] * (1. + m_xsi[13][2]), 2. * m_xsi[13][0] * m_xsi[13][1] * (1. + 2. * m_xsi[13][2]) },
	  { 0., (-0.5 + 2. * m_xsi[13][1]) * (1. + m_xsi[13][2]) * m_xsi[13][2], 0.5 * m_xsi[13][1] * (2. * m_xsi[13][1] - 1.) * (1. + 2. * m_xsi[13][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2], 0.5 * (3. - 4. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2], 0.5 * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (1. - 2. * m_xsi[14][2]) },
	  { (-2. + 4. * m_xsi[14][0] + 2. * m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2], 2. * m_xsi[14][0] * (1. - m_xsi[14][2]) * m_xsi[14][2], 2. * m_xsi[14][0] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][2]) },
	  { (0.5 - 2. * m_xsi[14][0]) * (1. - m_xsi[14][2]) * m_xsi[14][2], 0., -0.5 * m_xsi[14][0] * (2. * m_xsi[14][0] - 1.) * (1. - 2. * m_xsi[14][2]) },
	  { 2. * m_xsi[14][1] * (1. - m_xsi[14][2]) * m_xsi[14][2], (-2. + 2. * m_xsi[14][0] + 4. * m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2], 2. * m_xsi[14][1] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][2]) },
	  { -2. * m_xsi[14][1] * m_xsi[14][2] * (1. - m_xsi[14][2]), -2. * m_xsi[14][0] * m_xsi[14][2] * (1. - m_xsi[14][2]), -2. * m_xsi[14][0] * m_xsi[14][1] * (1. - 2. * m_xsi[14][2]) },
	  { 0., (0.5 - 2. * m_xsi[14][1]) * (1. - m_xsi[14][2]) * m_xsi[14][2], -0.5 * m_xsi[14][1] * (2. * m_xsi[14][1] - 1.) * (1. - 2. * m_xsi[14][2]) },
	  { (3. - 4. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]), (3. - 4. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]), (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (2. * m_xsi[14][2]) },
	  { (-4. + 8. * m_xsi[14][0] + 4. * m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]), 4. * m_xsi[14][0] * (-1. + m_xsi[14][2] * m_xsi[14][2]), 4. * m_xsi[14][0] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (2. * m_xsi[14][2]) },
	  { (1. - 4. * m_xsi[14][0]) * (-1. + m_xsi[14][2] * m_xsi[14][2]), 0., m_xsi[14][0] * (2. * m_xsi[14][0] - 1.) * (-2. * m_xsi[14][2]) },
	  { 4. * m_xsi[14][1] * (-1. + m_xsi[14][2] * m_xsi[14][2]), (-4. + 4. * m_xsi[14][0] + 8. * m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]), 4. * m_xsi[14][1] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (2. * m_xsi[14][2]) },
	  { -4. * m_xsi[14][1] * (-1. + m_xsi[14][2] * m_xsi[14][2]), -4. * m_xsi[14][0] * (-1. + m_xsi[14][2] * m_xsi[14][2]), -4. * m_xsi[14][0] * m_xsi[14][1] * (2. * m_xsi[14][2]) },
	  { 0., (1. - 4. * m_xsi[14][1]) * (-1. + m_xsi[14][2] * m_xsi[14][2]), m_xsi[14][1] * (2. * m_xsi[14][1] - 1.) * (-2. * m_xsi[14][2]) },
	  { -0.5 * (3. - 4. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2], -0.5 * (3. - 4. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2], -0.5 * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. - 2. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (1. + 2. * m_xsi[14][2]) },
	  { (2. - 4. * m_xsi[14][0] - 2. * m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2], -2. * m_xsi[14][0] * (1. + m_xsi[14][2]) * m_xsi[14][2], -2. * m_xsi[14][0] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. + 2. * m_xsi[14][2]) },
	  { (-0.5 + 2. * m_xsi[14][0]) * (1. + m_xsi[14][2]) * m_xsi[14][2], 0., 0.5 * m_xsi[14][0] * (2. * m_xsi[14][0] - 1.) * (1. + 2. * m_xsi[14][2]) },
	  { -2. * m_xsi[14][1] * (1. + m_xsi[14][2]) * m_xsi[14][2], (2. - 2. * m_xsi[14][0] - 4. * m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2], -2. * m_xsi[14][1] * (-1. + m_xsi[14][0] + m_xsi[14][1]) * (1. + 2. * m_xsi[14][2]) },
	  { 2. * m_xsi[14][1] * m_xsi[14][2] * (1. + m_xsi[14][2]), 2. * m_xsi[14][0] * m_xsi[14][2] * (1. + m_xsi[14][2]), 2. * m_xsi[14][0] * m_xsi[14][1] * (1. + 2. * m_xsi[14][2]) },
	  { 0., (-0.5 + 2. * m_xsi[14][1]) * (1. + m_xsi[14][2]) * m_xsi[14][2], 0.5 * m_xsi[14][1] * (2. * m_xsi[14][1] - 1.) * (1. + 2. * m_xsi[14][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2], 0.5 * (3. - 4. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2], 0.5 * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (1. - 2. * m_xsi[15][2]) },
	  { (-2. + 4. * m_xsi[15][0] + 2. * m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2], 2. * m_xsi[15][0] * (1. - m_xsi[15][2]) * m_xsi[15][2], 2. * m_xsi[15][0] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][2]) },
	  { (0.5 - 2. * m_xsi[15][0]) * (1. - m_xsi[15][2]) * m_xsi[15][2], 0., -0.5 * m_xsi[15][0] * (2. * m_xsi[15][0] - 1.) * (1. - 2. * m_xsi[15][2]) },
	  { 2. * m_xsi[15][1] * (1. - m_xsi[15][2]) * m_xsi[15][2], (-2. + 2. * m_xsi[15][0] + 4. * m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2], 2. * m_xsi[15][1] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][2]) },
	  { -2. * m_xsi[15][1] * m_xsi[15][2] * (1. - m_xsi[15][2]), -2. * m_xsi[15][0] * m_xsi[15][2] * (1. - m_xsi[15][2]), -2. * m_xsi[15][0] * m_xsi[15][1] * (1. - 2. * m_xsi[15][2]) },
	  { 0., (0.5 - 2. * m_xsi[15][1]) * (1. - m_xsi[15][2]) * m_xsi[15][2], -0.5 * m_xsi[15][1] * (2. * m_xsi[15][1] - 1.) * (1. - 2. * m_xsi[15][2]) },
	  { (3. - 4. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]), (3. - 4. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]), (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (2. * m_xsi[15][2]) },
	  { (-4. + 8. * m_xsi[15][0] + 4. * m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]), 4. * m_xsi[15][0] * (-1. + m_xsi[15][2] * m_xsi[15][2]), 4. * m_xsi[15][0] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (2. * m_xsi[15][2]) },
	  { (1. - 4. * m_xsi[15][0]) * (-1. + m_xsi[15][2] * m_xsi[15][2]), 0., m_xsi[15][0] * (2. * m_xsi[15][0] - 1.) * (-2. * m_xsi[15][2]) },
	  { 4. * m_xsi[15][1] * (-1. + m_xsi[15][2] * m_xsi[15][2]), (-4. + 4. * m_xsi[15][0] + 8. * m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]), 4. * m_xsi[15][1] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (2. * m_xsi[15][2]) },
	  { -4. * m_xsi[15][1] * (-1. + m_xsi[15][2] * m_xsi[15][2]), -4. * m_xsi[15][0] * (-1. + m_xsi[15][2] * m_xsi[15][2]), -4. * m_xsi[15][0] * m_xsi[15][1] * (2. * m_xsi[15][2]) },
	  { 0., (1. - 4. * m_xsi[15][1]) * (-1. + m_xsi[15][2] * m_xsi[15][2]), m_xsi[15][1] * (2. * m_xsi[15][1] - 1.) * (-2. * m_xsi[15][2]) },
	  { -0.5 * (3. - 4. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2], -0.5 * (3. - 4. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2], -0.5 * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. - 2. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (1. + 2. * m_xsi[15][2]) },
	  { (2. - 4. * m_xsi[15][0] - 2. * m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2], -2. * m_xsi[15][0] * (1. + m_xsi[15][2]) * m_xsi[15][2], -2. * m_xsi[15][0] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. + 2. * m_xsi[15][2]) },
	  { (-0.5 + 2. * m_xsi[15][0]) * (1. + m_xsi[15][2]) * m_xsi[15][2], 0., 0.5 * m_xsi[15][0] * (2. * m_xsi[15][0] - 1.) * (1. + 2. * m_xsi[15][2]) },
	  { -2. * m_xsi[15][1] * (1. + m_xsi[15][2]) * m_xsi[15][2], (2. - 2. * m_xsi[15][0] - 4. * m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2], -2. * m_xsi[15][1] * (-1. + m_xsi[15][0] + m_xsi[15][1]) * (1. + 2. * m_xsi[15][2]) },
	  { 2. * m_xsi[15][1] * m_xsi[15][2] * (1. + m_xsi[15][2]), 2. * m_xsi[15][0] * m_xsi[15][2] * (1. + m_xsi[15][2]), 2. * m_xsi[15][0] * m_xsi[15][1] * (1. + 2. * m_xsi[15][2]) },
	  { 0., (-0.5 + 2. * m_xsi[15][1]) * (1. + m_xsi[15][2]) * m_xsi[15][2], 0.5 * m_xsi[15][1] * (2. * m_xsi[15][1] - 1.) * (1. + 2. * m_xsi[15][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2], 0.5 * (3. - 4. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2], 0.5 * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (1. - 2. * m_xsi[16][2]) },
	  { (-2. + 4. * m_xsi[16][0] + 2. * m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2], 2. * m_xsi[16][0] * (1. - m_xsi[16][2]) * m_xsi[16][2], 2. * m_xsi[16][0] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][2]) },
	  { (0.5 - 2. * m_xsi[16][0]) * (1. - m_xsi[16][2]) * m_xsi[16][2], 0., -0.5 * m_xsi[16][0] * (2. * m_xsi[16][0] - 1.) * (1. - 2. * m_xsi[16][2]) },
	  { 2. * m_xsi[16][1] * (1. - m_xsi[16][2]) * m_xsi[16][2], (-2. + 2. * m_xsi[16][0] + 4. * m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2], 2. * m_xsi[16][1] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][2]) },
	  { -2. * m_xsi[16][1] * m_xsi[16][2] * (1. - m_xsi[16][2]), -2. * m_xsi[16][0] * m_xsi[16][2] * (1. - m_xsi[16][2]), -2. * m_xsi[16][0] * m_xsi[16][1] * (1. - 2. * m_xsi[16][2]) },
	  { 0., (0.5 - 2. * m_xsi[16][1]) * (1. - m_xsi[16][2]) * m_xsi[16][2], -0.5 * m_xsi[16][1] * (2. * m_xsi[16][1] - 1.) * (1. - 2. * m_xsi[16][2]) },
	  { (3. - 4. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]), (3. - 4. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]), (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (2. * m_xsi[16][2]) },
	  { (-4. + 8. * m_xsi[16][0] + 4. * m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]), 4. * m_xsi[16][0] * (-1. + m_xsi[16][2] * m_xsi[16][2]), 4. * m_xsi[16][0] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (2. * m_xsi[16][2]) },
	  { (1. - 4. * m_xsi[16][0]) * (-1. + m_xsi[16][2] * m_xsi[16][2]), 0., m_xsi[16][0] * (2. * m_xsi[16][0] - 1.) * (-2. * m_xsi[16][2]) },
	  { 4. * m_xsi[16][1] * (-1. + m_xsi[16][2] * m_xsi[16][2]), (-4. + 4. * m_xsi[16][0] + 8. * m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]), 4. * m_xsi[16][1] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (2. * m_xsi[16][2]) },
	  { -4. * m_xsi[16][1] * (-1. + m_xsi[16][2] * m_xsi[16][2]), -4. * m_xsi[16][0] * (-1. + m_xsi[16][2] * m_xsi[16][2]), -4. * m_xsi[16][0] * m_xsi[16][1] * (2. * m_xsi[16][2]) },
	  { 0., (1. - 4. * m_xsi[16][1]) * (-1. + m_xsi[16][2] * m_xsi[16][2]), m_xsi[16][1] * (2. * m_xsi[16][1] - 1.) * (-2. * m_xsi[16][2]) },
	  { -0.5 * (3. - 4. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2], -0.5 * (3. - 4. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2], -0.5 * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. - 2. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (1. + 2. * m_xsi[16][2]) },
	  { (2. - 4. * m_xsi[16][0] - 2. * m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2], -2. * m_xsi[16][0] * (1. + m_xsi[16][2]) * m_xsi[16][2], -2. * m_xsi[16][0] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. + 2. * m_xsi[16][2]) },
	  { (-0.5 + 2. * m_xsi[16][0]) * (1. + m_xsi[16][2]) * m_xsi[16][2], 0., 0.5 * m_xsi[16][0] * (2. * m_xsi[16][0] - 1.) * (1. + 2. * m_xsi[16][2]) },
	  { -2. * m_xsi[16][1] * (1. + m_xsi[16][2]) * m_xsi[16][2], (2. - 2. * m_xsi[16][0] - 4. * m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2], -2. * m_xsi[16][1] * (-1. + m_xsi[16][0] + m_xsi[16][1]) * (1. + 2. * m_xsi[16][2]) },
	  { 2. * m_xsi[16][1] * m_xsi[16][2] * (1. + m_xsi[16][2]), 2. * m_xsi[16][0] * m_xsi[16][2] * (1. + m_xsi[16][2]), 2. * m_xsi[16][0] * m_xsi[16][1] * (1. + 2. * m_xsi[16][2]) },
	  { 0., (-0.5 + 2. * m_xsi[16][1]) * (1. + m_xsi[16][2]) * m_xsi[16][2], 0.5 * m_xsi[16][1] * (2. * m_xsi[16][1] - 1.) * (1. + 2. * m_xsi[16][2]) } },

	{ { 0.5 * (3. - 4. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2], 0.5 * (3. - 4. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2], 0.5 * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (1. - 2. * m_xsi[17][2]) },
	  { (-2. + 4. * m_xsi[17][0] + 2. * m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2], 2. * m_xsi[17][0] * (1. - m_xsi[17][2]) * m_xsi[17][2], 2. * m_xsi[17][0] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][2]) },
	  { (0.5 - 2. * m_xsi[17][0]) * (1. - m_xsi[17][2]) * m_xsi[17][2], 0., -0.5 * m_xsi[17][0] * (2. * m_xsi[17][0] - 1.) * (1. - 2. * m_xsi[17][2]) },
	  { 2. * m_xsi[17][1] * (1. - m_xsi[17][2]) * m_xsi[17][2], (-2. + 2. * m_xsi[17][0] + 4. * m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2], 2. * m_xsi[17][1] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][2]) },
	  { -2. * m_xsi[17][1] * m_xsi[17][2] * (1. - m_xsi[17][2]), -2. * m_xsi[17][0] * m_xsi[17][2] * (1. - m_xsi[17][2]), -2. * m_xsi[17][0] * m_xsi[17][1] * (1. - 2. * m_xsi[17][2]) },
	  { 0., (0.5 - 2. * m_xsi[17][1]) * (1. - m_xsi[17][2]) * m_xsi[17][2], -0.5 * m_xsi[17][1] * (2. * m_xsi[17][1] - 1.) * (1. - 2. * m_xsi[17][2]) },
	  { (3. - 4. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]), (3. - 4. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]), (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (2. * m_xsi[17][2]) },
	  { (-4. + 8. * m_xsi[17][0] + 4. * m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]), 4. * m_xsi[17][0] * (-1. + m_xsi[17][2] * m_xsi[17][2]), 4. * m_xsi[17][0] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (2. * m_xsi[17][2]) },
	  { (1. - 4. * m_xsi[17][0]) * (-1. + m_xsi[17][2] * m_xsi[17][2]), 0., m_xsi[17][0] * (2. * m_xsi[17][0] - 1.) * (-2. * m_xsi[17][2]) },
	  { 4. * m_xsi[17][1] * (-1. + m_xsi[17][2] * m_xsi[17][2]), (-4. + 4. * m_xsi[17][0] + 8. * m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]), 4. * m_xsi[17][1] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (2. * m_xsi[17][2]) },
	  { -4. * m_xsi[17][1] * (-1. + m_xsi[17][2] * m_xsi[17][2]), -4. * m_xsi[17][0] * (-1. + m_xsi[17][2] * m_xsi[17][2]), -4. * m_xsi[17][0] * m_xsi[17][1] * (2. * m_xsi[17][2]) },
	  { 0., (1. - 4. * m_xsi[17][1]) * (-1. + m_xsi[17][2] * m_xsi[17][2]), m_xsi[17][1] * (2. * m_xsi[17][1] - 1.) * (-2. * m_xsi[17][2]) },
	  { -0.5 * (3. - 4. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2], -0.5 * (3. - 4. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2], -0.5 * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. - 2. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (1. + 2. * m_xsi[17][2]) },
	  { (2. - 4. * m_xsi[17][0] - 2. * m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2], -2. * m_xsi[17][0] * (1. + m_xsi[17][2]) * m_xsi[17][2], -2. * m_xsi[17][0] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. + 2. * m_xsi[17][2]) },
	  { (-0.5 + 2. * m_xsi[17][0]) * (1. + m_xsi[17][2]) * m_xsi[17][2], 0., 0.5 * m_xsi[17][0] * (2. * m_xsi[17][0] - 1.) * (1. + 2. * m_xsi[17][2]) },
	  { -2. * m_xsi[17][1] * (1. + m_xsi[17][2]) * m_xsi[17][2], (2. - 2. * m_xsi[17][0] - 4. * m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2], -2. * m_xsi[17][1] * (-1. + m_xsi[17][0] + m_xsi[17][1]) * (1. + 2. * m_xsi[17][2]) },
	  { 2. * m_xsi[17][1] * m_xsi[17][2] * (1. + m_xsi[17][2]), 2. * m_xsi[17][0] * m_xsi[17][2] * (1. + m_xsi[17][2]), 2. * m_xsi[17][0] * m_xsi[17][1] * (1. + 2. * m_xsi[17][2]) },
	  { 0., (-0.5 + 2. * m_xsi[17][1]) * (1. + m_xsi[17][2]) * m_xsi[17][2], 0.5 * m_xsi[17][1] * (2. * m_xsi[17][1] - 1.) * (1. + 2. * m_xsi[17][2]) } } };
