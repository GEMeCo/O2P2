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
// Plane element, with cubic interpolation functions, triangular shaped.
// 
// ================================================================================================
#pragma once

// Custom Header Files
#include "Element.h"
#include "IntegrationPoint.h"

/** @ingroup Elements
  * @class Elem_Tri10
  *
  * @brief Triangular cubic element with 10 nodes.
  * @details Plane element, with cubic interpolation functions, triangular shaped.
  * Options for integration points: 6, 7, 12 and 13.
  * @image html Elem_Tri10.png height=300
  * @note Minimum number of integration points: 6.
  */
class Elem_Tri10 : public ElementPlane
{
private:
	Elem_Tri10() = delete;

protected:
	/** Constructor for triangular cubic elements.
	  * @param Material Pointer to Material class.
	  * @param Section Pointer to Section class.
	  */
	Elem_Tri10(std::shared_ptr<Material>& Material, std::shared_ptr<Section>& Section)
		: ElementPlane(Material, Section) { };

public:
	// Output function for AcadView, based on element index.
	const std::string printByIndex_AV(const size_t add) const override {
		std::stringstream msg;
		msg << "2 3 " << this->v_Conect[0]->m_index + add << " " << this->v_Conect[1]->m_index + add << " "
			<< this->v_Conect[2]->m_index + add << " " << this->v_Conect[3]->m_index + add << " "
			<< this->v_Conect[4]->m_index + add << " " << this->v_Conect[5]->m_index + add << " "
			<< this->v_Conect[6]->m_index + add << " " << this->v_Conect[7]->m_index + add << " "
			<< this->v_Conect[8]->m_index + add << " " << this->v_Conect[9]->m_index + add << " "
			<< this->m_Mat->m_index << "\n";
		return msg.str();
	};

	// Output function for AcadView, based on element node number.
	const std::string printByAdder_AV(const size_t add) const override {
		std::stringstream msg;
		msg << "2 3 " << (1 + add) << " " << (2 + add) << " " << (3 + add) << " "
			<< (4 + add) << " " << (5 + add) << " " << (6 + add) << " "
			<< (7 + add) << " " << (8 + add) << " " << (9 + add) << " "
			<< (10 + add) << " " << this->m_Mat->m_index << "\n";
		return msg.str();
	};

	// Evaluates shape function in the point.
	Eigen::VectorXd getShapeFcOnPoint(const double* Point) override;

	// Evaluates the derivative of shape function in the point.
	Eigen::MatrixXd getShapeDerivOnPoint(const double* Point) override;

	// Returns the number of nodes of current element.
	int getNumNodes() override { return m_NumNodes; };

	// Returns the number of faces of current element.
	int getNumFaces() override { return m_NumFaces; };

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
	};

protected:
	// Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity.
	void setGeomProperties() override;

protected:
	/** @brief Number of Nodes */
	static const int m_NumNodes{ 10 };

	/** @brief Number of Faces */
	static const int m_NumFaces{ 1 };
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
	Elem_Tri10_IP() = delete;

public:
	/** Constructor for triangular cubic elements.
	  * @param Material Pointer to Material class.
	  * @param Section Pointer to Section class.
	  */
	explicit Elem_Tri10_IP(std::shared_ptr<Material>& Material, std::shared_ptr<Section>& Section)
		: Elem_Tri10(Material, Section) { };

	// Return a vector with values on the integration points currently known in the element' nodes.
	Eigen::VectorXd getValueOnIPs(const double* value) override;

	// Returns a pointer to the first element of the shape functions (with size [nIP][m_NumNodes]).
	double const* getShapeFc() const override { return &m_Psi[0][0]; };

	// Returns a pointer to the first element of the derivative of shape functions (with size [nIP][m_NumNodes][m_Dim]).
	double const* getShapeDerivative() const override { return &m_DPsi[0][0][0]; };

	// Returns a pointer to the weight of the integation points (with size [nIP]).
	double const* getWeight() const override { return m_weight; };

	// Returns the number of integration points of current element.
	int getNumIP() override { return nIP; };

private:
	// Weights for numerical integration
	static const double* m_weight;

	// Shape functions
	static const double m_Psi[nIP][m_NumNodes];

	// Shape functions derivative
	static const double m_DPsi[nIP][m_NumNodes][m_Dim];
};


// ================================================================================================
//
// Implementation of Member Function: getShapeFcOnPoint
// Shape functions evaluated on Point
// 
// ================================================================================================
inline Eigen::VectorXd Elem_Tri10::getShapeFcOnPoint(const double* Point) {
	Eigen::VectorXd Psi(10);

	Psi(0) = 1. - 5.5 * Point[0] - 5.5 * Point[1] + 9. * Point[0] * Point[0] + 18. * Point[0] * Point[1] + 9. * Point[1] * Point[1] - 4.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 4.5 * Point[1] * Point[1] * Point[1];
	Psi(1) = 9. * Point[0] - 22.5 * Point[0] * Point[0] - 22.5 * Point[0] * Point[1] + 13.5 * Point[0] * Point[0] * Point[0] + 27. * Point[0] * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1];
	Psi(2) = -4.5 * Point[0] + 18. * Point[0] * Point[0] + 4.5 * Point[0] * Point[1] - 13.5 * Point[0] * Point[0] * Point[0] - 13.5 * Point[0] * Point[0] * Point[1];
	Psi(3) = Point[0] - 4.5 * Point[0] * Point[0] + 4.5 * Point[0] * Point[0] * Point[0];
	Psi(4) = 9. * Point[1] - 22.5 * Point[0] * Point[1] - 22.5 * Point[1] * Point[1] + 13.5 * Point[0] * Point[0] * Point[1] + 27. * Point[0] * Point[1] * Point[1] + 13.5 * Point[1] * Point[1] * Point[1];
	Psi(5) = 27. * Point[0] * Point[1] - 27. * Point[0] * Point[0] * Point[1] - 27. * Point[0] * Point[1] * Point[1];
	Psi(6) = -4.5 * Point[0] * Point[1] + 13.5 * Point[0] * Point[0] * Point[1];
	Psi(7) = -4.5 * Point[1] + 18. * Point[1] * Point[1] + 4.5 * Point[0] * Point[1] - 13.5 * Point[0] * Point[1] * Point[1] - 13.5 * Point[1] * Point[1] * Point[1];
	Psi(8) = -4.5 * Point[0] * Point[1] + 13.5 * Point[0] * Point[1] * Point[1];
	Psi(9) = Point[1] - 4.5 * Point[1] * Point[1] + 4.5 * Point[1] * Point[1] * Point[1];

	return Psi;
};

// ================================================================================================
//
// Implementation of Member Function: getShapeDerivOnPoint
// Shape functions derivative evaluated on Point
// 
// ================================================================================================
inline Eigen::MatrixXd Elem_Tri10::getShapeDerivOnPoint(const double* Point) {
	Eigen::MatrixXd DPsi(10, 2);

	DPsi(0, 0) = -5.5 + 18. * Point[0] + 18. * Point[1] - 13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1];
	DPsi(1, 0) = 9. - 45. * Point[0] - 22.5 * Point[1] + 40.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 13.5 * Point[1] * Point[1];
	DPsi(2, 0) = -4.5 + 36. * Point[0] + 4.5 * Point[1] - 40.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1];
	DPsi(3, 0) = 1. - 9. * Point[0] + 13.5 * Point[0] * Point[0];
	DPsi(4, 0) = -22.5 * Point[1] + 27. * Point[0] * Point[1] + 27. * Point[1] * Point[1];
	DPsi(5, 0) = 27. * Point[1] - 54. * Point[0] * Point[1] - 27. * Point[1] * Point[1];
	DPsi(6, 0) = -4.5 * Point[1] + 27. * Point[0] * Point[1];
	DPsi(7, 0) = 4.5 * Point[1] - 13.5 * Point[1] * Point[1];
	DPsi(8, 0) = -4.5 * Point[1] + 13.5 * Point[1] * Point[1];
	DPsi(9, 0) = 0.;

	DPsi(0, 1) = -5.5 + 18. * Point[0] + 18. * Point[1] - 13.5 * Point[0] * Point[0] - 27. * Point[0] * Point[1] - 13.5 * Point[1] * Point[1];
	DPsi(1, 1) = -22.5 * Point[0] + 27. * Point[0] * Point[0] + 27. * Point[0] * Point[1];
	DPsi(2, 1) = 4.5 * Point[0] - 13.5 * Point[0] * Point[0];
	DPsi(3, 1) = 0.;
	DPsi(4, 1) = 9. - 22.5 * Point[0] - 45. * Point[1] + 13.5 * Point[0] * Point[0] + 54. * Point[0] * Point[1] + 40.5 * Point[1] * Point[1];
	DPsi(5, 1) = 27. * Point[0] - 27. * Point[0] * Point[0] - 54. * Point[0] * Point[1];
	DPsi(6, 1) = -4.5 * Point[0] + 13.5 * Point[0] * Point[0];
	DPsi(7, 1) = -4.5 + 36. * Point[1] + 4.5 * Point[0] - 27. * Point[0] * Point[1] - 40.5 * Point[1] * Point[1];
	DPsi(8, 1) = -4.5 * Point[0] + 27. * Point[0] * Point[1];
	DPsi(9, 1) = 1. - 9. * Point[1] + 13.5 * Point[1] * Point[1];

	return DPsi;
};

// ================================================================================================
//
// Implementation of Member Function: setGeomProperties
// Evaluate initial properties
// 
// ================================================================================================
inline void Elem_Tri10::setGeomProperties() {

	const int nVertices = 3;

	// Allocate an array with size m_Dim to which m_Centroid points to.
	m_Centroid = std::make_unique<double[]>(m_Dim);

	// Create a temporary array with the vertices of the polygon
	std::array<Node<m_Dim>*, nVertices> vertices;
	vertices[0] = v_Conect[0].get();
	vertices[1] = v_Conect[3].get();
	vertices[2] = v_Conect[9].get();

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
inline Eigen::VectorXd Elem_Tri10_IP<nIP>::getValueOnIPs(const double* value) {

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
template<> const double* Elem_Tri10_IP<6>::m_weight = &Hammer2D::Wg_6P[0];
template<> const double* Elem_Tri10_IP<7>::m_weight = &Hammer2D::Wg_7P[0];
template<> const double* Elem_Tri10_IP<12>::m_weight = &Hammer2D::Wg_12P[0];
template<> const double* Elem_Tri10_IP<13>::m_weight = &Hammer2D::Wg_13P[0];

// ================================================================================================
//
// Shape functions
//
// ================================================================================================
template<> const double Elem_Tri10_IP<6>::m_Psi[6][m_NumNodes] = {
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

template<> const double Elem_Tri10_IP<7>::m_Psi[7][m_NumNodes] = {
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

template<> const double Elem_Tri10_IP<12>::m_Psi[12][m_NumNodes] = {
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

template<> const double Elem_Tri10_IP<13>::m_Psi[13][m_NumNodes] = {
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
template<> const double Elem_Tri10_IP<6>::m_DPsi[6][m_NumNodes][m_Dim] = {
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

template<> const double Elem_Tri10_IP<7>::m_DPsi[7][m_NumNodes][m_Dim] = {
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

template<> const double Elem_Tri10_IP<12>::m_DPsi[12][m_NumNodes][m_Dim] = {
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

template<> const double Elem_Tri10_IP<13>::m_DPsi[13][m_NumNodes][m_Dim] = {
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
