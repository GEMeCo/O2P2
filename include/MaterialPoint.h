// ================================================================================================
//
// TODO:
// 1 - O ponto material está ligado ao tipo de material. Na classe base, apenas materiais elásticos.
// 2 - Classe derivada para material com danificação isotrópica.
// 3 - Classe derivada para material plástico.
// 4 - Distinguir 2D de 3D? Isso afeta a existência de matriz de rotação.
// 5 - Incluir temperatura como parâmetro? Ou extrai valor do elemento?
//
// ================================================================================================

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
// Material Point
// 
// Handler of stored information in the integration point, due to material properties.
// It stores previous damage, plastic strain, and any other information associated to the 
// integration point.
// 
// It also stores information that are often used, such the inverse of jacobian matrix.
// 
// ================================================================================================
#pragma once

// Eigen libraries
#include <Eigen/Dense>

/**
 * @brief Class for material information in the integration point.
 * @details Stores informations about damage, plastic strain, and such. There is one material point for each integration point.
*/
class MaterialPoint
{
public:
	MaterialPoint() = delete;

	/**
	* @brief Constructor for material point information recorded in the integration point.
	* @param Jacobian Reference jacobian matrix - A0 / F0.
	*/
	MaterialPoint(const Eigen::MatrixXd& Jacobian) {
		m_RefJacobian = Jacobian.inverse();
		m_Jacobian = Jacobian.determinant();
	};

	/**
	* @brief Constructor for material point information recorded in the integration point.
	* @param Jacobian Reference jacobian matrix - A0 / F0.
	* @param rot Rotation matrix.
	*/
	MaterialPoint(const Eigen::MatrixXd& Jacobian, const Eigen::MatrixXd& rot) {
		m_RefJacobian = Jacobian.inverse();
		m_Jacobian = Jacobian.determinant();
		m_Rot = rot;
	};

	/**
	* @brief Basic destructor
	*/
	~MaterialPoint() = default;

	/**
	 * @brief Returns the inverse of reference Jacobian matrix (A0i / F0i).
	 * @return Inverse of reference Jacobian matrix.
	*/
	Eigen::MatrixXd& getJacobianMatrix() { return m_RefJacobian; };

	/**
	 * @brief Returns the rotation matrix for 2D elements in 3D environments.
	 * @return Rotation matrix.
	*/
	Eigen::MatrixXd& getRotationMatrix() { return m_Rot; };

	/**
	 * @brief Returns the reference Jacobian.
	 * @return Reference Jacobian.
	*/
	double getJacobian() { return m_Jacobian; };

protected:
	/** @brief Reference Jacobian */
	double m_Jacobian{ 0. };

	/** @brief INVERSE of reference Jacobian Matrix / A0 or F0 - Mapping gradient from dimensionless coordinate system to initial position */
	Eigen::MatrixXd m_RefJacobian;

	/** @brief Rotation matrix for 2D inclusions in a 3D environment. Not used otherwise. */
	Eigen::MatrixXd m_Rot;
};


//std::vector<double> vMaterialPoint;		// relacionado ao material - plasticidade / dano / AAR / creep / etc
