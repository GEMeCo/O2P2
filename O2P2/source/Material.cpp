// ================================================================================================
// 
// This file is part of O2P2, an object oriented environment for the positional FEM
//
// Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved.
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
#include "Material.h"

// ================================================================================================
//
// Saint Venant Kirchhoff isotropic elastic material: constructor
// 
// ================================================================================================
O2P2::Prep::Mat_SVK_ISO::Mat_SVK_ISO(const size_t& index) : Material(index) {
	m_E11 = 1.;
	m_nu12 = 0.;
	m_Rho = 0.;
	m_Damp = 0.;

	m_G12 = 0.;
	m_Bulk = 0.;
	m_Lambda = 0.;
};

// ================================================================================================
//
// Saint Venant Kirchhoff isotropic elastic material: constructor
// 
// ================================================================================================
O2P2::Prep::Mat_SVK_ISO::Mat_SVK_ISO(const size_t& index, const std::vector<double>& Param) : Material(index) {

	if (Param.size() < 2 || Param.size() > 4) {
		LOG("\n\nMat_SVK_ISO.constructor: The number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(index) + "\n\n\n");
		throw std::invalid_argument("\n\n\nThe number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(index) + "\n\n\n");
	}
	m_E11 = Param.at(0);
	m_nu12 = Param.at(1);

	m_Rho = (Param.size() > 2) ? Param.at(2) : 0.;
	m_Damp = (Param.size() > 3) ? Param.at(3) : 0;

	m_G12 = m_E11 / (2. * (1. + m_nu12));
	m_Bulk = m_E11 / (3. * (1. - 2. * m_nu12));
	m_Lambda = m_E11 * m_nu12 / ((1. + m_nu12) * (1. - 2. * m_nu12));
};

// ================================================================================================
//
// Implementation of Member Function: setParameters
// Saint Venant Kirchhoff isotropic elastic material
// 
// ================================================================================================
void O2P2::Prep::Mat_SVK_ISO::setParameters(const std::vector<double>& Param) {

	if (Param.size() < 2 || Param.size() > 4) {
		LOG("\n\nMat_SVK_ISO.setParameters: The number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(this->m_index) + "\n\n\n");
		throw std::invalid_argument("\n\n\nThe number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(this->m_index) + "\n\n\n");
	}
	m_E11 = Param.at(0);
	m_nu12 = Param.at(1);

	m_Rho = (Param.size() > 2) ? Param.at(2) : 0.;
	m_Damp = (Param.size() > 3) ? Param.at(3) : 0;

	m_G12 = m_E11 / (2. * (1. + m_nu12));
	m_Bulk = m_E11 / (3. * (1. - 2. * m_nu12));
	m_Lambda = m_E11 * m_nu12 / ((1. + m_nu12) * (1. - 2. * m_nu12));
};
