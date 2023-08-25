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
#include "Material.h"

// ================================================================================================
//
// Saint Venant Kirchhoff isotropic elastic material: constructor
// 
// ================================================================================================
O2P2::Prep::Mat_SVK_ISO::Mat_SVK_ISO(const size_t& index) : Material(index)
{
	mv_E11 = 1.;
	mv_nu12 = 0.;
	mv_Rho = 0.;
	mv_Damp = 0.;

	mv_G12 = 0.;
	mv_Bulk = 0.;
	mv_Lambda = 0.;

}


// ================================================================================================
//
// Saint Venant Kirchhoff isotropic elastic material: constructor
// 
// ================================================================================================
O2P2::Prep::Mat_SVK_ISO::Mat_SVK_ISO(const size_t& index, const std::vector<double>& Param) : Material(index)
{
	if (Param.size() < 2 || Param.size() > 4) {
		LOG("\n\nMat_SVK_ISO.constructor: The number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(this->mv_index) + "\n\n\n");
		throw std::invalid_argument("\n\n\nThe number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(this->mv_index) + "\n\n\n");
	}

	mv_E11 = Param.at(0);
	mv_nu12 = Param.at(1);

	mv_Rho = (Param.size() > 2) ? Param.at(2) : 0.;
	mv_Damp = (Param.size() > 3) ? Param.at(3) : 0.;

	mv_G12 = mv_E11 / (2. * (1. + mv_nu12));
	mv_Bulk = mv_E11 / (3. * (1. - 2. * mv_nu12));
	mv_Lambda = mv_E11 * mv_nu12 / ((1. + mv_nu12) * (1. - 2. * mv_nu12));
}


// ================================================================================================
//
// Implementation of Member Function: setParameters
// Saint Venant Kirchhoff isotropic elastic material
// 
// ================================================================================================
void O2P2::Prep::Mat_SVK_ISO::setParameters(const std::vector<double>& Param)
{
	if (Param.size() < 2 || Param.size() > 4) {
		LOG("\n\nMat_SVK_ISO.setParameters: The number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(this->mv_index) + "\n\n\n");
		throw std::invalid_argument("\n\n\nThe number of input material parameters (" + std::to_string(Param.size()) + ") is wrong for material: " + std::to_string(this->mv_index) + "\n\n\n");
	}

	mv_E11 = Param.at(0);
	mv_nu12 = Param.at(1);

	mv_Rho = (Param.size() > 2) ? Param.at(2) : 0.;
	mv_Damp = (Param.size() > 3) ? Param.at(3) : 0.;

	mv_G12 = mv_E11 / (2. * (1. + mv_nu12));
	mv_Bulk = mv_E11 / (3. * (1. - 2. * mv_nu12));
	mv_Lambda = mv_E11 * mv_nu12 / ((1. + mv_nu12) * (1. - 2. * mv_nu12));
}

// ================================================================================================
//
// Implementation of Member Function: getConstitutiveMatrix
// Saint Venant Kirchhoff isotropic elastic material
// 
// ================================================================================================
M2S2::Dyadic4C O2P2::Prep::Mat_SVK_ISO::getConstitutiveMatrix(int nDim, PlaneStateType PS)
{
	// 3D Constitutive matrix
	M2S2::Dyadic4C mi_E(nDim);

	if (nDim == 3) {
		double mi_aux = mv_E11 / ((1. + mv_nu12) * (1. - 2. * mv_nu12));

		mi_E(0, 0) = mi_aux * (1. - mv_nu12);
		mi_E(1, 1) = mi_E(0, 0);
		mi_E(2, 2) = mi_E(0, 0);

		mi_E(0, 1) = mv_nu12 * mi_aux;
		mi_E(0, 2) = mi_E(0, 1);
		mi_E(1, 2) = mi_E(0, 1);

		mi_E(3, 3) = mv_G12 * 2.;
		mi_E(4, 4) = mv_G12 * 2.;
		mi_E(5, 5) = mv_G12 * 2.;
	}
	else if (PS == PlaneStateType::PLANE_STRESS) {
		double mi_aux = mv_E11 / (1. - mv_nu12 * mv_nu12);

		mi_E(0, 0) = mi_aux;
		mi_E(1, 1) = mi_aux;
		mi_E(0, 1) = mv_nu12 * mi_aux;
		mi_E(2, 2) = mv_G12 * 2.;
	}
	else {
		double mi_aux = mv_E11 / ((1. - 2. * mv_nu12) * (1. + mv_nu12));

		mi_E(0, 0) = (1. - mv_nu12) * mi_aux;
		mi_E(1, 1) = (1. - mv_nu12) * mi_aux;
		mi_E(0, 1) = mv_nu12 * mi_aux;
		mi_E(2, 2) = mv_G12 * 2.;
	}

	return mi_E;
}
