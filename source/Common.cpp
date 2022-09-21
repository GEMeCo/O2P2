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
#include "Common.h"

// Only create log file in debug mode
#ifdef _DEBUG
std::ofstream logFile;
#endif //_DEBUG

std::map<AnalysisType, std::string> analysisTypeNames =
{ {AnalysisType::STATIC,						"Quasi-static analysis"},
  {AnalysisType::TRANSIENT_1stORDER,			"First order transient analysis"},
  {AnalysisType::TRANSIENT_2ndORDER_NEWMARK,	"Newmark time step integration method for second order transient analysis"},
  {AnalysisType::TRANSIENT_2ndORDER_HHT_Alpha,	"HHT alpha time step integration method for second order transient analysis"},
  {AnalysisType::EIGENVALUE,					"Eingevalue and Eigenvector analysis"}
};

std::map<ProblemType, std::string> problemTypeNames =
{ {ProblemType::MECHANICAL, "Mechanical"},
  {ProblemType::THERMAL,	"Heat transfer"},
  {ProblemType::COUPLED,	"Coupled thermomechanical"}
};


std::map<NLSolverType, std::string> NLSolverTypeNames =
{ {NLSolverType::NEWTONRAPHSON,	"Newton-Raphson Method"},
  {NLSolverType::BFGS,			"Broyden-Fletcher-Goldfarb-Shanno Method"},
  {NLSolverType::NEWTONWLS,		"Newton based method with Line Search"},
  {NLSolverType::NEWTONWTR,		"Newton based method with Trust Region"}
};

std::map<MaterialType, std::string> MaterialTypeNames =
{ {MaterialType::SVK_ISO,		"Elastic isotropic, based on Saint-Venant-Kirchhof constitutive model"},
  {MaterialType::SVK_ORT,		"Elastic orthotropic, based on Saint-Venant-Kirchhof constitutive model"},
  {MaterialType::SVK_DAMAGE,	"Isotropic damage model, based on Saint-Venant-Kirchhof constitutive model"},
  {MaterialType::CONCRETE,		"Isotropic damage model, based on SVK constitutive model, including creep, shinkage, AAR and DEF"},
  {MaterialType::STEEL,			"SVK isotropic linear elastoplastic plus relaxation, presstress and losses"},
  {MaterialType::RSD_ISO,		"Elastic isotropic, based on Rivlin-Saunders-Duster constitutive model"},
  {MaterialType::RSD_PLASTIC,	"Elasto-Plastic isotropic, based on Rivlin-Saunders-Duster constitutive model"},
  {MaterialType::TEMP_ISO_LIN,	"Isotropic heat transfer, with linear properties"},
  {MaterialType::TEMP_ISO_NLIN,	"Isotropic heat transfer, with nonlinear properties"},
  {MaterialType::RSD_TP,		"For Thermomechanical problems - RSD_PLASTIC + TEMP_ISO_NLIN"}
};

std::map<PlaneStateType, std::string> PlaneStateTypeNames =
{ {PlaneStateType::PLANE_STRESS, "Stress Plane State"},
  {PlaneStateType::PLANE_STRAIN, "Strain Plane State"}
};

std::map<OutputType, std::string> OutputTypeExtension =
{ {OutputType::OGL, "_out.ogl"},
  {OutputType::VTU, "_out.vtu"}
};


// Format function definition
std::ostream& formatFixed(std::ostream& os) {
	os.fill(' ');
	return os << std::fixed << std::setprecision(6) << std::setw(10);
}

std::ostream& formatScien(std::ostream& os) {
	os.fill(' ');
	return os << std::scientific << std::setprecision(6) << std::setw(14);;
}
