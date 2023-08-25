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
#include "SolutionAlgorithm.h"


// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template bool O2P2::Proc::SolutionAlgorithm::initFEM<2>(O2P2::Prep::Domain<2>* theDomain, O2P2::Post::PostProcess* thePost);
template bool O2P2::Proc::SolutionAlgorithm::initFEM<3>(O2P2::Prep::Domain<3>* theDomain, O2P2::Post::PostProcess* thePost);

template void O2P2::Proc::SolutionAlgorithm::runSolution<2>(O2P2::Prep::Domain<2>* theDomain);
template void O2P2::Proc::SolutionAlgorithm::runSolution<3>(O2P2::Prep::Domain<3>* theDomain);


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initFEM
//
// ================================================================================================
template<int nDim> inline bool O2P2::Proc::SolutionAlgorithm::initFEM(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	LOG("\nSolutionAlgorithm.initFEModel: Basic Definitions");
	LOG("SolutionAlgorithm.initFEModel: Initiating DOF Mapping system");

	switch (mv_AnalysisType)
	{
	case AnalysisType::STATIC:
	{
		mv_theFEModel = std::make_unique<O2P2::Proc::Mesh_MQS<nDim>>(theDomain, thePost);
		O2P2::Proc::Mesh_MQS<nDim>* theModel = static_cast<O2P2::Proc::Mesh_MQS<nDim>*>(mv_theFEModel.get());
		theModel->initMesh(theDomain);
		break;
	}
	case AnalysisType::TRANSIENT_2ndORDER_NEWMARK:
	{
		mv_theFEModel = std::make_unique<O2P2::Proc::Mesh_MD<nDim>>(theDomain, thePost);
		O2P2::Proc::Mesh_MD<nDim>* theModel = static_cast<O2P2::Proc::Mesh_MD<nDim>*>(mv_theFEModel.get());
		theModel->initMesh(theDomain);
		break;
	}
	case AnalysisType::TRANSIENT_1stORDER:
	case AnalysisType::TRANSIENT_2ndORDER_HHT_Alpha:
	case AnalysisType::EIGENVALUE:
	default: { throw std::invalid_argument("\n\nUndefined type of analysis\nCheck input file\n\n\n"); break; }
	}
	LOG("SolutionAlgorithm.initFEModel: Total number of DOF: " << std::to_string(mv_theFEModel->getNumDof()));

	return true;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): runSolutionAlgorithm
//
// ================================================================================================
template<int nDim> inline void O2P2::Proc::SolutionAlgorithm::runSolution(O2P2::Prep::Domain<nDim>* theDomain)
{
	std::cout << "\n\n------------------------------------------------------------"
		<< "\nIniting Solution Algorithm"
		<< "\n------------------------------------------------------------";

	for (int loadStep = 0; loadStep < mv_numLoadSteps; ++loadStep)
	{
		LOG("\nSolutionAlgorithm.runSolutionAlgorithm: Current LoadStep: " << std::to_string(loadStep + 1));
		std::cout << "\nLoad Step = " << loadStep + 1
			<< "\n------------------------------------------------------------";

		// Setup Dirichlet boundary conditions
		mv_theFEModel->initDirichletBC(loadStep);

		// Second loop -> Time Stepping
		mv_theTimeStep->runTimeLoop(mv_theFEModel.get(), mv_theNLSolver.get());
	}
}
