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
// Implementation of Template Member Function (2D and 3D): initFEModel
//
// ================================================================================================
template<int nDim>
bool O2P2::Proc::SolutionAlgorithm::initFEM(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	LOG("\nSolutionAlgorithm.initFEModel: Basic Definitions");
	LOG("SolutionAlgorithm.initFEModel: Initiating DOF Mapping system");

	// Creates the container of Solution Components
	m_theFEModel = std::make_unique<O2P2::Proc::Mesh_Mec<nDim>>(theDomain, thePost);

	LOG("SolutionAlgorithm.initFEModel: Total number of DOF: " << std::to_string(m_theFEModel->getNumDof()));

	return true;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): runSolutionAlgorithm
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::SolutionAlgorithm::runSolution(O2P2::Prep::Domain<nDim>* theDomain)
{
	std::cout << "\n\n------------------------------------------------------------"
		<< "\nIniting Solution Algorithm"
		<< "\n------------------------------------------------------------";

	// First loop -> Load Steps
	for (int loadStep = 0; loadStep < m_numLoadSteps; ++loadStep)
	{
		LOG("\nSolutionAlgorithm.runSolutionAlgorithm: Current LoadStep: " << std::to_string(loadStep + 1));
		std::cout << "\nLoad Step = " << loadStep + 1
			<< "\n------------------------------------------------------------";

		// Setup Dirichlet boundary conditions
		m_theFEModel->initDirichletBC(loadStep);

		// Second loop -> Time Stepping
		m_theTimeStep->runTimeLoop(m_theFEModel.get(), m_theNLSolver.get());
	}
}