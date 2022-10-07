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
#include "TimeStepping.h"

// ================================================================================================
//
// Implementation of Member Function: runTimeLoop
//
// ================================================================================================
void O2P2::Proc::TimeStep_QsiStatic::runTimeLoop(O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) {
	PROFILE_FUNCTION();

	// Loop on time step
	for (int timeIt = 0; timeIt < theFEModel->getLoadStep()->m_NumSteps; ++timeIt) {

		// Update current time
		theFEModel->m_currentTime += theFEModel->getLoadStep()->m_TimeStep;

		LOG("TimeStep_QsiStatic.runTimeLoop: Current Time Step: " + std::to_string(timeIt) + "; Current Analysis Time: " + std::to_string(theFEModel->m_currentTime));

		std::cout << "\nCurrent analysis time = " << theFEModel->m_currentTime
			<< "\n------------------------------------------------------------";

		// Setup current analysis time step
		theFEModel->setTimeStep(timeIt);

		// Call loop of non-linear solver
		bool feasible = theSolver->runNLS(theFEModel, theFEModel->m_initialNorm, true);
		LOG("NLSolver.runNLS: Solution attained - feasibility: " << std::boolalpha << feasible << std::noboolalpha);
		std::cout << "\n------------------------------------------------------------";
	}
}


// ================================================================================================
//
// Implementation of Member Function: runTimeLoop
//
// ================================================================================================
void O2P2::Proc::TimeStep_2ndNew::runTimeLoop(O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) {
	PROFILE_FUNCTION();

	// At the beggining of the first load step, the initial acceleration must be evaluated


	// Loop on time step
	for (int timeIt = 0; timeIt < theFEModel->getLoadStep()->m_NumSteps; ++timeIt) {

		// Update current time
		theFEModel->m_currentTime += theFEModel->getLoadStep()->m_TimeStep;

		LOG("TimeStep_2ndNew.runTimeLoop: Current Time Step: " + std::to_string(timeIt) + "; Current Analysis Time: " + std::to_string(theFEModel->m_currentTime));

		std::cout << "\nCurrent analysis time = " << theFEModel->m_currentTime
			<< "\n------------------------------------------------------------";

		// Setup current analysis time step
		// Also evaluate the dynamic contribution from previous step (v_Qs and v_Rs)
		theFEModel->setTimeStep(timeIt, m_beta, m_gamma);

		// There is a inertia contribution in the external force, which uses Qs and Rs within runNLS
		

		// Call loop of non-linear solver
		bool feasible = theSolver->runNLS(theFEModel, theFEModel->m_initialNorm, true);
		LOG("NLSolver.runNLS: Solution attained - feasibility: " << std::boolalpha << feasible << std::noboolalpha);
		std::cout << "\n------------------------------------------------------------";
	}
}
