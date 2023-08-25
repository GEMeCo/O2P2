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
#include "TimeStepping.h"


// ================================================================================================
//
// Implementation of Member Function: runTimeLoop
//
// ================================================================================================
void O2P2::Proc::TimeStep_QsiStatic::runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver)
{
	PROFILE_FUNCTION();

	// Loop on time step
	for (int mi_timeIt = 0; mi_timeIt < theMesh->getLoadStep()->mv_numSteps; ++mi_timeIt) {
		// Update current time
		theMesh->mv_currentTime += theMesh->getLoadStep()->mv_timeStep;

		LOG("TimeStep_QsiStatic.runTimeLoop: Current Time Step: " + std::to_string(mi_timeIt) + "; Current Analysis Time: " + std::to_string(theMesh->mv_currentTime));
		std::cout << "\nCurrent analysis time = " << theMesh->mv_currentTime
			<< "\n------------------------------------------------------------";

		// Setup current analysis time step
		theMesh->setTimeStep(mi_timeIt);

		// Call loop of non-linear solver
		bool mi_feasible = theSolver->runNLS(theMesh, theMesh->mv_initialNorm, true);
		LOG("NLSolver.runNLS: Solution attained - feasibility: " << std::boolalpha << mi_feasible << std::noboolalpha);
		std::cout << "\n------------------------------------------------------------";
	}
}


// ================================================================================================
//
// Implementation of Member Function: runTimeLoop
//
// ================================================================================================
void O2P2::Proc::TimeStep_2ndNew::runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver)
{
	PROFILE_FUNCTION();

	// At the beggining of the first load step (at which currentTime is zero), the initial acceleration must be evaluated
	if ((int)(theMesh->mv_currentTime * 1000000) == 0) theMesh->setAccel();

	// Loop on time step
	for (int mi_timeIt = 0; mi_timeIt < theMesh->getLoadStep()->mv_numSteps; ++mi_timeIt) {
		// Update current time
		theMesh->mv_currentTime += theMesh->getLoadStep()->mv_timeStep;

		LOG("TimeStep_2ndNew.runTimeLoop: Current Time Step: " + std::to_string(mi_timeIt) + "; Current Analysis Time: " + std::to_string(theMesh->mv_currentTime));
		std::cout << "\nCurrent analysis time = " << theMesh->mv_currentTime
			<< "\n------------------------------------------------------------";

		// Setup current analysis time step
		// Also evaluate the dynamic contribution from previous step (v_Qs and v_Rs)
		theMesh->setTimeStep(mi_timeIt, mv_beta, mv_gamma);

		// Call loop of non-linear solver
		bool mi_feasible = theSolver->runNLS(theMesh, theMesh->mv_initialNorm, true);
		LOG("NLSolver.runNLS: Solution attained - feasibility: " << std::boolalpha << mi_feasible << std::noboolalpha);
		std::cout << "\n------------------------------------------------------------";
	}

}
