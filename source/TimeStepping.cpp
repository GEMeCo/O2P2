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
// TimeStepping
// 
// Time step integration schemes
// 
// ================================================================================================
#include "TimeStepping.h"

// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template void O2P2::Proc::TimeStep_QsiStatic::runTimeLooping<2>(O2P2::Prep::Domain<2>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver);
template void O2P2::Proc::TimeStep_QsiStatic::runTimeLooping<3>(O2P2::Prep::Domain<3>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver);

// ================================================================================================
//
// Implementation of Template Member Function: runTimeLoop
//
// ================================================================================================
template<int nDim> void O2P2::Proc::TimeStep_QsiStatic::runTimeLooping(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) {

    // Loop on time step
    for (int timeIt = 0; timeIt < theFEModel->getLoadStep()->m_NumSteps; ++timeIt) {

        // Update current time
        theFEModel->m_currentTime += theFEModel->getLoadStep()->m_TimeStep;

        LOG("TimeStep_QsiStatic.runTimeLoop: Current Time Step: " + std::to_string(timeIt) + "; Current Analysis Time: " + std::to_string(theFEModel->m_currentTime));

        std::cout << "\nCurrent analysis time = " << theFEModel->m_currentTime
            << "\n------------------------------------------------------------";

        // Setup current analysis time step
        theFEModel->setTimeStep(timeIt);

        // Initial norm
        double initialNorm = 0.;
        for (auto& x0 : theDomain->getNode()) {
            for (double coord : x0->getInitPos()) {
                initialNorm += coord * coord;
            }
        }
        initialNorm = std::sqrt(initialNorm);

        // Call loop of non-linear solver
        bool feasible = theSolver->runNLS(theFEModel, initialNorm, true);
        LOG("NLSolver.runNLS: Solution attained - feasibility: " << std::boolalpha << feasible << std::noboolalpha);
        std::cout << "\n------------------------------------------------------------";
    }
}
