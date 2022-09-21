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
#include "NonLinearSolver.h"

// ================================================================================================
//
// Explicit template member function instantiation
//
// ================================================================================================
template bool O2P2::Proc::NLS_NewtonRaphson::runNonLinearSolver<O2P2::Proc::Mesh>(O2P2::Proc::Mesh* theModel, const double initialNorm, const bool output);

// ================================================================================================
//
// Implementation of Template Member Function: runNonLinearSolver
//
// ================================================================================================
template<class T> bool O2P2::Proc::NLS_NewtonRaphson::runNonLinearSolver(T* theModel, const double initialNorm, const bool output)
{
	// Returns false if results are infeasible.
	bool feasible = false;

	// Hessian is a square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
	Eigen::SparseMatrix<double> Hessian(theModel->getNumDof(), theModel->getNumDof());

	// Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS
	Eigen::VectorXd RHS(theModel->getNumDof());

	// Vector of dependent terms, the left hand side, for system: Hessian.LHS = RHS
	Eigen::VectorXd LHS(theModel->getNumDof());

	// Load vector (Neumann Boundary conditions)
	Eigen::VectorXd FExt(theModel->getNumDof());

	// Sparse solver
	Eigen::PardisoLU<Eigen::SparseMatrix<double>> solver;

	// Assembly load vector (from input)
	theModel->imposeNeumannBC(FExt);

	// Nonlinear iteration (third nested loop)
	this->m_iteration = 0;

	// The first iteration is required - the difference here is that the Hessian pattern is stablished.
	{
		this->m_iteration++;

		// The external forces do not change inside this loop
		RHS = FExt;

		// Asks the template model to assemble the system of equations (RHS = FExt - FInt)
		theModel->assembleSOE(Hessian, RHS);

		{
			PROFILE_SCOPE("solver");
			// Solve the system of equation
			solver.analyzePattern(Hessian);
			solver.factorize(Hessian);
			if (solver.info() != Eigen::Success) {
				LOG("Decomposition failed in iteration " + std::to_string(m_iteration));
				throw std::invalid_argument("\n\n\nDecomposition failed!\n\n");
			}
			LHS = solver.solve(RHS);
		}

		// Update trial solution
		theModel->setTrial(LHS);

		// Evaluates the dot product of the solution (quadratic of variation).
		m_Norm = LHS.norm();
		m_Norm /= initialNorm;

		if (output) {
			std::cout << "\nIteration: " << this->m_iteration << " Norm: " << this->m_Norm;
			LOG("NLS_NewtonRaphson.runNonLinearSolver: Iteration: " + std::to_string(this->m_iteration) << "; Norm: " << std::scientific << this->m_Norm << std::fixed);
		}
	}

	while (((this->m_Norm > this->m_Tolerance) && (this->m_iteration < this->m_MaxIt)) || (this->m_iteration < this->m_MinIt))
	{
		this->m_iteration++;

		// The external forces do not change inside this loop
		RHS = FExt;

		// Asks the template model to assemble the system of equations (RHS = FExt - FInt)
		theModel->assembleSOE(Hessian, RHS);

		{
			PROFILE_SCOPE("solver");

			// Solve the system of equation
			//solver.analyzePattern(Hessian);
			solver.factorize(Hessian);
			if (solver.info() != Eigen::Success) {
				LOG("Decomposition failed in iteration " + std::to_string(m_iteration));
				throw std::invalid_argument("\n\n\nDecomposition failed!\n\n");
			}
			LHS = solver.solve(RHS);
		}

		// Update trial solution
		theModel->setTrial(LHS);

		// Evaluates the dot product of the solution (quadratic of variation).
		m_Norm = LHS.norm();
		m_Norm /= initialNorm;

		// Evaluate current norm
		if (output) {
			std::cout << "\nIteration: " << this->m_iteration << " Norm: " << this->m_Norm;
			LOG("NLS_NewtonRaphson.runNonLinearSolver: Iteration: " + std::to_string(this->m_iteration) << "; Norm: " << std::scientific << this->m_Norm << std::fixed);
		}
	};

	// If convergence was attained, the solution is feasible.
	if (this->m_Norm < this->m_Tolerance) feasible = true;

	// Update commit solution
	theModel->setCommit();

	return feasible;
}
