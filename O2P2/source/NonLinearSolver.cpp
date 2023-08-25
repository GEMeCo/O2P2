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
#include "NonLinearSolver.h"

// ================================================================================================
//
// Explicit template member function instantiation
//
// ================================================================================================
template bool O2P2::Proc::NLS_NewtonRaphson::runNonLinearSolver(O2P2::Proc::Mesh* theModel, const double initialNorm, const bool output);


// ================================================================================================
//
// Implementation of Template Member Function: runNonLinearSolver
//
// ================================================================================================
template<class T>
inline bool O2P2::Proc::NLS_NewtonRaphson::runNonLinearSolver(T* theModel, const double initialNorm, const bool output)
{
	// Returns false if results are infeasible.
	bool feasible = false;

	// Number of DOF
	int size = theModel->getNumDof();

	// Hessian is a symmetric square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
	M2S2::sparseMatrix mi_Hessian(size, true);

	// Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS
	std::vector<double> mi_RHS(size);

	// Vector of solution, the left hand side, for system: Hessian.LHS = RHS
	std::vector<double> mi_LHS(size);

	// Load vector (Neumann Boundary conditions)
	std::vector<double> mi_FExt(size);

	// Assembly load vector (from input)
	theModel->imposeNeumannBC(mi_FExt);

	// Nonlinear iteration (third nested loop)
	this->mv_iteration = 0;

	while (((this->mv_norm > this->mv_tolerance) && (this->mv_iteration < this->mv_maxIt)) || (this->mv_iteration < this->mv_minIt))
	{
		this->mv_iteration++;

		// The external forces do not change inside this loop
		mi_RHS = mi_FExt;

		// Asks the template model to assemble the system of equations (RHS = FExt - FInt)
		theModel->assembleSOE(mi_Hessian, mi_RHS);

		// Pardiso uses CSR matrices, thus needs conversion
		M2S2::CSR mi_csrMatrix;
		mi_Hessian.saveAsCSR(mi_csrMatrix);

		{
			PROFILE_SCOPE("solver");
			solve(size, mi_csrMatrix, mi_RHS, mi_LHS);
		}

		// Update trial solution
		theModel->setTrial(mi_LHS);

		// Evaluates the dot product of the solution (quadratic of variation).
		this->mv_norm = std::sqrt(std::inner_product(mi_LHS.begin(), mi_LHS.end(), mi_LHS.begin(), 0.)) / initialNorm;

		// Evaluate current norm
		if (output) {
			std::cout << "\nIteration: " << this->mv_iteration << " Norm: " << this->mv_norm;
			LOG("NLS_NewtonRaphson.runNonLinearSolver: Iteration: " + std::to_string(this->mv_iteration) << "; Norm: " << std::scientific << this->mv_norm << std::fixed);
		}
	}

	// If convergence was attained, the solution is feasible.
	if (this->mv_norm < this->mv_tolerance) feasible = true;

	// Update commit solution
	theModel->setCommit();

	return feasible;
}