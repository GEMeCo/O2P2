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
#pragma once

// Custom libraries
#include "Common.h"
#include "Mesh.h"
#include "Solver.h"

// M2S2 library
#include "M2S2/M2S2.h"

namespace O2P2 {
	namespace Proc {
		/** @ingroup Processor_Module
		  * @class NonLinearSolver
		  *
		  * @brief Nonlinear solvers
		  * @details This class manages the solution scheme for nonlinear problems.
		  */
		class NonLinearSolver
		{
		private:
			// Default constructor
			NonLinearSolver() = delete;

			// Copy constructor.
			NonLinearSolver(const NonLinearSolver& other) = delete;

			// Move constructor.
			NonLinearSolver(NonLinearSolver&& other) = delete;

		public:
			// Default destructor of private / protected pointers.
			virtual ~NonLinearSolver() = default;

			/** @return the number of iterations required to achieve solution in the last NL process.
			  */
			const int& getIteration() { return this->mv_iteration; }

			/** @return the norm in the last NL process.
			  */
			const double& getNorm() { return this->mv_norm; }

			/** Loop iterator for the selected NonLinear Solver.
			  * @return true if solution is feasible, false otherwise.
			  *
			  * @param theFEModel Container of solution components.
			  * @param initialNorm Initial norm to normalize the solution (norm = norm / initialNorm).
			  * @param symmetry Choose wheter the system is symmetric or not.
			  * @param output Choose whether to output results or not.
			  */
			virtual bool runNLS(O2P2::Proc::Mesh* theFEModel, const double initialNorm = 1., const bool symmetry = true, const bool output = false) = 0;

			/** Loop iterator for the selected NonLinear Solver.
			  * @return true if solution is feasible, false otherwise.
			  *
			  * @param theModel Container of solution components.
			  * @param initialNorm Initial norm to normalize the solution (norm = norm / initialNorm).
			  * @param symmetry Choose wheter the system is symmetric or not.
			  * @param output Choose whether to output results or not.
			  */
			virtual bool runNLS(O2P2::Proc::Comp::MeshPoint<2>* theModel, const double initialNorm = 1., const bool symmetry = false, const bool output = false) = 0;

			/** Loop iterator for the selected NonLinear Solver.
			  * @return true if solution is feasible, false otherwise.
			  *
			  * @param theModel Container of solution components.
			  * @param initialNorm Initial norm to normalize the solution (norm = norm / initialNorm).
			  * @param symmetry Choose wheter the system is symmetric or not.
			  * @param output Choose whether to output results or not.
			  */
			virtual bool runNLS(O2P2::Proc::Comp::MeshPoint<3>* theModel, const double initialNorm = 1., const bool symmetry = false, const bool output = false) = 0;

		protected:
			/** @brief Hidden Constructor - Implemented in derived classes.
			  * @param MaxIt Maximum number of iterations.
			  * @param Tol Tolerance.
			  */
			explicit NonLinearSolver(const int& MaxIt, const double& Tol)
				: mv_minIt(0), mv_maxIt(MaxIt), mv_tolerance(Tol) {};

			/** @brief Hidden Constructor - Implemented in derived classes.
			  * @param MinIt Minimum number of iterations.
			  * @param MaxIt Maximum number of iterations.
			  * @param Tol Tolerance.
			  */
			explicit NonLinearSolver(const int& MinIt, const int& MaxIt, const double& Tol)
				: mv_minIt(MinIt), mv_maxIt(MaxIt), mv_tolerance(Tol) {};

		protected:
			/** @brief Minimum number of iterations. */
			int mv_minIt;

			/** @brief Maximum number of iterations. */
			int mv_maxIt;

			/** @brief Iterations counter. */
			int mv_iteration{ 0 };

			/** @brief Tolerance. */
			double mv_tolerance;

			/** @brief Norm. */
			double mv_norm{ 0. };
		};


		/** @ingroup NLSolver
		  * @class NLS_NewtonRaphson
		  *
		  * @brief Newton-Raphson
		  * @details Iterative process to approach root of a nonlinear function.
		  * The concept is as follows:
		  * @code
		  * while (||x(i+1)|| / ||x(0)|| >= EPSILON) {
		  *	 x(i+1) = x(i) - f(x) / f'(x)
		  * }
		  * @endcode
		  *
		  * @note x(i+1) - x(i) -> LHS = Hessian^(-1) . RHS
		  * @note f(x) -> RHS = FExt - FInt
		  * @note f'(x) -> Hessian
		  */
		class NLS_NewtonRaphson : public NonLinearSolver
		{
		private:
			// Default constructor is deleted. Use explicit constructor only.
			NLS_NewtonRaphson() = delete;

		public:
			/** Constructor for Newton-Raphson solver.
			  * @param MinIt Minimum number of iterations.
			  * @param MaxIt Maximum number of iterations.
			  * @param Tol Tolerance.
			  */
			explicit NLS_NewtonRaphson(const int& MinIt, const int& MaxIt, const double& Tol)
				: NonLinearSolver(MinIt, MaxIt, Tol) {}

			/** Constructor for Newton-Raphson solver.
			  * @param MaxIt Maximum number of iterations.
			  * @param Tol Tolerance.
			  */
			explicit NLS_NewtonRaphson(const int& MaxIt, const double& Tol)
				: NonLinearSolver(MaxIt, Tol) {}

			// Default destructor of private / protected pointers.
			~NLS_NewtonRaphson() = default;

			// Loop iterator for the selected NonLinear Solver.
			bool runNLS(O2P2::Proc::Mesh* theFEModel, const double initialNorm = 1., const bool symmetry = true, const bool output = true) override {
				PROFILE_FUNCTION();
				return this->runNonLinearSolver(theFEModel, initialNorm, symmetry, output);
			};

			// Loop iterator for the selected NonLinear Solver.
			bool runNLS(O2P2::Proc::Comp::MeshPoint<2>* theModel, const double initialNorm = 1., const bool symmetry = false, const bool output = false) override {
				return this->runNonLinearSolver(theModel, initialNorm, symmetry, output);
			}

			// Loop iterator for the selected NonLinear Solver.
			bool runNLS(O2P2::Proc::Comp::MeshPoint<3>* theModel, const double initialNorm = 1., const bool symmetry = false, const bool output = false) override {
				return this->runNonLinearSolver(theModel, initialNorm, symmetry, output);
			}

		private:
			// Actual Implementation of non-linear solver.
			template<class T> bool runNonLinearSolver(T* theModel, const double initialNorm, const bool symmetry, const bool output);
		};

	} // End of Proc Namespace
} // End of O2P2 Namespace


