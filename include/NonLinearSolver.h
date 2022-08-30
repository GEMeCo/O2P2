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
#pragma once

// C++ standard libraries
#include <memory>				// required by std::shared_ptr
#include <math.h>               // required by std::sqrt
#include <iostream>             // standard input and output stream

// Eigen libraries
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

// Custom libraries
#include "Profiler.h"
#include "AnalysisComp.h"

/** @ingroup Processor_Module
  * @class NonLinearSolver
  *
  * @brief Nonlinear solvers
  * @details This class manages the solution scheme for nonlinear problems.
  *
  * @todo 1 - Funcionalidade para algorítmo BFGS.
  * @todo 2 - Funcionalidade para Newton with Line Search.
  * @todo 3 - Funcionalidade para Newton with Trust Region.
  */
class NonLinearSolver
{
private:
    NonLinearSolver() = delete;

public:
    /** @brief Destructor */
    virtual ~NonLinearSolver() = default;

	/** @return the number of iterations required to achieve solution in the last NL process.
      */
	const int& getIteration() { return this->m_iteration; };

    /** @return the norm in the last NL process.
      */
    const double& getNorm() { return this->m_Norm; };

    /** Loop iterator for the selected NonLinear Solver.
      * @return true if solution is feasible, false otherwise.
      * 
      * @param theFEModel Container of solution components.
      * @param initialNorm Initial norm to normalize the solution (norm = norm / initialNorm).
      * @param output Choose whether to output results or not.
      */
    virtual bool runNLS(AnalysisComp* theFEModel, const double initialNorm = 1., const bool output = false) = 0;

protected:
    /** @brief Minimum number of iterations. */
    int m_MinIt;

    /** @brief Maximum number of iterations. */
    int m_MaxIt;

    /** @brief Iterations counter. */
    int m_iteration{ 0 };

    /** @brief Tolerance. */
    double m_Tolerance;

    /** @brief Norm. */
    double m_Norm{ 0. };

protected:
    /**
      * @brief Hidden Constructor - Implemented in derived classes.
      * @param MaxIt Maximum number of iterations.
      * @param Tol Tolerance.
      */
    NonLinearSolver(const int& MaxIt, const double& Tol) : m_MinIt(0), m_MaxIt(MaxIt), m_Tolerance(Tol) {};

    /**
      * @brief Hidden Constructor - Implemented in derived classes.
      * @param MinIt Minimum number of iterations.
      * @param MaxIt Maximum number of iterations.
      * @param Tol Tolerance.
      */
    NonLinearSolver(const int& MinIt, const int& MaxIt, const double& Tol) : m_MinIt(MinIt), m_MaxIt(MaxIt), m_Tolerance(Tol) {};
};

/** @ingroup NLSolver
  * @class NLS_NewtonRaphson
  *
  * @brief Newton-Raphson
  * @details Iterative process to approach root of a nonlinear function.
  * The concept is as follows:
  * @code
  * while (||x(i+1)|| / ||x(0)|| >= EPSILON) {
  *     x(i+1) = x(i) - f(x) / f'(x) 
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
    NLS_NewtonRaphson() = delete;

public:
    /** Constructor for Newton-Raphson solver.
      * @param MinIt Minimum number of iterations.
      * @param MaxIt Maximum number of iterations.
      * @param Tol Tolerance.
      */
    NLS_NewtonRaphson(const int& MinIt, const int& MaxIt, const double& Tol) : NonLinearSolver(MinIt, MaxIt, Tol) {};

    /** Constructor for Newton-Raphson solver.
      * @param MaxIt Maximum number of iterations.
      * @param Tol Tolerance.
      */
    NLS_NewtonRaphson(const int& MaxIt, const double& Tol) : NonLinearSolver(MaxIt, Tol) {};

    // Destructor
    ~NLS_NewtonRaphson() = default;

    // Loop iterator for the selected NonLinear Solver.
    bool runNLS(AnalysisComp* theFEModel, const double initialNorm = 1., const bool output = false) override {
        PROFILE_FUNCTION();
        return this->runNonLinearSolver(theFEModel, initialNorm, output);
    };

private:
    // Actual Implementation of non-linear solver.
    template<class T> bool runNonLinearSolver(T* theModel, const double initialNorm, const bool output);
};
