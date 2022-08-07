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

// Custom Header Files
#include "Common.h"
#include "Domain.h"

#include "AnalysisComp.h"
#include "TimeStepping.h"
#include "NonLinearSolver.h"
#include "OutputSystem.h"

/** @ingroup Processor_Module
  * @class SolutionAlgorithm
  *
  * @brief Aggregation of classes to evaluate the analysis.
  * @details This class manages the entire processing, aggregating a nonlinear solver, a time step integration scheme and a solver for system of equations.
  * 
  * @todo Adicionar uma figura mostrando o fluxograma do algoritmo de solução (calls). SolAlg -> TimeStep -> NLSolver -> AnalysisComp -> ElemComp. E a LoadStep?
  *
  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
  */
template<int nDim>
class SolutionAlgorithm
{
private:
	SolutionAlgorithm() = delete;

	// SolutionAlgorithm is unique, therefore it should not be copied. Thus, copy constructor will be deleted.
	SolutionAlgorithm(const SolutionAlgorithm& other) = delete;

	// Move constructor will also be deleted.
	SolutionAlgorithm(SolutionAlgorithm&& other) = delete;

public:
	/** The only constructor that there is. It defines the type of analysis and solvers employed.
	  * @param Analysis Type of analysis - Quasi-static, dynamic or eigenvalue.
	  * @param Solver Type of nonlinear solver - i.e. Newton-Raphson.
	  * @param Problem Type of problem - mechanical, thermodynamics or coupled thermomechanical.
	  * @param NumLS Number of Load Steps.
	  * @param MinIt Minimum number of iterations.
	  * @param MaxIt Maximum number of iterations .
	  * @param Tolerance Tolerance.
	  */
	SolutionAlgorithm(const AnalysisType& Analysis, const NLSolverType& Solver, const ProblemType& Problem,
		              const int& NumLS, const int& MinIt, const int& MaxIt, const double& Tolerance) :
		m_AnalysisType(Analysis), m_SolverType(Solver), m_ProblemType(Problem), m_numLoadSteps(NumLS) {

		// Creates the object related to time stepping scheme
		std::cout << "\nAnalysis is set to ";
		switch (m_AnalysisType) {
		case AnalysisType::STATIC:
			std::cout << "Quasi-Static ";
			m_theTimeStep = std::make_unique<TimeStep_QsiStatic>();
			break;
		default:
			// if none, thows an error message
			throw std::invalid_argument("\n\n\nUndefined time integration scheme\nCheck input file\n\n\n");
		}

		switch (m_SolverType) {
		case NLSolverType::NEWTONRAPHSON:
			std::cout << "with Newton-Raphson solver.";
			m_theNLSolver = std::make_unique<NLS_NewtonRaphson>(MinIt, MaxIt, Tolerance);
			break;
		default:
			// if none, thows an error message
			throw std::invalid_argument("\n\n\nUndefined Solver\nCheck input file\n\n\n");
		}

		m_theFEModel = nullptr;
		m_SolutionAcquired = false;
	};

	// Default destructor.
	~SolutionAlgorithm() = default;

	/** Initiate the solution container.
	  * @return true if container is populated.
	  * 
	  * @param theDomain Container with nodal and elements information.
	  * @param thePost Container with solutions for post-process.
	  */
	bool initFEModel(Domain<nDim>* theDomain, PostProcess* thePost);

	/** Run the analysis, based on previous setup.
	  * @param theDomain Container with nodal and elements information.
	  */
	void runSolutionAlgorithm(Domain<nDim>* theDomain);

	/** @return the number of Load Steps. */
	int getNumLoadSteps() const { return m_numLoadSteps; };

	/** @return the selected analysis type. */
	AnalysisType getAnalysisType() const { return m_AnalysisType; };

	/** @return the selected problem type. */
	ProblemType getTypeOfAnalysis() const { return m_ProblemType; };

	/** @return the selected nonlinear solver. */
	NLSolverType getSolverType() const { return m_SolverType; };

	/** @return true if analysis succeeded. */
	bool isFinished() { return m_SolutionAcquired; };

	/** Add a load step to the AnalysisComp.
	  * @param NumSteps Number of time steps in the current load step.
	  * @param TimeStep Variation in time in the current load step.
	  * @param NumDBC Number of applied Dirichlet Boundary Conditions.
	  * @param NumNBC Number of applied Neumann Boundary Conditions.
	  * 
	  * @sa AnalysisComp
	  */
	void addLoadStep(const int& NumSteps, const double& TimeStep, const size_t& NumDBC, const size_t& NumNBC) {
		LOG("SolutionAlgorithm.addLoadStep: Number of time steps: " << std::to_string(NumSteps));
		LOG("SolutionAlgorithm.addLoadStep: Time Step: " << std::scientific << TimeStep << std::fixed);
		LOG("SolutionAlgorithm.addLoadStep: Number of Dirichlet Boundary Conditions: " << std::to_string(NumDBC));
		LOG("SolutionAlgorithm.addLoadStep: Number of Neumann Boundary Conditions: " << std::to_string(NumNBC));
		m_theFEModel->addLoadStep(NumSteps, TimeStep, NumDBC, NumNBC);
	};

	/** Add a Boundary Condition of Dirichlet type to a Load Step.
	  * @param nLS Number of the load step to recieve the BC.
	  * @param Dof Degree of freedom with imposed boundary condition.
	  * @param Value Value of the boundary condition.
	  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
	  * 
	  * @sa AnalysisComp
	  */
	void addDirichletBC(const int& nLS, const size_t& Dof, const double& Value, const double Var[]) {
		m_theFEModel->addDirichletBC(nLS, Dof, Value, Var);
	};

	/** Add a Boundary Condition of Neumann type to a Load Step.
	  * @param nLS Number of the load step to recieve the BC.
	  * @param Dof Degree of freedom with imposed boundary condition.
	  * @param Value Value of the boundary condition.
	  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
	  * 
	  * @sa AnalysisComp
	  */
	void addNeumannBC(const int& nLS, const size_t& Dof, const double& Value, const double Var[]) {
		m_theFEModel->addNeumannBC(nLS, Dof, Value, Var);
	};

private:
	/** Type of analysis */
	AnalysisType m_AnalysisType;

	/** Type of nonlinear solver */
	NLSolverType m_SolverType;

	/** Type of problem */
	ProblemType m_ProblemType;

	/** Number of load steps */
	int m_numLoadSteps;

	/** Boolean for the end of analysis */
	bool m_SolutionAcquired;

	/** Pointer to the nonlinear solver */
	std::unique_ptr<NonLinearSolver> m_theNLSolver;

	/** Pointer to the time integration step scheme */
	std::unique_ptr<TimeStepping> m_theTimeStep;

	/** Container of solution components */
	std::unique_ptr<AnalysisComp> m_theFEModel;
};
