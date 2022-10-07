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
#pragma once

// Custom Header Files
#include "Common.h"
#include "Domain.h"

#include "Mesh.h"
#include "TimeStepping.h"
#include "NonLinearSolver.h"
#include "OutputSystem.h"

namespace O2P2 {
	namespace Proc {
		/** @ingroup Processor_Module
		  * @class SolutionAlgorithm
		  *
		  * @brief Aggregation of classes to evaluate the analysis.
		  * @details This class manages the entire processing, aggregating a nonlinear solver, a time step integration scheme and a solver for system of equations.
		  */
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
			explicit SolutionAlgorithm(const AnalysisType& Analysis, const NLSolverType& Solver, const ProblemType& Problem,
				const int& NumLS, const int& MinIt, const int& MaxIt, const double& Tolerance) :
				m_AnalysisType(Analysis), m_SolverType(Solver), m_ProblemType(Problem), m_numLoadSteps(NumLS) {

				// Creates the object related to time stepping scheme
				std::cout << "\nAnalysis is set to ";
				switch (m_AnalysisType) {
				case AnalysisType::STATIC:
					std::cout << "Quasi-Static ";
					m_theTimeStep = std::make_unique<TimeStep_QsiStatic>();
					break;
				case AnalysisType::TRANSIENT_2ndORDER_NEWMARK:
					std::cout << "Dynamic ";
					m_theTimeStep = std::make_unique<TimeStep_2ndNew>();
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
			}

			// Default destructor of private / protected pointers.
			~SolutionAlgorithm() = default;

			/** Initiate the solution container.
			  * @return true if container is populated.
			  *
			  * @param theDomain Container with nodal and elements information.
			  * @param thePost Container with solutions for post-process.
			  */
			bool initFEModel(O2P2::Prep::Domain<2>* theDomain, O2P2::Post::PostProcess* thePost) {
				PROFILE_FUNCTION();
				return this->initFEM(theDomain, thePost);
			}

			/** Initiate the solution container.
			  * @return true if container is populated.
			  *
			  * @param theDomain Container with nodal and elements information.
			  * @param thePost Container with solutions for post-process.
			  */
			bool initFEModel(O2P2::Prep::Domain<3>* theDomain, O2P2::Post::PostProcess* thePost) {
				PROFILE_FUNCTION();
				return this->initFEM(theDomain, thePost);
			}

			/** Run the analysis, based on previous setup.
			  * @param theDomain Container with nodal and elements information.
			  */
			void runSolutionAlgorithm(O2P2::Prep::Domain<2>* theDomain) {
				PROFILE_FUNCTION();
				this->runSolution(theDomain);
			}

			/** Run the analysis, based on previous setup.
			  * @param theDomain Container with nodal and elements information.
			  */
			void runSolutionAlgorithm(O2P2::Prep::Domain<3>* theDomain) {
				PROFILE_FUNCTION();
				this->runSolution(theDomain);
			}

			/** @return the number of Load Steps. */
			int getNumLoadSteps() const { return m_numLoadSteps; }

			/** @return the selected analysis type. */
			AnalysisType getAnalysisType() const { return m_AnalysisType; }

			/** @return the selected problem type. */
			ProblemType getTypeOfAnalysis() const { return m_ProblemType; }

			/** @return the selected nonlinear solver. */
			NLSolverType getSolverType() const { return m_SolverType; }

			/** Add a load step to the Mesh.
			  * @param NumSteps Number of time steps in the current load step.
			  * @param TimeStep Variation in time in the current load step.
			  * @param NumDBC Number of applied Dirichlet Boundary Conditions.
			  * @param NumNBC Number of applied Neumann Boundary Conditions.
			  *
			  * @sa Mesh
			  */
			void addLoadStep(const int& NumSteps, const double& TimeStep, const size_t& NumDBC, const size_t& NumNBC) {
				LOG("SolutionAlgorithm.addLoadStep: Number of time steps: " << std::to_string(NumSteps));
				LOG("SolutionAlgorithm.addLoadStep: Time Step: " << std::scientific << TimeStep << std::fixed);
				LOG("SolutionAlgorithm.addLoadStep: Number of Dirichlet Boundary Conditions: " << std::to_string(NumDBC));
				LOG("SolutionAlgorithm.addLoadStep: Number of Neumann Boundary Conditions: " << std::to_string(NumNBC));
				m_theFEModel->addLoadStep(NumSteps, TimeStep, NumDBC, NumNBC);
			}

			/** Add a Boundary Condition of Dirichlet type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param Value Value of the boundary condition.
			  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  *
			  * @sa Mesh
			  */
			void addDirichletBC(const int& nLS, const size_t& index, const int& iDir, const double& Value, const double Var[]) {
				// Check node container for current DOF
				size_t Dof = m_theFEModel->m_meshNode[index]->m_DofIndex + iDir;

				m_theFEModel->addDirichletBC(nLS, Dof, Value, Var);
			}

			/** Add a Boundary Condition of Neumann type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param Value Value of the boundary condition.
			  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  *
			  * @sa Mesh
			  */
			void addNeumannBC(const int& nLS, const size_t& index, const int& iDir, const double& Value, const double Var[]) {
				// Check node container for current DOF
				size_t Dof = m_theFEModel->m_meshNode[index]->m_DofIndex + iDir;

				m_theFEModel->addNeumannBC(nLS, Dof, Value, Var);
			}

			/** Sets the time stepping parameters.
			  * @param vl1 Alfa for first order transient analysis.
			  * @param vl2 Beta for second order transient analysis.
			  * @param vl3 Gamma for second order transient analysis.
			  */
			void SetTSP(const double& vl1, const double& vl2, const double& vl3) { 
				m_theTimeStep->SetParameters(vl1, vl2, vl3);
			}

		private:
			template<int nDim>
			void runSolution(O2P2::Prep::Domain<nDim>* theDomain);

			template<int nDim>
			bool initFEM(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost);

		private:
			/** Type of analysis */
			AnalysisType m_AnalysisType;

			/** Type of nonlinear solver */
			NLSolverType m_SolverType;

			/** Type of problem */
			ProblemType m_ProblemType;

			/** Number of load steps */
			int m_numLoadSteps;

			/** Pointer to the nonlinear solver */
			std::unique_ptr<O2P2::Proc::NonLinearSolver> m_theNLSolver;

			/** Pointer to the time integration step scheme */
			std::unique_ptr<O2P2::Proc::TimeStepping> m_theTimeStep;

			/** Container of solution components */
			std::unique_ptr<O2P2::Proc::Mesh> m_theFEModel;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace
