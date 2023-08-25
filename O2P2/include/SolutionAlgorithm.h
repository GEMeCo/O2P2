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

// Custom Header Files
#include "Common.h"
#include "Domain.h"

#include "Mesh.h"
#include "TimeStepping.h"
#include "NonLinearSolver.h"

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
			// Default constructor will be deleted - use the explicit constructor below
			SolutionAlgorithm() = delete;

			// SolutionAlgorithm is unique, therefore it should not be copied. Thus, copy constructor will be deleted.
			SolutionAlgorithm(const SolutionAlgorithm& other) = delete;

			// Move constructor will also be deleted.
			SolutionAlgorithm(SolutionAlgorithm&& other) = delete;

		public:
			// Default destructor of private / protected pointers.
			~SolutionAlgorithm() = default;

			/** The only constructor that there is. It defines the type of analysis and solvers employed.
			  * @param Analysis Type of analysis - Quasi-static, dynamic or eigenvalue.
			  * @param Solver Type of nonlinear solver - i.e. Newton-Raphson.
			  * @param Problem Type of problem - mechanical, thermodynamics or coupled thermomechanical.
			  * @param numLS Number of Load Steps.
			  * @param minIt Minimum number of iterations.
			  * @param maxIt Maximum number of iterations .
			  * @param tolerance Tolerance.
			  */
			explicit SolutionAlgorithm(const AnalysisType& Analysis, const NLSolverType& Solver, const ProblemType& Problem,
				const int& numLS, const int& minIt, const int& maxIt, const double& tolerance) :
				mv_AnalysisType(Analysis), mv_SolverType(Solver), mv_ProblemType(Problem), mv_numLoadSteps(numLS) {

				std::cout << "\nAnalysis is set to ";
				switch (mv_AnalysisType)
				{
				case AnalysisType::STATIC:
					std::cout << "Quasi-static ";
					mv_theTimeStep = std::make_unique<TimeStep_QsiStatic>();
					break;
				case AnalysisType::TRANSIENT_2ndORDER_NEWMARK:
					mv_theTimeStep = std::make_unique<TimeStep_2ndNew>();
					std::cout << "Newmark dynamic ";
					break;
				case AnalysisType::TRANSIENT_1stORDER:
					std::cout << "1st order dynamic ";
					break;
				case AnalysisType::TRANSIENT_2ndORDER_HHT_Alpha:
					std::cout << "HHT Alpha dynamic  ";
					break;
				case AnalysisType::EIGENVALUE:
					std::cout << "Eigenvalue ";
					break;
				default:
					throw std::invalid_argument("\n\n\nUndefined time integration scheme\nCheck input file\n\n\n");
					break;
				}

				switch (mv_SolverType)
				{
				case NLSolverType::NEWTONRAPHSON:
					std::cout << "with Newton-Raphson solver.";
					mv_theNLSolver = std::make_unique<NLS_NewtonRaphson>(minIt, maxIt, tolerance);
					break;
				case NLSolverType::BFGS:
				case NLSolverType::NEWTONWLS:
				case NLSolverType::NEWTONWTR:
				default:
					throw std::invalid_argument("\n\n\nUndefined Solver\nCheck input file\n\n\n");
					break;
				}

				mv_theFEModel = nullptr;
			}

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
			int getNumLoadSteps() const { return mv_numLoadSteps; }

			/** @return the selected analysis type. */
			AnalysisType getAnalysisType() const { return mv_AnalysisType; }

			/** @return the selected problem type. */
			ProblemType getTypeOfAnalysis() const { return mv_ProblemType; }

			/** @return the selected nonlinear solver. */
			NLSolverType getSolverType() const { return mv_SolverType; }

			/** Add a load step to the Mesh.
			  * @param numSteps Number of time steps in the current load step.
			  * @param timeStep Variation in time in the current load step.
			  * @param numDBC Number of applied Dirichlet Boundary Conditions.
			  * @param numNBC Number of applied Neumann Boundary Conditions.
			  *
			  * @sa Mesh
			  */
			void addLoadStep(const int& numSteps, const double& timeStep, const size_t& numDBC, const size_t& numNBC) {
				LOG("SolutionAlgorithm.addLoadStep: Number of time steps: " << std::to_string(numSteps));
				LOG("SolutionAlgorithm.addLoadStep: Time Step: " << std::scientific << timeStep << std::fixed);
				LOG("SolutionAlgorithm.addLoadStep: Number of Dirichlet Boundary Conditions: " << std::to_string(numDBC));
				LOG("SolutionAlgorithm.addLoadStep: Number of Neumann Boundary Conditions: " << std::to_string(numNBC));
				mv_theFEModel->addLoadStep(numSteps, timeStep, numDBC, numNBC);
			}

			/** Add a Boundary Condition of Dirichlet type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param value Value of the boundary condition.
			  * @param var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  *
			  * @sa Mesh
			  */
			void addDirichletBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				mv_theFEModel->addDirichletBC(nLS, index, iDir, value, var);
			}

			/** Add a Boundary Condition of Neumann type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param value Value of the boundary condition.
			  * @param var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  *
			  * @sa Mesh
			  */
			void addNeumannBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				mv_theFEModel->addNeumannBC(nLS, index, iDir, value, var);
			}

			/** Sets the time stepping parameters.
			  * @param vl1 Alfa for first order transient analysis.
			  * @param vl2 Beta for second order transient analysis.
			  * @param vl3 Gamma for second order transient analysis.
			  */
			void SetTSP(const double& vl1, const double& vl2, const double& vl3) { 
				mv_theTimeStep->SetParameters(vl1, vl2, vl3);
			}

		private:
			template<int nDim>
			bool initFEM(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost);

			template<int nDim>
			void runSolution(O2P2::Prep::Domain<nDim>* theDomain);

		private:
			/** Type of analysis */
			AnalysisType mv_AnalysisType;

			/** Type of nonlinear solver */
			NLSolverType mv_SolverType;

			/** Type of problem */
			ProblemType mv_ProblemType;

			/** Number of load steps */
			int mv_numLoadSteps;

			/** Pointer to the nonlinear solver */
			std::unique_ptr<O2P2::Proc::NonLinearSolver> mv_theNLSolver;

			/** Pointer to the time integration step scheme */
			std::unique_ptr<O2P2::Proc::TimeStepping> mv_theTimeStep;

			/** Container of solution components */
			std::unique_ptr<O2P2::Proc::Mesh> mv_theFEModel;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace

