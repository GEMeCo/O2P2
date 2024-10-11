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

// custom library
#include "Common.h"
#include "Domain.h"

#include "Mesh Components\MeshNode.h"
#include "Mesh Components\MeshPoint.h"
#include "Mesh Components\MeshElem.h"
#include "Mesh Components\MeshBars.h"
#include "Mesh Components\MeshIncl.h"
#include "Mesh Components\MeshFibs.h"
#include "Mesh Components\MeshFace.h"

#include "LoadStep.h"
#include "PostProcess.h"

#include "Solver.h"

// M2S2 library
#include "M2S2/M2S2.h"

namespace O2P2 {
	namespace Proc {
		/** @ingroup Processor_Module
		  * @class Mesh
		  *
		  * @brief Container of solution components.
		  * @details This class holds vectors of solution and analysis components: Load Steps, DOF, System of equation, and such.
		  */
		class Mesh
		{
		private:
			// Default constructor is deleted. Use explicit constructor only.
			Mesh() = delete;

			// Mesh is unique. Copy constructor will be deleted.
			Mesh(const Mesh& other) = delete;

			// Mesh is unique. Move constructor will be deleted.
			Mesh(Mesh&& other) = delete;

		protected:
			/** Constructor for container of solution components.
			  * @param thePost Reference to post-process container.
			  */
			explicit Mesh(O2P2::Post::PostProcess* thePost) {
				m_PostPt = thePost;
			}

			/** Add new dof to the system. Used to calculate the total number of DOF.
			  * @param nDof Number of DOF to be added to the system of equation.
			  */
			void addDof(int nDof) { mv_TotalDof += nDof; }

			/** Stablish point references (which mesh element and what xsi values the point has).
			  */
			virtual bool evaluatePoints() = 0;

		public:
			// Default destructor of private / protected pointers.
			virtual ~Mesh() = default;

			/** Add a new load step.
			  * Any change in boundary conditions, time step and such, requires a new load step. All values must be included again.
			  * @param numSteps Number of time steps in the current load step.
			  * @param timeStep Variation in time in the current load step.
			  * @param numDBC Number of applied Dirichlet Boundary Conditions.
			  * @param numNBC Number of applied Neumann Boundary Conditions.
			  */
			void addLoadStep(const int& numSteps, const double& timeStep, const size_t& numDBC, const size_t& numNBC) {
				mv_LoadStep.push_back(std::make_unique<O2P2::Proc::Comp::LoadStep>(numSteps, timeStep, numDBC, numNBC));
			}

			/** Add a Boundary Condition of Dirichlet type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param value Value of the boundary condition.
			  * @param var Time behavior (var[0] + var[1].t + var[2].t^2 + var[3].sin(var[4] * dt) + var[5].cos(var[6] * dt))
			  */
			void addDirichletBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				// Check node container for current DOF
				size_t dof = mv_meshNode.at(index)->mv_dofList.at(0) + iDir;
				mv_LoadStep.at(nLS)->addDirichletBC(dof, value, var);
			}

			/** Add a Boundary Condition of Neumann type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param value Value of the boundary condition.
			  * @param var Time behavior (var[0] + var[1].t + var[2].t^2 + var[3].sin(var[4] * dt) + var[5].cos(var[6] * dt))
			  */
			void addNeumannBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				// Check node container for current DOF
				size_t dof = mv_meshNode.at(index)->mv_dofList.at(0) + iDir;
				mv_LoadStep.at(nLS)->addNeumannBC(dof, value, var);
			}

			/** Assemble the system of equation, made by the Hessian matrix and a right hand side vector.
			  * @param Hessian Square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
			  * @param RHS Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS.
			  */
			virtual void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) = 0;

			/** Initiates the auxiliary vector used to impose boundary conditions to the system of equations.
			  * @param loadStep Current load step under process.
			  */
			void initDirichletBC(const int& loadStep) {
				LOG("Mesh.initDirichletBC: Initiating Dirichlet BC system for load step " << std::to_string(loadStep));

				// Register the current load step for latter.
				mv_curLoadStep = loadStep;

				// Clear all boundary condition that was once here, and assign 1 to all.
				mv_BCindex.clear();
				mv_BCindex.assign(mv_TotalDof, 1);

				// If a Dirichlet boundary condition is imposed, m_BCIndex is changed to zero
				for (auto& DBC : mv_LoadStep.at(loadStep)->mv_DirichletBC) {
					mv_BCindex.at(DBC->mv_Dof) = 0;
				}
			}

			/** Sets the current timestep.
			  * @param timeStep Current timestep under process.
			  */
			void setTimeStep(const int& timeStep) { this->mv_curTimeStep = timeStep; }

			/** Sets the current timestep and Newmark-beta parameters.
			  * @param timeStep Current timestep under process.
			  * @param beta First Newmark-beta parameter, associated to displacement.
			  * @param gamma Second Newmark-beta parameter, associated to velocity.
			  */
			virtual void setTimeStep(const int& timeStep, const double& beta, const  double& gamma) = 0;

			/** Update trial solution.
			  * @param LHS Left hand side vector with current trial solution.
			  */
			virtual void setTrial(std::vector<double>& LHS) = 0;

			/** Once solution is achieved, commit it. */
			virtual void setCommit() = 0;

			/** Set initial acceleration. */
			virtual void setInitAccel() {};

			/** @return a pointer to the current load step.
			  */
			O2P2::Proc::Comp::LoadStep* getLoadStep() { return mv_LoadStep.at(mv_curLoadStep).get(); }

			/** @return a pointer to load step iLS.
			  * @param iLS Number of Load Step to be retrieved.
			  */
			O2P2::Proc::Comp::LoadStep* getLoadStep(int iLS) { return mv_LoadStep.at(iLS).get(); }

			/** @return the number of the degrees of freedom (size of the problem).
			  */
			const size_t& getNumDof() { return mv_TotalDof; }

			/** Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			  * @param RHS Right hand side vector to impose the Neumann boundary conditions.
			  */
			virtual void imposeNeumannBC(std::vector<double>& RHS) = 0;

		protected:
			/** Impose Inertial forces to the vector of independent terms (external forces).
			  * @param RHS Right hand side vector to impose the Neumann boundary conditions.
			  */
			virtual void imposeInertiaBC(std::vector<double>& RHS) = 0;

			/** Generates the system of equation for matrix elements.
			  * @param dt Time at current step. Only used in dynamic processes.
			  * @param elem Pointer to the mesh element. It is not a smart pointer since the container is made of unique_ptr.
  			  * @param Hessian Square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
			  * @param RHS Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS.
			  */
			virtual void ElemSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) = 0;

			/** Generates the system of equation for inclusion elements.
			  * @param dt Time at current step. Only used in dynamic processes.
			  * @param elem Pointer to the mesh element. It is not a smart pointer since the container is made of unique_ptr.
			  * @param Hessian Square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
			  * @param RHS Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS.
			  */
			virtual void InclSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) = 0;

		public:
			/** @brief Vector with Dirichlet boundary conditions - if there is a BC imposed, then it is equal to 0. Otherwise, it is 1. */
			std::vector<size_t> mv_BCindex;

			/** @brief Current time of anaysis. */
			double mv_currentTime{ 0. };

			/** @brief Initial norm - required by the NLSolver. */
			double mv_initialNorm{ 0 };

			/** @brief Container of element components associated to the analysis (Matrix). */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElemVI>> mv_meshElem;

			/** @brief Container of linear components associated to the analysis (Matrix). */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElemVI>> mv_meshBars;

			/** @brief Container of active face element components associated to the analysis (Matrix). */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElemVI>> mv_meshFace;

			/** @brief Container of inclusion components associated to the analysis (Immersed Particles). */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElemVI>> mv_meshIncl;

			/** @brief Container of fiber components associated to the analysis (Immersed Fibers). */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElemVI>> mv_meshFibs;

			/** @brief Container of node components associated to the analysis (Matrix). */
			std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> mv_meshNode;

			/** @brief Container of point components associated to the analysis (Immersed Elements). */
			std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> mv_meshPoint;

		protected:
			/** @brief Current load step under process. */
			int mv_curLoadStep{ 0 };

			/** @brief Current time step of current load step. */
			int mv_curTimeStep{ 0 };

			/** @brief Total number of DOF */
			size_t mv_TotalDof{ 0 };

			/** @brief Container of boundary conditions for load steps. */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::LoadStep>> mv_LoadStep;

			/** @brief Container of solution for post-process. */
			O2P2::Post::PostProcess* m_PostPt;

			/** @brief Mutex to avoid concurrency while assembling system of equation. */
			std::mutex mv_HesMutex;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace
