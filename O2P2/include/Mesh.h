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
#include "MeshElem.h"
#include "MeshNode.h"
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
			Mesh() = delete;

		protected:
			/** Constructor for container of solution components.
			  * @param vOut Reference to post-process container.
			  */
			Mesh(O2P2::Post::PostProcess* vOut) { m_PostPt = vOut; }

			/** Add new dof to the system. Used to calculate the total number of DOF.
			  * @param nDOF Number of DOF to be added to the system of equation.
			  */
			void addDof(int nDof) { mv_TotalDof += nDof; }

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
			  * @param var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  */
			void addDirichletBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				// Check node container for current DOF
				size_t dof = mv_meshNode[index]->mv_DofIndex + iDir;
				mv_LoadStep.at(nLS)->addDirichletBC(dof, value, var);
			}

			/** Add a Boundary Condition of Neumann type to a Load Step.
			  * @param nLS Number of the load step to recieve the BC.
			  * @param index Node container index with imposed boundary condition.
			  * @param iDir  Direction of the imposed boundary condition.
			  * @param value Value of the boundary condition.
			  * @param var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  */
			void addNeumannBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				// Check node container for current DOF
				size_t dof = mv_meshNode[index]->mv_DofIndex + iDir;
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
			void initDirichletBC(const int& loadStep);

			/** Sets the current timestep.
			  * @param timeStep Current timestep under process.
			  */
			void setTimeStep(const int& timeStep) { this->mv_curTimeStep = timeStep; }

			/** Sets the current timestep and Newmark-beta parameters.
			  * @param timeStep Current timestep under process.
			  * @param beta First Newmark-beta parameter, associated to displacement.
			  * @param gamm Second Newmark-beta parameter, associated to velocity.
			  */
			virtual void setTimeStep(const int& timeStep, const double& beta, const  double& gamma) = 0;

			/** Update trial solution.
			  * @param LHS Left hand side vector with current trial solution.
			  */
			virtual void setTrial(std::vector<double>& LHS) = 0;

			/** Once solution is achieved, commit it. */
			virtual void setCommit() = 0;

			/** Set initial acceleration. */
			virtual void setAccel() {};

			/** @return a pointer to the current load step.
			  */
			O2P2::Proc::Comp::LoadStep* getLoadStep() { return mv_LoadStep.at(mv_curLoadStep).get(); }

			/** Retrieve a pointer to a load step.
			  * @return A pointer to load step iLS.
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

		public:
			/** @brief Vector with Dirichlet boundary conditions - if there is a BC imposed, then it is equal to 0. Otherwise, it is 1. */
			std::vector<size_t> mv_BCindex;

			/** @brief Current time of anaysis. */
			double mv_currentTime{ 0. };

			/** @brief Initial norm - required by the NLSolver. */
			double mv_initialNorm{ 0 };

			/** @brief Container of element components associated to the analysis. */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElem>> mv_meshElem;

			/** @brief Container of node components associated to the analysis. */
			std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> mv_meshNode;

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
		};

		/**
		  * @class Mesh_MQS
		  *
		  * @brief Container of solution components for mechanical quasi-static problems.
		  * @details This class holds vectors of solution and analysis components: Load Steps, DOF, System of equation, and such.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  */
		template<int nDim>
		class Mesh_MQS : public Mesh
		{
		private:
			Mesh_MQS() = delete;

		public:
			// Default destructor of private / protected pointers.
			~Mesh_MQS() = default;

			/** Constructor for container of solution components.
			  * @param theDomain Container with nodal and elements information.
			  * @param vOut Reference to post-process container.
			  */
			explicit Mesh_MQS(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* vOut) : Mesh(vOut) {
				mv_meshNode.reserve(theDomain->mv_nNodes);
				mv_meshElem.reserve(theDomain->mv_nElem);

				// Initial norm (required by Non-linear Solver)
				for (std::shared_ptr<O2P2::Prep::Node<nDim>>& node : theDomain->getNode()) {
					for (double coord : node->getInitPos()) {
						this->mv_initialNorm += coord * coord;
					}
				}
				this->mv_initialNorm = std::sqrt(this->mv_initialNorm);
			}

			// Assemble the system of equation, made by the Hessian matrix and the right hand side vector.
			void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override;

			// Update trial solution.
			void setTrial(std::vector<double>& LHS) override;

			// Update commit solution.
			void setCommit() override;

			// Sets the current timestep and Newmark-beta parameters.
			void setTimeStep(const int& timeStep, const double& beta, const  double& gamma) override {
				O2P2::Proc::Mesh::setTimeStep(timeStep);
			}

			// Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			void imposeNeumannBC(std::vector<double>& RHS) override;

			// Initiates mesh elements and nodes.
			void initMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				this->DomainNodesToMesh(theDomain);
				this->DomainElemToMesh(theDomain);
			}

		protected:
			// Generates a mesh node from each domain node
			void DomainNodesToMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Prep::Node<nDim>>& node : theDomain->getNode()) {
					size_t dof = node->mv_index * nDim;
					this->mv_meshNode.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshNode_MQS<nDim>>(dof, node->getInitPos()));
					addDof(this->mv_meshNode.back()->getNumDOF());
				}
			}

			// Generates a mesh element for each domain element
			void DomainElemToMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Prep::Elem::Element<nDim>>& elem : theDomain->getElem()) {
					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

					for (int i = 0; i < elem->getNumNodes(); ++i) {
						conect[i] = this->mv_meshNode[elem->getConectivity(i)->mv_index];
					}
					this->mv_meshElem.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshElem_SVK<nDim>>(elem, conect));
				}
			}
		};


		/**
		  * @class Mesh_MD
		  *
		  * @brief Container of solution components for mechanical dynamic problems.
		  * @details This class holds vectors of solution and analysis components: Load Steps, DOF, System of equation, and such.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  */
		template<int nDim>
		class Mesh_MD : public Mesh_MQS<nDim>
		{
		private:
			Mesh_MD() = delete;

		public:
			// Default destructor of private / protected pointers.
			~Mesh_MD() = default;

			/** Constructor for container of solution components.
			  * @param theDomain Container with nodal and elements information.
			  * @param vOut Reference to post-process container.
			  */
			explicit Mesh_MD(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* vOut) : Mesh_MQS<nDim>(theDomain, vOut) {
				mv_beta = 0;
				mv_gamma = 0;
			}

			// Add a Boundary Condition of Dirichlet type to a Load Step.
			// Also imposes a initial velocity to the node.
			void addDirichletBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				if (nLS == 0) {
					// Get a downcasted pointer to the MD node
					O2P2::Proc::Comp::MeshNode_MD<nDim>* node = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(this->mv_meshNode[index].get());

					// The velocity is the derivative of the displacement
					node->setVel(iDir, value * var[1]);
				}
				O2P2::Proc::Mesh::addDirichletBC(nLS, index, iDir, value, var);
			}

			// Assemble the system of equation, made by the Hessian matrix and the right hand side vector.
			void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override;

			// Update trial solution.
			void setTrial(std::vector<double>& LHS) override;

			// Update commit solution.
			void setCommit() override;

			// Sets the current timestep and Newmark-beta parameters.
			void setTimeStep(const int& timeStep, const double& beta, const  double& gamma) override;

			// Set initial acceleration.
			void setAccel() override;

			// Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			void imposeNeumannBC(std::vector<double>& RHS) override {
				Mesh_MQS<nDim>::imposeNeumannBC(RHS);
			};

			// Initiates mesh elements and nodes.
			void initMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				this->DomainNodesToMesh(theDomain);
				Mesh_MQS<nDim>::DomainElemToMesh(theDomain);

				mv_Qs.resize(this->mv_TotalDof);
				mv_Rs.resize(this->mv_TotalDof);
			}

		protected:
			// Generates a mesh node for each domain node
			void DomainNodesToMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Prep::Node<nDim>>& node : theDomain->getNode()) {
					size_t dof = node->mv_index * nDim;
					this->mv_meshNode.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshNode_MD<nDim>>(dof, node->getInitPos()));
					this->addDof(this->mv_meshNode.back()->getNumDOF());
				}
			}

			// Impose Inertial forces to the vector of independent terms.
			void imposeInertiaBC(std::vector<double>& RHS);

		protected:
			// Dynamic contribution from previous step
			std::vector<double> mv_Qs;

			// Dynamic contribution from previous step
			std::vector<double> mv_Rs;

			// Time integration parameters
			double mv_beta, mv_gamma;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace