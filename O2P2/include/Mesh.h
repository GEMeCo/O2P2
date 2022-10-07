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

// Eigen libraries
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

// Custom Header Files
#include "Common.h"
#include "LoadStep.h"

#include "Domain.h"
#include "MeshElem.h"
#include "MeshNode.h"
#include "PostProcess.h"

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
			void addDOF(int nDOF) { m_TotalDof += nDOF; }

		public:
			// Default destructor of private / protected pointers.
			virtual ~Mesh() = default;

			/** Add a new load step.
			  * Any change in boundary conditions, time step and such, requires a new load step. All values must be included again.
			  * @param NumSteps Number of time steps in the current load step.
			  * @param TimeStep Variation in time in the current load step.
			  * @param NumDBC Number of applied Dirichlet Boundary Conditions.
			  * @param NumNBC Number of applied Neumann Boundary Conditions.
			  */
			void addLoadStep(const int& NumSteps, const double& TimeStep, const size_t& NumDBC, const size_t& NumNBC) {
				m_LoadStep.push_back(std::make_unique<O2P2::Proc::Comp::LoadStep>(NumSteps, TimeStep, NumDBC, NumNBC));
			}

			/** Add a Boundary Condition of Dirichlet type to a Load Step.
			  * @param nLS Number of the load step to receive the BC.
			  * @param Dof Degree of Freedom with imposed boundary condition.
			  * @param Value Value of the boundary condition.
			  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  */
			void addDirichletBC(const int& nLS, const size_t& Dof, const double& Value, const double Var[]) {
				m_LoadStep.at(nLS)->addDirichletBC(Dof, Value, Var);
			}

			/** Add a Boundary Condition of Neumann type to a Load Step.
			  * @param nLS Number of the load step to receive the BC.
			  * @param Dof Degree of Freedom with imposed boundary condition.
			  * @param Value Value of the boundary condition.
			  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
			  */
			void addNeumannBC(const int& nLS, const size_t& Dof, const double& Value, const double Var[]) {
				m_LoadStep.at(nLS)->addNeumannBC(Dof, Value, Var);
			}

			/** Assemble the system of equation, made by the Hessian matrix and a right hand side vector.
			  * @param Hessian Square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
			  * @param RHS Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS.
			  */
			virtual void assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS) = 0;

			/** Initiates the auxiliary vector used to impose boundary conditions to the system of equations.
			  * @param loadStep Current load step under process.
			  */
			void initDirichletBC(const int& loadStep);

			/** Sets the current timestep.
			  * @param timeStep Current timestep under process.
			  */
			void setTimeStep(const int& timeStep) { this->m_curTimeStep = timeStep; }

			/** Sets the current timestep and Newmark-beta parameters.
			  * @param timeStep Current timestep under process.
			  * @param beta First Newmark-beta parameter, associated to displacement.
			  * @param gamm Second Newmark-beta parameter, associated to velocity.
			  */
			virtual void setTimeStep(const int& timeStep, const double& beta, const double& gamma) = 0;

			/** Update trial solution.
			  * @param LHS Left hand side vector with current trial solution.
			  */
			virtual void setTrial(Eigen::VectorXd& LHS) = 0;

			/** Once solution is achieved, commit it. */
			virtual void setCommit() = 0;

			/** @return a pointer to the current load step.
			  */
			O2P2::Proc::Comp::LoadStep* getLoadStep() { return m_LoadStep.at(m_curLoadStep).get(); }

			/** Retrieve a pointer to a load step.
			  * @return A pointer to load step iLS.
			  * @param iLS Number of Load Step to be retrieved.
			  */
			O2P2::Proc::Comp::LoadStep* getLoadStep(int iLS) { return m_LoadStep.at(iLS).get(); }

			/** @return the number of the degrees of freedom (size of the problem).
			  */
			const size_t getNumDof() { return m_TotalDof; }

			/** Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			  * @param RHS Right hand side vector to impose the Neumann boundary conditions.
			  */
			virtual void imposeNeumannBC(Eigen::VectorXd& RHS) = 0;

		public:
			/** @brief Vector with Dirichlet boundary conditions - if there is a BC imposed, then it is equal to 0. Otherwise, it is 1. */
			std::vector<int> m_BCIndex;

			/** @brief Current time of anaysis. */
			double m_currentTime{ 0. };

			/** @brief Initial norm - required by the NLSolver. */
			double m_initialNorm{ 0 };

			/** @brief Container of element components associated to the analysis. */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::MeshElem>> m_meshElem;

			/** @brief Container of node components associated to the analysis. */
			std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> m_meshNode;

		protected:
			/** @brief Current load step under process. */
			int m_curLoadStep{ 0 };

			/** @brief Current time step of current load step. */
			int m_curTimeStep{ 0 };

			/** @brief Total number of DOF */
			size_t m_TotalDof{ 0 };

			/** @brief Container of boundary conditions for load steps. */
			std::vector<std::unique_ptr<O2P2::Proc::Comp::LoadStep>> m_LoadStep;

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
			/** Constructor for container of solution components.
			  * @param theDomain Container with nodal and elements information.
			  * @param vOut Reference to post-process container.
			  */
			explicit Mesh_MQS(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* vOut) : Mesh(vOut) {
				this->m_meshNode.reserve(theDomain->m_nNodes);
				this->m_meshElem.reserve(theDomain->m_nElem);

				this->DomainNodesToMesh(theDomain);
				this->DomainElemToMesh(theDomain);

				// Initial norm (required by Non-linear Solver)
				for (auto& x0 : theDomain->getNode()) {
					for (double coord : x0->getInitPos()) {
						this->m_initialNorm += coord * coord;
					}
				}
				this->m_initialNorm = std::sqrt(this->m_initialNorm);
			}

			// Default destructor of private / protected pointers.
			~Mesh_MQS() = default;

			// Assemble the system of equation, made by the Hessian matrix and the right hand side vector.
			void assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS) override;

			// Update trial solution.
			void setTrial(Eigen::VectorXd& LHS) override;

			// Update commit solution.
			void setCommit() override;

			// Sets the current timestep and Newmark-beta parameters.
			void setTimeStep(const int& timeStep, const double& beta, const double& gamma) override {
				O2P2::Proc::Mesh::setTimeStep(timeStep);
			};

			//Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			void imposeNeumannBC(Eigen::VectorXd& RHS) override;

		protected:
			// Generates a mesh node for each domain node
			void DomainNodesToMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Prep::Node<nDim>>& node : theDomain->getNode()) {
					// Index of DOF for current node
					// For now, each node has nDim DOF
					size_t iDof = node->m_index * nDim;
					this->m_meshNode.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshNode_MQS<nDim>>(iDof, node->getInitPos()));
					addDOF(this->m_meshNode.back()->getNumDOF());
				}
			}

			// Generates a mesh element for each domain element
			void DomainElemToMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Prep::Elem::Element<nDim>>& elem : theDomain->getElem()) {

					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

					for (int i = 0; i < elem->getNumNodes(); i++) {
						conect[i] = this->m_meshNode[elem->getConectivity(i)->m_index];
					}
					this->m_meshElem.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshElem_SVK<nDim>>(elem, conect));
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
			/** Constructor for container of solution components.
			  * @param theDomain Container with nodal and elements information.
			  * @param vOut Reference to post-process container.
			  */
			explicit Mesh_MD(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* vOut) : Mesh_MQS<nDim>(theDomain, vOut) { 
				v_Qs.resize(this->m_TotalDof);
				v_Rs.resize(this->m_TotalDof);

				m_beta = 0;
				m_gamma = 0;
			}

			// Default destructor of private / protected pointers.
			~Mesh_MD() = default;

			// Assemble the system of equation, made by the Hessian matrix and the right hand side vector.
			void assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS) override;

			// Update trial solution.
			void setTrial(Eigen::VectorXd& LHS) override;

			// Update commit solution.
			void setCommit() override;

			// Sets the current timestep and Newmark-beta parameters.
			void setTimeStep(const int& timeStep, const double& beta, const double& gamma) override;

			//Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			void imposeNeumannBC(Eigen::VectorXd& RHS) override;

		protected:
			// Generates a mesh node for each domain node
			void DomainNodesToMesh(O2P2::Prep::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Prep::Node<nDim>>& node : theDomain->getNode()) {
					// Index of DOF for current node
					// For now, each node has nDim DOF
					size_t iDof = node->m_index * nDim;
					this->m_meshNode.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshNode_MD<nDim>>(iDof, node->getInitPos()));
					this->addDOF(this->m_meshNode.back()->getNumDOF());
				}
			}

		protected:
			// Dynamic contribution from previous step
			Eigen::VectorXd v_Qs;

			// Dynamic contribution from previous step
			Eigen::VectorXd v_Rs;

			// Time integration parameters
			double m_beta, m_gamma;
		};

	} // End of Proc Namespace
} // End of O2P2 Namespace
