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
#include "Mesh_MQS.h"

namespace O2P2 {
	namespace Proc {
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
			// Default constructor is deleted. Use explicit constructor only.
			Mesh_MD() = delete;

			// Mesh is unique. Copy constructor will be deleted.
			Mesh_MD(const Mesh& other) = delete;

			// Mesh is unique. Move constructor will be deleted.
			Mesh_MD(Mesh&& other) = delete;

		public:
			// Default destructor of private / protected pointers.
			~Mesh_MD() = default;

			/** Constructor for container of solution components.
			  * @param theDomain Container with nodal and elements information.
			  * @param thePost Reference to post-process container.
			  */
			explicit Mesh_MD(O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost) : Mesh_MQS<nDim>(theDomain, thePost) {
				mv_beta = 0;
				mv_gamma = 0;
			}

			// Add a Boundary Condition of Dirichlet type to a Load Step.
			// Also imposes a initial velocity to the node.
			void addDirichletBC(const int& nLS, const size_t& index, const int& iDir, const double& value, const double var[]) {
				if (nLS == 0) {
					// Get a downcasted pointer to the MD node
					assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(this->mv_meshNode.at(index)));
					auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(this->mv_meshNode.at(index));

					// The velocity is the derivative of the displacement (at t = 0)
					pNode->setVel(iDir, value * (var[1] + var[3] * var[4]));
				}
				O2P2::Proc::Mesh::addDirichletBC(nLS, index, iDir, value, var);
			}

			// Assemble the system of equation, made by the Hessian matrix and the right hand side vector.
			void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override {
				Mesh_MQS<nDim>::assembleSOE(Hessian, RHS);
			};

			// Update trial solution.
			void setTrial(std::vector<double>& LHS) override;

			// Update commit solution.
			void setCommit() override;

			// Sets the current timestep and Newmark-beta parameters.
			void setTimeStep(const int& timeStep, const double& beta, const  double& gamma) override;

			// Set initial acceleration.
			void setInitAccel() override;

			// Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
			void imposeNeumannBC(std::vector<double>& RHS) override {
				Mesh_MQS<nDim>::imposeNeumannBC(RHS);
			};

			// Initiates mesh elements and nodes.
			void initMesh(O2P2::Geom::Domain<nDim>* theDomain) {
				this->domainNodesToMesh(theDomain);
				Mesh_MQS<nDim>::domainElemToMesh(theDomain);

				mv_Qs.resize(this->mv_TotalDof);
				mv_QsIm.resize(nDim * this->mv_meshPoint.size());
			}

		protected:
			// Generates a mesh node for each domain node
			void domainNodesToMesh(O2P2::Geom::Domain<nDim>* theDomain) {
				for (std::shared_ptr<O2P2::Geom::Node>& node : theDomain->getNode()) {
					size_t dof = node->mv_index * nDim;
					this->mv_meshNode.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshNode_MD<nDim>>(dof, node->getInitPos()));
					this->addDof(this->mv_meshNode.back()->getNumDOF());
				}

				for (std::shared_ptr<O2P2::Geom::Node>& node : theDomain->getPoint()) {
					this->mv_meshPoint.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(0, node->getInitPos()));
				}
			}

			// Impose Inertial forces to the vector of independent terms.
			void imposeInertiaBC(std::vector<double>& RHS) override;

			// Generates the system of equation for matrix elements
			// elem must be a pointer because the container is made of unique_ptr
			void ElemSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override;

			// Generates the system of equation for inclusion elements
			// elem must be a pointer because the container is made of unique_ptr
			void InclSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override;

		protected:
			// Generates element inertial forces
			void ElemIFO(O2P2::Proc::Comp::MeshElemVI* elem, std::vector<double>& RHS);

			// Generates inclusion inertial forces
			void InclIFO(O2P2::Proc::Comp::MeshElemVI* elem, std::vector<double>& RHS);

			// Generates initial acceleration for elements
			void ElemIAC(O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Mass, std::vector<double>& RHS);

			// Generates initial acceleration for inclusion elements
			void InclIAC(O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Mass, std::vector<double>& RHS);

		protected:
			// Dynamic contribution from previous step for matrix elements
			std::vector<double> mv_Qs;

			// Dynamic contribution from previous step for immersed / inclusion elements
			std::vector<double> mv_QsIm;

			// Time integration parameters
			double mv_beta, mv_gamma;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace