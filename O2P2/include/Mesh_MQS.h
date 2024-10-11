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
#include "Mesh.h"
#include "Domain.h"

namespace O2P2 {
	namespace Proc {
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
			// Default constructor is deleted. Use explicit constructor only.
			Mesh_MQS() = delete;

			// Mesh is unique. Copy constructor will be deleted.
			Mesh_MQS(const Mesh& other) = delete;

			// Mesh is unique. Move constructor will be deleted.
			Mesh_MQS(Mesh&& other) = delete;

		public:
			// Default destructor of private / protected pointers.
			~Mesh_MQS() = default;

			/** Constructor for container of solution components.
			  * @param theDomain Container with nodal and elements information.
			  * @param thePost Reference to post-process container.
			  */
			explicit Mesh_MQS(O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost) : Mesh(thePost) {
				mv_meshNode.reserve(theDomain->getNumNodes());
				mv_meshPoint.reserve(theDomain->getNumPoints());

				mv_meshElem.reserve(theDomain->getNumElems());
				mv_meshBars.reserve(theDomain->getNumBars());
				mv_meshFace.reserve(theDomain->getNumFace());
				mv_meshIncl.reserve(theDomain->getNumIncl());
				mv_meshFibs.reserve(theDomain->getNumFibs());

				// Initial norm (required by Nonlinear Solver)
				this->mv_initialNorm = 0.;
				for (std::shared_ptr<O2P2::Geom::Node>& node : theDomain->getNode()) {
					for (int i = 0; i < nDim; ++i) {
						this->mv_initialNorm += node->getInitPos()[i] * node->getInitPos()[i];
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
			void initMesh(O2P2::Geom::Domain<nDim>* theDomain) {
				domainNodesToMesh(theDomain);
				domainElemToMesh(theDomain);
			}

		protected:
			// Generates a mesh node from each domain node. Do the same for points (but do not increase the number of DOF)
			void domainNodesToMesh(O2P2::Geom::Domain<nDim>* theDomain);

			// Generates a mesh element for each domain element
			void domainElemToMesh(O2P2::Geom::Domain<nDim>* theDomain);

			// Stablish point references (which mesh element the point is immersed on and what xsi values the point has).
			bool evaluatePoints() override;

			// Impose Inertial forces to the vector of independent terms.
			void imposeInertiaBC(std::vector<double>& RHS) override {};

			// Generates the system of equation for matrix elements
			// dt is only used in dynamic processes
			// elem must be a pointer because the container is made of unique_ptr
			void ElemSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override;

			// Generates the system of equation for inclusion elements
			// dt is only used in dynamic processes
			// elem must be a pointer because the container is made of unique_ptr
			void InclSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace