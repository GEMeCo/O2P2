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
#include "MeshNode.h"
#include "MeshElem.h"

namespace O2P2 {
	namespace Proc {
		// Forward declaration of NLS class
		class NLS_NewtonRaphson;

		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshPoint
			  *
			  * @brief Mesh point, a Solution component. Nodes for inclusion elements.
			  * @details Adds to mesh node class the element in which the point is immersed and the dimensionless coordinates related to that element.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshPoint
			{
			protected:
				/** Constructor for point objects.
				  * @param index Node index.
				  * @param initPos Initial position of the geometry point.
				  * @param Tp Initial temperature.
				  */
				explicit MeshPoint() {
					mv_xsi.fill(0.);
				}

				// Default destructor of private / protected pointers.
				~MeshPoint() = default;

			public:
				/** Function to retrieve nodal dimensionless coordinates.
				  * @return Standard vector with nodal information.
				  */
				const std::array<double, nDim>& getXsi() { return this->mv_xsi; }

				/** Function to retrieve a pointer to the element in which the point is immersed.
				  * @return Pointer to the element in which the point is immersed on.
				  */
				O2P2::Geom::Elem::Element* getElement() { return this->mv_ptMeshElem->mv_ptElem.get(); }

				/** Function to retrieve a pointer to the mesh element in which the point is immersed.
				  * @return Pointer to the element in which the point is immersed on.
				  */
				O2P2::Proc::Comp::MeshElemVI* getMeshElement() { return this->mv_ptMeshElem; }

				/** Sets the dimensionless coordinates.
				  * @param xsi Vector with first guess for dimensionless coordinates.
				  */
				void setXsi(std::vector<double>& xsi) {
					for (int i = 0; i < nDim; ++i) {
						this->mv_xsi[i] = xsi.at(i);
					}
				}

				/** Sets the mesh element that the Point is immersed on.
				  * @param elem Pointer to the mesh element that the Point is immersed on.
				  */
				void setMeshElement(O2P2::Proc::Comp::MeshElemVI* elem) {
					this->mv_ptMeshElem = elem;
				}


			protected:
				// This section was created to stablish the dimensionless coordinates, using a non-linear solver.
				// Grant access to private member functions - friendship is not inherited :(
				friend class O2P2::Proc::NLS_NewtonRaphson;

				/** Required by NonLinearSolver to set matrices sizes.
				  * @return Size for matrices.
				  */
				int getNumDof() { return nDim; }

				/** Impose Neumann Boundary Conditions to the vector of independent terms (dimensionless nodal coordinates).
				  * @param RHS Right hand side vector to impose the Neumann boundary conditions - just an empty vector.
				  */
				void imposeNeumannBC(std::vector<double>& RHS) { 
					memset(&RHS[0], 0., RHS.size() * sizeof(double));
				}

				/** Assemble the system of equation, made by the Hessian matrix and a right hand side vector.
				  * @param Hessian Square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
				  * @param RHS Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS.
				  */
				virtual void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) = 0;

				/** Update trial solution.
				  * @param LHS Left hand side vector with current trial solution.
				  */
				void setTrial(std::vector<double>& LHS) {
					for (int i = 0; i < nDim; ++i) {
						this->mv_xsi.at(i) += LHS.at(i);
					}
				}

				/** Once solution is achieved, commit it. */
				void setCommit() { /* There is nothing to do here. */ }

			protected:
				/** @brief Dimensionless coordinates related to the domain element that contains it. */
				std::array<double, nDim> mv_xsi;

				/** @brief Pointer to the mesh element in which the node is immersed. */
				O2P2::Proc::Comp::MeshElemVI* mv_ptMeshElem{ nullptr };
			};


			/** @ingroup Processor_Module
			  * @class MeshPoint_MQS
			  *
			  * @brief Mesh point, a Solution component, for mechanical quasi-static problems.
			  * @details The solution point component for mechanical quasi-static problems includes current position (trial and commit).
			  * Also adds the element in which the point is immersed and dimensionless coordinates related to that element.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshPoint_MQS : public MeshPoint<nDim>, public MeshNode_MQS<nDim>
			{
			private:
				MeshPoint_MQS() = delete;

			public:
				/** Constructor for point objects of quasi-static mechanical problems.
				  * @param index Node index.
				  * @param initPos Initial position of the geometry point.
				  * @param Tp Initial temperature.
				  */
				explicit MeshPoint_MQS(const size_t& index, const double* initPos, const double Tp = 0.)
					: MeshPoint<nDim>(), MeshNode_MQS<nDim>(index, initPos, Tp) {
				}

				// Default destructor of private / protected pointers.
				~MeshPoint_MQS() = default;

				// Return a string with current position for printing.
				const std::string print() const override {
					return MeshNode_MQS<nDim>::print();
				}

				// Update trial position
				void updateTrial(const double dPos[]) override {
					for (int i = 0; i < nDim; ++i) {
						this->mv_Ytrial.at(i) = dPos[i];
					}
				};

			protected:
				// Assemble the system of equation, made by the Hessian matrix and a right hand side vector.
				void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override {
					Hessian.clear();

					auto mi_ptElem = this->getElement();

					std::vector<double> mi_pShapeDeriv = mi_ptElem->getShapeDerivOnPoint(this->mv_xsi.data());
					std::vector<double> mi_pShape = mi_ptElem->getShapeFcOnPoint(this->mv_xsi.data());

					// k is related to the conectivity
					int k = 0;
					int mi_nNodes = mi_pShape.size();
					double mi_value;
					std::vector<M2S2::triplet> mi_terms(mi_ptElem->getNumNodes() * mi_ptElem->getNumNodalDof());

					for (auto& node : mi_ptElem->getConectivity()) {
						for (int i = 0; i < nDim; ++i) {
							for (int j = 0; j < nDim; ++j) {
								mi_value = mi_pShapeDeriv.at(k + j * mi_nNodes) * node->getInitPos()[i];
								if ((int)(mi_value * 1000000) != 0) mi_terms.emplace_back(i, j, mi_value);
								//Hessian.coeffRef(i, j) += pShapeDeriv(k, j) * node->getInitPos()[i];
							}
							RHS.at(i) += node->getInitPos()[i] * mi_pShape.at(k);
						}
						k++;
					}
					Hessian.push(mi_terms);

					for (int i = 0; i < nDim; ++i) {
						RHS.at(i) = this->mv_Ycommit.at(i) - RHS.at(i);
					}
				}

				// From MeshPoint, we have:
				//std::array<double, nDim> mv_xsi;					// Dimensionless coordinates
				//O2P2::Proc::Comp::MeshElemVI* mv_ptMeshElem;	// Pointer to the element in which the node is immersed.

				// From MeshNode_MQS, we have:
				//std::vector<size_t> mv_dofList;		// Node DOF list (useless here)
				//std::array<double, nDim> mv_Ytrial;	// Current trial position
				//std::array<double, nDim> mv_Ycommit;	// Current commit position
				//double mv_Tp;							// Room / Commit temperature
			};


			/** @ingroup Processor_Module
			  * @class MeshPoint_MD
			  *
			  * @brief Mesh point, a Solution component, for mechanical dynamic problems.
			  * @details The solution point component for mechanical dynamic problems includes current position (trial and commit), and previous and current (timestep) velocity and accelarion.
			  * Also adds the element in which the point is immersed and dimensionless coordinates related to that element.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshPoint_MD : public MeshPoint<nDim>, public MeshNode_MD<nDim>
			{
			private:
				MeshPoint_MD() = delete;

			public:
				/** Constructor for point objects of dynamic mechanical problems.
				  * @param index Node index.
				  * @param initPos Initial position of the geometry point.
				  * @param Tp Initial temperature.
				  */
				explicit MeshPoint_MD(const size_t& index, const double* initPos, const double Tp = 0.)
					: MeshPoint<nDim>(), MeshNode_MD<nDim>(index, initPos, Tp) {
				}

				// Default destructor of private / protected pointers.
				~MeshPoint_MD() = default;

			protected:
				// Assemble the system of equation, made by the Hessian matrix and a right hand side vector.
				void assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS) override {
					Hessian.clear();

					auto mi_ptElem = this->getElement();

					std::vector<double> mi_pShapeDeriv = mi_ptElem->getShapeDerivOnPoint(this->mv_xsi.data());
					std::vector<double> mi_pShape = mi_ptElem->getShapeFcOnPoint(this->mv_xsi.data());

					// k is related to the conectivity
					int k = 0;
					int mi_nNodes = mi_pShape.size();
					double mi_value;
					std::vector<M2S2::triplet> mi_terms(mi_ptElem->getNumNodes() * mi_ptElem->getNumNodalDof());

					for (auto& node : mi_ptElem->getConectivity()) {
						for (int i = 0; i < nDim; ++i) {
							for (int j = 0; j < nDim; ++j) {
								mi_value = mi_pShapeDeriv.at(k + j * mi_nNodes) * node->getInitPos()[i];
								if ((int)(mi_value * 1000000) != 0) mi_terms.emplace_back(i, j, mi_value);
								//Hessian.coeffRef(i, j) += pShapeDeriv(k, j) * node->getInitPos()[i];
							}
							RHS.at(i) += node->getInitPos()[i] * mi_pShape.at(k);
						}
						k++;
					}
					Hessian.push(mi_terms);

					for (int i = 0; i < nDim; ++i) {
						RHS.at(i) = this->mv_Ycommit.at(i) - RHS.at(i);
					}
				}

				// From MeshPoint, we have:
				//std::array<double, nDim> mv_xsi;				// Dimensionless coordinates
				//O2P2::Proc::Comp::MeshElemVI* mv_ptMeshElem;	// Pointer to the element in which the node is immersed.

				// From MeshNode_MQS, we have:
				//std::vector<size_t> mv_dofList;		// Node DOF list (useless here)
				//std::array<double, nDim> mv_Ytrial;	// Current trial position
				//std::array<double, nDim> mv_Ycommit;	// Current commit position
				//double mv_Tp;							// Room / Commit temperature

				// From MeshNode_MD, we have:
				//std::array<double, nDim> mv_Vp;		// Previous velocity.
				//std::array<double, nDim> mv_Vc;		// Current velocity.
				//std::array<double, nDim> mv_Ap;		// Previous acceleration.
				//std::array<double, nDim> mv_Ac;		// Current acceleration.
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
