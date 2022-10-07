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

// C++ standard libraries
#include <vector>		// required by std::vector
#include <memory>		// required by std::shared_pointer
#include <utility>		// required by std::make_pair
#include <assert.h>		// required by assert

// Custom libraries
#include "Element.h"
#include "MeshNode.h"
#include "MaterialPoint.h"

// Eigen libraries
#include <Eigen/Dense>
#include <Eigen/Sparse>		// required by Triplets
#include <Eigen/Geometry>	// required by .cross (Cross product)

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshElem
			  *
			  * @brief Base class for elements and inclusion elements routines and local matrices.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  * It should also evaluate the internal stresses.
			  *
			  * @note For each Geometry Element, there is only one Mesh Element, according to material type.
			  */
			class MeshElem
			{
			public:
				// Default destructor of private / protected pointers.
				virtual ~MeshElem() = default;

				/** @return a reference to the elements Material object. */
				O2P2::Prep::Material* getMaterial() { return m_pElem->getMaterial(); }

				/** Prepare element contribution.
				  * @param Hessian Element contribution to the hessian matrix.
				  */
				virtual void getContribution(Eigen::MatrixXd& Hessian) = 0;

				/** Prepare element inertia contribution to hessian matrix.
				  * @param mult Multiplier of the mass matrix (such as Damping and Density, which are not included).
				  * @param Hessian Elemental hessian matrix.
				  */
				virtual void addMassContrib(const double& mult, Eigen::MatrixXd& Hessian) = 0;

			protected:
				/** Constructor for mechanical analysis elements, for SVK material model.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				MeshElem(std::shared_ptr<O2P2::Prep::Elem::BaseElement> pElem, std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect) {
					// Store pointer to geometry element
					m_pElem = pElem;
					// Sets the number of DOF for current element
					m_nDof = conect.size() * pElem->getNumNdDOF();

					// Set the size for the internal force vector and conectivity
					m_elFor.resize(m_nDof);
					v_Conect = std::move(conect);
				}

			protected:
				/** @brief Pointer to material point that holds information in the integration point. */
				std::vector<std::unique_ptr<O2P2::Proc::Comp::MaterialPoint>> m_matPoint;

				/** @brief Pointer to domain element. */
				std::shared_ptr<O2P2::Prep::Elem::BaseElement> m_pElem;

			public:
				/** @brief Number of DOF for current element. */
				int m_nDof;

				/** @brief Vector with element indexing. */
				std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> v_Conect;

				/** @brief Vector with element contribution to Hessian matrix (used for parallelism). */
				std::vector<Eigen::Triplet<double>> m_elHes;

				/** @brief Vector with element contribution to internal force (used for parallelism). */
				Eigen::VectorXd m_elFor;
			};


			/**
			  * @class MeshElem_SVK
			  *
			  * @brief Class for elements and inclusion elements routines and local matrices for SVK material model.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  * It should also evaluate the internal stresses (2nd Piola Kirchhoff).
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshElem_SVK : public MeshElem
			{
			private:
				MeshElem_SVK() = delete;

			public:
				/** Constructor for mechanical analysis elements, for SVK material model.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				explicit MeshElem_SVK(std::shared_ptr<O2P2::Prep::Elem::BaseElement> pElem, std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect)
					: MeshElem(pElem, conect) {
					this->m_matPoint.reserve(pElem->getNumIP());
					this->setMaterialPoint();
				}

				// Default destructor of private / protected pointers.
				virtual ~MeshElem_SVK() = default;

				// Prepare element contribution.
				void getContribution(Eigen::MatrixXd& Hessian) override;

				// Add element mass contribution to the Hessian.
				void addMassContrib(const double& mult, Eigen::MatrixXd& Hessian) override;

			protected:
				/** @brief Populates the material point vector. */
				void setMaterialPoint();

				/** @return the constitutive matrix (in Voigt notation).
				  */
				Eigen::MatrixXd getConstitutiveMatrix();

				/** @return Young modulus for linear elements.
				  */
				double getYoungModulus();

				/** Prepare element contribution of SVK_ISO material.
				  * @param Hessian Element contribution to the hessian matrix.
				  *
				  * @tparam nElDim The dimensionality of the element. It can be 1, 2 or 3 (linear, plane or solid).
				  */
				template<int nElDim>
				void getContribution_SVK_ISO(Eigen::MatrixXd& Hessian);
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
