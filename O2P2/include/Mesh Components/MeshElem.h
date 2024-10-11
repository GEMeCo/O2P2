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

// Custom libraries
#include "Element.h"
#include "MeshNode.h"
#include "MaterialPoint.h"

// M2S2 library
#include "M2S2/M2S2.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshElemVI
			  *
			  * @brief Base class for elements, inclusion elements, trusses and fibers. It handles routines and local matrices.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  * It also evaluates the internal stresses.
			  *
			  * @note For each Geometry Element (Domain), there is only one Mesh Element, according to the material type.
			  * @note VI stands for Virtual Interface.
			  */
			class MeshElemVI
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				MeshElemVI() = delete;

			public:
				// Default destructor of private / protected pointers.
				virtual ~MeshElemVI() = default;

				/** @return a reference to the elements Material object. */
				O2P2::Geom::Material* getMaterial() { return mv_ptElem->getMaterial(); };

				/** Prepare element contribution.
				  */
				virtual void getContribution() = 0;

				/** Prepare element inertia contribution to hessian matrix.
				  * @param mult Multiplier of the mass matrix (such as Damping and Density, which are not included).
				  */
				virtual void addMassContrib(const double& mult) = 0;

				/** @brief Populates the material point vector. */
				virtual void setMaterialPoint() = 0;

				/** @brief Populates the indexing (DOF) vector. */
				virtual void setIndexing() = 0;

				/** @return the spread matrix associated to inclusion. */
				virtual M2S2::MatrixX getSpreadMatrix() = 0;

			protected:
				/** Constructor for mesh components for elements.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				explicit MeshElemVI(std::shared_ptr<O2P2::Geom::Elem::Element> pElem,
					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect)
				{
					// Store pointer to geometry element
					mv_ptElem = pElem;
					
					// Sets the number of DOF for current element
					mv_nDof = conect.size() * pElem->getNumNodalDof();

					// Set the size for the internal force vector, Hessian and conectivity
					mv_elHes.resize(mv_nDof);
					mv_elFor.resize(mv_nDof);
					mv_dofIndex.reserve(mv_nDof);
					mv_conect = std::move(conect);

					// Reserves memory for material points
					this->mv_matPoint.reserve(pElem->getNumIP());
				}

				/** @brief Prepare element contribution of SVK isotropic material.
				  * @param section Only used for linear and plane elements. Value of cross section.
				  */
				virtual void getContribution_SVK(double section) = 0;

			public:
				/** @brief Number of DOF for current element. */
				int mv_nDof;

				/** @brief Dof numbering. */
				std::vector<int> mv_dofIndex;

				/** @brief Node indexing. */
				std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> mv_conect;

				/** @brief Element contribution to Hessian matrix (used for parallelism). */
				M2S2::MatrixS mv_elHes;

				/** @brief Element contribution to internal force (used for parallelism). */
				std::vector<double> mv_elFor;

				/** @brief Pointer to domain element. */
				std::shared_ptr<O2P2::Geom::Elem::Element> mv_ptElem;

			protected:
				/** @brief Pointer to material point that holds information in the integration point. */
				std::vector<std::unique_ptr<O2P2::Proc::Comp::MaterialPoint>> mv_matPoint;
			};


			/**
			  * @class MeshElem
			  *
			  * @brief Mesh elements routines and local matrices.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  * It also evaluates the internal stresses (2nd Piola Kirchhoff).
			  */
			class MeshElem : public MeshElemVI
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				MeshElem() = delete;

			public:
				/** Constructor for mesh elements.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				explicit MeshElem(std::shared_ptr<O2P2::Geom::Elem::Element> pElem,
					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect) :
					MeshElemVI(pElem, conect) {}

				// Default destructor of private / protected pointers.
				~MeshElem() = default;

				// Prepare element contribution.
				void getContribution() override;

				// Add element mass contribution to the Hessian.
				void addMassContrib(const double& mult) override;

				// Populates the material point vector.
				void setMaterialPoint() override;

				// Populates the indexing (DOF) vector.
				void setIndexing() override;

				// Return an empty spread matrix, since it is no inclusion
				M2S2::MatrixX getSpreadMatrix() override { return M2S2::MatrixX(); };

			protected:
				// Prepare element contribution of SVK isotropic material.
				void getContribution_SVK(double section) override;
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
