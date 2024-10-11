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

// Custom header files
#include "MeshElem.h"
#include "MeshPoint.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshIncl
			  *
			  * @brief Inclusion elements routines.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  * Expands results due immersion. It also evaluates the internal stresses (2nd Piola Kirchhoff).
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshIncl : public MeshElem
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				MeshIncl() = delete;

			public:
				/** Constructor for analysis of embedded elements.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				explicit MeshIncl(std::shared_ptr<O2P2::Geom::Elem::Element> pElem,
					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect) :
					MeshElem(pElem, conect) {}

				// Default destructor of private / protected pointers.
				~MeshIncl() = default;

				// Inherited from MeshElem:
				// 
				// Prepare element contribution.
				//void getContribution() override;

				// Inherited from MeshElem:
				// 
				// Add element mass contribution to the Hessian.
				//void addMassContrib(const double& mult) override;

				// Populates the material point vector.
				void setMaterialPoint() override;

				// Populates the indexing (DOF) vector.
				void setIndexing() override;

				// Return the spread matrix associated to current inclusion
				M2S2::MatrixX getSpreadMatrix() override;

			protected:
				// Prepare element contribution of SVK isotropic material.
				void getContribution_SVK(double section) override;

				/** @brief Prepare 2D element contribution of SVK isotropic material in 3D environment.
				  * @param section Only used for linear and plane elements. Value of cross section.
				  */
				void getRotatedContribution_SVK(double section);

				// Inherited from MeshElem:
				// 
				// int mv_nDof;							// Number of DOF for current element.
				// M2S2::MatrixS mv_elHes;				// Element contribution to Hessian matrix (used for parallelism).
				// std::vector<double> mv_elFor;		// Element contribution to internal force (used for parallelism).
				// std::vector<int> mv_dofIndex;		// Dof numbering.
				// std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> mv_conect;			// Node (Point) indexing.
				// std::shared_ptr<O2P2::Geom::Elem::Element> mv_ptElem;						// Pointer to domain element.
				// std::vector<std::unique_ptr<O2P2::Proc::Comp::MaterialPoint>> mv_matPoint;	// Pointers to material points.
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
