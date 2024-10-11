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
#include "MeshBars.h"
#include "MeshPoint.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshFibers
			  *
			  * @brief Fibers element routines.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshFibers : public MeshTruss<nDim>
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				MeshFibers() = delete;

			public:
				/** Constructor for analysis of immersed linear elements.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				explicit MeshFibers(std::shared_ptr<O2P2::Geom::Elem::Element> pElem,
					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect) :
					MeshTruss<nDim>(pElem, conect) { }

				// Default destructor of private / protected pointers.
				~MeshFibers() = default;

				// Inherited from MeshTruss:
				// 
				// Prepare element contribution.
				//void getContribution() override;

				// Inherited from MeshTruss:
				// 
				// Add element mass contribution to the Hessian.
				//void addMassContrib(const double& mult) override;

				// Inherited from MeshTruss:
				// 
				// Populates the material point vector.
				//void setMaterialPoint() override;

				// Populates the indexing (DOF) vector.
				void setIndexing() override;

				// Return the spread matrix associated to current inclusion
				M2S2::MatrixX getSpreadMatrix() override;

				// Inherited from MeshTruss:
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

