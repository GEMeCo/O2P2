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
#include "MeshIncl.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshFace
			  *
			  * @brief Active Face elements routines.
			  * @details Evaluates local matrices and handles material points elements, one for each integration point.
			  * Expands results due immersion. It also evaluates the internal stresses (2nd Piola Kirchhoff).
			  */
			class MeshFace : public MeshIncl<3>
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				MeshFace() = delete;

			public:
				/** Constructor for analysis of embedded elements.
				  * @param pElem Pointer to geometry element.
				  * @param conect Mesh nodes indexing.
				  */
				explicit MeshFace(std::shared_ptr<O2P2::Geom::Elem::Element> pElem,
					std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>>& conect) :
					MeshIncl<3>(pElem, conect) { }

				// Default destructor of private / protected pointers.
				~MeshFace() = default;

				// Inherited from MeshIncl<3>:
				// 
				// Prepare element contribution.
				//void getContribution() override;

				// Inherited from MeshIncl<3>:
				// 
				// Add element mass contribution to the Hessian.
				//void addMassContrib(const double& mult) override;

				// Inherited from MeshIncl<3>:
				// 
				// Populates the material point vector.
				//void setMaterialPoint() override;

				// Populates the indexing (DOF) vector.
				void setIndexing() override;

			//protected:
				// Inherited from MeshIncl<3>:
				// 
				// Prepare element contribution.
				//void getContribution_SVK(double section) override;

			public:
				// Inherited from MeshElemVI:
				// 
				// int mv_nDof;							// Number of DOF for current element.
				// M2S2::MatrixS mv_elHes;				// Element contribution to Hessian matrix (used for parallelism).
				// std::vector<double> mv_elFor;		// Element contribution to internal force (used for parallelism).
				// std::vector<size_t> mv_dofIndex;		// Dof numbering.
				// std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> mv_conect;			// Node indexing.
				// std::shared_ptr<O2P2::Geom::Elem::Element> mv_ptElem;						// Pointer to domain element.
				// std::vector<std::unique_ptr<O2P2::Proc::Comp::MaterialPoint>> mv_matPoint;	// Pointers to material points.
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
