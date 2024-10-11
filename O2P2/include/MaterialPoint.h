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

// Standard libraries
#include<vector>

// Custom libraries
#include "M2S2/M2S2.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MaterialPoint
			  *
			  * @brief Material information in the integration point.
			  * @details Stores informations about damage, plastic strain, and such. There is one material point for each integration point.
			  */
			class MaterialPoint
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				MaterialPoint() = delete;

			public:
				/** Constructor for material point information recorded in the integration point. Use for 2D and 3D elements.
				  * @param jacobian Reference jacobian matrix - F0.
				  */
				explicit MaterialPoint(const M2S2::Dyadic2N& jacobian) {
					mv_refJacobian = std::make_unique<M2S2::Dyadic2N>(jacobian.inverse());
					mv_Jacobian = jacobian.determinant();
				}

				/** Constructor for material point information recorded in the integration point. Use for linear elements (trusses and fibers).
				  * @param nDim Dimensionality of domain, not the element dimensionality.
				  * @param jacobian Reference jacobian matrix - F0.
				  */
				explicit MaterialPoint(const int nDim, const double jacobian[]) {
					for (int i = 0; i < nDim; ++i) {
						mv_Jacobian += jacobian[i] * jacobian[i];
					}
				}

				/** Constructor for material point information recorded in the integration point. Use only for 2D elements in 3D environment.
				  * @param jacobian Reference jacobian matrix - F0.
				  * @param rot Rotation matrix.
				  */
				explicit MaterialPoint(const M2S2::Dyadic2N& jacobian, const M2S2::Dyadic2N& rot) {
					mv_refJacobian = std::make_unique<M2S2::Dyadic2N>(jacobian.inverse());
					mv_Jacobian = jacobian.determinant();
					mv_rot = std::make_unique<M2S2::Dyadic2N>(rot); 
				}

				// Default destructor of private / protected pointers.
				~MaterialPoint() = default;

				/** @return the inverse of reference Jacobian matrix (F0i).
				  */
				M2S2::Dyadic2N& getJacobianMatrix() { return *this->mv_refJacobian; }

				/** @return the rotation matrix for 2D elements in 3D environments.
				  */
				M2S2::Dyadic2N& getRotationMatrix() { return *this->mv_rot; }

				/** @return the reference Jacobian.
				  */
				double getJacobian() { return mv_Jacobian; }

			protected:
				/** @brief Reference Jacobian */
				double mv_Jacobian{ 0. };

				/** @brief INVERSE of reference Jacobian Matrix / F0 - Mapping gradient from dimensionless coordinate system to initial position */
				std::unique_ptr<M2S2::Dyadic2N> mv_refJacobian;

				/** @brief Rotation matrix for 2D inclusions in a 3D environment. Not used otherwise. */
				std::unique_ptr<M2S2::Dyadic2N> mv_rot;
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace

