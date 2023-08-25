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
				MaterialPoint() = delete;

			public:
				/** Constructor for material point information recorded in the integration point.
				  * @param Jacobian Reference jacobian matrix - A0 / F0.
				  */
				MaterialPoint(const M2S2::Dyadic2N& jacobian) {
					if (jacobian.rows() == 2 || jacobian.rows() == 3)
					{
						mv_refJacobian = jacobian.inverse();
						mv_Jacobian = jacobian.determinant();
					}
					else
					{
						// Linear element
						//mv_refJacobian = jacobian;

						//for (int i = 0; i < jacobian.size(); ++i) {
						//	mv_Jacobian += jacobian(i) * jacobian(i);
						//}
						//mv_Jacobian = 1. / mv_Jacobian;
					}
				};

				// Default destructor of private / protected pointers.
				~MaterialPoint() = default;

				/** @return the inverse of reference Jacobian matrix (A0i / F0i).
				  */
				M2S2::Dyadic2N& getJacobianMatrix() { return this->mv_refJacobian; }

				/** @return the rotation matrix for 2D elements in 3D environments.
				  */
				M2S2::Dyadic2N& getRotationMatrix() { return this->mv_rot; }

				/** @return the reference Jacobian.
				  */
				double getJacobian() { return mv_Jacobian; }

			protected:
				/** @brief Reference Jacobian */
				double mv_Jacobian{ 0. };

				/** @brief INVERSE of reference Jacobian Matrix / A0 or F0 - Mapping gradient from dimensionless coordinate system to initial position */
				M2S2::Dyadic2N mv_refJacobian;

				/** @brief Rotation matrix for 2D inclusions in a 3D environment. Not used otherwise. */
				M2S2::Dyadic2N mv_rot;
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace

