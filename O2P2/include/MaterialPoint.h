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

// Eigen libraries
#include <Eigen/Dense>

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
				MaterialPoint(const Eigen::MatrixXd& Jacobian) {
					if (Jacobian.rows() == Jacobian.cols())
					{
						// Plane element in 2D, or solid elements in 3D
						m_RefJacobian = Jacobian.inverse();
						m_Jacobian = Jacobian.determinant();
					}
					else
					{
						// Linear element
						// m_RefJacobian should be the initial length in each direction
						m_RefJacobian = Jacobian;

						// and m_Jacobian is dA0
						for (int i = 0; i < Jacobian.rows(); ++i) {
							m_Jacobian += Jacobian(i, 0) * Jacobian(i, 0);
						}
						m_Jacobian = 1. / m_Jacobian;
					}
				}

				/** Constructor for material point information recorded in the integration point.
				  * @param Jacobian Reference jacobian matrix - A0 / F0.
				  * @param rot Rotation matrix.
				  */
				MaterialPoint(const Eigen::MatrixXd& Jacobian, const Eigen::MatrixXd& rot) {
					m_RefJacobian = Jacobian.inverse();
					m_Jacobian = Jacobian.determinant();
					m_Rot = rot;
				}

				// Default destructor of private / protected pointers.
				~MaterialPoint() = default;

				/** @return the inverse of reference Jacobian matrix (A0i / F0i).
				  */ 
				Eigen::MatrixXd& getJacobianMatrix() { return m_RefJacobian; }

				/** @return the rotation matrix for 2D elements in 3D environments.
				  */
				Eigen::MatrixXd& getRotationMatrix() { return m_Rot; }

				/** @return the reference Jacobian.
				  */
				double getJacobian() { return m_Jacobian; }

			protected:
				/** @brief Reference Jacobian */
				double m_Jacobian{ 0. };

				/** @brief INVERSE of reference Jacobian Matrix / A0 or F0 - Mapping gradient from dimensionless coordinate system to initial position */
				Eigen::MatrixXd m_RefJacobian;

				/** @brief Rotation matrix for 2D inclusions in a 3D environment. Not used otherwise. */
				Eigen::MatrixXd m_Rot;
			};

			//std::vector<double> vMaterialPoint;
			// relacionado ao material - plasticidade / dano / AAR / creep / etc

		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
