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

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class Constraint
			  *
			  * @brief DOF constraint
			  * @details Dirichlet and Neumann boundary conditions
			  */
			class Constraint
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Constraint() = delete;

			public:
				/** Constructor for both Dirichlet and Neumann boundary conditions.
				  * @param Dof Degree of Freedom with imposed boundary condition.
				  * @param value Value of the boundary condition.
				  * @param var Time behavior (var[0] + var[1].t + var[2].t^2 + var[3].sin(var[4] * dt) + var[5].cos(var[6] * dt))
				  */
				explicit Constraint(const size_t& Dof, const double& value, const double var[7]) {
					mv_Dof = Dof;
					mv_value = value;
					mv_timeVar[0] = var[0];
					mv_timeVar[1] = var[1];
					mv_timeVar[2] = var[2];
					mv_timeVar[3] = var[3];
					mv_timeVar[4] = var[4];
					mv_timeVar[5] = var[5];
					mv_timeVar[6] = var[6];
				}

				// Default destructor of private / protected pointers.
				~Constraint() = default;

				/** @brief Degree of Freedom with boundary condition */
				size_t mv_Dof{ 0 };

				/** @brief Value of the boundary condition, if apply */
				double mv_value{ 0. };

				/** @brief Time function: _value * (timeVar(1) + timeVar(2) * dT + timeVar(3) * dt^2 + timeVar[3].SIN(timeVar[4] * dt) + timeVar[5].COS(timeVar[6]) * dt) */
				double mv_timeVar[7]{ 0. };

			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
