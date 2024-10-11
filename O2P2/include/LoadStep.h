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

// C++ standard libraries
#include<vector>
#include<memory>

// Custom Header Files
#include "Constraints.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class LoadStep
			  *
			  * @brief Container of boundary conditions.
			  * @details This class holds Neumann and Dirichlet boundary conditions. It also keeps information about the current load step.
			  */
			class LoadStep
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				LoadStep() = delete;

			public:
				/** Constructor of the container of boundary conditions, that holds information about the current load step.
				  * @param numSteps Number of time steps in the current load step.
				  * @param timeStep Variation in time in the current load step.
				  * @param numDBC Number of applied Dirichlet Boundary Conditions.
				  * @param numNBC Number of applied Neumann Boundary Conditions.
				  */
				explicit LoadStep(const int& numSteps, const double& timeStep, const size_t& numDBC, const size_t& numNBC)
					: mv_numSteps(numSteps), mv_timeStep(timeStep) {
					mv_DirichletBC.reserve(numDBC);
					mv_NeumannBC.reserve(numNBC);
				}

				// Default destructor of private / protected pointers.
				~LoadStep() = default;

				/** Adds a Dirichlet Boundary Condition
				  * @param Dof Degree of Freedom with imposed boundary condition.
				  * @param value Value of the boundary condition
				  * @param var Time behavior (var[0] + var[1].t + var[2].t^2 + var[3].sin(var[4] * dt) + var[5].cos(var[6] * dt))
				  */
				void addDirichletBC(const size_t& Dof, const double& value, const double var[]) {
					mv_DirichletBC.push_back(std::make_unique<O2P2::Proc::Comp::Constraint>(Dof, value, var));
				}

				/** Adds a Neumann Boundary Condition on the flux variable
				  * @param Dof Degree of Freedom with imposed boundary condition.
				  * @param value Value of the boundary condition
				  * @param var Time behavior (var[0] + var[1].t + var[2].t^2 + var[3].sin(var[4] * dt) + var[5].cos(var[6] * dt))
				  */
				void addNeumannBC(const size_t& Dof, const double& value, const double var[]) {
					mv_NeumannBC.push_back(std::make_unique<O2P2::Proc::Comp::Constraint>(Dof, value, var));
				}

			public:
				/** @brief Number of time steps in the current load step */
				int mv_numSteps;

				/** @brief Time step */
				double mv_timeStep;

				/** @brief Container of Dirichlet Boundary Conditions (main variable, i.e. temperature) */
				std::vector<std::unique_ptr<O2P2::Proc::Comp::Constraint>> mv_DirichletBC;

				/** @brief Container of Neumann Boundary Conditions (flux variable, i.e. applied force) */
				std::vector<std::unique_ptr<O2P2::Proc::Comp::Constraint>> mv_NeumannBC;
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
