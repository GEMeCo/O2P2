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

// Custom Header Files
#include "Common.h"

namespace O2P2 {
	namespace Post {
		/** @ingroup PostProcessor_Module
		  * @class PostProcess
		  *
		  * @brief Container of solutions, for post-processing purposes.
		  * @details This class will manages data save for post-processing.
		  */
		class PostProcess
		{
		public:
			/** Constructor */
			PostProcess() {};

			// Default destructor of private / protected pointers.
			~PostProcess() = default;

			/** Add a new solution for post-processing.
			  * @param curTimeStep Current time step.
			  * @param sol Vector with the solution obtained.
			  */
			void addSolution(double curTimeStep, std::vector<double>& sol) {
				m_SolOnNode.emplace_back(std::make_pair(curTimeStep, sol));
			}

			/** Set the number of entries for post-processing.
			  * @param entries Number of entries for post-processing.
			  */
			void setNumSteps(int entries) {
				m_SolOnNode.reserve(entries);
			}

		public:
			/** @brief Solution on each node dof, for each time step. */
			std::vector<std::pair<double, std::vector<double>>> m_SolOnNode;

			/** @brief Solution on each material point dof, for each time step. */
			//std::vector<double> m_SolOnMatPoint;
		};
	} // End of Post Namespace
} // End of O2P2 Namespace
