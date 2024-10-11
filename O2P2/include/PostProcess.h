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
		private:
			// Default constructor will be deleted - use the explicit constructor below
			PostProcess() = delete;

			// Copy constructor. It must be deleted to avoid instantiation.
			PostProcess(const PostProcess&) = delete;

			// Move constructor will also be deleted.
			PostProcess(PostProcess&& other) = delete;

		public:
			/** The only constructor that there is. It defines the post processing system.
			  * @param outputFrequency Results output frequency.
			  * @param saveProcess Should intermediate results be saved?
			  * @param outputType Type of file for output.
			  */
			explicit PostProcess(int outputFrequency, int saveProcess, int outputType)
				: mv_OutputFreq(outputFrequency) {
				mv_Save = (saveProcess == 0) ? true : false;
				mv_OutputType = OutputType(outputType - 1);
			};

			// Default destructor of private / protected pointers.
			~PostProcess() = default;

			/** Add a new solution for post-processing.
			  * @param curTimeStep Current time step.
			  * @param sol Vector with the solution obtained.
			  */
			void addSolution(double curTimeStep, std::vector<double>& sol) {
				mv_SolOnNode.emplace_back(std::make_pair(curTimeStep, sol));
			}

			/** Set the number of entries for post-processing.
			  * @param entries Number of entries for post-processing.
			  */
			void setNumSteps(int entries) {
				mv_SolOnNode.reserve(entries);
			}

			/** @return the output file type.
			  */
			OutputType getOutputType() { return mv_OutputType; }

			/** @return the output frequency.
			  */
			inline int getOutputFrequency() { return mv_OutputFreq; }

			/** @return whether the process should be saved or not.
			  */
			inline bool getSaveStatus() { return mv_Save; }

		public:
			/** @brief Solution on each node dof, for each time step. */
			std::vector<std::pair<double, std::vector<double>>> mv_SolOnNode;

			/** @brief Solution on each material point dof, for each time step. */
			//std::vector<double> m_SolOnMatPoint;

		private:
			/** Results output frequency */
			int mv_OutputFreq;

			/** Save the entire process? */
			bool mv_Save;

			/** @brief Type of output */
			OutputType mv_OutputType;
		};
	} // End of Post Namespace
} // End of O2P2 Namespace
