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
#include "Domain.h"
#include "PostProcess.h"

namespace O2P2 {
	namespace Post {
		/** @ingroup PostProcessor_Module
		  * @class OutputSystem
		  *
		  * @brief Base class for output files
		  * @details Prepare files for AcadView (OGL), Paraview (VTU), etc.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  */
		template<int nDim>
		class OutputSystem
		{
		private:
			// Default constructor
			OutputSystem() = delete;

			// OutputSystem is unique. Copy constructor will be deleted.
			OutputSystem(const OutputSystem& other) = delete;

			// OutputSystem is unique. Move constructor will be deleted.
			OutputSystem(OutputSystem&& other) = delete;

		public:
			/** Constructor for OutputSystem object.
			  * @param projectName Name of the project, employed to give name to files.
			  */
			explicit OutputSystem(const std::string& projectName) {
				mv_Project = projectName;
			}

			// Default destructor
			~OutputSystem() = default;

			/** Prepare output.
			  * @param theDomain Reference to the Domain object.
			  * @param thePost Container with solutions for post-process.
			  */
			void draw(O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost) {
				std::ofstream file;

				if (thePost->getOutputType() == OutputType::OGL) {
					file.open(mv_Project + ".post.desloc" + OutputTypeExtension[OutputType::OGL], std::ios::trunc);

					draw_AcadView_Node(file, theDomain, thePost);
				}
			}

		private:
			// Write file for AcadView Visualizer, based on node information (such as displacement)
			void draw_AcadView_Node(std::ofstream& file, O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost);

			// Write file for AcadView Visualizer, based on element information (single value at nodes) 
			void draw_AcadView_Elem(std::ofstream& file, O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost);

		private:
			/** @brief Name of the Project */
			std::string mv_Project;
		};
	} // End of Post Namespace
} // End of O2P2 Namespace

