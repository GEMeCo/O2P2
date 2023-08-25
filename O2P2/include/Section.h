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
	namespace Prep {
		/** @ingroup PreProcessor_Module
		  * @class Section
		  *
		  * @brief Section definitions, a Domain component.
		  * @details Contains the basic definitions of the cross section of linear and plane elements.
		  */
		class Section
		{
		private:
			Section() = delete;

		public:
			/** Constructor for cross section objects.
			  * @param sec Cross sectional area / thickness.
			  */
			Section(const double& sec) : mv_section(sec) {
				mv_PlaneState = PlaneStateType::PLANE_STRESS;
			};

			/** Constructor for cross section objects.
			  * @param PS Plane state type.
			  * @param thick Thickness of the plane element.
			  */
			Section(const PlaneStateType& PS, const double& thick)
				: mv_section(thick), mv_PlaneState(PS) { }

			// Default destructor of private / protected pointers.
			~Section() = default;

			/** @return Cross sectional area / thickness. */
			double getSection() { return mv_section; }

			/** @return Plane state type. */
			PlaneStateType getPS() { return mv_PlaneState; }

		private:
			/** @brief Cross sectional area for truss elements and thickness for plane elements. */
			double mv_section;

			/** @brief Plane state for plane elements (Stress or Strain Plane State) */
			PlaneStateType mv_PlaneState;
		};
	} // End of Prep Namespace
} // End of O2P2 Namespace
