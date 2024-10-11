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
	namespace Geom {
		/** @ingroup PreProcessor_Module
		  * @class CrossSection
		  *
		  * @brief CrossSection definitions, a Domain component.
		  * @details Contains the basic definitions of the cross section of linear elements.
		  */
		class CrossSection
		{
		private:
			// Default constructor is deleted. Use explicit constructor only.
			CrossSection() = delete;

		public:
			/** Constructor for cross section objects.
			  * @param sec Cross sectional area.
			  */
			explicit CrossSection(const double& sec) : mv_section(sec) {	};

			// Default destructor of private / protected pointers.
			virtual ~CrossSection() = default;

			/** @return Cross sectional area / thickness. */
			double getSection() { return this->mv_section; }

		private:
			/** @brief Cross sectional area for truss elements and thickness for plane elements. */
			double mv_section;
		};

		/** @class PlaneSection
		  *
		  * @brief PlaneSection definitions, a Domain component.
		  * @details Contains the basic definitions of the cross section of plane elements.
		  */
		class PlaneSection : public CrossSection
		{
		private:
			// Default constructor is deleted. Use explicit constructor only.
			PlaneSection() = delete;

		public:
			/** Constructor for PlaneSection objects.
			  * @param PS Plane state type.
			  * @param thick Thickness of the plane element.
			  */
			explicit PlaneSection(const PlaneStateType& PS, const double& thick) : CrossSection(thick), mv_PlaneState(PS) { }

			/** @return Plane state type. */
			PlaneStateType getPS() { return mv_PlaneState; }

		private:
			/** @brief Plane state for plane elements (Stress or Strain Plane State) */
			PlaneStateType mv_PlaneState;
		};
	} // End of Geom Namespace
} // End of O2P2 Namespace
