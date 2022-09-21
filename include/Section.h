// ================================================================================================
//
// This file is part of O2P2, an object oriented environment for the positional FEM
//
// Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved.
// 
// This source code form is subject to the terms of 
// Creative Commons Attribution-NonCommerical 4.0 International license
// 
// ================================================================================================
#pragma once

// Standard libraries
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
			Section(const double& sec) : m_section(sec) {
				m_PlaneState = PlaneStateType::PLANE_STRESS;
			};

			/** Constructor for cross section objects.
			  * @param PS Plane state type.
			  * @param thick Thickness of the plane element.
			  */
			Section(const PlaneStateType& PS, const double& thick)
				: m_section(thick), m_PlaneState(PS) { }

			// Default destructor of private / protected pointers.
			virtual ~Section() = default;

			/** @return Cross sectional area / thickness. */
			double getSection() { return m_section; }

			/** @return Plane state type. */
			PlaneStateType getPS() { return m_PlaneState; }

		private:
			/** @brief Cross sectional area for truss elements and thickness for plane elements. */
			double m_section;

			/** @brief Plane state for plane elements (Stress or Strain Plane State) */
			PlaneStateType m_PlaneState;
		};
	} // End of Prep Namespace
} // End of O2P2 Namespace
