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
#include "M2S2/M2S2.h"

namespace O2P2 {
	namespace Prep {
		/** @ingroup PreProcessor_Module
		  * @class Material
		  *
		  * @brief Material properties, a Domain component.
		  * @details Contains the basic definitions for material properties related to the current analysis.
		  */
		class Material
		{
		private:
			Material() = delete;

		protected:
			/** Basic constructor for material objects.
			  * @param index Material number.
			  */
			Material(const size_t& index) : mv_index(index) {}

			/** Basic constructor for material objects.
			  * @param index Material number.
			  * @param Param Vector with the material properties.
			  */
			Material(const size_t& index, const std::vector<double>& Param) : mv_index(index) { }

			// Default destructor of private / protected pointers.
			virtual ~Material() = default;

		public:
			/** @return Density. */
			inline double getDensity() { return mv_Rho; }

			/** @return Damping ratio. */
			inline double getDamping() { return mv_Damp; }

			/** @param nDim Element dimensionality (2D or 3D).
			  * @return the constitutive matrix (row major Voigt notation).
			  */
			virtual M2S2::Dyadic4C getConstitutiveMatrix(int nDim, PlaneStateType PS = PlaneStateType::PLANE_STRESS) = 0;

			/** Set the material parameters.
			 * @param Param Vector with the material properties.
			*/
			virtual void setParameters(const std::vector<double>& Param) = 0;

			/** @return the material type associated to current object. */
			virtual MaterialType getMaterialType() const = 0;

		public:
			/** @brief Index of the material for output purposes */
			size_t mv_index;

		protected:
			/** @brief Density */
			double mv_Rho = 1.;

			/** @brief Damping ratio */
			double mv_Damp = 0.;
		};


		/** @ingroup Material
		  * @class Mat_SVK_ISO
		  *
		  * @brief Saint-Venant-Kirchhoff elastic isotropic constitutive model.
		  * @details Material properties for an elastic isotropic material, based on Saint-Venant-Kirchhoff constitutive model.
		  */
		class Mat_SVK_ISO : public Material
		{
		private:
			Mat_SVK_ISO() = delete;

		public:
			~Mat_SVK_ISO() = default;

			/** Constructor for Saint-Venant-Kirchhoff elastic isotropic constitutive model.
			  * @param index Material number.
			  */
			Mat_SVK_ISO(const size_t& index);

			/** Constructor for Saint-Venant-Kirchhoff elastic isotropic constitutive model.
			 * @param index Material number.
			 * @param Param Vector with the material properties  - Param[0] = Young's modulus; Param[1] = Poisson's ratio; Param[2] = Density; Param[3] = Damping.
			*/
			Mat_SVK_ISO(const size_t& index, const std::vector<double>& Param);

			// Set the material parameters: Param[0] = Young's modulus; Param[1] = Poisson's ratio; Param[2] = Density; Param[3] = Damping.
			void setParameters(const std::vector<double>& Param) override;

			/** @return Young's modulus / Longitudinal elastic modulus. */
			inline double getLongitudinalModulus() { return mv_E11; }

			/** @return Shear modulus / Transversal elastic modulus. */
			inline double getTransversalModulus() { return mv_G12; }

			/** @return Bulk modulus / Volumetric elastic modulus. */
			inline double getBulkModulus() { return mv_Bulk; }

			/** @return First Lamé parameter. */
			inline double getLameFirst() { return mv_Lambda; }

			/** @return Second Lamé parameter. */
			inline double getLameSecond() { return mv_nu12; }

			/** @return Poisson's ratio. */
			inline double getPoisson() { return mv_nu12; }

			// Returns the constitutive matrix - row major with Voigh notation
			M2S2::Dyadic4C getConstitutiveMatrix(int nDim, PlaneStateType PS = PlaneStateType::PLANE_STRESS) override;

			// Return the material type associated to current object.
			MaterialType getMaterialType() const override { return MaterialType::SVK_ISO; }

			/** Overloading operator << to stream the material properties. */
			friend std::ostream& operator<<(std::ostream& stream, Mat_SVK_ISO const& mat) {
				stream << " Young modulus: ";
				stream << formatScien << mat.mv_E11;
				stream << " Poisson ratio: ";
				stream << formatFixed << mat.mv_nu12;
				return stream;
			}

		private:
			/** @brief Young's modulus / Longitudinal elastic modulus */
			double mv_E11;

			/** @brief Poisson's ratio or Second Lamé parameter */
			double mv_nu12;

			/** @brief Shear modulus / Transversal elastic modulus */
			double mv_G12;

			/** @brief Bulk modulus / Volumetric elastic modulus */
			double mv_Bulk;

			/** @brief First Lamé parameter */
			double mv_Lambda;
		};

	} // End of Prep Namespace
} // End of O2P2 Namespace

