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

// Custom Header Files
#include "Common.h"

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
	Material(const size_t& index) : m_index(index) { };

	/** Basic constructor for material objects.
	  * @param index Material number.
	  * @param Param Vector with the material properties.
	  */
	Material(const size_t& index, const std::vector<double>& Param) : m_index(index) { };

	/** Defaul destructor. */
	virtual ~Material() = default;

public:
	/** Set the material parameters.
	 * @param Param Vector with the material properties.
	*/
	virtual void setParameters(const std::vector<double>& Param) = 0;

	/** @return the material type associated to current object. */
	virtual MaterialType getMaterialType() const = 0;

public:
	/** @brief Index of the material for output purposes */
	size_t m_index;
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
	double getLongitudinalModulus() { return m_E11; };

	/** @return Shear modulus / Transversal elastic modulus. */
	double getTransversalModulus() { return m_G12; };

	/** @return Bulk modulus / Volumetric elastic modulus. */
	double getBulkModulus() { return m_Bulk; };

	/** @return First Lamé parameter. */
	double getLameFirst() { return m_Lambda; };

	/** @return Second Lamé parameter. */
	double getLameSecond() { return m_nu12; };

	/** @return Poisson's ratio. */
	double getPoisson() { return m_nu12; };

	/** @return Density. */
	double getDensity() { return m_Rho; };

	/** @return Damping ratio. */
	double getDamping() { return m_Damp; };

	// Return the material type associated to current object.
	MaterialType getMaterialType() const override { return MaterialType::SVK_ISO; };

	/** Overloading operator << to stream the material properties. */
	friend std::ostream& operator<<(std::ostream& stream, Mat_SVK_ISO const& mat) {
		stream << " Young modulus: ";
		stream << formatScien << mat.m_E11;
		stream << " Poisson ratio: ";
		stream << formatFixed << mat.m_nu12;
		return stream;
	};

private:
	/** @brief Young's modulus / Longitudinal elastic modulus */
	double m_E11;

	/** @brief Poisson's ratio or Second Lamé parameter */
	double m_nu12;

	/** @brief Density */
	double m_Rho;

	/** @brief Damping ratio */
	double m_Damp;

	/** @brief Shear modulus / Transversal elastic modulus */
	double m_G12;

	/** @brief Bulk modulus / Volumetric elastic modulus */
	double m_Bulk;

	/** @brief First Lamé parameter */
	double m_Lambda;
};

