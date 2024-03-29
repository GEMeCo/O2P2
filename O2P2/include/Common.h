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
#include <string>			// standard string
#include <iostream>			// standard input and output stream | also required by operator << overloading
#include <sstream>			// required by istringstream
#include <fstream>			// standard input and output stream in files
#include <iomanip>			// required by std::setprecision / std::setw
#include <stdexcept>		// standard exceptions (try, throw, catch)
#include <memory>			// required by std::shared_ptr
#include <vector>			// required by std::vector
#include <array>			// required by std::array
#include <map>				// required by std::map
#include <algorithm>		// std::max / std::max_element
#include <chrono>			// standard timing functions
#include <limits>			// required by std::numeric_limits

#include <algorithm>		// for_each
#include <execution>		// parallel for_each

// Profiler
#include "Profiler.h"

// This function only exists in debug mode; Must set the _DEBUG pre-conditioner
#ifdef _DEBUG
extern std::ofstream logFile;

#define LOG(x) logFile << x << std::endl
#else
#define LOG(x)
#endif //_DEBUG

// Function to help output
std::ostream& formatFixed(std::ostream& os);
std::ostream& formatScien(std::ostream& os);

/**
  * @enum AnalysisType Type of analysis
  */
enum struct AnalysisType {
	STATIC,							/** @brief Quasi-estatic analysis */
	TRANSIENT_1stORDER,				/** @brief First order transient analysis, such as heat transfer */
	TRANSIENT_2ndORDER_NEWMARK,		/** @brief Second order transient analysis, using Newmark time step integration method */
	TRANSIENT_2ndORDER_HHT_Alpha,	/** @brief Second order transient analysis, using HHT alpha time step integration method */
	EIGENVALUE						/** @brief Eingevalue and Eigenvector analysis */
};

/**
  * @enum ProblemType Type of problem
  */
enum struct ProblemType {
	MECHANICAL,					/** @brief Only mechanical analysis */
	THERMAL,					/** @brief Only Heat transfer analysis */
	COUPLED						/** @brief Coupled thermomechanical analysis */
};

/**
  * @enum NLSolverType Type of Nonlinear solver
  */
enum struct NLSolverType {
	NEWTONRAPHSON,				/** @brief Newton-Raphson Method */
	BFGS,						/** @brief Broyden�Fletcher�Goldfarb�Shanno Method */
	NEWTONWLS,					/** @brief Newton based method with Line Search */
	NEWTONWTR					/** @brief Newton based method with Trust Region */
};

/**
  * @enum MaterialType Material Type
  */
enum struct MaterialType {
	SVK_ISO,					/** @brief Elastic isotropic, based on Saint-Venant-Kirchhof constitutive model */
	SVK_ORT,					/** @brief Elastic orthotropic, based on Saint-Venant-Kirchhof constitutive model */
	SVK_DAMAGE,					/** @brief Isotropic damage model, based on Saint-Venant-Kirchhof constitutive model */
	CONCRETE,					/** @brief Isotropic damage model, based on SVK constitutive model, including creep, shinkage, AAR and DEF */
	STEEL,						/** @brief SVK isotropic linear elastoplastic plus relaxation, presstress and losses */
	RSD_ISO,					/** @brief Elastic isotropic, based on Rivlin-Saunders-D�ster constitutive model */
	RSD_PLASTIC,				/** @brief Elasto-Plastic isotropic, based on Rivlin-Saunders-D�ster constitutive model */
	TEMP_ISO_LIN,				/** @brief Isotropic heat transfer, with linear properties */
	TEMP_ISO_NLIN,				/** @brief Isotropic heat transfer, with nonlinear properties */
	RSD_TP						/** @brief For Thermomechanical problems - RSD_PLASTIC + TEMP_ISO_NLIN */
};

/**
  * @enum PlaneStateType Plane State
  */
enum struct PlaneStateType {
	PLANE_STRESS,				/** @brief Stress Plane State */
	PLANE_STRAIN				/** @brief Strain Plane State */
};

enum struct OutputType {
	OGL,						/** @brief Acadview file format */
	VTU							/** @brief Paraview file format */
};

extern std::map<AnalysisType, std::string> analysisTypeNames;
extern std::map<ProblemType, std::string> problemTypeNames;
extern std::map<NLSolverType, std::string> NLSolverTypeNames;
extern std::map<MaterialType, std::string> MaterialTypeNames;
extern std::map<PlaneStateType, std::string> PlaneStateTypeNames;
extern std::map<OutputType, std::string> OutputTypeExtension;


/**
 * @brief Timer class for basic profiling - it outputs log file
*/
struct Timer
{
	/** @brief Register the timer start. */
	std::chrono::time_point<std::chrono::steady_clock> start;

	/** @brief Register the timer end. */
	std::chrono::time_point<std::chrono::steady_clock> end;

	/** @brief Function caller, just for output. */
	std::string caller;

	/**
	 * @brief Creates a timer. When destructed, outputs duration.
	 * @param text Output used with timer.
	*/
	Timer(std::string text) {
		caller = text;
		start = std::chrono::high_resolution_clock::now();
	}

	// Destructor of private / protected pointers.
	~Timer() {
		end = std::chrono::high_resolution_clock::now();
		auto var = end - start;

		// integral duration: requires duration_cast
		auto int_hs = std::chrono::duration_cast<std::chrono::hours>(var);
		auto int_mn = std::chrono::duration_cast<std::chrono::minutes>(var - int_hs);
		auto int_sc = std::chrono::duration_cast<std::chrono::seconds>(var - int_mn);
		auto int_ml = std::chrono::duration_cast<std::chrono::milliseconds>(var - int_sc);

		LOG("\n" << caller << " took (h:m:s.ms) "
			<< int_hs.count() << ":" << int_mn.count() << ":" << int_sc.count() << "." << int_ml.count());
		std::cout << "\n" << caller << " took (h:m:s.ms) "
			<< int_hs.count() << ":" << int_mn.count() << ":" << int_sc.count() << "." << int_ml.count();
		std::cout << "\n------------------------------------------------------------";
	}
};

// This was included here only for doxygen documentation :)
// 
//! Outer bounds of the implementation
namespace O2P2 {
	//! Preprocessor classes are gathered in Prep namespace
	namespace Prep {
		//! Domain elements are gathered in Elem namespace
		namespace Elem {}
	}
	//! Processor classes are gathered in Proc namespace
	namespace Proc {
		//! Processor component elements to build mesh are gathered in Comp namespace
		namespace Comp {}
	}
	//! Post-Processor classes are gathered in Post namespace
	namespace Post {}
}