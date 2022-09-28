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
#include "Domain.h"
#include "SolutionAlgorithm.h"
#include "NonLinearSolver.h"

namespace O2P2 {
	/** @ingroup PreProcessor_Module
	  * @class ModelBuilder
	  *
	  * @brief Reads project files and populates Domain object.
	  * @details This class manages the entire reading process, holding project files and asking Domain to create all components.
	  * Once reading is finished, it can be safely destroyed.
	  *
	  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
	  */
	template<int nDim>
	class ModelBuilder {

	public:
		// Default constructor.
		ModelBuilder() { };

		// Destructor. Should closes all files.
		~ModelBuilder() {
			m_fileNode.close();
			m_fileElem.close();
			m_fileMat.close();
			m_fileDOF.close();
			LOG("ModelBuilder.destructor: All project files were closed");
		};

		/** Reads the project file, defining the type of analysis and initiating Domain component.
		  * @return a pointer to a Domain object already initialized.
		  * @param fileProj Project file (Main file).
		  */
		std::unique_ptr<O2P2::Prep::Domain<nDim>> initMesh(std::istream& fileProj);

		/** Reads the project file, defines the type of analysis and initiates SolutionAlgorithm component.
		 * @return a pointer to a SolutionAlgorithm object, after initialization.
		 * @param fileProj Project file (Main file).
		*/
		std::unique_ptr<O2P2::Proc::SolutionAlgorithm> initAnalyzer(std::istream& fileProj);

		/** Populates nodes, materials, sections, elements, loads ands constrains.
		  * @return True if mesh components were correctly initiated.
		  * @param theDomain Domain container of all mesh components, holded by FEAnalysis.
		  *
		  * @sa FEAnalysis
		  * @sa Domain
		  */
		bool populateDomain(O2P2::Prep::Domain<nDim>* theDomain);

		/** Reads the Neumann and Dirichlet Boundary Conditions for all Load Steps.
		  * @return the total number of steps (load and time steps).
		  * @param theDomain Domain container of all mesh components, holded by FEAnalysis.
		  * @param theAnalyzer SolutionAlgorithm object, that manages the processing unity, holded by FEAnalysis.
		  *
		  * @sa FEAnalysis
		  * @sa Domain
		  * @sa SolutionAlgorithm
		  */
		int readBondaryConditions(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer);

	public:
		/** @brief Indicates whether input begins with 0 or 1. */
		int m_indexingInit{ 0 };

	private:
		// ModelBuilder is unique, therefore it must not be copied. Thus, copy constructor will be deleted.
		ModelBuilder(const ModelBuilder& other) = delete;

		// Move constructor will also be deleted.
		ModelBuilder(ModelBuilder&& other) = delete;

		// Reads node file and creates nodes and points.
		bool populateNodes(O2P2::Prep::Domain<nDim>* theDomain);

		// Reads materials file and creates materials and sections.
		bool populateMaterials(O2P2::Prep::Domain<nDim>* theDomain);

		// Reads elements file and creates elements and inclusions (particles and fibers).
		bool populateElements(O2P2::Prep::Domain<nDim>* theDomain);

	private:
		// File with nodes and points information.
		std::ifstream m_fileNode;

		// File with elements conectivity.
		std::ifstream m_fileElem;

		// File with material properties.
		std::ifstream m_fileMat;

		// File with boundary conditions.
		std::ifstream m_fileDOF;

		// Mesh boolean.
		bool m_IsDomainPopulated{ false };

		// Type of analysis
		AnalysisType m_AnalysisType{ AnalysisType::STATIC };

		// Type of problem
		ProblemType m_ProblemType{ ProblemType::MECHANICAL };
	};
} // End of O2P2 namespace