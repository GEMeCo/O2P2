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
#include "ModelBuilder.h"
#include "SolutionAlgorithm.h"
#include "OutputSystem.h"
#include "PostProcess.h"

namespace O2P2 {
	/** @ingroup Main_Module
	  * @class FEAnalysis
	  *
	  * @brief Singleton that manages the pre-processing, processing and post-processing.
	  * @details This singleton is the manager of the whole process (Pre-processing, processing and Pos-processing).
	  *
	  * It aggregates the following classes: Domain, SolutionAlgorithm, and PostProcessing.
	  * It creates a ModelBuilder object to read files and populates the domain.
	  *
	  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
	  */
	template<int nDim>
	class FEAnalysis
	{
	public:
		/** @return a reference to the singleton FEAnalysis. If it does not exists, creates an instance.
		  *
		  * A typical usage is:
		  * \code
		  * FEAnalysis<nDim>::initComponents(ifstream &file);
		  * FEAnalysis<nDim>::runAnalysis();
		  * FEAnalysis<nDim>::drawResults(string filename);
		  * \endcode
		  */
		static FEAnalysis<nDim>& instance() {
			static FEAnalysis<nDim> s_Instance;
			return s_Instance;
		}

		/** Create a ModelBuilder object to populate the Domain.
		 * @return True if Domain was populated.
		 * @param file Project File.
		 */
		static bool initComponents(std::istream& file) {
			PROFILE_FUNCTION();
			return instance().initComponentsImpl(file);
		}

		/** Asks the SolutionAlgorithm object to make the analysis.
		 * @return True if analysis is completed.
		 * @note There is no check if Domains was populated.
		 */
		static bool runAnalysis() {
			PROFILE_FUNCTION();
			return instance().runAnalysisImpl();
		}

		/** Asks the PostProcess object to output the results.
		 * @return True if results were drawn.
		 * @note It does not check if the analysis was feasible.
		*/
		static bool drawResults(const std::string& ProjectName) {
			PROFILE_FUNCTION();
			return instance().drawResultsImpl(ProjectName);
		}

	private:
		// Copy constructor. It must be deleted to avoid instantiation. Only references can be obtained.
		FEAnalysis(const FEAnalysis<nDim>&) = delete;

		// Move constructor will also be deleted.
		FEAnalysis(FEAnalysis&& other) = delete;

		// Private constructor. It is private because it is a singleton, thus it cannot be instantiated.
		FEAnalysis() { }

		// Default destructor of private / protected pointers.
		~FEAnalysis() = default;

		// Internal implementation of call to ModelBuilder.
		bool initComponentsImpl(std::istream& file);

		// Internal implementation of call to the solution algorithm.
		bool runAnalysisImpl();

		// Internal implementation of call to output system.
		bool drawResultsImpl(const std::string& ProjectName);

	private:
		// Container of Domain objects.
		std::unique_ptr<O2P2::Prep::Domain<nDim>> mv_theDomain;

		// Manages the processing unity.
		std::unique_ptr<O2P2::Proc::SolutionAlgorithm> mv_theAnalyzer;

		// Manages the pos-process. Also buffers solutions.
		std::unique_ptr<O2P2::Post::PostProcess> mv_thePost;
	};


	// ================================================================================================
	//
	// Implementation of Private Template Member Function (2D and 3D): initComponentsImpl
	//
	// ================================================================================================
	template<int nDim>
	inline bool FEAnalysis<nDim>::initComponentsImpl(std::istream& file)
	{
		LOG("\nFEAnalysis.initComponents: Reading project files");

		std::cout << "\n\n------------------------------------------------------------"
			      << "\nReading project files"
			      << "\n------------------------------------------------------------";

		O2P2::ModelBuilder<nDim> mi_theBuilder;

		mv_theDomain = mi_theBuilder.initMesh(file);
		mv_theAnalyzer = mi_theBuilder.initAnalyzer(file);
		mv_thePost = std::make_unique<O2P2::Post::PostProcess>();

		bool mi_initMesh = mi_theBuilder.populateDomain(file, mv_theDomain.get());

		// Once mesh was succefully initiated, initiate the solution algorithm elements
		if (mi_initMesh) mi_initMesh = mv_theAnalyzer->initFEModel(mv_theDomain.get(), mv_thePost.get());

		// Read boundary conditions
		int numSteps = mi_theBuilder.readBondaryConditions(file, mv_theDomain.get(), mv_theAnalyzer.get());
		mv_thePost->setNumSteps(numSteps);

		LOG("\nFEAnalysis.initComponents: Finished reading project files");

		return mi_initMesh;
	}

	// ================================================================================================
	//
	// Implementation of Private Template Member Function (2D and 3D): runAnalysisImpl
	//
	// ================================================================================================
	template<int nDim>
	inline bool FEAnalysis<nDim>::runAnalysisImpl() {
		LOG("\nFEAnalysis.runAnalysis: Starting the analysis");
		mv_theAnalyzer->runSolutionAlgorithm(mv_theDomain.get());
		return true;
	}

	// Internal implementation of call to output system.
	template<int nDim>
	inline bool FEAnalysis<nDim>::drawResultsImpl(const std::string& ProjectName) {
		LOG("\nFEAnalysis.drawResults: Outputting solution");

		O2P2::Post::OutputSystem<nDim> theSketcher(OutputType::OGL, ProjectName);
		theSketcher.draw(mv_theDomain.get(), mv_thePost.get());

		return true;
	}
} // End of O2P2 namespace
