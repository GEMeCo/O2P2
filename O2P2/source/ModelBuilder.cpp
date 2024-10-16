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
#include "ModelBuilder.h"

// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template std::unique_ptr<O2P2::Geom::Domain<2>> O2P2::ModelBuilder<2>::initMesh(std::istream& file);
template std::unique_ptr<O2P2::Geom::Domain<3>> O2P2::ModelBuilder<3>::initMesh(std::istream& file);

template std::unique_ptr<O2P2::Proc::SolutionAlgorithm> O2P2::ModelBuilder<2>::initAnalyzer(std::istream& file);
template std::unique_ptr<O2P2::Proc::SolutionAlgorithm> O2P2::ModelBuilder<3>::initAnalyzer(std::istream& file);

template std::unique_ptr<O2P2::Post::PostProcess> O2P2::ModelBuilder<2>::initPostProcess(std::istream& file);
template std::unique_ptr<O2P2::Post::PostProcess> O2P2::ModelBuilder<3>::initPostProcess(std::istream& file);

template int O2P2::ModelBuilder<2>::readBondaryConditions(std::istream& file, O2P2::Geom::Domain<2>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer);
template int O2P2::ModelBuilder<3>::readBondaryConditions(std::istream& file, O2P2::Geom::Domain<3>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer);

template bool O2P2::ModelBuilder<2>::populateDomain(std::istream& file, O2P2::Geom::Domain<2>* theDomain);
template bool O2P2::ModelBuilder<3>::populateDomain(std::istream& file, O2P2::Geom::Domain<3>* theDomain);

template bool O2P2::ModelBuilder<2>::populateNodes(std::istream& file, O2P2::Geom::Domain<2>* theDomain);
template bool O2P2::ModelBuilder<3>::populateNodes(std::istream& file, O2P2::Geom::Domain<3>* theDomain);

template bool O2P2::ModelBuilder<2>::populateMaterials(std::istream& file, O2P2::Geom::Domain<2>* theDomain);
template bool O2P2::ModelBuilder<3>::populateMaterials(std::istream& file, O2P2::Geom::Domain<3>* theDomain);

template bool O2P2::ModelBuilder<2>::populateElements(std::istream& file, O2P2::Geom::Domain<2>* theDomain);
template bool O2P2::ModelBuilder<3>::populateElements(std::istream& file, O2P2::Geom::Domain<3>* theDomain);

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initMesh
//
// ================================================================================================
template<int nDim>
inline std::unique_ptr<O2P2::Geom::Domain<nDim>> O2P2::ModelBuilder<nDim>::initMesh(std::istream& file)
{
	PROFILE_FUNCTION();
	LOG("\nModelBuilder.initMesh: Basic Definitions");

	std::string mi_stLine;
	std::string mi_stFlag = "#DATA#";

	// Basic definitions
	LOG("\nModelBuilder.initMesh: Reading flag " + mi_stFlag);

	while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
		std::getline(file, mi_stLine);
		if (file.eof()) {
			LOG("\n\nModelBuilder.initMesh: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
		}
	}
	size_t mi_nNodes, mi_nPoints, mi_nElem, mi_nIncl, mi_nMat, mi_nSec;
	file >> mi_nNodes >> mi_nPoints >> mi_nElem >> mi_nIncl >> mi_nMat >> mi_nSec;

	if (mi_nNodes < 1) throw std::invalid_argument("\n\n\nReading Error!\nAt least one node should exist!\n\n\n");
	if (mi_nElem < 1) throw std::invalid_argument("\n\n\nReading Error!\nAt least one element should exist!\n\n\n");
	if (mi_nMat < 1) throw std::invalid_argument("\n\n\nReading Error!\nAt least one material should exist!\n\n\n");

	LOG("ModelBuilder.initMesh: Number of Nodes: " << std::to_string(mi_nNodes));
	LOG("ModelBuilder.initMesh: Number of Points: " << std::to_string(mi_nPoints));
	LOG("ModelBuilder.initMesh: Number of Elements: " << std::to_string(mi_nElem));
	LOG("ModelBuilder.initMesh: Number of Inclusions: " << std::to_string(mi_nIncl));
	LOG("ModelBuilder.initMesh: Number of Materials: " << std::to_string(mi_nMat));
	LOG("ModelBuilder.initMesh: Number of Sections: " << std::to_string(mi_nSec));

	return std::make_unique<O2P2::Geom::Domain<nDim>>(mi_nNodes, mi_nPoints, mi_nElem, mi_nIncl, mi_nMat, mi_nSec);
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initAnalyzer
//
// ================================================================================================
template<int nDim>
std::unique_ptr<O2P2::Proc::SolutionAlgorithm> O2P2::ModelBuilder<nDim>::initAnalyzer(std::istream& file)
{
	std::string mi_stLine;
	std::string mi_stFlag = "#ANALYSIS#";

	LOG("\nModelBuilder.initAnalyzer: Defining Analysis");
	LOG("ModelBuilder.initAnalyzer: Reading flag " << mi_stFlag);

	while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
		std::getline(file, mi_stLine);
		if (file.eof()) {
			LOG("\n\nModelBuilder.initAnalyzer: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
		}
	}

	int mi_An, mi_Sl, mi_Pr, mi_MaxIt, mi_MinIt, mi_NumLS;
	double mi_Tol;
	file >> mi_An >> mi_Sl >> mi_Pr >> mi_Tol >> mi_MinIt >> mi_MaxIt >> mi_NumLS;

	if (mi_An > 3) {
		LOG("\n\nType of analysis unknown\nCheck problem input file\n\n\n");
		throw std::invalid_argument("\n\n\nType of analysis unknown\nCheck problem input file\n\n\n");
	}
	if (mi_Sl != 1) {
		LOG("\n\nMethod of solution unknown\nCheck problem input file\n\n\n");
		throw std::invalid_argument("\n\n\nMethod of solution unknown\nCheck problem input file\n\n\n");
	}

	mv_AnalysisType = AnalysisType(mi_An - 1);
	mv_ProblemType = ProblemType(mi_Pr - 1);

	LOG("ModelBuilder.initAnalyzer: Analysis Type: " << analysisTypeNames[mv_AnalysisType]);
	LOG("ModelBuilder.initAnalyzer: Problem Type: " << problemTypeNames[mv_ProblemType]);
	LOG("ModelBuilder.initAnalyzer: Solver: " << NLSolverTypeNames[NLSolverType(mi_Sl - 1)]);
	LOG("ModelBuilder.initAnalyzer: Number of Load Steps: " << std::to_string(mi_NumLS));
	LOG("ModelBuilder.initAnalyzer: Minimum number of iterations: " << std::to_string(mi_MinIt));
	LOG("ModelBuilder.initAnalyzer: Maximum number of iterations: " << std::to_string(mi_MaxIt));
	LOG("ModelBuilder.initAnalyzer: Tolerance for NonLinear Process: " << std::scientific << mi_Tol << std::fixed);

	// Stores a temporary unique_ptr of SolutionAlgorithm
	std::unique_ptr<O2P2::Proc::SolutionAlgorithm> theAnalyzer = std::make_unique<O2P2::Proc::SolutionAlgorithm>(mv_AnalysisType, NLSolverType(mi_Sl - 1), mv_ProblemType, mi_NumLS, mi_MinIt, mi_MaxIt, mi_Tol);

	if (mv_AnalysisType != AnalysisType::STATIC) {
		mi_stFlag = "#TIP#";

		LOG("ModelBuilder.initAnalyzer: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.initAnalyzer: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		double vl1, vl2, vl3;
		file >> vl1 >> vl2 >> vl3;
		theAnalyzer->SetTSP(vl1, vl2, vl3);

		LOG("ModelBuilder.initAnalyzer: Alfa (1st order): " << vl1);
		LOG("ModelBuilder.initAnalyzer: Beta (2nd order): " << vl2);
		LOG("ModelBuilder.initAnalyzer: Gamma (2nd order): " << vl3);
	}

	return std::move(theAnalyzer);
}

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initPostProcess
//
// ================================================================================================
template<int nDim>
inline std::unique_ptr<O2P2::Post::PostProcess> O2P2::ModelBuilder<nDim>::initPostProcess(std::istream& file)
{
	PROFILE_FUNCTION();
	LOG("\nModelBuilder.initPostProcess: Basic Definitions");

	std::string mi_stLine;
	std::string mi_stFlag = "#OUTPUT#";

	// Basic definitions
	LOG("\nModelBuilder.initPostProcess: Reading flag " + mi_stFlag);

	while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
		std::getline(file, mi_stLine);
		if (file.eof()) {
			LOG("\n\nModelBuilder.initPostProcess: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
		}
	}
	int mi_freq, mi_save, mi_fileType;
	file >> mi_freq >> mi_save >> mi_fileType;

	if (mi_freq < 1) throw std::invalid_argument("\n\n\nReading Error!\nUnexpected value for output frequency!\n\n\n");
	if (mi_save > 1 || mi_save < 0) throw std::invalid_argument("\n\n\nReading Error!\nUnexpected value for saving file!\n\n\n");
	if (mi_fileType != 1) throw std::invalid_argument("\n\n\nReading Error!\nOnly acadview (.ogl) is available!\n\n\n");

	LOG("ModelBuilder.PostProcess: Output Frequency: " << std::to_string(mi_freq));
	LOG("ModelBuilder.PostProcess: Save Process: " << std::to_string(mi_save));
	LOG("ModelBuilder.PostProcess: FileType: " << std::to_string(mi_fileType));

	return std::make_unique<O2P2::Post::PostProcess>(mi_freq, mi_save, mi_fileType);
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): populateDomain
//
// ================================================================================================
template<int nDim>
bool O2P2::ModelBuilder<nDim>::populateDomain(std::istream& file, O2P2::Geom::Domain<nDim>* theDomain)
{
	bool mi_isSorted;

	LOG("\nModelBuilder.populateDomain: Populating Geometry");

	if (!mv_IsDomainPopulated) {
		mi_isSorted = populateMaterials(file, theDomain);
		mv_IsDomainPopulated = mi_isSorted;

		mi_isSorted = populateNodes(file, theDomain);
		mv_IsDomainPopulated = mv_IsDomainPopulated && mi_isSorted;

		mi_isSorted = populateElements(file, theDomain);
		mv_IsDomainPopulated = mv_IsDomainPopulated && mi_isSorted;
	}

	if (mv_IsDomainPopulated) {
		std::cout << "\n\n---------------------------------------------"
			<< "\nInput domain data was successfull"
			<< "\n---------------------------------------------";
	}
	else {
		LOG("\nModelBuilder.populateDomain: Error while evaluating the Domain.");
		throw std::invalid_argument("\n\n\nSomething went wrong while populating the Domain\nCheck input files\n\n\n");
	}

	return mv_IsDomainPopulated;
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): populateMaterials
//
// ================================================================================================
template<int nDim>
bool O2P2::ModelBuilder<nDim>::populateMaterials(std::istream& file, O2P2::Geom::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	std::string mi_stLine;
	std::string mi_stFlag;
	size_t mi_index;
	int mi_iAux1, mi_iAux2;
	double mi_value;

	mi_stFlag = "#PARAMETERS#";

	LOG("\nModelBuilder.populateMaterials: Populating Materials");
	LOG("ModelBuilder.populateMaterials: Reading flag " << mi_stFlag);

	while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
		std::getline(file, mi_stLine);
		if (file.eof()) {
			LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
		}
	}

	file >> mv_indexingInit;
	file.unget();
	if(mv_indexingInit > 1) throw std::invalid_argument("\n\n\nReading Error!\nCheck materials numbering!\n\n\n");
	LOG("ModelBuilder.populatesMaterials: Index numbering begins with " << std::to_string(mv_indexingInit));

	for (size_t i = 0; i < theDomain->getNumMat(); ++i) {
		file >> mi_index >> mi_iAux1;

		if (mi_index != i + mv_indexingInit) {
			LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nMaterial numbering is not sorted. Last valid position: " << std::to_string(i - 1 + mv_indexingInit));
			throw std::invalid_argument("\n\n\nMaterial numbering is not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + mv_indexingInit));
		}

		std::vector<double> mi_param;

		mi_stLine = "";
		mi_stFlag = "#";

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nFlag " << mi_stFlag << " for input material parameters not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " for input material parameters not found\n\n\n");
			}

			std::istringstream mi_string(mi_stLine);

			while (mi_string >> mi_value) {
				mi_param.push_back(mi_value);
			}
		}

		mi_iAux1--;
		theDomain->addMaterial(mi_index, MaterialType(mi_iAux1), mi_param);
	}

	if (theDomain->getNumSec() > 0) {
		mi_stFlag = "#SECTION#";

		LOG("\nModelBuilder.populateMaterials: Populating Cross Sections");
		LOG("ModelBuilder.populateMaterials: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		for (size_t i = 0; i < theDomain->getNumSec(); ++i)
		{
			file >> mi_index >> mi_iAux1 >> mi_value;

			if (mi_index != i + mv_indexingInit) {
				LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nSection numbering is not sorted. Last valid position: " << std::to_string(i - 1 + mv_indexingInit));
				throw std::invalid_argument("\n\n\nSection numbering is not sorted or with empty position\nLast valid position: " + std::to_string(i - mv_indexingInit));
			}

			// mi_iAux1 holds the type of section
			// mi_value holds the cross section value (area or thickness)
			// mi_iAux2 holds the Plane State type, if any
			mi_iAux2 = 0;
			if (mi_value == 2) {
				file >> mi_iAux2;
			}

			LOG("ModelBuilder.populateMaterials: Adding Section " + std::to_string(mi_index));
			theDomain->addSection(mi_iAux1, mi_value, mi_iAux2);
		}
	}

	return true;
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): populateNodes
//
// ================================================================================================
template<int nDim>
bool O2P2::ModelBuilder<nDim>::populateNodes(std::istream& file, O2P2::Geom::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	std::string mi_stLine;
	std::string mi_stFlag;
	size_t mi_index;
	std::array<double, nDim> mi_x0;

	mi_stFlag = "#NODES#";

	LOG("ModelBuilder.populateNodes: Populating Nodes and Points");
	LOG("ModelBuilder.populateNodes: Reading flag " << mi_stFlag);

	while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
		std::getline(file, mi_stLine);
		if (file.eof()) {
			LOG("\n\nModelBuilder.populatesNodes: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
		}
	}

	for (size_t i = 0; i < theDomain->getNumNodes(); ++i) {
		file >> mi_index;

		if (mi_index != i + mv_indexingInit) {
			LOG("\n\nModelBuilder.populatesNodes: Reading Error!\nNodes are not sorted. Last valid position: " << std::to_string(i - 1 + mv_indexingInit));
			throw std::invalid_argument("\n\n\nNode numbering is not sorted or with empty position\nLast valid position: " + std::to_string(i - mv_indexingInit));
		}

		for (unsigned int j = 0; j < nDim; ++j) {
			file >> mi_x0[j];
		}

		// Index always is related to container index (which begins with 0)
		mi_index = mi_index - mv_indexingInit;
		theDomain->addNode(mi_index, mi_x0);
	}

	// Node for particles (points)
	if (theDomain->getNumPoints() != 0) {

		mi_stFlag = "#POINTS#";

		LOG("\nModelBuilder.populateNodes: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.populatesNodes: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		for (size_t i = 0; i < theDomain->getNumPoints(); i++)
		{
			file >> mi_index;

			if (mi_index != i + mv_indexingInit) {
				LOG("\n\nModelBuilder.populatesNodes: Reading Error!\nInclusion nodes are not sorted. Last valid position: " << std::to_string(i - 1 + mv_indexingInit));
				throw std::invalid_argument("\n\n\nInclusion node numbering is not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + mv_indexingInit));
			}

			for (int j = 0; j < nDim; ++j) {
				file >> mi_x0[j];
			}
			// Numbering in output files must begin with 1, thus indexing will be set accordingly
			mi_index = mi_index - mv_indexingInit;
			theDomain->addPoint(mi_index, mi_x0);
		}
	}
	return true;
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): populateElements
//
// ================================================================================================
template<int nDim>
bool O2P2::ModelBuilder<nDim>::populateElements(std::istream& file, O2P2::Geom::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	std::string mi_stLine;
	std::string mi_stFlag;
	size_t mi_index, mi_iAux4, mi_iAux5;
	int mi_iAux1, mi_iAux2, mi_iAux3;
	std::unordered_map<size_t, std::shared_ptr<O2P2::Geom::Elem::Element>> mi_prism;

	mi_stFlag = "#ELEMENTS#";

	LOG("\nModelBuilder.populateElements: Populating Elements");
	LOG("ModelBuilder.populateElements: Reading flag " << mi_stFlag);

	while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
		std::getline(file, mi_stLine);
		if (file.eof()) {
			LOG("\n\nModelBuilder.populateElements: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
		}
	}

	// Populating Elements and Bars (Trusses) Containers
	size_t mi_numElems = theDomain->getNumElems();
	for (size_t i = 0; i < mi_numElems; ++i) {
		file >> mi_index >> mi_iAux1 >> mi_iAux2 >> mi_iAux3 >> mi_iAux4;

		mi_iAux5 = 0;
		if (mi_iAux1 < 4) file >> mi_iAux5;
		if (mi_index != i + mv_indexingInit) {
			LOG("\n\nModelBuilder.populateElements: Reading Error!\nElement numbering is not sorted. Last valid position: " << std::to_string(i - 1 + mv_indexingInit));
			throw std::invalid_argument("\n\n\nElement numbering is not sorted or with empty position\nLast valid position: " + std::to_string(i - mv_indexingInit));
		}

		mi_iAux4 = mi_iAux4 - mv_indexingInit;
		mi_iAux5 = mi_iAux5 - mv_indexingInit;

		int mi_nParam = theDomain->addElement(mi_index, mi_iAux1, mi_iAux2, mi_iAux3, mi_iAux4, mi_iAux5);

		// The rest of current line is the element indexing
		std::string st;
		std::getline(file, st);
		std::istringstream mi_indexing(st);
		std::vector<size_t> mi_conect;
		mi_conect.reserve(mi_nParam);

		mi_iAux2 = 0;
		while (mi_indexing >> mi_iAux4) {
			mi_iAux2++;		// How many nodes are in the indexing list?
			mi_iAux4 = mi_iAux4 - mv_indexingInit;	// if indexing begins with 1, reduces 1 to match the container indexing
			mi_conect.push_back(mi_iAux4);
		}

		if (mi_iAux2 != mi_nParam) {
			LOG("\n\nModelBuilder.populateElements: Reading Error!\nWrong number of nodes in indexing of element " << std::to_string(mi_index));
			throw std::invalid_argument("\n\n\nReading error!\nWrong number of nodes in indexing of element " + std::to_string(mi_index) + "\n\n\n");
		}

		theDomain->addElementConectivity(i, mi_iAux1, mi_conect);

		// If current element is a prism, add to the map
		if(mi_iAux1 == 6) mi_prism[mi_index] = theDomain->getElem().back();
	}

	// Populating Active Faces container
	if (!mi_prism.empty())
	{
		mi_stFlag = "#FACES#";
		mi_iAux5 = 1;

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);

			if (file.eof()) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nFound End of File while looking for Flag " << mi_stFlag << "\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}

			if (!mi_stLine.compare(0, 13, "#INCLUSIONS#") || !mi_stLine.compare(0, 6, "#LS1#")) {
				LOG("\nModelBuilder.populateElements: Found no prism with Active Faces");
				mi_iAux5 = 0;
				break;
			}
		}

		// Found #FACES#
		if (mi_iAux5) {
			std::getline(file, mi_stLine);
			while (!mi_stLine.empty()) {
				if (!mi_stLine.compare(0, 13, "#INCLUSIONS#") || !mi_stLine.compare(0, 6, "#LS1#")) {
					LOG("\nModelBuilder.populateElements: Found an error during Active Faces Input!");
					break;
				}

				std::istringstream mi_stFace(mi_stLine);
				std::getline(file, mi_stLine);

				mi_stFace >> mi_index >> mi_iAux1 >> mi_iAux4 >> mi_iAux5;

				// Indexing correction
				mi_iAux1 -= mv_indexingInit;
				mi_iAux4 -= mv_indexingInit;
				mi_iAux5 -= mv_indexingInit;

				theDomain->addFaceElem(mi_prism[mi_index], mi_iAux1, mi_iAux4, mi_iAux5);
			}
		}
	}

	// Populating particles and fibers containers
	if (theDomain->getNumIncl() != 0) {

		mi_stFlag = "#INCLUSIONS#";

		LOG("\nModelBuilder.populateElements: Populating Inclusion Elements");
		LOG("ModelBuilder.populateElements: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		size_t mi_numIncl = theDomain->getNumIncl();
		for (size_t i = 0; i < mi_numIncl; i++)
		{
			file >> mi_index >> mi_iAux1 >> mi_iAux2 >> mi_iAux3 >> mi_iAux4;

			mi_iAux5 = 0;
			if (mi_iAux1 <= 4) file >> mi_iAux5;
			if (mi_index != i + mv_indexingInit) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nInclusions numbering is not sorted. Last valid position: " << std::to_string(i - 1 + mv_indexingInit));
				throw std::invalid_argument("\n\n\nInclusions element numbering are not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + mv_indexingInit));
			}

			// Does input files begins with 1 of 0? Containers always begins with 0.
			mi_iAux4 -= mv_indexingInit;
			mi_iAux5 -= mv_indexingInit;

			int iParam = theDomain->addInclusion(mi_index, mi_iAux1, mi_iAux2, mi_iAux3, mi_iAux4, mi_iAux5);

			// The rest of current line is the element indexing
			std::string st;
			std::getline(file, st);
			std::istringstream Indexing(st);
			std::vector<size_t> mi_conec;
			mi_conec.reserve(iParam);

			mi_iAux2 = 0;
			while (Indexing >> mi_iAux4) {
				mi_iAux2++;    // How many nodes in the indexing?
				mi_iAux4 -= mv_indexingInit;     // If index begins with 1, reduces 1 to match the container indexing
				mi_conec.push_back(mi_iAux4);
			}

			if (mi_iAux2 != iParam) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nWrong number of nodes in indexing of inclusion " << std::to_string(mi_index));
				throw std::invalid_argument("\n\n\nModelBuilder.populateElements: Reading Error!\nWrong number of nodes in indexing of inclusion " + std::to_string(mi_index) + "\n\n\n");
			}

			theDomain->addInclusionConectivity(i, mi_iAux1, mi_conec);
		}
	}

	// Reduce initial allocated size, since some may have be included as Trusses and Fibers
	theDomain->getElem().shrink_to_fit();
	theDomain->getIncl().shrink_to_fit();

	return true;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): readBoundaryConditions
//
// ================================================================================================
template<int nDim>
int O2P2::ModelBuilder<nDim>::readBondaryConditions(std::istream& file, O2P2::Geom::Domain<nDim>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer)
{
	PROFILE_FUNCTION();

	std::string mi_stLine;
	std::string mi_stFlag;
	int iDir, Aux1;
	int mi_totalSteps = 0;
	size_t index, mi_numDBC, mi_numNBC;
	double dbAux, var[7]{};

	int LS = theAnalyzer->getNumLoadSteps();

	LOG("\nModelBuilder.readBoundaryConditions: Reading Boundary Conditions");

	for (int i = 0; i < LS; ++i) {
		mi_stFlag = "#LS" + std::to_string(i + 1) + "#";
		LOG("ModelBuilder.readBoundaryConditions: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.readBoundaryConditions: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		LOG("ModelBuilder.readBoundaryConditions: Adding Load Step " << std::to_string(i + 1));

		file >> Aux1 >> dbAux >> mi_numDBC >> mi_numNBC;
		theAnalyzer->addLoadStep(Aux1, dbAux, mi_numDBC, mi_numNBC);
		mi_totalSteps += Aux1;

		mi_stFlag = "#DIR" + std::to_string(i + 1) + "#";
		LOG("\nModelBuilder.readBoundaryConditions: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		for (size_t j = 0; j < mi_numDBC; ++j) {
			file >> index >> iDir >> dbAux >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> var[5] >> var[6];

			index -= mv_indexingInit;
			if (mv_indexingInit != 0 && iDir != 0) iDir -= mv_indexingInit;

			LOG("SolutionAlgorithm.addDirichletBC: Node / Direction / Value: " << std::to_string(index) << " / " << std::to_string(iDir) << " / " << std::scientific << dbAux << std::fixed);
			theAnalyzer->addDirichletBC(i, index, iDir, dbAux, var);
		}

		mi_stFlag = "#NEU" + std::to_string(i + 1) + "#";
		LOG("\nModelBuilder.readBoundaryConditions: Reading flag " << mi_stFlag);

		while (mi_stLine.compare(0, mi_stFlag.size(), mi_stFlag)) {
			std::getline(file, mi_stLine);
			if (file.eof()) {
				LOG("\n\nModelBuilder.readBoundaryConditions: Reading Error!\nFlag " << mi_stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + mi_stFlag + " not found\n\n\n");
			}
		}

		for (size_t j = 0; j < mi_numNBC; ++j) {
			file >> index >> iDir >> dbAux >> var[0] >> var[1] >> var[2] >> var[3] >> var[4] >> var[5] >> var[6];

			index -= mv_indexingInit;
			if (mv_indexingInit != 0 && iDir != 0) iDir -= mv_indexingInit;

			LOG("SolutionAlgorithm.addNeumannBC: Node / Direction / Value: " << std::to_string(index) << " / " << std::to_string(iDir) << " / " << std::scientific << dbAux << std::fixed);
			theAnalyzer->addNeumannBC(i, index, iDir, dbAux, var);
		}
	}
	return mi_totalSteps;
}
