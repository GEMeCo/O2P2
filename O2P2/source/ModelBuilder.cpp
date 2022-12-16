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
#include "ModelBuilder.h"

// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template std::unique_ptr<O2P2::Prep::Domain<2>> O2P2::ModelBuilder<2>::initMesh(std::istream& fileProj);
template std::unique_ptr<O2P2::Prep::Domain<3>> O2P2::ModelBuilder<3>::initMesh(std::istream& fileProj);

template std::unique_ptr<O2P2::Proc::SolutionAlgorithm> O2P2::ModelBuilder<2>::initAnalyzer(std::istream& fileProj);
template std::unique_ptr<O2P2::Proc::SolutionAlgorithm> O2P2::ModelBuilder<3>::initAnalyzer(std::istream& fileProj);

template bool O2P2::ModelBuilder<2>::populateDomain(O2P2::Prep::Domain<2>* theDomain);
template bool O2P2::ModelBuilder<3>::populateDomain(O2P2::Prep::Domain<3>* theDomain);

template int O2P2::ModelBuilder<2>::readBondaryConditions(O2P2::Prep::Domain<2>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer);
template int O2P2::ModelBuilder<3>::readBondaryConditions(O2P2::Prep::Domain<3>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer);


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initMesh
//
// ================================================================================================
template<int nDim> std::unique_ptr<O2P2::Prep::Domain<nDim>> O2P2::ModelBuilder<nDim>::initMesh(std::istream& fileProj)
{
	PROFILE_FUNCTION();

	// Line from file
	std::string stLine;
	std::string stFlag = "#ARQUIVOS#";

	// Look up for flags
	LOG("\nModelBuilder.initMesh: Basic Definitions");
	LOG("ModelBuilder.initMesh: Reading flag " << stFlag);

	// Whatever comes before the flag is trashed away
	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(fileProj, stLine);
		if (fileProj.eof()) {
			LOG("\n\nModelBuilder.initMesh: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	// File 1: Nodes Information - For both elements
	std::getline(fileProj, stLine);
	LOG("ModelBuilder.initMesh: Node file: " << stLine);
	
	m_fileNode.open(stLine);
	if (!m_fileNode) {
		LOG("\n\nModelBuilder.initMesh: Error opening file " << stLine << "\n\n\n");
		throw std::invalid_argument("\n\n\nError opening file: " + stLine + "\n\n\n");
	}

	// File 2: Elements Information - For both elements
	std::getline(fileProj, stLine);
	LOG("ModelBuilder.initMesh: Element file: " << stLine);

	m_fileElem.open(stLine);
	if (!m_fileElem) {
		LOG("\n\nModelBuilder.initMesh: Error opening file " << stLine << "\n\n\n");
		throw std::invalid_argument("\n\n\nError opening file: " + stLine + "\n\n\n");
	}

	// File 3: Material Properties
	std::getline(fileProj, stLine);
	LOG("ModelBuilder.initMesh: Materials file: " << stLine);
	
	m_fileMat.open(stLine);
	if (!m_fileMat) {
		LOG("\n\nModelBuilder.initMesh: Error opening file " << stLine << "\n\n\n");
		throw std::invalid_argument("\n\n\nError opening file: " + stLine + "\n\n\n");
	}

	// File 4: Boundary Conditions
	std::getline(fileProj, stLine);
	LOG("ModelBuilder.initMesh: Boundary Conditions file: " << stLine);
	
	m_fileDOF.open(stLine);
	if (!m_fileDOF) {
		LOG("\n\nModelBuilder.initMesh: Error opening file " << stLine << "\n\n\n");
		throw std::invalid_argument("\n\n\nError opening file: " + stLine + "\n\n\n");
	}

	// Basic definitions from node file
	stFlag = "#DADOS_ND#";
	LOG("\nModelBuilder.initMesh: Reading flag " + stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(m_fileNode, stLine);
		if (m_fileNode.eof()) {
			LOG("\n\nModelBuilder.initMesh: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	size_t nNodes;
	m_fileNode >> nNodes;
	LOG("ModelBuilder.initMesh: Number of Nodes: " << std::to_string(nNodes));

	// Basic definitions from material file
	stFlag = "#MATERIAIS#";
	LOG("\nModelBuilder.initMesh: Reading flag " + stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(m_fileMat, stLine);
		if (m_fileMat.eof()) {
			LOG("\n\nModelBuilder.initMesh: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	size_t nMat;
	m_fileMat >> nMat;
	LOG("ModelBuilder.initMesh: Number of Materials: " << std::to_string(nMat));

	// Basic definitions from elements file
	stFlag = "#DADOS_EL#";
	LOG("\nModelBuilder.initMesh: Reading flag " << stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(m_fileElem, stLine);
		if (m_fileElem.eof()) {
			LOG("\n\nModelBuilder.initMesh: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	size_t nElem, nSec;
	m_fileElem >> nElem >> nSec;
	LOG("ModelBuilder.initMesh: Number of Sections: " << std::to_string(nSec));
	LOG("ModelBuilder.initMesh: Number of Elements: " << std::to_string(nElem));

	return std::make_unique<O2P2::Prep::Domain<nDim>>(nNodes, nElem, nMat, nSec);
}

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initAnalyzer
//
// ================================================================================================
template<int nDim> std::unique_ptr<O2P2::Proc::SolutionAlgorithm> O2P2::ModelBuilder<nDim>::initAnalyzer(std::istream& fileProj)
{
	PROFILE_FUNCTION();

	// Line from file
	std::string stLine;
	std::string stFlag;

	// Basic definitions from project file
	stFlag = "#ANALISE#";

	LOG("\nModelBuilder.initAnalyzer: Defining Analysis");
	LOG("ModelBuilder.initAnalyzer: Reading flag " << stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(fileProj, stLine);
		if (fileProj.eof()) {
			LOG("\n\nModelBuilder.initAnalyzer: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	int iAn, iSl, iPr, MaxIt, MinIt, NumLS;
	double Tol;
	fileProj >> iAn >> iSl >> iPr >> Tol >> MinIt >> MaxIt >> NumLS;

	if (iAn > 3) {
		LOG("\n\nType of analysis unknown\nCheck problem input file\n\n\n");
		throw std::invalid_argument("\n\n\nType of analysis unknown\nCheck problem input file\n\n\n");
	}
	if (iSl != 1) {
		LOG("\n\nMethod of solution unknown\nCheck problem input file\n\n\n");
		throw std::invalid_argument("\n\n\nMethod of solution unknown\nCheck problem input file\n\n\n");
	}

	m_AnalysisType = AnalysisType(iAn - 1);
	m_ProblemType = ProblemType(iPr - 1);

	LOG("ModelBuilder.initAnalyzer: Analysis Type: " << analysisTypeNames[m_AnalysisType]);
	LOG("ModelBuilder.initAnalyzer: Problem Type: " << problemTypeNames[m_ProblemType]);
	LOG("ModelBuilder.initAnalyzer: Solver: " << NLSolverTypeNames[NLSolverType(iSl - 1)]);
	LOG("ModelBuilder.initAnalyzer: Number of Load Steps: " << std::to_string(NumLS));
	LOG("ModelBuilder.initAnalyzer: Minimum number of iterations: " << std::to_string(MinIt));
	LOG("ModelBuilder.initAnalyzer: Maximum number of iterations: " << std::to_string(MaxIt));
	LOG("ModelBuilder.initAnalyzer: Tolerance for NonLinear Process: " << std::scientific << Tol << std::fixed);

	// Stores a temporary unique_ptr of SolutionAlgorithm
	std::unique_ptr<O2P2::Proc::SolutionAlgorithm> theAnalyzer = std::make_unique<O2P2::Proc::SolutionAlgorithm>(m_AnalysisType, NLSolverType(iSl - 1), m_ProblemType, NumLS, MinIt, MaxIt, Tol);

	// Reads Time integration parameters, if any
	stFlag = "#PIT#";

	LOG("\nModelBuilder.initAnalyzer: Time Integration Parameters");
	LOG("ModelBuilder.initAnalyzer: Reading flag " << stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(fileProj, stLine);
		if (fileProj.eof()) {
			LOG("\n\nModelBuilder.initAnalyzer: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	double vl1, vl2, vl3;
	fileProj >> vl1 >> vl2 >> vl3;
	theAnalyzer->SetTSP(vl1, vl2, vl3);

	LOG("ModelBuilder.initAnalyzer: Alfa (1st order): " << vl1);
	LOG("ModelBuilder.initAnalyzer: Beta (2nd order): " << vl2);
	LOG("ModelBuilder.initAnalyzer: Gamma (2nd order): " << vl3);

	return std::move(theAnalyzer);
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): populateDomain
//
// ================================================================================================
template<int nDim> bool O2P2::ModelBuilder<nDim>::populateDomain(O2P2::Prep::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	// The return value _isSorted indicates whether there are empty positions or numbering is not in sequence
	bool isSorted;

	LOG("\nModelBuilder.populateDomain: Populating Geometry");

	// Node file
	isSorted = populateNodes(theDomain);
	m_IsDomainPopulated = isSorted;

	// Material data file
	isSorted = populateMaterials(theDomain);
	m_IsDomainPopulated = m_IsDomainPopulated && isSorted;

	// Element information
	isSorted = populateElements(theDomain);
	m_IsDomainPopulated = m_IsDomainPopulated && isSorted;

	// If nothing went wrong, output a message
	if (m_IsDomainPopulated) {
		std::cout << "\n\n------------------------------------------------------------"
			<< "\nInput domain / mesh data was successfull"
			<< "\n------------------------------------------------------------";
	}
	else {
		LOG("\nModelBuilder.populateDomain: Error while evaluating the Domain.");
		throw std::invalid_argument("\n\n\nSomething went wrong while populating the Domain\nCheck input files\n\n\n");
	}

	return m_IsDomainPopulated;
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): populateNodes
//
// ================================================================================================
template<int nDim> bool O2P2::ModelBuilder<nDim>::populateNodes(O2P2::Prep::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	std::string stLine;
	std::string stFlag;
	size_t index;
	bool isSorted = true;
	std::array<double, nDim> x0;

	stFlag = "#NOS#";

	LOG("ModelBuilder.populateNodes: Populating Nodes");
	LOG("ModelBuilder.populateNodes: Reading flag " << stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(m_fileNode, stLine);
		if (m_fileNode.eof()) {
			LOG("\n\nModelBuilder.populatesNodes: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	m_fileNode >> m_indexingInit;
	LOG("ModelBuilder.populatesNodes: File numbering begins with " << std::to_string(m_indexingInit));
	m_fileNode.unget();

	for (size_t i = 0; i < theDomain->m_nNodes; ++i) {
		m_fileNode >> index;

		if (index != i + m_indexingInit) {
			LOG("\n\nModelBuilder.populatesNodes: Reading Error!\nNodes are not sorted. Last valid position: " << std::to_string(i - 1 + m_indexingInit));
			throw std::invalid_argument("\n\n\nNode numbering is not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + m_indexingInit));
		}

		for (int j = 0; j < nDim; ++j) {
			m_fileNode >> x0[j];
		}

		// Index always is related to container index (which begins with 0)
		index = index - m_indexingInit;
		theDomain->addGeomNode(index, x0);
	}

	m_fileNode.close();

	return isSorted;
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): populateMaterials
//
// ================================================================================================
template<int nDim> bool O2P2::ModelBuilder<nDim>::populateMaterials(O2P2::Prep::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	std::string stLine;
	std::string stFlag;
	size_t index;
	int iAux;
	bool isSorted = true;

	stFlag = "#PARAMETROS#";

	LOG("\nModelBuilder.populateMaterials: Populating Materials");
	LOG("ModelBuilder.populateMaterials: Reading flag " << stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(m_fileMat, stLine);
		if (m_fileMat.eof()) {
			LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		};
	}

	for (size_t i = 0; i < theDomain->m_nMat; ++i) {
		m_fileMat >> index >> iAux;

		if (index != i + m_indexingInit) {
			LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nMaterial numbering is not sorted. Last valid position: " << std::to_string(i - 1 + m_indexingInit));
			throw std::invalid_argument("\n\n\nMaterial numbering is not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + m_indexingInit));
		}

		std::vector<double> _param;

		stLine = "";
		stFlag = "#";
		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(m_fileMat, stLine);
			if (m_fileMat.eof()) {
				LOG("\n\nModelBuilder.populateMaterials: Reading Error!\nFlag " << stFlag << " for input material parameters not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " for input material parameters not found\n\n\n");
			};

			std::istringstream Str(stLine);
			double value;

			while (Str >> value) {
				_param.push_back(value);
			}
		}

		iAux--;
		theDomain->addMaterial(index, MaterialType(iAux), _param);
	}

	m_fileMat.close();

	return isSorted;
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): populateElements
//
// ================================================================================================
template<int nDim> bool O2P2::ModelBuilder<nDim>::populateElements(O2P2::Prep::Domain<nDim>* theDomain)
{
	PROFILE_FUNCTION();

	std::string stLine;
	std::string stFlag;
	int iAux1, iAux2, iAux3;
	size_t index, iAux4, iAux5;
	bool isSorted = true;

	if (theDomain->m_nSec > 0) {

		stFlag = "#SECOES#";

		LOG("\nModelBuilder.populateElements: Populating Cross Sections");
		LOG("ModelBuilder.populateElements: Reading flag " << stFlag);

		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(m_fileElem, stLine);
			if (m_fileElem.eof()) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
			}
		}

		for (size_t i = 0; i < theDomain->m_nSec; ++i)
		{
			m_fileElem >> index;

			if (index != i + m_indexingInit) {
				LOG("\n\nModelBuilder.populateElements: Reading Error!\nSection numbering is not sorted. Last valid position: " << std::to_string(i - 1 + m_indexingInit));
				throw std::invalid_argument("\n\n\nSection numbering are not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + m_indexingInit));
			}

			std::string st;
			std::getline(m_fileElem, st);
			std::istringstream Str(st);

			LOG("ModelBuilder.populateElements: Adding Section " + std::to_string(index));
			theDomain->addSection(Str);
		}
	}

	stFlag = "#ELEMENTOS#";

	LOG("\nModelBuilder.populateElements: Populating Elements");
	LOG("ModelBuilder.populateElements: Reading flag " << stFlag);

	while (stLine.compare(0, stFlag.size(), stFlag)) {
		std::getline(m_fileElem, stLine);
		if (m_fileElem.eof()) {
			LOG("\n\nModelBuilder.populateElements: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
			throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
		}
	}

	for (size_t i = 0; i < theDomain->m_nElem; ++i)
	{
		// Index, Type, Order, Integração, Material, [Section], Conectivity[]
		m_fileElem >> index >> iAux1 >> iAux2 >> iAux3 >> iAux4;

		iAux5 = 0;
		if (iAux1 < 4) m_fileElem >> iAux5;		// Elements of type 1 to 3 requires to input a section
		if (index != i + m_indexingInit) {
			LOG("\n\nModelBuilder.populateElements: Reading Error!\nElement numbering is not sorted. Last valid position: " << std::to_string(i - 1 + m_indexingInit));
			throw std::invalid_argument("\n\n\nElement numbering is not sorted or with empty positions\nLast valid position: " + std::to_string(i - 1 + m_indexingInit));
		}

		// Does input files begins with 1 of 0? Containers always begins with 0.
		iAux4 = iAux4 - m_indexingInit;
		iAux5 = iAux5 - m_indexingInit;

		int iParam = theDomain->addElement(index, iAux1, iAux2, iAux3, iAux4, iAux5);

		// The rest of current line is the element indexing
		std::string st;
		std::getline(m_fileElem, st);
		std::istringstream Indexing(st);
		std::vector<size_t> Conec;
		Conec.reserve(iParam);

		iAux1 = 0;
		while (Indexing >> iAux4) {
			iAux1++;	// How many nodes are in the indexing list?
			iAux4 = iAux4 - m_indexingInit;	 // If indexing begins with 1, reduces 1 to match the container indexing
			Conec.push_back(iAux4);
		}

		if (iAux1 != iParam) {
			LOG("\n\nModelBuilder.populateElements: Reading Error!\nWrong number of nodes in indexing of element " << std::to_string(index));
			throw std::invalid_argument("\n\n\nModelBuilder.populateElements: Reading Error!\nWrong number of nodes in indexing of element " + std::to_string(index) + "\n\n\n");
		}

		theDomain->addElementConect(i, Conec);
	}

	m_fileElem.close();

	return isSorted;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): readBoundaryConditions
//
// ================================================================================================
template<int nDim> int O2P2::ModelBuilder<nDim>::readBondaryConditions(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Proc::SolutionAlgorithm* theAnalyzer)
{
	PROFILE_FUNCTION();

	std::string stLine;
	std::string stFlag;
	int iDir, Aux1;
	int totalSteps = 0;
	size_t index, NumDBC, NumNBC;
	double dbAux, Var[3];

	int LS = theAnalyzer->getNumLoadSteps();

	LOG("\nModelBuilder.readBoundaryConditions: Reading Boundary Conditions");

	for (int i = 0; i < LS; ++i) {

		stFlag = "#FC" + std::to_string(i + 1) + "#";
		LOG("ModelBuilder.readBoundaryConditions: Reading flag " << stFlag);

		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(m_fileDOF, stLine);
			if (m_fileDOF.eof()) {
				LOG("\n\nModelBuilder.readBoundaryConditions: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
			}
		}

		LOG("ModelBuilder.readBoundaryConditions: Adding Load Step " << std::to_string(i+1));

		m_fileDOF >> Aux1 >> dbAux >> NumDBC >> NumNBC;
		theAnalyzer->addLoadStep(Aux1, dbAux, NumDBC, NumNBC);
		totalSteps += Aux1;

		stFlag = "#DIR" + std::to_string(i + 1) + "#";
		LOG("\nModelBuilder.readBoundaryConditions: Reading flag " << stFlag);

		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(m_fileDOF, stLine);
			if (m_fileDOF.eof()) {
				LOG("\n\nModelBuilder.readBoundaryConditions: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
			}
		}

		for (size_t j = 0; j < NumDBC; j++) {
			m_fileDOF >> index >> iDir >> dbAux >> Var[0] >> Var[1] >> Var[2];

			// If index begins with 1, reduces 1 to match the container indexing
			index = index - m_indexingInit;
			if (m_indexingInit != 0 && iDir != 0) iDir = iDir - m_indexingInit;

			LOG("SolutionAlgorithm.addDirichletBC: Node / Direction / Value: " << std::to_string(index) << " / " << std::to_string(iDir) << " / " << std::scientific << dbAux << std::fixed);
			theAnalyzer->addDirichletBC(i, index, iDir, dbAux, Var);
		}

		stFlag = "#NEU" + std::to_string(i + 1) + "#";
		LOG("\nModelBuilder.readBoundaryConditions: Reading flag " << stFlag);

		while (stLine.compare(0, stFlag.size(), stFlag)) {
			std::getline(m_fileDOF, stLine);
			if (m_fileDOF.eof()) {
				LOG("\n\nModelBuilder.readBoundaryConditions: Reading Error!\nFlag " << stFlag << " not found\n\n\n");
				throw std::invalid_argument("\n\n\nReading Error!\nFlag " + stFlag + " not found\n\n\n");
			}
		}

		for (size_t j = 0; j < NumNBC; j++) {
			m_fileDOF >> index >> iDir >> dbAux >> Var[0] >> Var[1] >> Var[2];

			// If index begins with 1, reduces 1 to match the container indexing
			index = index - m_indexingInit;
			if (m_indexingInit != 0 && iDir != 0) iDir = iDir - m_indexingInit;

			LOG("SolutionAlgorithm.addNeumannBC: Node / Direction / Value: " << std::to_string(index) << " / " << std::to_string(iDir) << " / " << std::scientific << dbAux << std::fixed);
			theAnalyzer->addNeumannBC(i, index, iDir, dbAux, Var);
		}
	}

	return totalSteps;
}