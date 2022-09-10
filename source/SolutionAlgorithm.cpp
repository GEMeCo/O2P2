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
//
// Solution Algorithm
// 
// Aggregation of classes to evaluate the analysis.
// 
// ================================================================================================
#include "SolutionAlgorithm.h"

// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template bool O2P2::Proc::SolutionAlgorithm::initFEM<2>(O2P2::Prep::Domain<2>* theDomain, O2P2::Post::PostProcess* thePost);
template bool O2P2::Proc::SolutionAlgorithm::initFEM<3>(O2P2::Prep::Domain<3>* theDomain, O2P2::Post::PostProcess* thePost);

template void O2P2::Proc::SolutionAlgorithm::runSolution<2>(O2P2::Prep::Domain<2>* theDomain);
template void O2P2::Proc::SolutionAlgorithm::runSolution<3>(O2P2::Prep::Domain<3>* theDomain);

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): initFEModel
//
// ================================================================================================
template<int nDim>
bool O2P2::Proc::SolutionAlgorithm::initFEM(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	LOG("\nSolutionAlgorithm.initFEModel: Basic Definitions");
	LOG("SolutionAlgorithm.initFEModel: Initiating DOF Mapping system");

	// Creates the container of Solution Components
	m_theFEModel = std::make_unique<O2P2::Proc::Mesh_Mec<nDim>>(theDomain, thePost);

	// Add the number of dof for every node to the system
	auto ptNode = theDomain->getNode();
	for (size_t i = 0; i < ptNode.size(); ++i) {

		// The number of DOF is associated to the problem and to type of element.
		// For now, the number of dof is the same as the dimensionality of the problem.
		ptNode[i]->v_DofIndex.reserve(nDim);

		// Should also look for hinge (pinned) connections, different number and type of DOF for each element, Lagrange, etc
		for (int j = 0; j < nDim; ++j) {
			ptNode[i]->v_DofIndex.emplace_back(i * nDim + j);
		}

		m_theFEModel->addDOF(nDim);
	}

	// Now assemble the elemental dof
	for (size_t i = 0; i < m_theFEModel->m_meshElem.size(); i++) {

		// Domain element associated to the current element component
		O2P2::Prep::Elem::Element<nDim>* pElem = theDomain->getElem(i);

		std::vector<size_t> elemDofIndex;
		elemDofIndex.reserve(pElem->getConectivity().size() * pElem->getNumNdDOF());

		for (std::shared_ptr<O2P2::Prep::Node<nDim>>& node : pElem->getConectivity()) {
			for (int i = 0; i < pElem->getNumNdDOF(); ++i) {
				elemDofIndex.emplace_back(node->v_DofIndex.at(i));
			}
		}

		// Transfer the elemental dof to the element component
		m_theFEModel->m_meshElem[i]->m_nDof = elemDofIndex.size();
		m_theFEModel->m_meshElem[i]->m_ElemDofIndex = std::move(elemDofIndex);

		// Set the size for the internal force vector
		m_theFEModel->m_meshElem[i]->m_elFor.resize(m_theFEModel->m_meshElem[i]->m_nDof);
	}

	LOG("SolutionAlgorithm.initFEModel: Total number of DOF: " << std::to_string(m_theFEModel->getNumDof()));

	return true;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): runSolutionAlgorithm
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::SolutionAlgorithm::runSolution(O2P2::Prep::Domain<nDim>* theDomain)
{
	std::cout << "\n\n------------------------------------------------------------"
		<< "\nIniting Solution Algorithm"
		<< "\n------------------------------------------------------------";

	// First loop -> Load Steps
	for (int loadStep = 0; loadStep < m_numLoadSteps; ++loadStep)
	{
		LOG("\nSolutionAlgorithm.runSolutionAlgorithm: Current LoadStep: " << std::to_string(loadStep + 1));
		std::cout << "\nLoad Step = " << loadStep + 1
			<< "\n------------------------------------------------------------";

		// Setup Dirichlet boundary conditions
		m_theFEModel->initDirichletBC(loadStep);

		// Second loop -> Time Stepping
		m_theTimeStep->runTimeLoop(m_theFEModel.get(), m_theNLSolver.get());
	}
}