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
// Mesh
// 
// Container of processing components: DOF, elements, particles and fibers
// 
// ================================================================================================
#include "Mesh.h"

#include <algorithm>			// for_each
#include <execution>			// parallel for_each

// ================================================================================================
//
// Implementation of Mesh Member Function: initDirichletBC
//
// ================================================================================================
void O2P2::Proc::Mesh::initDirichletBC(const int& loadStep)
{
	LOG("Mesh.initDirichletBC: Initiating Dirichlet BC system for load step " << std::to_string(loadStep));

	// Register the current load step for latter.
	m_curLoadStep = loadStep;

	// Clear all boundary that was once here.
	m_BCIndex.clear();

	// Assign 1 to all DOFs, resizing it.
	m_BCIndex.assign(m_TotalDof, 1);

	// If a Dirichlet boundary condition is imposed, m_BCIndex is changed to zero
	for (auto& DBC : m_LoadStep.at(loadStep)->v_DirichletBC) {
		m_BCIndex[DBC->m_Dof] = 0;
	}
}

// ================================================================================================
//
// Implementation of Mesh Member Function: imposeNeumannBC
//
// ================================================================================================
void O2P2::Proc::Mesh::imposeNeumannBC(Eigen::VectorXd& RHS)
{
	// A new step is beginning thus RHS vector is restarted
	RHS.setZero();

	LOG("Mesh.imposeNeumannBC: Initiating current load vector");

	// If a Dirichlet boundary vector is imposed (like a prescribed displacement), it is first imposed to the RHS 
	for (auto& DBC : m_LoadStep.at(m_curLoadStep)->v_DirichletBC) {
		double value = DBC->m_Value * (DBC->m_TimeVar[0] + DBC->m_TimeVar[1] * m_curTimeStep * m_LoadStep[m_curLoadStep]->m_TimeStep
			+ DBC->m_TimeVar[2] * (m_curTimeStep * m_LoadStep[m_curLoadStep]->m_TimeStep) * (m_curTimeStep * m_LoadStep[m_curLoadStep]->m_TimeStep));

		RHS(DBC->m_Dof) = value;
	}

	// Boundary conditions are saved directly on DOF number
	for (auto& NBC : m_LoadStep[m_curLoadStep]->v_NeumannBC) {
		RHS(NBC->m_Dof) += NBC->m_Value * (NBC->m_TimeVar[0] + NBC->m_TimeVar[1] * m_curTimeStep * m_LoadStep[m_curLoadStep]->m_TimeStep
			+ NBC->m_TimeVar[2] * (m_curTimeStep * m_LoadStep[m_curLoadStep]->m_TimeStep) * (m_curTimeStep * m_LoadStep[m_curLoadStep]->m_TimeStep));
	}

	// Even though there are BC imposed, the load vector must be set ot zero
	//for (size_t i = 0; i < RHS.size(); ++i) {
	//	RHS(i) *= this->m_BCIndex[i];
	//}
}

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Mesh_Mec<2>;
template class O2P2::Proc::Mesh_Mec<3>;

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): assembleSOE
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_Mec<nDim>::assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS)
{
	PROFILE_FUNCTION();

	// Every new iteration requires Hessian to be populated
	Hessian.setZero();

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, m_ElemComp.begin(), m_ElemComp.end(), [this](auto& elem)
		{
			// Matrices sizes
			int nDof = elem->m_ElemDofIndex.size();

			// Element Internal force and Hessian matrix
			Eigen::VectorXd elemVec = Eigen::VectorXd::Zero(nDof);
			Eigen::MatrixXd elemMat = Eigen::MatrixXd::Zero(nDof, nDof);

			// Get element contributions for internal force and hessian matrix
			elem->getContribution(elemVec, elemMat);

			// Impose Dirichlet boundary conditions on current element
			// Eigen matrices are stored in column-major order
			for (int i = 0; i < elemMat.rows(); ++i) {
				size_t dof = elem->m_ElemDofIndex.at(i);

				for (int j = 0; j < elemMat.cols(); ++j) {
					elemMat(j, i) *= this->m_BCIndex.at(dof);
					elemMat(i, j) *= this->m_BCIndex.at(dof);
				}

				elemMat(i, i) += (static_cast<size_t>(1) - this->m_BCIndex.at(dof));
				elemVec(i) *= this->m_BCIndex.at(dof);
			}

			// Register element contribution to local triplet
			for (int i = 0; i < elemMat.cols(); ++i) {
				for (int j = 0; j < elemMat.rows(); ++j) {
					elem->m_elHes.push_back(Eigen::Triplet<double>(static_cast<int>(elem->m_ElemDofIndex.at(j)), static_cast<int>(elem->m_ElemDofIndex.at(i)), elemMat(j, i)));
				}
			}

			// Register element contribution to the right hand side vector
			elemVec.swap(elem->m_elFor);
		}
	);

	// Add element contribution to global system
	std::vector<Eigen::Triplet<double>> gl_triplets;

	for (auto& elem : m_ElemComp) {
		std::move(elem->m_elHes.begin(), elem->m_elHes.end(), std::back_inserter(gl_triplets));
		elem->m_elHes.erase(elem->m_elHes.begin(), elem->m_elHes.end());

		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->m_elFor.size(); ++i) {
			RHS(elem->m_ElemDofIndex.at(i)) -= elem->m_elFor(i);
		}
	}


// I'm getting this way slower than for_each
//#pragma omp parallel for
//	//for (auto& elem : m_ElemComp) {
//	for (std::vector<std::unique_ptr<ElemComp>>::size_type i = 0; i != m_ElemComp.size(); i++)
//	{
//		auto& elem = m_ElemComp[i];
//
//		// Matrices sizes
//		int nDof = elem->m_ElemDofIndex.size();
//
//		// Element Internal force and Hessian matrix
//		Eigen::VectorXd elemVec = Eigen::VectorXd::Zero(nDof);
//		Eigen::MatrixXd elemMat = Eigen::MatrixXd::Zero(nDof, nDof);
//
//		// Get element contributions for internal force and hessian matrix
//		elem->getContribution(elemVec, elemMat);
//
//		// Impose Dirichlet boundary conditions on current element
//		// Eigen matrices are stored in column-major order
//		for (int i = 0; i < elemMat.rows(); ++i) {
//			size_t dof = elem->m_ElemDofIndex.at(i);
//
//			for (int j = 0; j < elemMat.cols(); ++j) {
//				elemMat(j, i) *= this->m_BCIndex.at(dof);
//				elemMat(i, j) *= this->m_BCIndex.at(dof);
//			}
//
//			elemMat(i, i) += (static_cast<size_t>(1) - this->m_BCIndex.at(dof));
//			elemVec(i) *= this->m_BCIndex.at(dof);
//		}
//
//		// Add element contribution to global system
//		for (int i = 0; i < elemMat.cols(); ++i) {
//			for (int j = 0; j < elemMat.rows(); ++j) {
//#pragma omp critical
//				gl_triplets.push_back(Eigen::Triplet<double>(static_cast<int>(elem->m_ElemDofIndex.at(j)), static_cast<int>(elem->m_ElemDofIndex.at(i)), elemMat(j, i)));
//			}
//		}
//
//		// Add element contribution to the right hand side vector
//		for (int i = 0; i < elemVec.size(); ++i) {
//#pragma omp critical
//			RHS(elem->m_ElemDofIndex.at(i)) -= elemVec(i);
//		}
//	}

	Hessian.setFromTriplets(gl_triplets.begin(), gl_triplets.end());
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setTrial
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_Mec<nDim>::setTrial(Eigen::VectorXd& LHS)
{
	//LOG("Mesh_Mec.setTrial: Updating trial solution to nodes");
	double* sol = new double[nDim];
	for (auto& node : (*m_NodePt)) {
		//double* sol = new double[node->v_DofIndex.size()];

		int i = 0;
		for (size_t dof : node->v_DofIndex) {
			sol[i] = LHS(dof);
			i++;
		}

		node->updateTrial(sol);
		//delete[] sol;
	}

	delete[] sol;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setCommit
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_Mec<nDim>::setCommit()
{
	LOG("Mesh_Mec.setCommit: Committing solution to nodes and saving for post-process");

	// Commit the solution and transfer it to Sol (for post-processing).

	// This was made only for displacement, nothing else.
	std::vector<double> Sol((*m_NodePt).size() * nDim);

	size_t i = 0;
	for (auto& node : (*m_NodePt)) {
		node->setCurrent();

		auto x = node->getInitPos();
		auto y = node->getCurrent();

		for (size_t j = 0; j < nDim; ++j) {
			Sol[i] = y[j] - x[j];
			i++;
		}
	}

	// Send the information to the post-process container.
	this->m_PostPt->addSolution(this->m_currentTime, Sol);
}
