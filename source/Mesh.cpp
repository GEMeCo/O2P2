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
	std::for_each(std::execution::par, m_meshElem.begin(), m_meshElem.end(), [&](auto& elem)
		{
			// Matrices sizes
			int nDof = elem->m_nDof;

			// Element Internal force and Hessian matrix
			Eigen::MatrixXd elemMat = Eigen::MatrixXd::Zero(nDof, nDof);

			// Get element contributions for internal force and hessian matrix
			elem->m_elFor.setZero();
			elem->getContribution(elemMat);

			// Impose Dirichlet boundary conditions on current element
			// Eigen matrices are stored in column-major order
			//for (int i = 0; i < elem->v_Conect.size(); ++i) {
			//	for (int j = 0; j < nDim; ++j) {
			//		// Global dof
			//		size_t dof = elem->v_Conect.at(i)->m_DofIndex + j;
			//		// Local dof
			//		size_t ldof = i * nDim + j;

			//		for (int k = 0; k < elem->m_nDof; ++k) {
			//			elemMat(k, ldof) *= this->m_BCIndex.at(dof);
			//			elemMat(ldof, k) *= this->m_BCIndex.at(dof);
			//		}

			//		elemMat(ldof, ldof) += (static_cast<size_t>(1) - this->m_BCIndex.at(dof));
			//		elem->m_elFor(ldof) *= this->m_BCIndex.at(dof);
			//	}
			//}

			//// Register element contribution to local triplet
			//for (int i = 0; i < elem->v_Conect.size(); ++i) {
			//	for (int j = 0; j < nDim; ++j) {
			//		// local dof 1
			//		size_t dof1 = i * nDim + j;

			//		for (int k = 0; k < elem->v_Conect.size(); ++k) {
			//			for (int l = 0; l < nDim; ++l) {
			//				size_t dof2 = k * nDim + l;

			//				elem->m_elHes.push_back(Eigen::Triplet<double>(
			//					static_cast<int>(elem->v_Conect.at(k)->m_DofIndex + l),
			//					static_cast<int>(elem->v_Conect.at(i)->m_DofIndex + j),
			//					elemMat(dof2, dof1)));
			//			}
			//		}
			//	}
			//}

			// Register element contribution to local triplet
			for (int i = 0; i < elem->v_Conect.size(); ++i) {
				for (int j = 0; j < nDim; ++j) {
					// Global dof 1
					size_t gdof1 = elem->v_Conect.at(i)->m_DofIndex + j;
					// local dof 1
					size_t ldof1 = i * nDim + j;

					if (this->m_BCIndex.at(gdof1)) {
						for (int k = 0; k < elem->v_Conect.size(); ++k) {
							for (int l = 0; l < nDim; ++l) {
								// Global dof 2
								size_t gdof2 = elem->v_Conect.at(k)->m_DofIndex + l;
								// local dof 1
								size_t ldof2 = k * nDim + l;

								if (this->m_BCIndex.at(gdof2)) {
									elem->m_elHes.push_back(Eigen::Triplet<double>(
										static_cast<int>(gdof2),
										static_cast<int>(gdof1),
										elemMat(ldof2, ldof1)));
								}
							}
						}
					}
					else
					{
						elem->m_elFor(ldof1) *= this->m_BCIndex.at(gdof1);

						elem->m_elHes.push_back(Eigen::Triplet<double>(
							static_cast<int>(gdof1), static_cast<int>(gdof1), 1.));
					}
				}
			}


			// Register element contribution to the right hand side vector
			//elemVec.swap(elem->m_elFor);
		}
	);

	// Add element contribution to global system
	std::vector<Eigen::Triplet<double>> gl_triplets;

	for (auto& elem : m_meshElem) {
		std::move(elem->m_elHes.begin(), elem->m_elHes.end(), std::back_inserter(gl_triplets));
		elem->m_elHes.clear();

		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->v_Conect.size(); ++i) {
			for (int j = 0; j < nDim; ++j) {
				// global dof
				size_t gdof = elem->v_Conect.at(i)->m_DofIndex + j;
				// Local dof
				size_t ldof = i * nDim + j;

				RHS(gdof) -= elem->m_elFor(ldof);
			}
		}
	}

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
	for (auto& node : m_meshNode) {
		// Number of DOF associated to node - since is fixed, using outside loop
		//double* sol = new double[node->v_DofIndex.size()];

		int i = 0;
		//for (int j = 0; j < node->m_nDOF; ++j) {
		for (int j = 0; j < nDim; ++j) {
			sol[i] = LHS(node->m_DofIndex + j);
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
	std::vector<double> Sol(m_TotalDof);

	size_t i = 0;
	for (auto& node : m_meshNode) {
		// Set trial position as current (commit position)
		node->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = node->getCurrent();

		for (size_t j = 0; j < node->getNumDOF(); ++j) {
			Sol[i] = *(y + j);
			i++;
		}
	}

	// Send the information to the post-process container.
	this->m_PostPt->addSolution(this->m_currentTime, Sol);
}
