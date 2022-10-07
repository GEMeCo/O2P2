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
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Mesh_MQS<2>;
template class O2P2::Proc::Mesh_MQS<3>;


// ================================================================================================
//
// Implementation of Mesh Member Function: imposeNeumannBC
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::imposeNeumannBC(Eigen::VectorXd& RHS)
{
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
// Implementation of Template Member Function (2D and 3D): assembleSOE
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS)
{
	PROFILE_FUNCTION();

	// Every new iteration requires Hessian to be populated
	Hessian.setZero();

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, this->m_meshElem.begin(), this->m_meshElem.end(), [&](auto& elem)
		{
			// Matrices sizes
			int nDof = elem->m_nDof;

			// Local element Hessian matrix
			Eigen::MatrixXd elemMat = Eigen::MatrixXd::Zero(nDof, nDof);

			// Get element contributions for internal force and hessian matrix
			elem->m_elFor.setZero();
			elem->getContribution(elemMat);

			// Register element contribution to local triplet
			// Dirichlet Boundary condition is applied as triplet is written (not writing it if BC is applied)
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
		}
	);

	// Add element contribution to global system
	std::vector<Eigen::Triplet<double>> gl_triplets;

	for (auto& elem : this->m_meshElem) {
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
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::setTrial(Eigen::VectorXd& LHS)
{
	//LOG("Mesh_MQS.setTrial: Updating trial solution to nodes");
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
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::setCommit()
{
	LOG("Mesh_MQS.setCommit: Committing solution to nodes and saving for post-process");

	// Commit the solution and transfer it to Sol (for post-processing).

	// This was made only for displacement, nothing else.
	std::vector<double> Sol(m_TotalDof);

	size_t i = 0;
	for (auto& node : m_meshNode) {
		// Set trial position as current (commit position)
		node->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = node->getCurPos();

		for (size_t j = 0; j < node->getNumDOF(); ++j) {
			Sol[i] = *(y + j);
			i++;
		}
	}

	// Send the information to the post-process container.
	this->m_PostPt->addSolution(this->m_currentTime, Sol);
}


// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Mesh_MD<2>;
template class O2P2::Proc::Mesh_MD<3>;


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setTimeStep
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::setTimeStep(const int& timeStep, const double& beta, const double& gamma)
{
	// Set the current step
	this->m_curTimeStep = timeStep;

	// Saves Newmark parameters
	m_beta = beta;
	m_gamma = gamma;

	// time step
	double dt = this->m_LoadStep[this->m_curLoadStep]->m_TimeStep;

	// Evalutes the dynamic contribution from previous step
	// downcast meshNode to dynamic type
	for (auto& node : this->m_meshNode) {
		O2P2::Proc::Comp::MeshNode_MD<nDim>* pNode = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(node.get());

		for (int j = 0; j < nDim; ++j) {
			v_Qs(pNode->m_DofIndex + j) = pNode->getCurPos()[j] / (beta * dt * dt) + pNode->getPrevVel()[j] / (beta * dt) + pNode->getPrevAcc()[j] * (0.5 / beta - 1.);
			v_Rs(pNode->m_DofIndex + j) = pNode->getPrevVel()[j] + pNode->getPrevAcc()[j] * dt * (1. - gamma);
		}
	}
}


// ================================================================================================
//
// Implementation of Mesh Member Function: imposeNeumannBC
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::imposeNeumannBC(Eigen::VectorXd& RHS)
{
	PROFILE_FUNCTION();

	// time step
	double dt = this->m_LoadStep[this->m_curLoadStep]->m_TimeStep;

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, this->m_meshElem.begin(), this->m_meshElem.end(), 
		// Execute this lambda function in parallel
		[&](auto& elem)
		{
			// Matrices sizes
			int nDof = elem->m_nDof;

			// Local element matrix
			Eigen::MatrixXd elemMat = Eigen::MatrixXd::Zero(nDof, nDof);

			// Pointer to the material, but downcast to SVK material
			O2P2::Prep::Mat_SVK_ISO* pMat = static_cast<O2P2::Prep::Mat_SVK_ISO*>(elem->getMaterial());

			double Rho = pMat->getDensity();
			double Cm = pMat->getDamping();

			// Evaluate element mass matrix
			elem->addMassContrib(Rho, elemMat);

			// Will save element inertia forces in m_elFor
			elem->m_elFor.setZero();

			for (int i = 0; i < elem->v_Conect.size(); ++i) {
				for (int j = 0; j < elem->v_Conect.size(); ++j) {
					
					// Number of DOF per node
					int nDof = elem->v_Conect[j]->getNumDOF();
					size_t nodeDof = elem->v_Conect[j]->m_DofIndex;
					
					for (int k = 0; k < nDof; ++k) {
						double mass = elemMat(i * nDof + k, j * nDof + k);
						double yt = elem->v_Conect[j]->getTrialPos()[3];
						elem->m_elFor(i * nDof + k) += mass * v_Qs(nodeDof + k) - mass * yt / (m_beta * dt * dt) - mass * Cm * v_Rs(nodeDof + k)
							- mass * Cm * yt * m_gamma / (m_beta * dt) + mass * Cm * v_Qs(nodeDof + k) * m_gamma * dt;
					}
				}
			}
		}
	);

	// Got the elemental inertial force, now transfer to RHS
	for (auto& elem : this->m_meshElem) {

		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->v_Conect.size(); ++i) {
			for (int j = 0; j < nDim; ++j) {
				// global dof
				size_t gdof = elem->v_Conect.at(i)->m_DofIndex + j;
				// Local dof
				size_t ldof = i * nDim + j;

				RHS(gdof) = elem->m_elFor(ldof);
			}
		}
	}

	// The rest of the external forces is just the same as Quasi-static
	Mesh_MQS<nDim>::imposeNeumannBC(RHS);
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): assembleSOE
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS)
{
	PROFILE_FUNCTION();

	// Every new iteration requires Hessian to be populated
	Hessian.setZero();

	// time step
	double dt = this->m_LoadStep[this->m_curLoadStep]->m_TimeStep;

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, this->m_meshElem.begin(), this->m_meshElem.end(),
		// Execute this lambda function in parallel
		[&](auto& elem)
		{
			// Matrices sizes
			int nDof = elem->m_nDof;

			// Local element Hessian matrix
			Eigen::MatrixXd elemMat = Eigen::MatrixXd::Zero(nDof, nDof);

			// Get element contributions for internal force and hessian matrix
			elem->m_elFor.setZero();
			elem->getContribution(elemMat);

			// -------------------------------------------------------------------------------------------
			// Add mass contribution (this is only needed in Dynamic processes)
			// This is the difference from Quasi-Static

			// Pointer to the material, but downcast to SVK material
			O2P2::Prep::Mat_SVK_ISO* pMat = static_cast<O2P2::Prep::Mat_SVK_ISO*>(elem->getMaterial());

			double Rho = pMat->getDensity();
			double Cm = pMat->getDamping();
			double mult = Rho * (1 + this->m_gamma * Cm * dt) / (this->m_beta * dt * dt);

			elem->addMassContrib(mult, elemMat);
			// -------------------------------------------------------------------------------------------

			// Register element contribution to local triplet
			// Dirichlet Boundary condition is applied as triplet is written (not writing it if BC is applied)
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
		}
	);

	// Add element contribution to global system
	std::vector<Eigen::Triplet<double>> gl_triplets;

	for (auto& elem : this->m_meshElem) {
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
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::setTrial(Eigen::VectorXd& LHS)
{
	//LOG("Mesh_MQS.setTrial: Updating trial solution to nodes");
	double* sol = new double[nDim*3];

	// time step
	double dt = this->m_LoadStep[this->m_curLoadStep]->m_TimeStep;

	for (auto& node : this->m_meshNode) {
		// Number of DOF associated to node - since is fixed, using outside loop
		//double* sol = new double[node->v_DofIndex.size()];
		O2P2::Proc::Comp::MeshNode_MD<nDim>* pNode = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(node.get());

		int i = 0;
		//for (int j = 0; j < node->m_nDOF; ++j) {
		for (int j = 0; j < nDim; ++j) {
			// New position
			sol[i] = pNode->getCurPos()[j] + LHS(pNode->m_DofIndex + j);
			// Update acceleration
			sol[2 * nDim + i] = sol[i] / (m_beta * dt * dt) - v_Qs(pNode->m_DofIndex + j);
			// Update velocity
			sol[nDim + i] = pNode->getPrevVel()[j] + dt * (1. - m_gamma) * pNode->getPrevAcc()[j] + dt * m_gamma * pNode->getCurAcc()[j];
			i++;
		}

		pNode->updateTrial(sol);
		//delete[] sol;
	}

	delete[] sol;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setCommit
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::setCommit()
{
	LOG("Mesh_MQS.setCommit: Committing solution to nodes and saving for post-process");

	// Commit the solution and transfer it to Sol (for post-processing).

	// This was made only for displacement, nothing else.
	std::vector<double> Sol(this->m_TotalDof);

	size_t i = 0;
	for (auto& node : this->m_meshNode) {
		// Set trial position as current (commit position)
		node->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = node->getCurPos();

		for (size_t j = 0; j < node->getNumDOF(); ++j) {
			Sol[i] = *(y + j);
			i++;
		}
	}

	// Send the information to the post-process container.
	this->m_PostPt->addSolution(this->m_currentTime, Sol);
}
