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
#include "Mesh.h"

// ================================================================================================
//
// Implementation of Mesh Member Function: initDirichletBC
//
// ================================================================================================
void O2P2::Proc::Mesh::initDirichletBC(const int& loadStep)
{
	LOG("Mesh.initDirichletBC: Initiating Dirichlet BC system for load step " << std::to_string(loadStep));

	// Register the current load step for latter.
	mv_curLoadStep = loadStep;

	// Clear all boundary that was once here, and assign 1 to all.
	mv_BCindex.clear();
	mv_BCindex.assign(mv_TotalDof, 1);

	// If a Dirichlet boundary condition is imposed, m_BCIndex is changed to zero
	for (auto& DBC : mv_LoadStep.at(loadStep)->mv_DirichletBC) {
		mv_BCindex.at(DBC->mv_Dof) = 0;
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
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::imposeNeumannBC(std::vector<double>& RHS)
{
	LOG("Mesh.imposeNeumannBC: Initiating current load vector");

	double mi_curTime = mv_curTimeStep * mv_LoadStep.at(mv_curLoadStep)->mv_timeStep;

	// If a Dirichlet boundary vector is imposed (like a prescribed displacement), it is first imposed to the RHS 
	for (auto& DBC : mv_LoadStep.at(mv_curLoadStep)->mv_DirichletBC) {
		RHS.at(DBC->mv_Dof) = DBC->mv_value * (DBC->mv_timeVar[0] + DBC->mv_timeVar[1] * mi_curTime + DBC->mv_timeVar[2] * mi_curTime * mi_curTime);
	}

	// Boundary conditions are saved directly on DOF number
	for (auto& NBC : mv_LoadStep.at(mv_curLoadStep)->mv_NeumannBC) {
		RHS.at(NBC->mv_Dof) = NBC->mv_value * (NBC->mv_timeVar[0] + NBC->mv_timeVar[1] * mi_curTime + NBC->mv_timeVar[2] * mi_curTime * mi_curTime);
	}

	// Even though there are BC imposed, the load vector must be set ot zero
	//for (size_t i = 0; i < RHS.size(); ++i) {
	//	RHS.at(i) *= this->m_BCIndex[i];
	//}
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): assembleSOE
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	PROFILE_FUNCTION();

	// Set Hessian to zero
	Hessian.clear();

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, this->mv_meshElem.begin(), this->mv_meshElem.end(), [&](auto& elem)
		{
			// Fastest way to set a vector to zero
			memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

			// The rest, is up to the library
			elem->mv_elHes.clear();
			elem->getContribution();

			// Impose restrictions in every element (faster)
			for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
				int mi_gdof = elem->mv_elIndex.at(i);

				if (!this->mv_BCindex.at(mi_gdof)) {
					for (int j = 0; j < elem->mv_elIndex.size(); ++j) {
						elem->mv_elHes.at(i, j) = 0.;
					}

					elem->mv_elHes.at(i, i) = 1.;
					elem->mv_elFor.at(i) *= this->mv_BCindex.at(mi_gdof);
				}
			}

			// Impose restriction directly in the sparse matrix (down below)
			//for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
			//	int mi_gdof = elem->mv_elIndex.at(i);
			//	elem->mv_elFor.at(i) *= this->mv_BCindex.at(mi_gdof);
			//}
		}
	);

	// Add element contribution to global system
	for (auto& elem : this->mv_meshElem) {
		Hessian.push(elem->mv_elHes, elem->mv_elIndex);

		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
			int mi_gdof = elem->mv_elIndex.at(i);
			RHS.at(mi_gdof) -= elem->mv_elFor.at(i);
		}

	}

	// Once every contribution was added, imposed boundary conditions
	//for (auto& DBC : mv_LoadStep.at(mv_curLoadStep)->mv_DirichletBC) {
	//	Hessian.zeroRowAndColumn(DBC->mv_Dof);
	//}

}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setTrial
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MQS<nDim>::setTrial(std::vector<double>& LHS)
{
	//LOG("Mesh_MQS.setTrial: Updating trial solution to nodes");
	double* sol = new double[nDim];

	for (int k = 0; k < LHS.size(); ++k) {
		LHS.at(k) *= this->mv_BCindex.at(k);
	}

	for (auto& node : mv_meshNode) {

		// Number of DOF associated to node - since it is fixed, using outside loop
		// double* sol = new double[node->v_DofIndex.size()];
		size_t i = 0;
		//for (int j = 0; j < node->m_nDOF; ++j) {
		for (int j = 0; j < nDim; ++j) {
			sol[i] = LHS.at(node->mv_DofIndex + j);
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

	// Commit the solution and transfer it to mi_sol (for post-processing).
	// This was made only for displacement, nothing else.

	std::vector<double> mi_sol(mv_TotalDof);

	size_t i = 0;
	for (auto& node : mv_meshNode) {
		// Set trial position as current (commit position)
		node->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = node->getCurPos();
		for (size_t j = 0; j < node->getNumDOF(); ++j) {
			mi_sol.at(i) = *(y + j);
			i++;
		}
	}

	// Send the information to the post-process container.
	this->m_PostPt->addSolution(this->mv_currentTime, mi_sol);
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
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::setTimeStep(const int& timeStep, const double& beta, const  double& gamma)
{
	// Set the current step
	this->mv_curTimeStep = timeStep;

	// Saves Newmark parameters
	this->mv_beta = beta;
	this->mv_gamma = gamma;

	// time step
	double dt = this->mv_LoadStep.at(this->mv_curLoadStep)->mv_timeStep;

	// Evalutes the dynamic contribution from previous step
	// downcast meshNode to dynamic type
	for (auto& node : this->mv_meshNode) {
		O2P2::Proc::Comp::MeshNode_MD<nDim>* pNode = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(node.get());

		for (int j = 0; j < nDim; j++) {
			mv_Qs.at(pNode->mv_DofIndex + j) = pNode->getCurPos()[j] / (beta * dt * dt) +
				pNode->getPrevVel()[j] / (beta * dt) +
				pNode->getPrevAcc()[j] * (0.5 / beta - 1.);
			mv_Rs.at(pNode->mv_DofIndex + j) = pNode->getPrevVel()[j] +
				pNode->getPrevAcc()[j] * dt * (1. - gamma);
		}
	}
}


// ================================================================================================
//
// Implementation of Mesh Member Function: imposeInertiaBC
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::imposeInertiaBC(std::vector<double>& RHS)
{
	PROFILE_FUNCTION();

	double dt = this->mv_LoadStep[this->mv_curLoadStep]->mv_timeStep;

	std::for_each(std::execution::par, this->mv_meshElem.begin(), this->mv_meshElem.end(),
		[&](auto& elem)
		{
			O2P2::Prep::Material* pMat = elem->getMaterial();

			double Rho = pMat->getDensity();
			double Cm = pMat->getDamping();

			elem->mv_elHes.clear();
			elem->addMassContrib(Rho);

			// Set element forces to zero and then evaluate inertia contribution
			memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

			for (int i = 0; i < elem->mv_conect.size(); ++i) {
				for (int j = 0; j < elem->mv_conect.size(); ++j) {
					int nDof = elem->mv_conect.at(j)->getNumDOF();
					size_t nodeDof = elem->mv_conect.at(j)->mv_DofIndex;

					double mass = elem->mv_elHes(i * nDof, j * nDof);

					for (int k = 0; k < nDof; ++k) {
						double yt = elem->mv_conect.at(j)->getTrialPos()[k];

						elem->mv_elFor.at(i * nDof + k) += mass * mv_Qs.at(nodeDof + k) - mass * yt / (mv_beta * dt * dt)
							- mass * Cm * mv_Rs.at(nodeDof + k) - mass * Cm * yt * mv_gamma / (mv_beta * dt)
							+ mass * Cm * mv_Qs.at(nodeDof + k) * mv_gamma * dt;
					}
				}
			}

		}
	);

	for (auto& elem : this->mv_meshElem) {
		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
			int mi_gdof = elem->mv_elIndex.at(i);
			RHS.at(mi_gdof) += elem->mv_elFor.at(i);
		}
	}
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): assembleSOE
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	PROFILE_FUNCTION();

	// Set Hessian to zero
	Hessian.clear();

	imposeInertiaBC(RHS);

	double dt = this->mv_LoadStep[this->mv_curLoadStep]->mv_timeStep;

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, this->mv_meshElem.begin(), this->mv_meshElem.end(), [&](auto& elem)
		{
			// Fastest way to set a vector to zero
			memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

			elem->mv_elHes.clear();
			elem->getContribution();

			// Add mass contribution
			O2P2::Prep::Material* pMat = elem->getMaterial();

			double Rho = pMat->getDensity();
			double Cm = pMat->getDamping();
			double mult = Rho * (1 + this->mv_gamma * Cm * dt) / (this->mv_beta * dt * dt);

			elem->addMassContrib(mult);

			for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
				int mi_gdof = elem->mv_elIndex.at(i);

				if (!this->mv_BCindex.at(mi_gdof)) {
					for (int j = 0; j < elem->mv_elIndex.size(); ++j) {
						elem->mv_elHes.at(i, j) = 0.;
					}

					elem->mv_elHes.at(i, i) = 1.;
					elem->mv_elFor.at(i) *= this->mv_BCindex.at(mi_gdof);
				}
			}
		}
	);

	// Add element contribution to global system
	for (auto& elem : this->mv_meshElem) {
		Hessian.push(elem->mv_elHes, elem->mv_elIndex);

		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
			int mi_gdof = elem->mv_elIndex.at(i);
			RHS.at(mi_gdof) -= elem->mv_elFor.at(i);
		}
	}
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setTrial
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::setTrial(std::vector<double>& LHS)
{
	//LOG("Mesh_MD.setTrial: Updating trial solution to nodes");
	double* sol = new double[nDim * 3];
	double dt = this->mv_LoadStep[this->mv_curLoadStep]->mv_timeStep;

	for (int k = 0; k < LHS.size(); ++k) {
		LHS.at(k) *= this->mv_BCindex.at(k);
	}

	for (auto& node : this->mv_meshNode) {
		O2P2::Proc::Comp::MeshNode_MD<nDim>* pNode = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(node.get());
		size_t i = 0;
		for (int j = 0; j < nDim; ++j) {
			// position variation
			sol[i] = LHS.at(node->mv_DofIndex + j);
			// updated acceleration
			sol[2 * nDim + i] = (pNode->getTrialPos()[j] + sol[i]) / (mv_beta * dt * dt) - mv_Qs.at(pNode->mv_DofIndex + j);
			// updated velocity
			sol[nDim + i] = pNode->getPrevVel()[j] + dt * (1. - mv_gamma) * pNode->getPrevAcc()[j] + dt * mv_gamma * sol[2 * nDim + i];
			i++;
		}
		pNode->updateTrial(sol);
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
	LOG("Mesh_MD.setCommit: Committing solution to nodes and saving for post-process");

	// Commit the solution and transfer it to mi_sol (for post-processing).
	// This was made only for displacement, nothing else.
	std::vector<double> mi_sol(this->mv_TotalDof);

	size_t i = 0;
	for (auto& node : this->mv_meshNode) {
		// Set trial position as current (commit position)
		node->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = node->getCurPos();
		for (size_t j = 0; j < node->getNumDOF(); ++j) {
			mi_sol.at(i) = *(y + j);
			i++;
		}
	}
	// Send the information to the post-process container.
	this->m_PostPt->addSolution(this->mv_currentTime, mi_sol);
}

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setAccel
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Mesh_MD<nDim>::setAccel()
{
	PROFILE_FUNCTION();
	LOG("Mesh_MD.setAccel: Set initial acceleration");

	// Fourth nested loop - Elements
	std::for_each(std::execution::par, this->mv_meshElem.begin(), this->mv_meshElem.end(), [&](auto& elem)
		{
			// Fastest way to set a vector to zero
			memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

			elem->mv_elHes.clear();
			elem->getContribution();

			// Get element contribution for the initial force
			// Although the local hessiam matrix is evaluated, is not required - Only the internal force is needed (saved in elem->m_elFor)
			elem->mv_elHes.clear();

			// Add mass contribution
			O2P2::Prep::Material* pMat = elem->getMaterial();

			double Rho = pMat->getDensity();
			double Cm = pMat->getDamping();

			elem->addMassContrib(Rho);

			// Evaluate initial velocity field
			std::vector<double> va(elem->mv_nDof);

			int i = 0;
			for (auto& node : elem->mv_conect) {
				O2P2::Proc::Comp::MeshNode_MD<nDim>* pNode = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(node.get());

				for (int j = 0; j < nDim; ++j) {
					va.at(i + j) = pNode->getPrevVel()[j];
				}
				i++;
			}

			elem->mv_elFor = Cm * elem->mv_elHes * va;
		}
	);

	// System of equations to be solved
	int size = this->mv_TotalDof;

	M2S2::sparseMatrix GlMass(size, true);
	std::vector<double> mi_RHS(size);
	std::vector<double> mi_LHS(size);

	// Add element contribution to global system
	for (auto& elem : this->mv_meshElem) {
		//GlMass.push(elem->mv_elHes, elem->mv_elIndex);
		std::vector<M2S2::triplet> mi_triplets;

		for (int i = 0; i < elem->mv_elHes.rows(); i++) {
			for (int j = i; j < elem->mv_elHes.cols(); j++) {
				if (std::abs(elem->mv_elHes.at(i, j)) > std::numeric_limits<double>::epsilon()) {
					if (elem->mv_elIndex.at(i) >= elem->mv_elIndex.at(j)) {
						mi_triplets.emplace_back(elem->mv_elIndex.at(j), elem->mv_elIndex.at(i), elem->mv_elHes.at(i, j));
					}
					else {
						mi_triplets.emplace_back(elem->mv_elIndex.at(i), elem->mv_elIndex.at(j), elem->mv_elHes.at(i, j));
					}
				}
			}
		}
		// Push triplet to the sparse matrix
		GlMass.push(mi_triplets);

		// Add element contribution to the right hand side vector
		for (int i = 0; i < elem->mv_elIndex.size(); ++i) {
			int mi_gdof = elem->mv_elIndex.at(i);
			mi_RHS.at(mi_gdof) += elem->mv_elFor.at(i);
		}

		// Fastest way to set a vector to zero
		memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));
	}

	for (auto& NBC : this->mv_LoadStep.at(0)->mv_NeumannBC) {
		mi_RHS.at(NBC->mv_Dof) += NBC->mv_value * (NBC->mv_timeVar[0]);
	}

	// LHS = GlMass^-1 * RHS;
	M2S2::CSR mi_csrMatrix;
	GlMass.saveAsCSR(mi_csrMatrix);

	{
		PROFILE_SCOPE("setAccel");
		solve(size, mi_csrMatrix, mi_RHS, mi_LHS);
	}

	for (int i = 0; i < this->mv_meshNode.size(); ++i) {
		for (int j = 0; j < this->mv_meshNode.at(i)->getNumDOF(); ++j) {
			O2P2::Proc::Comp::MeshNode_MD<nDim>* pNode = static_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>*>(this->mv_meshNode.at(i).get());

			pNode->setAcc(j, mi_LHS.at(pNode->mv_DofIndex + j) * (this->mv_BCindex.at(pNode->mv_DofIndex + j)));
		}
	}
}
