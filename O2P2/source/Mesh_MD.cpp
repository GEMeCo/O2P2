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
#include "Mesh_MD.h"
#include "NonLinearSolver.h"

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
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::setTimeStep(const int& timeStep, const double& beta, const  double& gamma)
{
	// Set the current step
	this->mv_curTimeStep = timeStep;

	// Saves Newmark parameters
	this->mv_beta = beta;
	this->mv_gamma = gamma;

	// time step
	double dt = this->mv_LoadStep.at(this->mv_curLoadStep)->mv_timeStep;

	// Evaluated acceleration
	std::array<double, nDim> mi_Acc;

	// Evalutes the dynamic contribution from previous step
	// downcast meshNode to dynamic type
	for (auto& node : this->mv_meshNode) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node));
		auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node);

		for (int j = 0; j < nDim; j++) {
			mi_Acc.at(j) = -pNode->getPrevVel()[j] / (beta * dt) - pNode->getPrevAcc()[j] * (0.5 / beta - 1.);
			mv_Qs.at(pNode->mv_dofList.at(j)) = pNode->getCurPos()[j] / (beta * dt * dt) - mi_Acc.at(j);
		}
		pNode->setCurAcc(mi_Acc);
	}

	// Evalutes the dynamic contribution from previous step
	// downcast meshPoint to dynamic type
	for (size_t i = 0; i < this->mv_meshPoint.size(); i++) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(this->mv_meshPoint.at(i)));
		auto pPoint = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(this->mv_meshPoint.at(i));

		for (int j = 0; j < nDim; j++) {
			mi_Acc.at(j) = -pPoint->getPrevVel()[j] / (beta * dt) - pPoint->getPrevAcc()[j] * (0.5 / beta - 1.);
			mv_QsIm.at(i * nDim + j) = pPoint->getCurPos()[j] / (beta * dt * dt) - mi_Acc.at(j);
		}
		pPoint->setCurAcc(mi_Acc);
	}
}


// ================================================================================================
//
// Implementation of Mesh Member Function: imposeInertiaBC
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::imposeInertiaBC(std::vector<double>& RHS)
{
	PROFILE_FUNCTION();

	// --------------------------------------------------------------------------------------------
	// Matrix Elements
	std::for_each(std::execution::par_unseq, this->mv_meshElem.begin(), this->mv_meshElem.end(), [&](auto& elem)
		{ this->ElemIFO(elem.get(), RHS); } );

	// --------------------------------------------------------------------------------------------
	// Bars Elements (Matrix)
	std::for_each(std::execution::par_unseq, this->mv_meshBars.begin(), this->mv_meshBars.end(), [&](auto& elem)
		{ this->ElemIFO(elem.get(), RHS); });

	// --------------------------------------------------------------------------------------------
	// Particle Elements (Inclusions)
	std::for_each(std::execution::par_unseq, this->mv_meshIncl.begin(), this->mv_meshIncl.end(), [&](auto& elem)
		{ this->InclIFO(elem.get(), RHS); });

	// --------------------------------------------------------------------------------------------
	// Fibers Elements (Inclusions)
	std::for_each(std::execution::par_unseq, this->mv_meshFibs.begin(), this->mv_meshFibs.end(), [&](auto& elem)
		{ this->InclIFO(elem.get(), RHS); });

	// --------------------------------------------------------------------------------------------
	// Active Face Elements
	std::for_each(std::execution::par_unseq, this->mv_meshFace.begin(), this->mv_meshFace.end(), [&](auto& elem)
		{ this->ElemIFO(elem.get(), RHS); });
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setTrial
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::setTrial(std::vector<double>& LHS)
{
	//LOG("Mesh_MD.setTrial: Updating trial solution to nodes");
	double* sol = new double[nDim * 3];
	double dt = this->mv_LoadStep[this->mv_curLoadStep]->mv_timeStep;

	for (int k = 0; k < LHS.size(); ++k) {
		LHS.at(k) *= this->mv_BCindex.at(k);
	}

	for (auto& node : this->mv_meshNode) {
		// Downcast the MeshNode
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node));
		auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node);

		for (size_t j = 0; j < nDim; ++j) {
			// position variation
			sol[j] = LHS.at(pNode->mv_dofList.at(j));
			// updated acceleration
			sol[2 * nDim + j] = (pNode->getTrialPos()[j] + sol[j]) / (mv_beta * dt * dt) - mv_Qs.at(pNode->mv_dofList.at(j));
			// updated velocity
			sol[nDim + j] = pNode->getPrevVel()[j] + dt * (1. - mv_gamma) * pNode->getPrevAcc()[j] + dt * mv_gamma * sol[2 * nDim + j];
		}
		pNode->updateTrial(sol);
	}

	//LOG("Mesh_MD.setTrial: Updating trial solution to points");
	for (size_t i = 0; i < this->mv_meshPoint.size(); i++) {
		// Fastest way to set a vector to zero
		memset(&sol[0], 0., nDim * sizeof(double));

		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(this->mv_meshPoint.at(i)));
		auto pPoint = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(this->mv_meshPoint.at(i));

		// Value of shape functions
		std::vector<double> mi_shape = pPoint->getElement()->getShapeFcOnPoint(pPoint->getXsi().data());

		int k = 0;
		for (auto& node : pPoint->getMeshElement()->mv_conect) {
			// Trial position of element node
			const auto Pos = node->getTrialPos();

			// Evaluate contribution following incidence
			for (int i = 0; i < nDim; i++) {
				sol[i] += Pos[i] * mi_shape.at(k);
			}
			k++;
		}

		for (int j = 0; j < nDim; ++j) {
			// updated acceleration
			sol[2 * nDim + j] = (pPoint->getTrialPos()[j] + sol[j]) / (mv_beta * dt * dt) - mv_QsIm.at(i * nDim + j);
			// updated velocity
			sol[nDim + j] = pPoint->getPrevVel()[j] + dt * (1. - mv_gamma) * pPoint->getPrevAcc()[j] + dt * mv_gamma * sol[2 * nDim + j];
		}

		pPoint->updateTrial(sol);
	}

	delete[] sol;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setCommit
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::setCommit()
{
	LOG("Mesh_MD.setCommit: Committing solution to nodes and saving for post-process");

	// Commit the solution and transfer it to mi_sol (for post-processing).
	// This was made only for displacement, nothing else.
	std::vector<double> mi_sol(this->mv_TotalDof + this->mv_meshPoint.size() * nDim);

	size_t i = 0;
	for (auto& node : this->mv_meshNode) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node));
		auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node);

		// Set trial position as current (commit position)
		pNode->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = node->getCurPos();
		for (size_t j = 0; j < node->getNumDOF(); ++j) {
			mi_sol.at(i) = *(y + j);
			i++;
		}
	}

	for (auto& point : this->mv_meshPoint) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(point));
		auto pPoint = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(point);

		// Set trial position as current (commit position)
		pPoint->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = point->getCurPos();
		for (size_t j = 0; j < point->getNumDOF(); ++j) {
			mi_sol.at(i) = *(y + j);
			i++;
		}
	}

	if ((this->mv_curLoadStep == 0) && double(this->mv_curLoadStep / this->m_PostPt->getOutputFrequency()) - int(this->mv_curLoadStep / this->m_PostPt->getOutputFrequency()) < 1.e-6) {
		// Send the information to the post-process container.
		this->m_PostPt->addSolution(this->mv_currentTime, mi_sol);
	}
}

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): setInitAccel
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::setInitAccel()
{
	PROFILE_FUNCTION();
	LOG("Mesh_MD.setInitAccel: Set initial acceleration");

	// System of equations to be solved
	int size = this->mv_TotalDof;

	M2S2::sparseMatrix GlMass(size, true);
	std::vector<double> mi_RHS(size);
	std::vector<double> mi_LHS(size);

	// --------------------------------------------------------------------------------------------
	// Matrix Elements
	std::for_each(std::execution::par_unseq, this->mv_meshElem.begin(), this->mv_meshElem.end(), [&](auto& elem)
		{ this->ElemIAC(elem.get(), GlMass, mi_RHS); });

	// --------------------------------------------------------------------------------------------
	// Bars Elements (Matrix)
	std::for_each(std::execution::par_unseq, this->mv_meshBars.begin(), this->mv_meshBars.end(), [&](auto& elem)
		{ this->ElemIAC(elem.get(), GlMass, mi_RHS); });

	// --------------------------------------------------------------------------------------------
	// Particle Elements (Inclusions)
	std::for_each(std::execution::par_unseq, this->mv_meshIncl.begin(), this->mv_meshIncl.end(), [&](auto& elem)
		{ this->InclIAC(elem.get(), GlMass, mi_RHS); });

	// --------------------------------------------------------------------------------------------
	// Fibers Elements (Inclusions)
	std::for_each(std::execution::par_unseq, this->mv_meshFibs.begin(), this->mv_meshFibs.end(), [&](auto& elem)
		{ this->InclIAC(elem.get(), GlMass, mi_RHS); });

	// --------------------------------------------------------------------------------------------
	// Active Face Elements
	std::for_each(std::execution::par_unseq, this->mv_meshFace.begin(), this->mv_meshFace.end(), [&](auto& elem)
		{ this->ElemIAC(elem.get(), GlMass, mi_RHS); });

	// --------------------------------------------------------------------------------------------
	// Once assembled, impose Neumann BC and solve for initial acceleration
	for (auto& NBC : this->mv_LoadStep.at(0)->mv_NeumannBC) {
		mi_RHS.at(NBC->mv_Dof) += NBC->mv_value * (NBC->mv_timeVar[0]);
	}

	// LHS = GlMass^-1 * RHS;
	M2S2::CSR mi_csrMatrix;
	GlMass.saveAsCSR(mi_csrMatrix);

	{
		PROFILE_SCOPE("setInitAccel");
		solve(size, mi_csrMatrix, mi_RHS, mi_LHS);
	}

	// --------------------------------------------------------------------------------------------
	// Now register the initial acceleration to nodes
	std::array<double, nDim> mi_accNd;

	for (auto& node : this->mv_meshNode) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node));
		auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node);

		for (int j = 0; j < nDim; ++j) {
			mi_accNd.at(j) = mi_LHS.at(pNode->mv_dofList.at(j)) * (this->mv_BCindex.at(pNode->mv_dofList.at(j)));
		}
		pNode->setPrevAcc(mi_accNd);
	}

	// Then register the initial acceleration to points
	std::array<double, nDim> mi_accPt;

	for (size_t i = 0; i < this->mv_meshPoint.size(); ++i) {

		// Fastest way to set a vector to zero
		memset(&mi_accPt[0], 0., nDim * sizeof(double));

		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(this->mv_meshPoint.at(i)));
		auto pPoint = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(this->mv_meshPoint.at(i));

		// Value of shape functions
		std::vector<double> mi_shape = pPoint->getElement()->getShapeFcOnPoint(pPoint->getXsi().data());

		int k = 0;
		for (auto& node : pPoint->getMeshElement()->mv_conect) {
			// Shades previous declaration of mi_accNd
			assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node));
			auto mi_accNd = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node)->getPrevAcc();

			// Evaluate contribution following incidence
			for (int i = 0; i < nDim; i++) {
				mi_accPt[i] += mi_accNd[i] * mi_shape.at(k);
			}
			k++;
		}
		pPoint->setPrevAcc(mi_accPt);
	}
}


// ================================================================================================
//
// Implementation of Mesh Member Function: ElemIFO
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::ElemIFO(O2P2::Proc::Comp::MeshElemVI* elem, std::vector<double>& RHS)
{
	O2P2::Geom::Material* pMat = elem->getMaterial();

	double Rho = pMat->getDensity();
	double Cm = pMat->getDamping();

	elem->mv_elHes.clear();
	elem->addMassContrib(Rho);

	// Set element forces to zero and then evaluate inertia contribution
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	for (int i = 0; i < elem->mv_conect.size(); ++i) {
		for (int j = 0; j < elem->mv_conect.size(); ++j) {
			assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(elem->mv_conect.at(j)));
			auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(elem->mv_conect.at(j));

			int nDof = pNode->getNumDOF();
			double mass = elem->mv_elHes(i * nDof, j * nDof);

			const auto vt = pNode->getCurVel();
			const auto at = pNode->getCurAcc();

			for (int k = 0; k < nDof; ++k) {
				elem->mv_elFor.at(i * nDof + k) -= mass * at[k] + mass * Cm * vt[k];
			}
		}
	}

	// Add element contribution to the right hand side vector
	this->mv_HesMutex.lock();
	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		RHS.at(mi_gdof) += elem->mv_elFor.at(i);
	}
	this->mv_HesMutex.unlock();
}


// ================================================================================================
//
// Implementation of Mesh Member Function: InclIFO
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::InclIFO(O2P2::Proc::Comp::MeshElemVI* elem, std::vector<double>& RHS)
{
	O2P2::Geom::Material* pMat = elem->getMaterial();

	double Rho = pMat->getDensity();
	double Cm = pMat->getDamping();

	elem->mv_elHes.clear();
	elem->addMassContrib(Rho);

	// Set element forces to zero and then evaluate inertia contribution
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	for (int i = 0; i < elem->mv_conect.size(); ++i) {
		for (int j = 0; j < elem->mv_conect.size(); ++j) {
			assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(elem->mv_conect.at(j)));
			auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(elem->mv_conect.at(j));

			int nDof = pNode->getNumDOF();
			double mass = elem->mv_elHes(i * nDof, j * nDof);

			const auto vt = pNode->getCurVel();
			const auto at = pNode->getCurAcc();

			for (int k = 0; k < nDof; ++k) {
				elem->mv_elFor.at(i * nDof + k) -= mass * at[k] + mass * Cm * vt[k];
			}
		}
	}

	M2S2::MatrixX mi_spread = elem->getSpreadMatrix();
	std::vector<double> expFor = mi_spread * elem->mv_elFor;

	// Add element contribution to the right hand side vector
	this->mv_HesMutex.lock();

	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		RHS.at(mi_gdof) -= expFor.at(i);
	}
	this->mv_HesMutex.unlock();
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): ElemSOE
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::ElemSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	// Fastest way to set a vector to zero
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	// The rest, is up to the library
	elem->mv_elHes.clear();
	elem->getContribution();

	// Add mass contribution
	O2P2::Geom::Material* pMat = elem->getMaterial();

	double Rho = pMat->getDensity();
	double Cm = pMat->getDamping();
	double mult = Rho * (1 + this->mv_gamma * Cm * dt) / (this->mv_beta * dt * dt);

	elem->addMassContrib(mult);

	// Impose restrictions directly in the element (faster)
	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);

		if (!this->mv_BCindex.at(mi_gdof)) {
			for (int j = 0; j < elem->mv_dofIndex.size(); ++j) {
				elem->mv_elHes.at(i, j) = 0.;
			}

			elem->mv_elHes.at(i, i) = 1.;
			elem->mv_elFor.at(i) *= this->mv_BCindex.at(mi_gdof);
		}
	}

	// Add element contribution to global system
	this->mv_HesMutex.lock();
	Hessian.push(elem->mv_elHes, elem->mv_dofIndex);

	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		RHS.at(mi_gdof) -= elem->mv_elFor.at(i);
	}
	this->mv_HesMutex.unlock();
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): InclSOE
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::InclSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	// Fastest way to set a vector to zero
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	// The rest, is up to the library
	elem->mv_elHes.clear();
	elem->getContribution();

	// Add mass contribution
	O2P2::Geom::Material* pMat = elem->getMaterial();

	double Rho = pMat->getDensity();
	double Cm = pMat->getDamping();
	double mult = Rho * (1 + this->mv_gamma * Cm * dt) / (this->mv_beta * dt * dt);

	elem->addMassContrib(mult);

	M2S2::MatrixX mi_spread = elem->getSpreadMatrix();

	// exp stands for expanded
	std::vector<double> expFor = mi_spread * elem->mv_elFor;
	M2S2::MatrixS expHes = (mi_spread * elem->mv_elHes * mi_spread.transpose()).SaveAsMatrixS();

	// Impose restrictions directly in every element (faster)
	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		if (!this->mv_BCindex.at(mi_gdof)) {
			for (int j = 0; j < elem->mv_dofIndex.size(); ++j) {
				expHes.at(i, j) = 0.;
			}

			expHes.at(i, i) = 1.;
			expFor.at(i) *= this->mv_BCindex.at(mi_gdof);
		}
	}

	// Add element contribution to global system
	this->mv_HesMutex.lock();
	Hessian.push(expHes, elem->mv_dofIndex);

	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		RHS.at(mi_gdof) -= expFor.at(i);
	}
	this->mv_HesMutex.unlock();
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): ElemIAC
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::ElemIAC(O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Mass, std::vector<double>& RHS)
{
	// Fastest way to set a vector to zero
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	// Get element contribution for the initial force
	elem->mv_elHes.clear();
	elem->getContribution();

	// Although the local hessian matrix is evaluated, is not required
	// Only the internal force is needed (saved in elem->m_elFor)
	elem->mv_elHes.clear();

	// Material properties
	O2P2::Geom::Material* pMat = elem->getMaterial();
	double Rho = pMat->getDensity();
	double Cm = pMat->getDamping();

	// Add mass contribution
	elem->addMassContrib(Rho);

	// Evaluate initial (known) velocity field
	std::vector<double> va(elem->mv_nDof);

	int i = 0;
	for (auto& node : elem->mv_conect) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node));
		auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshNode_MD<nDim>>(node);

		for (int j = 0; j < nDim; ++j) {
			va.at(i + j) = pNode->getPrevVel()[j];
		}
		i++;
	}

	// Inertial forces
	elem->mv_elFor = Cm * elem->mv_elHes * va;

	// Add element contribution to global system
	this->mv_HesMutex.lock();
	Mass.push(elem->mv_elHes, elem->mv_dofIndex);

	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		RHS.at(mi_gdof) -= elem->mv_elFor.at(i);
	}
	this->mv_HesMutex.unlock();
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): InclIAC
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MD<nDim>::InclIAC(O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Mass, std::vector<double>& RHS)
{
	// Fastest way to set a vector to zero
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	// Get element contribution for the initial force
	elem->mv_elHes.clear();
	elem->getContribution();

	// Although the local hessian matrix is evaluated, is not required
	// Only the internal force is needed (saved in elem->m_elFor)
	elem->mv_elHes.clear();

	// Material properties
	O2P2::Geom::Material* pMat = elem->getMaterial();
	double Rho = pMat->getDensity();
	double Cm = pMat->getDamping();

	// Add mass contribution
	elem->addMassContrib(Rho);

	// Evaluate initial (known) velocity field
	std::vector<double> va(elem->mv_nDof);

	int i = 0;
	for (auto& node : elem->mv_conect) {
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(node));
		auto pNode = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MD<nDim>>(node);

		for (int j = 0; j < nDim; ++j) {
			va.at(i + j) = pNode->getPrevVel()[j];
		}
		i++;
	}

	// Inertial forces
	elem->mv_elFor = Cm * elem->mv_elHes * va;

	// Expand inclusion contribution to dof's
	M2S2::MatrixX mi_spread = elem->getSpreadMatrix();

	// exp stands for expanded
	std::vector<double> expFor = mi_spread * elem->mv_elFor;
	M2S2::MatrixS expMas = (mi_spread * elem->mv_elHes * mi_spread.transpose()).SaveAsMatrixS();

	// Add element contribution to global system
	this->mv_HesMutex.lock();
	Mass.push(expMas, elem->mv_dofIndex);

	for (int i = 0; i < elem->mv_dofIndex.size(); ++i) {
		int mi_gdof = elem->mv_dofIndex.at(i);
		RHS.at(mi_gdof) -= expFor.at(i);
	}
	this->mv_HesMutex.unlock();
}
