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
#include "Mesh_MQS.h"
#include "NonLinearSolver.h"

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Mesh_MQS<2>;
template class O2P2::Proc::Mesh_MQS<3>;

template void O2P2::Proc::Mesh_MQS<2>::domainNodesToMesh(O2P2::Geom::Domain<2>* theDomain);
template void O2P2::Proc::Mesh_MQS<3>::domainNodesToMesh(O2P2::Geom::Domain<3>* theDomain);

// ================================================================================================
//
// Implementation of Protected Template Member Function: DomainNodesToMesh
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MQS<nDim>::domainNodesToMesh(O2P2::Geom::Domain<nDim>* theDomain) {

	// Generates a mesh node from each domain node. Do the same for points (but do not increase the number of DOF)
	for (std::shared_ptr<O2P2::Geom::Node>& node : theDomain->getNode()) {
		size_t dof = node->mv_index * nDim;
		this->mv_meshNode.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshNode_MQS<nDim>>(dof, node->getInitPos()));
		addDof(this->mv_meshNode.back()->getNumDOF());
	}

	for (std::shared_ptr<O2P2::Geom::Node>& node : theDomain->getPoint()) {
		this->mv_meshPoint.emplace_back(std::make_shared<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(0, node->getInitPos()));
	}
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): DomainElemToMesh
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MQS<nDim>::domainElemToMesh(O2P2::Geom::Domain<nDim>* theDomain) {

	// Matrix elements
	for (std::shared_ptr<O2P2::Geom::Elem::Element>& elem : theDomain->getElem()) {
		std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

		for (int i = 0; i < elem->getNumNodes(); ++i) {
			conect[i] = this->mv_meshNode[elem->getConectivity(i)->mv_index];
		}
		this->mv_meshElem.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshElem>(elem, conect));
		this->mv_meshElem.back()->setMaterialPoint();
		this->mv_meshElem.back()->setIndexing();
	}

	// Linear elements (Trusses and Rods)
	for (std::shared_ptr<O2P2::Geom::Elem::Element>& elem : theDomain->getBars()) {
		std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

		for (int i = 0; i < elem->getNumNodes(); ++i) {
			conect[i] = this->mv_meshNode[elem->getConectivity(i)->mv_index];
		}
		this->mv_meshBars.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshTruss<nDim>>(elem, conect));
		this->mv_meshBars.back()->setMaterialPoint();
		this->mv_meshBars.back()->setIndexing();
	}

	// Points should be evaluated after the MeshElements are created
	if (!this->evaluatePoints()) throw std::invalid_argument("\n\n\nPoint with undefined dimensionless coordinates\nCheck log file\n\n\n");

	// Particle elements
	for (std::shared_ptr<O2P2::Geom::Elem::Element>& elem : theDomain->getIncl()) {
		std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

		for (int i = 0; i < elem->getNumNodes(); ++i) {
			conect[i] = this->mv_meshPoint[elem->getConectivity(i)->mv_index];
		}
		this->mv_meshIncl.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshIncl<nDim>>(elem, conect));
		this->mv_meshIncl.back()->setMaterialPoint();
		this->mv_meshIncl.back()->setIndexing();
	}

	// Fiber elements
	for (std::shared_ptr<O2P2::Geom::Elem::Element>& elem : theDomain->getFibs()) {
		std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

		for (int i = 0; i < elem->getNumNodes(); ++i) {
			conect[i] = this->mv_meshPoint[elem->getConectivity(i)->mv_index];
		}
		this->mv_meshFibs.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshFibers<nDim>>(elem, conect));
		this->mv_meshFibs.back()->setMaterialPoint();
		this->mv_meshFibs.back()->setIndexing();
	}

	// Active face elements
	for (std::shared_ptr<O2P2::Geom::Elem::Element>& elem : theDomain->getFaces()) {
		std::vector<std::shared_ptr<O2P2::Proc::Comp::MeshNode>> conect(elem->getNumNodes());

		for (int i = 0; i < elem->getNumNodes(); ++i) {
			conect[i] = this->mv_meshNode[elem->getConectivity(i)->mv_index];
		}
		this->mv_meshFace.emplace_back(std::make_unique<O2P2::Proc::Comp::MeshFace>(elem, conect));
		this->mv_meshFace.back()->setMaterialPoint();
		this->mv_meshFace.back()->setIndexing();
	}
}


// ================================================================================================
//
// Implementation of Mesh Member Function: imposeNeumannBC
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MQS<nDim>::imposeNeumannBC(std::vector<double>& RHS)
{
	LOG("Mesh.imposeNeumannBC: Initiating current load vector");

	double mi_curTime = mv_curTimeStep * mv_LoadStep.at(mv_curLoadStep)->mv_timeStep;

	// If a Dirichlet boundary vector is imposed (like a prescribed displacement), it is first imposed to the RHS 
	for (auto& DBC : mv_LoadStep.at(mv_curLoadStep)->mv_DirichletBC) {
		RHS.at(DBC->mv_Dof) = DBC->mv_value * (DBC->mv_timeVar[0] + DBC->mv_timeVar[1] * mi_curTime + DBC->mv_timeVar[2] * mi_curTime * mi_curTime +
			DBC->mv_timeVar[3] * sin(DBC->mv_timeVar[4] * mi_curTime) + DBC->mv_timeVar[5] * cos(DBC->mv_timeVar[6] * mi_curTime));
	}

	// Boundary conditions are saved directly on DOF number
	for (auto& NBC : mv_LoadStep.at(mv_curLoadStep)->mv_NeumannBC) {
		RHS.at(NBC->mv_Dof) = NBC->mv_value * (NBC->mv_timeVar[0] + NBC->mv_timeVar[1] * mi_curTime + NBC->mv_timeVar[2] * mi_curTime * mi_curTime +
			NBC->mv_timeVar[3] * sin(NBC->mv_timeVar[4] * mi_curTime) + NBC->mv_timeVar[5] * cos(NBC->mv_timeVar[6] * mi_curTime));
	}

	// Even though there are BC imposed, the load vector must be set ot zero
	//for (size_t i = 0; i < RHS.size(); ++i) {
	//	RHS.at(i) *= this->m_BCIndex[i];
	//}
}


// ================================================================================================
//
// Implementation of Private Template Member Function (2D and 3D): evaluatePoints
//
// ================================================================================================
template<int nDim>
bool O2P2::Proc::Mesh_MQS<nDim>::evaluatePoints()
{
	PROFILE_FUNCTION();

	if(!this->mv_meshPoint.empty()) LOG("Mesh.evaluatePoints: Stablishing dimensionless coordinates for immersed points");

	// If no point exist, return true
	bool all_feasible = true;
#ifdef _DEBUG
	int mi_curPoint = 0;
#endif

	// Every Point is attached to one (and only one) element
	// A Nonlinear system of equation is required to find this out
	O2P2::Proc::NLS_NewtonRaphson theNLSolver(1, 50, 1.e-6);

	for (auto& mi_point : this->mv_meshPoint) {
		// Downcast the point to MeshPoint_MQS
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(mi_point));
		auto point = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(mi_point);

		bool mi_feasible = false;
#ifdef _DEBUG
		int mi_elIndex = 0;
#endif

		for (auto& mi_elem : mv_meshElem) {
			auto& elem = mi_elem->mv_ptElem;

			// Check the distance between element center and current point.
			double mi_distance = 0.;
			for (int i = 0; i < nDim; i++) {
				mi_distance += std::pow(point->getCurPos()[i] - *(elem->getCentroid() + i), 2);
			}
			mi_distance = std::sqrt(mi_distance);

			// If this distance is less than element minimum bounding circle radius, try to find feasible xsi.
			if (mi_distance < 1.05 * elem->getRadius()) {

				// Set current element as trial
				point->setMeshElement(mi_elem.get());

				std::vector<double> mi_trialXsi(nDim);

				// Initial xsi trial
				for (int i = 0; i < nDim; i++) {
					if (*(elem->getCentroid() + i) > 0.) {
						mi_trialXsi.at(i) = point->getCurPos()[i] / *(elem->getCentroid() + i);
					}
					else {
						mi_trialXsi.at(i) = 0.;
					}

				}
				point->setXsi(mi_trialXsi);

				// If feasible is true, values for xsi were found
				mi_feasible = theNLSolver.runNLS(point.get());
			}

			// Ask the current element to verify trial xsi
			if (mi_feasible) {
				mi_feasible = elem->evaluateXsi(point->getXsi().data());
				if (mi_feasible) break;
			}
#ifdef _DEBUG
			mi_elIndex++;
#endif
		}
		if (!mi_feasible) {
			LOG("Mesh.evaluatePoints: Unable to find Xsi for point " << std::to_string(mi_curPoint));
			all_feasible = mi_feasible;
		}
		else {
			LOG("Mesh.evaluatePoints: Point " << std::to_string(mi_curPoint) << " is immersed on element " << std::to_string(mi_elIndex));
		}
#ifdef _DEBUG
		mi_curPoint++;
#endif
	}

	return all_feasible;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): assembleSOE
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MQS<nDim>::assembleSOE(M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	PROFILE_FUNCTION();

	// Set Hessian to zero
	Hessian.clear();

	// imposeInertia will only have effect in Mesh_MD
	this->imposeInertiaBC(RHS);
	double dt = this->mv_LoadStep[this->mv_curLoadStep]->mv_timeStep;

	// --------------------------------------------------------------------------------------------
	// Fourth nested loop - Matrix Elements
	std::for_each(std::execution::par_unseq, this->mv_meshElem.begin(), this->mv_meshElem.end(), [&](auto& elem)
		{ this->ElemSOE(dt, elem.get(), Hessian, RHS); } );

	// --------------------------------------------------------------------------------------------
	// Another fourth nested loop - Bars Elements (Matrix)
	std::for_each(std::execution::par_unseq, this->mv_meshBars.begin(), this->mv_meshBars.end(), [&](auto& elem)
		{ this->ElemSOE(dt, elem.get(), Hessian, RHS); } );

	// --------------------------------------------------------------------------------------------
	// Yet fourth nested loop - Particle Elements (Inclusions)
	std::for_each(std::execution::par_unseq, this->mv_meshIncl.begin(), this->mv_meshIncl.end(), [&](auto& elem)
		{ this->InclSOE(dt, elem.get(), Hessian, RHS); } );

	// --------------------------------------------------------------------------------------------
	// More of the fourth nested loop - Fibers Elements (Inclusions)
	std::for_each(std::execution::par_unseq, this->mv_meshFibs.begin(), this->mv_meshFibs.end(), [&](auto& elem)
		{ this->InclSOE(dt, elem.get(), Hessian, RHS); } );

	// --------------------------------------------------------------------------------------------
	// Even more of the fourth nested loop - Active Face Elements
	std::for_each(std::execution::par_unseq, this->mv_meshFace.begin(), this->mv_meshFace.end(), [&](auto& elem)
		{ this->ElemSOE(dt, elem.get(), Hessian, RHS); });
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): ElemSOE
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MQS<nDim>::ElemSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	// Fastest way to set a vector to zero
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	// The rest, is up to the library
	elem->mv_elHes.clear();
	elem->getContribution();

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
void O2P2::Proc::Mesh_MQS<nDim>::InclSOE(double dt, O2P2::Proc::Comp::MeshElemVI* elem, M2S2::sparseMatrix& Hessian, std::vector<double>& RHS)
{
	// Fastest way to set a vector to zero
	memset(&elem->mv_elFor[0], 0., elem->mv_elFor.size() * sizeof(double));

	// The rest, is up to the library
	elem->mv_elHes.clear();
	elem->getContribution();

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
// Implementation of Template Member Function (2D and 3D): setTrial
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Mesh_MQS<nDim>::setTrial(std::vector<double>& LHS)
{
	//LOG("Mesh_MQS.setTrial: Updating trial solution to nodes");
	double* sol = new double[nDim];

	for (int k = 0; k < LHS.size(); ++k) {
		LHS.at(k) *= this->mv_BCindex.at(k);
	}

	for (auto& node : mv_meshNode) {
		for (int j = 0; j < nDim; ++j) {
			sol[j] = LHS.at(node->mv_dofList.at(j));
		}
		node->updateTrial(sol);
	}

	//LOG("Mesh_MQS.setTrial: Updating trial solution to points");
	for (auto& point : mv_meshPoint) {

		// Fastest way to set a vector to zero
		memset(&sol[0], 0., nDim * sizeof(double));

		// Downcast do MeshPoint
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(point));
		auto pPoint = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(point);

		// Value of shape functions
		std::vector<double> mi_shape = pPoint->getElement()->getShapeFcOnPoint(pPoint->getXsi().data());

		int j = 0;
		for (auto& node : pPoint->getMeshElement()->mv_conect) {
			// Trial position of element node
			const auto Pos = node->getTrialPos();

			// Evaluate contribution following incidence
			for (int i = 0; i < nDim; i++) {
				sol[i] += Pos[i] * mi_shape.at(j);
			}
			j++;
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
void O2P2::Proc::Mesh_MQS<nDim>::setCommit()
{
	LOG("Mesh_MQS.setCommit: Committing solution to nodes and saving for post-process");

	// Commit the solution and transfer it to mi_sol (for post-processing).
	// This was made only for displacement, nothing else.

	// The number of DOF of Points is still associated to the nDim.
	std::vector<double> mi_sol(mv_TotalDof + this->mv_meshPoint.size() * nDim);

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

	for (auto& point : mv_meshPoint) {
		// Set trial position as current (commit position)
		point->setCurrent();

		// Save current position to Sol vector for post-processing
		auto y = point->getCurPos();
		for (size_t j = 0; j < point->getNumDOF(); ++j) {
			mi_sol.at(i) = *(y + j);
			i++;
		}
	}

	if ((mv_curLoadStep == 0) && double(mv_curLoadStep / this->m_PostPt->getOutputFrequency()) - int(mv_curLoadStep / this->m_PostPt->getOutputFrequency()) < 1.e-6) {
		// Send the information to the post-process container.
		this->m_PostPt->addSolution(this->mv_currentTime, mi_sol);
	}
}
