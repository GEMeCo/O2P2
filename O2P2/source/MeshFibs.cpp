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
#include "Mesh Components\MeshFibs.h"

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Comp::MeshFibers<2>;
template class O2P2::Proc::Comp::MeshFibers<3>;

template void O2P2::Proc::Comp::MeshFibers<2>::setIndexing();
template void O2P2::Proc::Comp::MeshFibers<3>::setIndexing();

template M2S2::MatrixX O2P2::Proc::Comp::MeshFibers<2>::getSpreadMatrix();
template M2S2::MatrixX O2P2::Proc::Comp::MeshFibers<3>::getSpreadMatrix();


// ================================================================================================
//
// Implementation of MeshFibers Member Function: setIndexing
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshFibers<nDim>::setIndexing()
{
	// Calculate the number of DOF for expanded indexing
	int mi_size = 0;
	for (auto& node : this->mv_conect) {
		// Downcasting to inclusion node (point)
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(node));
		auto mi_point = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(node);
		mi_size += mi_point->getElement()->getNumNodes() * mi_point->getElement()->getNumNodalDof();
	}
	this->mv_dofIndex.clear();
	this->mv_dofIndex.reserve(mi_size);

	// Let's assemble the inclusion dof
	for (auto& node : this->mv_conect) {
		// Downcasting to inclusion node (point)
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(node));
		auto mi_point = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(node);

		// Copy dof indexing
		std::copy(mi_point->getMeshElement()->mv_dofIndex.begin(), mi_point->getMeshElement()->mv_dofIndex.end(), std::back_inserter(this->mv_dofIndex));
	}
}


// ================================================================================================
//
// Implementation of MeshFibers Member Function: getSpreadMatrix
//
// ================================================================================================
template<int nDim>
M2S2::MatrixX O2P2::Proc::Comp::MeshFibers<nDim>::getSpreadMatrix()
{
	M2S2::MatrixX mi_spread(this->mv_dofIndex.size(), this->mv_nDof);

	int m = 0;
	// Point indexing
	for (auto& node : this->mv_conect) {
		// Downcasting to inclusion node (point)
		assert(std::dynamic_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(node));
		auto mi_point = std::static_pointer_cast<O2P2::Proc::Comp::MeshPoint_MQS<nDim>>(node);

		// Element that the inclusion node points to
		auto mi_ptElem = mi_point->getElement();

		// Element information
		int mi_elDof = mi_ptElem->getNumNodalDof();
		int mi_nNodes = mi_ptElem->getNumNodes();

		// Element shape functions in current coordinates
		std::vector<double> mi_shape = mi_ptElem->getShapeFcOnPoint(mi_point->getXsi().data());

		// Loop over the number of degrees of freedom of current element (incidence)
		for (int i = 0; i < mi_nNodes; i++) {
			int k = m * mi_nNodes + i * mi_elDof;
			for (int j = 0; j < mi_elDof; j++) {
				mi_spread(k + j, m + j) += mi_shape.at(i);
			}
		}
		m = m + mi_elDof;
	}
	return mi_spread;
}
