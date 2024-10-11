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
#include "Mesh Components\MeshBars.h"

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Comp::MeshTruss<2>;
template class O2P2::Proc::Comp::MeshTruss<3>;

template void O2P2::Proc::Comp::MeshTruss<2>::setMaterialPoint();
template void O2P2::Proc::Comp::MeshTruss<3>::setMaterialPoint();

template void O2P2::Proc::Comp::MeshTruss<2>::setIndexing();
template void O2P2::Proc::Comp::MeshTruss<3>::setIndexing();

template void O2P2::Proc::Comp::MeshTruss<2>::getContribution_SVK(double section);
template void O2P2::Proc::Comp::MeshTruss<3>::getContribution_SVK(double section);

template void O2P2::Proc::Comp::MeshTruss<2>::addMassContrib(const double& mult);
template void O2P2::Proc::Comp::MeshTruss<3>::addMassContrib(const double& mult);

// ================================================================================================
//
// Implementation of Template Member Function: setMaterialPoint
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshTruss<nDim>::setMaterialPoint()
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_ptElem->getShapeDerivative();

	// Number of integration points
	int mi_nIP = mv_ptElem->getNumIP();

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {

		// Once created, evaluates F0 (Ref. Jacobian)
		// For fiber elements, F0 holds the initial length in each direction
		std::array<double, nDim> mi_F0;
		mi_F0.fill(0.);

		// Some pointer arithmetic is required
		int n = i * mv_ptElem->getNumNodes();

		for (int j = 0; j < nDim; ++j) {
			for (int m = 0; m < mv_ptElem->getNumNodes(); ++m) {
				auto mi_pNode = mv_ptElem->getConectivity(m);

				mi_F0[j] += mi_pNode->getInitPos()[j] * *(mi_pShape + n + m);
			}
		}

		// Record the Jacobian Matrix on every Material (Integration) point
		this->mv_matPoint.emplace_back(std::make_unique<O2P2::Proc::Comp::MaterialPoint>(nDim, mi_F0.data()));
	}
};


// ================================================================================================
//
// Implementation of MeshTruss Member Function: setIndexing
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshTruss<nDim>::setIndexing()
{
	for (int i = 0; i < mv_conect.size(); ++i) {
		for (int j = 0; j < nDim; ++j) {
			mv_dofIndex.push_back(mv_conect.at(i)->mv_dofList.at(j));
		}
	}
}

// ================================================================================================
//
// Implementation of MeshTruss Member Function: getContribution
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshTruss<nDim>::getContribution()
{
	auto mi_pMat = mv_ptElem->getMaterial();

	// Access cross section through linear element
	assert(std::dynamic_pointer_cast<O2P2::Geom::Elem::ElemLinear<nDim>>(mv_ptElem));
	double mi_section = std::static_pointer_cast<O2P2::Geom::Elem::ElemLinear<nDim>>(mv_ptElem)->getSection()->getSection();

	switch (mi_pMat->getMaterialType())
	{
	case MaterialType::SVK_ISO:
	{
		this->getContribution_SVK(mi_section);
		break;
	}
	default: { throw std::invalid_argument("\n\n\nMeshTruss::getContribution - The material type was not defined.\n\n\n"); break; }
	}
}


// ================================================================================================
//
// Implementation of MeshTruss Member Function: getContribution_SVK
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshTruss<nDim>::getContribution_SVK(double section)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_ptElem->getShapeDerivative();
	auto mi_pWeight = mv_ptElem->getWeight();

	auto mi_nNodes = mv_ptElem->getNumNodes();
	auto mi_nIP = mv_ptElem->getNumIP();

	auto mi_pMat = dynamic_cast<O2P2::Geom::Mat_SVK_ISO*>(mv_ptElem->getMaterial());

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {

		// Weight for Gauss point and Jacobian
		double numInt = *(mi_pWeight + i) * std::sqrt(this->mv_matPoint.at(i)->getJacobian()) * section;

		// For truss and fiber elements, F0 and F1 holds the initial and current length in each direction
		std::array<double, nDim> F1;
		F1.fill(0.);

		// Some pointer arithmetic is required
		int n = i * mi_nNodes;

		for (int j = 0; j < nDim; ++j) {
			int m = 0;
			for (auto& node : mv_conect) {
				F1.at(j) += node->getTrialPos()[j] * *(mi_pShape + n + m);
				m++;
			}
		}

		double dA0 = this->mv_matPoint.at(i)->getJacobian();;
		double dA1 = 0.;

		for (int j = 0; j < nDim; ++j) {
			dA1 += F1.at(j) * F1.at(j);
		}

		// Green-Lagrange Strain
		double L = 0.5 * (dA1 - dA0) / dA0;

		// Young modulus
		double E = mi_pMat->getLongitudinalModulus();

		// Second Piola-Kirchhoff Stress
		double S = E * L;

		// Green-Lagrange Strain Derivative
		std::array<double, nDim> dLdy;
		for (int j = 0; j < nDim; ++j) {
			dLdy.at(j) = F1.at(j) / dA0;
		}

		// Green-Lagrange Strain Second Derivative
		// It is the same in every direction, and equal to dA0
		double d2Ldy2 = dA0;

		// Specific energy derivative
		std::array<double, nDim> dUdy;
		for (int j = 0; j < nDim; ++j) {
			dUdy.at(j) = S * dLdy.at(j);
		}

		// Auxiliary for Hessian matrix
		M2S2::MatrixS mi_Ha(nDim);

		for (int j = 0; j < nDim; ++j) {
			for (int k = j; k < nDim; ++k) {
				mi_Ha(j, k) = E * dLdy.at(j) * dLdy.at(k);
			}
			mi_Ha(j, j) += S * d2Ldy2;
		}

		// Integration of the specific energy derivative = internal force
		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = 0; k < nDim; ++k) {
				mv_elFor.at(j * nDim + k) += dUdy.at(k) * *(mi_pShape + n + j) * numInt;
			}
		}

		// Elemental Hessian matrix
		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = 0; k < nDim; ++k) {
				for (int o = k; o < nDim; ++o) {
					mv_elHes.at(j * nDim + k, j * nDim + o) += *(mi_pShape + n + j) * mi_Ha(k, o) * *(mi_pShape + n + j) * numInt;
				}
				for (int m = j + 1; m < mi_nNodes; ++m) {
					for (int o = 0; o < nDim; ++o) {
						mv_elHes.at(j * nDim + k, m * nDim + o) += *(mi_pShape + n + j) * mi_Ha(k, o) * *(mi_pShape + n + m) * numInt;
					}
				}
			}
		}
	}
}


// ================================================================================================
//
// Implementation of MeshTruss Member Function: addMassContrib
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshTruss<nDim>::addMassContrib(const double& mult)
{
	auto mi_pShape = mv_ptElem->getShapeFc();
	auto mi_pWeight = mv_ptElem->getWeight();

	auto mi_nNodes = mv_ptElem->getNumNodes();
	auto mi_nIP = mv_ptElem->getNumIP();

	assert(std::dynamic_pointer_cast<O2P2::Geom::Elem::ElemLinear<nDim>>(mv_ptElem));
	double mi_section = std::static_pointer_cast<O2P2::Geom::Elem::ElemLinear<nDim>>(mv_ptElem)->getSection()->getSection();

	for (int i = 0; i < mi_nIP; ++i) {
		int p = i * mi_nNodes;
		double numInt = mult * *(mi_pWeight + i) * mv_matPoint[i]->getJacobian() * mi_section;

		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = j; k < mi_nNodes; ++k) {
				for (int m = 0; m < nDim; ++m) {
					mv_elHes.at(j * nDim + m, k * nDim + m) += *(mi_pShape + p + j) * *(mi_pShape + p + k) * numInt;
				}
			}
		}
	}
}
