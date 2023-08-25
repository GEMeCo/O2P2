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
#include "MeshElem.h"

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Comp::MeshElem_SVK<2>;
template class O2P2::Proc::Comp::MeshElem_SVK<3>;

template void O2P2::Proc::Comp::MeshElem_SVK<2>::setMaterialPoint();
template void O2P2::Proc::Comp::MeshElem_SVK<3>::setMaterialPoint();

template void O2P2::Proc::Comp::MeshElem_SVK<2>::setIndexing();
template void O2P2::Proc::Comp::MeshElem_SVK<3>::setIndexing();

// ================================================================================================
//
// Implementation of Template Member Function: setMaterialPoint
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Comp::MeshElem_SVK<nDim>::setMaterialPoint()
{
	int mi_nIP = mv_pElem->getNumIP();
	int mi_nElDim = mv_pElem->getDIM();

	// The BaseElement pointer does not have access to getConectivity
	// thus, must cast it to a derived class
	O2P2::Prep::Elem::Element<nDim>* mi_pElem = static_cast<O2P2::Prep::Elem::Element<nDim>*>(mv_pElem.get());

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {
		// Once created, evaluates F0 (Ref. Jacobian)
		M2S2::Dyadic2N mi_F0(nDim);

		// Some pointer arithmetic is required
		int n = i * mv_pElem->getNumNodes() * nDim;

		for (int j = 0; j < nDim; ++j) {
			for (int m = 0; m < mv_pElem->getNumNodes(); ++m) {
				auto mi_pNode = mi_pElem->getConectivity(m);

				for (int k = 0; k < mi_nElDim; ++k) {
					mi_F0(j, k) += mi_pNode->getInitPos()[j] * *(mv_pElem->getShapeDerivative() + n + m * mi_nElDim + k);
				}
			}
		}

		// Record the Jacobian Matrix on every Material (Integration) point
		this->mv_matPoint.emplace_back(std::make_unique<O2P2::Proc::Comp::MaterialPoint>(mi_F0));
	}
}

// ================================================================================================
//
// Implementation of Template Member Function: setIndexing
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshElem_SVK<nDim>::setIndexing()
{
	for (int i = 0; i < mv_conect.size(); ++i) {
		for (int j = 0; j < nDim; ++j) {
			mv_elIndex.push_back(mv_conect.at(i)->mv_DofIndex + j);
		}
	}
}


// ================================================================================================
//
// Implementation of MeshElem_SVK Member Function: getYoungModulus
//
// ================================================================================================
template<int nDim> double O2P2::Proc::Comp::MeshElem_SVK<nDim>::getYoungModulus()
{
	assert((mv_pElem->getDIM() == 1) && "MeshElem_SVK<nDim>::getYoungModulus() should only be used with linear elements");

	// Downcasting to linear element.
	O2P2::Prep::Elem::ElementLinear<nDim>* mi_pElem = static_cast<O2P2::Prep::Elem::ElementLinear<nDim>*>(mv_pElem.get());

	// Pointer to section
	auto mi_pSec = mi_pElem->getSection();

	// Pointer to the material
	auto mi_pMat = mi_pElem->getMaterial();

	// Should also downcast to SVK material
	O2P2::Prep::Mat_SVK_ISO* mi_pSVK = static_cast<O2P2::Prep::Mat_SVK_ISO*>(mi_pMat);

	return mi_pSVK->getLongitudinalModulus() * mi_pSec->getSection();
}


// ================================================================================================
//
// Implementation of MeshElem_SVK Member Function: getContribution
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshElem_SVK<nDim>::getContribution()
{
	auto mi_pMat = mv_pElem->getMaterial();
	auto mi_nElDim = mv_pElem->getDIM();

	if (mi_nElDim == 1) {
		if (mi_pMat->getMaterialType() == MaterialType::SVK_ISO)
			this->getContribution_SVK<1>();
	}

	if (mi_nElDim == 2) {
		if (mi_pMat->getMaterialType() == MaterialType::SVK_ISO)
			this->getContribution_SVK<2>();
	}

	if (mi_nElDim == 3) {
		if (mi_pMat->getMaterialType() == MaterialType::SVK_ISO)
			this->getContribution_SVK<3>();
	}
}


// ================================================================================================
//
// Specialization of MeshElem_SVK Member Function for 1D elem in 2D environment: getContribution 
//
// ================================================================================================
template<> template<>
void O2P2::Proc::Comp::MeshElem_SVK<2>::getContribution_SVK<1>()
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_pElem->getShapeDerivative();
	auto mi_pWeight = mv_pElem->getWeight();

	auto mi_nNodes = mv_pElem->getNumNodes();
	auto mi_nIP = mv_pElem->getNumIP();

	const int nDim = 2;

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {
		
		// For fiber elements, F0 and F1 holds the initial and current length in each direction
		M2S2::Dyadic2N mi_F1(nDim);

		// Some pointer arithmetic is required
		int n = i * mi_nNodes * nDim;

		for (int j = 0; j < nDim; ++j) {
			int m = 0;

			for (auto& node : mv_conect) {
				mi_F1(j, 0) += node->getTrialPos()[j] * *(mi_pShape + n + m);
				m++;
			}
		}

		double dA0 = mv_matPoint.at(i)->getJacobian();
		double dA1 = 0.;

		for (int j = 0; j < nDim; ++j) {
			dA1 += mi_F1(j, 0) * mi_F1(j, 0);
		}
		
		// Green-Lagrange Strain
		double L = 0.5 * (dA1 - dA0) * dA0;

		// Elasticity properties
		double E = this->getYoungModulus();
		
		// Second Piola-Kirchhoff Stress
		double S = E * L;

		// Green-Lagrange Strain Derivative
		std::vector<double> dLdy(nDim);
		for (int j = 0; j < nDim; ++j) {
			dLdy.at(j) = mi_F1(j, 0) * mi_F1(j, 0) * dA0;
		}

		// Green-Lagrange Strain Second Derivative
		// It is the same in every direction, and equal to dA0
		double d2Ldy2 = dA0;

		// Specific energy derivative
		double dU2dy = S * d2Ldy2;

		// Weight for Gauss point and Jacobian
		double numInt = *(mi_pWeight + i) * this->mv_matPoint.at(i)->getJacobian();

		// Integration of the specific energy derivative = internal force
		for (int j = 0; j < mi_nNodes; ++j) {
			mv_elFor.at(j) += dU2dy * *(mi_pShape + n + j) * numInt;
		}

		// Auxiliary for Hessian matrix
		M2S2::MatrixS mi_Ha(nDim);

		for (int j = 0; j < nDim; ++j) {
			for (int k = j; k < nDim; ++k) {
				mi_Ha(j, k) = E * dLdy.at(j) * dLdy.at(k);
			}
			mi_Ha(j, j) += S * d2Ldy2;
		}

		// First part of Hessian
		M2S2::MatrixX H1(nDim, nDim * mi_nNodes);

		// Kronecker delta
		M2S2::Dyadic2S I = M2S2::Dyadic2S::identity(nDim);

		for (int j = 0; j < nDim; ++j) {
			for (int k = 0; k < nDim; ++k) {
				for (int l = 0; l < nDim; ++l) {
					for (int m = 0; m < nDim; ++m) {
						H1(l, j * nDim + k) += mi_Ha(l, m) * I(m, k) * *(mi_pShape + n + j);
					}
				}
			}
		}

		// Local Hessian
		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = 0; k < nDim; ++k) {

				for (int l = 0; l < mi_nNodes; ++l) {
					for (int m = 0; m < nDim; ++m) {

						for (int o = 0; o < nDim; ++o) {
							mv_elHes.at(j * mi_nNodes + k, l * mi_nNodes + m) += *(mi_pShape + n + j) * I(o, k) * H1(o, l * mi_nNodes + m) * numInt;
						}
					}
				}
			}
		}
	}
}


// ================================================================================================
//
// Specialization of MeshElem_SVK Member Function for 1D elem in 2D environment: getContribution 
//
// ================================================================================================
template<> template<>
void O2P2::Proc::Comp::MeshElem_SVK<3>::getContribution_SVK<1>()
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_pElem->getShapeDerivative();
	auto mi_pWeight = mv_pElem->getWeight();

	auto mi_nNodes = mv_pElem->getNumNodes();
	auto mi_nIP = mv_pElem->getNumIP();

	const int nDim = 2;

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {

		// For fiber elements, F0 and F1 holds the initial and current length in each direction
		M2S2::Dyadic2N mi_F1(nDim);

		// Some pointer arithmetic is required
		int n = i * mi_nNodes * nDim;

		for (int j = 0; j < nDim; ++j) {
			int m = 0;

			for (auto& node : mv_conect) {
				mi_F1(j, 0) += node->getTrialPos()[j] * *(mi_pShape + n + m);
				m++;
			}
		}

		double dA0 = mv_matPoint.at(i)->getJacobian();
		double dA1 = 0.;

		for (int j = 0; j < nDim; ++j) {
			dA1 += mi_F1(j, 0) * mi_F1(j, 0);
		}

		// Green-Lagrange Strain
		double L = 0.5 * (dA1 - dA0) * dA0;

		// Elasticity properties
		double E = this->getYoungModulus();

		// Second Piola-Kirchhoff Stress
		double S = E * L;

		// Green-Lagrange Strain Derivative
		std::vector<double> dLdy(nDim);
		for (int j = 0; j < nDim; ++j) {
			dLdy.at(j) = mi_F1(j, 0) * mi_F1(j, 0) * dA0;
		}

		// Green-Lagrange Strain Second Derivative
		// It is the same in every direction, and equal to dA0
		double d2Ldy2 = dA0;

		// Specific energy derivative
		double dU2dy = S * d2Ldy2;

		// Weight for Gauss point and Jacobian
		double numInt = *(mi_pWeight + i) * this->mv_matPoint.at(i)->getJacobian();

		// Integration of the specific energy derivative = internal force
		for (int j = 0; j < mi_nNodes; ++j) {
			mv_elFor.at(j) += dU2dy * *(mi_pShape + n + j) * numInt;
		}

		// Auxiliary for Hessian matrix
		M2S2::MatrixS mi_Ha(nDim);

		for (int j = 0; j < nDim; ++j) {
			for (int k = j; k < nDim; ++k) {
				mi_Ha(j, k) = E * dLdy.at(j) * dLdy.at(k);
			}
			mi_Ha(j, j) += S * d2Ldy2;
		}

		// First part of Hessian
		M2S2::MatrixX H1(nDim, nDim * mi_nNodes);

		// Kronecker delta
		M2S2::Dyadic2S I = M2S2::Dyadic2S::identity(nDim);

		for (int j = 0; j < nDim; ++j) {
			for (int k = 0; k < nDim; ++k) {
				for (int l = 0; l < nDim; ++l) {
					for (int m = 0; m < nDim; ++m) {
						H1(l, j * nDim + k) += mi_Ha(l, m) * I(m, k) * *(mi_pShape + n + j);
					}
				}
			}
		}

		// Local Hessian
		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = 0; k < nDim; ++k) {

				for (int l = 0; l < mi_nNodes; ++l) {
					for (int m = 0; m < nDim; ++m) {

						for (int o = 0; o < nDim; ++o) {
							mv_elHes.at(j * mi_nNodes + k, l * mi_nNodes + m) += *(mi_pShape + n + j) * I(o, k) * H1(o, l * mi_nNodes + m) * numInt;
						}
					}
				}
			}
		}
	}
}


// ================================================================================================
//
// Implementation of MeshElem_SVK Member Function: getContribution_SVK_ISO
//
// ================================================================================================
template<int nDim> template<int nElDim>
void O2P2::Proc::Comp::MeshElem_SVK<nDim>::getContribution_SVK()
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_pElem->getShapeDerivative();
	auto mi_pWeight = mv_pElem->getWeight();

	auto mi_nNodes = mv_pElem->getNumNodes();
	auto mi_nIP = mv_pElem->getNumIP();

	auto mi_pMat = mv_pElem->getMaterial();

	// Elasticity tensor in material description
	M2S2::Dyadic4C E(nDim);
	E = mi_pMat->getConstitutiveMatrix(nDim);

	if (nElDim == 2) {
		// Downcasting to plane element.
		O2P2::Prep::Elem::ElementPlane* mi_pElem = static_cast<O2P2::Prep::Elem::ElementPlane*>(mv_pElem.get());

		// Pointer to section
		auto mi_pSec = mi_pElem->getSection();

		// Thickness
		auto mi_tk = mi_pSec->getSection();

		E *= mi_tk;
	}

	const int mi_nDof = mi_nNodes * nDim;
	const int mi_nVoigt = 3 * nDim - 3;

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {

		// Deformation gradient from undeformed to non-dimensional space
		M2S2::Dyadic2N F0(nDim);
		F0 = mv_matPoint.at(i)->getJacobianMatrix();

		// Deformation gradient from deformed to non-dimensional space
		M2S2::Dyadic2N F1(nDim);

		// Some pointer arithmetic is required
		// p and numInt are related to current integration point and weighting numerical integration
		int p = i * mi_nDof;
		double numInt = *(mi_pWeight + i) * this->mv_matPoint.at(i)->getJacobian();

		for (int j = 0; j < nDim; ++j) {
			for (int m = 0; m < mi_nNodes; ++m) {
				auto node = mv_conect.at(m);

				for (int k = 0; k < nDim; ++k) {
					F1.at(j, k) += node->getTrialPos()[j] * *(mi_pShape + p + m * nDim + k);
				}
			}
		}

		// Deformation gradient
		M2S2::Dyadic2N F(nDim);

		// evaluate F = F1 * F0i (F0i -> inverse of jacobian matrix)
		F = F1 * F0;

		// Right Cauchy-Green (or Green deformation tensor)
		M2S2::Dyadic2S C(nDim);
		C = F.getATA();

		// Green-Lagrange Strain
		M2S2::Dyadic2S L(nDim);
		L = 0.5 * (C - M2S2::Dyadic2S::identity(nDim));

		// Second Piola-Kirchhoff Stress
		M2S2::Dyadic2S S(nDim);
		S = E * L;
		auto SM = S.getVoigtMnenomics();

		// Green-Lagrange Strain Derivative
		// Notice that dLdy is written in Voigt notation with mnemonic rule
		M2S2::MatrixX dLdy(mi_nVoigt, mi_nDof);

		for (int r = 0; r < mi_nNodes; ++r) {
			for (int d = 0; d < nDim; ++d) {
				int mi_curDof = r * nDim + d;

				for (int m = 0; m < nDim; ++m) {
					int n = m;

					for (int j = 0; j < nDim; ++j) {
						dLdy(m, mi_curDof) += 0.5 * (F(d, m) * F0(j, n) * *(mi_pShape + p + r * nDim + j) + F0(j, m) * F(d, n) * *(mi_pShape + p + r * nDim + j));
					}

					for (n = m + 1; n < nDim; ++n) {
						int l = mi_nVoigt - m - n;

						for (int j = 0; j < nDim; ++j) {
							dLdy(l, mi_curDof) += 0.5 * (F(d, m) * F0(j, n) * *(mi_pShape + p + r * nDim + j) + F0(j, m) * F(d, n) * *(mi_pShape + p + r * nDim + j));
						}
					}
				}

				// Specific energy derivative - it is a vector, but we are not going to use it again.
				double dUedy = 0;

				for (int j = 0; j < mi_nVoigt; ++j) {
					dUedy += SM.at(j) * dLdy(j, mi_curDof);
				}

				for (int j = nDim; j < mi_nVoigt; ++j) {
					dUedy += SM.at(j) * dLdy(j, mi_curDof);
				}

				// Integration of the specific energy derivative = internal force
				mv_elFor.at(mi_curDof) += dUedy * numInt;
			}
		}

		//Hessian += dLdy.transpose() * E * dLdy * numInt;
		for (int j = 0; j < mi_nDof; ++j) {
			for (int k = j; k < mi_nDof; ++k) {
				for (int l = 0; l < nDim; ++l) {
					for (int m = l; m < nDim; ++m) {
						mv_elHes.at(j, k) += dLdy.at(m, j) * E.at(m, l) * dLdy.at(l, k) * numInt;
					}
				}
				for (int l = nDim; l < mi_nVoigt; ++l) {
					for (int m = l; m < mi_nVoigt; ++m) {
						mv_elHes.at(j, k) += 2. * dLdy.at(m, j) * E.at(m, l) * dLdy.at(l, k) * numInt;
					}
				}
			}
		}

		// Green-Lagrange Strain Second Derivative
		for (int j = 0; j < mi_nDof; j = j + nDim) {
			for (int k = 0; k < mi_nDof; k = k + nDim) {
				double d2Ldy2[mi_nVoigt];

				for (int m = 0; m < nDim; ++m) {
					d2Ldy2[m] = 0.;

					for (int n = 0; n < nDim; ++n) {
						for (int o = 0; o < nDim; ++o) {
							d2Ldy2[m] += (F0(n, m) * *(mi_pShape + p + j + n)) * (F0(o, m) * *(mi_pShape + p + k + o));
						}
					}

					for (int n = m + 1; n < nDim; ++n) {
						int  v = mi_nVoigt - m - n;
						d2Ldy2[v] = 0.;

						for (int l = 0; l < nDim; ++l) {
							for (int o = 0; o < nDim; ++o) {
								d2Ldy2[v] += (F0(l, m) * *(mi_pShape + p + j + l)) * (F0(o, n) * *(mi_pShape + p + k + o));
							}
						}
					}
				}

				for (int m = 0; m < nDim; ++m) {
					for (int n = 0; n < nDim; ++n) {
						if (j <= k) mv_elHes.at(j + n, k + n) += SM.at(m) * d2Ldy2[m] * numInt;
					}
				}

				for (int m = nDim; m < mi_nVoigt; ++m) {
					for (int n = 0; n < nDim; ++n) {
						if (j <= k) mv_elHes.at(j + n, k + n) += 2. * SM.at(m) * d2Ldy2[m] * numInt;
					}
				}
			}
		}
	}
}

template<int nDim>
void O2P2::Proc::Comp::MeshElem_SVK<nDim>::addMassContrib(const double& mult)
{
	auto mi_pShape = mv_pElem->getShapeFc();
	auto mi_pWeight = mv_pElem->getWeight();

	auto mi_nNodes = mv_pElem->getNumNodes();
	auto mi_nIP = mv_pElem->getNumIP();

	for (int i = 0; i < mi_nIP; ++i) {
		int p = i * mi_nNodes;
		double numInt = mult * *(mi_pWeight + i) * mv_matPoint[i]->getJacobian();

		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = j; k < mi_nNodes; ++k) {
				for (int m = 0; m < nDim; ++m) {
					mv_elHes.at(j * nDim + m, k * nDim + m) += *(mi_pShape + p + j) * *(mi_pShape + p + k) * numInt;
				}
			}
		}
	}
}