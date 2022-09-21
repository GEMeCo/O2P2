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
// Mesh Element Component
// 
// Handler of element routines and local (element) matrices
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

// ================================================================================================
//
// Implementation of Template Member Function: setMaterialPoint
//
// ================================================================================================
template<int nDim> void O2P2::Proc::Comp::MeshElem_SVK<nDim>::setMaterialPoint()
{
	int nIP = m_pElem->getNumIP();
	int nElDim = m_pElem->getDIM();
	
	// The BaseElement pointer does not have access to getConectivity
	// thus, must cast it to a derived class
	O2P2::Prep::Elem::Element<nDim>* pElem = static_cast<O2P2::Prep::Elem::Element<nDim>*> (m_pElem.get());

	// One material point for each integration point
	for (int i = 0; i < nIP; ++i) {
		// Once created, evaluates F0 (Ref. Jacobian)
		Eigen::MatrixXd F0 = Eigen::MatrixXd::Zero(nDim, nElDim);

		// Some pointer arithmetic is required
		int n = i * m_pElem->getNumNodes() * nDim;

		for (int j = 0; j < nDim; ++j) {
			for (int m = 0; m < m_pElem->getNumNodes(); ++m) {
				auto node = pElem->getConectivity(m);

				for (int k = 0; k < nElDim; ++k) {
					F0(j, k) += node->getInitPos()[j] * *(m_pElem->getShapeDerivative() + n + m*nElDim + k);
				}
			}
		}

		// Record the Jacobian Matrix on every Material (Integration) point
		this->m_matPoint.emplace_back(std::make_unique<O2P2::Proc::Comp::MaterialPoint>(F0));
	}
};


// ================================================================================================
//
// Implementation of MeshElem_SVK Member Function: getYoungModulus
//
// ================================================================================================
template<int nDim>
double O2P2::Proc::Comp::MeshElem_SVK<nDim>::getYoungModulus()
{
	assert((m_pElem->getDIM() == 1) && "MeshElem_SVK<nDim>::getYoungModulus() should only be used with linear elements");

	// Downcasting to linear element.
	O2P2::Prep::Elem::ElementLinear<nDim>* pElem = static_cast<O2P2::Prep::Elem::ElementLinear<nDim>*> (m_pElem.get());

	// Pointer to section
	auto pSec = pElem->getSection();

	// Pointer to the material
	auto pMat = pElem->getMaterial();

	// Should also downcast to SVK material
	O2P2::Prep::Mat_SVK_ISO* pSVK_Mat = static_cast<O2P2::Prep::Mat_SVK_ISO*> (pMat);

	double E = pSVK_Mat->getLongitudinalModulus() * pSec->getSection();

	return E;
}


// ================================================================================================
//
// Specialization of MeshElem_SVK Member Function (2D only): getConstitutiveMatrix
//
// ================================================================================================
template<>
Eigen::MatrixXd O2P2::Proc::Comp::MeshElem_SVK<2>::getConstitutiveMatrix()
{
	// Dimensionality of the element
	auto nElDim = m_pElem->getDIM();

	// Downcasting to plane element.
	O2P2::Prep::Elem::ElementPlane* pElem = static_cast<O2P2::Prep::Elem::ElementPlane*> (m_pElem.get());

	// Pointer to section
	auto pSec = pElem->getSection();

	// Plane State Type
	auto PS = pSec->getPS();

	// Thickness
	auto tk = pSec->getSection();

	// Pointer to the material
	auto pMat = pElem->getMaterial();

	// Should also downcast to SVK material
	O2P2::Prep::Mat_SVK_ISO* pSVK_Mat = static_cast<O2P2::Prep::Mat_SVK_ISO*> (pMat);

	// Size of Constitutive matrix
	int nVoigt = 3 * nElDim - 3;

	double Yg = pSVK_Mat->getLongitudinalModulus();
	double nu = pSVK_Mat->getPoisson();
	double G = pSVK_Mat->getTransversalModulus();

	Eigen::MatrixXd E(nVoigt, nVoigt);
	E.setZero();

	if (PS == PlaneStateType::PLANE_STRESS) {

		double Aux = Yg / (1. - nu * nu);

		E(0, 0) = Aux;
		E(1, 1) = Aux;

		E(0, 1) = nu * Aux;
		E(1, 0) = E(0, 1);

		E(2, 2) = G * 2.;
	}
	else {
		double Aux = Yg / ((1. - 2. * nu) * (1. + nu));

		E(0, 0) = (1. - nu) * Aux;
		E(1, 1) = (1. - nu) * Aux;

		E(0, 1) = nu * Aux;
		E(1, 0) = E(0, 1);

		E(2, 2) = G * 2.;
	}

	return E * tk;
};


// ================================================================================================
//
// Specialization of MeshElem_SVK Member Function (3D only): getConstitutiveMatrix
//
// ================================================================================================
template<>
Eigen::MatrixXd O2P2::Proc::Comp::MeshElem_SVK<3>::getConstitutiveMatrix()
{
	// Dimensionality of the element
	auto nElDim = m_pElem->getDIM();

	// Downcasting to solid element.
	O2P2::Prep::Elem::ElementSolid* pElem = static_cast<O2P2::Prep::Elem::ElementSolid*>(m_pElem.get());

	// Pointer to the material
	auto pMat = pElem->getMaterial();

	// Should also downcast to SVK material
	O2P2::Prep::Mat_SVK_ISO* pSVK_Mat = static_cast<O2P2::Prep::Mat_SVK_ISO*>(pMat);

	// Size of Constitutive matrix
	int nVoigt = 3 * nElDim - 3;

	double Yg = pSVK_Mat->getLongitudinalModulus();
	double nu = pSVK_Mat->getPoisson();
	double G = pSVK_Mat->getTransversalModulus();

	Eigen::MatrixXd E(nVoigt, nVoigt);
	E.setZero();

	double Aux = Yg / (1. + nu) / (1. - 2 * nu);

	E(0, 0) = Aux * (1 - nu);
	E(1, 1) = E(0, 0);
	E(2, 2) = E(0, 0);

	E(0, 1) = nu * Aux;
	E(0, 2) = E(0, 1);

	E(1, 0) = E(0, 1);
	E(1, 2) = E(0, 1);

	E(2, 0) = E(0, 1);
	E(2, 1) = E(0, 1);

	E(3, 3) = G * 2.;
	E(4, 4) = G * 2.;
	E(5, 5) = G * 2.;

	return E;
};


// ================================================================================================
//
// Implementation of MeshElem_SVK Member Function: getContribution
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshElem_SVK<nDim>::getContribution(Eigen::MatrixXd& Hessian)
{
	auto pMat = m_pElem->getMaterial();
	auto nElDim = m_pElem->getDIM();

	if (nElDim == 1) {
		if (pMat->getMaterialType() == MaterialType::SVK_ISO)
			this->getContribution_SVK_ISO<1>(Hessian);
	}

	if (nElDim == 2) {
		if (pMat->getMaterialType() == MaterialType::SVK_ISO)
			this->getContribution_SVK_ISO<2>(Hessian);
	}

	if (nElDim == 3) {
		if (pMat->getMaterialType() == MaterialType::SVK_ISO)
			this->getContribution_SVK_ISO<3>(Hessian);
	}
};


// ================================================================================================
//
// Specialization of MeshElem_SVK Member Function: getContribution_SVK_ISO
//
// ================================================================================================
template<> template<>
void O2P2::Proc::Comp::MeshElem_SVK<2>::getContribution_SVK_ISO<1>(Eigen::MatrixXd& Hessian)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto pShape = m_pElem->getShapeDerivative();
	auto pWeight = m_pElem->getWeight();

	auto nNodes = m_pElem->getNumNodes();
	auto nIP = m_pElem->getNumIP();

	const int nDim = 2;

	// One material point for each integration point
	for (int i = 0; i < nIP; ++i) {

		// For fiber elements, F0 and F1 holds the initial and current length in each direction
		Eigen::MatrixXd F1 = Eigen::MatrixXd::Zero(nDim, 1);

		// Some pointer arithmetic is required
		int n = i * m_pElem->getNumNodes() * nDim;

		for (int j = 0; j < nDim; ++j) {
			int m = 0;

			for (auto& node : v_Conect) {
				F1(j, 0) += node->getTrial()[j] * *(m_pElem->getShapeDerivative() + n + m);
				m++;
			}
		}

		double dA0 = m_matPoint[i]->getJacobian();
		double dA1 = 0.;

		for (int j = 0; j < nDim; ++j) {
			dA1 += F1(j, 0) * F1(j, 0);
		}

		// Green-Lagrange Strain
		double L = 0.5 * (dA1 - dA0) * dA0;

		// Elasticity properties
		double E = this->getYoungModulus();

		// Second Piola-Kirchhoff Stress
		double S = E * L;

		// Green-Lagrange Strain Derivative
		Eigen::VectorXd dLdy = Eigen::VectorXd::Zero(nDim);

		for (int j = 0; j < nDim; ++j) {
			dLdy(j) = F1(j, 0) * F1(j, 0) * dA0;
		}

		// Green-Lagrange Strain Second Derivative
		// It is the same in every direction, and equal to dA0
		double d2Ldy2 = dA0;

		// Specific energy derivative
		double dUedy = S * d2Ldy2;

		// Weight for Gauss point and Jacobian
		double numInt = *(pWeight + i) * this->m_matPoint[i]->getJacobian();

		// Integration of the specific energy derivative = internal force
		for (int j = 0; j < nNodes; ++j) {
			m_elFor(j) += dUedy * *(pShape + n + j) * numInt;
		}

		// Auxiliary for Hessian matrix
		Eigen::MatrixXd auxH = Eigen::MatrixXd::Zero(nDim, nDim);

		for (int j = 0; j < nDim; ++j) {
			for (int k = 0; k < nDim; ++k) {
				auxH(j, k) = E * dLdy(j) * dLdy(k);
			}
			auxH(j, j) += S * d2Ldy2;
		}

		// First part of Hessian
		Eigen::MatrixXd H1 = Eigen::MatrixXd::Zero(nDim, nDim * nNodes);

		// Kronecker delta
		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nDim, nDim);

		for (int j = 0; j < nNodes; ++j) {
			for (int k = 0; k < nDim; ++k) {
				for (int l = 0; l < nDim; ++l) {
					for (int m = 0; m < nDim; ++m) {
						H1(l, j * nDim + k) += auxH(l, m) * I(m, k) * *(pShape + n + j);
					}
				}
			}
		}

		// Local Hessian
		for (int j = 0; j < nNodes; ++j) {
			for (int k = 0; k < nDim; ++k) {

				for (int l = 0; l < nNodes; ++l) {
					for (int m = 0; m < nDim; ++m) {

						for (int o = 0; o < nDim; ++o) {
							Hessian(j * nNodes + k, l * nNodes + m) += *(pShape + n + j) * I(o, k) * H1(o, l * nNodes + m) * numInt;
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
void O2P2::Proc::Comp::MeshElem_SVK<nDim>::getContribution_SVK_ISO(Eigen::MatrixXd& Hessian)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto pShape = m_pElem->getShapeDerivative();
	auto pWeight = m_pElem->getWeight();

	auto nNodes = m_pElem->getNumNodes();
	auto nIP = m_pElem->getNumIP();

	// Matrices sizes
	const int nVoigt = 3 * nDim - 3;
	const int nDof = nNodes * nDim;

	// Elasticity tensor in material description
	Eigen::MatrixXd E(nVoigt, nVoigt);
	E = this->getConstitutiveMatrix();

	// Contribution from each material point
	for (int i = 0; i < nIP; ++i) {
		// Deformation gradient from undeformed to non-dimensional space
		Eigen::MatrixXd F0i(nDim, nDim);

		// The reference jacobian matrix is already inversed in the material point.
		F0i = m_matPoint[i]->getJacobianMatrix();

		// Deformation gradient from deformed to non-dimensional space
		Eigen::MatrixXd F1 = Eigen::MatrixXd::Zero(nDim, nDim);

		// Some pointer arithmetic is required
		// p and numInt are related to current integration point and weighting numerical integration
		int p = i * nDof;
		double numInt = *(pWeight + i) * m_matPoint[i]->getJacobian();

		for (int j = 0; j < nDim; ++j) {

			for (int m = 0; m < m_pElem->getNumNodes(); ++m) {
				auto node = v_Conect[m];

				for (int k = 0; k < nDim; ++k) {
					F1(j, k) += node->getTrial()[j] * *(pShape + p + m * nDim + k);
				}
			}
		}

		// Deformation gradient
		Eigen::MatrixXd F(nDim, nDim);

		// evaluate F = F1 * F0i (F0i -> inverse of jacobian matrix)
		F = F1 * F0i;

		// Right Cauchy-Green (or Green deformation tensor)
		Eigen::MatrixXd C(nDim, nDim);
		C = F.transpose() * F;

		// Green-Lagrange Strain in Voigt notation
		Eigen::VectorXd L(nVoigt);

		for (int j = 0; j < nDim; ++j) {
			L(j) = 0.5 * (C(j, j) - 1.);

			for (int k = j + 1; k < nDim; ++k) {
				int l = nVoigt - j - k;
				L(l) = 0.5 * C(j, k);
			}
		}

		// Second Piola-Kirchhoff Stress
		Eigen::VectorXd S(nVoigt);
		S = E * L;

		// Green-Lagrange Strain Derivative
		Eigen::MatrixXd dLdy = Eigen::MatrixXd::Zero(nVoigt, nDof);

		// Some pointer arithmetic is required
		for (int r = 0; r < nNodes; ++r) {
			for (int d = 0; d < nDim; ++d) {
				int curDof = r * nDim + d;

				for (int m = 0; m < nDim; ++m) {
					int n = m;

					for (int j = 0; j < nDim; ++j) {
						dLdy(m, curDof) += 0.5 * (F(d, m) * F0i(j, n) * *(pShape + p + r * nDim + j) + F0i(j, m) * F(d, n) * *(pShape + p + r * nDim + j));
					}

					for (n = m + 1; n < nDim; ++n) {
						int l = nVoigt - m - n;

						for (int j = 0; j < nDim; ++j) {
							dLdy(l, curDof) += 0.5 * (F(d, m) * F0i(j, n) * *(pShape + p + r * nDim + j) + F0i(j, m) * F(d, n) * *(pShape + p + r * nDim + j));
						}
					}
				}

				// Specific energy derivative - it is a vector, but we are not going to use it again.
				double dUedy = 0;

				for (int j = 0; j < nVoigt; ++j) {
					dUedy += S(j) * dLdy(j, curDof);
				}

				for (int j = nDim; j < nVoigt; ++j) {
					dUedy += S(j) * dLdy(j, curDof);
				}

				// Integration of the specific energy derivative = internal force
				m_elFor(curDof) += dUedy * numInt;
			}
		}

		Eigen::MatrixXd temporary(nVoigt, nDof);

		// The transversal modulus must be multiplied by 2 only here for Hessian
		// since I am disregarding the last terms (thus reducing all vector/matrices sizes).
		for (int j = nDim; j < nVoigt; ++j) {
			E(j, j) = E(j, j) * 2.;
		}

		temporary = E * dLdy;
		temporary *= numInt;
		Hessian += dLdy.transpose() * temporary;

		// Returning the transversal modulus for its original value
		for (int j = nDim; j < nVoigt; ++j) {
			E(j, j) = E(j, j) * 0.5;
		}

		// Green-Lagrange Strain Second Derivative
		for (int j = 0; j < nDof; j = j + nDim) {
			for (int k = 0; k < nDof; k = k + nDim) {
				//Eigen::VectorXd d2Ldy2 = Eigen::VectorXd::Zero(nVoigt);
				double d2Ldy2[nVoigt];

				for (int m = 0; m < nDim; ++m) {
					d2Ldy2[m] = 0.;

					for (int n = 0; n < nDim; ++n) {
						for (int o = 0; o < nDim; ++o) {
							d2Ldy2[m] += (F0i(n, m) * *(pShape + p + j + n)) * (F0i(o, m) * *(pShape + p + k + o));
						}
					}

					for (int n = m + 1; n < nDim; ++n) {
						int v = nVoigt - m - n;
						d2Ldy2[v] = 0.;

						for (int l = 0; l < nDim; ++l) {
							for (int o = 0; o < nDim; ++o) {
								d2Ldy2[v] += (F0i(l, m) * *(pShape + p + j + l)) * (F0i(o, n) * *(pShape + p + k + o));
							}
						}
					}
				}

				for (int m = 0; m < nDim; ++m) {
					for (int n = 0; n < nDim; ++n) {
						Hessian(j + n, k + n) += S(m) * d2Ldy2[m] * numInt;
					}
				}

				for (int m = nDim; m < nVoigt; ++m) {
					for (int n = 0; n < nDim; ++n) {
						Hessian(j + n, k + n) += 2. * S(m) * d2Ldy2[m] * numInt;
					}
				}
			}
		}
	}

	return;
};
