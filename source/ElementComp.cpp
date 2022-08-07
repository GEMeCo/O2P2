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
// Element Component
// 
// Handler of element routines and local (element) matrices
// 
// ================================================================================================
#include "ElementComp.h"

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class ElemComponent<2>;

// ================================================================================================
//
// Implementation of Template Member Function: setMaterialPoint
//
// ================================================================================================
template<int nDim> void ElemComponent<nDim>::setMaterialPoint()
{
	int nIP = m_pElem->getNumIP();

	// One material point for each integration point
	for (int i = 0; i < nIP; ++i) {
		// Once created, evaluates F0 (Ref. Jacobian)
		Eigen::MatrixXd F0 = Eigen::MatrixXd::Zero(nDim, nDim);
		// Once created, evaluates F0 (Ref. Jacobian)
		//Eigen::MatrixXd F02 = Eigen::MatrixXd::Zero(nDim, nDim);

		// Some pointer arithmetic is required
		int n = i * m_pElem->getNumNodes() * nDim;

		for (int j = 0; j < nDim; ++j) {
			for (int m = 0; m < m_pElem->getNumNodes(); ++m) {
				auto node = getConectivity(m);

				for (int k = 0; k < nDim; ++k) {
					F0(j, k) += node->getInitPos()[j] * *(m_pElem->getShapeDerivative() + n + m*nDim + k);
				}
			}

			//int m = 0;
			//for (auto& node : getConectivity()) {
			//	for (int k = 0; k < nDim; k++) {
			//		F02(k, j) += node->getInitPos()[j] * *(m_pElem->getShapeDerivative() + n + m);
			//		m++;
			//	}
			//}
		}

		// Record the Jacobian Matrix on every Material (Integration) point
		this->m_MatPoint.emplace_back(std::make_unique<MaterialPoint>(F0));
	}
};


// ================================================================================================
//
// Implementation of ElemComp Member Function (2D only): getConstitutiveMatrix
//
// ================================================================================================
template<>
Eigen::MatrixXd ElemComponent<2>::getConstitutiveMatrix()
{
	// Dimensionality of the element
	auto nDim = m_pElem->getDIM();

	// Downcasting to plane element.
	ElementPlane* pElem = static_cast<ElementPlane*> (m_pElem.get());

	// Pointer to section
	auto pSec = pElem->getSection();

	// Plane State Type
	auto PS = pSec->getPS();

	// Pointer to the material
	auto pMat = pElem->getMaterial();

	// Should also downcast to SVK material
	Mat_SVK_ISO* pSVK_Mat = static_cast<Mat_SVK_ISO*> (pMat);

	// Size of Constitutive matrix
	int nVoigt = 3 * nDim - 3;

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

	return E;
};


// ================================================================================================
//
// Implementation of ElemComp Member Function: getContribution
//
// ================================================================================================
template<int nDim>
void ElemComponent<nDim>::getContribution(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian)
{
	auto pMat = m_pElem->getMaterial();

	if (pMat->getMaterialType() == MaterialType::SVK_ISO)
		this->getContribution_SVK_ISO(FInt, Hessian);
};


// ================================================================================================
//
// Specialization of ElemComp Member Function: getContribution_SVK_ISO
//
// ================================================================================================
template<int nDim>
void ElemComponent<nDim>::getContribution_SVK_ISO(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto pShape = m_pElem->getShapeDerivative();
	auto pWeight = m_pElem->getWeight();

	auto nNodes = m_pElem->getNumNodes();
	auto nIP = m_pElem->getNumIP();

	// Matrices sizes
	int nVoigt = 3 * nDim - 3;
	int nDof = nNodes * nDim;

	// Elasticity tensor in material description
	Eigen::MatrixXd E(nVoigt, nVoigt);
	E = this->getConstitutiveMatrix();

	// Contribution from each material point
	for (int i = 0; i < nIP; ++i) {
		// Deformation gradient from undeformed to non-dimensional space
		Eigen::MatrixXd F0i(nDim, nDim);

		// The reference jacobian matrix is already inversed in the material point.
		F0i = m_MatPoint[i]->getJacobianMatrix();

		// Deformation gradient from deformed to non-dimensional space
		Eigen::MatrixXd F1 = Eigen::MatrixXd::Zero(nDim, nDim);

		// Some pointer arithmetic is required
		// p and numInt are related to current integration point and weighting numerical integration
		int p = i * nDof;
		double numInt = *(pWeight + i) * m_MatPoint[i]->getJacobian();

		for (int j = 0; j < nDim; ++j) {

			for (int m = 0; m < m_pElem->getNumNodes(); ++m) {
				auto node = getConectivity(m);

				for (int k = 0; k < nDim; ++k) {
					F1(j, k) += node->getTrial()[j] * *(pShape + p + m*nDim + k);
					//F1(k, j) += node->getTrial()[j] * *(pShape + p + m * nDim + k);
				}
			}

			// m is related to the nDof -> nNodes * nDim
			//int m = 0;
			//for (auto& node : getConectivity()) {
			//	for (int k = 0; k < nDim; ++k) {
			//		F1(j, k) += node->getTrial()[j] * *(pShape + p + m);
			//		m++;
			//	}
			//}
		}

		// Deformation gradient
		Eigen::MatrixXd F(nDim, nDim);

		// evaluate F = F1 * F0i (F0i -> inverse of jacobian matrix)
		F = F1 * F0i;

		// Right Cauchy-Green (or Green deformation tensor)
		Eigen::MatrixXd C(nDim, nDim);
		C = F.transpose() * F;

		// Thermal, AAR, Damage and other stuff should be implemented here

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
				FInt(curDof) += dUedy * numInt;
			}
		}

		Eigen::MatrixXd temporary(nVoigt, nDof);

		// The transversal modulus must be multiplied by 2 only here for Hessian
		// since I am disregarding the last term (thus reducing all vector/matrices sizes).
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

		// We still have to add the Green-Lagrange strain second derivative to the Hessian Matrix
		// As it is the same in every direction, only one must be evaluated

		// Green-Lagrange Strain Second Derivative
		Eigen::VectorXd d2Ldy2(nVoigt);

		for (int r = 0; r < nNodes; ++r) {
			int curDof_1 = r * nDim;

			for (int s = 0; s < nNodes; ++s) {
				int curDof_2 = s * nDim;

				d2Ldy2.setZero();

				for (int m = 0; m < nDim; ++m) {
					for (int l = 0; l < nDim; ++l) {
						for (int o = 0; o < nDim; ++o) {
							d2Ldy2(m) += (F0i(l, m) * *(pShape + p + curDof_1 + l)) * (F0i(o, m) * *(pShape + p + curDof_2 + o));
						}
					}

					for (int n = m + 1; n < nDim; ++n) {
						int v = nVoigt - m - n;

						for (int l = 0; l < nDim; ++l) {
							for (int o = 0; o < nDim; ++o) {
								d2Ldy2(v) += (F0i(l, m) * *(pShape + p + curDof_1 + l)) * (F0i(o, n) * *(pShape + p + curDof_2 + o));
							}
						}
					}

					for (int v = 0; v < nDim; ++v) {
						Hessian(curDof_1 + m, curDof_2 + m) += S(v) * d2Ldy2(v) * numInt;
					}

					for (int v = nDim; v < nVoigt; ++v) {
						Hessian(curDof_1 + m, curDof_2 + m) += 2. * S(v) * d2Ldy2(v) * numInt;
					}
				}
			}
		}
	}

	return;
};
