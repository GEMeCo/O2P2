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
template class ElemComponent<3>;

template void ElemComponent<2>::setMaterialPoint();
template void ElemComponent<3>::setMaterialPoint();

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
			//		F0(k, j) += node->getInitPos()[j] * *(m_pElem->getShapeDerivative() + n + m);
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
// Implementation of ElemComp Member Function (3D only): getConstitutiveMatrix
//
// ================================================================================================
template<>
Eigen::MatrixXd ElemComponent<3>::getConstitutiveMatrix()
{
	// Dimensionality of the element
	auto nDim = m_pElem->getDIM();

	// Downcasting to solid element.
	ElementSolid* pElem = static_cast<ElementSolid*>(m_pElem.get());

	// Pointer to the material
	auto pMat = pElem->getMaterial();

	// Should also downcast to SVK material
	Mat_SVK_ISO* pSVK_Mat = static_cast<Mat_SVK_ISO*>(pMat);

	// Size of Constitutive matrix
	int nVoigt = 3 * nDim - 3;

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
				FInt(curDof) += dUedy * numInt;
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


// ================================================================================================
//
// Implementation of ElemComp Member Function (2D only): getContribution
//
// ================================================================================================
/*template<> void ElemComponent<2>::getContribution_SVK_ISO(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto pShape = m_pElem->getShapeDerivative();
	auto pWeight = m_pElem->getWeight();

	auto nNodes = m_pElem->getNumNodes();
	auto nIP = m_pElem->getNumIP();
	auto nDim = m_pElem->getDIM();

	// Matrices sizes
	int nVoigt = 3 * nDim - 3;
	int nDof = nNodes * nDim;

	// Elasticity tensor in material description
	Eigen::MatrixXd E(nVoigt, nVoigt);
	E = this->getConstitutiveMatrix();

	// Contribution from each material point
	for (int i = 0; i < nIP; ++i) {
		// Deformation gradient from undeformed to non-dimensional space
		Eigen::Matrix2d F0i;

		// The reference jacobian matrix is already inversed in the material point.
		F0i = m_MatPoint[i]->getJacobianMatrix();

		// Deformation gradient from deformed to non-dimensional space
		Eigen::Matrix2d F1 = Eigen::Matrix2d::Zero();

		// Some pointer arithmetic is required
		// p is related to current integration point
		int p = i * nDof;

		for (int j = 0; j < nDim; ++j) {
			// m is related to the nDof -> nNodes * nDim
			int m = 0;

			for (auto& node : getConectivity()) {
				for (int k = 0; k < nDim; ++k) {
					F1(j, k) += node->getTrial()[j] * *(pShape + p + m);
					m++;
				}
			}
		}
		// Deformation gradient
		Eigen::Matrix2d F;

		// evaluate F = F1 * F0i (F0i -> inverse of jacobian matrix)
		F = F1 * F0i;

		// Right Cauchy-Green (or Green deformation tensor)
		Eigen::Matrix2d C;
		C = F.transpose() * F;

		// Green-Lagrange Strain in Voigt notation
		Eigen::Vector3d L;

		for (int j = 0; j < nDim; ++j) {
			L(j) = 0.5 * (C(j, j) - 1.);

			for (int k = j + 1; k < nDim; ++k) {
				int l = 3 * nDim - 3 - j - k;
				L(l) = 0.5 * C(j, k);
			}
		}

		// Second Piola-Kirchhoff Stress
		Eigen::Vector3d S;
		S = E * L;

		// Green-Lagrange Strain Derivative
		Eigen::MatrixXd dLdy = Eigen::MatrixXd::Zero(nVoigt, nDof);

		// Specif energy derivative
		double dUedy[2];

		// Some pointer arithmetic is required
		for (int j = 0; j < nDof; j = j + nDim) {
			dLdy(0, j) = 0.5 * (F(0, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)) + F(0, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)));
			dLdy(1, j) = 0.5 * (F(0, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)) + F(0, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)));
			dLdy(2, j) = 0.5 * (F(0, 0) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)) + F(0, 1) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)));

			dLdy(0, j + 1) = 0.5 * (F(1, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)) + F(1, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)));
			dLdy(1, j + 1) = 0.5 * (F(1, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)) + F(1, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)));
			dLdy(2, j + 1) = 0.5 * (F(1, 0) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)) + F(1, 1) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)));

			dUedy[0] = S(0) * dLdy(0, j) + S(1) * dLdy(1, j) + 2. * (S(2) * dLdy(2, j));
			dUedy[1] = S(0) * dLdy(0, j + 1) + S(1) * dLdy(1, j + 1) + 2. * (S(2) * dLdy(2, j + 1));

			// Integration of the specific energy derivative = internal force
			FInt(j) += dUedy[0] * *(pWeight + i) * m_MatPoint[i]->getJacobian();
			FInt(j + 1) += dUedy[1] * *(pWeight + i) * m_MatPoint[i]->getJacobian();
		}

		Eigen::MatrixXd temporary(nVoigt, nDof);

		// The transversal modulus must be multiplied by 2 only here for Hessian
		// since I am disregarding the last term (thus reducing all vector/matrices sizes).
		for (int i = nDim; i < nVoigt; ++i) {
			E(i, i) = E(i, i) * 2.;
		}

		temporary = E * dLdy;
		temporary *= *(pWeight + i) * m_MatPoint[i]->getJacobian();
		Hessian += dLdy.transpose() * temporary;

		// Returning the transversal modulus for its original value
		for (int i = nDim; i < nVoigt; ++i) {
			E(i, i) = E(i, i) * 0.5;
		}

		// Green-Lagrange Strain Second Derivative
		Eigen::MatrixXd d2Ldy2(nVoigt, nDim);

		for (int j = 0; j < nDof; j = j + nDim) {
			for (int k = 0; k < nDof; k = k + nDim) {
				d2Ldy2(0, 0) = (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)) * (F0i(0, 0) * *(pShape + p + k) + F0i(1, 0) * *(pShape + p + k + 1));
				d2Ldy2(1, 0) = (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1)) * (F0i(0, 1) * *(pShape + p + k) + F0i(1, 1) * *(pShape + p + k + 1));
				d2Ldy2(2, 0) = (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1)) * (F0i(0, 1) * *(pShape + p + k) + F0i(1, 1) * *(pShape + p + k + 1));

				d2Ldy2(0, 1) = d2Ldy2(0, 0);
				d2Ldy2(1, 1) = d2Ldy2(1, 0);
				d2Ldy2(2, 1) = d2Ldy2(2, 0);

				Hessian(j, k) += (S(0) * d2Ldy2(0, 0) + S(1) * d2Ldy2(1, 0) + 2. * S(2) * d2Ldy2(2, 0)) * *(pWeight + i) * m_MatPoint[i]->getJacobian();
				Hessian(j + 1, k + 1) += (S(0) * d2Ldy2(0, 1) + S(1) * d2Ldy2(1, 1) + 2. * S(2) * d2Ldy2(2, 1)) * *(pWeight + i) * m_MatPoint[i]->getJacobian();
			}
		}
	}
};*/

// ================================================================================================
//
// Implementation of ElemComp Member Function (3D only): getContribution
//
// ================================================================================================
/*template<> void ElemComponent<3>::getContribution_SVK_ISO(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto pShape = m_pElem->getShapeDerivative();
	auto pWeight = m_pElem->getWeight();

	auto nNodes = m_pElem->getNumNodes();
	auto nIP = m_pElem->getNumIP();
	auto nDim = m_pElem->getDIM();

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
			for (int m = 0; m < nNodes; ++m) {
				auto node = getConectivity(m);

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

		// Specif energy derivative
		double dUedy[3];

		for (int j = 0; j < nDof; j = j + nDim) {
			dLdy(0, j) = 0.5 * (F(0, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)) + F(0, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));
			dLdy(1, j) = 0.5 * (F(0, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) + F(0, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)));
			dLdy(2, j) = 0.5 * (F(0, 2) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(0, 2) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)));
			dLdy(3, j) = 0.5 * (F(0, 1) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(0, 2) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)));
			dLdy(4, j) = 0.5 * (F(0, 0) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(0, 2) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));
			dLdy(5, j) = 0.5 * (F(0, 0) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) + F(0, 1) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));

			dLdy(0, j + 1) = 0.5 * (F(1, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)) + F(1, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));
			dLdy(1, j + 1) = 0.5 * (F(1, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) + F(1, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)));
			dLdy(2, j + 1) = 0.5 * (F(1, 2) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(1, 2) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)));
			dLdy(3, j + 1) = 0.5 * (F(1, 1) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(1, 2) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)));
			dLdy(4, j + 1) = 0.5 * (F(1, 0) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(1, 2) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));
			dLdy(5, j + 1) = 0.5 * (F(1, 0) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) + F(1, 1) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));

			dLdy(0, j + 2) = 0.5 * (F(2, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)) + F(2, 0) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));
			dLdy(1, j + 2) = 0.5 * (F(2, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) + F(2, 1) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)));
			dLdy(2, j + 2) = 0.5 * (F(2, 2) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(2, 2) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)));
			dLdy(3, j + 2) = 0.5 * (F(2, 1) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(2, 2) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)));
			dLdy(4, j + 2) = 0.5 * (F(2, 0) * (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) + F(2, 2) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));
			dLdy(5, j + 2) = 0.5 * (F(2, 0) * (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) + F(2, 1) * (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)));

			dUedy[0] = S(0) * dLdy(0, j) + S(1) * dLdy(1, j) + S(2) * dLdy(2, j) + 2. * (S(3) * dLdy(3, j) + S(4) * dLdy(4, j) + S(5) * dLdy(5, j));
			dUedy[1] = S(0) * dLdy(0, j + 1) + S(1) * dLdy(1, j + 1) + S(2) * dLdy(2, j + 1) + 2. * (S(3) * dLdy(3, j + 1) + S(4) * dLdy(4, j + 1) + S(5) * dLdy(5, j + 1));
			dUedy[2] = S(0) * dLdy(0, j + 2) + S(1) * dLdy(1, j + 2) + S(2) * dLdy(2, j + 2) + 2. * (S(3) * dLdy(3, j + 2) + S(4) * dLdy(4, j + 2) + S(5) * dLdy(5, j + 2));

			//	// Integration of the specific energy derivative = internal force
			FInt(j) += dUedy[0] * *(pWeight + i) * m_MatPoint[i]->getJacobian();
			FInt(j + 1) += dUedy[1] * *(pWeight + i) * m_MatPoint[i]->getJacobian();
			FInt(j + 2) += dUedy[2] * *(pWeight + i) * m_MatPoint[i]->getJacobian();
		}

		Eigen::MatrixXd temporary(nVoigt, nDof);

		// The transversal modulus must be multiplied by 2 only here for Hessian
		// since I am disregarding the last term (thus reducing all vector/matrices sizes).
		for (int i = nDim; i < nVoigt; ++i) {
			E(i, i) = E(i, i) * 2.;
		}

		temporary = E * dLdy;
		temporary *= numInt;
		Hessian += dLdy.transpose() * temporary;

		// Returning the transversal modulus for its original value
		for (int i = nDim; i < nVoigt; ++i) {
			E(i, i) = E(i, i) * 0.5;
		}

		// Green-Lagrange Strain Second Derivative
		Eigen::VectorXd d2Ldy2(nVoigt);
		
		for (int j = 0; j < nDof; j = j + nDim) {
			for (int k = 0; k < nDof; k = k + nDim) {
				d2Ldy2(0) = (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)) * (F0i(0, 0) * *(pShape + p + k) + F0i(1, 0) * *(pShape + p + k + 1) + F0i(2, 0) * *(pShape + p + k + 2));
				d2Ldy2(1) = (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) * (F0i(0, 1) * *(pShape + p + k) + F0i(1, 1) * *(pShape + p + k + 1) + F0i(2, 1) * *(pShape + p + k + 2));
				d2Ldy2(2) = (F0i(0, 2) * *(pShape + p + j) + F0i(1, 2) * *(pShape + p + j + 1) + F0i(2, 2) * *(pShape + p + j + 2)) * (F0i(0, 2) * *(pShape + p + k) + F0i(1, 2) * *(pShape + p + k + 1) + F0i(2, 2) * *(pShape + p + k + 2));
				d2Ldy2(3) = (F0i(0, 1) * *(pShape + p + j) + F0i(1, 1) * *(pShape + p + j + 1) + F0i(2, 1) * *(pShape + p + j + 2)) * (F0i(0, 2) * *(pShape + p + k) + F0i(1, 2) * *(pShape + p + k + 1) + F0i(2, 2) * *(pShape + p + k + 2));
				d2Ldy2(4) = (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)) * (F0i(0, 2) * *(pShape + p + k) + F0i(1, 2) * *(pShape + p + k + 1) + F0i(2, 2) * *(pShape + p + k + 2));
				d2Ldy2(5) = (F0i(0, 0) * *(pShape + p + j) + F0i(1, 0) * *(pShape + p + j + 1) + F0i(2, 0) * *(pShape + p + j + 2)) * (F0i(0, 1) * *(pShape + p + k) + F0i(1, 1) * *(pShape + p + k + 1) + F0i(2, 1) * *(pShape + p + k + 2));

				Hessian(j, k) += (S(0) * d2Ldy2(0) + S(1) * d2Ldy2(1) + S(2) * d2Ldy2(2) + 2. * (S(3) * d2Ldy2(3) + S(4) * d2Ldy2(4) + S(5) * d2Ldy2(5))) * numInt;
				Hessian(j + 1, k + 1) += (S(0) * d2Ldy2(0) + S(1) * d2Ldy2(1) + S(2) * d2Ldy2(2) + 2. * (S(3) * d2Ldy2(3) + S(4) * d2Ldy2(4) + S(5) * d2Ldy2(5))) * numInt;
				Hessian(j + 2, k + 2) += (S(0) * d2Ldy2(0) + S(1) * d2Ldy2(1) + S(2) * d2Ldy2(2) + 2. * (S(3) * d2Ldy2(3) + S(4) * d2Ldy2(4) + S(5) * d2Ldy2(5))) * numInt;
			}
		}
	}
};*/