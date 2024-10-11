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
#include "Mesh Components\MeshElem.h"

// ================================================================================================
//
// Implementation of MeshElem Member Function: setMaterialPoint
//
// ================================================================================================
void O2P2::Proc::Comp::MeshElem::setMaterialPoint()
{
	int mi_nIP = mv_ptElem->getNumIP();
	int mi_ElDim = mv_ptElem->getDIM();

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {
		// Once created, evaluates F0 (Ref. Jacobian)
		M2S2::Dyadic2N mi_F0(mi_ElDim);

		// Some pointer arithmetic is required
		int n = i * mv_ptElem->getNumNodes() * mi_ElDim;

		for (int j = 0; j < mi_ElDim; ++j) {
			for (int m = 0; m < mv_ptElem->getNumNodes(); ++m) {
				auto mi_pNode = mv_ptElem->getConectivity(m);

				for (int k = 0; k < mi_ElDim; ++k) {
					mi_F0(j, k) += mi_pNode->getInitPos()[j] * *(mv_ptElem->getShapeDerivative() + n + m * mi_ElDim + k);
				}
			}
		}

		// Record the Jacobian Matrix on every Material (Integration) point
		this->mv_matPoint.emplace_back(std::make_unique<O2P2::Proc::Comp::MaterialPoint>(mi_F0));
	}
}


// ================================================================================================
//
// Implementation of MeshElem Member Function: setIndexing
//
// ================================================================================================
void O2P2::Proc::Comp::MeshElem::setIndexing()
{
	int mi_ElDim = mv_ptElem->getDIM();
	for (int i = 0; i < mv_conect.size(); ++i) {
		for (int j = 0; j < mi_ElDim; ++j) {
			mv_dofIndex.push_back(mv_conect.at(i)->mv_dofList.at(j));
		}
	}
}


// ================================================================================================
//
// Implementation of MeshElem Member Function: getContribution
//
// ================================================================================================
void O2P2::Proc::Comp::MeshElem::getContribution()
{
	auto mi_pMat = mv_ptElem->getMaterial();
	auto mi_ElDim = mv_ptElem->getDIM();

	switch (mi_pMat->getMaterialType())
	{
	case MaterialType::SVK_ISO:
	{
		double mi_section = 1.;

		if (mi_ElDim == 2) {
			// Downcasting element pointer to plane element pointer.
			assert(std::dynamic_pointer_cast<O2P2::Geom::Elem::ElemPlane>(mv_ptElem));
			mi_section = std::static_pointer_cast<O2P2::Geom::Elem::ElemPlane>(mv_ptElem)->getSection()->getSection();
		}

		this->getContribution_SVK(mi_section);
		break;
	}
	default: { throw std::invalid_argument("\n\n\nMeshElem::getContribution - The selected material type was not defined yet.\n\n\n"); break; }
	}
}


// ================================================================================================
//
// Implementation of MeshElem Member Function: getContribution_SVK
//
// ================================================================================================
void O2P2::Proc::Comp::MeshElem::getContribution_SVK(double section)
{
	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_ptElem->getShapeDerivative();
	auto mi_pWeight = mv_ptElem->getWeight();

	auto mi_ElDim = mv_ptElem->getDIM();
	auto mi_nNodes = mv_ptElem->getNumNodes();
	auto mi_nIP = mv_ptElem->getNumIP();
	auto mi_pMat = mv_ptElem->getMaterial();

	// Elasticity tensor in material description
	M2S2::Dyadic4C E(mi_ElDim);
	E = mi_pMat->getConstitutiveMatrix(mi_ElDim);

	const int mi_nVoigt = 3 * mi_ElDim - 3;
	const int mi_nSize = mi_ElDim * mi_ElDim;

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {

		// Deformation gradient from undeformed to non-dimensional space
		M2S2::Dyadic2N F0(mi_ElDim);
		F0 = mv_matPoint.at(i)->getJacobianMatrix();

		// Deformation gradient from deformed to non-dimensional space
		M2S2::Dyadic2N F1(mi_ElDim);

		// Some pointer arithmetic is required
		// p and numInt are related to current integration point and weighting numerical integration
		int p = i * mv_nDof;
		double numInt = *(mi_pWeight + i) * this->mv_matPoint.at(i)->getJacobian() * section;

		for (int j = 0; j < mi_ElDim; ++j) {
			for (int m = 0; m < mi_nNodes; ++m) {
				auto& node = mv_conect.at(m);

				for (int k = 0; k < mi_ElDim; ++k) {
					F1.at(j, k) += node->getTrialPos()[j] * *(mi_pShape + p + m * mi_ElDim + k);
				}
			}
		}

		// Deformation gradient
		M2S2::Dyadic2N F(mi_ElDim);

		// evaluate F = F1 * F0i (F0i -> inverse of jacobian matrix)
		F = F1 * F0;

		// Right Cauchy-Green (or Green deformation tensor)
		M2S2::Dyadic2S C(mi_ElDim);
		C = F.getATA();

		// Green-Lagrange Strain
		M2S2::Dyadic2S L(mi_ElDim);
		L = 0.5 * (C - M2S2::Dyadic2S::identity(mi_ElDim));

		// Second Piola-Kirchhoff Stress
		M2S2::Dyadic2S S(mi_ElDim);
		S = E * L;
		auto SM = S.getVoigtMnemonics();

		// Green-Lagrange Strain Derivative
		// Notice that dLdy is written in Voigt notation with mnemonic rule
		M2S2::MatrixX dLdy(mi_nVoigt, mv_nDof);

		for (int r = 0; r < mi_nNodes; ++r) {
			for (int d = 0; d < mi_ElDim; ++d) {
				int mi_curDof = r * mi_ElDim + d;

				for (int m = 0; m < mi_ElDim; ++m) {
					for (int j = 0; j < mi_ElDim; ++j) {
						dLdy(m, mi_curDof) += 0.5 * (F(d, m) * F0(j, m) * *(mi_pShape + p + r * mi_ElDim + j) + F0(j, m) * F(d, m) * *(mi_pShape + p + r * mi_ElDim + j));
					}

					for (int n = m + 1; n < mi_ElDim; ++n) {
						int l = mi_nVoigt - m - n;

						for (int j = 0; j < mi_ElDim; ++j) {
							dLdy(l, mi_curDof) += 0.5 * (F(d, m) * F0(j, n) * *(mi_pShape + p + r * mi_ElDim + j) + F0(j, m) * F(d, n) * *(mi_pShape + p + r * mi_ElDim + j));
						}
					}
				}

				// Specific energy derivative - it is a vector, but we are not going to use it again.
				double dUedy = 0.;

				for (int j = 0; j < mi_nVoigt; ++j) {
					dUedy += SM.at(j) * dLdy(j, mi_curDof);
				}

				for (int j = mi_ElDim; j < mi_nVoigt; ++j) {
					dUedy += SM.at(j) * dLdy(j, mi_curDof);
				}

				// Integration of the specific energy derivative = internal force
				mv_elFor.at(mi_curDof) += dUedy * numInt;
			}
		}

		//Hessian += dLdy.transpose() * E * dLdy * numInt;
		for (int j = 0; j < mv_nDof; ++j) {
			for (int k = j; k < mv_nDof; ++k) {
				for (int l = 0; l < mi_ElDim; ++l) {
					for (int m = l; m < mi_ElDim; ++m) {
						mv_elHes.at(j, k) += dLdy.at(m, j) * E.at(m, l) * dLdy.at(l, k) * numInt;
					}
				}
				for (int l = mi_ElDim; l < mi_nVoigt; ++l) {
					for (int m = l; m < mi_nVoigt; ++m) {
						mv_elHes.at(j, k) += 2. * dLdy.at(m, j) * E.at(m, l) * dLdy.at(l, k) * numInt;
					}
				}
			}
		}

		// Green-Lagrange Strain Second Derivative
		std::vector<double> d2Ldy2(mi_nSize);		// Non-symmetric

		for (int j = 0; j < mi_nNodes; ++j) {
			int p1 = p + j * mi_ElDim;
			for (int k = 0; k < mi_nNodes; ++k) {
				int p2 = p + k * mi_ElDim;

				// Fastest way to set a vector to zero
				memset(&d2Ldy2[0], 0., mi_nSize * sizeof(double));

				for (int o = 0; o < mi_ElDim; ++o) {

					for (int m = 0; m < mi_ElDim; ++m) {
						for (int n = 0; n < mi_ElDim; ++n) {
							d2Ldy2.at(o) += (F0(m, o) * *(mi_pShape + p1 + m)) * (F0(n, o) * *(mi_pShape + p2 + n));
						}
					}

					for (int p = o + 1; p < mi_ElDim; ++p) {
						for (int m = 0; m < mi_ElDim; ++m) {
							for (int n = 0; n < mi_ElDim; ++n) {
								d2Ldy2.at(mi_nVoigt - o - p) += (F0(m, o) * *(mi_pShape + p1 + m)) * (F0(n, p) * *(mi_pShape + p2 + n));
								d2Ldy2.at(mi_nSize - o - p) += (F0(m, p) * *(mi_pShape + p1 + m)) * (F0(n, o) * *(mi_pShape + p2 + n));
							}
						}
					}
				}

				for (int m = 0; m < mi_ElDim; ++m) {
					for (int n = 0; n < mi_nVoigt; ++n) {
						if (j * mi_ElDim <= k * mi_ElDim) mv_elHes.at(j* mi_ElDim + m, k * mi_ElDim + m) += SM.at(n) * d2Ldy2.at(n) * numInt;
					}

					for (int n = mi_nVoigt; n < mi_nSize; ++n) {
						if (j * mi_ElDim <= k * mi_ElDim) mv_elHes.at(j * mi_ElDim + m, k * mi_ElDim + m) += SM.at(n + mi_nVoigt - mi_nSize) * d2Ldy2.at(n) * numInt;
					}
				}
			}
		}
	}
}


// ================================================================================================
//
// Implementation of MeshElem Member Function: addMassContrib
//
// ================================================================================================
void O2P2::Proc::Comp::MeshElem::addMassContrib(const double& mult)
{
	auto mi_pShape = mv_ptElem->getShapeFc();
	auto mi_pWeight = mv_ptElem->getWeight();

	auto mi_ElDim = mv_ptElem->getDIM();
	auto mi_nNodes = mv_ptElem->getNumNodes();
	auto mi_nIP = mv_ptElem->getNumIP();
	auto mi_tk = 1.;

	int p;
	double numInt;

	if (mi_ElDim == 2) {
		// Downcasting to plane element.
		assert(std::dynamic_pointer_cast<O2P2::Geom::Elem::ElemPlane>(mv_ptElem));
		mi_tk = std::static_pointer_cast<O2P2::Geom::Elem::ElemPlane>(mv_ptElem)->getSection()->getSection();
	}

	for (int i = 0; i < mi_nIP; ++i) {
		p = i * mi_nNodes;
		numInt = *(mi_pWeight + i) * this->mv_matPoint.at(i)->getJacobian() * mult * mi_tk;

		for (int j = 0; j < mi_nNodes; ++j) {
			for (int k = j; k < mi_nNodes; ++k) {
				for (int m = 0; m < mi_ElDim; ++m) {
					this->mv_elHes.at(j * mi_ElDim + m, k * mi_ElDim + m) += *(mi_pShape + p + j) * *(mi_pShape + p + k) * numInt;
				}
			}
		}
	}
}