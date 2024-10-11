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
#include "Mesh Components\MeshIncl.h"

// ================================================================================================
//
// Explicit template class instantiation
//
// ================================================================================================
template class O2P2::Proc::Comp::MeshIncl<2>;
template class O2P2::Proc::Comp::MeshIncl<3>;

template void O2P2::Proc::Comp::MeshIncl<2>::setMaterialPoint();
template void O2P2::Proc::Comp::MeshIncl<3>::setMaterialPoint();

template void O2P2::Proc::Comp::MeshIncl<2>::setIndexing();
template void O2P2::Proc::Comp::MeshIncl<3>::setIndexing();

template void O2P2::Proc::Comp::MeshIncl<2>::getContribution_SVK(double section);
template void O2P2::Proc::Comp::MeshIncl<3>::getContribution_SVK(double section);

template M2S2::MatrixX O2P2::Proc::Comp::MeshIncl<2>::getSpreadMatrix();
template M2S2::MatrixX O2P2::Proc::Comp::MeshIncl<3>::getSpreadMatrix();


// ================================================================================================
//
// Implementation of MeshIncl Member Function: setMaterialPoint
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshIncl<nDim>::setMaterialPoint()
{
	int mi_ElDim = mv_ptElem->getDIM();

	if (nDim == mi_ElDim) {
		O2P2::Proc::Comp::MeshElem::setMaterialPoint();
	}
	else {
		auto mi_pShape = mv_ptElem->getShapeDerivative();
		int mi_nIP = mv_ptElem->getNumIP();

		// One material point for each integration point
		for (int i = 0; i < mi_nIP; ++i) {
			// Once created, evaluates F0 (Ref. Jacobian)
			M2S2::Dyadic2N mi_F0(mi_ElDim);
			M2S2::Dyadic2N mi_rot(nDim);

			std::array<double, 2> mi_cd{ 0, 0 };
			std::array<double, 3> mi_vectX{ 0, 0, 0 }, mi_vectY{ 0, 0, 0 }, mi_vectZ{ 0, 0, 0 };

			// Some pointer arithmetic is required
			int n = i * mv_ptElem->getNumNodes() * mi_ElDim;
			int m = 0;

			double mi_normX = 0;
			double mi_normY = 0;
			double mi_normZ = 0;

			for (auto& node : mv_ptElem->getConectivity()) {
				for (int j = 0; j < 3; ++j) {
					mi_vectX.at(j) += node->getInitPos()[j] * *(mi_pShape + n + m);
					mi_vectY.at(j) += node->getInitPos()[j] * *(mi_pShape + n + m + 1);
				}
				m = m + 2;
			}

			// Norm of vectors in X and Y
			for (int j = 0; j < 3; ++j) {
				mi_normX += mi_vectX.at(j) * mi_vectX.at(j);
				mi_normY += mi_vectY.at(j) * mi_vectY.at(j);
			}

			mi_normX = sqrt(mi_normX);
			mi_normY = sqrt(mi_normY);

			// Unity vector in X and Y
			for (int j = 0; j < 3; ++j) {
				mi_vectX.at(j) = mi_vectX.at(j) / mi_normX;
				mi_vectY.at(j) = mi_vectY.at(j) / mi_normY;
			}

			// Cross product -> VectZ
			mi_vectZ.at(0) = mi_vectX.at(1) * mi_vectY.at(2) - mi_vectX.at(2) * mi_vectY.at(1);
			mi_vectZ.at(1) = mi_vectX.at(2) * mi_vectY.at(0) - mi_vectX.at(0) * mi_vectY.at(2);
			mi_vectZ.at(2) = mi_vectX.at(0) * mi_vectY.at(1) - mi_vectX.at(1) * mi_vectY.at(0);

			mi_normZ = sqrt(mi_vectZ.at(0) * mi_vectZ.at(0) + mi_vectZ.at(1) * mi_vectZ.at(1) + mi_vectZ.at(2) * mi_vectZ.at(2));

			// Unity vector in Z
			for (int j = 0; j < 3; ++j) {
				mi_vectZ.at(j) = mi_vectZ.at(j) / mi_normZ;
			}

			// Correction of vector in X direction with cross product of Y and Z
			mi_vectX.at(0) = mi_vectY.at(1) * mi_vectZ.at(2) - mi_vectY.at(2) * mi_vectZ.at(1);
			mi_vectX.at(1) = mi_vectY.at(2) * mi_vectZ.at(0) - mi_vectY.at(0) * mi_vectZ.at(2);
			mi_vectX.at(2) = mi_vectY.at(0) * mi_vectZ.at(1) - mi_vectY.at(1) * mi_vectZ.at(0);

			mi_normX = sqrt(mi_vectX.at(0) * mi_vectX.at(0) + mi_vectX.at(1) * mi_vectX.at(1) + mi_vectX.at(2) * mi_vectX.at(2));

			// Unity vector in X
			for (int j = 0; j < 3; ++j) {
				mi_vectX.at(j) = mi_vectX.at(j) / mi_normX;
			}

			// Now lets create the rotation matrix
			for (int j = 0; j < 3; ++j) {
				mi_rot.at(0, j) = mi_vectX.at(j);
				mi_rot.at(1, j) = mi_vectY.at(j);
				mi_rot.at(2, j) = mi_vectZ.at(j);
			}

			// Finally evaluates the Jacobian matrix
			for (int j = 0; j < 2; ++j) {
				m = 0;

				for (auto node : mv_ptElem->getConectivity()) {
					// Auxiliary coordinates
					mi_cd.at(0) = mi_rot(0, 0) * node->getInitPos()[0] + mi_rot(0, 1) * node->getInitPos()[1] + mi_rot(0, 2) * node->getInitPos()[2];
					mi_cd.at(1) = mi_rot(1, 0) * node->getInitPos()[0] + mi_rot(1, 1) * node->getInitPos()[1] + mi_rot(1, 2) * node->getInitPos()[2];

					for (int k = 0; k < 2; ++k) {
						mi_F0(j, k) += mi_cd.at(j) * *(mi_pShape + n + m);
						m++;
					}
				}
			}

			// Record the Jacobian Matrix and the Rotation Matrix on every Material (Integration) point
			this->mv_matPoint.emplace_back(std::make_unique<O2P2::Proc::Comp::MaterialPoint>(mi_F0, mi_rot));
		}
	}
}


// ================================================================================================
//
// Implementation of Template Member Function: setIndexing
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshIncl<nDim>::setIndexing()
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

	// When 2D elements are immersed in 3D elements, elemental contribution must be resized
	if (nDim != mv_ptElem->getDIM()) {
		mv_nDof = mv_conect.size() * nDim;
		mv_elHes.resize(mv_nDof);
		mv_elFor.resize(mv_nDof);
	}

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
// Implementation of MeshIncl Member Function: getContribution_SVK
//
// ================================================================================================
template<int nDim>
void O2P2::Proc::Comp::MeshIncl<nDim>::getContribution_SVK(double section)
{
	int mi_ElDim = mv_ptElem->getDIM();

	if (nDim == mi_ElDim) {
		O2P2::Proc::Comp::MeshElem::getContribution_SVK(section);
	}
	else {
		// Must to rotate elements
		getRotatedContribution_SVK(section);
	}
}


// ================================================================================================
//
// Implementation of 3D MeshIncl Member Function: getRotatedContribution_SVK
// Only called when 2D elements are immersed in 3D environment
//
// ================================================================================================
template<>
void O2P2::Proc::Comp::MeshIncl<3>::getRotatedContribution_SVK(double section)
{
	// 2D element in 3D environment
	const int mi_ElDim = 2;
	const int mi_Dim = 3;

	// Pointer to the first Shape Fuctions Derivative (const static)
	auto mi_pShape = mv_ptElem->getShapeDerivative();
	auto mi_pWeight = mv_ptElem->getWeight();

	auto mi_nNodes = mv_ptElem->getNumNodes();
	auto mi_nIP = mv_ptElem->getNumIP();
	auto mi_pMat = mv_ptElem->getMaterial();

	// Elasticity tensor in material description
	M2S2::Dyadic4C E(mi_ElDim);
	E = mi_pMat->getConstitutiveMatrix(mi_ElDim);

	// 2D element constants
	const int mi_nVoigt = 3 * mi_ElDim - 3;
	const int mi_nSize = mi_ElDim * mi_ElDim;

	// Local element matrix and vector. These will be rotated and then stored
	M2S2::MatrixS mi_elHes(mv_nDof);
	std::vector<double> mi_elFor(mv_nDof);

	// Rotation matrix
	M2S2::MatrixX mi_R(mv_nDof, mv_nDof);

	// One material point for each integration point
	for (int i = 0; i < mi_nIP; ++i) {

		// Clear matrices for current integration point
		mi_R.clear();
		mi_elHes.clear();
		memset(&mi_elFor[0], 0., mv_nDof * sizeof(double));

		// Deformation gradient from undeformed to non-dimensional space
		M2S2::Dyadic2N F0(mi_ElDim);
		F0 = mv_matPoint.at(i)->getJacobianMatrix();

		// Rotation matrix associated to the integration point
		M2S2::Dyadic2N mi_rot(mi_Dim);
		mi_rot = mv_matPoint.at(i)->getRotationMatrix();

		// Coordinates in rotated space
		std::array<double, 2> mi_cd{ 0, 0 };

		for (int j = 0; j < mv_nDof; j = j + 3) {
			for (int k = 0; k < 3; ++k) {
				for (int m = 0; m < 3; ++m) {
					mi_R(j + k, j + m) = mi_rot(k, m);
				}
			}
		}

		// Deformation gradient from deformed to non-dimensional space
		M2S2::Dyadic2N F1(mi_ElDim);

		// Some pointer arithmetic is required
		// p and numInt are related to current integration point and weighting numerical integration
		int p = i * mi_nNodes * mi_ElDim;
		double numInt = *(mi_pWeight + i) * this->mv_matPoint.at(i)->getJacobian() * section;

		for (int j = 0; j < mi_ElDim; ++j) {
			int m = 0;

			for (auto& node : mv_conect) {
				// Auxiliary coordinates
				mi_cd.at(0) = mi_rot(0, 0) * node->getTrialPos()[0] + mi_rot(0, 1) * node->getTrialPos()[1] + mi_rot(0, 2) * node->getTrialPos()[2];
				mi_cd.at(1) = mi_rot(1, 0) * node->getTrialPos()[0] + mi_rot(1, 1) * node->getTrialPos()[1] + mi_rot(1, 2) * node->getTrialPos()[2];

				for (int k = 0; k < mi_ElDim; ++k) {
					F1(j, k) += mi_cd.at(j) * *(mi_pShape + p + m);
					m++;
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
				int mi_curDof = r * mi_Dim + d;

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
				mi_elFor.at(mi_curDof) += dUedy * numInt;
			}
		}

		//Hessian += dLdy.transpose() * E * dLdy * numInt;
		for (int j = 0; j < mv_nDof; ++j) {
			for (int k = j; k < mv_nDof; ++k) {
				for (int l = 0; l < mi_ElDim; ++l) {
					for (int m = l; m < mi_ElDim; ++m) {
						mi_elHes.at(j, k) += dLdy.at(m, j) * E.at(m, l) * dLdy.at(l, k) * numInt;
					}
				}
				for (int l = mi_ElDim; l < mi_nVoigt; ++l) {
					for (int m = l; m < mi_nVoigt; ++m) {
						mi_elHes.at(j, k) += 2. * dLdy.at(m, j) * E.at(m, l) * dLdy.at(l, k) * numInt;
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
						if (j * mi_Dim <= k * mi_Dim) mi_elHes.at(j * mi_Dim + m, k * mi_Dim + m) += SM.at(n) * d2Ldy2.at(n) * numInt;
					}

					for (int n = mi_nVoigt; n < mi_nSize; ++n) {
						if (j * mi_Dim <= k * mi_Dim) mi_elHes.at(j * mi_Dim + m, k * mi_Dim + m) += SM.at(n + mi_nVoigt - mi_nSize) * d2Ldy2.at(n) * numInt;
					}
				}
			}
		}

		// Rotation of mi_elHes e mi_elFor
		mv_elHes += (mi_R.transpose() * mi_elHes * mi_R).SaveAsMatrixS();

		for (int m = 0; m < mv_nDof; ++m) {
			for (int n = 0; n < mv_nDof; ++n) {
				mv_elFor.at(m) += mi_R.at(n, m) * mi_elFor.at(n);
			}
		}
	}
}


// ================================================================================================
//
// Implementation of 2D MeshIncl Member Function: getRotatedContribution_SVK
// Only here for completeness. It will never be called, but it is here for completeness.
//
// ================================================================================================
template<>
void O2P2::Proc::Comp::MeshIncl<2>::getRotatedContribution_SVK(double section)
{
	O2P2::Proc::Comp::MeshElem::getContribution_SVK(section);
}


// ================================================================================================
//
// Implementation of MeshIncl Member Function: getSpreadMatrix
//
// ================================================================================================
template<int nDim>
M2S2::MatrixX O2P2::Proc::Comp::MeshIncl<nDim>::getSpreadMatrix()
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
