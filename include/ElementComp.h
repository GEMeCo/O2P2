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
#pragma once

// C++ standard libraries
#include <vector>		// required by std::vector
#include <memory>		// required by std::shared_pointer
#include <utility>		// required by std::make_pair

// Custom libraries
#include "Element.h"
#include "MaterialPoint.h"

// Eigen libraries
#include <Eigen/Dense>
#include <Eigen/Geometry>	// required by .cross (Cross product)


/** @ingroup Processor_Module
  * @class ElemComp
  *
  * @brief Base class for elements and inclusion elements routines and local matrices.
  * @details Evaluates local matrices and handles material points elements, one for each integration point.
  * It should also evaluate the internal stresses.
  * 
  * @note For each Geometry Element, there is only one Element Component, according to material and section types
  *
  * @todo 1 - ElemComponent deve estar relacionada ao material (MaterialType). A implementação atual se refere apenas à SVK isotrópico.
  * @todo 2 - Seria interessante implementar uma versão teste com menos temporários e mais funções, para auxiliar novas derivações.
  * Por exemplo, em cada ponto de integração, calcula-se F1, F, C, L, S, dLdy, dUedy, e d2Ldy2. 
  * Uma depende da outra, e tem notações distintas (matricial e Voigt). Poderia registrar só Voigt e retornar conforme solicitado.
  * Será que o overhead seria muito grande? E o uso de memória, mantido (afinal está usando temporárias)?
  * @todo 3 - getConstitutiveMatrix não depende da temperatura?
  * @todo 4 - Não deveria estar em função da seção transversal em algum lugar (2D)?
  */
class ElemComp
{
public:
	/** @brief Basic destructor. */
	virtual ~ElemComp() = default;

	/** Prepare element contribution.
	  * @param FInt Element contribution to the internal force.
	  * @param Hessian Element contribution to the hessian matrix.
	  */
	virtual void getContribution(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian) = 0;

protected:
	/** @brief Basic constructor. */ 
	ElemComp() { m_nDof = 0; };

	/** @return a pointer to the element indexing, after casting.
	  */
	template<int nDim>
	std::vector<std::shared_ptr<Node<nDim>>> getConectivity() { return nullptr; };

	/** @return a pointer to a node.
	  * @param index Element convectivity container index.
	  */
	template<int nDim>
	Node<nDim>* getConectivity(const int& index) { return nullptr; };

protected:
	/** @brief Pointer to material point that holds information in the integration point. */
	std::vector<std::unique_ptr<MaterialPoint>> m_MatPoint;

public:
	/** @brief Number of DOF for current element. */
	int m_nDof;

	/** @brief Vector with element DOF index. */
	std::vector<size_t> m_ElemDofIndex;
};


/** 
  * @class ElemComponent
  *
  * @brief Class for elements and inclusion elements routines and local matrices for SVK material model.
  * @details Evaluates local matrices and handles material points elements, one for each integration point.
  * It should also evaluate the internal stresses (2nd Piola Kirchhoff).
  *
  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
  */
template<int nDim>
class ElemComponent : public ElemComp
{
private:
	ElemComponent() = delete;

public:
	/** Constructor for mechanical analysis elements, for SVK material model.
	  * @param pElmt Pointer to geometry element.
	  */
	explicit ElemComponent(std::shared_ptr<BaseElement> pElmt) : ElemComp() {
		m_pElem = pElmt;
		this->m_MatPoint.reserve(pElmt->getNumIP());

		this->setMaterialPoint();
	};

	/** @brief Basic destructor. */
	virtual ~ElemComponent() = default;

	// Prepare element contribution.
	void getContribution(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian) override;

protected:
	/** @brief Populates the material point vector. */
	void setMaterialPoint();

	/** @return the constitutive matrix (in Voigt notation).
	  */
	Eigen::MatrixXd getConstitutiveMatrix();

	/** Prepare element contribution fo SVK_ISO material.
	  * @param FInt Element contribution to the internal force.
	  * @param Hessian Element contribution to the hessian matrix.
	  */
	void getContribution_SVK_ISO(Eigen::VectorXd& FInt, Eigen::MatrixXd& Hessian);

	/** @return a pointer to the element indexing, after casting.
	  */
	std::vector<std::shared_ptr<Node<nDim>>> getConectivity() {
		Element<nDim>* pElem = static_cast<Element<nDim>*> (m_pElem.get());
		return pElem->getConectivity();
	};

	/** @return a pointer to a node.
	  * @param index Element convectivity container index.
	  */
	Node<nDim>* getConectivity(const int& index) {
		Element<nDim>* pElem = static_cast<Element<nDim>*> (m_pElem.get());
		return pElem->getConectivity(index);
	};

protected:
	/** @brief Pointer to domain element. */
	std::shared_ptr<BaseElement> m_pElem;
};
