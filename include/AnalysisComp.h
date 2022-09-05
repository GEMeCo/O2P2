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

// Eigen libraries
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

// Custom Header Files
#include "Common.h"
#include "LoadStep.h"

#include "Domain.h"
#include "ElementComp.h"
#include "PostProcess.h"

/** @ingroup Processor_Module
  * @class AnalysisComp
  *
  * @brief Container of solution components.
  * @details This class holds vectors of solution and analysis components: Load Steps, DOF, System of equation, and such.
  *
  * @todo 1 - Funcionalidade para problemas térmicos.
  * @todo 2 - Geração de elemComp não está associada ao problema (térmico, mecânico, etc), afinal só há uma implementação de elemComp.
  * @todo 3 - Não seria mais interessante levar os container de elementos e nós para classe base (atualmente na derivada _MEC)?
  * Desta forma, as classes derivadas conteriam apenas implementações específicas, sem novos dados...
  */
class AnalysisComp
{
private:
	AnalysisComp() = delete;

protected:
	/** Constructor for container of solution components.
	  * @param vOut Reference to post-process container.
	  */
	AnalysisComp(PostProcess* vOut) { m_PostPt = vOut; };

public:
	/** Default destructor. */
	virtual ~AnalysisComp() = default;

	/** Add a new load step.
	  * Any change in boundary conditions, time step and such, requires a new load step. All values must be included again.
	  * @param NumSteps Number of time steps in the current load step.
	  * @param TimeStep Variation in time in the current load step.
	  * @param NumDBC Number of applied Dirichlet Boundary Conditions.
	  * @param NumNBC Number of applied Neumann Boundary Conditions.
	  */
	void addLoadStep(const int& NumSteps, const double& TimeStep, const size_t& NumDBC, const size_t& NumNBC) {
		m_LoadStep.push_back(std::make_unique<LoadStep>(NumSteps, TimeStep, NumDBC, NumNBC));
	};

	/** Add a Boundary Condition of Dirichlet type to a Load Step.
	  * @param nLS Number of the load step to receive the BC.
	  * @param Dof Degree of Freedom with imposed boundary condition.
	  * @param Value Value of the boundary condition.
	  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
	  */
	void addDirichletBC(const int& nLS, const size_t& Dof, const double& Value, const double Var[]) {
		m_LoadStep.at(nLS)->addDirichletBC(Dof, Value, Var);
	};

	/** Add a Boundary Condition of Neumann type to a Load Step.
	  * @param nLS Number of the load step to receive the BC.
	  * @param Dof Degree of Freedom with imposed boundary condition.
	  * @param Value Value of the boundary condition.
	  * @param Var Time behavior, for a quadratic polinomial(var[0] + var[1].t + var[2].t ^ 2).
	  */
	void addNeumannBC(const int& nLS, const size_t& Dof, const double& Value, const double Var[]) {
		m_LoadStep.at(nLS)->addNeumannBC(Dof, Value, Var);
	};

	/** Add new dof to the system. Used to calculate the total number of DOF.
	  * @param nDOF Number of DOF to be added to the system of equation.
	  */
	void addDOF(int nDOF) { m_TotalDof += nDOF; };

	/** Assemble the system of equation, made by the Hessian matrix and a right hand side vector.
	  * @param Hessian Square matrix of second-order partial derivatives of a scalar-valued function or scalar field.
	  * @param RHS Vector of independent terms, the right hand side, for system: Hessian.LHS = RHS.
	  */
	virtual void assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS) = 0;

	/** Initiates the auxiliary vector used to impose boundary conditions to the system of equations.
	  * @param loadStep Current load step under process.
	  */
	void initDirichletBC(const int& loadStep);

	/** Sets the current timestep.
	  * @param timeStep Current timestep under process.
	  */
	void setTimeStep(const int& timeStep) { this->m_curTimeStep = timeStep; };

	/** Update trial solution.
	  * @param LHS Left hand side vector with current trial solution.
	  */
	virtual void setTrial(Eigen::VectorXd& LHS) = 0;

	/** Once solution is achieved, commit it. */
	virtual void setCommit() = 0;

	/** @return a pointer to the current load step.
	  */
	LoadStep* getLoadStep() { return m_LoadStep.at(m_curLoadStep).get(); };

	/** Retrieve a pointer to a load step.
	  * @return A pointer to load step iLS.
	  * @param iLS Number of Load Step to be retrieved.
	  */
	LoadStep* getLoadStep(int iLS) { return m_LoadStep.at(iLS).get(); };

	/** @return the number of the degrees of freedom (size of the problem). 
	  */
	const size_t getNumDof() { return m_TotalDof; };

	/** Impose Neumann Boundary Conditions to the vector of independent terms (external forces).
	  * @param RHS Right hand side vector to impose the Neumann boundary conditions.
	  */
	void imposeNeumannBC(Eigen::VectorXd& RHS);

public:
	/** @brief Vector with Dirichlet boundary conditions - if there is a BC imposed, then it is equal to 0. Otherwise, it is 1. */
	std::vector<int> m_BCIndex;

	/** @brief Current time of anaysis. */
	double m_currentTime{ 0. };

	/** @brief Container of element components associated to the analysis. */
	std::vector<std::unique_ptr<ElemComp>> m_ElemComp;

protected:
	/** @brief Current load step under process. */
	int m_curLoadStep{ 0 };

	/** @brief Current time step of current load step. */
	int m_curTimeStep{ 0 };

	/** @brief Total number of DOF */
	size_t m_TotalDof{ 0 };

	/** @brief Container of boundary conditions for load steps. */
	std::vector<std::unique_ptr<LoadStep>> m_LoadStep;

	/** @brief Container of solution for post-process. */
	PostProcess* m_PostPt;
};


/** 
  * @class AnalysisComp_Mec
  *
  * @brief Container of solution components for mechanical problems.
  * @details This class holds vectors of solution and analysis components: Load Steps, DOF, System of equation, and such.
  * 
  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
  */
template<int nDim>
class AnalysisComp_Mec : public AnalysisComp
{
private:
	AnalysisComp_Mec() = delete;

public:
	/** Constructor for container of solution components.
	  * @param theDomain Container with nodaland elements information.
	  * @param vOut  Reference to post-process container.
	  */
	explicit AnalysisComp_Mec(Domain<nDim>* theDomain, PostProcess* vOut) : AnalysisComp(vOut) {
		m_NodePt = &theDomain->getNode();
		m_ElemComp.reserve(theDomain->getElem().size());

		// Generates element components for each domain element
		for (std::shared_ptr<Element<nDim>> elem : theDomain->getElem()) {
			m_ElemComp.emplace_back(std::make_unique<ElemComponent<nDim>>(elem));
		}
	};

	/** @brief Default destructor. */
	~AnalysisComp_Mec() = default;

	// Assemble the system of equation, made by the Hessian matrix and the right hand side vector.
	void assembleSOE(Eigen::SparseMatrix<double>& Hessian, Eigen::VectorXd& RHS) override;

	// Update trial solution.
	void setTrial(Eigen::VectorXd& LHS) override;

	// Update commit solution.
	void setCommit() override;

public:
	/** @brief Pointer to the node container. The pointer is required to directly access and update its contents. */
	std::vector<std::shared_ptr<Node<nDim>>>* m_NodePt;
};