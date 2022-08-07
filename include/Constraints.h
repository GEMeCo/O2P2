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

/** @ingroup Processor_Module
  * @class Constraint
  *
  * @brief DOF constraint
  * @details Dirichlet and Neumann boundary conditions
  */
class Constraint
{
private:
    Constraint() = delete;

public:
    /** Constructor for both Dirichlet and Neumann boundary conditions.
      * @param Dof Degree of Freedom with imposed boundary condition.
      * @param Value Value of the boundary condition.
      * @param Var Time behavior, for a quadratic polinomial (var[0] + var[1].t + var[2].t^2).
      */
    Constraint(const size_t& Dof, const double& Value, const double Var[3]) {
        m_Dof = Dof;
        m_Value = Value;
        m_TimeVar[0] = Var[0];
        m_TimeVar[1] = Var[1];
        m_TimeVar[2] = Var[2];
    };

    /** Destructor of private / protected pointers. */
    ~Constraint() = default;

    /** @brief Degree of Freedom with boundary condition */
    size_t m_Dof{ 0 };

    /** @brief Value of the boundary condition, if apply */
    double m_Value{ 0 };

    /** @brief Time function: _Value * (TimeVar(1) + TimeVar(2)*dT + TimeVar(3)*dT^2 */
    double m_TimeVar[3]{ 0 }; 
};
