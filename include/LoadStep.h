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
#include <memory>				// required by std::shared_ptr
#include <vector>				// required by std::vector

// Custom Header Files
#include "Constraints.h"

namespace O2P2 {
    namespace Proc {
        namespace Comp {
            /** @ingroup Processor_Module
              * @class LoadStep
              *
              * @brief Container of boundary conditions.
              * @details This class holds Neumann and Dirichlet boundary conditions. It also keeps information about the current load step.
              */
            class LoadStep
            {
            private:
                LoadStep() = delete;

            public:
                /** Constructor of the container of boundary conditions, that holds information about the current load step.
                  * @param NumSteps Number of time steps in the current load step.
                  * @param TimeStep Variation in time in the current load step.
                  * @param NumDBC Number of applied Dirichlet Boundary Conditions.
                  * @param NumNBC Number of applied Neumann Boundary Conditions.
                  */
                explicit LoadStep(const int& NumSteps, const double& TimeStep, const size_t& NumDBC, const size_t& NumNBC) :
                    m_NumSteps(NumSteps), m_TimeStep(TimeStep) {
                    v_DirichletBC.reserve(NumDBC);
                    v_NeumannBC.reserve(NumNBC);
                }

                // Default destructor of private / protected pointers.
                ~LoadStep() = default;

                /** Adds a Dirichlet Boundary Condition
                  * @param Dof Degree of Freedom with imposed boundary condition.
                  * @param Value Value of the boundary condition
                  * @param Var Time behavior, for a quadratic polinomial (var[0] + var[1].t + var[2].t^2)
                  */
                void addDirichletBC(const size_t& Dof, const double& Value, const double Var[]) {
                    v_DirichletBC.push_back(std::make_unique<O2P2::Proc::Comp::Constraint>(Dof, Value, Var));
                }

                /** Adds a Neumann Boundary Condition on the flux variable
                  * @param Dof Degree of Freedom with imposed boundary condition.
                  * @param Value Value of the boundary condition
                  * @param Var Time behavior, for a quadratic polinomial (var[0] + var[1].t + var[2].t^2)
                  */
                void addNeumannBC(const size_t& Dof, const double& Value, const double Var[]) {
                    v_NeumannBC.push_back(std::make_unique<O2P2::Proc::Comp::Constraint>(Dof, Value, Var));
                }

            public:
                /** @brief Number of time steps in the current load step */
                int m_NumSteps;

                /** @brief Time step */
                double m_TimeStep;

                /** @brief Container of Dirichlet Boundary Conditions (main variable, i.e. temperature) */
                std::vector<std::unique_ptr<O2P2::Proc::Comp::Constraint>> v_DirichletBC;

                /** @brief Container of Neumann Boundary Conditions (flux variable, i.e. applied force) */
                std::vector<std::unique_ptr<O2P2::Proc::Comp::Constraint>> v_NeumannBC;
            };
        } // End of Comp Namespace
    } // End of Proc Namespace
} // End of O2P2 Namespace
