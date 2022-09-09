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

// Custom Header Files
#include "Common.h"
#include "Domain.h"
#include "Mesh.h"
#include "NonLinearSolver.h"

namespace O2P2 {
	namespace Proc {
		/** @ingroup Processor_Module
		  * @class TimeStepping
		  *
		  * @brief Time step integration schemes.
		  * @details This class manages the time step integration scheme and handles variables for quasi-static / parabolic / hyperbolic / eigenvalue analysis.
		  * The values required for the time step integration scheme are available at LoadStep.
		  *
		  * @todo 1 - Funcionalidade para problemas parabólicos (temperatura).
		  * @todo 2 - Funcionalidade para problemas dinâmicos hiperbólicos (Newmark).
		  * @todo 3 - Funcionalidade para problemas dinâmicos hiperbólicos (HHT alpha).
		  * @todo 4 - Funcionalidade para problemas de autovalor.
		  *
		  * @sa Mesh
		  * @sa LoadStep
		  */
		class TimeStepping {

		public:
			// Default destructor of private / protected pointers.
			virtual ~TimeStepping() = default;

			/** Manages the time stepping integration scheme.
			  * @param theDomain Container with nodal and elements information.
			  * @param theFEModel Container of solution components.
			  * @param theSolver Pointer to the non-linear solver.
			  */
			virtual void runTimeLoop(O2P2::Prep::Domain<2>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) = 0;

			/** Manages the time stepping integration scheme.
			  * @param theDomain Container with nodal and elements information.
			  * @param theFEModel Container of solution components.
			  * @param theSolver Pointer to the non-linear solver.
			  */
			virtual void runTimeLoop(O2P2::Prep::Domain<3>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) = 0;

		protected:
			/** @brief Default construtor. Implemented in derived classes.  */
			TimeStepping() { };
		};


		/** @ingroup TimeStep
		  * @class TimeStep_QsiStatic
		  * @brief Time step integration schemes for quasi-static problems.
		  *
		  * @todo Representação do esquema de integração quase-estático.
		  */
		class TimeStep_QsiStatic : public TimeStepping
		{
		public:
			/** Constructor for quasi-static time stepping integration scheme. */
			TimeStep_QsiStatic() {};

			// Work-around to implement Quasi-Static time stepping integration for 2D and 3D problems.
			void runTimeLoop(O2P2::Prep::Domain<2>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) override {
				PROFILE_FUNCTION();
				this->runTimeLooping(theDomain, theFEModel, theSolver);
			}

			// Work-around to implement Quasi-Static time stepping integration for 2D and 3D problems.
			void runTimeLoop(O2P2::Prep::Domain<3>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) override {
				PROFILE_FUNCTION();
				this->runTimeLooping(theDomain, theFEModel, theSolver);
			}

		private:
			// Actual implementation of Quasi-Static time stepping integration. 2D AND 3D.
			template<int nDim>
			void runTimeLooping(O2P2::Prep::Domain<nDim>* theDomain, O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver);
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace
