// ================================================================================================
// 
// This file is part of O2P2, an object oriented environment for the positional FEM
//
// Copyright(C) 2022 Rogerio Carrazedo - All Rights Reserved.
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
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
		  * @sa Mesh
		  * @sa LoadStep
		  */
		class TimeStepping {

		public:
			// Default destructor of private / protected pointers.
			virtual ~TimeStepping() = default;

			/** Manages the time stepping integration scheme.
			  * @param theFEModel Container of solution components.
			  * @param theSolver Pointer to the non-linear solver.
			  */
			virtual void runTimeLoop(O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) = 0;

		protected:
			/** @brief Default construtor. Implemented in derived classes.  */
			TimeStepping() { };
		};


		/** @ingroup TimeStep
		  * @class TimeStep_QsiStatic
		  * @brief Time step integration schemes for quasi-static problems.
		  *
		  */
		class TimeStep_QsiStatic : public TimeStepping
		{
		public:
			/** Constructor for quasi-static time stepping integration scheme. */
			TimeStep_QsiStatic() {};

			// Work-around to implement Quasi-Static time stepping integration for 2D and 3D problems.
			void runTimeLoop(O2P2::Proc::Mesh* theFEModel, O2P2::Proc::NonLinearSolver* theSolver) override;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace
