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
		class TimeStepping
		{
		public:
			// Default destructor of private / protected pointers.
			virtual ~TimeStepping() = default;

			/** Manages the time stepping integration scheme.
			  * @param theMesh Container of solution components.
			  * @param theSolver Pointer to the non-linear solver.
			  */
			virtual void runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver) = 0;

			/** Set time stepping parameters
			  * @param alfa Alfa for first order transient analysis (thermal or coupled).
			  * @param beta Beta for second order transient analysis (mechanical).
			  * @param gamma Gamma for second order transient analysis (mechanical).
			  */
			virtual void SetParameters(const double& alfa, const double& beta, const double& gamma) = 0;

		protected:
			/** @brief Default construtor. Implemented in derived classes. */
			TimeStepping() = default;

			// TimeStepping is unique. Copy constructor will be deleted.
			TimeStepping(const TimeStepping& other) = delete;

			// TimeStepping is unique. Move constructor will be deleted.
			TimeStepping(TimeStepping&& other) = delete;
		};


		/** @ingroup TimeStep
		  * @class TimeStep_QsiStatic
		  * @brief Time step integration schemes for quasi-static problems.
		  */
		class TimeStep_QsiStatic : public TimeStepping
		{	
		public:
			/** Constructor for quasi-static time stepping integration scheme. */
			explicit TimeStep_QsiStatic() = default;

			// Default destructor of private / protected pointers.
			~TimeStep_QsiStatic() = default;

			// Non used parameters.
			void SetParameters(const double& alfa, const double& beta, const double& gamma) override {};

			// Implement Quasi-Static time stepping integration for 2D and 3D problems.
			void runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver) override;

		private:
			// TimeStep_QsiStatic is unique. Copy constructor will be deleted.
			TimeStep_QsiStatic(const TimeStep_QsiStatic& other) = delete;

			// TimeStep_QsiStatic is unique. Move constructor will be deleted.
			TimeStep_QsiStatic(TimeStep_QsiStatic&& other) = delete;
		};


		/** @ingroup TimeStep
		  * @class TimeStep_1stEul
		  * @brief Time step integration schemes for dynamic problems, employing Euler time step integration method for first order transient analysis.
		  */
		class TimeStep_1stEul : public TimeStepping
		{
		public:
			/** Constructor for first order time stepping integration scheme (flux), using Euler scheme. */
			explicit TimeStep_1stEul() = default;

			// Default destructor of private / protected pointers.
			~TimeStep_1stEul() = default;

			// Sets time step integration parameter
			void SetParameters(const double& alfa, const double& beta, const double& gamma) override {
				mv_alfa = alfa;
			};

			// Implement Euler time stepping integration for 2D and 3D problems.
			void runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver) override { 
					throw std::invalid_argument("\n\nUndefined type of analysis\nCheck input file\n\n\n");
			};

		private:
			double mv_alfa = 0.5;		// associated to the flux

			// TimeStep_1stEul is unique. Copy constructor will be deleted.
			TimeStep_1stEul(const TimeStep_1stEul& other) = delete;

			// TimeStep_1stEul is unique. Move constructor will be deleted.
			TimeStep_1stEul(TimeStep_1stEul&& other) = delete;
		};


		/** @ingroup TimeStep
		  * @class TimeStep_2ndNew
		  * @brief Time step integration scheme for dynamic problems, employing Newmark-beta time step integration method for second order transient analysis.
		  */
		class TimeStep_2ndNew : public TimeStepping
		{
		public:
			/** Constructor for second order time stepping integration scheme (velocity + acceleration), using Newmark-beta method. */
			explicit TimeStep_2ndNew() = default;

			// Default destructor of private / protected pointers.
			~TimeStep_2ndNew() = default;

			// Sets the Newmark time step integration parameters
			void SetParameters(const double& alfa, const double& beta, const double& gamma) override {
				// mv_alfa is not used in Newmark-beta
				mv_beta = beta;
				mv_gamma = gamma;
			};

			// Implement Newmark time stepping integration for 2D and 3D problems.
			void runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver) override;

		private:
			// Newmark parameters
			double mv_beta = 0.5;		// associated to the displacement
			double mv_gamma = 0.25;		// associated to the velocity

			// TimeStep_2ndNew is unique. Copy constructor will be deleted.
			TimeStep_2ndNew(const TimeStep_2ndNew& other) = delete;

			// TimeStep_2ndNew is unique. Move constructor will be deleted.
			TimeStep_2ndNew(TimeStep_2ndNew&& other) = delete;
		};


		/** @ingroup TimeStep
		  * @class TimeStep_Eigen
		  * @brief Time step integration scheme for eigenvalue problems.
		  */
		class TimeStep_Eigen : public TimeStepping
		{
		public:
			/** Constructor for quasi-static time stepping integration scheme. */
			/** Constructor for first order time stepping integration scheme (flux), using Euler scheme. */
			explicit TimeStep_Eigen() = default;

			// Default destructor of private / protected pointers.
			~TimeStep_Eigen() = default;

			// Sets time step integration parameter
			void SetParameters(const double& alfa, const double& beta, const double& gamma) override { };

			// Implement Eigenvalue integration scheme for 2D and 3D problems.
			void runTimeLoop(O2P2::Proc::Mesh* theMesh, O2P2::Proc::NonLinearSolver* theSolver) override {
				throw std::invalid_argument("\n\nUndefined type of analysis\nCheck input file\n\n\n");
			};

		private:
			// TimeStep_Eigen is unique. Copy constructor will be deleted.
			TimeStep_Eigen(const TimeStep_Eigen& other) = delete;

			// TimeStep_Eigen is unique. Move constructor will be deleted.
			TimeStep_Eigen(TimeStep_Eigen&& other) = delete;
		};
	} // End of Proc Namespace
} // End of O2P2 Namespace
