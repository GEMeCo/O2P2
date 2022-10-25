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

// C++ standard libraries
#include <array>				// Standard array manipulation (fixed size)
#include <vector>				// Standard vector manipulation
#include <iomanip>				// Required by ios (stream)

#include "Common.h"

namespace O2P2 {
	namespace Proc {
		namespace Comp {
			/** @ingroup Processor_Module
			  * @class MeshNode
			  *
			  * @brief Mesh node, a Solution component.
			  * @details This base class only provides interface for derived classes, which contains current trial and commit information on nodes.
			  */
			class MeshNode
			{
			private:
				MeshNode() = delete;

			protected:
				/** Constructor for mesh node objects.
				  * @param index Node DOF index (the first one).
				  * @param Tp Temperature.
				  */
				MeshNode(const size_t& index, const double Tp = 0.) : m_DofIndex(index), m_Tp(Tp) {	}

			public:
				// Default destructor of private / protected pointers.
				virtual ~MeshNode() = default;

				/** Overloading operator << to stream the node current position. */
				friend std::ostream& operator<<(std::ostream& stream, MeshNode const& node) {
					stream << node.print();
					return stream;
				}

				/** @return string with current position for printing. */
				virtual const std::string print() const = 0;

				/** Retrieve current commited position.
				  * @return pointer to nodal commited position (with size nDim).
				  */
				virtual const double* getCurPos() = 0;

				/** Retrieve current trial position.
				  * @return pointer to nodal trial position (with size nDim).
				  */
				virtual const double* getTrialPos() = 0;

				/** Retrieve current nodal temperature.
				  * @return double with nodal temperature.
				  */
				const double& getCurTemp() { return m_Tp; }

				/** Retrieve current trial nodal temperature.
				  * @return double with nodal trial temperature.
				  */
				const double& getTrialTemp() { return m_Tp; }

				/** Update trial position / velocity / acceleration.
				  * @param dPos Increase in the trial position and new values for velociy and acceleration (in a single vector, and only if applied).
				  * @warning The size of dPos is not validated anywhere. Undoubtedly will lead to access error, if you don't know what you are doing. 
				  */
				virtual void updateTrial(const double dPos[]) = 0;

				/** Update trial temperature.
				  * @param dTp New value for trial temperature.
				  */
				void updateTrial(const double dTp) { m_Tp = dTp; }

				/** Commit trial to current. */
				virtual void setCurrent() = 0;

				/** @return the number of DOF per node.
				  */
				virtual const int getNumDOF() = 0;

				/** @brief Node DOF index. */
				size_t m_DofIndex = 0;

			protected:
				/** @brief Room / Commit temperature. */
				double m_Tp;
			};


			/**
			  * @class MeshNode_MQS
			  *
			  * @brief Mesh node, a Solution component, for mechanical quasi-static problems.
			  * @details The solution node component for mechanical quasi-static problems includes current position (trial and commit).
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshNode_MQS : public MeshNode
			{
			private:
				MeshNode_MQS() = delete;

			public:
				/** Constructor for node objects of quasi-static mechanical problems.
				  * @param index Node DOF index (the first one).
				  * @param InitPos Initial position of the corresponding geometry node. Used to initialize trial position.
				  * @param Tp Initial temperature.
				  */
				explicit MeshNode_MQS(const size_t& index, const std::array<double, nDim>& InitPos, const double Tp = 0.)
					: MeshNode(index, Tp), v_Ytrial(InitPos), v_Ycommit(InitPos) { }

				// Default destructor of private / protected pointers.
				~MeshNode_MQS() = default;

				// Return a string with current position for printing.
				const std::string print() const override {
					std::string str;
					for (auto& pos : v_Ycommit) {
						str.append(std::to_string(pos));
						str.append(" ");
					}
					return str;
				}

				// Function to retrieve current nodal coordinates.
				const double* getCurPos() override { return this->v_Ycommit.data(); }

				// Function to retrieve current trial nodal coordinates.
				const double* getTrialPos() override { return this->v_Ytrial.data(); }

				// Function to increase trial position.
				void updateTrial(const double dPos[]) override {
					for (int i = 0; i < nDim; ++i) {
						this->v_Ytrial[i] += dPos[i];
					}
				}

				// Function to commit trial position to current.
				void setCurrent() override { this->v_Ycommit = this->v_Ytrial; }

				/** @return the number of DOF per node.
				  */
				const int getNumDOF() override { return nDim; }

			protected:
				/** @brief Trial current position. */
				std::array<double, nDim> v_Ytrial;

				/** @brief Commit current position. */
				std::array<double, nDim> v_Ycommit;
			};


			/**
			  * @class MeshNode_MD
			  *
			  * @brief Mesh node, a Solution component, for mechanical dynamic problems.
			  * @details The solution node component for mechanical dynamic problems includes current position (trial and commit), and previous and current (timestep) velocity and accelarion.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshNode_MD : public MeshNode_MQS<nDim>
			{
			private:
				MeshNode_MD() = delete;

			public:
				/** Constructor for node objects of dynamic mechanical problems.
				  * @param index Node DOF index (only the first one).
				  * @param InitPos Initial position of the corresponding geometry node. Used to initialize trial position.
				  * @param Tp Initial temperature.
				  */
				explicit MeshNode_MD(const size_t& index, const std::array<double, nDim>& InitPos, const double Tp = 0.)
					: MeshNode_MQS<nDim>(index, InitPos, Tp) {
					v_Vp.fill(0);
					v_Vc.fill(0);
					v_Ap.fill(0);
					v_Ac.fill(0);
				}

				// Default destructor of private / protected pointers.
				~MeshNode_MD() = default;

				// Function to retrieve current nodal coordinates.
				const double* getCurPos() override { return this->v_Ycommit.data(); }

				// Function to retrieve current trial nodal coordinates.
				const double* getTrialPos() override { return this->v_Ytrial.data(); }

				/** Retrieve previous velocity.
				  * @return pointer to nodal previous velocity (with size nDim).
				  */
				const double* getPrevVel() { return v_Vp.data(); }

				/** Retrieve current velocity.
				  * @return pointer to nodal current velocity (with size nDim).
				  */
				const double* getCurVel() { return v_Vc.data(); }

				/** Retrieve previous acceleration.
				  * @return pointer to nodal previous acceleration (with size nDim).
				  */
				const double* getPrevAcc() { return v_Ap.data(); }

				/** Retrieve current acceleration.
				  * @return pointer to nodal current acceleration (with size nDim).
				  */
				const double* getCurAcc() { return v_Ac.data(); }

				/** Initiates velocity.
				  * @param Dir Direction of the imposed velocity.
				  * @param value Imposed velocity.
				  */
				void setVel(int Dir, double value) { 
					this->v_Vp[Dir] = value;
					this->v_Vc[Dir] = value;
				}

				/** Initiates acceleration.
				  * @param Dir Direction of the imposed acceleration.
				  * @param value Imposed acceleration.
				  */
				void setAcc(int Dir, double value) {
					this->v_Ap[Dir] = value;
					this->v_Ac[Dir] = value;
				}


				// Update trial position, velocity and acceleration.
				void updateTrial(const double dPos[]) override {
					for (int i = 0; i < nDim; ++i) {
						this->v_Ytrial[i] = dPos[i];
						this->v_Vc[i] = dPos[nDim + i];
						this->v_Ac[i] = dPos[2 * nDim + i];
					}
				}

				// Function to commit trial position, velocity and acceleration to current.
				void setCurrent() override { 
					this->v_Ycommit = this->v_Ytrial;
					this->v_Vp = this->v_Vc;
					this->v_Ap = this->v_Ac;
				}

				/** @return the number of DOF per node.
				  */
				const int getNumDOF() override { return nDim; }

			protected:
				/** @brief Previous velocity. */
				std::array<double, nDim> v_Vp;

				/** @brief Current velocity. */
				std::array<double, nDim> v_Vc;

				/** @brief Previous acceleration. */
				std::array<double, nDim> v_Ap;

				/** @brief Current acceleration. */
				std::array<double, nDim> v_Ac;
			};

		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
