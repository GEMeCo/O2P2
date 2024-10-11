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
				  * @param Tp Temperature.
				  */
				explicit MeshNode(const double Tp = 0.) : mv_Tp(Tp) { }

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
				const double& getCurTemp() { return mv_Tp; }

				/** Retrieve current trial nodal temperature.
				  * @return double with nodal trial temperature.
				  */
				const double& getTrialTemp() { return mv_Tp; }

				/** Update trial position / velocity / acceleration.
				  * @param dPos Increase in the trial position and new values for velociy and acceleration (in a single vector, and only if applied).
				  * @warning The size of dPos is not validated anywhere. Undoubtedly will lead to access error, if you don't know what you are doing.
				  */
				virtual void updateTrial(const double dPos[]) = 0;

				/** Update trial temperature.
				  * @param dTp New value for trial temperature.
				  */
				void updateTrial(const double dTp) { mv_Tp = dTp; }

				/** Commit trial to current. */
				virtual void setCurrent() = 0;

				/** @return the number of DOF per node.
				  */
				virtual const int getNumDOF() = 0;

				/** @brief Node DOF list. */
				std::vector<size_t> mv_dofList;

			protected:
				/** @brief Room / Commit temperature. */
				double mv_Tp;
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
				  * @param initPos Initial position of the corresponding geometry node. Used to initialize trial position.
				  * @param Tp Initial temperature.
				  */
				explicit MeshNode_MQS(const size_t& index, const double* initPos, const double Tp = 0.)
					: MeshNode(Tp) {
					for (int i = 0; i < nDim; i++) {
						this->mv_dofList.push_back(index + i);		// WARNING: Here we are assuming that the number of a node dof is fixed and equal to nDim
						mv_Ytrial.at(i) = initPos[i];
						mv_Ycommit.at(i) = initPos[i];
					}
				}

				// Default destructor of private / protected pointers.
				~MeshNode_MQS() = default;

				// Return a string with current position for printing.
				const std::string print() const override {
					std::string str;
					for (auto& pos : mv_Ycommit) {
						str.append(std::to_string(pos));
						str.append(" ");
					}
					return str;
				}

				// Function to retrieve current nodal coordinates.
				const double* getCurPos() override { return this->mv_Ycommit.data(); }

				// Function to retrieve current trial nodal coordinates.
				const double* getTrialPos() override { return this->mv_Ytrial.data(); }

				// Function to increase trial position.
				void updateTrial(const double dPos[]) override {
					for (int i = 0; i < nDim; ++i) {
						this->mv_Ytrial.at(i) += dPos[i];
					}
				}

				// Function to commit trial position to current.
				void setCurrent() override { this->mv_Ycommit = this->mv_Ytrial; }

				/** @return the number of DOF per node.
				  */
				const int getNumDOF() override { return nDim; }

			protected:
				/** @brief Trial current position. */
				std::array<double, nDim> mv_Ytrial;

				/** @brief Commit current position. */
				std::array<double, nDim> mv_Ycommit;
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
				  * @param initPos Initial position of the corresponding geometry node. Used to initialize trial position.
				  * @param Tp Initial temperature.
				  */
				explicit MeshNode_MD(const size_t& index, const double* initPos, const double Tp = 0.)
					: MeshNode_MQS<nDim>(index, initPos, Tp) {
					this->mv_Vp.fill(0.);
					this->mv_Vc.fill(0.);
					this->mv_Ap.fill(0.);
					this->mv_Ac.fill(0.);
				}

				// Default destructor of private / protected pointers.
				~MeshNode_MD() = default;

				// Function to retrieve current nodal coordinates.
				const double* getCurPos() override { return this->mv_Ycommit.data(); }

				// Function to retrieve current trial nodal coordinates.
				const double* getTrialPos() override { return this->mv_Ytrial.data(); }

				/** Retrieve previous velocity.
				  * @return pointer to nodal previous velocity (with size nDim).
				  */
				const double* getPrevVel() { return this->mv_Vp.data(); }

				/** Retrieve current velocity.
				  * @return pointer to nodal current velocity (with size nDim).
				  */
				const double* getCurVel() { return this->mv_Vc.data(); }

				/** Retrieve previous acceleration.
				  * @return pointer to nodal previous acceleration (with size nDim).
				  */
				const double* getPrevAcc() { return this->mv_Ap.data(); }

				/** Retrieve current acceleration.
				  * @return pointer to nodal current acceleration (with size nDim).
				  */
				const double* getCurAcc() { return this->mv_Ac.data(); }

				/** Initiates velocity.
				  * @param Dir Direction of the imposed velocity.
				  * @param value Imposed velocity.
				  */
				void setVel(int Dir, double value) {
					this->mv_Vp.at(Dir) = value;
					this->mv_Vc.at(Dir) = value;
				}

				/** Initiates velocity.
				  * @param value Imposed velocity.
				  */
				void setVel(const std::array<double, nDim> value) {
					for (int i = 0; i < nDim; ++i) {
						this->mv_Vp.at(i) = value.at(i);
						this->mv_Vc.at(i) = value.at(i);
					}
				}

				/** Sets previous acceleration.
				  * @param value Imposed acceleration.
				  */
				void setPrevAcc(const std::array<double, nDim> value) {
					for (int i = 0; i < nDim; ++i) {
						mv_Ap.at(i) = value.at(i);
					}
				}

				/** Sets current acceleration.
				  * @param value Imposed acceleration.
				  */
				void setCurAcc(const std::array<double, nDim> value) {
					for (int i = 0; i < nDim; ++i) {
						mv_Ac.at(i) = value.at(i);
					}
				}

				// Update trial position, velocity and acceleration.
				void updateTrial(const double dPos[]) override {
					for (int i = 0; i < nDim; ++i) {
						this->mv_Ytrial.at(i) += dPos[i];
						this->mv_Vc.at(i) = dPos[nDim + i];
						this->mv_Ac.at(i) = dPos[2 * nDim + i];
					}
				}

				// Function to commit trial position, velocity and acceleration to current.
				void setCurrent() override { 
					this->mv_Ycommit = this->mv_Ytrial;
					this->mv_Vp = this->mv_Vc;
					this->mv_Ap = this->mv_Ac;
				}

			protected:
				/** @brief Previous velocity. */
				std::array<double, nDim> mv_Vp;

				/** @brief Current velocity. */
				std::array<double, nDim> mv_Vc;

				/** @brief Previous acceleration. */
				std::array<double, nDim> mv_Ap;

				/** @brief Current acceleration. */
				std::array<double, nDim> mv_Ac;
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
