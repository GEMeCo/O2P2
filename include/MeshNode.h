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
				virtual const double* getCurrent() = 0;

				/** Retrieve current trial position.
				  * @return pointer to nodal trial position (with size nDim).
				  */
				virtual const double* getTrial() = 0;

				/** Retrieve current nodal temperature.
				  * @return double with nodal temperature.
				  */
				const double& getCurrentTemp() { return m_Tp; }

				/** Retrieve current trial nodal temperature.
				  * @return double with nodal trial temperature.
				  */
				const double& getTrialTemp() { return m_Tp; }

				/** Update trial position.
				  * @param dPos Increase in the trial position (with size nDim).
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
			  * @class MeshNode_MQ
			  *
			  * @brief Mesh node, a Solution component, for mechanical quasi-static problems.
			  * @details The solution node component for mechanical quasi-static problems includes current position (trial and commit).
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class MeshNode_MQ : public MeshNode
			{
			private:
				MeshNode_MQ() = delete;

			public:
				/** Constructor for node objects for quasi-static mechanical problems.
				  * @param index Node DOF index (the first one).
				  * @param InitPos Initial position of the corresponding geometry node. Use to initialize trial position.
				  * @param Tp Initial temperature.
				  */
				explicit MeshNode_MQ(const size_t& index, const std::array<double, nDim>& InitPos, const double Tp = 0.)
					: MeshNode(index, Tp), v_Ytrial(InitPos), v_Ycommit(InitPos) { }

				// Default destructor of private / protected pointers.
				~MeshNode_MQ() = default;

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
				const double* getCurrent() override { return this->v_Ycommit.data(); }

				// Function to retrieve current trial nodal coordinates.
				const double* getTrial() override { return this->v_Ytrial.data(); }

				// Function to increase trial position.
				void updateTrial(const double dPos[]) override {
					for (int i = 0; i < nDim; ++i) {
						this->v_Ytrial[i] += dPos[i];
					}
				}

				// Function to commit trial position to commit position (current).
				void setCurrent() override { this->v_Ycommit = this->v_Ytrial; }

				/** @return the number of DOF per node.
				  */
				const int getNumDOF() override { return nDim; };

			protected:
				/** @brief Trial current position. */
				std::array<double, nDim> v_Ytrial;

				/** @brief Commit current position. */
				std::array<double, nDim> v_Ycommit;
			};
		} // End of Comp Namespace
	} // End of Proc Namespace
} // End of O2P2 Namespace
