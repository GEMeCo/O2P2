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
#include <iostream>				// Required by cout
#include <tuple>				// Required by tuple

namespace O2P2 {
	namespace Prep {
		/** @ingroup PreProcessor_Module
		  * @class Node
		  *
		  * @brief Geometry node, a Domain component.
		  * @details Contains the initial position and an inverse indexing matrix, with the elements connected to the node.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  */
		template<int nDim>
		class Node
		{
		private:
			Node() = delete;

		protected:
			/** Constructor for node objects.
			  * @param index Node index.
			  * @param InitPos Initial position of the geometry node.
			  * @param Tp Initial temperature.
			  */
			Node(const size_t& index, const std::array<double, nDim>& InitPos, const double Tp = 0.) : m_index(index), m_Tp(Tp), v_X(InitPos) {	}

			// Default destructor of private / protected pointers.
			virtual ~Node() = default;

		public:
			/** Overloading operator << to stream the node coordinates. */
			friend std::ostream& operator<<(std::ostream& stream, Node<nDim> const& node) {
				for (auto& x : node.v_X) {
					stream << formatFixed << x;
				}
				return stream;
			}

			/** Whenever an element points to the node, the inverse indexing (v_IncInv) registers the element index container.
			  * @param idx Container index number of the element that is pointing to the node.
			  */
			void addConecToElem(const size_t& idx) {
				this->v_InvInc.push_back(idx);
			}

			/** Retrieve initial nodal coordinates.
			  * @return standard array with nodal initial coordinates.
			  */
			const std::array<double, nDim>& getInitPos() { return this->v_X; }


			/** Retrieve the inverse indexing - elements that are connected to the node.
			  * @return standard vector with elements connected to the node.
			  */
			const std::vector<size_t>& getInvInc() { return this->v_InvInc; }

			/** Retrieve current commited position.
			  * @return standard vector with nodal commited position.
			  */
			virtual const std::array<double, nDim>& getCurrent() = 0;

			/** Retrieve current trial position.
			  * @return standard vector with nodal trial position.
			  */
			virtual const std::array<double, nDim>& getTrial() = 0;

			/** Retrieve current nodal temperature.
			  * @return double with nodal temperature.
			  */
			const double& getCurrentTemp() { return m_Tp; }

			/** Retrieve current trial nodal temperature.
			  * @return double with nodal trial temperature.
			  */
			const double& getTrialTemp() { return m_Tp; }

			/** Update trial position.
			  * @param dPos Increase in the trial position.
			  */
			virtual void updateTrial(const double dPos[]) = 0;

			/** Update trial temperature.
			  * @param dTp New value for trial temperature.
			  */
			void updateTrial(const double dTp) { m_Tp = dTp; }

			/** Commit trial to current. */
			virtual void setCurrent() = 0;

			/** @brief Index of DOF that the node is associated to. */
			std::vector<size_t> v_DofIndex;

			/** @brief Node index. */
			size_t m_index;

		protected:
			/** @brief Room / Commit temperature. */
			double m_Tp;

			/** @brief Initial coordinates. */
			std::array<double, nDim> v_X;

			/** @brief Inverse indexing (elements index). */
			std::vector<size_t> v_InvInc;
		};


		/**
		  * @class Node_Mech_Qse
		  *
		  * @brief Geometry node, a Domain component, for quasi-static mechanical problems.
		  * @details Includes current position (trial and commit). Developed for mechanical quasi-static problems.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  */
		template<int nDim>
		class Node_Mech_Qse : public Node<nDim>
		{
		public:
			Node_Mech_Qse() = delete;

			/** Constructor for node objects for quasi-static mechanical problems.
			  * @param index Node index.
			  * @param InitPos Initial position of the geometry node.
			  * @param Tp Initial temperature.
			  */
			explicit Node_Mech_Qse(const size_t& index, const std::array<double, nDim>& InitPos, const double Tp = 0.)
				: Node<nDim>(index, InitPos, Tp), v_Ytrial(InitPos), v_Ycommit(InitPos) { }

			// Default destructor of private / protected pointers.
			~Node_Mech_Qse() = default;

			// Function to retrieve current nodal coordinates.
			const std::array<double, nDim>& getCurrent() override { return this->v_Ycommit; }

			// Function to retrieve current trial nodal coordinates.
			const std::array<double, nDim>& getTrial() override { return this->v_Ytrial; }

			// Function to increase trial position.
			void updateTrial(const double dPos[]) override {
				for (int i = 0; i < nDim; ++i) {
					this->v_Ytrial[i] += dPos[i];
				}
			};

			// Function to commit trial position to commit position (current).
			void setCurrent() override { this->v_Ycommit = this->v_Ytrial; }

		protected:
			/** @brief Trial current position. */
			std::array<double, nDim> v_Ytrial;

			/** @brief Commit current position. */
			std::array<double, nDim> v_Ycommit;
		};
	} // End of Prep Namespace
} // End of O2P2 Namespace