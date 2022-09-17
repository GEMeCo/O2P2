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

		public:
			/** Constructor for node objects.
			  * @param index Node index.
			  * @param InitPos Initial position of the geometry node.
			  * @param Tp Initial temperature.
			  */
			explicit Node(const size_t& index, const std::array<double, nDim>& InitPos, const double Tp = 0.)
				: m_index(index), m_Tp(Tp), v_X(InitPos) {	}

			// Default destructor of private / protected pointers.
			~Node() = default;

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
			const std::array<double, nDim>& getInitPos() const { return this->v_X; }


			/** Retrieve the inverse indexing - elements that are connected to the node.
			  * @return standard vector with elements connected to the node.
			  */
			const std::vector<size_t>& getInvInc() { return this->v_InvInc; }

			/** Retrieve initial nodal temperature.
			  * @return double with nodal temperature.
			  */
			const double& getInitialTemp() { return m_Tp; }

			/** @brief Node index. */
			size_t m_index;

		private:
			/** @brief Initial temperature. */
			double m_Tp;

			/** @brief Initial coordinates. */
			std::array<double, nDim> v_X;

			/** @brief Inverse indexing (elements index). */
			std::vector<size_t> v_InvInc;
		};
	} // End of Prep Namespace
} // End of O2P2 Namespace