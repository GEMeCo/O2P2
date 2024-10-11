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
	namespace Geom {
		/** @ingroup PreProcessor_Module
		  * @class Node
		  *
		  * @brief Geometry node, a Domain component.
		  * @details This base class only provides interface for derived classes, which contains initial information on nodes.
		  */
		class Node
		{
		private:
			// Default constructor is deleted. Use explicit constructor only.
			Node() = delete;

		protected:
			/** Constructor for geometry node objects.
			  * @param index Node index.
			  * @param Tp Initial temperature.
			  */
			explicit Node(const size_t& index, const double Tp = 0.)
				: mv_index(index), mv_Tp(Tp) {}

		public:
			// Default destructor of private / protected pointers.
			virtual ~Node() = default;

			/** Overloading operator << to stream the node coordinates. */
			friend std::ostream& operator<<(std::ostream& stream, Node const& node) {
				stream << node.print();
				return stream;
			}

			/** @return string with node coordinates. */
			virtual const std::string print() const = 0;

			/** Whenever an element points to the node, the inverse indexing registers the element index container.
			  * @param idx Container index number of the element that is pointing to the node.
			  */
			void addConectToElem(const size_t& idx) {
				this->mv_invInc.push_back(idx);
			}

			/** Retrieve initial nodal coordinates.
			  * @return standard array with nodal initial coordinates.
			  */
			virtual const double* getInitPos() const = 0;

			/** Retrieve the dimensionality of the vector space.
			  * @return int with the dimensionality.
			  */
			virtual const int getDIM() const = 0;

			/** Retrieve the inverse indexing - elements that are connected to the node.
			  * @return standard vector with elements connected to the node.
			  */
			const std::vector<size_t>& getInvInc() const { return this->mv_invInc; }

			/** Retrieve initial temperature.
			  * @return double with nodal temperature.
			  */
			const double& getInitialTemp() const { return this->mv_Tp; }

		public:
			/** @brief Node index */
			size_t mv_index;

		protected:
			/** @brief Initial temperature. */
			double mv_Tp;

			/** @brief Inverse indexing. */
			std::vector<size_t> mv_invInc;
		};
		  
		/** @ingroup PreProcessor_Module
		  * @class DomainNode
		  * 
		  * @brief Geometry node, a Domain component.
		  * @details Contains the initial position and an inverse indexing matrix, with the elements connected to the node.
		  * 
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  */
		template<int nDim>
		class DomainNode : public Node
		{
		private:
			// Default constructor is deleted. Use explicit constructor only.
			DomainNode() = delete;

		public:
			/** Constructor for node objects.
			* @param index Node index.
			* @param initPos Initial position of the geometry node.
			* @param Tp Initial temperature.
			*/
			explicit DomainNode(const size_t& index, const std::array<double, nDim>& initPos, const double Tp = 0.)
				: Node(index, Tp), mv_x(initPos) {}

			// Default destructor of private / protected pointers.
			~DomainNode() = default;

			// Return a string with initial position for printing.
			const std::string print() const override {
				std::string str;
				for (auto& pos : this->mv_x) {
					str.append(" ");
					str.append(std::to_string(pos));
				}
				return str;
			}

			// Retrieve iniital nodal coordinates
			const double* getInitPos() const override { return this->mv_x.data(); }

			// Retrieve the dimensionality of the vector space.
			const int getDIM() const override { return nDim; };

		protected:
			/** @brief Initial coordinates. */
			std::array<double, nDim> mv_x;
		};
	} // End of Geom Namespace
} // End of O2P2 Namespace