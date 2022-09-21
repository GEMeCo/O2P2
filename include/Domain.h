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

#include "Node.h"
#include "Element.h"
#include "Material.h"
#include "Section.h"

namespace O2P2 {
	namespace Prep {
		/** @ingroup PreProcessor_Module
		  * @class Domain
		  *
		  * @brief Container of geometry components: nodes, elements, sections and materials.
		  * @details This class holds vectors of geometry components. These components are objects of the following classes:
		  * Node, Element, Material and Section.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  * 
		  * @sa Node
		  * @sa Element
		  * @sa Material
		  * @sa Section
		  */
		template<int nDim>
		class Domain
		{
		private:
			// Default constructor
			Domain() = delete;

			// Domain is unique, therefore it should not be copied. Thus, copy constructor will be deleted.
			Domain(const Domain& other) = delete;

			// Move constructor will also be deleted.
			Domain(Domain&& other) = delete;

		public:
			/** The only constructor that there is. It resizes mesh containers (Node, Element, Section, Material)
			  * @param nNodes Number of matrix nodes.
			  * @param nElem Number of matrix elements.
			  * @param nMat Number of materials.
			  * @param nSec Number of cross sections.
			  */
			explicit Domain(const size_t& nNodes, const size_t& nElem, const size_t& nMat, const size_t& nSec) :
				m_nNodes(nNodes), m_nElem(nElem), m_nMat(nMat), m_nSec(nSec) {

				// Resize containers
				v_Node.reserve(nNodes);
				v_Mat.reserve(nMat);
				v_Sect.reserve(nSec);
				v_Elem.reserve(nElem);
			}

			// Default destructor.
			~Domain() = default;

			/** Add a new node to the matrix node vector.
			  * @param index Node indexing number.
			  * @param x0 Node initial position coordinates.
			  *
			  * @sa Node
			  */
			void addGeomNode(const size_t& index, const std::array<double, nDim>& x0);

			/** Add a new material to the container of materials.
			  * @param index Material number.
			  * @param matType Type of material, see also MaterialType struct.
			  * @param matParam Material properties vector.
			  *
			  * @sa Material
			  */
			void addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);

			/** Sets the material properties to the material with the input index. It throws an out_of_range exception if index is out of bounds.
			  * @note There is no validation if container index and material index are the same.
			  *
			  * @param index Material container index number - Notice that this is NOT the material index.
			  * @param Parameters Material parameters.
			  *
			  * @sa Material
			  */
			void addMaterialParam(const size_t& index, const std::vector<double>& Parameters) {
				v_Mat.at(index)->setParameters(Parameters);
			}

			/** Add a new section to the container of cross sections.
			  * @param data Section parameters.
			  */
			void addSection(std::istringstream& data);

			/** Add a new element to the container vector.
			  * @return the number of incidence nodes to be read, after adding a new element to the container.
			  * @note Only bidimensional triangular and rectangular elements are accepted so far.
			  *
			  * @param index Indexing number of current element.
			  * @param Type Type of element.
			  * @param Order Order of the element.
			  * @param numIP Number of integration points.
			  * @param Material Material associated to the element.
			  * @param Section Section, if there is, associated to the element.
			  */
			int addElement(const size_t& index, const int& Type, const int& Order, const int& numIP, const size_t& Material, const size_t& Section);

			/** Add the incidence of the element.
			  * @note There is no validation if container index and element index are the same.
			  *
			  * @param index Element container index number - Notice that this is NOT the element index.
			  * @param Conectivity Node incidence.
			  */
			void addElementConect(const size_t& index, const std::vector<size_t>& Conectivity);

			/** Get a single node from the Node container index.
			  * @return a pointer to a single node.
			  * @param index Node container index.
			  */
			const O2P2::Prep::Node<nDim>* getNode(const size_t& index) {
				return v_Node.at(index).get();
			}

			/** Get access to the container of nodes.
			  * @return a pointer to the node's container.
			  * @note Direct access to node's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>>& getNode() { return v_Node; }

			/** Get a single element from the Element container index.
			  * @return a pointer to a single element.
			  * @param index Element container index.
			  */
			O2P2::Prep::Elem::Element<nDim>* getElem(const size_t& index) { return v_Elem.at(index).get(); }

			/** Get access to the container of elements.
			  * @return a pointer to the container of elements.
			  * @note Direct access to element's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Prep::Elem::Element<nDim>>>& getElem() { return v_Elem; }

			/** @return Number of element' nodes.
			  * @param index Element container index.
			  */
			int getNumNodes(const size_t& index) { return v_Elem.at(index)->getNumNodes(); }

		public:
			/** @brief Number of geometry nodes. */
			size_t m_nNodes;

			/** @brief Number of geometry elements. */
			size_t m_nElem;

			/** @brief Number of materials. */
			size_t m_nMat;

			/** @brief Number of cross sections. */
			size_t m_nSec;

		private:
			// Container of Node objects.
			std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>> v_Node;

			// Container of Material objects.
			std::vector<std::shared_ptr<O2P2::Prep::Material>> v_Mat;

			// Container of Section objects.
			std::vector<std::shared_ptr<O2P2::Prep::Section>> v_Sect;

			// Container of Element objects.
			std::vector<std::shared_ptr<O2P2::Prep::Elem::Element<nDim>>> v_Elem;
		};
	} // End of Prep Namespace
} // End of O2P2 Namespace