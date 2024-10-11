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

#include "Node.h"
#include "Material.h"
#include "Section.h"
#include "Element.h"

namespace O2P2 {
	namespace Geom {
		/** @ingroup PreProcessor_Module
		  * @class Domain
		  *
		  * @brief Container of geometry components: nodes, elements, sections and materials.
		  * @details This class holds vectors of geometry components. These components are objects of the following classes:
		  * Node, Point, Element, Inclusion, Material and Section.
		  *
		  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
		  *
		  * @sa Node
		  * @sa Element
		  * @sa Material
		  * @sa CrossSection
		  */
		template<int nDim>
		class Domain
		{
		private:
			// Default constructor
			Domain() = delete;

			// Domain is unique. Copy constructor will be deleted.
			Domain(const Domain& other) = delete;

			// Domain is unique. Move constructor will be deleted.
			Domain(Domain&& other) = delete;

		public:
			/** The only constructor that there is. It resizes mesh containers (Node, Element, Cross Section, Material)
			  * @param nNodes Number of nodes.
			  * @param nPoints Number of particles and fibers nodes.
			  * @param nElem Number of elements.
			  * @param nIncl Number of embedded elements.
			  * @param nMat Number of materials.
			  * @param nSec Number of cross sections.
			  */
			explicit Domain(const size_t& nNodes, const size_t& nPoints, const size_t& nElem, const size_t& nIncl, const size_t& nMat, const size_t& nSec) :
				mv_numNodes(nNodes), mv_numPoints(nPoints), mv_numElem(nElem), mv_numIncl(nIncl), mv_numMat(nMat), mv_numSec(nSec)
			{
				mv_Node.reserve(nNodes);
				mv_Point.reserve(nPoints);
				mv_Mat.reserve(nMat);
				mv_Sect.reserve(nSec);
				mv_Elem.reserve(nElem);
				mv_Incl.reserve(nIncl);

				mv_numBars = 0;
				mv_numFibs = 0;
				mv_numFace = 0;
			}

			// Default destructor
			~Domain() = default;

			/** Add a new node to the geometry node vector.
			  * @param index Node indexing number.
			  * @param x0 Node initial position coordinates.
			  *
			  * @sa Node
			  */
			void addNode(const size_t& index, const std::array<double, nDim>& x0);

			/** Add a new inclusion node to the point vector.
			  * @param index Point indexing number.
			  * @param x0 Point initial position coordinates.
			  *
			  * @sa Node
			  */
			void addPoint(const size_t& index, const std::array<double, nDim>& x0);

			/** Add a new material to the container of materials.
			  * @param index Material number.
			  * @param matType Type of material, see also MaterialType struct.
			  * @param matParam Material properties vector.
			  *
			  * @sa Material
			  */
			void addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);

			/** Add a new section to the container of cross sections.
			  * @param type Type of section (cross section of linear or plane elements).
			  * @param value Area or Thickness value.
			  * @param PSType If element is plane, indicates the plane state.
			  */
			void addSection(const int& type, const double& value, const int& PSType = 1);

			/** Add a new element to the container vector.
			  * @return the number of incidence nodes to be read, after adding a new element to the container.
			  *
			  * @param index Indexing number of current element.
			  * @param type Type of element.
			  * @param order Order of the element.
			  * @param numIP Number of integration points.
			  * @param material Material associated to the element.
			  * @param section Cross section information, if there is, associated to the element.
			  *
			  * @sa Element
			  */
			int addElement(const size_t& index, const int& type, const int& order, const int& numIP, const size_t& material, const size_t& section);

			/** Add the incidence of the element.
			  * @note There is no validation if container index and element index are the same.
			  *
			  * @param index Element container index number - Notice that this is NOT the element index.
			  * @param Type Type of element.
			  * @param conectivity Node incidence.
			  */
			void addElementConectivity(const size_t& index, const int& Type, const std::vector<size_t>& conectivity);

			/** Add an active face element to the container vector.
			  * @note Every active face must retain the pointer to the prism element. This is saved in mv_
			  * 
			  * @param element A pointer to element which face was Actived.
			  * @param face Index of the active face.
			  * @param material Material associated to the element.
			  * @param section Cross section information, if there is, associated to the inclusion element.
			  *
			  * @sa Element
			  */
			void addFaceElem(std::shared_ptr<O2P2::Geom::Elem::Element> element, const int& face, const size_t& material, const size_t& section);

			/** Add a new particle or fiber element element to the container vector.
			  * @return the number of incidence nodes to be read, after adding the new element to the container.
			  *
			  * @param index Indexing number of current element.
			  * @param type Type of element.
			  * @param order Order of the element.
			  * @param numIP Number of integration points.
			  * @param material Material associated to the element.
			  * @param section Cross section information, if there is, associated to the inclusion element.
			  *
			  * @sa Element
			  */
			int addInclusion(const size_t& index, const int& type, const int& order, const int& numIP, const size_t& material, const size_t& section);

			/** Add the incidence of the inclusion element.
			  * @note There is no validation if container index and element index are the same.
			  *
			  * @param index Inclusion container index number - Notice that this is NOT the element index.
			  * @param Type Type of element.
			  * @param Conectivity Point incidence.
			  */
			void addInclusionConectivity(const size_t& index, const int& Type, const std::vector<size_t>& Conectivity);

			/** Get a pointer to a single node from the Node container.
			  * @return a pointer to a single node.
			  * @param index Node container index.
			  */
			const O2P2::Geom::Node* getNode(const size_t& index) {
				return mv_Node.at(index).get();
			}

			/** Get a pointer to a single point from the Point container.
			  * @param index Point container index.
			  * @return a pointer to a inclusion node.
			  */
			const O2P2::Geom::Node* getPoint(const size_t& index) {
				return mv_Point.at(index).get();
			}

			/** Get access to the container of nodes.
			  * @return a pointer to the node's container.
			  * @note Direct access to node's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Node>>& getNode() { return mv_Node; }

			/** Get access to the container of points.
			  * @return a pointer to the inclusion node's container.
			  * @note Direct access to point's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Node>>& getPoint() { return mv_Point; }

			/** Get a pointer to a single element from the Element container index.
			  * @return a pointer to a single element.
			  * @param index Element container index.
			  */
			const O2P2::Geom::Elem::Element* getElem(const size_t& index) {
				return mv_Elem.at(index).get();
			}

			/** Get a pointer to a single element from the Linear Element container index.
			  * @return a pointer to a single element.
			  * @param index Linear Element container index.
			  */
			const O2P2::Geom::Elem::Element* getBars(const size_t& index) {
				return mv_Bars.at(index).get();
			}

			/** Get a pointer to a single element from the Active Face Element container index.
			  * @return a pointer to a single element.
			  * @param index Active Face Element container index.
			  */
			const O2P2::Geom::Elem::Element* getFaces(const size_t& index) {
				return mv_Face.at(index).get();
			}

			/** Get a pointer to a single inclusion element from the Inclusion container index.
			  * @return a pointer to a single inclusion element.
			  * @param index Inclusion container index.
			  */
			const O2P2::Geom::Elem::Element* getIncl(const size_t& index) {
				return mv_Incl.at(index).get();
			}

			/** Get a pointer to a single fiber element from the Fiber Inclusion container index.
			  * @return a pointer to a single fiber element.
			  * @param index Fiber Inclusion container index.
			  */
			const O2P2::Geom::Elem::Element* getFibs(const size_t& index) {
				return mv_Fibs.at(index).get();
			}

			/** Get access to the container of elements.
			  * @return a pointer to the container of elements.
			  * @note Direct access to element's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>>& getElem() { return mv_Elem; }

			/** Get access to the container of linear elements.
			  * @return a pointer to the container of linear elements.
			  * @note Direct access to linear element's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>>& getBars() { return mv_Bars; }

			/** Get access to the container of active face elements.
			  * @return a pointer to the container of active face elements.
			  * @note Direct access to active face element's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>>& getFaces() { return mv_Face; }

			/** Get access to the container of embedded elements.
			  * @return a pointer to the container of inclusion elements.
			  * @note Direct access to inclusion's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>>& getIncl() { return mv_Incl; }

			/** Get access to the container of fiber elements.
			  * @return a pointer to the container of fiber elements.
			  * @note Direct access to fiber element's container should be avoided. This was allowed for fasten up process.
			  */
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>>& getFibs() { return mv_Fibs; }

			/** @return Number of element' nodes.
			  * @param index Element container index.
			  */
			const int getNumElemNodes(const size_t& index) const noexcept { return mv_Elem.at(index)->getNumNodes(); }

			/** @return Number of inclusion' nodes.
			  * @param index Inclusion container index.
			  */
			const int getNumInclNodes(const size_t& index) const noexcept { return mv_Incl.at(index)->getNumNodes(); }

			/** @return Number of domain nodes.
			  */
			inline const size_t getNumNodes() const noexcept { return mv_numNodes; }

			/** @return Number of domain immersed nodes.
			  */
			inline const size_t getNumPoints() const noexcept { return mv_numPoints; }

			/** @return Number of domain elements.
			  */
			inline const size_t getNumElems() const noexcept { return mv_numElem; }

			/** @return Number of domain linear elements (Trusses).
			  */
			inline const size_t getNumBars() const noexcept { return mv_numBars; }

			/** @return Number of domain active face elements (2D in 3D).
			  */
			inline const size_t getNumFace() const noexcept { return mv_numFace; }

			/** @return Number of domain immersed elements.
			  */
			inline const size_t getNumIncl() const noexcept { return mv_numIncl; }

			/** @return Number of domain immersed linear elements (Fibers).
			  */
			inline const size_t getNumFibs() const noexcept { return mv_numFibs; }

			/** @return Number of materials.
			  */
			inline const size_t getNumMat() const noexcept { return mv_numMat; }

			/** @return Number of sections.
			  */
			inline const size_t getNumSec() const noexcept { return mv_numSec; }

		private:
			/** @brief Number of geometry nodes. */
			size_t mv_numNodes;

			/** @brief Number of particles and fibers nodes. */
			size_t mv_numPoints;

			/** @brief Number of geometry elements (total). */
			size_t mv_numElem;

			/** @brief Number of linear elements. */
			size_t mv_numBars;

			/** @brief Number of inclusion elements (total). */
			size_t mv_numIncl;

			/** @brief Number of fiber elements. */
			size_t mv_numFibs;

			/** @brief Number of active face elements. */
			size_t mv_numFace;

			/** @brief Number of materials. */
			size_t mv_numMat;

			/** @brief Number of cross sections. */
			size_t mv_numSec;

			// Container of Node objects.
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mv_Node;

			// Container of inclusion node objects.
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mv_Point;

			// Container of Material objects.
			std::vector<std::shared_ptr<O2P2::Geom::Material>> mv_Mat;

			// Container of Cross section objects.
			std::vector<std::shared_ptr<O2P2::Geom::CrossSection>> mv_Sect;

			// Container of Element objects.
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>> mv_Elem;

			// Container of Bar element objects.
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>> mv_Bars;

			// Container of Active Face element objects.
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>> mv_Face;

			// Container of Inclusion / Embedded element objects.
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>> mv_Incl;

			// Container of Fiber element objects.
			std::vector<std::shared_ptr<O2P2::Geom::Elem::Element>> mv_Fibs;
		};
	} // End of Geom Namespace
} // End of O2P2 Namespace