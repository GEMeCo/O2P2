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

// Custom header files
#include "Common.h"

#include "Material.h"
#include "Section.h"
#include "Node.h"

namespace O2P2 {
	namespace Prep {
		namespace Elem {
			/** @ingroup PreProcessor_Module
			  * @class BaseElement
			  *
			  * @brief Base geometry element, a Domain component.
			  * @details Contains the basic definitions for geometry elements (Domain).
			  * It is also a container for the integration points, weights and functions.
			  */
			class BaseElement
			{
			private:
				BaseElement() = delete;

			protected:
				/** Constructor for element objects.
				  * @param Material Pointer to Material class.
				  */
				explicit BaseElement(const std::shared_ptr<O2P2::Prep::Material>& Material) : mv_Mat(Material) {}

				// Default destructor of private / protected pointers.
				virtual ~BaseElement() = default;

			public:
				/** Output function for AcadView, based on element index.
				  * @return a string with incidence for printing in AcadView.
				  * @param add Adder to number indexing.
				  */
				virtual const std::string printByIndex_AV(const size_t add) const = 0;

				/** Output function for AcadView, based on element node number.
				  * @return a string with incidence for printing in AcadView.
				  * @param add Adder to number incidence.
				  */
				virtual const std::string printByAdder_AV(const size_t add) const = 0;

				/** Evaluates shape function in the point.
				  * @return a vector with the shape function.
				  * @param Point Dimensionless coordinates of the point.
				  */
				virtual std::vector<double> getShapeFcOnPoint(const double* Point) = 0;

				/** Evaluates the derivative of shape function in the point.
				  * @return a matrix with the value of the derivative of shape functions.
				  * @param Point Dimensionless coordinates of the point.
				  */
				virtual std::vector<double> getShapeDerivOnPoint(const double* Point) = 0;

				/** @return a pointer to the shape functions (with size [mv_numIP][mv_numNodes]). */
				virtual double const* getShapeFc() const = 0;

				/** @return a pointer to the derivative of shape functions (with size [mv_numIP][mv_Dim][mv_numNodes]). */
				virtual double const* getShapeDerivative() const = 0;

				/** @return a pointer to the weight of the integation points (with size [mv_numIP]). */
				virtual double const* getWeight() const = 0;

				/** @return the number of nodes of current element. */
				virtual int getNumNodes() = 0;

				/** @return the number of faces of current element. */
				virtual int getNumFaces() = 0;

				/** @return the number of degrees of freedom per node for current element. */
				virtual int getNumNdDOF() = 0;

				/** @return the number of integration points for current element. */
				virtual int getNumIP() = 0;

				/** @return the dimensionality of current element. */
				virtual int getDIM() = 0;

				/** @return a reference to the elements Material object. */
				O2P2::Prep::Material* getMaterial() { return mv_Mat.get(); }

				virtual std::vector<double> getValueOnIPs(const double* value) = 0;

			protected:
				/** @brief Pointer to the Material. */
				std::shared_ptr<O2P2::Prep::Material> mv_Mat;
			};


			/**
			  * @class Element
			  *
			  * @brief Geometry element for derived classes.
			  * @details Includes element indexing and element centroid.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class Element : public BaseElement
			{
			private:
				Element() = delete;

			protected:
				/** Constructor for element objects.
				  * @param Material Pointer to Material class.
				  */
				explicit Element(const std::shared_ptr<O2P2::Prep::Material>& Material)
					: BaseElement(Material) {
					mv_Centroid = nullptr;
				};

				// Default destructor of private / protected pointers.
				~Element() = default;

				/** Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity. */
				virtual void setGeomProperties() {}

			public:
				/** Overloading operator << to stream the node indexing. */
				friend std::ostream& operator<<(std::ostream& stream, Element<nDim> const& elem) {
					stream << elem.print();
					return stream;
				}

				/** @return string with incidence for printing. */
				const std::string print() const {
					std::string str;
					for (auto& node : mv_Conect) {
						str.append(std::to_string(node->mv_index + 1));
						str.append(" ");
					}
					return str;
				}

				/** @return String with geometric properties. */
				const std::string printGeom() const {
					std::string str;
					str.append("Radius: ");
					str.append(std::to_string(mv_Radius));
					str.append("; Centroid: ");
					for (int i = 0; i < nDim; i++) {
						str.append(std::to_string(mv_Centroid[i]));
						str.append(" ");
					}
					return str;
				}

				/** Populates the indexing matrix and evaluates initial properties.
				  * @param Conect The indexing matrix.
				  */
				void setConectivity(const std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>>& Conect) {
					mv_Conect = std::move(Conect);
					this->setGeomProperties();
				}

				/** @return a reference to the element indexing. */
				std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>>& getConectivity() { return mv_Conect; }

				/** Get a pointer to a specific node from the conectivity.
				  * @return a pointer to a node.
				  * @param index Element convectivity container index.
				  */
				O2P2::Prep::Node<nDim>* getConectivity(const int& index) { return mv_Conect.at(index).get(); }

				/** @return the radius of the element' minimum bounding circle. */
				double getRadius() { return this->mv_Radius; }

				/** @return a pointer to the centroid of the element (with size [nDim]). */
				double const* getCentroid() const { return &mv_Centroid[0]; }

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates.
				  */
				virtual inline bool evaluateXsi(const std::array<double, nDim> xsi) { return false; }

			protected:
				/** @brief Vector with element indexing. */
				std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>> mv_Conect;

				/** @brief Centroid of the polygon (It is not the circumcenter). */
				std::unique_ptr<double[]> mv_Centroid;

				/** @brief Radius of the minimum bounding circle. */
				double mv_Radius{ 0. };
			};


			/** @class ElementLinear
			  *
			  * @brief Base geometry class for linear elements.
			  * @details Holds basic definitios for geometry linear elements.
			  *
			  */
			template<int nDim>
			class ElementLinear : public Element<nDim>
			{
			private:
				ElementLinear() = delete;

			protected:
				/** Constructor for linear element objects.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				ElementLinear(std::shared_ptr<O2P2::Prep::Material>& Material, std::shared_ptr<O2P2::Prep::Section>& Section)
					: Element<nDim>(Material), mv_Section(Section) {}

			protected:
				/** @brief Pointer to the Section. */
				std::shared_ptr<O2P2::Prep::Section> mv_Section;

				/** @brief Element dimensionality */
				static inline int const mv_Dim{ 1 };

			public:
				// Returns the number of DOF per node of current element.
				int getNumNdDOF() override { return nDim; }

				// Returns the dimensionality of current element.
				int getDIM() override { return mv_Dim; }

				/** @return a reference to a section object. */
				O2P2::Prep::Section* getSection() { return mv_Section.get(); }
			};

			/**
			  * @class ElementPlane
			  *
			  * @brief Base geometry class for plane elements.
			  * @details Includes basic definitions for plane elements.
			  *
			  */
			class ElementPlane : public Element<2>
			{
			private:
				ElementPlane() = delete;

			protected:
				/** Constructor for plane element objects.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				ElementPlane(std::shared_ptr<O2P2::Prep::Material>& Material, std::shared_ptr<O2P2::Prep::Section>& Section)
					: Element<2>(Material), mv_Section(Section) {}

			protected:
				/** @brief Pointer to Section. */
				std::shared_ptr<O2P2::Prep::Section> mv_Section;

				/** @brief Element dimensionality */
				static inline int const mv_Dim{ 2 };

			public:
				// Returns the number of DOF per node of current element.
				int getNumNdDOF() override { return 2; }

				// Returns the dimensionality of current element.
				int getDIM() override { return mv_Dim; }

				/** @return a reference to a section object. */
				O2P2::Prep::Section* getSection() { return mv_Section.get(); }
			};


			/**
			  * @class ElementSolid
			  *
			  * @brief Base geometry class for solid elements.
			  * @details Includes basic definitions for solid elements.
			  *
			  */
			class ElementSolid : public Element<3>
			{
			private:
				ElementSolid() = delete;

			protected:
				/** Constructor for plane solid objects.
				  * @param Material Pointer to Material class.
				  */
				ElementSolid(std::shared_ptr<O2P2::Prep::Material>& Material)
					: Element<3>(Material) {}

			protected:
				/** @brief Element dimensionality */
				static inline int const mv_Dim{ 3 };

			public:
				// Returns the number of DOF per node of current element.
				int getNumNdDOF() override { return 3; }

				// Returns the dimensionality of current element.
				int getDIM() override { return mv_Dim; }
			};
		} // End of Elem Namespace
	} // End of Prep Namespace
} // End of O2P2 Namespace
