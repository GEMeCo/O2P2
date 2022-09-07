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
#include <vector>		// required by std::vector
#include <memory>		// required by std::shared_pointer
#include <iostream>		// required by operator << overloading
#include <algorithm>    // std::max / std::max_element
#include <math.h>		// std::sqrt

// Eigen libraries
#include <Eigen/Dense>

// Custom Header Files
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
				explicit BaseElement(const std::shared_ptr<O2P2::Prep::Material>& Material) : m_Mat(Material) { }

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
				virtual Eigen::VectorXd getShapeFcOnPoint(const double* Point) = 0;

				/** Evaluates the derivative of shape function in the point.
				  * @return a matrix with the value of the derivative of shape functions.
				  * @param Point Dimensionless coordinates of the point.
				  */
				virtual Eigen::MatrixXd getShapeDerivOnPoint(const double* Point) = 0;

				/** @return a pointer to the shape functions (with size [m_NumIP][m_NumNodes]). */
				virtual double const* getShapeFc() const = 0;

				/** @return a pointer to the derivative of shape functions (with size [m_NumIP][m_Dim][m_NumNodes]). */
				virtual double const* getShapeDerivative() const = 0;

				/** @return a pointer to the weight of the integation points (with size [m_NumIP]). */
				virtual double const* getWeight() const = 0;

				/** return a vector with values on the integration points currently known in the element' nodes (such as temperature).
				  * @param value element nodal values of the required parameter. There must be the same number of values than there are nodes in the element.
				  */
				virtual Eigen::VectorXd getValueOnIPs(const double* value) = 0;

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
				O2P2::Prep::Material* getMaterial() { return m_Mat.get(); }

			protected:
				/** @brief Pointer to the Material. */
				std::shared_ptr<O2P2::Prep::Material> m_Mat;
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
				explicit Element(const std::shared_ptr<O2P2::Prep::Material>& Material) : BaseElement(Material) { m_Centroid = nullptr; }

				// Default destructor of private / protected pointers.
				~Element() = default;

			public:
				/** Overloading operator << to stream the node indexing. */
				friend std::ostream& operator<<(std::ostream& stream, Element<nDim> const& elem) {
					stream << elem.print();
					return stream;
				}

				/** @return string with incidence for printing. */
				const std::string print() const {
					std::string str;
					for (auto& node : v_Conect) {
						str.append(std::to_string(node->m_index));
						str.append(" ");
					}
					return str;
				}

				/** @return String with geometric properties. */
				const std::string printGeom() const {
					std::string str;
					str.append("Radius: ");
					str.append(std::to_string(m_Radius));
					str.append("; Centroid: ");
					for (int i = 0; i < nDim; i++) {
						str.append(std::to_string(m_Centroid[i]));
						str.append(" ");
					}
					return str;
				}

				/** Populates the indexing matrix and evaluates initial properties.
				  * @param Conect The indexing matrix.
				  */
				void setConectivity(const std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>>& Conect) {
					v_Conect = std::move(Conect);
					this->setGeomProperties();
				}

				/** @return a reference to the element indexing. */
				std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>>& getConectivity() { return v_Conect; }

				/** Get a pointer to a specific node from the conectivity.
				  * @return a pointer to a node.
				  * @param index Element convectivity container index.
				  */
				O2P2::Prep::Node<nDim>* getConectivity(const int& index) { return v_Conect.at(index).get(); }

				/** @return the radius of the element' minimum bounding circle. */
				double getRadius() { return this->m_Radius; }

				/** @return a pointer to the centroid of the element (with size [nDim]). */
				double const* getCentroid() const { return &m_Centroid[0]; }

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates.
				  */
				virtual bool evaluateXsi(const std::array<double, nDim> xsi) { return false; }

			protected:
				/** Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity. */
				virtual void setGeomProperties() {}
				//virtual void setGeomProperties() = 0;

			protected:
				/** @brief Vector with element indexing. */
				std::vector<std::shared_ptr<O2P2::Prep::Node<nDim>>> v_Conect;

				/** @brief Centroid of the polygon (It is not the circumcenter). */
				std::unique_ptr<double[]> m_Centroid;

				/** @brief Radius of the minimum bounding circle. */
				double m_Radius{ 0. };
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
					: Element<nDim>(Material), m_Section(Section) { }

			protected:
				/** @brief Pointer to the Section. */
				std::shared_ptr<O2P2::Prep::Section> m_Section;

				/** @brief Element dimensionality */
				static inline int const m_Dim{ 1 };

			public:
				// Returns the number of DOF per node of current element.
				int getNumNdDOF() override { return nDim; }

				// Returns the dimensionality of current element.
				int getDIM() override { return m_Dim; }

				/** @return a reference to a section object. */
				O2P2::Prep::Section* getSection() { return m_Section.get(); }
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
					: Element<2>(Material), m_Section(Section) { }

			protected:
				/** @brief Pointer to Section. */
				std::shared_ptr<O2P2::Prep::Section> m_Section;

				/** @brief Element dimensionality */
				static inline int const m_Dim{ 2 };

			public:
				// Returns the number of DOF per node of current element.
				int getNumNdDOF() override { return 2; }

				// Returns the dimensionality of current element.
				int getDIM() override { return m_Dim; }

				/** @return a reference to a section object. */
				O2P2::Prep::Section* getSection() { return m_Section.get(); }
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
					: Element<3>(Material) { }

			public:
				// Returns the number of DOF per node of current element.
				int getNumNdDOF() { return 3; }

				// Returns the dimensionality of current element.
				int getDIM() { return m_Dim; }

			protected:
				/** @brief Element dimensionality */
				static inline int const m_Dim{ 3 };
			};
		} // End of Elem Namespace
	} // End of Prep Namespace
} // End of O2P2 Namespace
