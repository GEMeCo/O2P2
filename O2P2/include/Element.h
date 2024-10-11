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
	namespace Geom {
		namespace Elem {
			/** @ingroup PreProcessor_Module
			  * @class Element
			  *
			  * @brief Geometry element for derived classes.
			  * @details Contains the basic definitions for geometry elements (Domain), including element indexing and element centroid.
			  * It is also a container for the integration points, weights and functions.
			  */
			class Element
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				Element() = delete;

			protected:
				/** Constructor for element objects.
				  * @param Material Pointer to Material class.
				  */
				explicit Element(const std::shared_ptr<O2P2::Geom::Material>& Material)
					: mv_Mat(Material) {
					mv_Centroid = nullptr;
				}

				// Default destructor of private / protected pointers.
				virtual ~Element() = default;

				/** Evaluate centroid and circumsphere Radius. Must be called after setting the conectivity. */
				virtual void setGeomProperties() {};

			public:
				/** Overloading operator << to stream the node indexing. */
				friend std::ostream& operator<<(std::ostream& stream, Element const& elem) {
					stream << elem.print();
					return stream;
				}

				/** @return string with incidence for printing. */
				const std::string print() const {
					std::string str;
					for (auto& node : mv_Conect) {
						str.append(" ");
						str.append(std::to_string(node->mv_index));
					}
					return str;
				}

				/** @return String with geometric properties. */
				virtual const std::string printGeom() const = 0;

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

				/** @return a pointer to the derivative of shape functions (with size [mv_numIP][mv_ElDim][mv_numNodes]). */
				virtual double const* getShapeDerivative() const = 0;

				/** @return a pointer to the weight of the integation points (with size [mv_numIP]). */
				virtual double const* getWeight() const = 0;

				/** @return the number of nodes of current element. */
				virtual const int getNumNodes() const = 0;

				/** @return the number of faces of current element. */
				virtual const int getNumFaces() const = 0;

				/** @return the number of degrees of freedom per node for current element. */
				virtual const int getNumNodalDof() const = 0;

				/** @return the number of integration points for current element. */
				virtual const int getNumIP() const = 0;

				/** @return the dimensionality of current element. */
				virtual const int getDIM() const = 0;

				/** Evaluates values on integration point based on input.
				  * @return a vector with the approximated value on integration points.
				  * @param value Vector of values currently known on nodes.
				  */
				virtual std::vector<double> getValueOnIPs(const double* value) = 0;

				/** @return a reference to the elements Material object. */
				O2P2::Geom::Material* getMaterial() { return mv_Mat.get(); }

				/** Get a pointer to a specific node from the conectivity.
				  * @return a pointer to a node.
				  * @param index Element convectivity container index.
				  */
				O2P2::Geom::Node* getConectivity(const int& index) { return mv_Conect.at(index).get(); }

				/** @return a reference to the element indexing. */
				std::vector<std::shared_ptr<O2P2::Geom::Node>>& getConectivity() { return mv_Conect; }

				/** @return the radius of the element' minimum bounding circle. */
				double getRadius() { return this->mv_Radius; }

				/** @return a pointer to the centroid of the element (with size [nDim]). */
				double const* getCentroid() const { return &mv_Centroid[0]; }

				/** Populates the indexing matrix and evaluates initial properties.
				  * @param Conect The indexing matrix.
				  */
				void setConectivity(std::vector<std::shared_ptr<O2P2::Geom::Node>>& Conect) {
					mv_Conect = std::move(Conect);
					this->setGeomProperties();
				}

				/** Verifies dimensionless coordinates from input - if it is immersed on the element.
				  * @return True if input falls within the element.
				  * @param xsi Trial dimensionless coordinates (with size [nDim]).
				  */
				virtual inline bool evaluateXsi(const double* xsi) = 0;

			protected:
				/** @brief Pointer to the Material. */
				std::shared_ptr<O2P2::Geom::Material> mv_Mat;

				/** @brief Vector with element indexing. */
				std::vector<std::shared_ptr<O2P2::Geom::Node>> mv_Conect;

				/** @brief Centroid of the polygon (It is not the circumcenter). */
				std::unique_ptr<double[]> mv_Centroid;

				/** @brief Radius of the minimum bounding circle. */
				double mv_Radius{ 0. };
			};


			/** @class ElemLinear
			  *
			  * @brief Base geometry class for linear elements.
			  * @details Holds basic definitios for geometry linear elements.
			  *
			  * @tparam nDim The dimensionality of the problem. It is either 2 or 3 (bidimensional or tridimensional).
			  */
			template<int nDim>
			class ElemLinear : public Element
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				ElemLinear() = delete;

			protected:
				/** Constructor for linear element objects.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to Section class.
				  */
				explicit ElemLinear(std::shared_ptr<O2P2::Geom::Material>& Material, std::shared_ptr<O2P2::Geom::CrossSection>& Section)
					: Element(Material), mv_Section(Section) {}

			protected:
				/** @brief Pointer to the Section. */
				std::shared_ptr<O2P2::Geom::CrossSection> mv_Section;

				/** @brief Element dimensionality */
				static inline int const mv_ElDim{ 1 };

			private:
				// There is no need to evaluateXsi in Linear Elements. It is an empty function
				inline bool evaluateXsi(const double* xsi) override { return false; }

			public:
				// Returns string with geometric properties.
				const std::string printGeom() const override {
					std::string str;
					str.append("Approx. Length: ");
					str.append(std::to_string(this->mv_Radius));
					str.append("; Centroid:");
					for (int i = 0; i < nDim; i++) {
						str.append(" ");
						str.append(std::to_string(this->mv_Centroid[i]));
					}
					return str;
				}

				// Returns the number of DOF per node of current element.
				const int getNumNodalDof() const override { return nDim; }

				// Returns the dimensionality of current element.
				const int getDIM() const override { return mv_ElDim; }

				/** @return a reference to a section object. */
				O2P2::Geom::CrossSection* getSection() { return mv_Section.get(); }
			};

			/** @class ElemPlane
			  *
			  * @brief Base geometry class for plane elements.
			  * @details Includes basic definitions for plane elements.
			  */
			class ElemPlane : public Element
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				ElemPlane() = delete;

			protected:
				/** Constructor for plane element objects.
				  * @param Material Pointer to Material class.
				  * @param Section Pointer to PlaneSection class.
				  */
				explicit ElemPlane(std::shared_ptr<O2P2::Geom::Material>& Material, std::shared_ptr<O2P2::Geom::PlaneSection>& Section)
					: Element(Material), mv_Section(Section) {}

			protected:
				/** @brief Pointer to PlaneSection. */
				std::shared_ptr<O2P2::Geom::PlaneSection> mv_Section;

				/** @brief Element dimensionality */
				static inline int const mv_ElDim{ 2 };

			public:
				// Returns string with geometric properties.
				const std::string printGeom() const override {
					std::string str;
					const int mi_Dim = mv_Conect.at(0)->getDIM();	// Dimensionality of vector space (2D or 3D)

					str.append("Radius: ");
					str.append(std::to_string(this->mv_Radius));
					str.append("; Centroid:");
					for (int i = 0; i < mi_Dim; i++) {
						str.append(" ");
						str.append(std::to_string(this->mv_Centroid[i]));
					}
					return str;
				}

				// Returns the number of DOF per node of current element.
				const int getNumNodalDof() const override { return 2; }

				// Returns the dimensionality of current element.
				const int getDIM() const override { return mv_ElDim; }

				/** @return a reference to a section object. */
				O2P2::Geom::PlaneSection* getSection() { return mv_Section.get(); }
			};


			/** @class ElemSolid
			  *
			  * @brief Base geometry class for solid elements.
			  * @details Includes basic definitions for solid elements.
			  */
			class ElemSolid : public Element
			{
			private:
				// Default constructor is deleted. Use explicit constructor only.
				ElemSolid() = delete;

			protected:
				/** Constructor for plane solid objects.
				  * @param Material Pointer to Material class.
				  */
				explicit ElemSolid(std::shared_ptr<O2P2::Geom::Material>& Material)
					: Element(Material) {}

			protected:
				/** @brief Element dimensionality */
				static inline int const mv_ElDim{ 3 };

			public:
				// Returns string with geometric properties.
				const std::string printGeom() const override {
					std::string str;
					str.append("Radius: ");
					str.append(std::to_string(this->mv_Radius));
					str.append("; Centroid:");
					for (int i = 0; i < mv_ElDim; i++) {
						str.append(" ");
						str.append(std::to_string(this->mv_Centroid[i]));
					}
					return str;
				}

				// Returns the number of DOF per node of current element.
				const int getNumNodalDof() const override { return 3; }

				// Returns the dimensionality of current element.
				const int getDIM() const override { return mv_ElDim; }
			};
		} // End of Elem Namespace
	} // End of Geom Namespace
} // End of O2P2 Namespace
