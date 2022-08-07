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
//
// Domain
// 
// Container of geometry components: nodes, points, elements, particles, fibers and materials
// 
// ================================================================================================
#include "Domain.h"
#include "Elements/Elem_Tri3.h"
#include "Elements/Elem_Tri6.h"
#include "Elements/Elem_Tri10.h"

#include "Elements/Elem_Rect4.h"
#include "Elements/Elem_Rect8.h"
#include "Elements/Elem_Rect9.h"
#include "Elements/Elem_Rect12.h"
#include "Elements/Elem_Rect16.h"
#include "Elements/Elem_Rect20.h"


// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template void Domain<2>::addMatrixNode(const size_t& index, const AnalysisType& AnType, const ProblemType& PrType, const std::array<double, 2>& x0);

template void Domain<2>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);

template void Domain<2>::addSection(std::istringstream& data);

template void Domain<2>::addElementConect(const size_t& index, const std::vector<size_t>& Conectivity);

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addMatrixNode
//
// ================================================================================================
template<int nDim>
void Domain<nDim>::addMatrixNode(const size_t& index, const AnalysisType& AnType, const ProblemType& PrType, const std::array<double, nDim>& x0)
{
	// Node is related to the type of problem and analysis
	if (PrType == ProblemType::MECHANICAL) {
		switch (AnType)
		{
		case AnalysisType::STATIC:
			v_Node.emplace_back(std::make_shared<Node_Mech_Qse<nDim>>(index, x0));
			break;
		default:
			LOG("\n\nDomain.addMatrixNode: Node creation error - Could not find the correct AnalysisType for Mechanical Problem\n\n\n");
			throw std::invalid_argument("\n\n\nNode creation error - Could not find the correct AnalysisType for Mechanical Problem\n\n\n");
			break;
		}
	}
	else {
		LOG("\n\nDomain.addMatrixNode: Node creation error - Coupled analysis not yet implemented\n\n\n");
		throw std::invalid_argument("\n\n\nNode creation error - Coupled analysis not yet implemented\n\n\n");
	};

	LOG("Domain.addMatrixNode: Adding Node: " << std::to_string(index) << *v_Node.back());
};


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addMaterial
//
// ================================================================================================
template<int nDim>
void Domain<nDim>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam)
{
	LOG("Domain.addMaterial: Adding Material " << MaterialTypeNames[matType] << " with index " << std::to_string(index));

	// The type of material defines which class Material
	switch (matType)
	{
	case MaterialType::SVK_ISO:
	{
		//1 - Elastic, Isotropic, Saint - Venant - Kirchhoff
		v_Mat.emplace_back(std::make_shared<Mat_SVK_ISO>(index, matParam));
		LOG("Domain.addMaterial: added Material: " << std::to_string(index) << *(std::dynamic_pointer_cast<Mat_SVK_ISO>(v_Mat.back())));
		break;
	}
	default:
		LOG("\n\nDomain.addMaterial: Undefined material. Check input file.\n\n\n");
		throw std::invalid_argument("\n\n\n***** Error in material definition, check input *****\n\n\n");
		break;
	};
};


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addSection
//
// ================================================================================================
template<int nDim>
void Domain<nDim>::addSection(std::istringstream& data)
{
	size_t iAux;
	double dAux;

	// For now, only two options are provided
	// Cross section (for truss elements) 
	// Thickness and Plane State (for plane elements)
	data >> dAux;
	if (!data.eof()) {
		data >> iAux;
		iAux--;
		LOG("Domain.addSection: Plane Thickness: " << std::to_string(dAux) << "; " << PlaneStateTypeNames[PlaneStateType(iAux)]);
		v_Sect.emplace_back(std::make_shared<Section>(PlaneStateType(iAux), dAux));
	}
	else {
		LOG("Domain.addSection: Bar Cross Section: " << std::to_string(dAux));
		v_Sect.emplace_back(std::make_shared<Section>(dAux));
	}
}


// ================================================================================================
//
// Specialized implementation of Template Member Function (only 2D): addElement
//
// ================================================================================================
template<> int Domain<2>::addElement(const size_t& index, const int& Type, const int& Order, const int& numIP, const size_t& Material, const size_t& Section)
{
	// Return value - number of incidence elements
	int Param = 0;

	// The type of element defines which class Element
	switch (Type)
	{
	case 1:
	{
		//1 - Bar / Truss element
		LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. There is no implementation for Trusses. Element used only as immersed fibers.\n\n\n");
		throw std::invalid_argument("\n\n\nThere is no implementation for Trusses so far. Linear element used only as immersed fibers.\n\n\n");
		break;
	}
	case 2:
	{
		//2 - Triangular element

		// Number of nodes for incidence
		Param = (Order + 1) * (Order + 2) / 2;

		switch (Order)
		{
		case 1:
		{
			//1	- Linear interpolation (3 nodes)
			switch (numIP) {
			case 1:  { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<1>>(v_Mat[Material], v_Sect[Section])); break; }
			case 3:  { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<3>>(v_Mat[Material], v_Sect[Section])); break; }
			case 4:  { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<4>>(v_Mat[Material], v_Sect[Section])); break; }
			case 6:  { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<6>>(v_Mat[Material], v_Sect[Section])); break; }
			case 7:  { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<7>>(v_Mat[Material], v_Sect[Section])); break; }
			case 12: { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<12>>(v_Mat[Material], v_Sect[Section])); break; }
			case 13: { v_Elem.emplace_back(std::make_shared<Elem_Tri3_IP<13>>(v_Mat[Material], v_Sect[Section])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (6 nodes)
			switch (numIP) {
			case 1:  { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<1>>(v_Mat[Material], v_Sect[Section])); break; }
			case 3:  { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<3>>(v_Mat[Material], v_Sect[Section])); break; }
			case 4:  { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<4>>(v_Mat[Material], v_Sect[Section])); break; }
			case 6:  { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<6>>(v_Mat[Material], v_Sect[Section])); break; }
			case 7:  { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<7>>(v_Mat[Material], v_Sect[Section])); break; }
			case 12: { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<12>>(v_Mat[Material], v_Sect[Section])); break; }
			case 13: { v_Elem.emplace_back(std::make_shared<Elem_Tri6_IP<13>>(v_Mat[Material], v_Sect[Section])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (10 nodes)
			switch (numIP) {
			case 1:  { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<1>>(v_Mat[Material], v_Sect[Section])); break; }
			case 3:  { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<3>>(v_Mat[Material], v_Sect[Section])); break; }
			case 4:  { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<4>>(v_Mat[Material], v_Sect[Section])); break; }
			case 6:  { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<6>>(v_Mat[Material], v_Sect[Section])); break; }
			case 7:  { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<7>>(v_Mat[Material], v_Sect[Section])); break; }
			case 12: { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<12>>(v_Mat[Material], v_Sect[Section])); break; }
			case 13: { v_Elem.emplace_back(std::make_shared<Elem_Tri10_IP<13>>(v_Mat[Material], v_Sect[Section])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		default:
			// Higher order
			LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Only up to cubic order is provided.\n\n\n");
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	case 3:
	{
		//3 - Rectangular element

		// Number of nodes for incidence
		if (Order < 4) {
			Param = (Order + 1) * (Order + 1);
		}
		else if (Order < 6) {
			Param = (4) * (Order - 2);
		}
		else Param = 20;

		switch (Order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
			case 4: { v_Elem.emplace_back(std::make_shared<Elem_Rect4_IP<4>>(v_Mat[Material], v_Sect[Section])); break; }
			case 9: { v_Elem.emplace_back(std::make_shared<Elem_Rect4_IP<9>>(v_Mat[Material], v_Sect[Section])); break; }
			case 16: { v_Elem.emplace_back(std::make_shared<Elem_Rect4_IP<16>>(v_Mat[Material], v_Sect[Section])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (9 nodes)
			switch (numIP) {
			case 4: { v_Elem.emplace_back(std::make_shared<Elem_Rect9_IP<4>>(v_Mat[Material], v_Sect[Section])); break; }
			case 9: { v_Elem.emplace_back(std::make_shared<Elem_Rect9_IP<9>>(v_Mat[Material], v_Sect[Section])); break; }
			case 16: { v_Elem.emplace_back(std::make_shared<Elem_Rect9_IP<16>>(v_Mat[Material], v_Sect[Section])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (16 nodes)
			switch (numIP) {
			case 4: { v_Elem.emplace_back(std::make_shared<Elem_Rect16_IP<4>>(v_Mat[Material], v_Sect[Section])); break; }
			case 9: { v_Elem.emplace_back(std::make_shared<Elem_Rect16_IP<9>>(v_Mat[Material], v_Sect[Section])); break; }
			case 16: { v_Elem.emplace_back(std::make_shared<Elem_Rect16_IP<16>>(v_Mat[Material], v_Sect[Section])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 4:
		{
			//4	- Cubic / linear interpolation (8 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Rect8>(v_Mat[Material], v_Sect[Section]));
			break;
		}
		case 5:
		{
			//5 - Cubic / quadratic interpolation (12 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Rect12>(v_Mat[Material], v_Sect[Section]));
			break;
		}
		case 6:
		{
			//6	- Cubic / quartic interpolation (20 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Rect20>(v_Mat[Material], v_Sect[Section]));
			break;
		}
		default:
			// Higher order
			LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Only up to cubic order is provided.\n\n\n");
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	default:
		LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Invalid element type.\n\n\n");
		throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Invalid Element Type\n\n\n");
		break;
	}

	return Param;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addElementConect
//
// ================================================================================================
template<int nDim>
void Domain<nDim>::addElementConect(const size_t& index, const std::vector<size_t>& Conectivity)
{
	std::vector<std::shared_ptr<Node<nDim>>> Conect(Conectivity.size());

	for (size_t i = 0; i < Conectivity.size(); i++) {
		Conect[i] = v_Node[Conectivity[i]];

		// Whenever a element points to a node, it is added to the inverse incidence
		Conect[i]->addConecToElem(index);
	}
	v_Elem.at(index)->setConectivity(Conect);

	LOG("Domain.addElement: Element " << std::to_string(index+1) << " indexing: " << *v_Elem.back());
	LOG("Domain.addElement: " << v_Elem.back()->printGeom());

	return;
}
