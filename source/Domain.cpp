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

#include "Elements/Elem_Tet4.h"
#include "Elements/Elem_Tet10.h"
#include "Elements/Elem_Tet20.h"

#include "Elements/Elem_Hex8.h"
#include "Elements/Elem_Hex27.h"
#include "Elements/Elem_Hex64.h"

#include "Elements/Elem_Pri6.h"
#include "Elements/Elem_Pri18.h"
#include "Elements/Elem_Pri20.h"
#include "Elements/Elem_Pri30.h"
#include "Elements/Elem_Pri40.h"


// ================================================================================================
//
// Explicit template member functions instantiation
//
// ================================================================================================
template void Domain<2>::addMatrixNode(const size_t& index, const AnalysisType& AnType, const ProblemType& PrType, const std::array<double, 2>& x0);
template void Domain<3>::addMatrixNode(const size_t& index, const AnalysisType& AnType, const ProblemType& PrType, const std::array<double, 3>& x0);

template void Domain<2>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);
template void Domain<3>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);

template void Domain<2>::addSection(std::istringstream& data);
template void Domain<3>::addSection(std::istringstream& data);

template void Domain<2>::addElementConect(const size_t& index, const std::vector<size_t>& Conectivity);
template void Domain<3>::addElementConect(const size_t& index, const std::vector<size_t>& Conectivity);

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
// Specialized implementation of Template Member Function (only 3D): addElement
//
// ================================================================================================
template<> int Domain<3>::addElement(const size_t& index, const int& Type, const int& Order, const int& numIP, const size_t& Material, const size_t& Section)
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
		//2 - Triangular element
	case 3:
	{
		//3 - Rectangular element
		LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Invalid 2D element in 3D domain.\n\n\n");
		throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Invalid 2D Element\n\n\n");
		break;
	}
	case 4:
	{
		//4 - Tetrahedral Element

		// Number of nodes for incidence
		Param = (Order + 1) * (Order + 2) * (Order + 3) / 6;

		switch (Order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
			case  1: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP< 1>>(v_Mat[Material])); break; }
			case  4: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP< 4>>(v_Mat[Material])); break; }
			case 10: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP<10>>(v_Mat[Material])); break; }
			case 11: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP<11>>(v_Mat[Material])); break; }
			case 14: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP<14>>(v_Mat[Material])); break; }
			case 15: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP<15>>(v_Mat[Material])); break; }
			case 24: { v_Elem.emplace_back(std::make_shared<Elem_Tet4_IP<24>>(v_Mat[Material])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (10 nodes)
			switch (numIP) {
			case  4: { v_Elem.emplace_back(std::make_shared<Elem_Tet10_IP< 4>>(v_Mat[Material])); break; }
			case 10: { v_Elem.emplace_back(std::make_shared<Elem_Tet10_IP<10>>(v_Mat[Material])); break; }
			case 11: { v_Elem.emplace_back(std::make_shared<Elem_Tet10_IP<11>>(v_Mat[Material])); break; }
			case 14: { v_Elem.emplace_back(std::make_shared<Elem_Tet10_IP<14>>(v_Mat[Material])); break; }
			case 15: { v_Elem.emplace_back(std::make_shared<Elem_Tet10_IP<15>>(v_Mat[Material])); break; }
			case 24: { v_Elem.emplace_back(std::make_shared<Elem_Tet10_IP<24>>(v_Mat[Material])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (20 nodes)
			switch (numIP) {
			case 10: { v_Elem.emplace_back(std::make_shared<Elem_Tet20_IP<10>>(v_Mat[Material])); break; }
			case 11: { v_Elem.emplace_back(std::make_shared<Elem_Tet20_IP<11>>(v_Mat[Material])); break; }
			case 14: { v_Elem.emplace_back(std::make_shared<Elem_Tet20_IP<14>>(v_Mat[Material])); break; }
			case 15: { v_Elem.emplace_back(std::make_shared<Elem_Tet20_IP<15>>(v_Mat[Material])); break; }
			case 24: { v_Elem.emplace_back(std::make_shared<Elem_Tet20_IP<24>>(v_Mat[Material])); break; }
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
	case 5:
	{
		//5 - Hexahedral Element

		// Number of nodes for incidence
		Param = (Order + 1) * (Order + 1) * (Order + 1);

		switch (Order)
		{
		case 1:
		{
			//1	- Linear interpolation (8 nodes)
			switch (numIP) {
			case  8: { v_Elem.emplace_back(std::make_shared<Elem_Hex8_IP< 8>>(v_Mat[Material])); break; }
			case 27: { v_Elem.emplace_back(std::make_shared<Elem_Hex8_IP<27>>(v_Mat[Material])); break; }
			case 64: { v_Elem.emplace_back(std::make_shared<Elem_Hex8_IP<64>>(v_Mat[Material])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (27 nodes)
			switch (numIP) {
			case 27: { v_Elem.emplace_back(std::make_shared<Elem_Hex27_IP<27>>(v_Mat[Material])); break; }
			case 64: { v_Elem.emplace_back(std::make_shared<Elem_Hex27_IP<64>>(v_Mat[Material])); break; }
			default:
				LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error. Wrong number of integration points.\n\n\n");
				throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n");
				break;
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (64 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Hex64>(v_Mat[Material]));
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
	case 6:
	{
		//6 - Prism Element

		// Number of nodes for incidence
		Param = ((Order + 1) * (Order + 2) / 2) * (Order + 1);

		switch (Order)
		{
		case 1:
		{
			//1	- Linear interpolation (6 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Pri6>(v_Mat[Material]));
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (18 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Pri18>(v_Mat[Material]));
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (40 nodes)
			v_Elem.emplace_back(std::make_shared<Elem_Pri40>(v_Mat[Material]));
			break;
		}
		case 4:
		{
			//4	- Cubic / Linear interpolation (20 nodes)
			Param = 20;
			v_Elem.emplace_back(std::make_shared<Elem_Pri20>(v_Mat[Material]));
			break;
		}
		case 5:
		{
			//5	- Cubic / Quadratic interpolation (30 nodes)
			Param = 30;
			v_Elem.emplace_back(std::make_shared<Elem_Pri30>(v_Mat[Material]));
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
		LOG("\n\nDomain.addElement: Element " << std::to_string(index) << " creation error.\nCheck type: " << std::to_string(Type) << ", order: " << std::to_string(Order) << " and number of IP: " << std::to_string(numIP) << "\n\n\n");
		throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Invalid Element\n\n\n");
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
