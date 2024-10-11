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
#include "Domain.h"
#include "Elements/Elem_Lin.h"

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
template void O2P2::Geom::Domain<2>::addNode(const size_t& index, const std::array<double, 2>& x0);
template void O2P2::Geom::Domain<3>::addNode(const size_t& index, const std::array<double, 3>& x0);

template void O2P2::Geom::Domain<2>::addPoint(const size_t& index, const std::array<double, 2>& x0);
template void O2P2::Geom::Domain<3>::addPoint(const size_t& index, const std::array<double, 3>& x0);

template void O2P2::Geom::Domain<2>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);
template void O2P2::Geom::Domain<3>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam);

template void O2P2::Geom::Domain<2>::addSection(const int& type, const double& value, const int& PSType);
template void O2P2::Geom::Domain<3>::addSection(const int& type, const double& value, const int& PSType);

template void O2P2::Geom::Domain<2>::addElementConectivity(const size_t& index, const int& Type, const std::vector<size_t>& conectivity);
template void O2P2::Geom::Domain<3>::addElementConectivity(const size_t& index, const int& Type, const std::vector<size_t>& conectivity);

template void O2P2::Geom::Domain<2>::addInclusionConectivity(const size_t& index, const int& Type, const std::vector<size_t>& conectivity);
template void O2P2::Geom::Domain<3>::addInclusionConectivity(const size_t& index, const int& Type, const std::vector<size_t>& conectivity);

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addNode
//
// ================================================================================================
template<int nDim>
void O2P2::Geom::Domain<nDim>::addNode(const size_t& index, const std::array<double, nDim>& x0)
{
	mv_Node.emplace_back(std::make_shared<DomainNode<nDim>>(index, x0));
	LOG("Domain.addNode: Adding Node: " << std::to_string(index) << *mv_Node.back());
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addPoint
//
// ================================================================================================
template<int nDim>
void O2P2::Geom::Domain<nDim>::addPoint(const size_t& index, const std::array<double, nDim>& x0)
{
	mv_Point.emplace_back(std::make_shared<DomainNode<nDim>>(index, x0));
	LOG("Domain.addPoint: Adding Point: " << std::to_string(index) << *mv_Point.back());
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addMaterial
//
// ================================================================================================
template<int nDim>
void O2P2::Geom::Domain<nDim>::addMaterial(const size_t& index, const MaterialType& matType, const std::vector<double> matParam)
{
	LOG("Domain.addMaterial: Adding Material " << MaterialTypeNames[matType] << " with index " << std::to_string(index));

	switch (matType)
	{
	case MaterialType::SVK_ISO:
	{
		// 1 - Elastic, Isotropic, Saint-Venant-Kirchhoff
		mv_Mat.emplace_back(std::make_shared<Mat_SVK_ISO>(index, matParam));
		assert(std::dynamic_pointer_cast<O2P2::Geom::Mat_SVK_ISO>(mv_Mat.back()));
		LOG("Domain.addMaterial: added Material: " << std::to_string(index) << *(std::static_pointer_cast<O2P2::Geom::Mat_SVK_ISO>(mv_Mat.back())));
		break;
	}
	case MaterialType::SVK_ORT:
	case MaterialType::SVK_DAMAGE:
	case MaterialType::CONCRETE:
	case MaterialType::STEEL:
	case MaterialType::RSD_ISO:
	case MaterialType::RSD_PLASTIC:
	case MaterialType::TEMP_ISO_LIN:
	case MaterialType::TEMP_ISO_NLIN:
	case MaterialType::RSD_TP:
	default:
		LOG("\n\nDomain.addMaterial: Undefined material. Check input file.\n\n\n");
		throw std::invalid_argument("\n\n\n******** Error in material definition, check input ********\n\n\n");
		break;
	}
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addCrossSection
//
// ================================================================================================
template<int nDim>
void O2P2::Geom::Domain<nDim>::addSection(const int& type, const double& value, const int& PSType)
{
	// For now, only two options are provided
	// 1 - Cross section of linear elements (trusses)
	// 2 - Thickness of plane elements (membrane)
	if (type == 1) {
		LOG("Domain.addSection: Bar Cross Section: " << std::to_string(value));
		mv_Sect.emplace_back(std::make_shared<CrossSection>(value));
	}
	else {
		LOG("Domain.addSection: Plane Thickness: " << std::to_string(value) << "; " << PlaneStateTypeNames[PlaneStateType(PSType)]);
		mv_Sect.emplace_back(std::make_shared<PlaneSection>(PlaneStateType(PSType), value));
	}
}


// ================================================================================================
//
// Specialized implementation of Template Member Function (only 2D): addElement
//
// ================================================================================================
template<>
int O2P2::Geom::Domain<2>::addElement(const size_t& index, const int& type, const int& order, const int& numIP, const size_t& material, const size_t& section)
{
	// Return value - number of nodes in incidence
	int mi_param = 0;

	// The type of element defines which class Element
	switch (type)
	{
	case 1:
	{
		// 1 - Truss element
		mi_param = order + 1;
		mv_numElem--;
		mv_numBars++;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (2 nodes)
			switch (numIP) {
				case 2: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 2, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 2, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 2, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (3 nodes)
			switch (numIP) {
				case 2: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 3, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 3, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 3, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (4 nodes)
			switch (numIP) {
				case 2: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 4, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 4, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 4, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided.\n\n\n");
			break;
		}

		break;
	}
	case 2:
	{
		//2 - Triangular element
		assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
		auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 2) / 2;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (3 nodes)
			switch (numIP) {
				case  1: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 1>>(mv_Mat[material], mi_Sect)); break; }
				case  3: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 3>>(mv_Mat[material], mi_Sect)); break; }
				case  4: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  6: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (6 nodes)
			switch (numIP) {
				case  3: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 3>>(mv_Mat[material], mi_Sect)); break; }
				case  4: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  6: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (10 nodes)
			switch (numIP) {
				case  6: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}

		break;
	}
	case 3:
	{
		//3 - Rectangular element
		assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
		auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

		// Number of nodes for incidence
		if (order < 4) {
			mi_param = (order + 1) * (order + 1);
		}
		else if (order < 6) {
			mi_param = (4) * (order - 2);
		}
		else mi_param = 20;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
				case  4: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  9: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (9 nodes)
			switch (numIP) {
				case  9: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (16 nodes)
			switch (numIP) {
				case  9: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 4:
		{
			//4	- Cubic / linear interpolation (8 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect8>(mv_Mat[material], mi_Sect));
			break;
		}
		case 5:
		{
			//5 - Cubic / quadratic interpolation (12 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect12>(mv_Mat[material], mi_Sect));
			break;
		}
		case 6:
		{
			//6	- Cubic / quartic interpolation (20 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect20>(mv_Mat[material], mi_Sect));
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	default:
		throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Invalid Element Type\n\n\n");
		break;
	}

	return mi_param;
}


// ================================================================================================
//
// Specialized implementation of Template Member Function (only 3D): addElement
//
// ================================================================================================
template<>
int O2P2::Geom::Domain<3>::addElement(const size_t& index, const int& type, const int& order, const int& numIP, const size_t& material, const size_t& section)
{
	int mi_param = 0;

	// The type of element defines which class Element
	switch (type)
	{
	case 1:
	{
		//1 - Bar / Truss element

		// Number of nodes for incidence
		mi_param = (order + 1);
		mv_numElem--;
		mv_numBars++;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (2 nodes)
			switch (numIP) {
				case 2: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 2, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 2, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 2, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (3 nodes)
			switch (numIP) {
				case 2: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 3, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 3, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 3, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (4 nodes)
			switch (numIP) {
				case 2: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 4, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 4, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Bars.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 4, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided.\n\n\n");
			break;
		}
		break;
	}
	case 2:
		//2 - Triangular element
	case 3:
	{
		//3 - Rectangular element
		throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Invalid 2D Element in 3D environment. Don't you mean an immersed element instead?\n\n\n");
		break;
	}
	case 4:
	{
		//4 - Tetrahedral Element

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 2) * (order + 3) / 6;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
				case  1: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP< 1>>(mv_Mat[material])); break; }
				case  4: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP< 4>>(mv_Mat[material])); break; }
				case 10: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<10>>(mv_Mat[material])); break; }
				case 11: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<11>>(mv_Mat[material])); break; }
				case 14: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<14>>(mv_Mat[material])); break; }
				case 15: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<15>>(mv_Mat[material])); break; }
				case 24: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<24>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (10 nodes)
			switch (numIP) {
				case  4: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP< 4>>(mv_Mat[material])); break; }
				case 10: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<10>>(mv_Mat[material])); break; }
				case 11: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<11>>(mv_Mat[material])); break; }
				case 14: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<14>>(mv_Mat[material])); break; }
				case 15: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<15>>(mv_Mat[material])); break; }
				case 24: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<24>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (20 nodes)
			switch (numIP) {
				case 10: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<10>>(mv_Mat[material])); break; }
				case 11: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<11>>(mv_Mat[material])); break; }
				case 14: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<14>>(mv_Mat[material])); break; }
				case 15: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<15>>(mv_Mat[material])); break; }
				case 24: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<24>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	case 5:
	{
		//5 - Hexahedral Element

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 1) * (order + 1);

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (8 nodes)
			switch (numIP) {
				case  8: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex8_IP< 8>>(mv_Mat[material])); break; }
				case 27: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex8_IP<27>>(mv_Mat[material])); break; }
				case 64: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex8_IP<64>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (27 nodes)
			switch (numIP) {
				case 27: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex27_IP<27>>(mv_Mat[material])); break; }
				case 64: { mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex27_IP<64>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (64 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex64>(mv_Mat[material]));
			break;
		}
		default:
			// Higher order
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	case 6:
	{
		//6 - Prism Element

		// Number of nodes for incidence
		mi_param = ((order + 1) * (order + 2) / 2) * (order + 1);

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (6 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri6>(mv_Mat[material]));
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (18 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri18>(mv_Mat[material]));
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (40 nodes)
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri40>(mv_Mat[material]));
			break;
		}
		case 4:
		{
			//4	- Cubic / Linear interpolation (20 nodes)
			mi_param = 20;
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri20>(mv_Mat[material]));
			break;
		}
		case 5:
		{
			//5	- Cubic / Quadratic interpolation (30 nodes)
			mi_param = 30;
			mv_Elem.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri30>(mv_Mat[material]));
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	default:
		throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Invalid Element\n\n\n");
		break;
	}

	return mi_param;
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addElementConectivity
//
// ================================================================================================
template<int nDim>
void O2P2::Geom::Domain<nDim>::addElementConectivity(const size_t& index, const int& Type, const std::vector<size_t>& conectivity)
{
	std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(conectivity.size());

	for (int i = 0; i < conectivity.size(); i++)
	{
		mi_conect[i] = mv_Node.at(conectivity[i]);

		// Whenever a element points to a node, it is added to the inverse indexing
		mi_conect[i]->addConectToElem(index);
	}
	if (Type == 1) {
		mv_Bars.back()->setConectivity(mi_conect);

		LOG("Domain.addElement: Element index " << std::to_string(index) << " incidence: " << *mv_Bars.back());
		LOG("Domain.addElement: " << mv_Bars.back()->printGeom());
	}
	else {
		mv_Elem.back()->setConectivity(mi_conect);

		LOG("Domain.addElement: Element index " << std::to_string(index) << " incidence: " << *mv_Elem.back());
		LOG("Domain.addElement: " << mv_Elem.back()->printGeom());
	}
}


// ================================================================================================
//
// Specialized implementation of Template Member Function (only 2D): addInclusion
//
// ================================================================================================
template<>
int O2P2::Geom::Domain<2>::addInclusion(const size_t& index, const int& type, const int& order, const int& numIP, const size_t& material, const size_t& section)
{
	// Return value - number of nodes in incidence
	int mi_param = 0;

	// The type of element defines which class Element
	switch (type)
	{
	case 1:
	{
		// 1 - Truss element
		mi_param = order + 1;
		mv_numIncl--;
		mv_numFibs++;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (2 nodes)
			switch (numIP) {
				case 2: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 2, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 2, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 2, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (3 nodes)
			switch (numIP) {
				case 2: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 3, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 3, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 3, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (4 nodes)
			switch (numIP) {
				case 2: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 4, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 4, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<2, 4, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided.\n\n\n");
			break;
		}

		break;
	}
	case 2:
	{
		//2 - Triangular element
		assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
		auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 2) / 2;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (3 nodes)
			switch (numIP) {
				case  1: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 1>>(mv_Mat[material], mi_Sect)); break; }
				case  3: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 3>>(mv_Mat[material], mi_Sect)); break; }
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  6: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (6 nodes)
			switch (numIP) {
				case  3: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 3>>(mv_Mat[material], mi_Sect)); break; }
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  6: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (10 nodes)
			switch (numIP) {
				case  6: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}

		break;
	}
	case 3:
	{
		//3 - Rectangular element
		assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
		auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

		// Number of nodes for incidence
		if (order < 4) {
			mi_param = (order + 1) * (order + 1);
		}
		else if (order < 6) {
			mi_param = (4) * (order - 2);
		}
		else mi_param = 20;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  9: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (9 nodes)
			switch (numIP) {
				case  9: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (16 nodes)
			switch (numIP) {
				case  9: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 4:
		{
			//4	- Cubic / linear interpolation (8 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect8>(mv_Mat[material], mi_Sect));
			break;
		}
		case 5:
		{
			//5 - Cubic / quadratic interpolation (12 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect12>(mv_Mat[material], mi_Sect));
			break;
		}
		case 6:
		{
			//6	- Cubic / quartic interpolation (20 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect20>(mv_Mat[material], mi_Sect));
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	default:
		throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Invalid Element Type\n\n\n");
		break;
	}

	return mi_param;
}

// ================================================================================================
//
// Specialized implementation of Template Member Function (only 3D): addInclusion
//
// ================================================================================================
template<>
int O2P2::Geom::Domain<3>::addInclusion(const size_t& index, const int& type, const int& order, const int& numIP, const size_t& material, const size_t& section)
{
	int mi_param = 0;

	// The type of element defines which class Element
	switch (type)
	{
	case 1:
	{
		//1 - Bar / Truss element

		// Number of nodes for incidence
		mi_param = (order + 1);
		mv_numIncl--;
		mv_numFibs++;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (2 nodes)
			switch (numIP) {
				case 2: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 2, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 2, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 2, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (3 nodes)
			switch (numIP) {
				case 2: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 3, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 3, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 3, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (4 nodes)
			switch (numIP) {
				case 2: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 4, 2>>(mv_Mat[material], mv_Sect[section])); break; }
				case 3: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 4, 3>>(mv_Mat[material], mv_Sect[section])); break; }
				case 4: { mv_Fibs.push_back(std::make_shared<O2P2::Geom::Elem::Elem_Lin<3, 4, 4>>(mv_Mat[material], mv_Sect[section])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided.\n\n\n");
			break;
		}
		break;
	}
	case 2:
	{
		//2 - Triangular element
		assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
		auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 2) / 2;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (3 nodes)
			switch (numIP) {
				case  1: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 1>>(mv_Mat[material], mi_Sect)); break; }
				case  3: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 3>>(mv_Mat[material], mi_Sect)); break; }
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  6: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (6 nodes)
			switch (numIP) {
				case  3: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 3>>(mv_Mat[material], mi_Sect)); break; }
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  6: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (10 nodes)
			switch (numIP) {
				case  6: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP< 6>>(mv_Mat[material], mi_Sect)); break; }
				case  7: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP< 7>>(mv_Mat[material], mi_Sect)); break; }
				case 12: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect)); break; }
				case 13: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<13>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nElement " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}

		break;
	}
	case 3:
	{
		//3 - Rectangular element
		assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
		auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

		// Number of nodes for incidence
		if (order < 4) {
			mi_param = (order + 1) * (order + 1);
		}
		else if (order < 6) {
			mi_param = (4) * (order - 2);
		}
		else mi_param = 20;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 4>>(mv_Mat[material], mi_Sect)); break; }
				case  9: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (9 nodes)
			switch (numIP) {
				case  9: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (16 nodes)
			switch (numIP) {
				case  9: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP< 9>>(mv_Mat[material], mi_Sect)); break; }
				case 16: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP<16>>(mv_Mat[material], mi_Sect)); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 4:
		{
			//4	- Cubic / linear interpolation (8 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect8>(mv_Mat[material], mi_Sect));
			break;
		}
		case 5:
		{
			//5 - Cubic / quadratic interpolation (12 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect12>(mv_Mat[material], mi_Sect));
			break;
		}
		case 6:
		{
			//6	- Cubic / quartic interpolation (20 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect20>(mv_Mat[material], mi_Sect));
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	case 4:
	{
		//4 - Tetrahedral Element

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 2) * (order + 3) / 6;

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (4 nodes)
			switch (numIP) {
				case  1: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP< 1>>(mv_Mat[material])); break; }
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP< 4>>(mv_Mat[material])); break; }
				case 10: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<10>>(mv_Mat[material])); break; }
				case 11: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<11>>(mv_Mat[material])); break; }
				case 14: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<14>>(mv_Mat[material])); break; }
				case 15: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<15>>(mv_Mat[material])); break; }
				case 24: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet4_IP<24>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (10 nodes)
			switch (numIP) {
				case  4: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP< 4>>(mv_Mat[material])); break; }
				case 10: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<10>>(mv_Mat[material])); break; }
				case 11: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<11>>(mv_Mat[material])); break; }
				case 14: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<14>>(mv_Mat[material])); break; }
				case 15: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<15>>(mv_Mat[material])); break; }
				case 24: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet10_IP<24>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (20 nodes)
			switch (numIP) {
				case 10: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<10>>(mv_Mat[material])); break; }
				case 11: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<11>>(mv_Mat[material])); break; }
				case 14: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<14>>(mv_Mat[material])); break; }
				case 15: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<15>>(mv_Mat[material])); break; }
				case 24: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tet20_IP<24>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusoin element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	case 5:
	{
		//5 - Hexahedral Element

		// Number of nodes for incidence
		mi_param = (order + 1) * (order + 1) * (order + 1);

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (8 nodes)
			switch (numIP) {
				case  8: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex8_IP< 8>>(mv_Mat[material])); break; }
				case 27: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex8_IP<27>>(mv_Mat[material])); break; }
				case 64: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex8_IP<64>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (27 nodes)
			switch (numIP) {
				case 27: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex27_IP<27>>(mv_Mat[material])); break; }
				case 64: { mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex27_IP<64>>(mv_Mat[material])); break; }
				default: { throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Wrong number of integration points.\n\n\n"); break; }
			}
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (64 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Hex64>(mv_Mat[material]));
			break;
		}
		default:
			// Higher order
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	case 6:
	{
		//6 - Prism Element

		// Number of nodes for incidence
		mi_param = ((order + 1) * (order + 2) / 2) * (order + 1);

		switch (order)
		{
		case 1:
		{
			//1	- Linear interpolation (6 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri6>(mv_Mat[material]));
			break;
		}
		case 2:
		{
			//2	- Quadratic interpolation (18 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri18>(mv_Mat[material]));
			break;
		}
		case 3:
		{
			//3	- Cubic interpolation (40 nodes)
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri40>(mv_Mat[material]));
			break;
		}
		case 4:
		{
			//4	- Cubic / Linear interpolation (20 nodes)
			mi_param = 20;
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri20>(mv_Mat[material]));
			break;
		}
		case 5:
		{
			//5	- Cubic / Quadratic interpolation (30 nodes)
			mi_param = 30;
			mv_Incl.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Pri30>(mv_Mat[material]));
			break;
		}
		default:
			throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Only up to cubic order is provided\n\n\n");
			break;
		}
		break;
	}
	default:
		throw std::invalid_argument("\n\n\nInclusion element " + std::to_string(index) + " creation error. Invalid Element\n\n\n");
		break;
	}

	return mi_param;
}

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): addInclusionConectivity
//
// ================================================================================================
template<int nDim>
void O2P2::Geom::Domain<nDim>::addInclusionConectivity(const size_t& index, const int& Type, const std::vector<size_t>& Conectivity)
{
	std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(Conectivity.size());

	for (int i = 0; i < Conectivity.size(); i++)
	{
		mi_conect[i] = mv_Point.at(Conectivity[i]);
	}

	if (Type == 1) {
		mv_Fibs.back()->setConectivity(mi_conect);

		LOG("Domain.addInclusion: Element index " << std::to_string(index) << " incidence: " << *mv_Fibs.back());
		LOG("Domain.addInclusion: " << mv_Fibs.back()->printGeom());
	}
	else {
		mv_Incl.back()->setConectivity(mi_conect);

		LOG("Domain.addInclusion: Element index " << std::to_string(index) << " incidence: " << *mv_Incl.back());
		LOG("Domain.addInclusion: " << mv_Incl.back()->printGeom());
	}
}


// ================================================================================================
//
// Specialized implementation of Template Member Function (only 2D): addFaceElem
//
// ================================================================================================
template<>
void O2P2::Geom::Domain<2>::addFaceElem(std::shared_ptr<O2P2::Geom::Elem::Element> element, const int& face, const size_t& material, const size_t& section)
{
	// This function cannot be called. The prism is always 3D
	throw std::invalid_argument("\n\n\nDomain.addFaceElem: How do you get in here?\n\n\n");
}

// ================================================================================================
//
// Specialized implementation of Template Member Function (only 3D): addFaceElem
//
// ================================================================================================
template<>
void O2P2::Geom::Domain<3>::addFaceElem(std::shared_ptr<O2P2::Geom::Elem::Element> element, const int& face, const size_t& material, const size_t& section)
{
	assert(std::dynamic_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]));
	auto mi_Sect = std::static_pointer_cast<O2P2::Geom::PlaneSection>(mv_Sect[section]);

	// Active face element is created upon number of nodes of the prismatic element and the selected face
	switch (element->getNumNodes()) {
	case 6:
	{
		switch (face)
		{
		case 0:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 3>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 3> mi_conectivity = { 0, 2, 1 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 1:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri3_IP< 3>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 3> mi_conectivity = { 3, 4, 5 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 2:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 4>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 4> mi_conectivity = { 0, 1, 3, 4 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 3:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 4>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 4> mi_conectivity = { 1, 2, 4, 5 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 4:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect4_IP< 4>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 4> mi_conectivity = { 2, 0, 5, 3 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		default: { throw std::invalid_argument("\n\n\nDomain.addFaceElem: Invalid input. Check face number\n\n\n"); break; }
		}
		break;
	}
	case 18:
	{
		switch (face)
		{
		case 0:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 7>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 6> mi_conectivity = { 0, 1, 2, 3, 4, 5 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 1:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri6_IP< 7>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 6> mi_conectivity = { 12, 13, 14, 15, 16, 17 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 2:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP< 9>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 9> mi_conectivity = { 0, 1, 2, 6, 7, 8, 12, 13, 14 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 3:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP< 9>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 9> mi_conectivity = { 2, 4, 5, 8, 10, 11, 14, 16, 17 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 4:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect9_IP< 9>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 9> mi_conectivity = { 5, 3, 0, 11, 9, 6, 17, 15, 12 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		default: { throw std::invalid_argument("\n\n\nDomain.addFaceElem: Invalid input. Check face number\n\n\n"); break; }
		}
		break;
	}
	case 20:
	{
		switch (face)
		{
		case 0:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 10> mi_conectivity = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 1:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 10> mi_conectivity = { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 2:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect8>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 8> mi_conectivity = { 0, 1, 2, 3, 10, 11, 12, 13 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 3:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect8>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 8> mi_conectivity = { 3, 6, 8, 9, 13, 16, 18, 19 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 4:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect8>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 8> mi_conectivity = { 9, 7, 4, 0, 19, 17, 14, 10 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		default: { throw std::invalid_argument("\n\n\nDomain.addFaceElem: Invalid input. Check face number\n\n\n"); break; }
		}
		break;
	}
	case 30:
	{
		switch (face)
		{
		case 0:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 10> mi_conectivity = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 1:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 10> mi_conectivity = { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 2:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect12>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 12> mi_conectivity = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 3:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect12>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 12> mi_conectivity = { 3, 6, 8, 9, 13, 16, 18, 19, 23, 26, 28, 29 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 4:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect12>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 12> mi_conectivity = { 9, 7, 4, 0, 19, 17, 14, 10, 29, 27, 24, 20 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		default: { throw std::invalid_argument("\n\n\nDomain.addFaceElem: Invalid input. Check face number\n\n\n"); break; }
		}
		break;
	}
	case 40:
	{
		switch (face)
		{
		case 0:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 10> mi_conectivity = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 1:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Tri10_IP<12>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 10> mi_conectivity = { 30, 31, 32, 33, 34, 35, 36, 37, 38, 39 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 2:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP<16>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 16> mi_conectivity = { 0, 1, 2, 3, 10, 11, 12, 13, 20, 21, 22, 23, 30, 31, 32, 33 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 3:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP<16>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 16> mi_conectivity = { 3, 6, 8, 9, 13, 16, 18, 19, 23, 26, 28, 29, 33, 36, 38, 39 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		case 4:
		{
			mv_Face.emplace_back(std::make_shared<O2P2::Geom::Elem::Elem_Rect16_IP<16>>(mv_Mat[material], mi_Sect));
			std::vector<std::shared_ptr<O2P2::Geom::Node>> mi_conect(mv_Face.back()->getNumNodes());
			std::array<int, 16> mi_conectivity = { 9, 7, 4, 0, 19, 17, 14, 10, 29, 27, 24, 20, 39, 37, 34, 30 };
			for (int i = 0; i < mi_conectivity.size(); i++)
				mi_conect[i] = element->getConectivity().at(mi_conectivity.at(i));
			mv_Face.back()->setConectivity(mi_conect);
			LOG("Domain.addFaceElem: Face Element with incidence: " << *mv_Face.back());
			break;
		}
		default: { throw std::invalid_argument("\n\n\nDomain.addFaceElem: Invalid input. Check face number\n\n\n"); break; }
		}
		break;
	}
	default: { throw std::invalid_argument("\n\n\nDomain.addFaceElem: How do you get in here?\n\n\n"); break; }
	}

	mv_numFace++;
}
