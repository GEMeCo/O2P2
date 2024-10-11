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
#include "OutputSystem.h"

// ================================================================================================
//
// Explicit template instantiation
//
// ================================================================================================
template class O2P2::Post::OutputSystem<2>;
template class O2P2::Post::OutputSystem<3>;

template void O2P2::Post::OutputSystem<2>::draw_AcadView_Node(std::ofstream& file, O2P2::Geom::Domain<2>* theDomain, O2P2::Post::PostProcess* thePost);
template void O2P2::Post::OutputSystem<3>::draw_AcadView_Node(std::ofstream& file, O2P2::Geom::Domain<3>* theDomain, O2P2::Post::PostProcess* thePost);

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): draw_AcadView_Node
// 
// ================================================================================================
template<int nDim>
void O2P2::Post::OutputSystem<nDim>::draw_AcadView_Node(std::ofstream& file, O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	file << "Project Description.\n";

	// Acadview uses # flag to begin a new section / list
	file << "#\n";

	// As only plane elements are printed, the number of faces must be evaluated
	size_t nFaces = 0;
	for (auto& elem : theDomain->getElem()) {
		nFaces += elem->getNumFaces();
	}

	// Add Truss elements
	nFaces += theDomain->getNumBars();

	// Add inclusion elements
	for (auto& elem : theDomain->getIncl()) {
		nFaces += elem->getNumFaces();
	}

	// Add immersed fibers
	nFaces += theDomain->getNumFibs();

	// Add active face elements
	for (auto& elem : theDomain->getFaces()) {
		nFaces += elem->getNumFaces();
	}

	// Total number of nodes and elements
	file << theDomain->getNumNodes() + theDomain->getNumPoints() << " "
		<< nFaces << " " << thePost->mv_SolOnNode.size() * nDim << "\n#\n";

	for (auto& node : theDomain->getNode()) {
		file << *node;

		for (int k = nDim; k < 6; k++) {
			// format is a function defined in common.h
			file << formatFixed << 0.F;
		}
		file << "\n";
	}

	// Add points
	for (auto& node : theDomain->getPoint()) {
		file << *node;

		for (int k = nDim; k < 6; k++) {
			// format is a function defined in common.h
			file << formatFixed << 0.F;
		}
		file << "\n";
	}

	file << "#\n";

	// Element conectivity
	for (auto& elem : theDomain->getElem()) {
		file << elem->printByIndex_AV(1);
	}

	// Truss conectivity
	for (auto& elem : theDomain->getBars()) {
		file << elem->printByIndex_AV(1);
	}

	// Inclusion conectivity
	for (auto& elem : theDomain->getIncl()) {
		file << elem->printByIndex_AV(theDomain->getNumNodes() + 1);
	}

	// Immersed fibers conectivity
	for (auto& elem : theDomain->getFibs()) {
		file << elem->printByIndex_AV(theDomain->getNumNodes() + 1);
	}

	// Active face conectivity
	for (auto& elem : theDomain->getFaces()) {
		file << elem->printByIndex_AV(1);
	}

	for (int i = 0; i < thePost->mv_SolOnNode.size(); ++i) {
		for (int j = 0; j < nDim; ++j) {
			file << "#\n" << "Desl_" << j << "_t_" << std::defaultfloat << std::get<0>(thePost->mv_SolOnNode[i]) << "\n";

			// Solution on nodes
			for (int k = 0; k < theDomain->getNumNodes(); ++k) {
				// Must write four data here -> disp x, y, z, "color" - value to be ploted
				const auto& node = theDomain->getNode(k);
				const auto& x = node->getInitPos();

				for (int l = 0; l < nDim; l++) {
					file << formatScien << std::get<1>(thePost->mv_SolOnNode[i])[k * nDim + l] - x[l];
				}
				if (nDim == 2) file << formatScien << 0.F;
				file << formatScien << std::get<1>(thePost->mv_SolOnNode[i])[k * nDim + j] - x[j] << "\n";
			}

			// Solution on points
			for (int l = 0; l < theDomain->getNumPoints(); ++l) {
				const auto& point = theDomain->getPoint(l);
				const auto& x = point->getInitPos();

				// Must write four data here -> disp x, y, z, color
				for (int k = 0; k < nDim; ++k) {
					file << formatScien << std::get<1>(thePost->mv_SolOnNode[i])[theDomain->getNumNodes() * nDim + l * nDim + k] - x[k];
				}
				if (nDim == 2) file << formatScien << 0.F;
				file << formatScien << std::get<1>(thePost->mv_SolOnNode[i])[theDomain->getNumNodes() * nDim + l * nDim + j] - x[j] << "\n";
			}
		}
	}
	file.close();
}


// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): draw_AcadView_Elem
// 
// In this output, Each and every element are independent, thus they have their own nodes.
// It works for element information (such as stresses), but only displacement has been made so far.
// Thus, it was kept only for future reference.
// 
// ================================================================================================
template<int nDim>
void O2P2::Post::OutputSystem<nDim>::draw_AcadView_Elem(std::ofstream& file, O2P2::Geom::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	file << "Project Description.\n";

	// Acadview uses # flag to begin a new section / list
	file << "#\n";

	// Total number of nodes
	int iNd = 0;
	for (auto& elem : theDomain->getElem()) {
		iNd += elem->getNumNodes();
	}
	file << iNd << " " << theDomain->getNumElems() << " " << thePost->mv_SolOnNode.size() * nDim << "\n#\n";

	for (auto& elem : theDomain->getElem()) {
		for (auto& node : elem->getConectivity()) {
			file << *node;

			for (int k = nDim; k < 6; k++) {
				// format is a function defined in common.h
				file << formatFixed << 0.F;
			}
			file << "\n";
		}
	}

	file << "#\n";

	size_t iS = 0;
	for (auto& elem : theDomain->getElem()) {
		file << elem->printByAdder_AV(iS);
		iS += elem->getNumNodes();
	}

	for (int i = 0; i < thePost->mv_SolOnNode.size(); ++i) {
		for (int j = 0; j < nDim; ++j) {
			file << "#\n" << "Desl_" << j << "_t_" << std::defaultfloat << std::get<0>(thePost->mv_SolOnNode[i]) << "\n";

			for (auto& elem : theDomain->getElem()) {
				for (auto& node : elem->getConectivity()) {
					size_t pos = node->mv_index * nDim;

					// Must write four data here -> disp x, y, z, color
					for (int k = 0; k < nDim; ++k) {
						file << formatScien << std::get<1>(thePost->mv_SolOnNode[i])[pos + k];
					}
					if (nDim == 2) file << formatScien << 0.F;

					file << formatScien << std::get<1>(thePost->mv_SolOnNode[i])[pos + j] << "\n";
				}
			}
		}
	}
	file.close();
}
