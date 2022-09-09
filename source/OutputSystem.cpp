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
// OutputSystem
// 
// Base class for writting files for post-processing, such as AcadView (OGL), Paraview (VTU), etc
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

template void O2P2::Post::OutputSystem<2>::draw_AcadView_Node(std::ofstream& file, O2P2::Prep::Domain<2>* theDomain, O2P2::Post::PostProcess* thePost);
template void O2P2::Post::OutputSystem<3>::draw_AcadView_Node(std::ofstream& file, O2P2::Prep::Domain<3>* theDomain, O2P2::Post::PostProcess* thePost);

// ================================================================================================
//
// Implementation of Template Member Function (2D and 3D): draw_AcadView_Node
// 
// ================================================================================================
template<int nDim>
void O2P2::Post::OutputSystem<nDim>::draw_AcadView_Node(std::ofstream& file, O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	file << "Project Description.\n";

	// Acadview uses # flag to begin a new section / list
	file << "#\n";

	// As only plane elements are printed, the number of faces must be evaluated
	size_t nFaces = 0;
	for (auto& elem : theDomain->getElem()) {
		nFaces += elem->getNumFaces();
	}

	// Total number of nodes and elements
	file << theDomain->m_nNodes << " "
		 << nFaces << " " << thePost->m_SolOnNode.size() * nDim << "\n#\n";

	for (auto& node : theDomain->getNode()) {
		file << *node;

		for (int k = nDim; k < 6; k++) {
			// format is a function defined in common.h
			file << formatFixed << 0.F;
		}
		file << "\n";
	}

	file << "#\n";

	// Conectivity and stuff
	for (auto& elem : theDomain->getElem()) {
		file << elem->printByIndex_AV(0);
	}

	for (int i = 0; i < thePost->m_SolOnNode.size(); ++i) {
		for (int j = 0; j < nDim; ++j) {
			file << "#\n" << "Desl_" << j << "_t_" << std::defaultfloat << std::get<0>(thePost->m_SolOnNode[i]) << "\n";

			for (auto& node : theDomain->getNode()) {
				// Must write four data here -> disp x, y, z, color
				// For now, we assume that node->v_DofIndex is made only by displacements
				for (auto dof : node->v_DofIndex) {
					file << formatScien << std::get<1>(thePost->m_SolOnNode[i])[dof];
				}
				if (node->v_DofIndex.size() == 2) file << formatScien << 0.F;
				file << formatScien << std::get<1>(thePost->m_SolOnNode[i])[node->v_DofIndex[j]] << "\n";
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
void O2P2::Post::OutputSystem<nDim>::draw_AcadView_Elem(std::ofstream& file, O2P2::Prep::Domain<nDim>* theDomain, O2P2::Post::PostProcess* thePost)
{
	file << "Project Description.\n";

	// Acadview uses # flag to begin a new section / list
	file << "#\n";

	// Total number of nodes
	int iNd = 0;
	for (auto& elem : theDomain->getElem()) {
		iNd += elem->getNumNodes();
	}
	file << iNd << " " << theDomain->m_nElem << " " << thePost->m_SolOnNode.size()*nDim << "\n#\n";
	
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

	for (int i = 0; i < thePost->m_SolOnNode.size(); ++i) {
		for (int j = 0; j < nDim; ++j) {
			file << "#\n" << "Desl_" << j << "_t_" << std::defaultfloat << std::get<0>(thePost->m_SolOnNode[i]) << "\n";

			for (auto& elem : theDomain->getElem()) {
				for (auto& node : elem->getConectivity()) {
					size_t pos = (node->m_index - 1) * nDim;

					// Must write four data here -> disp x, y, z, color
					for (int k = 0; k < nDim; ++k) {
						file << formatScien << std::get<1>(thePost->m_SolOnNode[i])[pos + k];
					}
					if (nDim == 2) file << formatScien << 0.F;

					file << formatScien << std::get<1>(thePost->m_SolOnNode[i])[pos + j] << "\n";
				}
			}
		}
	}
	file.close();
}

