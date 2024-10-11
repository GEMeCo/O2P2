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
#include "Mesh Components\MeshFace.h"

// ================================================================================================
//
// Implementation of Template Member Function: setIndexing
//
// ================================================================================================
void O2P2::Proc::Comp::MeshFace::setIndexing()
{
	for (int i = 0; i < mv_conect.size(); ++i) {
		for (int j = 0; j < 3; ++j) {
			mv_dofIndex.push_back(mv_conect.at(i)->mv_dofList.at(j));
		}
	}

	// Since these are 2D elements immersed in 3D elements, elemental contribution must be resized
	mv_nDof = mv_conect.size() * 3;
	mv_elHes.resize(mv_nDof);
	mv_elFor.resize(mv_nDof);
}
