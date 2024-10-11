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

// MKL library
#include "mkl_pardiso.h"

// M2S2 library
#include "M2S2/M2S2.h"

/** Call to MKL Pardiso Solver
  * @param size Total number of DOF.
  * @param csrMatrix Symmetric sparse matrix in CSR format.
  * @param RHS Vector of independent terms, the right hand side.
  * @param LHS Vector of dependent terms, the left hand side vector.
  */
inline void solve(int size, M2S2::CSR& csrMatrix, std::vector<double>& RHS, std::vector<double>& LHS)
{
	int mtype;
	// 2: real and symmetric positive definite matrix
	// -2: real and symmetric indefinite matrix
	// 11: real and nonsymmetric
	mtype = csrMatrix.mv_sym ? 2 : 11;

	int nrhs = 1;         /* Number of right hand sides. */

	int idum;             /* Integer dummy. */
	double ddum;          /* Double dummy */

	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void* pt[64];

	/* Pardiso control parameters. */
	int iparm[64];
	int maxfct, mnum, phase, error, msglvl;

	/* -------------------------------------*/
	/* .. Setup Pardiso control parameters. */
	/* -------------------------------------*/
	memset(&iparm[0], 0., 64 * sizeof(int));

	iparm[0] = 1;         /* No solver default */
	iparm[1] = 2;         /* Fill-in reordering from METIS */
	iparm[3] = 0;         /* No iterative-direct algorithm */
	iparm[4] = 0;         /* No user fill-in reducing permutation */
	iparm[5] = 0;         /* Write solution into x */
	iparm[7] = 3;         /* Max numbers of iterative refinement steps */
	iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
	iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
	iparm[13] = 0;        /* Output: Number of perturbed pivots */
	iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1;       /* Output: Mflops for LU factorization */
	iparm[19] = 0;        /* Output: Numbers of CG Iterations */
	iparm[34] = 1;        /* PARDISO use C-style indexing for ia and ja arrays */

	maxfct = 1;           /* Maximum number of numerical factorizations. */
	mnum = 1;             /* Which factorization to use. */
	msglvl = 0;           /* Do not print statistical information */
	error = 0;            /* Initialize error flag */

	/* ----------------------------------------------------------------*/
	/* .. Initialize the internal solver memory pointer. This is only  */
	/*   necessary for the FIRST call of the PARDISO solver.           */
	/* ----------------------------------------------------------------*/
	for (int i = 0; i < 64; i++)
	{
		pt[i] = 0;
	}

	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates  */
	/*    all memory that is necessary for the factorization.              */
	/* --------------------------------------------------------------------*/
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, csrMatrix.mv_value.data(), csrMatrix.mv_rowIndex.data(), csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		std::cout << "\nERROR during symbolic factorization: " << error;
		std::cout << "\nClick a button to exit!";
		std::cin.get();
		exit(error);
	}
	//std::cout << "\nReordering completed ... ";
	//std::cout << "\nNumber of nonzeros in factors = " << iparm[17];
	//std::cout << "\nNumber of factorization MFLOPS = " << iparm[18];

	/* ----------------------------*/
	/* .. Numerical factorization. */
	/* ----------------------------*/
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, csrMatrix.mv_value.data(), csrMatrix.mv_rowIndex.data(), csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0)
	{
		std::cout << "\nERROR during numerical factorization: " << error;
		std::cout << "\n\nClick a button to exit!";
		std::cin.get();
		exit(error);
	}
	//std::cout << "\nFactorization completed ... ";

	/* -----------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* -----------------------------------------------*/
	phase = 33;
	iparm[7] = 2;         /* Max numbers of iterative refinement steps. */

	/* Set right hand side to one. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, csrMatrix.mv_value.data(), csrMatrix.mv_rowIndex.data(), csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, RHS.data(), LHS.data(), &error);
	if (error != 0)
	{
		std::cout << "\nERROR during solution: " << error;
		std::cout << "\n\nClick a button to exit!";
		std::cin.get();
		exit(error);
	}
	//std::cout << "\nSolve completed ... ";

	/* --------------------------------------*/
	/* .. Termination and release of memory. */
	/* --------------------------------------*/
	phase = -1;           /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&size, &ddum, csrMatrix.mv_rowIndex.data(), csrMatrix.mv_colIndex.data(),
		&idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	// ------------------------------------------------------------------------------------
}