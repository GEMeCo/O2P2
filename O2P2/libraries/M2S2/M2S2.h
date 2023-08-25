// ================================================================================================
// 
// This file is part of M2S2 - Matrices for Mechanices of Solids and Structures
//
// Copyright(C) 2023 
//		Dorival Piedade Neto &
//		Rodrigo Ribeiro Paccola &
//		Rogério Carrazedo
// 
// This source code form is subject to the terms of the Apache License 2.0.
// If a copy of Apache License 2.0 was not distributed with this file, you can obtain one at
// https://www.apache.org/licenses/LICENSE-2.0
// 
// ================================================================================================
#pragma once

// M2S2 libraries
#include "Dyadic2S.h"
#include "Dyadic2N.h"
#include "Dyadic4S.h"
#include "Dyadic4C.h"
#include "MatrixS.h"
#include "MatrixX.h"
#include "MatrixSparse.h"

namespace M2S2 {
	// ================================================================================================
	//
	// Helper functions
	//
	// ================================================================================================
	//
	// Operators for Symmetric and Asymmetric dyadics
	//
	// ================================================================================================
	//
	// Symmetric + Asymmetric dyadics
	//
	// ================================================================================================
	inline Dyadic2N operator+(const Dyadic2S& first, const Dyadic2N& second)
	{
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < result.cols(); ++j) {
				result.at(i, j) = first.at(i, j) + second.at(i, j);
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Symmetric - Asymmetric dyadics
	//
	// ================================================================================================
	inline Dyadic2N operator-(const Dyadic2S& first, const Dyadic2N& second)
	{
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < result.cols(); ++j) {
				result.at(i, j) = first.at(i, j) - second.at(i, j);
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Symmetric * Asymmetric dyadics
	//
	// ================================================================================================
	inline Dyadic2N operator*(const Dyadic2S& first, const Dyadic2N& second)
	{
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < first.cols(); ++j) {
				for (unsigned int k = 0; k < second.rows(); ++k) {
					result.at(i, j) = first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	// ================================================================================================
	//
	// Asymmetric + Symmetric dyadics
	//
	// ================================================================================================
	inline Dyadic2N operator+(const Dyadic2N& first, const Dyadic2S& second)
	{
		return (second + first);
	}

	// ================================================================================================
	//
	// Asymmetric - Symmetric dyadics
	//
	// ================================================================================================
	inline Dyadic2N operator-(const Dyadic2N& first, const Dyadic2S& second)
	{
		return (second - first) * (-1.);
	}

	// ================================================================================================
	//
	// Asymmetric * Symmetric dyadics
	//
	// ================================================================================================
	inline Dyadic2N operator*(const Dyadic2N& first, const Dyadic2S& second)
	{
		assert(first.rows() == second.rows());		// Size of dyadics does not correspond!

		Dyadic2N result(second.rows());
		for (unsigned int i = 0; i < result.rows(); ++i) {
			for (unsigned int j = 0; j < first.cols(); ++j) {
				for (unsigned int k = 0; k < second.rows(); ++k) {
					result.at(i, j) = first.at(i, k) * second.at(k, j);
				}
			}
		}
		return result;
	}

	inline Dyadic2N operator+(const double& alfa, const Dyadic2N& dyadic)
	{
		return dyadic + alfa;
	}

	inline Dyadic2N operator-(const double& alfa, const Dyadic2N& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	inline Dyadic2N operator*(const double& alfa, const Dyadic2N& dyadic)
	{
		return dyadic * alfa;
	}

	inline Dyadic2S operator+(const double& alfa, const Dyadic2S& dyadic)
	{
		return dyadic + alfa;
	}

	inline Dyadic2S operator-(const double& alfa, const Dyadic2S& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	inline Dyadic2S operator*(const double& alfa, const Dyadic2S& dyadic)
	{
		return dyadic * alfa;
	}

	inline MatrixX operator+(const double& alfa, const MatrixX& dyadic)
	{
		return dyadic + alfa;
	}

	inline MatrixX operator-(const double& alfa, const MatrixX& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	inline MatrixX operator*(const double& alfa, const MatrixX& dyadic)
	{
		return dyadic * alfa;
	}

	inline MatrixS operator+(const double& alfa, const MatrixS& dyadic)
	{
		return dyadic + alfa;
	}

	inline MatrixS operator-(const double& alfa, const MatrixS& dyadic)
	{
		return (dyadic - alfa) * (-1.);
	}

	inline MatrixS operator*(const double& alfa, const MatrixS& dyadic)
	{
		return dyadic * alfa;
	}


	// ================================================================================================
	//
	// Functions
	//
	// ================================================================================================

	/** @return the determinant of the dyadic.
	  */
	static inline double determinant(const Dyadic2N& T) {
		return T.determinant();
	}

	/** @return the determinant of the dyadic.
	  */
	static inline double determinant(const Dyadic2S& T) {
		return T.determinant();
	}

	/** @return the norm of the dyadic.
	  */
	static inline double norm(const Dyadic2N& T) {
		return T.norm();
	}

	/** @return the norm of the dyadic.
	  */
	static inline double norm(const Dyadic2S& T) {
		return T.norm();
	}

	/** @return the trace of the dyadic.
	  */
	static inline double trace(const Dyadic2N& T) {
		return T.trace();
	}

	/** @return the trace of the dyadic.
	  */
	static inline double trace(const Dyadic2S& T) {
		return T.trace();
	}

	/** @return a vector with the dyadic' eigenvalues.
	  */
	static inline std::vector<double> eigenvalues(const Dyadic2N& T) {
		return T.eigenvalues();
	}

	/** @return a vector with the dyadic' eigenvalues.
	  */
	static inline std::vector<double> eigenvalues(const Dyadic2S& T) {
		return T.eigenvalues();
	}

	/** @return a vector with the dyadic' inverse.
	  */
	static inline std::vector<double> invariants(const Dyadic2N& T) {
		return T.invariants();
	}

	/** @return a vector with the dyadic' inverse.
	  */
	static inline std::vector<double> invariants(const Dyadic2S& T) {
		return T.invariants();
	}

	/** @return the transpose of the dyadic.
	  */
	static inline Dyadic2N transpose(const Dyadic2N& T) {
		return T.transpose();
	}

	/** @return the inverse of the dyadic.
	  */
	static inline Dyadic2N inverse(const Dyadic2N& T) {
		return T.inverse();
	}

	/** @return the inverse of the dyadic.
	  */
	static inline Dyadic2S inverse(const Dyadic2S& T) {
		return T.inverse();
	}

	/** @return the symmetric part of Tensor.
	  */
	static inline Dyadic2S getSymmetric(const Dyadic2N& T) {
		return T.getSymmetric();
	}

	/** @return the antisymmetric part of Tensor.
	  */
	static inline Dyadic2N getAsymmetric(const Dyadic2N& T) {
		return T.getAsymmetric();
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2N& T1, const Dyadic2N& T2) {
		return T1.contraction(T2);
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2S& T1, const Dyadic2S& T2) {
		return T1.contraction(T2);
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2S& T1, const Dyadic2N& T2) {
		assert(T1.rows() == T2.rows());		// Size of dyadics does not correspond!

		double result = 0.;
		auto& V1 = T1.getVector();
		auto& V2 = T2.getVector();

		for (unsigned int i = 0; i < T1.rows(); ++i) {
			result += V1.at(i) * V2.at(i);
		}
		for (unsigned int i = T1.rows(); i < T1.size(); ++i) {
			result += 2. * V1.at(i) * V2.at(i);
		}
		return result;
	}

	/** @return the double-dot product of the dyadics (contraction).
	  * @param T1 First dyadic to be used on the double-dot product.
	  * @param T2 First dyadic to be used on the double-dot product.
	  */
	static inline double contraction(const Dyadic2N& T1, const Dyadic2S& T2) {
		return contraction(T2, T1);
	}
}  // End of M2S2 namespace
