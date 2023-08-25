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

// Standard libraries
#include <vector>
#include <iostream>		// required by std::cout
#include <iomanip>		// Required by ios manipulations
#include <sstream>		// required by std::ostringstream
#include <cmath>		// required by std::sqrt / std::acos
#include <cassert>		// required by assert

// M2S2 libraries
#include "Dyadic2S.h"

// ================================================================================================
//
// Dyadic2N class
//
// ================================================================================================
namespace M2S2 {
	/** @class Dyadic2N
	 * @brief NON-Symmetric 2nd order tensor.
	 * @details Asymmetric 2nd order (rank) tensors of 2 or 3 dimensional vector space (saved in Voigt notation using row major).
	 */
	class Dyadic2N {
	public:
		/** Asymmetric 2nd order tensor.
		  */
		Dyadic2N()
		{
			mv_nDim = 0;
			mv_nSize = 0;
			mv_Values.clear();
		}

		/** Asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space
		  */
		Dyadic2N(unsigned int nDim) : mv_nDim(nDim)
		{
			mv_nSize = nDim * nDim;
			mv_Values.resize(mv_nSize);
		}

		/** Asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space.
		  * @param value Value to initiate the entire dyadic.
		  */
		Dyadic2N(unsigned int nDim, const double& value) : mv_nDim(nDim)
		{
			mv_nSize = nDim * nDim;
			mv_Values.resize(mv_nSize, value);
		}

		/** Asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param value Vector with either 4 or 9 values.
		  */
		Dyadic2N(const std::vector<double>& value)
		{
			mv_nDim = (int)(std::sqrt(value.size()));
			mv_nSize = mv_nDim * mv_nDim;

			assert(mv_nDim == 2 || mv_nDim == 3);	// Size of input vector is not from 2D or 3D
			mv_Values = { value.begin(), value.begin() + mv_nSize };
		}

		/** Copy constructor for asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param other Asymmetric dyadic to be copied.
		  */
		Dyadic2N(const Dyadic2N& other) {
			mv_nDim = other.mv_nDim;
			mv_nSize = other.mv_nSize;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for asymmetric 2nd order tensor, for 2 or 3 dimensional vector space.
		  * @param other Dyadic to be moved.
		  */
		Dyadic2N(Dyadic2N&& other) noexcept : mv_nDim(other.mv_nDim), mv_nSize(other.mv_nSize), mv_Values(std::move(other.mv_Values)) { }

		/** Destructor.
		  */
		~Dyadic2N() { }

		/** Generate an identity with the required dimensionality.
		  * @param nDim Dyadic dimensionality.
		  * @param value Diagonal value. Default is 1.
		  */
		static Dyadic2N identity(unsigned int nDim, const double& value = double(1))
		{
			Dyadic2N result(nDim);
			for (unsigned int i = 0; i < result.mv_nDim; ++i) {
				result.at(i, i) = value;
			}
			return result;
		}

		/** Overloads operator << to stream the dyadic. */
		friend std::ostream& operator<<(std::ostream& output, const Dyadic2N& tensor)
		{
			output << tensor.print();
			return output;
		}

		/** Overloads operator >> to stream the dyadic. */
		friend std::istream& operator>>(std::istream& input, Dyadic2N& tensor)
		{
			for (unsigned int i = 0; i < tensor.rows(); i++) {
				for (unsigned int j = 0; j < tensor.cols(); j++) {
					input >> tensor.at(i, j);
				}
			}
			return input;
		}

		/** Prepare a string to print (to file or screen)
		  * @param precision Number of decimal digits after the decimal point (default is 4)
		  * @param width Minimum number of characters to be written (default is 8)
		  */
		const std::string print(const int precision = 4, const int width = 8) const
		{
			std::ostringstream output;
			output << std::endl;
			for (unsigned int i = 0; i < mv_nDim; i++) {
				output << "\t" << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, 0);

				for (unsigned int j = 1; j < mv_nDim; j++) {
					output << " " << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, j);
				}
				output << "\n";
			}
			output << std::endl;
			return output.str();
		}

		/** Swap dyadics
		  * @param other Dyadic to be swaped.
		  */
		void swap(Dyadic2N& other) noexcept
		{
			std::swap(mv_nDim, other.mv_nDim);
			std::swap(mv_nSize, other.mv_nSize);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(unsigned int i, unsigned int j) { return mv_Values.at(i * mv_nDim + j); }

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(unsigned int i, unsigned int j) const { return mv_Values.at(i * mv_nDim + j); }

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator()(unsigned int i, unsigned int j) { return mv_Values.at(i * mv_nDim + j); }

		/** Access specified element - returns Tensor_ij
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator()(unsigned int i, unsigned int j) const { return mv_Values.at(i * mv_nDim + j); }

		/** @return the row size.
		  */
		unsigned int rows() const noexcept { return mv_nDim; }

		/** @return the column size.
		  */
		unsigned int cols() const noexcept { return mv_nDim; }

		/** @return the number of itens.
		  */
		unsigned int size() const noexcept { return mv_nSize; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		std::vector<double>& getVector() { return mv_Values; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		const std::vector<double>& getVector() const { return mv_Values; }

		/** Returns a vector with tensor written in Voigt notation with mnemonics rule: 
		  * V(1) = T(1,1); V(2) = T(2,2); V(3) = T(3,3); V(4) = T(2,3); V(5) = T(1,3); V(6) = T(1,2); V(7) = T(3,2); V(8) = T(3,1); V(9) = T(2,1)
		  * @return a vector with tensor written in Voigt notation with mnemonics rule.
		  */
		std::vector<double> getVoigtMnenomics()
		{
			unsigned int mi_nVoigt = 3 * mv_nDim - 3;
			std::vector<double> result(mv_nSize);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.at(i) = at(i, i);
				for (unsigned int j = i + 1; j < mv_nDim; ++j) {
					result.at(mi_nVoigt - i - j) = at(i, j);
					result.at(mv_nSize - i - j) = at(j, i);
				}
			}
			return result;
		}

		/** @return a iterator to the beginning of the storing vector, providing direct access to the values.
		  */
		std::vector<double>::const_iterator begin() const { return mv_Values.begin(); }

		/** Overloads operator = for copy assignment operations -> T = O
		  * @param other Dyadic to be copied.
		  */
		Dyadic2N& operator=(const Dyadic2N& other)
		{
			if (this == &other) {
				return *this;
			}

			if (mv_nDim == other.mv_nDim)
			{
				mv_Values = other.mv_Values;
			}
			else
			{
				mv_nDim = other.mv_nDim;
				mv_nSize = other.mv_nSize;
				mv_Values.resize(mv_nSize);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations -> T = &O
		  * @param other Dyadic to be moved.
		  */
		Dyadic2N& operator=(Dyadic2N&& other) noexcept
		{
			if (this != &other) {
				mv_nDim = other.mv_nDim;
				mv_nSize = other.mv_nSize;
				mv_Values = std::move(other.mv_Values);

				other.mv_nDim = 0;
				other.mv_nSize = 0;
			}
			return *this;
		}

		/** Overloads operator += for cumulative addition -> T += O -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic2N& operator+=(const Dyadic2N& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) += other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator -= for cumulative substraction -> T -= O -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2N& operator-=(const Dyadic2N& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator *= for cumulative multiplication -> T *= O -> T_ij = T_ik * O_kj
		  * @param other Dyadic to be multiplied with.
		  */
		Dyadic2N& operator*=(const Dyadic2N& other)
		{
			check_order(other);
			Dyadic2N result(mv_nDim);
			auto& values = result.getVector();
			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					for (unsigned int k = 0; k < mv_nDim; k++) {
						values.at(i * mv_nDim + j) += mv_Values.at(i * mv_nDim + k) * other.mv_Values.at(k * mv_nDim + j);
					}
				}
			}
			swap(result);
			return *this;
		}

		/** Overloads operator + for addition -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic2N operator+(const Dyadic2N& other) const
		{
			check_order(other);
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) += other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic2N operator-(const Dyadic2N& other) const
		{
			check_order(other);
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
		  * @param other Dyadic to be multiplied with.
		  */
		Dyadic2N operator*(const Dyadic2N& other) const
		{
			check_order(other);
			Dyadic2N result(mv_nDim);
			auto& values = result.getVector();
			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					for (unsigned int k = 0; k < mv_nDim; k++) {
						values.at(i * mv_nDim + j) += mv_Values.at(i * mv_nDim + k) * other.mv_Values.at(k * mv_nDim + j);
					}
				}
			}
			return result;
		}

		/** Overloads operator * for Dot product -> V_i = T_ij * R_j
		  * @param R Vector to be multiplied with.
		  * @return a vector with the inner product.
		  */
		std::vector<double> operator*(const std::vector<double>& R)
		{
			assert(R.size() == mv_nDim);		// Order of tensor and vector differ
			std::vector<double> result(mv_nDim);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				for (unsigned int j = 0; j < mv_nDim; j++) {
					result.at(i) += at(i, j) * R.at(j);
				}
			}
			return result;
		}

		/** Overloads operator + to sum a scalar -> T_ij = T_ij + alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator+(const double& alfa) const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator-(const double& alfa) const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nDim; i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator*(const double& alfa) const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N operator/(const double& alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument("\nDivision by zero!\n");
			}

			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nSize; i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator+=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_nDim; i++) {
				mv_Values.at(i) += alfa;
			}
			return *this;
		}

		/** Overloads operator -= to substract a scalar -> T_ij -= alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator-=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_nDim; i++) {
				mv_Values.at(i) -= alfa;
			}
			return *this;
		}

		/** Overloads operator *= to multiply with a scalar -> T_ij *= alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator*=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) *= alfa;
			}
			return *this;
		}

		/** Overloads operator /= to divide with a non-zero scalar -> T_ij /= alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic2N& operator/=(const double& alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument("\nDivision by zero!\n");
			}

			for (unsigned int i = 0; i < mv_nSize; i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

		/** @return the determinant of the dyadic.
		  */
		double determinant() const
		{
			if (mv_nDim == 2) {
				return mv_Values.at(0) * mv_Values.at(3) - mv_Values.at(1) * mv_Values.at(2);
			}
			else {
				return mv_Values.at(0) * mv_Values.at(4) * mv_Values.at(8) +
					mv_Values.at(1) * mv_Values.at(5) * mv_Values.at(6) +
					mv_Values.at(2) * mv_Values.at(3) * mv_Values.at(7) -
					mv_Values.at(0) * mv_Values.at(5) * mv_Values.at(7) -
					mv_Values.at(1) * mv_Values.at(3) * mv_Values.at(8) -
					mv_Values.at(2) * mv_Values.at(4) * mv_Values.at(6);
			}
		}

		/** @return the norm of the dyadic -> square root of T:T
		  */
		double norm() const
		{
			return std::sqrt(contraction(*this));
		}

		/** @return the trace of the dyadic -> T_ii
		  */
		double trace() const
		{
			double result = 0;
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				result += at(i, i);
			}
			return result;
		}

		/** @return a vector with the dyadic' eigenvalues.
		  */
		std::vector<double> eigenvalues() const
		{
			std::cout << "\n\nWARNING: Eigenvalues of symmetric part only.\nOnly those will be always real." << std::endl << std::endl;
			auto TS = getSymmetric();
			std::vector<double> result = TS.eigenvalues();
			return result;
		}

		/** @return a vector with the dyadic' invariants.
		  */
		std::vector<double> invariants() const
		{
			std::vector<double> result(mv_nDim);
			if (mv_nDim == 2) {
				result.at(0) = trace();
				result.at(1) = determinant();
			}
			else {
				result.at(0) = trace();
				result.at(1) = mv_Values.at(0) * mv_Values.at(4) + mv_Values.at(4) * mv_Values.at(8) + mv_Values.at(8) * mv_Values.at(0)
					- mv_Values.at(1) * mv_Values.at(3) - mv_Values.at(5) * mv_Values.at(7) - mv_Values.at(2) * mv_Values.at(6);
				result.at(2) = determinant();
			}
			return result;
		}

		/** @return the transpose of the dyadic -> T^t
		  */
		Dyadic2N transpose() const
		{
			Dyadic2N result(*this);
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				for (unsigned int j = i + 1; j < mv_nDim; ++j) {
					result.at(i, j) = at(j, i);
					result.at(j, i) = at(i, j);
				}
			}
			return result;
		}

		/** @return the inverse of the dyadic -> T^-1
		  */
		Dyadic2N inverse() const
		{
			Dyadic2N result(mv_nDim);
			double det = determinant();
			if ((int)(det * 100000) == 0) {
				throw std::invalid_argument("\nTensor is singular!\n");
			}
			det = 1. / det;

			if (mv_nDim == 2) {
				result.at(0, 0) = det * at(1, 1);
				result.at(0, 1) = -1. * det * at(0, 1);
				result.at(1, 0) = -1. * det * at(1, 0);
				result.at(1, 1) = det * at(0, 0);
			}
			else
			{
				// For non-symmetric tensors
				result.at(0, 0) = at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1);
				result.at(0, 1) = at(0, 2) * at(2, 1) - at(0, 1) * at(2, 2);
				result.at(0, 2) = at(0, 1) * at(1, 2) - at(0, 2) * at(1, 1);
				result.at(1, 0) = at(1, 2) * at(2, 0) - at(1, 0) * at(2, 2);
				result.at(1, 1) = at(0, 0) * at(2, 2) - at(0, 2) * at(2, 0);
				result.at(1, 2) = at(0, 2) * at(1, 0) - at(0, 0) * at(1, 2);
				result.at(2, 0) = at(1, 0) * at(2, 1) - at(1, 1) * at(2, 0);
				result.at(2, 1) = at(0, 1) * at(2, 0) - at(0, 0) * at(2, 1);
				result.at(2, 2) = at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);

				result *= det;
			}
			return result;
		}

		/** @return the symmetric part of Tensor -> Sym(T)
		  */
		Dyadic2S getSymmetric() const
		{
			Dyadic2S result(mv_nDim);
			auto& values = result.getVector();

			if (mv_nDim == 2) {
				values.at(0) = mv_Values.at(0);
				values.at(1) = 0.5 * (mv_Values.at(1) + mv_Values.at(2));
				values.at(2) = mv_Values.at(3);
			}
			else
			{
				values.at(0) = mv_Values.at(0);
				values.at(1) = 0.5 * (mv_Values.at(1) + mv_Values.at(3));
				values.at(2) = 0.5 * (mv_Values.at(2) + mv_Values.at(6));
				values.at(3) = mv_Values.at(4);
				values.at(4) = 0.5 * (mv_Values.at(5) + mv_Values.at(7));
				values.at(5) = mv_Values.at(8);
			}
			return result;
		}

		/** @return the antisymmetric part of Tensor -> Asym(T)
		  */
		Dyadic2N getAsymmetric() const
		{
			Dyadic2S Sym = this->getSymmetric();
			Dyadic2N result(this->mv_nDim);
			for (unsigned int i = 0; i < result.rows(); ++i) {
				for (unsigned int j = 0; j < result.cols(); ++j) {
					result.at(i, j) = at(i, j) - Sym.at(i, j);
				}
			}
			return result;
		}

		/** @return the product of the transposed tensor and itself (always symmetric) -> T1^T * T1
		  */
		Dyadic2S getATA() const
		{
			Dyadic2S result(mv_nDim);
			for (unsigned int i = 0; i < mv_nDim; ++i) {
				for (unsigned int j = 0; j < mv_nDim; ++j) {
					for (unsigned int k = 0; k < mv_nDim; ++k) {
						result.at(i, j) += mv_Values.at(i * mv_nDim + k) * mv_Values.at(k * mv_nDim + j);
					}
				}
			}
			return result;
		}

		/** @return the double-dot product of the dyadics (contraction) -> T1:T2.
		  * @param other Dyadic to be used on the double-dot product.
		  */
		double contraction(const Dyadic2N& other) const
		{
			check_order(other);
			double result = 0.;
			for (unsigned int i = 0; i < mv_nSize; ++i) {
				result += mv_Values.at(i) * other.mv_Values.at(i);
			}
			return result;
		}

	private:
		unsigned int mv_nDim;
		unsigned int mv_nSize;
		std::vector<double> mv_Values;

		inline void check_order(const Dyadic2N& other) const
		{
			assert(mv_nDim == other.mv_nDim);	// Order of tensors differ
		}
	};
}  // End of M2S2 namespace
