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
#include "Dyadic2N.h"

// ================================================================================================
//
// Dyadic4S class
//
// ================================================================================================
namespace M2S2 {
	/** @class Dyadic4S
	 * @brief Symmetric 4th order tensor.
	 * @details Symmetric 4th order (rank) tensors of 2 or 3 dimensional vector space (saved in Voigt notation using row major).
	 */
	class Dyadic4S {
	public:
		/** Symmetric 4th order tensor.
		  */
		Dyadic4S()
		{
			mv_nDim = 0;
			mv_nSize = 0;
			mv_nVoigt = 0;
			mv_Values.clear();
		}

		/** Symmetric 4th order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space
		  */
		Dyadic4S(unsigned int nDim) : mv_nDim(nDim)
		{
			mv_nVoigt = 3 * nDim - 3;
			mv_nSize = (unsigned int)(0.5 * mv_nVoigt * mv_nVoigt + 0.5 * mv_nVoigt);
			mv_Values.resize(mv_nSize);
		}

		/** Symmetric 4th order tensor, for 2 or 3 dimensional vector space.
		  * @param nDim Dimensionality of vector space.
		  * @param value Value to initiate the entire dyadic.
		  */
		Dyadic4S(unsigned int nDim, const double& value) : mv_nDim(nDim)
		{
			mv_nVoigt = 3 * nDim - 3;
			mv_nSize = (unsigned int)(0.5 * mv_nVoigt * mv_nVoigt + 0.5 * mv_nVoigt);
			mv_Values.resize(mv_nSize, value);
		}

		/** Symmetric 4th order tensor, for 2 or 3 dimensional vector space.
		  * @param value Vector with either 6 or 21 values.
		  */
		Dyadic4S(const std::vector<double>& value)
		{
			mv_nSize = (unsigned int)(value.size());
			assert(mv_nSize == 6 || mv_nSize == 21);	// Size of vector does not correspond to 2nd order tensor from 2 or 3 dimensional vector space!

			mv_nDim = (mv_nSize == 6) ? 2 : 3;
			mv_nVoigt = 3 * mv_nDim - 3;
			mv_Values = { value.begin(), value.begin() + mv_nSize };
		}

		/** Copy constructor for symmetric 4th order tensor, for 2 or 3 dimensional vector space.
		  * @param other Dyadic to be copied.
		  */
		Dyadic4S(const Dyadic4S& other) {
			mv_nDim = other.mv_nDim;
			mv_nSize = other.mv_nSize;
			mv_nVoigt = other.mv_nVoigt;
			mv_Values = other.mv_Values;
		}

		/** Move constructor for symmetric 4th order tensor, for 2 or 3 dimensional vector space.
		  * @param other Dyadic to be moved.
		  */
		Dyadic4S(Dyadic4S&& other) noexcept : mv_nDim(other.mv_nDim), mv_nSize(other.mv_nSize), mv_nVoigt(other.mv_nVoigt), mv_Values(std::move(other.mv_Values)) { }

		/** Destructor.
		  */
		~Dyadic4S() { }

		/** Generate an identity with the required dimensionality.
		  * @param nDim Dyadic dimensionality.
		  * @param value Diagonal value. Default is 1.
		  */
		static Dyadic4S identity(unsigned int nDim, const double& value = double(1))
		{
			Dyadic4S result(nDim);
			auto& R = result.getVector();

			for (unsigned int i = 0; i < result.mv_nVoigt; ++i) {
				R.at((unsigned int)(i * (result.mv_nVoigt - i * 0.5 + 0.5))) = value;
			}
			return result;
		}

		/** Overloads operator << to stream the dyadic. */
		friend std::ostream& operator<<(std::ostream& output, const Dyadic4S& tensor)
		{
			output << tensor.print();
			return output;
		}

		/** Overloads operator >> to stream the dyadic. */
		friend std::istream& operator>>(std::istream& input, Dyadic4S& tensor)
		{
			for (unsigned int i = 0; i < tensor.rows(); i++) {
				for (unsigned int j = i; j < tensor.cols(); j++) {
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
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				output << "\t" << std::fixed << std::setprecision(precision) << std::setw(width) << at(i, 0);

				for (unsigned int j = 1; j < mv_nVoigt; j++) {
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
		void swap(Dyadic4S& other) noexcept
		{
			std::swap(mv_nDim, other.mv_nDim);
			std::swap(mv_nSize, other.mv_nSize);
			std::swap(mv_nVoigt, other.mv_nVoigt);
			mv_Values.swap(other.mv_Values);
		}

		/** Set all values to zero. Size remains unchanged.
		  */
		void clear()
		{
			memset(&mv_Values[0], 0., mv_Values.size() * sizeof(double));
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& at(unsigned int i, unsigned int j) {
			unsigned int pos = (i < j) ? (unsigned int)(i * (mv_nVoigt - i * 0.5 - 0.5) + j) : (unsigned int)(j * (mv_nVoigt - j * 0.5 - 0.5) + i);
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& at(unsigned int i, unsigned int j) const {
			unsigned int pos = (i < j) ? (unsigned int)(i * (mv_nVoigt - i * 0.5 - 0.5) + j) : (unsigned int)(j * (mv_nVoigt - j * 0.5 - 0.5) + i);
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline double& operator()(unsigned int i, unsigned int j) {
			unsigned int pos = (i < j) ? (unsigned int)(i * (mv_nVoigt - i * 0.5 - 0.5) + j) : (unsigned int)(j * (mv_nVoigt - j * 0.5 - 0.5) + i);
			return mv_Values.at(pos);
		}

		/** Access specified element - returns Tensor_ij (in Voigt notation)
		  * @param i First component.
		  * @param j Second component.
		  */
		inline const double& operator()(unsigned int i, unsigned int j) const {
			unsigned int pos = (i < j) ? (unsigned int)(i * (mv_nVoigt - i * 0.5 - 0.5) + j) : (unsigned int)(j * (mv_nVoigt - j * 0.5 - 0.5) + i);
			return mv_Values.at(pos);
		}

		/** @return the row size.
		  */
		unsigned int rows() const noexcept { return mv_nVoigt; }

		/** @return the column size.
		  */
		unsigned int cols() const noexcept { return mv_nVoigt; }

		/** @return the number of itens.
		  */
		unsigned int size() const noexcept { return mv_nSize; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		std::vector<double>& getVector() { return mv_Values; }

		/** @return a reference to the storing vector, providing direct access to the values.
		  */
		const std::vector<double>& getVector() const { return mv_Values; }

		/** @return a iterator to the beginning of the storing vector, providing direct access to the values.
		  */
		std::vector<double>::const_iterator begin() const { return mv_Values.begin(); }

		/** Overloads operator = for copy assignment operations -> T = O
		  * @param other Dyadic to be copied.
		  */
		Dyadic4S& operator=(const Dyadic4S& other)
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
				mv_nVoigt = other.mv_nVoigt;
				mv_Values.resize(mv_nVoigt * mv_nVoigt);
				mv_Values = other.mv_Values;
			}
			return *this;
		}

		/** Overloads operator = for move assignment operations -> T = &O
		  * @param other Dyadic to be moved.
		  */
		Dyadic4S& operator=(Dyadic4S&& other) noexcept
		{
			if (this != &other) {
				mv_nDim = other.mv_nDim;
				mv_nSize = other.mv_nSize;
				mv_nVoigt = other.mv_nVoigt;
				mv_Values = std::move(other.mv_Values);

				other.mv_nDim = 0;
				other.mv_nSize = 0;
				other.mv_nVoigt = 0;
			}
			return *this;
		}

		/** Overloads operator += for cumulative addition -> T += O -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic4S& operator+=(const Dyadic4S& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) += other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator -= for cumulative substraction -> T -= O -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic4S& operator-=(const Dyadic4S& other)
		{
			check_order(other);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return *this;
		}

		/** Overloads operator *= for cumulative multiplication -> T *= O -> T_ij = T_ik * O_kj
		  * @param other Dyadic to be multiplied with.
		  */
		Dyadic4S& operator*=(const Dyadic4S& other)
		{
			check_order(other);
			Dyadic4S result(mv_nDim);
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				for (unsigned int j = i; j < mv_nVoigt; j++) {
					for (unsigned int k = 0; k < mv_nVoigt; k++) {
						result.at(i, j) += at(i, k) * other.at(k, j);
					}
				}
			}
			swap(result);
			return *this;
		}

		/** Overloads operator + for addition -> T = T + O
		  * @param other Dyadic to be added.
		  */
		Dyadic4S operator+(const Dyadic4S& other) const
		{
			check_order(other);
			Dyadic4S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) += other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator - for substraction -> T = T - O
		  * @param other Dyadic to be substracted.
		  */
		Dyadic4S operator-(const Dyadic4S& other) const
		{
			check_order(other);
			Dyadic4S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); ++i) {
				result.mv_Values.at(i) -= other.mv_Values.at(i);
			}
			return result;
		}

		/** Overloads operator * for multiplication -> T_ij = T_ik * O_kj
		  * @param other Dyadic to be multiplied with.
		  */
		Dyadic4S operator*(const Dyadic4S& other) const
		{
			check_order(other);
			Dyadic4S result(mv_nDim);
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				for (unsigned int j = i; j < mv_nVoigt; j++) {
					for (unsigned int k = 0; k < mv_nVoigt; k++) {
						result.at(i, j) += at(i, k) * other.at(k, j);
					}
				}
			}
			return result;
		}

		/** Overloads operator * for Dot product -> T1_ij = T_ijkl * T_kl
		  * @param T Symmetric dyadic to be multiplied with.
		  * @return a symmetric dyadic with the inner product.
		  */
		Dyadic2S operator*(const Dyadic2S& T)
		{
			assert(T.rows() == mv_nVoigt); // Order of tensors differ!

			Dyadic2S result(mv_nDim);
			auto& R = T.getVector();
			auto& L = result.getVector();
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				for (unsigned int j = 0; j < mv_nVoigt; j++) {
					L.at(i) += at(i, j) * R.at(j);
				}
			}
			return result;
		}

		/** Overloads operator * for Dot product -> T1_ij = T_ijkl * T_kl
		  * @param T Asymmetric dyadic to be multiplied with.
		  * @return a asymmetric dyadic with the inner product.
		  */
		Dyadic2N operator*(const Dyadic2N& T)
		{
			assert(T.rows() == mv_nVoigt); // Order of tensors differ!

			Dyadic2N result(mv_nDim);
			auto& R = T.getVector();
			auto& L = result.getVector();
			for (unsigned int i = 0; i < mv_nVoigt; i++) {
				for (unsigned int j = 0; j < mv_nVoigt; j++) {
					L.at(i) += at(i, j) * R.at(j);
				}
			}
			return result;
		}

		/** Overloads operator + to sum a scalar -> T_ij = T_ij + alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic4S operator+(const double& alfa) const
		{
			Dyadic4S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) += alfa;
			}
			return result;
		}

		/** Overloads operator - to substract a scalar -> T_ij = T_ij - alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic4S operator-(const double& alfa) const
		{
			Dyadic4S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) -= alfa;
			}
			return result;
		}

		/** Overloads operator * to multiply with a scalar -> T_ij = T_ij * alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic4S operator*(const double& alfa) const
		{
			Dyadic4S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) *= alfa;
			}
			return result;
		}

		/** Overloads operator / to divide with a non-zero scalar -> T_ij = T_ij / alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic4S operator/(const double& alfa) const
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument("\nDivision by zero!\n");
			}

			Dyadic4S result(*this);
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				result.mv_Values.at(i) /= alfa;
			}
			return result;
		}

		/** Overloads operator += to sum a scalar -> T_ij += alfa
		  * @param alfa Scalar to be added.
		  * @return a dyadic with the result.
		  */
		Dyadic4S& operator+=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) += alfa;
			}
			return *this;
		}

		/** Overloads operator -= to substract a scalar -> T_ij -= alfa
		  * @param alfa Scalar to be substracted.
		  * @return a dyadic with the result.
		  */
		Dyadic4S& operator-=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) -= alfa;
			}
			return *this;
		}

		/** Overloads operator *= to multiply with a scalar -> T_ij *= alfa
		  * @param alfa Scalar to be multiplied with.
		  * @return a dyadic with the result.
		  */
		Dyadic4S& operator*=(const double& alfa)
		{
			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) *= alfa;
			}
			return *this;
		}

		/** Overloads operator /= to divide with a non-zero scalar -> T_ij /= alfa
		  * @param alfa Scalar to be divided with.
		  * @return a dyadic with the result.
		  */
		Dyadic4S& operator/=(const double& alfa)
		{
			if ((int)(alfa * 100000) == 0) {
				throw std::invalid_argument("\nDivision by zero!\n");
			}

			for (unsigned int i = 0; i < mv_Values.size(); i++) {
				mv_Values.at(i) /= alfa;
			}
			return *this;
		}

	private:
		unsigned int mv_nDim;
		unsigned int mv_nSize;
		unsigned int mv_nVoigt;
		std::vector<double> mv_Values;

		inline void check_order(const Dyadic4S& other) const
		{
			assert(mv_nVoigt == other.mv_nVoigt);	// Order of tensors differ
		}
	};
}  // End of M2S2 namespace
