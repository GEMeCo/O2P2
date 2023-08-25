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
#include <utility>      // std::move
#include <algorithm>    // std::copy() and std::assign()
#include <iostream>		// std::cout
#include <iomanip>		// Required by ios manipulations (std::setprecision)
#include <sstream>		// std::ostringstream
#include <cassert>      // assert

namespace M2S2 {
    /** @class triplet
      * @brief A triplet is a record for each non-zero entry in the matrix (row number, column number, and value)
      * @details These triplet format is auxiliary only for sparseMatrix. Use with caution.
      */
    class triplet {
    public:
        /** M2S2 triplet implementation. A triplet corresponds to the non-zero values of a sparse matrix along with their row and column index values
          */
        triplet() {
            mv_row = 0;
            mv_col = 0;
            mv_value = 0.;
        }

        /** M2S2 triplet implementation. A triplet corresponds to the non-zero values of a sparse matrix along with their row and column index values
          * @param row Row number
          * @param col Column number
          * @param value Value of the entry
          */
        triplet(int row, int col, double value) {
            mv_row = row;
            mv_col = col;
            mv_value = value;
        }

        /** Copy constructor for M2S2 triplet.
          * @param other Triplet to be copied.
          */
        triplet(const triplet& other) {
            mv_row = other.mv_row;
            mv_col = other.mv_col;
            mv_value = other.mv_value;
        }

        /** Move constructor for M2S2 triplet.
          * @param other Triplet to be moved.
          */
        triplet(triplet&& other) noexcept : mv_row(other.mv_row), mv_col(other.mv_col), mv_value(other.mv_value) { }

        /** Destructor.
          */
        ~triplet() { }

    public:
        /** @brief row index */
        int mv_row;

        /** @brief column index */
        int mv_col;

        /** @brief Value of the matrix element. */
        double mv_value;
    };

    /** @class CSR
      * @brief CSR stands for Compressed Sparse Row, a Sparse BLAS format for matrices.
      * @details These CSR class was created only to output data. All operations must be performed in sparseMatrix class.
      */
    class CSR {
    public:
        /** M2S2 Compressed Sparse Row matrices.
          */
        CSR() {
            mv_sym = false;
            mv_size = 0;
        }

        /** Copy constructor for M2S2 CSR.
          * @param other CSR to be copied.
          */
        CSR(const CSR& other) {
            mv_sym = other.mv_sym;
            mv_size = other.mv_size;

            mv_nz.assign(other.mv_nz.begin(), other.mv_nz.end());
            mv_col.assign(other.mv_col.begin(), other.mv_col.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());
        }

        /** Move constructor for M2S2 CSR.
          * @param other CSR to be moved.
          */
        CSR(CSR&& other) noexcept 
            : mv_sym(other.mv_sym), mv_size(other.mv_size), mv_nz(other.mv_nz), mv_col(other.mv_col), mv_value(other.mv_value) { }

        /** Destructor.
          */
        ~CSR() { }

        /** Removes all elements from the matrix, leaving the container empty
          */
        void destroy()
        {
            mv_sym = false;
            mv_size = 0;

            std::vector<int>().swap(mv_nz);
            std::vector<int>().swap(mv_col);
            std::vector<double>().swap(mv_value);
        }

    public:
        /** @brief If not symmetric, sym = false */
        bool mv_sym = false;

        /** @brief Number of rows / columns */
        int mv_size;

        /** @brief nz[i] - start of line i; nz[i+1] - end of line i. */
        std::vector<int> mv_nz;

        /** @brief Column index. */
        std::vector<int> mv_col;

        /** @brief Value of the matrix elements. */
        std::vector<double> mv_value;
    };


    /** @class CSC
      * @brief CSC stands for Compressed Sparse Column, a Sparse BLAS format for matrices.
      * @details These CSC class was created only to output data. All operations must be performed in sparseMatrix class.
      */
    class CSC {
    public:
        /** M2S2 Compressed Sparse Row matrices.
          */
        CSC() {
            mv_sym = false;
            mv_size = 0;
        }

        /** Copy constructor for M2S2 CSC.
          * @param other CSC to be copied.
          */
        CSC(const CSC& other) {
            mv_sym = other.mv_sym;
            mv_size = other.mv_size;

            mv_nz.assign(other.mv_nz.begin(), other.mv_nz.end());
            mv_row.assign(other.mv_row.begin(), other.mv_row.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());
        }

        /** Move constructor for M2S2 CSC.
          * @param other CSC to be moved.
          */
        CSC(CSC&& other) noexcept
            : mv_sym(other.mv_sym), mv_size(other.mv_size), mv_nz(other.mv_nz), mv_row(other.mv_row), mv_value(other.mv_value) { }

        /** Destructor.
          */
        ~CSC() { }

        /** Removes all elements from the matrix, leaving the container empty
          */
        void destroy()
        {
            mv_sym = false;
            mv_size = 0;

            std::vector<int>().swap(mv_nz);
            std::vector<int>().swap(mv_row);
            std::vector<double>().swap(mv_value);

            mv_nz.clear();
            mv_row.clear();
            mv_value.clear();
        }

    public:
        /** @brief If not symmetric, sym = false */
        bool mv_sym = false;

        /** @brief Number of rows / columns */
        int mv_size;

        /** @brief nz[i] - start of line i; nz[i+1] - end of line i. */
        std::vector<int> mv_nz;

        /** @brief Row index. */
        std::vector<int> mv_row;

        /** @brief Value of the matrix elements. */
        std::vector<double> mv_value;
    };


    /** @class line
      * @brief Auxiliary class for sparseMatrix.
      * @details Holds information about each line of the sparseMatrix.
      *
      * Comment (1)
      * resize variable is used to control the resize, according to the
      * adopted policy
      * First attempt to create a policy:
      * 1st resize -> new size is 2x initial size;
      * 2nd resize -> new size is 3x initial size;
      * 3rd resize -> new size is 4x initial size;
      * 4th resize -> new size is 2x last size;
      * 5th resize -> new size is 2x last size;
      * 6th resize -> stop program and show message!
      *
      * Comment (2)
      * Each time new values are pushed in the line, assembled is set to 0, i.e. false.
      * Each time the line is sorted and equal terms are summed, assembled is set to 1.
      * In other words, this variables is used to verify if it is necessary to assemble
      * the line. It is not necessary to assemble it to push new terms, but to use the
      * sparse matrix, all lines must be assembled.
      */
    class line {
    public:
        /** M2S2 row / column constructor.
          */
        line() {}

        /** M2S2 row / column constructor.
          * @param index vector of indexes of the values in the line.
          * @param value vector of values in the line.
          */
        line(const std::vector<int>& index, const std::vector<double>& value) {
            assert(index.size() == value.size()); // Size of input vectors are not of the same size

            mv_index.assign(index.begin(), index.end());
            mv_value.assign(value.begin(), value.end());
        }

        /** M2S2 row / column constructor.
          * @param size Size to be allocated.
          */
        line(const int& size) {
            reserve(size);
        }

        /** Copy constructor for M2S2 line.
          * @param other Line to be copied.
          */
        line(const line& other) {
            mv_index.reserve(other.mv_index.capacity());
            mv_value.reserve(other.mv_value.capacity());

            mv_index.assign(other.mv_index.begin(), other.mv_index.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());

            mv_assembled = other.mv_assembled;
        }

        /** Move constructor for M2S2 line.
          * @param other Line to be moved.
          */
        line(line&& other) noexcept
        {
            mv_index = std::move(other.mv_index);
            mv_value = std::move(other.mv_value);

            mv_assembled = other.mv_assembled;
        }

        /** Destructor.
          */
        ~line() { }

        /** Move assignment operator for M2S2 line.
          * @param other Line to be moved.
          */
        line& operator=(line&& other) noexcept
        {
            mv_index = std::move(other.mv_index);
            mv_value = std::move(other.mv_value);

            mv_assembled = other.mv_assembled;
            return *this;
        }

        /** Copy assignment operator for M2S2 line.
          * @param other Line to be copied.
          */
        line& operator=(const line& other)
        {
            mv_index.reserve(other.mv_index.capacity());
            mv_value.reserve(other.mv_value.capacity());

            mv_index.assign(other.mv_index.begin(), other.mv_index.end());
            mv_value.assign(other.mv_value.begin(), other.mv_value.end());

            mv_assembled = other.mv_assembled;

            return *this;
        }

    public:
        /** @brief index of the values in the line */
        std::vector<int> mv_index;

        /** @brief value of the line item */
        std::vector<double> mv_value;

        /** @brief number of times line was resized */
        int mv_resizeCount = 0;

        /** @brief false -> not assembled; true -> assembled */
        bool mv_assembled = false;

    public:
        /** Reserve capacity to the line. Notice that size is not checked (if it is smaller than current).
          * @param size New size to be allocated.
          */
        void reserve(const int& size)
        {
            mv_index.reserve(size);
            mv_value.reserve(size);

            mv_resizeCount = 0;
            mv_assembled = false;
        }

        /** Removes all elements from the line, leaving the line empty
          * Notice that std::vector::clear does not change capacity, but here it does.
          */
        void destroy()
        {
            if (mv_index.size()) {
                std::vector<int>().swap(mv_index);
                std::vector<double>().swap(mv_value);

                mv_resizeCount = 0;
                mv_assembled = false;
            }
        }

        /** Deletes all elements in the line, but does not change its capacity. Like std::vector::clear()
          */
        void clear()
        {
            mv_index.clear();
            mv_value.clear();
            mv_resizeCount = 0;
            mv_assembled = false;
        }

        /** Overloads operator << to stream the line. */
        friend std::ostream& operator<<(std::ostream& output, const line& line)
        {
            output << line.print();
            return output;
        }

        /** Overloads operator >> to stream the line. */
        friend std::istream& operator>>(std::istream& input, line& line)
        {
            for (int i = 0; i < line.mv_index.size(); ++i) {
                input >> line.mv_index.at(i);
            }
            for (int i = 0; i < line.mv_index.size(); ++i) {
                input >> line.mv_value.at(i);
            }
            return input;
        }

        /** Prepare a string to print (to file or screen) line
          * @param precision Number of decimal digits after the decimal point (default is 4)
          * @param width Minimum number of characters to be written (default is 8)
          */
        const std::string print(const int precision = 4, const int width = 8) const
        {
            std::ostringstream output;
            output << "Line size: " << mv_index.size() << "\t" << std::fixed << std::setprecision(precision) << std::setw(width);
            for (int i = 0; i < mv_index.size(); ++i) {
                output << mv_index.at(i) << ":" << mv_value.at(i) << " ";
            }
            std::cout << std::endl;
            return output.str();
        }

        /** If there are terms with repeated index, these are summed
          */
        void addEqualTerms() {
            if (mv_index.size() && !mv_assembled) {
                sort();
                /* Adding equal terms will only work if line is alredy sorted */
                int pos = 0;
                for (int i = 1; i < mv_index.size(); i++) {
                    if (mv_index.at(i) == mv_index.at(pos))
                        mv_value.at(pos) += mv_value.at(i);
                    else {
                        pos++;
                        mv_index.at(pos) = mv_index.at(i);
                        mv_value.at(pos) = mv_value.at(i);
                    }
                }
                mv_index.resize(++pos); /* new size, but keeping memory (reserve) */
                mv_value.resize(pos); /* new size, but keeping memory (reserve) */
                mv_assembled = true;
            }
        }

        /** @return the position in line of index. If index is not found return -1.
          * @param index Index to be found.
          */
        int search(int& index)
        {
            // Trying to enhance performance by using a bissection search algorithm
            if (mv_index.size()) {
                if (!mv_assembled) addEqualTerms();
                int ib, im, ie;
                ib = 0;
                ie = mv_index.size() - 1;
                if ((mv_index.at(ib) > index) || (mv_index.at(ie) < index))
                    return -1;
                if (mv_index.at(ib) == index)
                    return ib;
                if (mv_index.at(ie) == index)
                    return ie;

                /* int iter = 0;*/
                while ((ie - ib) > 1) {
                    if (mv_index.at(ib) == index)
                        return ib;
                    if (mv_index.at(ie) == index)
                        return ie;
                    im = ib + (ie - ib) / 2;
                    if (mv_index.at(im) == index)
                        return im;
                    if (index <= mv_index.at(im)) {
                        ie = im;
                    }
                    else {
                        ib = im;
                    }
                }
                return -1;
            }
            else {
                return -1;
            }
        }

        /** @return the position in line of index. If index is not found return -1.
          * @param index Index to be found.
          */
        int search(int index) const
        {
            assert(mv_assembled); // Const version can only be used with already sorted line.

            // Trying to enhance performance by using a bissection search algorithm
            if (mv_index.size()) {
                int ib, im, ie;
                ib = 0;
                ie = mv_index.size() - 1;
                if ((mv_index.at(ib) > index) || (mv_index.at(ie) < index))
                    return -1;
                if (mv_index.at(ib) == index)
                    return ib;
                if (mv_index.at(ie) == index)
                    return ie;

                /* int iter = 0;*/
                while ((ie - ib) > 1) {
                    if (mv_index.at(ib) == index)
                        return ib;
                    if (mv_index.at(ie) == index)
                        return ie;
                    im = ib + (ie - ib) / 2;
                    if (mv_index.at(im) == index)
                        return im;
                    if (index <= mv_index.at(im)) {
                        ie = im;
                    }
                    else {
                        ib = im;
                    }
                }
                return -1;
            }
            else {
                return -1;
            }
        }

        /* Increase initial capacity up to five times */
        void resize() {
            /* To be used for automatic resize */
            int initSize, lastSize, newSize;
            switch (mv_resizeCount) {
            case 0:
                /* First resize */
                initSize = mv_index.capacity();
                newSize = 2 * initSize;
                break;
            case 1:
                /* Second resize */
                initSize = mv_index.capacity() / 2;
                newSize = 3 * initSize;
                break;
            case 2:
                /* Third resize */
                initSize = mv_index.capacity() / 3;
                newSize = 4 * initSize;
                break;
            case 3:
                /* Fourth resize */
                lastSize = mv_index.capacity();
                newSize = 2 * lastSize;
                break;
            case 4:
                /* Fifth resize */
                lastSize = mv_index.capacity();
                newSize = 2 * lastSize;
                break;
            default:
                throw std::invalid_argument("\n\nLine resize failed!\nVerify initial size!");
            }
            mv_resizeCount++;
            reserve(newSize);
        }

    private:

        /* Intended just to test search */
        const int naiveSearch(const int& index)
        {
            if (mv_index.size()) {
                if (!mv_assembled) addEqualTerms();
                if ((mv_index.at(0) > index) || (mv_index.at(mv_index.size() - 1) < index))
                    return -1;
                for (int i = 0; i < mv_index.size(); i++) {
                    if (index == mv_index.at(i))
                        return i;
                    else if (index < mv_index.at(i))
                        return -1;
                }
            }
            else {
                return -1;
            }
            return -1;
        }

        /* Sort line by indexes */
        void quicksort(int* index, double* value, const int& begin, const int& end)
        {
            int i, j, pivo, aux;
            double daux;
            i = begin;
            j = end - 1;
            pivo = index[(begin + end) / 2];
            while (i <= j)
            {
                while (index[i] < pivo && i < end)
                {
                    i++;
                }
                while (index[j] > pivo && j > begin)
                {
                    j--;
                }
                if (i <= j)
                {
                    aux = index[i];
                    index[i] = index[j];
                    index[j] = aux;

                    daux = value[i];
                    value[i] = value[j];
                    value[j] = daux;

                    i++;
                    j--;
                }
            }
            if (j > begin)
                quicksort(index, value, begin, j + 1);
            if (i < end)
                quicksort(index, value, i, end);
        }

        /* Sort line by indexes */
        void quicksort(std::vector<int>& index, std::vector<double>& value, const int& begin, const int& end)
        {
            int i, j, pivo, aux;
            double daux;
            i = begin;
            j = end - 1;
            pivo = index[(begin + end) / 2];
            while (i <= j)
            {
                while (index[i] < pivo && i < end)
                {
                    i++;
                }
                while (index[j] > pivo && j > begin)
                {
                    j--;
                }
                if (i <= j)
                {
                    aux = index[i];
                    index[i] = index[j];
                    index[j] = aux;

                    daux = value[i];
                    value[i] = value[j];
                    value[j] = daux;

                    i++;
                    j--;
                }
            }
            if (j > begin)
                quicksort(index, value, begin, j + 1);
            if (i < end)
                quicksort(index, value, i, end);
        }

        /* Sort line by indexes */
        void sort() {
            quicksort(mv_index, mv_value, 0, mv_index.size());
        }
    };
}  // End of M2S2 namespace
