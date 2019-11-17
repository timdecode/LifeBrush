// Changes Copyright (c) 2019, Timothy Davison.
// Original version by emmarose88:
// https://github.com/emmarose88/HungarianAlgorithm/blob/master/matching.cpp

#pragma once


#include <cmath>

#include <algorithm>

#include <Eigen/Dense>

using namespace std;

struct Hungarian
{
protected:
	enum class MinOrMax
	{
		Min,
		Max
	};

	static void square(Eigen::MatrixXf& m)
	{
		auto nRows = m.rows();
		auto nCols = m.cols();


		// if the dimensions are equal (square matrix), we're done
		// else, have to figure out larger dimension and pad with matrix max
		if (nRows != nCols)
		{
			float max_elem = m.maxCoeff(); // find the max element
			auto dim = std::max(nRows, nCols); // find the dimension for the new, square matrix

			Eigen::MatrixXf newm(dim, dim);


			// fill the matrix with the elements from temp and pad with max element
			for (int i = 0; i < dim; i++)
			{
				for (int j = 0; j < dim; j++)
				{
					if (i >= nRows || j >= nCols)
						newm(i, j) = max_elem;
					else
						newm(i, j) = m(i, j);
				}
			}

			m = std::move(newm);
		}
	}

	/*
	* reduce
	* reduces matrix based on row and column minimums
	*/
	static void reduce(Eigen::MatrixXf& m)
	{
		// subtract row minimum from each row
		for (int i = 0; i<m.rows(); i++) {
			float minElement = m.row(i).minCoeff();
			Eigen::VectorXf rMinusMin(m.rows());
			rMinusMin.fill(-minElement);
			m.row(i) += rMinusMin;
		}
	}

	/*
	* hasMark
	* if there is a starred/primed zero in the given row/col, returns it's index
	* else, returns -1
	*/
	static int hasMark(Eigen::VectorXf& v)
	{
		for (int i = 0; i<v.size(); i++) {
			if (v(i)) {
				return i;
			}
		}
		return -1;
	}

	/*
	* swapStarsAndPrimes
	* Swap stars and primes based on step 5 of Hungarian algorithm
	* Z0 is uncovered primed zero we've found
	* Z1 is the stared zero in the column of Z0 (if any)
	* Z2 is the primed zero in the row of Z1 (will always be one)
	* ...continue series until we reach a primed zero with no starred zero in its column
	* Unstar each starred zero, star each primed zero, erase all primes and uncover every line in the matrix
	*/
	static void swapStarsAndPrimes(int i, int j, Eigen::MatrixXf& stars, Eigen::MatrixXf& primes)
	{
		int primeRow = i;
		int primeCol = j;

		bool done = false;
		while (!done) {
			// find row index of row that has a 0* in the same col as the current 0'
			Eigen::VectorXf col = stars.col(primeCol);
			int starInPrimeColRow = hasMark(col);

			if (starInPrimeColRow < 0) {
				// star the prime we're looking at
				primes(primeRow, primeCol) = 0;
				stars(primeRow, primeCol) = 1;
				done = true;
			}
			else {
				// find which col has a 0' in the same row as z1
				Eigen::VectorXf row = primes.row(starInPrimeColRow);
				int primeInStarRowCol = hasMark(row);

				// star first primed zero
				primes(primeRow, primeCol) = 0;
				stars(primeRow, primeCol) = 1;
				//primes(starInPrimeColRow, primeInStarRowCol) = 0;
				//stars(starInPrimeColRow, primeInStarRowCol) = 1;

				// unstar starred zero
				stars(starInPrimeColRow, primeCol) = 0;

				// set index of last prime, will check it's column for 0*s next
				primeRow = starInPrimeColRow;
				primeCol = primeInStarRowCol;
			}
		}
		// clear primes
		primes.fill(0);
	}

	/*
	* findMatching
	* implementation of the Hungarian matching algorithm
	* referenced from: http://csclab.murraystate.edu/bob.pilgrim/445/munkres.html
	*/
	static void findMatching(Eigen::MatrixXf& m, Eigen::MatrixXf& result, MinOrMax minOrMax)
	{
		Eigen::MatrixXf n = m; // make a copy of m for reducing
		int dim = n.rows(); // dimension of matrix, used for checking if we've reduced
							// the matrix enough yet

		Eigen::MatrixXf stars(m.rows(), m.cols()); // matrix for storing our "starred" 0s (0*)
		stars.fill(0);
		Eigen::MatrixXf primes(m.rows(), m.cols()); // matrix for storing our "primed" 0s (0')
		primes.fill(0);
		Eigen::VectorXf rowCover(m.rows()); // keep track of which rows are "covered"
		rowCover.fill(0);
		Eigen::VectorXf colCover(m.cols()); // keep track of which columns are "covered"
		colCover.fill(0);

		// to do maximization rather than minimization, we have to
		// transform the matrix by subtracting every value from the maximum
		if (minOrMax == MinOrMax::Max) {
			float max = n.maxCoeff();
			Eigen::MatrixXf maxMat(n.rows(), n.cols());
			maxMat.fill(max);
			n = maxMat - n;
		}

		// Step 1 
		// Reduce matrix
		reduce(n);

		// Step 2
		// Find a zero in the matrix. If there is no starred zero in 
		// its row or column, star Z. Repeat for each element in the matrix.
		for (int i = 0; i<n.rows(); i++) {
			for (int j = 0; j<n.cols(); j++) {
				if (n(i, j) == 0 && !rowCover(i) && !colCover(j)) {
					stars(i, j) = 1;
					rowCover(i) = 1;
					colCover(j) = 1;
				}
			}
		}
		// covers need to be cleared for following steps
		rowCover.fill(0);
		colCover.fill(0);

		while (true) {
			// Step 3
			// Cover all columns that have a starred zero
			// If the number of columns with starred zeroes equals the matrix
			// dimensions, we are done! Otherwise, move on to step 4.
		step3:
			for (int j = 0; j<n.cols(); j++) {
				Eigen::VectorXf col = stars.col(j);
				if (hasMark(col) >= 0) {
					colCover(j) = 1;
				}
			}
			if (colCover.sum() == dim) {
				result = stars;
				return;
			}

			// Step 4
			// Find a non-covered zero and prime it
		step4:
			for (int i = 0; i<n.rows(); i++) {
				for (int j = 0; j<n.cols(); j++) {
					if (n(i, j) == 0 && !rowCover(i) && !colCover(j)) {
						primes(i, j) = 1;
						// if no starred zero in the row...
						Eigen::VectorXf row = stars.row(i);
						if (hasMark(row) < 0) {
							// Step 5
							// swap stars and primes            
							swapStarsAndPrimes(i, j, stars, primes);

							// clear lines
							rowCover.fill(0);
							colCover.fill(0);

							goto step3;
						}
						else {
							// cover row
							rowCover(i) = 1;

							// uncover column of the starred zero in the same row
							int col = hasMark(row);
							colCover(col) = 0;
						}
					}
				}
			}

			// Step 6
			// Should now be no more uncovered zeroes
			// Get the minimum uncovered element
			float min = 1000000;
			for (int i = 0; i<n.rows(); i++) {
				for (int j = 0; j<n.cols(); j++) {
					if (!rowCover(i) && !colCover(j) && n(i, j) < min) {
						min = n(i, j);
					}
				}
			}

			// Subtract minimum from uncovered elements, add it to elements covered twice
			for (int i = 0; i<n.rows(); i++) {
				for (int j = 0; j<n.cols(); j++) {
					if (!rowCover(i) && !colCover(j)) {
						n(i, j) -= min;
					}
					else if (rowCover(i) && colCover(j)) {
						n(i, j) += min;
					}
				}
			}

			goto step4;
		}
	}

public:
	// m, the cost matrix, must be square.
	static auto solve(Eigen::MatrixXf m) -> std::pair<float, std::vector<uint16>>
	{
		const auto rows = m.rows();
		const auto cols = m.cols();

		checkfSlow(rows == cols, TEXT("The matrix must be square."));

		// create an empty matrix to put the result in
		Eigen::MatrixXf result(rows, cols);
		result.fill(0);

		// run the Hungarian (Munkres) algorithm to find the maximal matching
		findMatching(m, result, MinOrMax::Min);

		float cost = 0.0f;

		auto pairs = std::make_pair(1.0f, std::vector<uint16>(rows));

		for (int c = 0; c < cols; ++c)
		{
			for( int r = 0; r < rows; ++r )
			{
				if( result(r,c) >= 1 )
				{
					cost += m(r,c);
					pairs.second[r] = c;
				}
			}
		}

		return std::move(pairs);
	}
};

