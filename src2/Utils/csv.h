#ifndef CSV_H_
#define CSV_H_

#include <armadillo>
#include "Vector.h"

const char* fileName(const int t) {
	std::stringstream stream;
	stream << "particles/week-";
	stream << t;
	stream << ".csv";

	return stream.str().c_str();
}

const char* fileNameSim(const int t) {
	std::stringstream stream;
	stream << "particles/sim-week-";
	stream << t;
	stream << ".csv";

	return stream.str().c_str();
}

const char* matchFile(const double s) {
	std::stringstream stream;
	stream << "particles/sigma-";
	stream << s;
	stream << ".csv";

	return stream.str().c_str();
}

void writeCSV(const DoubleMatrix &mat, const char* name) {

	std::ofstream file;

	file.open(name);

	int nRows = mat.size();
	int nCols = mat[0].size();

	for (int row = 0; row < nRows; ++row) {
		for (int col = 0; col < nCols; ++col) {
			file << mat[row][col];
			if (col != nCols - 1) {
				file << ", ";
			}
		}
		file << "\n";
	}

	file.close();
}

void writeCSV(const arma::mat &mat, const char* name) {

	std::ofstream file;

	file.open(name);

	int nRows = mat.n_rows;
	int nCols = mat.n_cols;

	for (int row = 0; row < nRows; ++row) {
		for (int col = 0; col < nCols; ++col) {
			file << mat.at(row, col);
			if (col != nCols - 1) {
				file << ", ";
			}
		}
		file << "\n";
	}

	file.close();
}

/*void writeCSV(const Eigen::MatrixXd &mat, const char* name) {

	std::ofstream file;

	file.open(name);

	int nRows = mat.rows();
	int nCols = mat.cols();

	for (int row = 0; row < nRows; ++row) {
		for (int col = 0; col < nCols; ++col) {
			file << mat(row, col);
			if (col != nCols - 1) {
				file << ", ";
			}
		}
		file << "\n";
	}

	file.close();
}*/
#endif /* CSV_H_ */
