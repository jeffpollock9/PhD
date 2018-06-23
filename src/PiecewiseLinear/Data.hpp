#ifndef PIECEWISELINEARDATA_HPP_
#define PIECEWISELINEARDATA_HPP_

#include <armadillo>
#include "PiecewiseLinearModel.hpp"

using namespace std;

class PiecewiseLinearData {
public:
	PiecewiseLinearData() {
		n = int();
	}
	PiecewiseLinearData(const string& fileNameX, const string& fileNameY, const bool hasHeader);

	void readX(const string& fileNameX, const bool hasHeader);
	void readY(const string& fileNameY, const bool hasHeader);

	const arma::mat& getX() const {
		return X;
	}

	const arma::colvec& getY() const {
		return Y;
	}

	const arma::rowvec atX(const int row) const {
		return X.row(row);
	}

	const double atY(const int i) const {
		return Y.at(i);
	}

	const int size() const {
		return n;
	}

	const Knot randomKnot() const;

private:
	int n;
	arma::colvec Y;
	arma::mat X;
};

inline PiecewiseLinearData::PiecewiseLinearData(const string& fileNameX, const string& fileNameY,
		const bool hasHeader) {
	readX(fileNameX, hasHeader);
	readY(fileNameY, hasHeader);
	if (X.n_rows != Y.n_rows) {
		throw runtime_error("rows in data X != rows in data Y");
	}
	n = X.n_rows;
}

inline void PiecewiseLinearData::readX(const string& fileName, const bool hasHeader) {

	string line, field;

	//determine rows
	ifstream file1(fileName.c_str());
	if (hasHeader) {
		getline(file1, line);
	}
	int nrow = 0;
	while (getline(file1, line)) {
		++nrow;
	}

	//determine cols
	ifstream file2(fileName.c_str());
	int ncol = 0;
	getline(file2, line);
	stringstream ss(line);
	while (getline(ss, field, ',')) {
		++ncol;
	}

	X = arma::mat(nrow, ncol);

	ifstream file3(fileName.c_str());
	if (hasHeader) {
		getline(file3, line);
	}
	int i = 0;
	while (getline(file3, line)) {
		stringstream ss(line);
		int j = 0;
		while (getline(ss, field, ',')) {
			X.at(i, j) = atof(field.c_str());
			++j;
		}
		++i;
	}
}

inline void PiecewiseLinearData::readY(const string& fileName, const bool hasHeader) {
	ifstream file(fileName.c_str());

	string line, field;
	vector<string> lineVector;

	vector<double> data;

	if (hasHeader) {
		getline(file, line);
	}

	while (getline(file, line)) {
		data.push_back(atof(line.c_str()));
	}

	int n = data.size();

	Y.resize(n);

	for (int i = 0; i < n; ++i) {
		Y[i] = data[i];
	}
}

inline const Knot PiecewiseLinearData::randomKnot() const {
	double u = arma::randu();
	double sum = 0.0;
	for (int i = 0; i < n; ++i) {
		sum += 1.0 / n;
		if (u < sum) {
			//TODO: temp while X is a column
			double tmp = X.row(i)[0];
			return Knot(tmp, Y.at(i));
		}
	}
	throw runtime_error("error finding random Knot");
}

#endif /* PIECEWISELINEARDATA_HPP_ */
