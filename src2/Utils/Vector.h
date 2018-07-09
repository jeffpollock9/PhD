#ifndef VECTOR_H_
#define VECTOR_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <functional>
#include <armadillo>

#ifdef myvector
class Vector {
public:
	Vector(int length);
	Vector(int length, double x);

	double& at(int index);
	const double& at(int index) const;

	void fill(double x);
	int size();

	double mean();
	double sum();
	double var();
	double min();
	double max();
private:
	std::vector<double> data;
	int length;
};

inline Vector::Vector(int length) {
	this->length = length;
	data.resize(length, 0.0);
}

inline Vector::Vector(int length, double x) {
	this->length = length;
	data.resize(length, x);
}

inline double& Vector::at(int index) {
	return data[index];
}

inline const double& Vector::at(int index) const {
	return data[index];
}

inline void Vector::fill(double x) {
	std::fill(data.begin(), data.end(), x);
}

inline int Vector::size() {
	return data.size();
}

inline double Vector::mean() {
	return sum() / data.size();
}

inline double Vector::sum() {
	return std::accumulate(data.begin(), data.end(), 0.0);
}

inline double Vector::var() {

	double m = mean();
	int n = data.size();

	std::vector<double> diff(n);
	std::transform(data.begin(), data.end(), diff.begin(),
			std::bind2nd(std::minus<double>(), m));

	return std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0)
			/ (n - 1);
}

inline double Vector::min() {
	return *std::min_element(data.begin(), data.end());
}

inline double Vector::max() {
	return *std::max_element(data.begin(), data.end());
}

class Matrix {
public:
	Matrix(int nrow, int ncol);
	Matrix(int nrow, int ncol, double x);

	double& at(int row, int col);
	const double& at(int row, int col) const;

	const int nrow() const;
	const int ncol() const;

	const Vector row(int row) const;
	const Vector col(int col) const;

	void setRow(int row, Vector data);
	void setCol(int col, Vector data);

//	Vector& col(int col);
//	const Vector& col(int col) const;

//	Vector rowMeans();
//	Vector rowSums();
//	Vector rowVars();
//	Vector rowMins();
//	Vector rowMaxs();
private:
	const int index(int row, int col) const;
	std::vector<double> _data;
	int _nrow, _ncol;
};

inline Matrix::Matrix(int nrow, int ncol) {
	this->_nrow = nrow;
	this->_ncol = ncol;
	_data.resize(nrow * ncol, 0.0);
}

inline Matrix::Matrix(int nrow, int ncol, double x) {
	this->_nrow = nrow;
	this->_ncol = ncol;
	_data.resize(nrow * ncol, x);
}

inline const int Matrix::nrow() const {
	return this->_nrow;
}

inline const int Matrix::ncol() const {
	return this->_ncol;
}

inline const int Matrix::index(int row, int col) {
	return row * this->_ncol + col;
}

inline double& Matrix::at(int row, int col) {
	return _data.at(index(row, col));
}

inline const double& Matrix::at(int row, int col) const {
	return _data.at(index(row, col));
}

inline void Matrix::setRow(int row, Vector data) {
	for (int col = 0; col < this->_ncol; ++col) {
		_data.at(index(row, col)) = data.at(col);
	}
}

inline void Matrix::setCol(int col, Vector data) {
	for (int row = 0; row < this->_nrow; ++row) {
		_data.at(index(row, col)) = data.at(row);
	}
}

inline const Vector Matrix::row(int row) const {
	int length = this->_ncol;
	Vector result(length);
	for (int col = 0; col < length; ++col) {
		result.at(col) = _data.at(index(row, col));
	}
	return result;
}

inline const Vector Matrix::col(int row) const {
	int length = this->_nrow;
	Vector result(length);
	for (int row = 0; row < length; ++row) {
		result.at(row) = _data.at(index(row, col));
	}
	return result;
}

//inline const Vector& Matrix::col(int col) const {
//	int nrows = nrow();
//	Vector vector(nrows);
//	for (int row = 0; row < nrows; ++row) {
//		vector.at(row) = matrix.at(row).at(col);
//	}
//	return vector;
//}

//inline Vector Matrix::rowMeans() {
//	Vector means(this->nrow);
//	for (int i = 0; i < this->nrow; ++i) {
//		means.at(i) = data.at(i).mean();
//	}
//	return means;
//}
#endif

typedef std::vector<double> DoubleVector;
typedef std::vector<int> IntegerVector;
typedef std::vector<std::vector<double> > DoubleMatrix;

DoubleMatrix createDoubleMatrix(const int &nrow, const int &ncol) {
	DoubleMatrix mat(nrow);
	for (int i = 0; i < nrow; ++i) {
		mat[i].resize(ncol);
	}
	return mat;
}

DoubleVector seq(const int &size, const double &from, const double &to) {
	DoubleVector vec(size);
	double by = (to - from) / (size - 1.0);
	double addition = 0.0;
	for (int i = 0; i < size; ++i) {
		vec[i] = from + addition;
		addition += by;
	}
	return vec;
}

arma::rowvec arma_seq(const int &size, const double &from, const double &to) {
	arma::rowvec vec(size);
	double by = (to - from) / (size - 1.0);
	double addition = 0.0;
	for (int i = 0; i < size; ++i) {
		vec[i] = from + addition;
		addition += by;
	}
	return vec;
}

#endif /* VECTOR_H_ */
