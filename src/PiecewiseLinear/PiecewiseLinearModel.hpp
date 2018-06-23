#ifndef PIECEWISELINEARMODEL_HPP_
#define PIECEWISELINEARMODEL_HPP_

#include <vector>
#include "../Math/Stats.hpp"

class Knot {
public:
	Knot() {
		this->x = double();
		this->y = double();
	}

	Knot(const double x, const double y) {
		this->x = x;
		this->y = y;
	}

	friend ostream &operator<<(ostream &os, const Knot &knot) {
		return os << "{" << knot.x << ", " << knot.y << "}";
	}

	double x;
	double y;
};

inline bool operator<(const Knot& lhs, const Knot& rhs) {
	return lhs.x < rhs.x;
}

class PiecewiseLinearModel {
public:
	PiecewiseLinearModel() {
		this->k = int();
	}

	PiecewiseLinearModel(const std::vector<Knot>& knots) {
		setKnots(knots);
	}

	int getK() const {
		return k;
	}

	const std::vector<Knot>& getKnots() const {
		return knots;
	}

	void setKnots(const std::vector<Knot>& knots) {
		this->k = knots.size();
		this->knots = knots;
		sort();
	}

	void sort() {
		std::sort(knots.begin(), knots.end());
	}

	double at(const double x) const;

	void addKnot(const Knot& knot);

	void removeKnot(const int idx, Knot& removedKnot);
	void removeRandomKnot(Knot& removedKnot);

	void moveRandomKnot(Knot& oldKnot, Knot& newKnot, const double sigma);
	void moveKnot(const int idx, Knot& oldKnot, Knot& newKnot, const double sigma);

private:
	int k;
	std::vector<Knot> knots;
};

inline double PiecewiseLinearModel::at(const double x) const {
	if (k == 1) {
		return knots[0].y;
	}

	for (int i = 1; i < k; ++i) {
		if (x <= knots[i].x || i == k - 1) {
			double gradient = (knots[i].y - knots[i - 1].y) / (knots[i].x - knots[i - 1].x);
			double intercept = knots[i].y - gradient * knots[i].x;
			return intercept + x * gradient;
		}
	}
	throw runtime_error("error calculating PiecewiseLinearModel::at");
}

inline void PiecewiseLinearModel::addKnot(const Knot& knot) {
	++k;
	knots.push_back(knot);
	sort();
}

inline void PiecewiseLinearModel::removeKnot(const int idx, Knot& removedKnot) {
	if (idx >= k) {
		throw runtime_error("tried to remove knot out of bounds");
	} else {
		removedKnot = knots[idx];
		knots.erase(knots.begin() + idx);
		--k;
	}
}

inline void PiecewiseLinearModel::removeRandomKnot(Knot& removedKnot) {
	double u = arma::randu();
	double sum = 0.0;
	for (int i = 0; i < k; ++i) {
		sum += 1.0 / k;
		if (u < sum) {
			removeKnot(i, removedKnot);
			return;
		}
	}
}

inline void PiecewiseLinearModel::moveKnot(const int idx, Knot& oldKnot, Knot& newKnot, const double sigma) {
	oldKnot = knots[idx];
	knots[idx].y += rnorm(0.0, sigma);
	newKnot = knots[idx];
}

inline void PiecewiseLinearModel::moveRandomKnot(Knot& oldKnot, Knot& newKnot, const double sigma) {
	double u = arma::randu();
	double sum = 0.0;
	for (int i = 0; i < k; ++i) {
		sum += 1.0 / k;
		if (u < sum) {
			moveKnot(i, oldKnot, newKnot, sigma);
			return;
		}
	}
}

#endif /* PIECEWISELINEARMODEL_HPP_ */
