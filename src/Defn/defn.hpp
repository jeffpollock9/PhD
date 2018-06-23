#ifndef DEFN_HPP_
#define DEFN_HPP_

#include <armadillo>
#include <string>
#include <iomanip>
#include <chrono>

using namespace std;

typedef pair<int, int> ScoreLine;

typedef arma::mat Scores;

typedef vector<double> DoubleVector;
typedef vector<int> IntegerVector;

typedef string Date;
typedef string Season;
typedef string Team;

int iteration = 0;
void trace(int &iteration, const double logLikelihood, const DoubleVector& theta) {
	if (iteration % 50 == 0) {
		std::cout << "iteration: " << iteration << "\t logLikelihood: "
				<< logLikelihood << std::endl;
		int n = theta.size();
		for (int i = 0; i < n; ++i) {
			std::cout << theta[i] << ", ";
		}
		std::cout << "\n" << std::endl;
	}
	++iteration;
}

class DoublePair: public pair<double, double> {
public:
	DoublePair() : pair<double, double>() {};
	DoublePair(const double& x1, const double& x2) {
		this->first = x1;
		this->second = x2;
	}
	DoublePair& operator+=(const DoublePair& rhs) {
		this->first += rhs.first;
		this->second += rhs.second;
		return *this;
	}
	DoublePair& operator-=(const DoublePair& rhs) {
		this->first -= rhs.first;
		this->second -= rhs.second;
		return *this;
	}
	DoublePair& operator/=(const double& x) {
		this->first /= x;
		this->second /= x;
		return *this;
	}
	DoublePair& operator*=(const double& x) {
		this->first *= x;
		this->second *= x;
		return *this;
	}
	friend ostream &operator<<(ostream &os, const DoublePair &x) {
		return os << "{" << x.first << ", " << x.second << "}";
	}
};

inline DoublePair operator+(DoublePair lhs, const DoublePair& rhs) {
	lhs += rhs;
	return lhs;
}
inline DoublePair operator-(DoublePair lhs, const DoublePair& rhs) {
	lhs -= rhs;
	return lhs;
}
inline DoublePair operator/(DoublePair lhs, const double& x) {
	lhs /= x;
	return lhs;
}
inline DoublePair operator*(DoublePair lhs, const double& x) {
	lhs *= x;
	return lhs;
}

#define INIT_TIMER auto start = std::chrono::high_resolution_clock::now();
#define START_TIMER start = std::chrono::high_resolution_clock::now();
#define STOP_TIMER(name) std::cout << "RUNTIME of " << name << ": " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()-start).count() / 1000.0 << " seconds " << std::endl;

#endif /* DEFN_HPP_ */
