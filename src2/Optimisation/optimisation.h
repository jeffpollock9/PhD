#ifndef OPTIMISATION_H_
#define OPTIMISATION_H_

#include "../Utils/Vector.h"
#include <nlopt.hpp>

int iteration = 0;

struct Data {
	bool constrained;
	int team;
	int upToWeek;
	bool printTrace;
};

void trace(int &iteration, double logLikelihood, DoubleVector theta) {
	if (iteration % 500 == 0) {
		std::cout << "iteration: " << iteration << "\t logLikelihood: "
				<< logLikelihood << std::endl;
		int n = theta.size();
		for (int i = 0; i < n; ++i) {
			std::cout << theta[i] << ", ";
		}
		std::cout << "\n" << std::endl;
	}
	iteration++;

}

void trace(int &iteration, double logLikelihood, double penalty,
		DoubleVector theta) {
	if (iteration % 25 == 0) {
		std::cout << "iteration: " << iteration << ", logLikelihood: "
				<< logLikelihood << ", penalty: " << penalty
				<< ", penalisedLogLikelihood: " << logLikelihood - penalty
				<< std::endl;
		int n = theta.size();
		for (int i = 0; i < n; ++i) {
			std::cout << theta[i] << ", ";
		}
		std::cout << "\n" << std::endl;
	}
	iteration++;

}
#endif /* OPTIMISATION_H_ */
