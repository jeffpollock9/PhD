#ifndef THREADING_H_
#define THREADING_H_

#include <pthread.h>
#include <cmath>
#include <ctime>
#include "../Bayesian/simulateBT.h"

struct Gamma {
	double gamma_1, gamma_2, metric;
};

void *threadFunction(void *arg) {
	Gamma *gamma = static_cast<Gamma*>(arg);

	int gameIdx;
	double metric = 0.0;
	static const int nSamplesMCMC = 10000;
	static const int burn = 2000;
	static const int thin = 1;

	for (int week = 0; week < NUM_WEEKS; ++week) {

		DoubleVector theta(startBT, startBT + nParametersBT);
		DoubleMatrix probs = getProbsMatrix(theta, week, nSamplesMCMC, burn,
				thin, gamma->gamma_1, gamma->gamma_2);

		for (int g = 0; g < 10; ++g) {
			gameIdx = week * 10 + g;
			if (MATCH_HOME_SCORE[gameIdx] > MATCH_AWAY_SCORE[gameIdx]) {
				metric += log(probs[g][2]);
			} else if (MATCH_HOME_SCORE[gameIdx] == MATCH_AWAY_SCORE[gameIdx]) {
				metric += log(probs[g][1]);
			} else {
				metric += log(probs[g][0]);
			}
		}
	}

	std::cout << "thread running: " << gamma->gamma_1 << ", " << gamma->gamma_2 << " done." << std::endl;

	gamma->metric = metric;

	pthread_exit(NULL);

	return 0;
}

//#define getGamma
//#define thread
#ifdef thread
#include "Utils/threading.h"
#include "Bayesian/simulateBT.h"
int main(int argc, char *argvs[]) {

	clock_t startTime = clock();

	static const double rGamma1[] = {1.0};//{log(3.5), log(4), log(4.5)};
	static const double rGamma2[] = {0.5, 0.25};//{log(1.5), log(1.75), log(2)};

	static const int nRGamma1 = sizeof(rGamma1) / sizeof(double);
	static const int nRGamma2 = sizeof(rGamma2) / sizeof(double);

	pthread_t threads[nRGamma1][nRGamma2];
	Gamma data[nRGamma1][nRGamma2];

	for (int i = 0; i < nRGamma1; i++) {
		for (int j = 0; j < nRGamma2; ++j) {
			data[i][j].gamma_1 = rGamma1[i];
			data[i][j].gamma_2 = rGamma2[j];
			pthread_create(&threads[i][j], NULL, threadFunction,
					static_cast<void*>(&data[i][j]));
		}
	}

	for (int i = 0; i < nRGamma1; i++) {
		for (int j = 0; j < nRGamma2; ++j) {
			pthread_join(threads[i][j], NULL);
		}
	}

	DoubleMatrix metricsMat = createDoubleMatrix(nRGamma1, nRGamma2);
	for (int i = 0; i < nRGamma1; i++) {
		for (int j = 0; j < nRGamma2; ++j) {
			metricsMat[i][j] = data[i][j].metric;
		}
	}

	writeCSV(metricsMat, "metricsBT3.csv");

	std::cout << "time: " << (double) (clock() - startTime) / CLOCKS_PER_SEC
	<< " seconds" << std::endl;
}
//example of threading
#include <thread>
#include <iostream>
#include <cmath>
#include <chrono>

static const int n = 10000;
static double ans[n];

void threadFunction(int i) {
	double x = 0.0;
	for (int j = 0; j < 100000; ++j) {
		x += log(j + 3.14);
	}
	ans[i] = x;
}

int main(int argc, char **argv) {

	auto start = std::chrono::high_resolution_clock::now();

	std::thread t[n];

	for (int i = 0; i < n; ++i) {
		t[i] = std::thread(threadFunction, i);
	}

	for (int i = 0; i < n; ++i) {
		t[i].join();
	}
	auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Time for threads: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";

    start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < n; ++i) {
		threadFunction(i);
	}

    end = std::chrono::high_resolution_clock::now();

    std::cout << "Time for normal: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";
}
#endif

#endif /* THREADING_H_ */
