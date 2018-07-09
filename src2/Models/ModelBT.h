#ifndef MODELBT_H_
#define MODELBT_H_

#include "../Utils/Stats.h"
#include <nlopt.hpp>

class ModelBT {
public:
	ModelBT(double delta, double h);
	DoubleVector getHomeDrawAway();
	double getHome();
	double getDraw();
	double getAway();

	double logLikelihood(const double homeScore, const double awayScore);

	void setAwayR(double awayR) {
		this->awayR = awayR;
	}
	void setH(double h) {
		this->h = h;
	}
	void setHomeR(double homeR) {
		this->homeR = homeR;
	}
	void setDelta(double delta) {
		this->delta[0] = -delta;
		this->delta[1] = delta;
	}
private:
	double delta[2], h, homeR, awayR;
};

inline ModelBT::ModelBT(double delta, double h) {
	setH(h);
	setDelta(delta);

	setHomeR(0.0);
	setAwayR(0.0);
}

inline DoubleVector ModelBT::getHomeDrawAway() {

	double linPred = h + homeR - awayR;

	double awayProb = plogis(delta[0] - linPred);
	double drawProb = plogis(delta[1] - linPred) - awayProb;
	double homeProb = 1.0 - (drawProb + awayProb);

	DoubleVector out(3);
	out[0] = awayProb;
	out[1] = drawProb;
	out[2] = homeProb;

	return out;
}

inline double ModelBT::logLikelihood(const double homeScore,
		const double awayScore) {

	if (homeScore > awayScore) {
		return log(getHome());
	} else if (homeScore == awayScore) {
		return log(getDraw());
	} else {
		return log(getAway());
	}
}

inline double ModelBT::getHome() {
	return 1.0 - plogis(delta[1] - (h + homeR - awayR));
}

inline double ModelBT::getDraw() {
	double linPred =  h + homeR - awayR;
	double awayProb = plogis(delta[0] - linPred);

	return plogis(delta[1] - linPred) - awayProb;
}

inline double ModelBT::getAway() {
	return plogis(delta[0] - (h + homeR - awayR));
}

#endif /* MODELBT_H_ */
