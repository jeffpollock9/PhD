#ifndef SCORESMATRIX_H_
#define SCORESMATRIX_H_

#include "../Utils/Vector.h"
#include "Line.h"

const static int DEFAULT_DIM = 15;

class ScoresMatrix {
public:
	ScoresMatrix();
	ScoresMatrix(int dim);

	int getDim();

	const double& at(int homeScore, int awayScore) const;
	double& at(int homeScore, int awayScore);

	double getHomeWin();
	double getAwayWin();
	double getDraw();

	double getGoal(int goal);
	double getUnderGoal(int goal);

	double getSupremacy(int supremacy);
	double getHomeWinByMoreThan(int supremacy);

	double getHomeHandicap(double line);
	double getUnderLine(double line);
private:
	DoubleMatrix values;
	int dim;
};

inline ScoresMatrix::ScoresMatrix() {
	this->dim = DEFAULT_DIM;
	this->values.resize(this->dim);
	for (int i = 0; i < this->dim; ++i) {
		this->values[i].resize(this->dim);
	}
}

inline ScoresMatrix::ScoresMatrix(int dim) {
	this->dim = dim;
	this->values.resize(this->dim);
	for (int i = 0; i < this->dim; ++i) {
		this->values[i].resize(this->dim);
	}
}

inline int ScoresMatrix::getDim() {
	return this->dim;
}

inline const double& ScoresMatrix::at(int homeScore, int awayScore) const {
	return this->values[homeScore][awayScore];
}

inline double& ScoresMatrix::at(int homeScore, int awayScore) {
	return this->values[homeScore][awayScore];
}

inline double ScoresMatrix::getHomeWin() {

	int homeScore, awayScore;
	double prob = 0.0;

	for (homeScore = 1; homeScore < this->dim; ++homeScore) {
		for (awayScore = 0; awayScore < homeScore; ++awayScore) {
			prob += this->values[homeScore][awayScore];
		}
	}

	return prob;
}

inline double ScoresMatrix::getAwayWin() {

	int homeScore, awayScore;
	double prob = 0.0;

	for (awayScore = 1; awayScore < this->dim; ++awayScore) {
		for (homeScore = 0; homeScore < awayScore; ++homeScore) {
			prob += this->values[homeScore][awayScore];
		}
	}

	return prob;
}

inline double ScoresMatrix::getDraw() {

	int score;
	double prob = 0.0;

	for (score = 0; score < this->dim; ++score) {
		prob += this->values[score][score];
	}

	return prob;
}

inline double ScoresMatrix::getGoal(int goal) {

	double prob = 0.0;

	for (int i = 0; i <= goal; ++i) {
		prob += this->values[i][goal - i];
	}

	return prob;
}

inline double ScoresMatrix::getUnderGoal(int goal) {

	double prob = 0.0;

	for (int g = 0; g < goal; ++g) {
		prob += getGoal(g);
	}

	return prob;
}

inline double ScoresMatrix::getUnderLine(double line) {

	int integerLine = (int) line;
	double underProb = getUnderGoal(integerLine);

	if (Line::isWholeLine(line)) {
		return underProb / (1.0 - getGoal(integerLine));
	} else if (Line::isOneQuarterLine(line)) {
		double exactProb = getGoal(integerLine);
		return (underProb + exactProb / 2.0) / (1.0 - exactProb / 2.0);
	} else if (Line::isHalfLine(line)) {
		return underProb + getGoal(integerLine);
	} else {
		return (underProb + getGoal(integerLine))
				/ (1.0 - getGoal(integerLine + 1) / 2.0);
	}
}

inline double ScoresMatrix::getSupremacy(int supremacy) {

	double prob = 0.0;
	int homeScore, awayScore;

	if (supremacy > 0) {
		for (awayScore = 0; awayScore < (this->dim - supremacy); ++awayScore) {
			homeScore = awayScore + supremacy;
			prob += this->values[homeScore][awayScore];
		}
	} else {
		for (homeScore = 0; homeScore < (this->dim + supremacy); ++homeScore) {
			awayScore = homeScore - supremacy;
			prob += this->values[homeScore][awayScore];
		}
	}

	return prob;
}

inline double ScoresMatrix::getHomeWinByMoreThan(int supremacy) {

	double prob = 0.0;

	for (int sup = (supremacy + 1); sup < this->dim; ++sup) {
		prob += getSupremacy(sup);
	}

	return prob;
}

inline double ScoresMatrix::getHomeHandicap(double line) {

	int integerLine = (int) line;
	double homeProb = getHomeWinByMoreThan(-integerLine);

	if (Line::isWholeLine(line)) {
		return homeProb / (1.0 - getSupremacy(-integerLine));
	}

	if (line > 0.0) {
		if (Line::isOneQuarterLine(line)) {
			double exactProb = getSupremacy(-integerLine);
			return (homeProb + exactProb / 2.0) / (1.0 - exactProb / 2.0);
		} else if (Line::isHalfLine(line)) {
			return homeProb + getSupremacy(-integerLine);
		} else {
			return (homeProb + getSupremacy(integerLine))
					/ (1.0 - getSupremacy(integerLine - 1) / 2.0);
		}
	} else {
		if (Line::isOneQuarterLine(line)) {
			return homeProb / (1.0 - getSupremacy(-integerLine) / 2.0);
		} else if (Line::isHalfLine(line)) {
			return homeProb;
		} else {
			homeProb -= getSupremacy(integerLine + 1);
			double exactProb = getSupremacy(integerLine + 1);
			return (homeProb + exactProb / 2.0) / (1.0 - exactProb / 2.0);
		}
	}
}

#endif /* SCORESMATRIX_H_ */
