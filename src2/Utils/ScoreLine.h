
#ifndef SCORELINE_H_
#define SCORELINE_H_

#include <sstream>
#include <string>

class ScoreLine {
public:
	ScoreLine();
	ScoreLine(int homeScore, int awayScore);
	std::string toString();
    int getAwayScore();
    int getHomeScore();
    void setAwayScore(int awayScore);
    void setHomeScore(int homeScore);
private:
	int homeScore;
	int awayScore;
};

inline ScoreLine::ScoreLine() {
	setHomeScore(0);
	setAwayScore(0);
}

inline ScoreLine::ScoreLine(int homeScore, int awayScore) {
	setHomeScore(homeScore);
	setAwayScore(awayScore);
}

inline std::string ScoreLine::toString() {
	std::stringstream stringStream;
	stringStream << this->homeScore << " - " << this->awayScore;
	return stringStream.str();
}

inline int ScoreLine::getAwayScore() {
    return this->awayScore;
}

inline int ScoreLine::getHomeScore() {
    return this->homeScore;
}

inline void ScoreLine::setAwayScore(int awayScore) {
    this->awayScore = awayScore;
}

inline void ScoreLine::setHomeScore(int homeScore) {
    this->homeScore = homeScore;
}

#endif /* SCORELINE_H_ */
