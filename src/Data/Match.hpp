#ifndef MATCH_HPP_
#define MATCH_HPP_

#include <iostream>

#include "../Defn/defn.hpp"

class Match {
public:
	Match(const Team& homeTeam, const Team& awayTeam) {
		setHomeTeam(homeTeam);
		setAwayTeam(awayTeam);
	}

	friend ostream& operator<<(ostream &os, const Match &m) {
		return os << m.getHomeTeam() << " - " << m.getAwayTeam();
	}

	bool operator<(const Match& rhs) const {
		return getHomeTeam() < rhs.getHomeTeam();
	}

	const Team& getAwayTeam() const {
		return awayTeam;
	}

	void setAwayTeam(const Team& awayTeam) {
		this->awayTeam = awayTeam;
	}

	const string& getHomeTeam() const {
		return homeTeam;
	}

	void setHomeTeam(const Team& homeTeam) {
		this->homeTeam = homeTeam;
	}

private:
	Team homeTeam;
	Team awayTeam;
};

//template<> struct less<Match> {
//	bool operator()(const Match& lhs, const Match& rhs) {
//		return lhs.getHomeTeam() < rhs.getHomeTeam();
//	}
//};

#endif /* MATCH_HPP_ */
