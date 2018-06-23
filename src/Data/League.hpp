#ifndef LEAGUE_HPP_
#define LEAGUE_HPP_

#include <vector>
#include "../Defn/defn.hpp"
#include "Team.hpp"

class League {
public:
	League(const Season& season);
	friend ostream &operator<<(ostream &os, const League &l);
	void currentMatch(const Team& homeTeam, const Team& awayTeam, const int homeScore, const int awayScore);
	void endMatches(const map<Match, ScoreLine>& scores);
	void updateLeague();
	int matchesPlayedForTeam(const Team& team) const;
	int positionOfTeam(const Team& team) const;

	const IntegerVector& getCurrentGoalDifference() const {
		return currentGoalDifference;
	}

	const IntegerVector& getCurrentGoalsFor() const {
		return currentGoalsFor;
	}

	const IntegerVector& getCurrentPoints() const {
		return currentPoints;
	}

	const vector<Team>& getTeams() const {
		return teams;
	}

	const IntegerVector& getMatchesPlayed() const {
		return matchesPlayed;
	}

private:
	IntegerVector currentPoints, currentGoalsFor, currentGoalDifference;
	IntegerVector endMatchPoints, endMatchGoalsFor, endMatchGoalDifference;
	IntegerVector matchesPlayed;
	vector<Team> teams;
	bool teamsWrongWayRound(const int pos, const int posBelow);
	void swapTeams(const int pos, const int posBelow);
};

inline League::League(const Season& season) {
	currentPoints.resize(20, 0);
	currentGoalsFor.resize(20, 0);
	currentGoalDifference.resize(20, 0);

	endMatchPoints.resize(20, 0);
	endMatchGoalsFor.resize(20, 0);
	endMatchGoalDifference.resize(20, 0);

	matchesPlayed.resize(20, 0);

	teams.reserve(20);

	map<Team, int> teamMap = TEAM_MAP_MAP.at(season);

	for (const auto element : teamMap) {
		teams.push_back(element.first);
	}
}

inline ostream &operator<<(ostream &os, const League &l) {

	static const int w1 = 5;
	static const int w2 = 25;

	os << setw(w1) << left << "pos" << setw(w2) << "team"  << setw(w1) << "pl" << setw(w1) << "pts" << setw(w1) << "gd" << setw(w1)
			<< "gf" << std::endl;
	for (int pos = 0; pos < 20; ++pos) {
		os << setw(w1) << pos + 1 << setw(w2) << l.getTeams()[pos] << setw(w1) <<l. getMatchesPlayed()[pos] << setw(w1) << l.getCurrentPoints()[pos]
				<< setw(w1) << l.getCurrentGoalDifference()[pos] << setw(w1) << l.getCurrentGoalsFor()[pos]
				<< std::endl;
	}
	return os;
}

inline int League::positionOfTeam(const Team& team) const {

	vector<Team>::const_iterator iter = find(teams.begin(), teams.end(), team);

	if (iter != teams.end()) {
		return distance(teams.begin(), iter) + 1;
	} else {
		cout << *this << endl;
		throw runtime_error("Couldn't find position of team " + team);
	}

}

inline int League::matchesPlayedForTeam(const Team& team) const {

	vector<Team>::const_iterator iter = find(teams.begin(), teams.end(), team);

	if (iter != teams.end()) {
		int idx = distance(teams.begin(), iter);
		return matchesPlayed[idx];
	} else {
		cout << *this << endl;
		throw runtime_error("Couldn't find matches played for team " + team);
	}

}

inline void League::currentMatch(const Team& homeTeam, const Team& awayTeam, int homeScore, int awayScore) {

	int homePos = -1, awayPos = -1;
	for (int pos = 0; pos < 20; ++pos) {
		if (homeTeam == teams[pos]) {
			homePos = pos;
		}
		if (awayTeam == teams[pos]) {
			awayPos = pos;
		}
	}

	if (homeScore > awayScore) {
		currentPoints[homePos] = endMatchPoints[homePos] + 3;
		currentPoints[awayPos] = endMatchPoints[awayPos];
	} else if (homeScore == awayScore) {
		currentPoints[homePos] = endMatchPoints[homePos] + 1;
		currentPoints[awayPos] = endMatchPoints[awayPos] + 1;
	} else /*(homeScore < awayScore)*/ {
		currentPoints[homePos] = endMatchPoints[homePos];
		currentPoints[awayPos] = endMatchPoints[awayPos] + 3;
	}
	currentGoalsFor[homePos] = endMatchGoalsFor[homePos] + homeScore;
	currentGoalsFor[awayPos] = endMatchGoalsFor[awayPos] + awayScore;

	currentGoalDifference[homePos] = endMatchGoalDifference[homePos] + (homeScore - awayScore);
	currentGoalDifference[awayPos] = endMatchGoalDifference[awayPos] + (awayScore - homeScore);

	updateLeague();
}

inline void League::endMatches(const map<Match, ScoreLine>& matchScoresMap) {

	for (const auto match : matchScoresMap) {

		Team homeTeam = match.first.getHomeTeam();
		Team awayTeam = match.first.getAwayTeam();

		int homeScore = match.second.first;
		int awayScore = match.second.second;

		int homePos = -1, awayPos = -1;
		for (int pos = 0; pos < 20; ++pos) {
			if (homeTeam == teams[pos]) {
				homePos = pos;
			}
			if (awayTeam == teams[pos]) {
				awayPos = pos;
			}
		}

		if (homeScore > awayScore) {
			endMatchPoints[homePos] += 3;
		} else if (homeScore == awayScore) {
			endMatchPoints[homePos] += 1;
			endMatchPoints[awayPos] += 1;
		} else {
			endMatchPoints[awayPos] += 3;
		}
		endMatchGoalsFor[homePos] += homeScore;
		endMatchGoalsFor[awayPos] += awayScore;

		endMatchGoalDifference[homePos] += (homeScore - awayScore);
		endMatchGoalDifference[awayPos] += (awayScore - homeScore);

		++matchesPlayed[homePos];
		++matchesPlayed[awayPos];
	}
}

inline bool League::teamsWrongWayRound(const int pos, const int posBelow) {
	if (currentPoints[pos] > currentPoints[posBelow]) {
		return false;
	} else if (currentPoints[pos] == currentPoints[posBelow]) {
		//tied on points
		if (currentGoalDifference[pos] > currentGoalDifference[posBelow]) {
			return false;
		} else if (currentGoalDifference[pos] == currentGoalDifference[posBelow]) {
			//tied on points and goal difference
			//don't swap if teams are tied on points, goal difference, and goals for
			//does this ever happen anyway?
			return currentGoalsFor[pos] < currentGoalsFor[posBelow];
		} else /*(currentGoalDifference[pos] < currentGoalDifference[posBelow])*/ {
			return true;
		}
	} else /*(currentPoints[pos] < currentPoints[posBelow])*/ {
		return true;
	}
}

inline void League::swapTeams(const int pos, const int posBelow) {
	swap(teams[pos], teams[posBelow]);
	swap(matchesPlayed[pos], matchesPlayed[posBelow]);

	swap(currentPoints[pos], currentPoints[posBelow]);
	swap(currentGoalsFor[pos], currentGoalsFor[posBelow]);
	swap(currentGoalDifference[pos], currentGoalDifference[posBelow]);

	swap(endMatchPoints[pos], endMatchPoints[posBelow]);
	swap(endMatchGoalsFor[pos], endMatchGoalsFor[posBelow]);
	swap(endMatchGoalDifference[pos], endMatchGoalDifference[posBelow]);
}

inline void League::updateLeague() {
	bool madeSwaps = true;
	while (madeSwaps) {
		madeSwaps = false;
		for (int pos = 0; pos < 19; ++pos) {
			if (teamsWrongWayRound(pos, pos + 1)) {
				madeSwaps = true;
				swapTeams(pos, pos + 1);
			}
		}
	}
}

#endif /* LEAGUE_HPP_ */
