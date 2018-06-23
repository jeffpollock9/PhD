#ifndef EVENT_HPP_
#define EVENT_HPP_

#include <map>
#include <string>

#include "../Defn/defn.hpp"

map<int, string> KEY = {
		{ -1, string("Nothing") },
		{ 1, string("Home Goal") },
		{ 2, string("Away Goal") },
		{ 3, string("Home Red Card") },
		{ 4, string("Away Red Card") }
};

class Event {
public:
	Event(const Team& homeTeam, const Team& awayTeam, const int number, const double time) {
		setHomeTeam(homeTeam);
		setAwayTeam(awayTeam);
		setNumber(number);
		setTime(time);
	}

	friend ostream &operator<<(ostream &os, const Event &e) {
		return os << e.getHomeTeam() << " - " << e.getAwayTeam() << ", Event: " << e.getNumber() << ", Time: " << e.getTime();
	}

	const Team& getAwayTeam() const {
		return awayTeam;
	}

	void setAwayTeam(const Team& awayTeam) {
		this->awayTeam = awayTeam;
	}

	const Team& getHomeTeam() const {
		return homeTeam;
	}

	void setHomeTeam(const Team& homeTeam) {
		this->homeTeam = homeTeam;
	}

	double getTime() const {
		return time;
	}

	void setTime(const double time) {
		this->time = time;
	}

	int getNumber() const {
		return number;
	}

	void setNumber(const int number) {
		this->number = number;
	}

private:
	Team homeTeam;
	Team awayTeam;
	int number;
	double time;
};

#endif /* EVENT_HPP_ */
