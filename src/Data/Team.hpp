#ifndef TEAM_HPP_
#define TEAM_HPP_

#include <map>

#include "../Defn/defn.hpp"

static const vector<Season> ALL_SEASONS = {
		Season("2007-2008"),
		Season("2008-2009"),
		Season("2009-2010"),
		Season("2010-2011"),
		Season("2011-2012")
};

static const map<Team, int> TEAMS_2007_2012 = {
		{make_pair(Team("Arsenal"), 0)},
		{make_pair(Team("Aston Villa"), 1)},
		{make_pair(Team("Birmingham City"), 2)},
		{make_pair(Team("Blackburn Rovers"), 3)},
		{make_pair(Team("Blackpool"), 4)},
		{make_pair(Team("Bolton Wanderers"), 5)},
		{make_pair(Team("Burnley"), 6)},
		{make_pair(Team("Chelsea"), 7)},
		{make_pair(Team("Derby County"), 8)},
		{make_pair(Team("Everton"), 9)},
		{make_pair(Team("Fulham"), 10)},
		{make_pair(Team("Hull City"), 11)},
		{make_pair(Team("Liverpool"), 12)},
		{make_pair(Team("Manchester City"), 13)},
		{make_pair(Team("Manchester United"), 14)},
		{make_pair(Team("Middlesbrough"), 15)},
		{make_pair(Team("Newcastle United"), 16)},
		{make_pair(Team("Norwich City"), 17)},
		{make_pair(Team("Portsmouth"), 18)},
		{make_pair(Team("Queens Park Rangers"), 19)},
		{make_pair(Team("Reading"), 20)},
		{make_pair(Team("Stoke City"), 21)},
		{make_pair(Team("Sunderland"), 22)},
		{make_pair(Team("Swansea City"), 23)},
		{make_pair(Team("Tottenham Hotspur"), 24)},
		{make_pair(Team("West Bromwich Albion"), 25)},
		{make_pair(Team("West Ham United"), 26)},
		{make_pair(Team("Wigan Athletic"), 27)},
		{make_pair(Team("Wolverhampton Wanderers"), 28)}
};

static const map<Team, int> TEAMS_2007_2008 = {
		{make_pair(Team("Arsenal"), 0)},
		{make_pair(Team("Aston Villa"), 1)},
		{make_pair(Team("Birmingham City"), 2)},
		{make_pair(Team("Blackburn Rovers"), 3)},
		{make_pair(Team("Bolton Wanderers"), 4)},
		{make_pair(Team("Chelsea"), 5)},
		{make_pair(Team("Derby County"), 6)},
		{make_pair(Team("Everton"), 7)},
		{make_pair(Team("Fulham"), 8)},
		{make_pair(Team("Liverpool"), 9)},
		{make_pair(Team("Manchester City"), 10)},
		{make_pair(Team("Manchester United"), 11)},
		{make_pair(Team("Middlesbrough"), 12)},
		{make_pair(Team("Newcastle United"), 13)},
		{make_pair(Team("Portsmouth"), 14)},
		{make_pair(Team("Reading"), 15)},
		{make_pair(Team("Sunderland"), 16)},
		{make_pair(Team("Tottenham Hotspur"), 17)},
		{make_pair(Team("West Ham United"), 18)},
		{make_pair(Team("Wigan Athletic"), 19)}
};
static const map<Team, int> TEAMS_2008_2009 = {
		{make_pair(Team("Arsenal"), 0)},
		{make_pair(Team("Aston Villa"), 1)},
		{make_pair(Team("Blackburn Rovers"), 2)},
		{make_pair(Team("Bolton Wanderers"), 3)},
		{make_pair(Team("Chelsea"), 4)},
		{make_pair(Team("Everton"), 5)},
		{make_pair(Team("Fulham"), 6)},
		{make_pair(Team("Hull City"), 7)},
		{make_pair(Team("Liverpool"), 8)},
		{make_pair(Team("Manchester City"), 9)},
		{make_pair(Team("Manchester United"), 10)},
		{make_pair(Team("Middlesbrough"), 11)},
		{make_pair(Team("Newcastle United"), 12)},
		{make_pair(Team("Portsmouth"), 13)},
		{make_pair(Team("Stoke City"), 14)},
		{make_pair(Team("Sunderland"), 15)},
		{make_pair(Team("Tottenham Hotspur"), 16)},
		{make_pair(Team("West Bromwich Albion"), 17)},
		{make_pair(Team("West Ham United"), 18)},
		{make_pair(Team("Wigan Athletic"), 19)}
};
static const map<Team, int> TEAMS_2009_2010 = {
		{make_pair(Team("Arsenal"), 0)},
		{make_pair(Team("Aston Villa"), 1)},
		{make_pair(Team("Birmingham City"), 2)},
		{make_pair(Team("Blackburn Rovers"), 3)},
		{make_pair(Team("Bolton Wanderers"), 4)},
		{make_pair(Team("Burnley"), 5)},
		{make_pair(Team("Chelsea"), 6)},
		{make_pair(Team("Everton"), 7)},
		{make_pair(Team("Fulham"), 8)},
		{make_pair(Team("Hull City"), 9)},
		{make_pair(Team("Liverpool"), 10)},
		{make_pair(Team("Manchester City"), 11)},
		{make_pair(Team("Manchester United"), 12)},
		{make_pair(Team("Portsmouth"), 13)},
		{make_pair(Team("Stoke City"), 14)},
		{make_pair(Team("Sunderland"), 15)},
		{make_pair(Team("Tottenham Hotspur"), 16)},
		{make_pair(Team("West Ham United"), 17)},
		{make_pair(Team("Wigan Athletic"), 18)},
		{make_pair(Team("Wolverhampton Wanderers"), 19)}
};
static const map<Team, int> TEAMS_2010_2011 = {
		{make_pair(Team("Arsenal"), 0)},
		{make_pair(Team("Aston Villa"), 1)},
		{make_pair(Team("Birmingham City"), 2)},
		{make_pair(Team("Blackburn Rovers"), 3)},
		{make_pair(Team("Blackpool"), 4)},
		{make_pair(Team("Bolton Wanderers"), 5)},
		{make_pair(Team("Chelsea"), 6)},
		{make_pair(Team("Everton"), 7)},
		{make_pair(Team("Fulham"), 8)},
		{make_pair(Team("Liverpool"), 9)},
		{make_pair(Team("Manchester City"), 10)},
		{make_pair(Team("Manchester United"), 11)},
		{make_pair(Team("Newcastle United"), 12)},
		{make_pair(Team("Stoke City"), 13)},
		{make_pair(Team("Sunderland"), 14)},
		{make_pair(Team("Tottenham Hotspur"), 15)},
		{make_pair(Team("West Bromwich Albion"), 16)},
		{make_pair(Team("West Ham United"), 17)},
		{make_pair(Team("Wigan Athletic"), 18)},
		{make_pair(Team("Wolverhampton Wanderers"), 19)}
};
static const map<Team, int> TEAMS_2011_2012 = {
		{make_pair(Team("Arsenal"), 0)},
		{make_pair(Team("Aston Villa"), 1)},
		{make_pair(Team("Blackburn Rovers"), 2)},
		{make_pair(Team("Bolton Wanderers"), 3)},
		{make_pair(Team("Chelsea"), 4)},
		{make_pair(Team("Everton"), 5)},
		{make_pair(Team("Fulham"), 6)},
		{make_pair(Team("Liverpool"), 7)},
		{make_pair(Team("Manchester City"), 8)},
		{make_pair(Team("Manchester United"), 9)},
		{make_pair(Team("Newcastle United"), 10)},
		{make_pair(Team("Norwich City"), 11)},
		{make_pair(Team("Queens Park Rangers"), 12)},
		{make_pair(Team("Stoke City"), 13)},
		{make_pair(Team("Sunderland"), 14)},
		{make_pair(Team("Swansea City"), 15)},
		{make_pair(Team("Tottenham Hotspur"), 16)},
		{make_pair(Team("West Bromwich Albion"), 17)},
		{make_pair(Team("Wigan Athletic"), 18)},
		{make_pair(Team("Wolverhampton Wanderers"), 19)}
};
static const map<Season, map<Team, int>> TEAM_MAP_MAP = {
		{make_pair(Season("2007-2008"), TEAMS_2007_2008)},
		{make_pair(Season("2008-2009"), TEAMS_2008_2009)},
		{make_pair(Season("2009-2010"), TEAMS_2009_2010)},
		{make_pair(Season("2010-2011"), TEAMS_2010_2011)},
		{make_pair(Season("2011-2012"), TEAMS_2011_2012)}
};

#endif /* TEAM_HPP_ */
