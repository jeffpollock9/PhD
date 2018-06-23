
matches <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/matches.csv", stringsAsFactors = FALSE)
events <- read.csv("/home/jeff/workspace/Latex-Particle-Filter/data/results-simulated.csv", stringsAsFactors = FALSE)

matchIDs <- character(380)
g <- 1L
for (h in 0:19) {
	for (a in 0:19) {
		if (h != a) {
			matchIDs[g] <- paste(h, a, sep = "-")
			g <- g + 1L
		}
	}
}

minBeforeHalfTime <- data.frame(matchID = matchIDs, 
		time = 44.0 / 90.0, 
		event = "break", 
		stringsAsFactors = FALSE)
halfTime <- data.frame(matchID = matchIDs, 
		time = 45.0 / 90.0, 
		event = "break", 
		stringsAsFactors = FALSE)
minBeforeEndMatch <- data.frame(matchID = matchIDs, 
		time = 89.0 / 90.0, 
		event = "break", 
		stringsAsFactors = FALSE)
endMatch <- data.frame(matchID = matchIDs, 
		time = 90.0 / 90.0, 
		event = "break", 
		stringsAsFactors = FALSE)

goalAndRedCardTimes <- rbind(events, minBeforeHalfTime, halfTime, minBeforeEndMatch, endMatch)
goalAndRedCardTimes <- goalAndRedCardTimes[order(goalAndRedCardTimes$matchID, goalAndRedCardTimes$time), ]
rownames(goalAndRedCardTimes) <- NULL

goalAndRedCardTimes$start<- unlist(aggregate(time ~ matchID, goalAndRedCardTimes, na.action = na.pass,
				function(t) {
					c(0L, t[1:(length(t) - 1L)])
				}
		)$time)

countEvent <- function(df, eventID) {
	unlist(aggregate(event ~ matchID, df, 
					function(x) {
						c(0L, cumsum(x[-length(x)] == eventID))
					}
			)$event)
}

goalAndRedCardTimes <- within(goalAndRedCardTimes, c(
				homeScore <- countEvent(goalAndRedCardTimes, "home"),
				awayScore <- countEvent(goalAndRedCardTimes, "away"))
)

colnames(goalAndRedCardTimes)[2L] <- "end"

goalAndRedCardTimes$isHomeGoal <- (goalAndRedCardTimes$event == "home") * 1
goalAndRedCardTimes$isAwayGoal <- (goalAndRedCardTimes$event == "away") * 1

library(sqldf)

matches <- matches[matches$competitionid == 8L, ]
matches <- matches[order(matches$matchid), ]
matches$homeID <- unclass(as.factor(matches$hometeamname)) - 1L
matches$awayID <- unclass(as.factor(matches$awayteamname)) - 1L
matches$matchID <- paste0(matches$homeID, "-", matches$awayID)
matches$week <- rep(1:38, each = 10L)
matches$matchid <- NULL

data <- sqldf("
				SELECT g.matchID, g.start, g.end, g.homeScore, g.awayScore, 
				g.isHomeGoal, g.isAwayGoal, m.homeID, m.awayID, m.week
				FROM matches m 
				JOIN goalAndRedCardTimes g
				ON g.matchID = m.matchID
				")

makeVector <- function(name, type, x) {
	
	n <- length(x)
	if (type == "double") {
		d <- paste(sprintf("%0.1f / 90.0", x), collapse = ", ")
	} else {
		d <- paste(x, collapse = ", ")
	}
	string <- sprintf("static const %s %s[%i] = {%s};\n", type, name, n, d)
	
	invisible(string)
}

{
	cat(makeVector("START_SIM", "doublem", data$start))
	cat(makeVector("END_SIM", "doublem", data$end))
	cat(makeVector("IS_HOME_GOAL_SIM", "bool", data$isHomeGoal))
	cat(makeVector("IS_AWAY_GOAL_SIM", "bool", data$isAwayGoal))
	cat(makeVector("HOME_SCORE_SIM", "int", data$homeScore))
	cat(makeVector("AWAY_SCORE_SIM", "int", data$awayScore))
	cat(makeVector("HOME_ID_SIM", "int", data$homeID))
	cat(makeVector("AWAY_ID_SIM", "int", data$awayID))
	cat(makeVector("WEEK_SIM", "int", data$week))
	cat("\n")
}

simScores <- aggregate(cbind(isHomeGoal, isAwayGoal) ~ matchID, data, sum)

simMatches <- sqldf("
				SELECT m.matchID, m.week, m.homeID, m.awayID, 
				s.isHomeGoal AS homeScore, s.isAwayGoal AS awayScore
				FROM matches m
				JOIN simScores s
				ON s.matchID = m.matchID
				ORDER BY m.week
				")

{
	cat(makeVector("MATCH_HOME_SCORE_SIM", "int", simMatches$homeScore))
	cat(makeVector("MATCH_AWAY_SCORE_SIM", "int", simMatches$awayScore))
	cat(makeVector("MATCH_HOME_ID_SIM", "int", simMatches$homeID))
	cat(makeVector("MATCH_AWAY_ID_SIM", "int", simMatches$awayID))
	cat(makeVector("MATCH_WEEK_SIM", "int", simMatches$week))
	cat("\n")
}

### team resource priors
R <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/R-simulated.csv", header = FALSE)
cat(sprintf("double Rprior[20] = { %s };\n", paste(R[1L, ], collapse = ", ")))
