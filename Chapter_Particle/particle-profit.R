library(sqldf)

matches <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/matches.csv", stringsAsFactors = FALSE)
matches <- subset(matches, competitionid == 8L)
matches <- matches[sort.list(matches$matchid), ]
finalScoresList <- strsplit(matches$score, " - ")
matches$finalHomeScore <- sapply(finalScoresList, function(x) as.integer(x[1]))
matches$finalAwayScore <- sapply(finalScoresList, function(x) as.integer(x[2]))
matches$competition <- matches$competitionid <- matches$score <- NULL
matches$homeID <- as.integer(unclass(factor(matches$hometeamname))) - 1L
matches$awayID <- as.integer(unclass(factor(matches$awayteamname))) - 1L
matches$week <- rep(1:38, each = 10)


odds <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/E0-2011-2012.csv", stringsAsFactors = FALSE)
namesLink <- data.frame(matches = sort(unique(matches$hometeamname)), odds = sort(unique(odds$HomeTeam)))

joinedNamesData <- sqldf("
				SELECT m.*, l1.odds AS oddsHome, l2.odds AS oddsAway 
				FROM matches m 
				JOIN namesLink l1
				ON m.hometeamname = l1.matches
				JOIN namesLink l2
				ON m.awayteamname = l2.matches
				")

fullData <- sqldf("
				SELECT t1.matchid AS matchID, t1.hometeamname AS homeTeam, t1.awayteamname AS awayTeam, 
				t1.finalHomeScore, t1.finalAwayScore,
				t2.B365H, t2.B365D, t2.B365A
				FROM joinedNamesData t1
				JOIN odds t2 
				ON t1.oddsHome = t2.HomeTeam AND t1.oddsAway = t2.AwayTeam 
				ORDER BY t1.matchid
				")

modelProbsDyn <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/dynamic-probs-openmp.csv", header = FALSE)
modelProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsM.csv", header = FALSE)

colnames(modelProbsDyn)[1:3] <- c("home", "draw", "away")
colnames(modelProbs)[1:3] <- c("home", "draw", "away")

bookieOdds <- cbind(fullData, modelProbsDyn[ , 1:3])

bookieOdds <- tail(bookieOdds, 330)

d <- within(bookieOdds, c(
				betHome <- home > 1 / B365H,
				betDraw <- draw > 1 / B365D,
				betAway <- away > 1 / B365A,
				profitHome <- ifelse(betHome, (finalHomeScore > finalAwayScore) * B365H - 1.0, 0.0),
				profitDraw <- ifelse(betDraw, (finalHomeScore == finalAwayScore) * B365D - 1.0, 0.0),
				profitAway <- ifelse(betAway, (finalHomeScore < finalAwayScore) * B365A - 1.0, 0.0)
		)
)

profitBettingHome <- sum(d$profitHome, na.rm = TRUE)
profitBettingDraw <- sum(d$profitDraw, na.rm = TRUE)
profitBettingAway <- sum(d$profitAway, na.rm = TRUE)

nBetsHome <- sum(d$betHome, na.rm = TRUE)
nBetsDraw <- sum(d$betDraw, na.rm = TRUE)
nBetsAway <- sum(d$betAway, na.rm = TRUE)

nBets <- nBetsHome + nBetsDraw + nBetsAway
profit <- profitBettingHome + profitBettingDraw + profitBettingAway

LSR <- 0.0

for (i in 1:nrow(d)) {
	if (d$finalHomeScore[i] > d$finalAwayScore[i]) {
		LSR <- LSR + log(d$home[i])
	} else if (d$finalHomeScore[i] == d$finalAwayScore[i]) {
		LSR <- LSR + log(d$draw[i])
	} else {
		LSR <- LSR + log(d$away[i])
	}
}

print(nBets)
print(profit)
print(LSR)

par(mfrow = c(3, 1))
plot(modelProbs[, 1], modelProbsDyn[, 1])
abline(0, 1, col = "red", lty = 2)
plot(modelProbs[, 2], modelProbsDyn[, 2])
abline(0, 1, col = "red", lty = 2)
plot(modelProbs[, 3], modelProbsDyn[, 3])
abline(0, 1, col = "red", lty = 2)


