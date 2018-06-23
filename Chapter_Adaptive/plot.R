library(ggplot2)
library(gridExtra)
library(reshape2)
library(sqldf)
library(RColorBrewer)

rm(list = ls())

##########
##########

priorSamples <- read.csv("/home/jeff/workspace/C-Football/csv/priorSamples-1.0-0.5-thinned.csv", 
		header = FALSE)

names <- c("h", "a", "c_-1", "c_0", "c_1", "d_-1", "d_0", 
		"d_1", "rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Bromwich Albion", "Wigan Athletic", "Wolverhampton Wanderers")

colnames(priorSamples) <- names[11:30]

pl1 <- ggplot(priorSamples) + geom_histogram(aes(`Manchester United`, ..density..), binwidth = 0.05) + 
		xlab(expression(R["Manchester United"])) + xlim(0, 3) + ylim(0, 5) + theme_bw()
pl2 <- ggplot(priorSamples) + geom_histogram(aes(`Arsenal`, ..density..), binwidth = 0.05) + 
		xlab(expression(R[Arsenal])) + xlim(0, 3) + ylim(0, 5) + theme_bw()
pl3 <- ggplot(priorSamples) + geom_histogram(aes(`Everton`, ..density..), binwidth = 0.05) + 
		xlab(expression(R[Everton])) + xlim(0, 3) + ylim(0, 5) + theme_bw()
pl4 <- ggplot(priorSamples) + geom_histogram(aes(`Wolverhampton Wanderers`, ..density..), binwidth = 0.05) + 
		xlab(expression(R["Wolverhampton Wanderers"])) + xlim(0, 3) + ylim(0, 5) + theme_bw()

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/marginal-prior.pdf")
grid.arrange(pl1, pl2, pl3, pl4, ncol = 2, nrow = 2)
dev.off()

##########
##########

cols <- colorRampPalette(c("white", "red"))(16)

g1 <- ggplot(priorSamples, aes(x = `Chelsea`, y = `Manchester United`)) + stat_binhex(bins = 50) + 
		theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = 3) + 
		scale_fill_gradient(low = cols[2], high = cols[16], na.value = "white", guide = "none") + 
		xlab(expression(R["Chelsea"])) + ylab(expression(R["Manchester United"])) + 
		xlim(c(0, 3)) + ylim(c(0, 3))

g2 <- ggplot(priorSamples, aes(x = `Wolverhampton Wanderers`, y = `Everton`)) + stat_binhex(bins = 50) + 
		theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = 3) + 
		scale_fill_gradient(low = cols[2], high = cols[16], na.value = "white", guide = "none") + 
		xlab(expression(R["Wolverhampton Wanderers"])) + ylab(expression(R["Everton"])) + 
		xlim(c(0, 0.8)) + ylim(c(0, 0.8))

g3 <- ggplot(priorSamples, aes(x = `Bolton Wanderers`, y = `Newcastle United`)) + stat_binhex(bins = 50) + 
		theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = 3) + 
		scale_fill_gradient(low = cols[2], high = cols[16], na.value = "white", guide = "none") + 
		xlab(expression(R["Bolton Wanderers"])) + ylab(expression(R["Newcastle United"])) + 
		xlim(c(0, 0.8)) + ylim(c(0, 0.8))

g4 <- ggplot(priorSamples, aes(x = `Swansea City`, y = `Queens Park Rangers`)) + stat_binhex(bins = 50) + 
		theme_bw() + geom_abline(intercept = 0, slope = 1, linetype = 3) + 
		scale_fill_gradient(low = cols[2], high = cols[16], na.value = "white", guide = "none") + 
		xlab(expression(R["Swansea City"])) + ylab(expression(R["Queens Park Rangers"])) + 
		xlim(c(0, 0.8)) + ylim(c(0, 0.8))

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/bivariate-prior.pdf")
grid.arrange(g1, g2, g3, g4, ncol = 2, nrow = 2)
dev.off()

##########
##########

posteriorSamplesNonLinear <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/R/MCMC.csv",
		header = FALSE)[-(1:5000), ]

namesNonLinear <- c("h", "a", "c_-1", "c_0", "c_1", "e_-1", "e_0", "e_1", "d_-1", "d_0", 
		"d_1", "rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Bromwich Albion", "Wigan Athletic", "Wolverhampton Wanderers")

colnames(posteriorSamplesNonLinear) <- namesNonLinear

xlims <- c(-1, 1)
ylims <- c(0, 3.5)
breaks <- seq(-2, 2, by = 0.05)

nonLinearHist1 <- ggplot(posteriorSamplesNonLinear, aes(x = `e_-1` - 0.5 * (`c_-1` + `d_-1`), y = ..density..)) + 
		geom_histogram(breaks = breaks) + theme_bw() + xlab(expression(e[-1] - 0.5 (c[-1] + d[-1]))) + xlim(xlims) + 
		ylim(ylims)

nonLinearHist2 <- ggplot(posteriorSamplesNonLinear, aes(x = `e_0` - 0.5 * (`c_0` + `d_0`), y = ..density..)) + 
		geom_histogram(breaks = breaks) + theme_bw() + xlab(expression(e[0] - 0.5 (c[0] + d[0]))) + xlim(xlims) + 
		ylim(ylims) 

nonLinearHist3 <- ggplot(posteriorSamplesNonLinear, aes(x = `e_1` - 0.5 * (`c_1` + `d_1`), y = ..density..)) + 
		geom_histogram(breaks = breaks) + theme_bw() + xlab(expression(e[1] - 0.5 (c[1] + d[1]))) + xlim(xlims) + 
		ylim(ylims)

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/non-linear.pdf", height = 9)
grid.arrange(nonLinearHist1, nonLinearHist2, nonLinearHist3, ncol = 1)
dev.off()

##########
##########

posteriorSamplesPiecewiseLinear <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/data/MCMC_piecewise_linear.csv",
		header = FALSE, colClasses = "numeric")[-(1:5000), ]

namesPiecewiseLinear <- c("h", "a", 
		paste0("K_-1_", 0:6 * 15), paste0("K_0_", 0:6 * 15), paste0("K_1_", 0:6 * 15), 
		"rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Bromwich Albion", "Wigan Athletic", "Wolverhampton Wanderers")

colnames(posteriorSamplesPiecewiseLinear) <- namesPiecewiseLinear

meanData <- colMeans(posteriorSamplesPiecewiseLinear)
upperData <- apply(posteriorSamplesPiecewiseLinear, 2L, function(x) quantile(x, 0.025))
lowerData <- apply(posteriorSamplesPiecewiseLinear, 2L, function(x) quantile(x, 0.975))

t <- 0:6 * 15
nIterations <- nrow(posteriorSamplesPiecewiseLinear)
posteriorSamplesPiecewiseLinear$iteration <- 1:nIterations

nSamples <- 1500L
s <- sort(sample.int(nIterations, nSamples, FALSE))

d1 <- melt(posteriorSamplesPiecewiseLinear[s , c(3:9, 46)], id.vars = "iteration")
d2 <- melt(posteriorSamplesPiecewiseLinear[s , c(10:16, 46)], id.vars = "iteration")
d3 <- melt(posteriorSamplesPiecewiseLinear[s , c(17:23, 46)], id.vars = "iteration")

d1$x <- rep(t / 90.0, each = nSamples)
d2$x <- rep(t / 90.0, each = nSamples)
d3$x <- rep(t / 90.0, each = nSamples)

m1 <- data.frame(t = t / 90.0, value = meanData[3:9])
m2 <- data.frame(t = t / 90.0, value = meanData[10:16])
m3 <- data.frame(t = t / 90.0, value = meanData[17:23])

l1 <- data.frame(t = t / 90.0, value = lowerData[3:9])
l2 <- data.frame(t = t / 90.0, value = lowerData[10:16])
l3 <- data.frame(t = t / 90.0, value = lowerData[17:23])

u1 <- data.frame(t = t / 90.0, value = upperData[3:9])
u2 <- data.frame(t = t / 90.0, value = upperData[10:16])
u3 <- data.frame(t = t / 90.0, value = upperData[17:23])

piecewiseLinearHist1 <- ggplot() + 
		geom_line(aes(x, value, group = iteration), d1, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(t, value), m1, col = "red") + 
		geom_line(aes(t, value), l1, col = "black", linetype = 2) + 
		geom_line(aes(t, value), u1, col = "black", linetype = 2) + 
		theme_bw() + ylab(expression(alpha[k](t))) + xlab("t") + xlim(c(0, 1)) + ylim(c(0, 1))

piecewiseLinearHist2 <- ggplot() + 
		geom_line(aes(x, value, group = iteration), d2, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(t, value), m2, col = "red") + 
		geom_line(aes(t, value), l2, col = "black", linetype = 2) + 
		geom_line(aes(t, value), u2, col = "black", linetype = 2) + 
		theme_bw() + ylab(expression(alpha[k](t))) + xlab("t") + xlim(c(0, 1)) + ylim(c(0, 1))

piecewiseLinearHist3 <- ggplot() + 
		geom_line(aes(x, value, group = iteration), d3, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(t, value), m3, col = "red") + 
		geom_line(aes(t, value), l3, col = "black", linetype = 2) + 
		geom_line(aes(t, value), u3, col = "black", linetype = 2) + 
		theme_bw() +ylab(expression(alpha[k](t))) + xlab("t") + xlim(c(0, 1)) + ylim(c(0, 1))

#pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/piecewise-linear.pdf", height = 9)
png("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/piecewise-linear.png")
grid.arrange(piecewiseLinearHist1, piecewiseLinearHist2, piecewiseLinearHist3, ncol = 1)
dev.off()

##########
##########

DrProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsDR.csv", header = FALSE)
BtProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsBT.csv", header = FALSE)
BBtProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsBBT.csv", header = FALSE)
MProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsM.csv", header = FALSE)

DrProbs <- tail(DrProbs, 330)
BtProbs <- tail(BtProbs, 330)
MProbs <- tail(MProbs, 330)
BBtProbs <- tail(BBtProbs, 330)

colnames(DrProbs)[1:3] <- paste0("DR_", c("home", "draw", "away"))
colnames(BtProbs)[1:3] <- paste0("BT_", c("home", "draw", "away"))
colnames(BBtProbs)[1:3] <- paste0("BBT_", c("home", "draw", "away"))
colnames(MProbs)<- c(paste0("M_", c("home", "draw", "away")), "home", "away")

probs <- cbind(DrProbs[1:3], BtProbs[1:3], BBtProbs[1:3], MProbs)
odds <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/E0-2011-2012.csv", stringsAsFactors = FALSE)
probs$order <- seq_len(nrow(probs))

odds$home <- unclass(as.factor(odds$HomeTeam)) - 1L
odds$away <- unclass(as.factor(odds$AwayTeam)) - 1L

library(sqldf)

data <- sqldf("
				SELECT p.*, o.FTHG AS homeScore, o.FTAG AS awayScore
				FROM probs p
				JOIN odds o 
				ON o.home = p.home AND o.away = p.away 
				")

data$DRLogLik <- with(data, log(DR_home) * (homeScore > awayScore) + 
				log(DR_draw) * (homeScore == awayScore) + 
				log(DR_away) * (homeScore < awayScore))
#has a log(1) should be 0
data$DRLogLik[is.nan(data$DRLogLik)] <- log(1.0)

data$MLogLik <- with(data, log(M_home) * (homeScore > awayScore) + 
				log(M_draw) * (homeScore == awayScore) + 
				log(M_away) * (homeScore < awayScore))

data$BTLogLik <- with(data, log(BT_home) * (homeScore > awayScore) + 
				log(BT_draw) * (homeScore == awayScore) + 
				log(BT_away) * (homeScore < awayScore))

data$BBTLogLik <- with(data, log(BBT_home) * (homeScore > awayScore) + 
				log(BBT_draw) * (homeScore == awayScore) + 
				log(BBT_away) * (homeScore < awayScore))

weeks <- (39 - nrow(data) / 10):38
data$week <- rep(weeks, each = 10)

startWeek <- 6

agg <- aggregate(cbind(DRLogLik, MLogLik, BTLogLik, BBTLogLik) ~ week, data, sum)
#agg <- aggregate(cbind(DRLogLik - MLogLik) ~ week, data, sum)
agg1 <- agg[agg$week >= startWeek, ]

n <- 100000
sums <- numeric(n)

for (i in 1:n) {
	samp <- sample.int(4, 33, TRUE)
	for (j in 1:33) {
		sums[i] <- sums[i] + agg1[j, 1L + samp[j]] 
	}
}

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/scoring-rule-sample.pdf")
dat <- colSums(agg1)[-1][c(2, 1, 3, 4)] # so in right order with thesis
datCols <- c("blue", "red", "green", "black")
p <- sapply(dat, function(x) mean(x < sums))
ggplot() + geom_histogram(aes(sums, ..density..), binwidth = 0.5) +
		scale_x_continuous(name = expression(Sigma["w=6"]^"38"*LSR[w]), breaks = seq(-360, -315, by = 5), limits = c(-360, -320)) + 
		geom_line(aes(x = rep(dat[1], 2)), y = c(0.0, 10.0), size = 1, col = datCols[1]) + 
		geom_line(aes(x = rep(dat[2], 2)), y = c(0.0, 10.0), size = 1, col = datCols[2]) + 
		geom_line(aes(x = rep(dat[3], 2)), y = c(0.0, 10.0), size = 1, col = datCols[3]) + 
		geom_line(aes(x = rep(dat[4], 2)), y = c(0.0, 10.0), size = 1, col = datCols[4]) + 
		theme_bw()
dev.off()

sums2 <- numeric(n)

for (i in 1:n) {
	samp <- sample(c(1L, 2L, 4L), 33, TRUE)
	for (j in 1:33) {
		sums2[i] <- sums2[i] + agg1[j, 1L + samp[j]] 
	}
}

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/scoring-rule-sample2.pdf")
dat <- colSums(agg1)[-1][c(2, 1, 3, 4)] # so in right order with thesis
datCols <- c("blue", "red", "green", "black")
p2 <- sapply(dat, function(x) mean(x < sums2))
ggplot() + geom_histogram(aes(sums2, ..density..), binwidth = 0.5) +
		scale_x_continuous(name = expression(Sigma["w=6"]^"38"*LSR[w]), breaks = seq(-345, -320, by = 5), limits = c(-345, -320)) + 
		geom_line(aes(x = rep(dat[1], 2)), y = c(0.0, 10.0), size = 1, col = datCols[1]) + 
		geom_line(aes(x = rep(dat[2], 2)), y = c(0.0, 10.0), size = 1, col = datCols[2]) + 
		geom_line(aes(x = rep(dat[4], 2)), y = c(0.0, 10.0), size = 1, col = datCols[4]) + 
		theme_bw()
dev.off()

##########
##########

posteriorSamples <- read.csv("/home/jeff/workspace/C-Football/csv/history.csv", 
		header = FALSE, colClasses = "numeric")

#temp
apply(posteriorSamples, 2, function(x) {
			n <- length(x)
			corr <- acf(x, plot = FALSE, lag.max = 100)
			corr1 <- as.numeric(corr$acf)
			n / (1 + 2 * sum(corr1))
		})

names <- c("h", "a", "c_-1", "c_0", "c_1", "d_-1", "d_0", 
		"d_1", "rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Bromwich Albion", "Wigan Athletic", "Wolverhampton Wanderers")

colnames(posteriorSamples) <- names

teams <- names[11:30]
teamData <- posteriorSamples[ , 11:30]
d <- data.frame(R = unlist(teamData), 
		team = factor(rep(teams, each = nrow(teamData)), 
				levels = teams[order(colMeans(teamData), decreasing = FALSE)]))

Rposterior <- ggplot(d, aes(team, R)) + geom_violin(fill = "red") + coord_flip() + theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/Rposterior.pdf", Rposterior)

##########
##########

rankings <- matrix(NA_integer_, nrow(teamData), ncol(teamData))
colnames(rankings) <- colnames(teamData)

for (i in seq_len(nrow(teamData))) {
	rankings[i, ] <- 21L - rank(teamData[i,])
}

tab <- apply(rankings, 2L, function(x) table(factor(x, levels = 1:20)))

positions <- c(1, 2, 6, 3, 4, 8, 7, 9, 13, 16, 5, 10, 14, 18, 19, 15, 11, 12, 17, 20)
pos <- data.frame(x = positions, y = 1:20)

grid1 <- expand.grid(1:21 - 0.5, c(0.5, 20.5))
grid1$g <- grid1$Var1

grid2 <- expand.grid(c(0.5, 20.5), 1:21 - 0.5)
grid2$g <- grid2$Var2 + 100

grid <- rbind(grid1, grid2)

rankData <- as.data.frame(tab[ , order(colMeans(rankings))] / nrow(rankings))
rankData$position <- 1:20

rankDataMelted <- melt(rankData, id.vars = "position", variable.name = "team", value.name = "prob")

cols <- colorRampPalette(c("white", "red"))(16)

rankDataMelted <- melt(rankData, id.vars = "position", variable.name = "team", value.name = "prob")
rankDataMelted$prob[rankDataMelted$prob < 0.001] <- NA

rankingsPlot <- ggplot(rankDataMelted) + 
		geom_tile(aes(x = position, y = team, fill = prob), stat = "identity") +
		scale_fill_gradient(low = cols[2], high = cols[16], na.value = "white", guide = "none") + 
		geom_point(aes(x, y), pos, size = 3) + 
		geom_line(aes(Var1, Var2, group = g), grid, col = "lightgrey", size = 0.1) + 
		theme_classic() + 
		scale_x_discrete(breaks = 1:20) + 
		theme(panel.grid.major = element_blank())

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/rankings.pdf", rankingsPlot)

##########
##########

graphNames <- as.list(names)
graphNames[3:5] <- list(expression(c[-1]), expression(c[0]), expression(c[1]))
graphNames[6:8] <- list(expression(d[-1]), expression(d[0]), expression(d[1]))
graphNames[9:10] <- list(expression(rho[1]), expression(rho[2]))


colnames(posteriorSamples)[seq_along(names)] <- paste0("col", seq_along(names))
posteriorSamples$iteration <- seq_len(nrow(posteriorSamples))

#1
plotList <- vector("list", 10L)
for (i in 1:5) {
	plotList[[i * 2 - 1]] <- ggplot(posteriorSamples, aes(x = iteration)) + 
			geom_line(aes_string(y = paste0("col", i))) + 
			ylab(graphNames[[i]]) + theme_bw()
	
	plotList[[i * 2]] <- ggplot(posteriorSamples, aes(y = ..density..)) + 
			geom_histogram(aes_string(x = paste0("col", i)), binwidth = 0.02) + 
			xlab(graphNames[[i]]) + theme_bw()
}
pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/MCMCplots1.pdf", height = 9)
do.call(grid.arrange, c(plotList, ncol = 2))
dev.off()

#2
plotList <- vector("list", 10L)
for (i in 1:5) {
	j <- i + 5L
	plotList[[i * 2 - 1]] <- ggplot(posteriorSamples, aes(x = iteration)) + 
			geom_line(aes_string(y = paste0("col", j))) + 
			ylab(graphNames[[j]]) + theme_bw()
	
	plotList[[i * 2]] <- ggplot(posteriorSamples, aes(y = ..density..)) + 
			geom_histogram(aes_string(x = paste0("col", j)), binwidth = 0.02) + 
			xlab(graphNames[[j]]) + theme_bw()
}
pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/MCMCplots2.pdf", height = 9)
do.call(grid.arrange, c(plotList, ncol = 2))
dev.off()

posteriorSamples$iteration <- NULL
colnames(posteriorSamples) <- names

##########
##########

alphaStar <- function(t) {
	t^2 * (2*c + 2 * d - 4 * e) + t * (4 * e - d - 3 * c) + c
}

alpha <- function(t) {
	c + (d - c) * t
}

n <- 1500L
m <- 100L
s <- sample.int(nrow(posteriorSamples), n, FALSE)
t <- seq(0.0, 1.0, len = m)

plotList <- vector("list", 9)

#alpha
for (j in 1:3) {
	
	l <- vector("list", n)
	
	for (i in seq_len(n)) {
		c <- posteriorSamples[s[i], j + 2]
		d <- posteriorSamples[s[i], j + 5]
		l[[i]] <- data.frame(t, alpha = alpha(t), sample = i)
	}
	
	alphaData <- do.call(rbind, l)
	
	c <- mean(posteriorSamples[ , j + 2])
	d <- mean(posteriorSamples[ , j + 5])
	meanData <- data.frame(t = c(0.0, 1.0), alpha = c(alpha(0.0), alpha(1.0)))
	
	plotList[[j]] <- ggplot(alphaData) + 
			geom_line(aes(x = t, y = alpha, group = sample), col = rgb(0, 0, 1, 0.1)) + 
			geom_line(aes(x = t, y = alpha), meanData, col = "red", size = 1.0) + 
			xlim(c(0, 1)) + ylim(c(0, 1)) + 
			labs(y = expression(alpha[k](t)), title = switch(j, "losing", "drawing", "winning")) + 
			theme_bw()
}

#alpha(1)
for (j in 1:3) {
	
	l <- vector("list", n)
	
	for (i in seq_len(n)) {
		c <- posteriorSamplesNonLinear[s[i], j + 2]
		e <- posteriorSamplesNonLinear[s[i], j + 5]
		d <- posteriorSamplesNonLinear[s[i], j + 8]
		l[[i]] <- data.frame(t, alpha = alphaStar(t), sample = i)
	}
	
	alphaData <- do.call(rbind, l)
	
	c <- mean(posteriorSamplesNonLinear[ , j + 2])
	e <- mean(posteriorSamplesNonLinear[ , j + 5])
	d <- mean(posteriorSamplesNonLinear[ , j + 8])
	meanData <- data.frame(t, alpha = alphaStar(t))
	
	plotList[[j + 3L]] <- ggplot(alphaData) + 
			geom_line(aes(x = t, y = alpha, group = sample), col = rgb(0, 0, 1, 0.1)) + 
			geom_line(aes(x = t, y = alpha), meanData, col = "red", size = 1.0) + 
			xlim(c(0, 1)) + ylim(c(0, 1)) + 
			labs(y = expression(alpha[k]^(1)*(t)), title = switch(j, "losing", "drawing", "winning")) + 
			theme_bw()
}

#alpha(2)
posteriorSamplesPiecewiseLinear2 <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/data/MCMC_piecewise_linear2.csv",
		header = FALSE, colClasses = "numeric")[-(1:5000), ]

namesPiecewiseLinear2 <- c("h", "a", 
		paste0("K_-1_", 0:3 * 30), paste0("K_0_", 0:3 * 15), paste0("K_1_", 0:3 * 15), 
		"rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Bromwich Albion", "Wigan Athletic", "Wolverhampton Wanderers")

colnames(posteriorSamplesPiecewiseLinear2) <- namesPiecewiseLinear2

meanData <- colMeans(posteriorSamplesPiecewiseLinear2)
upperData <- apply(posteriorSamplesPiecewiseLinear2, 2L, function(x) quantile(x, 0.025))
lowerData <- apply(posteriorSamplesPiecewiseLinear2, 2L, function(x) quantile(x, 0.975))

t <- 0:3 * 30
nIterations <- nrow(posteriorSamplesPiecewiseLinear2)
posteriorSamplesPiecewiseLinear2$iteration <- 1:nIterations

nSamples <- 1500L
s <- sort(sample.int(nIterations, nSamples, FALSE))

d1 <- melt(posteriorSamplesPiecewiseLinear2[s , c(3:6, 37)], id.vars = "iteration")
d2 <- melt(posteriorSamplesPiecewiseLinear2[s , c(7:10, 37)], id.vars = "iteration")
d3 <- melt(posteriorSamplesPiecewiseLinear2[s , c(11:14, 37)], id.vars = "iteration")

d1$x <- rep(t / 90.0, each = nSamples)
d2$x <- rep(t / 90.0, each = nSamples)
d3$x <- rep(t / 90.0, each = nSamples)

m1 <- data.frame(t = t / 90.0, value = meanData[3:6])
m2 <- data.frame(t = t / 90.0, value = meanData[7:10])
m3 <- data.frame(t = t / 90.0, value = meanData[11:14])

l1 <- data.frame(t = t / 90.0, value = lowerData[3:6])
l2 <- data.frame(t = t / 90.0, value = lowerData[7:10])
l3 <- data.frame(t = t / 90.0, value = lowerData[11:14])

u1 <- data.frame(t = t / 90.0, value = upperData[3:6])
u2 <- data.frame(t = t / 90.0, value = upperData[7:10])
u3 <- data.frame(t = t / 90.0, value = upperData[11:14])

piecewiseLinearHist21 <- ggplot() + 
		geom_line(aes(x, value, group = iteration), d1, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(t, value), m1, col = "red", size = 1.0) + ggtitle("losing") + 
		theme_bw() + ylab(expression(alpha[k]^(2)*(t))) + xlab("t") + xlim(c(0, 1)) + ylim(c(0, 1))

piecewiseLinearHist22 <- ggplot() + 
		geom_line(aes(x, value, group = iteration), d2, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(t, value), m2, col = "red", size = 1.0) + ggtitle("drawing") + 
		theme_bw() + ylab(expression(alpha[k]^(2)*(t))) + xlab("t") + xlim(c(0, 1)) + ylim(c(0, 1))

piecewiseLinearHist23 <- ggplot() + 
		geom_line(aes(x, value, group = iteration), d3, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(t, value), m3, col = "red", size = 1.0) + ggtitle("winning") + 
		theme_bw() +ylab(expression(alpha[k]^(2)*(t))) + xlab("t") + xlim(c(0, 1)) + ylim(c(0, 1))

plotList[[7]] <- piecewiseLinearHist21
plotList[[8]] <- piecewiseLinearHist22
plotList[[9]] <- piecewiseLinearHist23

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/alphaBayesian3.pdf", height = 9)
do.call(grid.arrange, c(plotList, ncol = 3))
dev.off()

##########
##########

posteriorSamples <- read.csv("/home/jeff/workspace/C-Football/csv/history.csv", header = FALSE)
events <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/data/2011-2012-data.csv", stringsAsFactors = FALSE)
manCityQPR <- events[events$matchID == 1117411, ] #man city - qpr

names <- c("h", "a", "c_-1", "c_0", "c_1", "d_-1", "d_0", 
		"d_1", "rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Bromwich Albion", "Wigan Athletic", "Wolverhampton Wanderers")

posteriorMean <- colMeans(posteriorSamples)

h <- posteriorMean[1]
a <- posteriorMean[2]
c <- posteriorMean[3:5]
d <- posteriorMean[6:8]
rho <- c(0.0, posteriorMean[9:10])
homeR <- posteriorMean[names == "Manchester City"]
awayR <- posteriorMean[names == "Queens Park Rangers"]

getRhoIdx <- function(t) {
	if (t <= 44.0 / 90.0) {
		return(1)
	} else if (t <= 45.0 / 90) {
		return(2);
	} else if (t <= 89.0 / 90.0) {
		return(1);
	} else {
		return(3);
	}
}

rate <- function(t, homeScore, awayScore, isHome) {
	
	if (homeScore > awayScore) {
		homeIdx <- 3
	} else if (homeScore == awayScore) {
		homeIdx <- 2
	} else {
		homeIdx <- 1
	}
	awayIdx <- 4 - homeIdx
	
	if (isHome) {
		return(exp(h + rho[getRhoIdx(t)] + (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t) * homeR - 
								(1.0 - (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t)) * awayR))
	} else {
		return(exp(a + rho[getRhoIdx(t)] + (c[homeIdx] + (d[homeIdx] - c[homeIdx]) * t) * homeR - 
								(1.0 - (c[awayIdx] + (d[awayIdx] - c[awayIdx]) * t)) * awayR))
	}
}

n <- 100L
m <- nrow(manCityQPR)

l <- vector("list", m)
lam <- double(n)
mu <- double(n)

for (i in seq_len(m)) {
	
	start <- manCityQPR$start[i] / 90.0
	end <- manCityQPR$end[i] / 90.0
	homeScore <- manCityQPR$homeScore[i]
	awayScore <- manCityQPR$awayScore[i]
	event <- manCityQPR$event[i]
	
	s <- seq(start, end, length = 100)
	
	for (j in seq_len(n)) {
		lam[j] <- rate(s[j], homeScore, awayScore, TRUE)
		mu[j] <- rate(s[j], homeScore, awayScore, FALSE)
	}
	
	l[[i]] <- data.frame(lam, mu, t = s)
}

inplayData <- do.call(rbind, l)
goalTimes <- manCityQPR[manCityQPR$event %in% 1:2, ]

arrowData <- data.frame(
		x0 = c(34, 43, 61, 85, 85), 
		y0 = c(13, 13, 13, 14.5, 11.5), 
		x1 = c(38.5, 47.5, 65.5, 89.5, 89.5), 
		y1 = c(11.5, 11.5, 11.5, 13, 10),
		score = c("1-0", "1-1", "1-2", "3-2", "2-2"))

inplay <- ggplot(inplayData) + geom_vline(aes(xintercept = end), goalTimes, size = 1.0, col = "grey") +
		geom_line(aes(x = t * 90.0, y = lam), col = "red", size = 1.0) + 
		geom_line(aes(x = t * 90.0, y = mu), col = "blue", size = 1.0, linetype = 5) + 
		scale_x_continuous(breaks = seq(0, 90, by = 5)) +
		labs(x = "time (minutes)", y = "rate") +
		geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), arrowData, arrow = arrow()) + 
		geom_text(aes(x = x0, y = y0 + 0.5, label = score), arrowData) + 
		theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/inplay.pdf", inplay)

##########
##########

stats <- read.csv("/home/jeff/workspace/C-Football/csv/chi-sq.csv", header = FALSE)
chiSqStat <- 16.4905
df <- 17
pVals <- pchisq(stats[ , 2], stats[ , 1], lower.tail = FALSE)

chisqPlot <- ggplot() + geom_histogram(aes(stats[, 2], ..density..), binwidth = 1) + 
		geom_line(aes(x = rep(chiSqStat, 2)), y = c(0.0, 10), size = 1, col = "red") + 
		xlab(expression(paste(chi^2, " statistic"))) + 
		theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/chiSqStats.pdf", chisqPlot)

##########
##########

DrProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsDR.csv", header = FALSE)
BtProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsBT.csv", header = FALSE)
BBtProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsBBT.csv", header = FALSE)
MProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsM.csv", header = FALSE)

DrProbs <- tail(DrProbs, 330)
BtProbs <- tail(BtProbs, 330)
MProbs <- tail(MProbs, 330)
BBtProbs <- tail(BBtProbs, 330)

colnames(DrProbs)[1:3] <- paste0("DR_", c("home", "draw", "away"))
colnames(BtProbs)[1:3] <- paste0("BT_", c("home", "draw", "away"))
colnames(BBtProbs)[1:3] <- paste0("BBT_", c("home", "draw", "away"))
colnames(MProbs)<- c(paste0("M_", c("home", "draw", "away")), "home", "away")

probs <- cbind(DrProbs[1:3], BtProbs[1:3], BBtProbs[1:3], MProbs)
odds <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/E0-2011-2012.csv", 
		stringsAsFactors = FALSE)
probs$order <- seq_len(nrow(probs))

odds$home <- unclass(as.factor(odds$HomeTeam)) - 1L
odds$away <- unclass(as.factor(odds$AwayTeam)) - 1L

data <- sqldf("
				SELECT p.*, o.FTHG AS homeScore, o.FTAG AS awayScore
				FROM probs p
				JOIN odds o 
				ON o.home = p.home AND o.away = p.away 
				")

observedHome <- with(data, homeScore > awayScore)
observedDraw <- with(data, homeScore == awayScore)
observedAway <- with(data, awayScore > homeScore)

expectedMHome <- data$M_home
expectedMDraw <- data$M_draw
expectedMAway <- data$M_away

expectedBTHome <- data$BT_home
expectedBTDraw <- data$BT_draw
expectedBTAway <- data$BT_away

expectedBBTHome <- data$BBT_home
expectedBBTDraw <- data$BBT_draw
expectedBBTAway <- data$BBT_away

expectedDRHome <- data$DR_home
expectedDRDraw <- data$DR_draw
expectedDRAway <- data$DR_away

cutPointsHome <- c(0.0,0.3,0.4,0.5,1.0)
cutPointsDraw <- c(0.0,0.25,0.28,1.0)
cutPointsAway <- c(0.0,0.2,0.3,0.4,1.0)

f1 <- function(x) {
	mu <- sum(x)
	sigma2 <- sum(x * (1.0 - x))
	return(c(mean = mu, sd = sqrt(sigma2)))
}

#home
binnedHomeExpM <- aggregate(expectedMHome ~ cut(expectedMHome, breaks = cutPointsHome), , f1)
binnedHomeObsM <- aggregate(observedHome ~ cut(expectedMHome, breaks = cutPointsHome), , sum)

binnedHomeExpBT <- aggregate(expectedBTHome ~ cut(expectedBTHome, breaks = cutPointsHome), , f1)
binnedHomeObsBT <- aggregate(observedHome ~ cut(expectedBTHome, breaks = cutPointsHome), , sum)

binnedHomeExpBBT <- aggregate(expectedBBTHome ~ cut(expectedBBTHome, breaks = cutPointsHome), , f1)
binnedHomeObsBBT <- aggregate(observedHome ~ cut(expectedBBTHome, breaks = cutPointsHome), , sum)

binnedHomeExpDR <- aggregate(expectedDRHome ~ cut(expectedDRHome, breaks = cutPointsHome), , f1)
binnedHomeObsDR <- aggregate(observedHome ~ cut(expectedDRHome, breaks = cutPointsHome), , sum)

ggBinnedHomeExpM <- data.frame(binnedHomeExpM$expectedMHome, x = seq(0.5, 12.5, by = 4), group = "M")
ggBinnedHomeExpDR <- data.frame(binnedHomeExpDR$expectedDRHome, x = seq(0.5, 12.5, by = 4) + 1, group = "DR")
ggBinnedHomeExpBT <- data.frame(binnedHomeExpBT$expectedBTHome, x = seq(0.5, 12.5, by = 4) + 2, group = "BT")
ggBinnedHomeExpBBT <- data.frame(binnedHomeExpBBT$expectedBBTHome, x = seq(0.5, 12.5, by = 4) + 3, group = "BBT")

ggBinnedHomeObsM <- data.frame(obs = binnedHomeObsM$observedHome, x = seq(0.5, 12.5, by = 4), group = "M")
ggBinnedHomeObsDR <- data.frame(obs = binnedHomeObsDR$observedHome, x = seq(0.5, 12.5, by = 4) + 1, group = "DR")
ggBinnedHomeObsBT <- data.frame(obs = binnedHomeObsBT$observedHome, x = seq(0.5, 12.5, by = 4) + 2, group = "BT")
ggBinnedHomeObsBBT <- data.frame(obs = binnedHomeObsBBT$observedHome, x = seq(0.5, 12.5, by = 4) + 3, group = "BBT")

ggHomeExp <- rbind(ggBinnedHomeExpM, ggBinnedHomeExpDR, ggBinnedHomeExpBT, ggBinnedHomeExpBBT)
ggHomeObs <- rbind(ggBinnedHomeObsM, ggBinnedHomeObsDR, ggBinnedHomeObsBT, ggBinnedHomeObsBBT)

home <- ggplot() + theme_bw() + labs(x = "home win probability interval (I)", y = "frequency") + 
		geom_errorbar(aes(x, ymin = qnorm(0.025, mean, sd), ymax = qnorm(0.975, mean, sd)), ggHomeExp, width = 0.3) + 
		geom_point(aes(x, y = obs, col = group), ggHomeObs, size = 3) + 
		scale_color_manual(values = c("blue", "red", "green", "black"), guide = "none") + 
		scale_y_continuous(breaks = seq(0, 100, by = 20)) + 
		scale_x_continuous(breaks = seq(2, 14, by = 4), labels = gsub(",", ", ", binnedHomeExpM[ , 1])) + 
		geom_vline(aes(xintercept = seq(0, 16, by = 4)), linetype = 3) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#draw
binnedDrawExpM <- aggregate(expectedMDraw ~ cut(expectedMDraw, breaks = cutPointsDraw), , f1)
binnedDrawObsM <- aggregate(observedDraw ~ cut(expectedMDraw, breaks = cutPointsDraw), , sum)

binnedDrawExpBT <- aggregate(expectedBTDraw ~ cut(expectedBTDraw, breaks = cutPointsDraw), , f1)
binnedDrawObsBT <- aggregate(observedDraw ~ cut(expectedBTDraw, breaks = cutPointsDraw), , sum)

binnedDrawExpBBT <- aggregate(expectedBBTDraw ~ cut(expectedBBTDraw, breaks = cutPointsDraw), , f1)
binnedDrawObsBBT <- aggregate(observedDraw ~ cut(expectedBBTDraw, breaks = cutPointsDraw), , sum)

binnedDrawExpDR <- aggregate(expectedDRDraw ~ cut(expectedDRDraw, breaks = cutPointsDraw), , f1)
binnedDrawObsDR <- aggregate(observedDraw ~ cut(expectedDRDraw, breaks = cutPointsDraw), , sum)

ggBinnedDrawExpM <- data.frame(binnedDrawExpM$expectedMDraw, x = seq(0.5, 8.5, by = 4), group = "M")
ggBinnedDrawExpDR <- data.frame(binnedDrawExpDR$expectedDRDraw, x = seq(0.5, 8.5, by = 4) + 1, group = "DR")
ggBinnedDrawExpBT <- data.frame(binnedDrawExpBT$expectedBTDraw, x = seq(0.5, 8.5, by = 4) + 2, group = "BT")
ggBinnedDrawExpBBT <- data.frame(binnedDrawExpBBT$expectedBBTDraw, x = seq(0.5, 8.5, by = 4) + 3, group = "BBT")

ggBinnedDrawObsM <- data.frame(obs = binnedDrawObsM$observedDraw, x = seq(0.5, 8.5, by = 4), group = "M")
ggBinnedDrawObsDR <- data.frame(obs = binnedDrawObsDR$observedDraw, x = seq(0.5, 8.5, by = 4) + 1, group = "DR")
ggBinnedDrawObsBT <- data.frame(obs = binnedDrawObsBT$observedDraw, x = seq(0.5, 8.5, by = 4) + 2, group = "BT")
ggBinnedDrawObsBBT <- data.frame(obs = binnedDrawObsBBT$observedDraw, x = seq(0.5, 8.5, by = 4) + 3, group = "BBT")

ggDrawExp <- rbind(ggBinnedDrawExpM, ggBinnedDrawExpDR, ggBinnedDrawExpBT, ggBinnedDrawExpBBT)
ggDrawObs <- rbind(ggBinnedDrawObsM, ggBinnedDrawObsDR, ggBinnedDrawObsBT, ggBinnedDrawObsBBT)

draw <- ggplot() + theme_bw() + labs(x = "draw win probability interval (I)", y = "frequency") + 
		geom_errorbar(aes(x, ymin = qnorm(0.025, mean, sd), ymax = qnorm(0.975, mean, sd)), ggDrawExp, width = 0.3) + 
		geom_point(aes(x, y = obs, col = group), ggDrawObs, size = 3) + 
		scale_color_manual(values = c("blue", "red", "green", "black"), guide = "none") + 
		scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) + 
		scale_x_continuous(breaks = seq(2, 10, by = 4), labels = gsub(",", ", ", binnedDrawExpM[ , 1])) + 
		geom_vline(aes(xintercept = seq(0, 12, by = 4)), linetype = 3) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#away
binnedAwayExpM <- aggregate(expectedMAway ~ cut(expectedMAway, breaks = cutPointsAway), , f1)
binnedAwayObsM <- aggregate(observedAway ~ cut(expectedMAway, breaks = cutPointsAway), , sum)

binnedAwayExpBT <- aggregate(expectedBTAway ~ cut(expectedBTAway, breaks = cutPointsAway), , f1)
binnedAwayObsBT <- aggregate(observedAway ~ cut(expectedBTAway, breaks = cutPointsAway), , sum)

binnedAwayExpBBT <- aggregate(expectedBBTAway ~ cut(expectedBBTAway, breaks = cutPointsAway), , f1)
binnedAwayObsBBT <- aggregate(observedAway ~ cut(expectedBBTAway, breaks = cutPointsAway), , sum)

binnedAwayExpDR <- aggregate(expectedDRAway ~ cut(expectedDRAway, breaks = cutPointsAway), , f1)
binnedAwayObsDR <- aggregate(observedAway ~ cut(expectedDRAway, breaks = cutPointsAway), , sum)

ggBinnedAwayExpM <- data.frame(binnedAwayExpM$expectedMAway, x = seq(0.5, 12.5, by = 4), group = "M")
ggBinnedAwayExpDR <- data.frame(binnedAwayExpDR$expectedDRAway, x = seq(0.5, 12.5, by = 4) + 1, group = "DR")
ggBinnedAwayExpBT <- data.frame(binnedAwayExpBT$expectedBTAway, x = seq(0.5, 12.5, by = 4) + 2, group = "BT")
ggBinnedAwayExpBBT <- data.frame(binnedAwayExpBBT$expectedBBTAway, x = seq(0.5, 12.5, by = 4) + 3, group = "BBT")

ggBinnedAwayObsM <- data.frame(obs = binnedAwayObsM$observedAway, x = seq(0.5, 12.5, by = 4), group = "M")
ggBinnedAwayObsDR <- data.frame(obs = binnedAwayObsDR$observedAway, x = seq(0.5, 12.5, by = 4) + 1, group = "DR")
ggBinnedAwayObsBT <- data.frame(obs = binnedAwayObsBT$observedAway, x = seq(0.5, 12.5, by = 4) + 2, group = "BT")
ggBinnedAwayObsBBT <- data.frame(obs = binnedAwayObsBBT$observedAway, x = seq(0.5, 12.5, by = 4) + 3, group = "BBT")

ggAwayExp <- rbind(ggBinnedAwayExpM, ggBinnedAwayExpDR, ggBinnedAwayExpBT, ggBinnedAwayExpBBT)
ggAwayObs <- rbind(ggBinnedAwayObsM, ggBinnedAwayObsDR, ggBinnedAwayObsBT, ggBinnedAwayObsBBT)

away <- ggplot() + theme_bw() + labs(x = "away win probability interval (I)", y = "frequency") + 
		geom_errorbar(aes(x, ymin = qnorm(0.025, mean, sd), ymax = qnorm(0.975, mean, sd)), ggAwayExp, width = 0.3) + 
		geom_point(aes(x, y = obs, col = group), ggAwayObs, size = 3) + 
		scale_color_manual(values = c("blue", "red", "green", "black"), guide = "none") + 
		scale_y_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100)) + 
		scale_x_continuous(breaks = seq(2, 14, by = 4), labels = gsub(",", ", ", binnedAwayExpM[ , 1])) + 
		geom_vline(aes(xintercept = seq(0, 16, by = 4)), linetype = 3) + 
		theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/bins.pdf", height = 9)
grid.arrange(home, draw, away, ncol = 1)
dev.off()

##########
##########

DrProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsDR.csv", header = FALSE)
BtProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsBT.csv", header = FALSE)
BBtProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsBBT.csv", header = FALSE)
MProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsM.csv", header = FALSE)

DrProbs <- tail(DrProbs, 330)
BtProbs <- tail(BtProbs, 330)
MProbs <- tail(MProbs, 330)
BBtProbs <- tail(BBtProbs, 330)

colnames(DrProbs)[1:3] <- paste0("DR_", c("home", "draw", "away"))
colnames(BtProbs)[1:3] <- paste0("BT_", c("home", "draw", "away"))
colnames(BBtProbs)[1:3] <- paste0("BBT_", c("home", "draw", "away"))
colnames(MProbs)<- c(paste0("M_", c("home", "draw", "away")), "home", "away")

probs <- cbind(DrProbs[1:3], BtProbs[1:3], BBtProbs[1:3], MProbs)
odds <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/E0-2011-2012.csv", 
		stringsAsFactors = FALSE)
probs$order <- seq_len(nrow(probs))

odds$home <- unclass(as.factor(odds$HomeTeam)) - 1L
odds$away <- unclass(as.factor(odds$AwayTeam)) - 1L

data <- sqldf("
				SELECT p.*, o.FTHG AS homeScore, o.FTAG AS awayScore, o.B365H, o.B365D, o.B365A
				FROM probs p
				JOIN odds o 
				ON o.home = p.home AND o.away = p.away 
				")

r <- seq(from = 0.0, to = 2.5, by = 0.01)
nBets <- integer(length(r))
profit <- double(length(r))
tiny <- .Machine$double.eps # cause DR is stupid

profitFunc <- function(home, draw, away) {
	
	for (i in seq_along(r)) {
		
		betHome <- home * data$B365H > r[i] - tiny
		betDraw <- draw *data$B365D > r[i] - tiny
		betAway <- away * data$B365A > r[i] - tiny
		
		profitHomeM <- ifelse(betHome, (data$homeScore > data$awayScore) * data$B365H - 1.0, 0.0)
		profitDrawM <- ifelse(betDraw, (data$homeScore == data$awayScore) * data$B365D - 1.0, 0.0)
		profitAwayM <- ifelse(betAway, (data$homeScore < data$awayScore) * data$B365A - 1.0, 0.0)
		
		profit[i] <- sum(profitHomeM + profitDrawM + profitAwayM, na.rm = TRUE)
	}
	return(profit)
}

plotData <- melt(
		data.frame(r, 
				M = profitFunc(data$M_home, data$M_draw, data$M_away),
				DR = profitFunc(data$DR_home, data$DR_draw, data$DR_away),
				BT = profitFunc(data$BT_home, data$BT_draw, data$BT_away),
				BBT = profitFunc(data$BBT_home, data$BBT_draw, data$BBT_away)
		), id.vars = "r")

pl <- ggplot(plotData) + geom_line(aes(r, value, col = variable), size = 0.75) + 
		scale_color_manual(values = c("blue", "red", "green", "black"), guide = "none") + 
		scale_y_continuous(breaks = seq(-50, 70, by = 10)) + 
		geom_hline(aes(yintercept = 0), lty = 3) + ylab(expression(O["P,r"])) + theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/ggFigures/profit.pdf", pl)

# a check (again!)
sum(
(data$homeScore > data$awayScore) * log(data$BT_home) +
(data$homeScore == data$awayScore) * log(data$BT_draw) +
(data$homeScore < data$awayScore) * log(data$BT_away)
)

##########
##########

kelly <- function(home, draw, away, bank) {
	
	fHome <- (home * data$B365H - 1.0) / (data$B365H - 1.0)
	fdraw <- (home * data$B365D - 1.0) / (data$B365D - 1.0)
	faway <- (home * data$B365A - 1.0) / (data$B365A - 1.0)

	for (i in seq_along(home)) {
		bank <- bank + (fHome > 0.0) * fHome * bank * ((data$homeScore > data$awayScore) * data$B365H - 1.0)
		bank <- bank + (fDraw > 0.0) * fDraw * bank * ((data$homeScore == data$awayScore) * data$B365D - 1.0)
		bank <- bank + (fAway > 0.0) * fAway * bank * ((data$homeScore < data$awayScore) * data$B365A - 1.0)
	}
	
	return(bank)
}

