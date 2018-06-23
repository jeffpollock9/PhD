library(ggplot2)
library(gridExtra)

x <- 1:20
y <- 20:1
y[18:20] <- y[18:20]*0.2 - 2.0
y[1] <- y[1]*1.5
y[2:4] <- y[2:4]*1.25
y[5] <- y[5]*1.15
y <- y + 1.8

knots <- c(1, 2, 4, 5, 6, 17, 18, 20)
plotData <- data.frame(x, y, knot = 1:20 %in% knots)
arrowData <- data.frame(
		x0 = c(x[19], x[19], x[19], x[6], x[4], x[4], x[4], x[4], x[3]),
		y0 = c(10, 10, 10, 21, 29, 29, 29, 29, 34),
		x1 = c(x[20], x[19], x[18], x[5]+0.1, x[1]+0.2, x[2]+0.1, x[3]+0.1, x[4]+0.1, x[1]+0.1),
		y1 = c(y[20]+0.4, y[19]+0.4, y[18]+0.4, y[5]+0.1, y[1], y[2]+0.1, y[3]+0.1, y[4]+0.1, y[1]+0.1)
)
textData <- data.frame(
		x = c(x[8]+0.25, x[7], x[19], x[5]+0.75),
		y = c(22, 29.5, 10.9, 35),
		text = c("Europa League", "Champions League", "Relegation", "League Champion")
)

leaguePlot <- ggplot(plotData) + geom_line(aes(x, y)) + geom_point(aes(x, y, fill = knot), size = 3, shape = 21) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), arrowData, arrow = arrow(length = unit(0.25, "cm"))) + 
		geom_text(aes(x, y, label = text), textData) + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "league position", y = "utility")

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/utility_league.pdf", leaguePlot)

##########
##########

u <- c(10.0, 8.0, 7.75, 7.5, 7.0, 6.5, 6.27273, 6.04545, 5.81818, 5.59091, 5.36364, 5.13636, 
		4.90909, 4.68182, 4.45455, 4.22727, 4.0, 1.0, 0.5)
alpha <- 2.5
beta <- alpha / u
cols <- c("red", "green", "blue", "pink", "black", "grey", "orange")

utilityPrior <- ggplot(NULL, aes(x, colour = g)) + 
		stat_function(data = data.frame(x = 0:20, g = factor(1)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[1])) + 
		stat_function(data = data.frame(x = 0:20, g = factor(2)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[2])) + 
		stat_function(data = data.frame(x = 0:20, g = factor(3)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[4])) + 
		stat_function(data = data.frame(x = 0:20, g = factor(4)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[5])) + 
		stat_function(data = data.frame(x = 0:20, g = factor(5)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[6])) + 
		stat_function(data = data.frame(x = 0:20, g = factor(6)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[17])) + 
		stat_function(data = data.frame(x = 0:20, g = factor(7)), fun = dgamma, group = 1, args = list(shape = alpha, rate = beta[18])) + 
		scale_colour_manual(values = cols, guide = "none") + theme_bw() + labs(x = expression(U[p]), y = expression(p(U[p])))

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/utility_prior.pdf", utilityPrior)

##########
##########
#temp added /new into file path for updated code

historyX <- c(0, 1, 2, 3)

piecewiseLinear <- function(p, x, y) {
	
	stopifnot(length(x) == length(y))
	
	n <- length(x)
	if (n == 1L) {
		return(y)
	}
	
	for (i in 2:n) {
		if (p <= x[i] | i == n) {
			gradient <- (y[i] - y[i-1]) / (x[i] - x[i-1])
			intercept <- y[i] - gradient * x[i]
			return(intercept + p * gradient)
		}
	}
}

#when estimating parameters, used matchesPlayed, not week
beta <- function(w) {
	matchesPlayed <- (w - 1L)
	return(plogis(A + B * matchesPlayed))
}

names0 <- c("h", "a", "c_-1", "c_0", "c_1", "d_-1", "d_0", "d_1", "rho1", "rho2",
		"Arsenal", "Aston Villa", "Birmingham City", "Blackburn Rovers", 
		"Blackpool", "Bolton Wanderers", "Burnley", "Chelsea", "Derby County", 
		"Everton", "Fulham", "Hull City", "Liverpool", "Manchester City", 
		"Manchester United", "Middlesbrough", "Newcastle United", "Norwich City", 
		"Portsmouth", "Queens Park Rangers", "Reading", "Stoke City", 
		"Sunderland", "Swansea City", "Tottenham Hotspur", "West Bromwich Albion", 
		"West Ham United", "Wigan Athletic", "Wolverhampton Wanderers", "A", "B")

for (history in historyX) {
	
	if (history == 0) {
		x <- c(1.0, 17.0, 18.0, 20.0)
	} else if (history == 1) {
		x <- c(1.0, 5.0, 17.0, 18.0, 20.0)
	} else if (history == 2) {
		x <- c(1.0, 2.0, 5.0, 17.0, 18.0, 20.0)
	} else if (history == 3) {
		x <- c(1.0, 2.0, 4.0, 5.0, 6.0, 17.0, 18.0, 20.0)
	}
	
	thetaHistory <- read.csv(sprintf("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/new/MCMC%s/thetaHistory%s.txt", 
					history, history), header = FALSE, colClasses = "numeric")
	colnames(thetaHistory) <- names0
	colMeans(thetaHistory)
	idx <- (length(names0) + 1):(length(names0) + length(x) - 1)
	colnames(thetaHistory)[idx] <- paste0("knot", x[-length(x)])
	
	f <- numeric(20)
	nIterations <- nrow(thetaHistory)
	burn <- 1L
	thin <- 1L
	I <- (burn:nIterations)[(burn:nIterations) %% thin == 0L]
	nSamples <- 1500
	
	yHistory <- vector("list", nIterations)
	for (i in 1:nIterations) {
		yHistory[[i]] <- c(as.numeric(thetaHistory[i, idx]), 0.0)
	}
	
	samples <- sample.int(nIterations, nSamples, FALSE)
	utilityValues <- vector("list", nSamples)
	for (i in seq_along(samples)) {
		for (p in 1:20) {
			f[p] <- piecewiseLinear(p, x, yHistory[[samples[i]]])
		}
		utilityValues[[i]] <- data.frame(p = 1:20, utility = f, sample = i)	
	}
	
	utilityValueData <- do.call(rbind, utilityValues)
	
	firstKnotIdx <- length(names0) + 1
	lastKnotIdx <- ncol(thetaHistory)
	knots <- thetaHistory[ , firstKnotIdx:lastKnotIdx]
	meanData <- data.frame(p = x, utility = c(colMeans(knots), 0))		

	if (history == 3) {
		cat("MCMC3 P(U(6) > U(5)) = ", mean(knots$knot6 > knots$knot5), "\n")
	}
	
	pl1 <- ggplot() + theme_bw() + labs(x = "p", y = "U(p)") + 
			geom_line(aes(p, utility, group = sample), utilityValueData, col = rgb(0, 0, 1, 0.1)) + 
			geom_line(aes(p, utility), meanData, col = "red", size = 1) +
			scale_x_discrete(breaks = 1:20)
	
	samples <- sample.int(nIterations, nSamples, FALSE)
	betaValues <- vector("list", nSamples)
	w <- seq(1, 38, by = 0.1)
	for (i in seq_along(samples)) {
		A <- thetaHistory[samples[i], "A"]
		B <- thetaHistory[samples[i], "B"]
		betaValues[[i]] <- data.frame(w, beta = beta(w), sample = i)
	}
	
	betaData <- do.call(rbind, betaValues)
	
	A <- mean(thetaHistory[ , "A"])
	B <- mean(thetaHistory[ , "B"])
	meanData <- data.frame(w, beta = beta(w))
	
	pl2 <- ggplot() + theme_bw() + labs(x = "w", y = expression(beta(w))) + 
			geom_line(aes(w, beta, group = sample), betaData, col = rgb(0, 0, 1, 0.1)) + 
			geom_line(aes(w, beta), meanData, col = "red", size = 1) + ylim(c(0, 1)) + 
			scale_x_discrete(breaks = seq(0, 38, by = 2))
	
	pdf(sprintf("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/new/utility_%s.pdf", history), height = 9)
	grid.arrange(pl1, pl2, ncol = 1)
	dev.off()
	gc()
}

##########
##########

knotsM <- c(TRUE, rep(FALSE, 8), FALSE, rep(FALSE, 9), TRUE)
dataM <- data.frame(p = 1:20, utility = 20:1, knot = knotsM)

utilityPrime <- numeric(20)
for (i in 1:20) {
	utilityPrime[i] <- piecewiseLinear(i, c(1, 10, 20), c(20, 15, 1))	
}
knotsPrime <- c(TRUE, rep(FALSE, 8), TRUE, rep(FALSE, 9), TRUE)
dataMPrime <- data.frame(p = 1:20, utility = utilityPrime, knot = knotsPrime)

pl1 <- ggplot(dataM) + geom_line(aes(p, utility)) + geom_point(aes(p, utility, fill = factor(knot)), size = 3, shape = 21) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "p", y = "U(p)")

y <- seq(0, 20, by = 0.1)
x <- dnorm(y, mean = 11, 3)*20 + 10
dataDensity <- data.frame(x, y)

pl2 <- ggplot(dataM) + geom_line(aes(p, utility), col = "grey") + 
		geom_vline(xintercept = 10, linetype = 3) + 
		geom_point(aes(p, utility, fill = knot), col = "grey", size = 3, shape = 21) + 
		geom_line(aes(p, utility), dataMPrime) +
		geom_point(aes(p, utility, fill = knot), dataMPrime, size = 3, shape = 21) + 
		geom_segment(aes(x = 10, y = 11.5, xend = 10, yend = 14.3), size = 0.5, arrow = arrow(length = unit(0.25, "cm"))) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "p", y = "U(p)") + 
		geom_path(aes(x, y), dataDensity, linetype = 3)

pl3 <- ggplot(dataMPrime) + geom_line(aes(p, utility)) + geom_point(aes(p, utility, fill = factor(knot)), size = 3, shape = 21) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "p", y = "U(p)")

pdf("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/knotBirthExample.pdf", height = 9)
grid.arrange(pl1, pl2, pl3, ncol = 1)
dev.off()

##########
##########

x <- c(1, 3, 5, 10, 20)
y <- c(20, 21, 16, 7, 0)		
knotsM <- 1:20 %in% x

utility <- numeric(20)
for (i in 1:20) {
	utility[i] <- piecewiseLinear(i, x, y)	
}
dataM <- data.frame(p = 1:20, utility = utility, knot = knotsM)

x <- c(1, 5, 10, 20)
y <- c(20, 16, 7, 0)		
knotsPrime <- 1:20 %in% x
utilityPrime <- numeric(20)
for (i in 1:20) {
	utilityPrime[i] <- piecewiseLinear(i, x, y)	
}
dataMPrime <- data.frame(p = 1:20, utility = utilityPrime, knot = knotsPrime)

pl1 <- ggplot(dataM) + geom_line(aes(p, utility)) + geom_point(aes(p, utility, fill = factor(knot)), size = 3, shape = 21) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "p", y = "U(p)")

pl2 <- ggplot(dataM) + geom_line(aes(p, utility), col = "grey") + 
		geom_point(aes(p, utility, fill = knot), col = "grey", size = 3, shape = 21) + 
		geom_line(aes(p, utility), dataMPrime) +
		geom_point(aes(p, utility, fill = knot), dataMPrime, size = 3, shape = 21) + 
		geom_segment(aes(x = 3, y = 20.8, xend = 3, yend = 18.25), size = 0.5, arrow = arrow(length = unit(0.25, "cm"))) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "p", y = "U(p)")

pl3 <- ggplot(dataMPrime) + geom_line(aes(p, utility)) + geom_point(aes(p, utility, fill = factor(knot)), size = 3, shape = 21) + 
		scale_fill_manual(values = c("white", "black"), guide = "none") + 
		scale_x_discrete(breaks = 1:20) + 
		scale_y_continuous(breaks = NULL) + 
		theme_bw() + labs(x = "p", y = "U(p)")

pdf("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/knotDeathExample.pdf", height = 9)
grid.arrange(pl1, pl2, pl3, ncol = 1)
dev.off()

##########
##########

kHistory <- as.integer(readLines(sprintf("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/RJMCMC/kHistory.txt", history)))
thetaHistory <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/RJMCMC/thetaHistory.txt", 
		header = FALSE, colClasses = "numeric")
colnames(thetaHistory) <- names0

xHistory <- as.list(readLines("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/RJMCMC/knotXHistory.txt", ))
yHistory <- as.list(readLines("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/RJMCMC/knotYHistory.txt", ))

xHistory <- lapply(xHistory, function(x) as.double(unlist(strsplit(x, ","))))
yHistory <- lapply(yHistory, function(x) as.double(unlist(strsplit(x, ","))))

f <- numeric(20)
nIterations <- nrow(thetaHistory)
burn <- 1L
thin <- 1L
nSamples <- 1500
I <- (burn:nIterations)[(burn:nIterations) %% thin == 0L]

samples <- sample.int(nIterations, nSamples, FALSE)
utilityValues <- vector("list", nSamples)
for (i in seq_along(samples)) {
	for (p in 1:20) {
		f[p] <- piecewiseLinear(p, xHistory[[samples[i]]], yHistory[[samples[i]]])
	}
	utilityValues[[i]] <- data.frame(p = 1:20, utility = f, sample = i)	
}

knots <- matrix(NA_real_, length(I), 20)
colnames(knots) <- paste0("knot", 1:20)

for (i in I) {
	for (p in seq_len(20)) {
		f[p] <- piecewiseLinear(p, xHistory[[i]], yHistory[[i]])
	}
	knots[i, ] <- f
}

utilityValueData <- do.call(rbind, utilityValues)
meanData <- data.frame(p = 1:20, utility = colMeans(knots))

pl1 <- ggplot() + theme_bw() + labs(x = "p", y = "U(p)") + 
		geom_line(aes(p, utility, group = sample), utilityValueData, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(p, utility), meanData, col = "red", size = 1) +
		scale_x_discrete(breaks = 1:20)

samples <- sample.int(nIterations, nSamples, FALSE)
betaValues <- vector("list", nSamples)
w <- seq(1, 38, by = 0.1)
for (i in seq_along(samples)) {
	A <- thetaHistory[samples[i], "A"]
	B <- thetaHistory[samples[i], "B"]
	betaValues[[i]] <- data.frame(w, beta = beta(w), sample = i)
}

betaData <- do.call(rbind, betaValues)

A <- mean(thetaHistory[ , "A"])
B <- mean(thetaHistory[ , "B"])
meanData <- data.frame(w, beta = beta(w))

pl2 <- ggplot() + theme_bw() + labs(x = "w", y = expression(beta(w))) + 
		geom_line(aes(w, beta, group = sample), betaData, col = rgb(0, 0, 1, 0.1)) + 
		geom_line(aes(w, beta), meanData, col = "red", size = 1) + ylim(c(0, 1)) + 
		scale_x_discrete(breaks = seq(0, 38, by = 2))

pdf("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/utility_RJ.pdf", height = 9)
grid.arrange(pl1, pl2, ncol = 1)
dev.off()

##########
##########

prior <- function(p) {
	dbinom(p, 20, priorProb) / (1 - pbinom(1, 20, priorProb))
}

priorProb <- 0.05

posteriorTable <- table(factor(kHistory[I])) / nIterations
priorValues <- prior(sort(unique(kHistory)))

knotNumberData <- rbind(
		data.frame(
				number = sort(unique(kHistory)), 
				prob = as.vector(posteriorTable),
				type = "posterior"),
		data.frame(	
				number = sort(unique(kHistory)),
				prob = priorValues,
				type = "prior")
)

pl1 <- ggplot(knotNumberData) + geom_bar(aes(number, prob, fill = type), stat = "identity", position = "dodge") + 
		scale_x_discrete() + theme_bw() + scale_fill_manual(values = c("blue", "red"), guide = "none") + 
		labs(x = "K", y = "probability")

knotPositionData <- data.frame(
		p = 1:20,
		prob = as.vector(table(factor(unlist(xHistory[I]), levels = 1:20)) / nIterations)
)

pl2 <- ggplot(knotPositionData) + geom_bar(aes(p, prob), stat = "identity", fill = "blue") + 
		scale_x_discrete(breaks = 1:20) + theme_bw() + labs(x = "p", y = "probability")

pdf("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/number_of_knots_RJ.pdf", height = 9)
grid.arrange(pl1, pl2, ncol = 1)
dev.off()

cat("RJMCMC P(U(6) > U(5)) = ", mean(knots[ , 6] > knots[ , 5]), "\n")

##########
##########

events <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/data/2011-2012-data.csv", stringsAsFactors = FALSE)
matches <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/data/2011-2012-matches.csv", stringsAsFactors = FALSE)

manCityQPR <- events[events$matchID == 1117411, ] 
sunderlandManUtd <- events[events$matchID == 1117414, ] 
StokeBolton <- events[events$matchID == 1117413, ] 
norwichAstonVilla <- events[events$matchID == 1117412, ]
WiganWolves <- events[events$matchID == 1117418, ]

rates <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/rates.csv", header = FALSE)
alpha <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Utility/R/alpha.csv", header = FALSE)

colnames(rates) <- c("t", "home", "away")
colnames(alpha) <- c("t", "home", "away")

x <- rnorm(100)

acf(x, plot = FALSE)

#calc lag 10 acf
a <- 0.0
for (i in 1:(100 - 10)) {
	a <- a + (x[i] - mean(x)) * (x[i+10] - mean(x))
}
#a <- a / (100 - 10)
a <- a / 100

a / (1/100 * sum((x - mean(x))^2))
acf(x, plot = FALSE)$acf[11]

#rates[rates$t < 21 / 90 & rates$t > 19 / 90, ]
#before - 230 0.221778 1.59094 0.475236
#after  - 231 0.222222 2.37187 0.708512

#the same (lol)
1.59094 / (1.59094 + 0.475236)
2.37187 / (2.37187 + 0.708512)
1.59094 + 0.475236
2.37187 + 0.708512

arrowData1 <- data.frame(
		x0 = c(34, 43, 61, 85, 85), 
		y0 = c(13, 13, 13, 14.5, 11.5), 
		x1 = c(38.5, 47.5, 65.5, 89.5, 89.5), 
		y1 = c(11.5, 11.5, 11.5, 13, 10),
		score = c("1-0", "1-1", "1-2", "3-2", "2-2"))

## arrowData1 <- data.frame(
##         x0 = c(34, 43, 61, 85, 85), 
##         y0 = c(13, 13, 13, 14.5, 11.5), 
##         x1 = c(38.5, 47.5, 65.5, 89.5, 89.5), 
##         y1 = c(11.5, 11.5, 11.5, 13, 10),
##         score = c("1-0", "1-1", "1-2", "3-2", "2-2"))

#arrowData2 <- data.frame(
#		x0 = c(34, 43, 61, 85, 85), 
#		y0 = c(13, 13, 13, 14.5, 11.5) / 15.0, 
#		x1 = c(38.5, 47.5, 65.5, 89.5, 89.5), 
#		y1 = c(11.5, 11.5, 11.5, 13, 10) / 15.0,
#		score = c("1-0", "1-1", "1-2", "3-2", "2-2"))

goalTimes <- rbind(manCityQPR, sunderlandManUtd, StokeBolton, norwichAstonVilla, WiganWolves)
goalTimes <- goalTimes[goalTimes$event %in% 1:2, ]
goalTimes <- goalTimes[!(goalTimes$end %in% c(8, 13, 14, 22, 79, 86)), ]

inplay1 <- ggplot() + geom_vline(aes(xintercept = end, linetype = !(matchID ==1117411)), goalTimes, size = 1.0, col = "grey") +
		geom_line(aes(t * 90.0, home), rates, col = "red", size = 1.0) + 
		geom_line(aes(t * 90.0, away), rates, col = "blue", size = 1.0, linetype = 5) + 
		scale_x_continuous(breaks = seq(0, 90, by = 5)) +
		labs(x = "time (minutes)", y = "rate") +
		## geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), arrowData1, arrow = arrow()) + 
		## geom_text(aes(x = x0, y = y0 + 0.5, label = score), arrowData1) + 
		theme_bw()

inplay2 <- ggplot() + geom_vline(aes(xintercept = end, linetype = !(matchID ==1117411)), goalTimes, size = 1.0, col = "grey") +
		geom_line(aes(t * 90.0, home), alpha, col = "red", size = 1.0) + 
		geom_line(aes(t * 90.0, away), alpha, col = "blue", size = 1.0, linetype = 5) + 
		scale_x_continuous(breaks = seq(0, 90, by = 5)) +
		scale_y_continuous(limits = c(0, 1)) +
		labs(x = "time (minutes)", y = expression(alpha[k](t, w))) +
		theme_bw()

#geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), arrowData2, arrow = arrow()) + 
#geom_text(aes(x = x0, y = y0 + 0.5 / 15, label = score), arrowData2) + 

pdf("/home/jeff/workspace/Latex-Main/Chapter_Utility/ggFigures/inplay2.pdf", height = 9)
grid.arrange(inplay1, inplay2, ncol = 1)
dev.off()







