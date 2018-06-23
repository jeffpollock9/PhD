library(ggplot2)
library(reshape2)
library(gridExtra)

rm(list = ls())

nParticles <- 100000

particles <- runif(nParticles, -5, 5)
weights <- dnorm(particles)
weights <- weights / sum(weights)

expectedCounts <- nParticles * weights
floorExpectedCounts <- floor(expectedCounts)

residualSamples <- nParticles - sum(floorExpectedCounts)
residualWeights <- (expectedCounts - floorExpectedCounts) / residualSamples
residualCounts <- as.vector(rmultinom(1, residualSamples, residualWeights))

newParticles <- rep(particles, times = floorExpectedCounts + residualCounts)

expectedParticles <- rep(particles, times = floorExpectedCounts)
residualParticles <- rep(particles, times = residualCounts)

data <- rbind(
		data.frame(x = expectedParticles, particle = "expected"), 
		data.frame(x = residualParticles, particle = "residual")) 

residualSampling <- ggplot(data, aes(x = x, fill = particle)) + 
		geom_histogram(binwidth = 0.2) + labs(x = "value", y = "frequency") +
		scale_fill_manual(values = c("blue", "red"), guide = "none") +
		theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/residual-sampling.pdf", residualSampling)

##########
##########

samples <- rnorm(5000)

smooth1 <- function(x, samples, h) {
	
	var <- var(samples)
	n <- length(x)
	m <- length(samples)
	density <- numeric(n)
	
	for (i in seq_len(n)) {
		density[i] <- sum(dnorm(x[i], samples, sqrt(h^2 * var))) / m
	}	
	
	return(density)
}

smooth2 <- function(x, samples, h) {
	
	var <- var(samples)
	n <- length(x)
	m <- length(samples)
	density <- numeric(n)
	a <- sqrt(1 - h^2)
	loc <- a * samples + (1 - a) * mean(samples) 
	
	for (i in seq_len(n)) {
		density[i] <- sum(dnorm(x[i], loc, sqrt(h^2 * var))) / m
	}	
	
	return(density)
}

x <- seq(-5, 5, length = 1000)
smoothData <- data.frame(x, s1 = smooth1(x, samples, 0.1), s2 = smooth2(x, samples, 0.1))
smoothDataMelt <- melt(smoothData, id.vars = "x", variable.name = "smooth")

smooth01 <- ggplot() + geom_histogram(aes(x = samples, y = ..density..), binwidth = 0.25) + 
		geom_line(aes(x, value, col = smooth, linetype = smooth), smoothDataMelt, size = 1) + 
		scale_color_manual(values = c("red", "blue"), guide = "none") +
		scale_linetype_manual(values = c("solid", "dashed"), guide = "none") + 
		xlab("x") + theme_bw() 

smoothData <- data.frame(x, s1 = smooth1(x, samples, 0.5), s2 = smooth2(x, samples, 0.5))
smoothDataMelt <- melt(smoothData, id.vars = "x", variable.name = "smooth")

smooth05 <- ggplot() + geom_histogram(aes(x = samples, y = ..density..), binwidth = 0.25) + 
		geom_line(aes(x, value, col = smooth, linetype = smooth), smoothDataMelt, size = 1) + 
		scale_color_manual(values = c("red", "blue"), guide = "none") +
		scale_linetype_manual(values = c("solid", "dashed"), guide = "none") + 
		xlab("x") + theme_bw() 

smoothData <- data.frame(x, s1 = smooth1(x, samples, 0.9), s2 = smooth2(x, samples, 0.9))
smoothDataMelt <- melt(smoothData, id.vars = "x", variable.name = "smooth")

smooth09 <- ggplot() + geom_histogram(aes(x = samples, y = ..density..), binwidth = 0.25) + 
		geom_line(aes(x, value, col = smooth, linetype = smooth), smoothDataMelt, size = 1) + 
		scale_color_manual(values = c("red", "blue"), guide = "none") +
		scale_linetype_manual(values = c("solid", "dashed"), guide = "none") + 
		xlab("x") + theme_bw() 

pdf("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/kernel-smooth.pdf", height = 9)
grid.arrange(smooth01, smooth05, smooth09, ncol = 1)
dev.off()

##########
##########

#particle filter for fixed parameter
#corrected variance over artificial evolution
# data follows N(muTrue, sigmaTrue^2)
# prior for mu is N(5, 1^2)
# estimate mu

set.seed(100)
T <- 25
nParticles <- 50000
x <- matrix(NA_real_, T + 1L, nParticles)
w1 <- numeric(nParticles)
w2 <- numeric(nParticles)
mu <- numeric(T)
sigma2 <- numeric(T)
ybar <- numeric(T)
muTrue <- 1.0
sigmaTrue <- 2.0
sigma02 <- 1.0
mu0 <- 5.0
delta <- 0.95
a <- (3 * delta - 1) / (2 * delta)
h2 <- 1 - a^2

##initialisation
x[1L, ] <- rnorm(nParticles, mu0, sigma02) # samples from prior

n <- 2L
yT <- NULL

for (t in 1:T) {
	#random static data
	y <- rnorm(n, muTrue, sigmaTrue)
	ybar[t] <- mean(y)
	yT <- c(yT, y)
	b <- length(yT)
	
	#theoretical posterior
	mu[t] <- sigma02 / (sigmaTrue^2 / b + sigma02) * mean(yT) + (sigmaTrue^2 / b) / (sigmaTrue^2 / b + sigma02) * mu0
	sigma2[t] <- 1 / (1 / sigma02 + b / sigmaTrue^2)
	
	#first stage weights
	kernalLocations <- a * x[t, ] + (1 - a) * mean(x[t, ]) 
	for (i in seq_len(nParticles)) {
		w1[i] <- sum(dnorm(y, kernalLocations[i], sigmaTrue, TRUE))
	}
	w1 <- exp(w1 - max(w1))
	
	#aux. resample based on first stage weights
	aux <- sample(kernalLocations, nParticles, TRUE, w1)
	
	#artificial evolution
	x[t + 1L, ] <- rnorm(nParticles, aux, sqrt(h2 * var(x[t, ]))) # should be var of aux sample??
	
	#second stage weights
	for (i in seq_len(nParticles)) {
		w2[i] <- sum(dnorm(y, x[t + 1L, i], sigmaTrue, TRUE)) - 
				sum(dnorm(y, aux[i], sigmaTrue, TRUE))
	}
	w2 <- exp(w2 - max(w2))
	w2 <- w2 / sum(w2)
	
	#resample
	x[t + 1L, ] <- sample(x[t + 1L, ], nParticles, TRUE, w2)
}

plotData <- as.data.frame(t(x))
colnames(plotData) <- paste0("time", seq_len(T + 1) - 1)

t <- c(1, 10, 20)
plotList <- vector("list", length(t))

for (i in seq_along(t)) {
	plotList[[i]] <- ggplot(plotData) + geom_histogram(aes_string(paste0("time", t[i]), "..density.."), binwidth = 0.1) + 
			xlim(c(0, 8)) + ylim(c(0, 1.5)) + theme_bw() + xlab(expression(mu)) + 
			stat_function(fun = dnorm, args = list(mean = mu0, sd = sqrt(sigma02)), size = 1, col = "blue") + 
			stat_function(fun = dnorm, args = list(mean = mu[t[i]], sd = sqrt(sigma2[t[i]])), size = 1, col = "green")
}

pdf("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/example-filter.pdf", height = 9)
do.call(grid.arrange, c(plotList, ncol = 1))
dev.off()

##########
##########

#metrics1 <- read.csv("/home/jeff/workspace/C-Football/metrics-real.csv", header = FALSE)
#metrics2 <- read.csv("/home/jeff/workspace/C-Football/metrics-real-full.csv", header = FALSE)
#
#metrics <- rbind(metrics1[-1L, ], metrics2[-1L, ])
#colnames(metrics) <- c("sigma2", "metric")
#
#duplicates <- metrics$metric[metrics$sigma2 == 0.05]
#metrics <- metrics[!metrics$sigma2 == 0.05, ]
#metrics <- rbind(metrics, data.frame(sigma2 = 0.05, metric = mean(duplicates)))
#
#metrics <- metrics[order(metrics[ , 1]), ]

metrics <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/metrics-real-openmp.csv")
metricsSmall <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/metrics-real-openmp-small.csv")
colnames(metrics) <- c("sigma2", "metric")
colnames(metricsSmall) <- c("sigma2", "metric")

#metricPlot <- ggplot(metrics) + geom_line(aes(sigma2, metric), size = 1, col = "red") + 
#		labs(x = expression(sigma^2), y = expression(sum(LSR[w], "w=1", 38))) + theme_bw()

metricPlot <- ggplot(metrics, aes(sigma2, exp(metric / 380))) + geom_line(size = 1, col = "black") + 
		labs(x = expression(sigma^2), y = expression(GM["1, 38"])) + 
		stat_smooth(method = "loess", col = "red", lty = 2, size = 1, se = FALSE) + theme_bw()

#metricPlotSmall <- ggplot(metricsSmall) + geom_line(aes(sigma2, metric), size = 1, col = "red") + 
#		labs(x = expression(sigma^2), y = expression(sum(LSR[w], "w=1", 38))) + theme_bw()

metricPlotSmall <- ggplot(metricsSmall, aes(sigma2, exp(metric / 380))) + geom_line(size = 1, col = "black") + 
		labs(x = expression(sigma^2), y = expression(GM["1, 38"])) + 
		stat_smooth(method = "loess", col = "red", lty = 2, size = 1, se = FALSE) + theme_bw()

metrics[which.max(metrics$metric), ]
metricsSmall[which.max(metricsSmall$metric), ]

#ggsave("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/metricPlot.pdf", metricPlot)
#ggsave("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/metricPlotSmall.pdf", metricPlotSmall)

pdf("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/metricPlotComb.pdf", height = 9)
grid.arrange(metricPlot, metricPlotSmall, ncol = 1)
dev.off()

##########
##########

names <- c("h", "a", "c_-1", "c_0", "c_1", "d_-1", "d_0", 
		"d_1", "rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Brom", "Wigan Athletic", "Wolverhampton")

graphNames <- as.list(names)
graphNames[3:5] <- list(expression(c[-1]), expression(c[0]), expression(c[1]))
graphNames[6:8] <- list(expression(d[-1]), expression(d[0]), expression(d[1]))
graphNames[9:10] <- list(expression(rho[1]), expression(rho[2]))

particles <- vector("list", 39L)
for (t in 1:39) {
	particles[[t]] <- read.csv(paste0("/home/jeff/workspace/Latex-Main/Chapter_Particle/particles/week-", t-1, ".csv"), 
			header = FALSE, colClasses = "numeric")
	colnames(particles[[t]]) <- paste0("col", 1:30)
	particles[[t]]$week <- t - 1L
}

particlesData <- do.call(rbind, particles)

for (i in 0:4) {
	
	plotList <- vector("list", 6)
	for (j in 1:6) {
		k <- j + i * 6
		agg <- aggregate(as.formula(paste0("col", k ," ~ week")), particlesData, function(x) {
					c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
				})
		
		aggData <- as.data.frame(agg[ , 2])
		aggData$week <- 0:38
		aggDataMelt <- melt(aggData, id.vars = "week")
		
		plotList[[j]] <- ggplot(aggDataMelt) + 
				geom_line(aes(week, value, group = variable, col = variable, linetype = variable)) + 
				scale_color_manual(values = c("black", "red", "black"), guide = "none") + 
				scale_linetype_manual(values = c(3, 1, 3), guide = "none") + 
				scale_x_discrete(breaks = seq(0, 38, by = 4)) + ylab(graphNames[[k]]) + 
				theme_bw()
	}
	
	pdf(paste0("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/filter-", i,".pdf"), height = 9)
	do.call(grid.arrange, c(plotList, ncol = 2))
	dev.off()
}

##########
##########

modelProbsDyn <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/dynamic-probs-openmp.csv", header = FALSE)
modelProbs <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/output/modelProbsM.csv", header = FALSE)[ , 1:3]

colnames(modelProbsDyn) <- paste0("DM", c("home", "draw", "away"))
colnames(modelProbs) <- paste0("M", c("home", "draw", "away"))

data <- cbind(modelProbs, modelProbsDyn)
data$match <- 1:380
data <- tail(data, 330)

pl1 <- ggplot(data, aes(Mhome, DMhome)) + geom_point() + geom_abline(intercept = 0, gradient = 1, col = "red", linetype = 2, size = 1) + 
		labs(x = "model M home win probability", y = "model DM home win probability") + 
		scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
		scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
		stat_smooth(method = "loess", col = "blue", lty = 1, size = 1, se = FALSE) + theme_bw()

pl2 <- ggplot(data, aes(Mdraw, DMdraw)) + geom_point() + geom_abline(intercept = 0, gradient = 1, col = "red", linetype = 2, size = 1) + 
		labs(x = "model M draw probability", y = "model DM draw probability") + 
		scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 0.5)) + 
		scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 0.5)) + 
		stat_smooth(method = "loess", col = "blue", lty = 1, size = 1, se = FALSE) + theme_bw()

pl3 <- ggplot(data, aes(Maway, DMaway)) + geom_point() + geom_abline(intercept = 0, gradient = 1, col = "red", linetype = 2, size = 1) + 
		labs(x = "model M away win probability", y = "model DM away win probability") + 
		scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) + 
		scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
		stat_smooth(method = "loess", col = "blue", lty = 1, size = 1, se = FALSE) + theme_bw()

pdf("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/MandDM.pdf", height = 9)
grid.arrange(pl1, pl2, pl3, ncol = 1)
dev.off()

##########
##########

mh <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/mh.txt", header = FALSE)
pf <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/pf.txt", header = FALSE)

colnames(mh) <- c("method", "time")
colnames(pf) <- c("method", "time")

mh[ , 1L] <- "mh"
pf[ , 1L] <- "pf"

mh$week <- 1:38
pf$week <- 1:38

times <- rbind(mh, pf)

pl1 <- ggplot() + geom_line(aes(week, time, col = method), times, size = 1) +
		scale_x_discrete(breaks = seq(0, 38, by = 2)) + 
		scale_y_continuous(breaks = c(seq(0, 13000, by = 2500)), name = "cumulative time (seconds)") +
		scale_color_manual(values = c("black", "red"), guide = "none") + 
		theme_bw()

pl2 <- ggplot() + geom_line(aes(week, time, col = method), times, size = 1) + 
		scale_x_discrete(breaks = seq(0, 38, by = 2)) + 
		scale_y_continuous(name = "cumulative time (seconds)", limits = c(0, 300)) +
		scale_color_manual(values = c("black", "red"), guide = "none") + 
		theme_bw()

pl2 <- pl1 + coord_cartesian(ylim = c(-15, 315)) + 
		scale_y_continuous(breaks = c(seq(0, 500, by = 100)), name = "cumulative time (seconds)")

pdf("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/runtimes.pdf", height = 9)
grid.arrange(pl1, pl2, ncol = 1)
dev.off()

##########
##########

#metrics <- read.csv("/home/jeff/workspace/C-Football/metrics-simulated.csv", header = FALSE)[-1 , ]
metrics <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/metrics-simulated-openmp.csv", header = FALSE)[-1 , ]
metrics <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Particle/data/metrics-simulated-openmp.csv", header = FALSE)[-1 , ]

colnames(metrics) <- c("sigma2", "metric")

metricPlotSim <- ggplot(metrics) + geom_line(aes(sigma2, exp(metric / 380)), size = 1, col = "red") + 
		labs(x = expression(sigma^2), y = "GM") + theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/metricPlotSim.pdf", metricPlotSim)

metrics[which.max(metrics$metric), ]
exp(-304.522/380)
##########
##########

R <- read.csv("/home/jeff/workspace/Latex-Particle-Filter/data/R-simulated.csv", header = FALSE)
colnames(R) <- paste0("team", 1:20)
truePars <- c(0.4, 0.1, 0.7, 0.4, 0.6, 0.9, 0.5, 0.3, 1.0, 1.5)
R <- rbind(R[1 , ] , R)
R$week <- 0:38

particles <- vector("list", 39L)
for (t in 1:39) {
	particles[[t]] <- read.csv(paste0("/home/jeff/workspace/Latex-Main/Chapter_Particle/particles/sim-week-", t-1, ".csv"), 
			header = FALSE, colClasses = "numeric")
	colnames(particles[[t]]) <- paste0("col", 1:30)
	particles[[t]]$week <- t - 1L
}

particlesData <- do.call(rbind, particles)

for (i in 0:4) {
	
	plotList <- vector("list", 6)
	for (j in 1:6) {
		k <- j + i * 6
		agg <- aggregate(as.formula(paste0("col", k ," ~ week")), particlesData, function(x) {
					c(quantile(x, 0.025), mean(x), quantile(x, 0.975))
				})
		
		aggData <- as.data.frame(agg[ , 2])
		aggData$week <- 0:38
		aggDataMelt <- melt(aggData, id.vars = "week")
		
		plotList[[j]] <- ggplot(aggDataMelt) + 
				geom_line(aes(week, value, group = variable, col = variable, linetype = variable)) + 
				scale_color_manual(values = c("black", "red", "black")) + 
				scale_linetype_manual(values = c(3, 1, 3)) + 
				scale_x_discrete(breaks = seq(0, 38, by = 4)) + ylab(graphNames[[k]]) + 
				theme_bw() + theme(legend.position = "none")
		if (k > 10) {
			plotList[[j]] <- plotList[[j]] + geom_line(aes_string("week", paste0("team", k-10)), R, col = "green", linetype = 1)
		} else {
			plotList[[j]] <- plotList[[j]] + geom_hline(yintercept = truePars[k], col = "green", linetype = 1)
		}
	}
	
	pdf(paste0("/home/jeff/workspace/Latex-Main/Chapter_Particle/ggFigures/sim-filter-", i,".pdf"), height = 9)
	do.call(grid.arrange, c(plotList, ncol = 2))
	dev.off()
}

