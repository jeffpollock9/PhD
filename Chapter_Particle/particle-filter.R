library(ggplot2)
setwd("/home/jeff/workspace/Latex-Main/Chapter_Particle/Figures")


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
		geom_histogram(breaks = seq(-4, 4, by = 0.2)) + 
		scale_fill_manual(values = c("black", "grey")) +
		theme_classic() + xlab("value") + 
		theme(legend.position = "none")
ggsave("residual-sampling.pdf", residualSampling)
###

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

pdf("kernel-smooth-h01.pdf")
par(las = 1)
h <- 0.1
hist(samples, breaks = "fd", prob = TRUE, main = "")
#legend("topright", 
#		legend = c("no shrinkage", "shrinkage", "theoretical"), 
#		col = c("red", "blue", "black"), 
#		lty = c(1, 2, 2))
curve(smooth1(x, samples, h), add = TRUE, col = "red", lty = 1, lwd = 2)
curve(smooth2(x, samples, h), add = TRUE, col = "blue", lty = 2, lwd = 2)
curve(dnorm(x), add = TRUE, col = "black", lty = 2, lwd = 2)
dev.off()

###
pdf("kernel-smooth-h05.pdf")
par(las = 1)
h <- 0.5
par(las = 1)
hist(samples, breaks = "fd", prob = TRUE, main = "")
#legend("topright", 
#		legend = c("no shrinkage", "shrinkage", "theoretical"), 
#		col = c("red", "blue", "black"), 
#		lty = c(1, 2, 2))
curve(smooth1(x, samples, h), add = TRUE, col = "red", lty = 1, lwd = 2)
curve(smooth2(x, samples, h), add = TRUE, col = "blue", lty = 2, lwd = 2)
curve(dnorm(x), add = TRUE, col = "black", lty = 2, lwd = 2)
dev.off()

###

pdf("kernel-smooth-h09.pdf")
par(las = 1)
h <- 0.9
par(las = 1)
hist(samples, breaks = "fd", prob = TRUE, main = "")
#legend("topright", 
#		legend = c("no shrinkage", "shrinkage", "theoretical"), 
#		col = c("red", "blue", "black"), 
#		lty = c(1, 2, 2))
curve(smooth1(x, samples, h), add = TRUE, col = "red", lty = 1, lwd = 2)
curve(smooth2(x, samples, h), add = TRUE, col = "blue", lty = 2, lwd = 2)
curve(dnorm(x), add = TRUE, col = "black", lty = 2, lwd = 2)
dev.off()
###

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

logLikelihood <- function(x, y) {
	len <- length(x)
	out <- numeric(len)
	for (i in seq_len(len)) {
		out[i] <- sum(dnorm(y, x[i], 2, TRUE))
	}
	return(exp(out))
}

###animation
if (FALSE) {
	library(animation)
	saveLatex({
				for (i in 0:T) {
					hist(x[i + 1 , ], 
							xlim = c(-2, 10), ylim = c(0, 1.8), xlab = "x",
							main = paste("t =", i), breaks = 25, probability = TRUE)
					legend("topright", 
							legend = c("prior\n", "likelihood of\n2 new points\n", "theoretical\nposterior\n"), 
							col = c("blue", "red", "green"), lty = c(1, 1, 1))
					curve(dnorm(x, 5, sqrt(sigma02)), add = TRUE, col = "blue", lwd = 2)
					if (i > 0) {
						curve(dnorm(x, mu[i], sqrt(sigma2[i])), add = TRUE, col = "green", lwd = 2)
						curve(logLikelihood(x, yT[1:n + (i-1)*n]), add = TRUE, col = "red", lwd = 2)
					}
					
					ani.pause()
				}
			}, outdir = getwd(), img.name = "particle-filter", 
			htmlfile = "particle-filter.html", autobrowse = TRUE, 
			title = "aux particle filter")
}

pdf("example-filter-t1.pdf")
par(las = 1)
i <- 1
hist(x[i + 1 , ], 
		xlim = c(-2, 10), ylim = c(0, 1.8), xlab = "x",
		main = "", breaks = 25, probability = TRUE)
#legend("topright", 
#		legend = c("prior\n", "likelihood of\n2 new points\n", "theoretical\nposterior\n"), 
#		col = c("blue", "red", "green"), lty = c(1, 1, 1))
curve(dnorm(x, 5, sqrt(sigma02)), add = TRUE, col = "blue", lwd = 2)
if (i > 0) {
	curve(dnorm(x, mu[i], sqrt(sigma2[i])), add = TRUE, col = "green", lwd = 2)
	curve(logLikelihood(x, yT[1:n + (i-1)*n]), add = TRUE, col = "red", lwd = 2)
}
dev.off()

###

pdf("example-filter-t10.pdf")
par(las = 1)
i <- 10
hist(x[i + 1 , ], 
		xlim = c(-2, 10), ylim = c(0, 1.8), xlab = "x",
		main = "", breaks = 25, probability = TRUE)
#legend("topright", 
#		legend = c("prior\n", "likelihood of\n2 new points\n", "theoretical\nposterior\n"), 
#		col = c("blue", "red", "green"), lty = c(1, 1, 1))
curve(dnorm(x, 5, sqrt(sigma02)), add = TRUE, col = "blue", lwd = 2)
if (i > 0) {
	curve(dnorm(x, mu[i], sqrt(sigma2[i])), add = TRUE, col = "green", lwd = 2)
	curve(logLikelihood(x, yT[1:n + (i-1)*n]), add = TRUE, col = "red", lwd = 2)
}
dev.off()

###
pdf("example-filter-t20.pdf")
par(las = 1)
i <- 20
hist(x[i + 1 , ], 
		xlim = c(-2, 10), ylim = c(0, 1.8), xlab = "x",
		main = "", breaks = 25, probability = TRUE)
#legend("topright", 
#		legend = c("prior\n", "likelihood of\n2 new points\n", "theoretical\nposterior\n"), 
#		col = c("blue", "red", "green"), lty = c(1, 1, 1))
curve(dnorm(x, 5, sqrt(sigma02)), add = TRUE, col = "blue", lwd = 2)
if (i > 0) {
	curve(dnorm(x, mu[i], sqrt(sigma2[i])), add = TRUE, col = "green", lwd = 2)
	curve(logLikelihood(x, yT[1:n + (i-1)*n]), add = TRUE, col = "red", lwd = 2)
}
dev.off()

###
metrics <- as.matrix(read.csv("/home/jeff/workspace/C-Football/metrics-real.csv", header = FALSE))
metrics1 <- as.matrix(read.csv("/home/jeff/workspace/C-Football/metrics-real-full.csv", header = FALSE))

#TODO: combine metrics and metrics1 when you can be bothered (plot already done)
#pdf("metrics.pdf")
#par(las = 1)
#plot(metrics[2:nrow(metrics), 1], metrics[2:nrow(metrics), 2], 
#		type = "l", lwd = 2, col = "red", 
#		xlab = expression(sigma^2), ylab = "M")
#grid()
#dev.off()

metrics <- as.matrix(read.csv("/home/jeff/workspace/C-Football/metrics-simulated.csv", header = FALSE))

###

pdf("metrics-simulated.pdf")
par(las = 1)
plot(metrics[2:nrow(metrics), 1], metrics[2:nrow(metrics), 2], 
		type = "l", lwd = 2, col = "red", 
		xlab = expression(sigma^2), ylab = "M")
grid()
dev.off()
###

names1 <- c("h", "a", "c_-1", "c_0", "c_1", "d_-1", "d_0", 
		"d_1", "rho1", "rho2", "Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", 
		"Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", 
		"Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", 
		"Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", 
		"West Brom", "Wigan Athletic", "Wolverhampton")
graphNames <- as.list(names1)
graphNames[3:5] <- list(expression(c[-1]), expression(c[0]), expression(c[1]))
graphNames[6:8] <- list(expression(d[-1]), expression(d[0]), expression(d[1]))
graphNames[9:10] <- list(expression(rho[1]), expression(rho[2]))

##particle output
particles <- list()
for (t in 1:39) {
	particles[[t]] <- read.csv(paste0("/home/jeff/Documents/particles/week-", t-1, ".csv"), 
			header = FALSE, colClasses = "numeric")
}

plotTimeSeries <- function(par) {
	points <- matrix(NA_real_, 39, 3)
	for (t in 0:38) {
		points[t+1 , ] <- quantile(particles[[t+1]][, par], c(0.025, 0.5, 0.975))
	}
#	if (par > 10) {
#		points <- exp(points)
#	}
	par(las = 1)
	plot(c(0, 38), range(points), main = "", type = "n", 
			xlab = "T", ylab = graphNames[[par]])
	grid()
	
	lines(0:38, points[ , 1], lty = 2)
	lines(0:38, points[ , 2], col = "red", lwd = 2)
	lines(0:38, points[ , 3], lty = 2)
	
	return(invisible(points))
}

###
for (i in 0:4) {
	k <- 1L + i * 6L
	pdf(paste0("filter-", k, "-", k + 5L, ".pdf"), width = 12, height = 18)
	par(mfrow = c(3, 2), cex = 1.2)
	for (j in k:(k+5)) {
		plotTimeSeries(j)
	}
	dev.off()
}

plotBivariateDensity2 <- function(par1, par2, t) {
	
	require(ggplot2)
	require(RColorBrewer)
	require(gridExtra)
	
	xlim <- range(particles[[1]][ , par1])
	ylim <- range(particles[[1]][ , par2])
	
	dat <- particles[[t+1]]
	
	p <- ggplot(dat, aes(dat[, par1], dat[, par2]), environment = environment()) + 
			stat_binhex(bins = 25) + 
			theme_classic() + 
			scale_fill_gradientn(colours = brewer.pal(9, "Reds")[2:9]) + 
			theme(legend.position="none") + 
			xlab(graphNames[[par1]]) + ylab(graphNames[[par2]]) + 
			xlim(xlim) + ylim(ylim) + 
			ggtitle(paste("t =", t))
	
	return(p)
}

pairs <- list(a = c(1, 3), b = c(2, 4), c = c(5, 8), d = c(6, 10))
for (p in pairs) {
	times <- c(2, 5, 10, 15, 30, 35)
	plots <- list()
	for (i in seq_along(times)) {
		plots[[i]] <- plotBivariateDensity2(p[1], p[2], times[i])
	}
	multi <- do.call(arrangeGrob, plots)
	ggsave(paste0("filter-bivariate-", p[1], "-", p[2], ".png"), multi, width = 12, height = 18)
}

### simulated data ### 

R <- read.csv("/home/jeff/workspace/Latex-Particle-Filter/data/R-simulated.csv", header = FALSE)

colnames(R) <- names1[11:30]

#par(mfrow = c(2, 3))
#k <- sample(1:20, 6, FALSE)
#for (i in seq_along(k)) {
#	week <- 1:38
#	par(las =1)
#	plot(week, R[ , k[i]][week], type = "l", col = "green", lwd = 2,
#			xlab = "T", ylab = names1[11:30][k[i]], main = "")
#	grid()
#}

#particle output
particles <- list()
for (t in 1:39) {
	particles[[t]] <- read.csv(paste0("/home/jeff/Documents/particles/sim-week-", t-1, ".csv"), 
			header = FALSE, colClasses = "numeric")
	colnames(particles[[t]]) <- names1
}

plotTimeSeries <- function(par) {
	points <- matrix(NA_real_, 39, 3)
	for (t in 1:39) {
		points[t , ] <- quantile(particles[[t]][, par], c(0.025, 0.5, 0.975))
	}
	
#	if (par > 10) {
#		points <- exp(points)
#	}
	
	par(las = 1)
	plot(c(0, 38), range(points), ylab = graphNames[[par]], type = "n", 
			main = "", xlab = "T")
	grid()
	
	lines(0:38, points[ , 1], lty = 3)
	lines(0:38, points[ , 2], col = "black", lwd = 2)
	lines(0:38, points[ , 3], lty = 3)
	if (par > 10) {
#		lines(1:38, exp(R[ , par - 10]), col = "green", lwd = 2, pch = 16, type = "o")
		lines(0:38, c(R[1, par - 10], R[ , par - 10]), col = "green", lwd = 2)
	} else {
		truePars <- c(0.4, 0.1, 0.7, 0.4, 0.6, 0.9, 0.5, 0.3, 1.0, 1.5)
		abline(h = truePars[par], col = "green", lwd = 2)
	}
	
	return(invisible(points))
}

#true values
#double c[3] = { 0.7, 0.4, 0.6 };
#double d[3] = { 0.9, 0.5, 0.3 };
#double rho[2] = { 1.0, 1.5 };
#
#Model9 model(0.4, 0.1, c, d, rho);

#k <- 5:10
#k <- sample(1:30, 6, FALSE)
#par(mfrow = c(2, length(k) / 2))
#for (i in seq_along(k)) {
#	plotTimeSeries(k[i])
#}

for (i in 0:4) {
	k <- 1L + i * 6L
	pdf(paste0("sim-filter-", k, "-", k + 5L, ".pdf"), width = 12, height = 18)
	par(mfrow = c(3, 2), cex = 1.2)
	for (j in k:(k+5)) {
		plotTimeSeries(j)
	}
	dev.off()
}

##plots for talk
if (FALSE) {
	pdf(paste0("sim-filter-talk", 1, ".pdf"), width = 18)
	par(mfrow = c(1, 2), cex = 1.2)
	for (j in c(2, 8)) {
		plotTimeSeries(j)
	}
	dev.off()

	pdf(paste0("sim-filter-talk", 2, ".pdf"), width = 18)
	par(mfrow = c(1, 2), cex = 1.2)
	for (j in c(13, 18)) {
		plotTimeSeries(j)
	}
	dev.off()
}

### bivariate plots ###

if (FALSE) {
	plotDensity <- function(par, t, xlim, ylim) {
		points <- particles[[t + 1]][ , par]
		if (par > 10) {
			points <- exp(points)
		}
		hist(points, breaks = "fd", main = paste0("parameter: ", graphNames[par], "\ntime: ", t), 
				probability = TRUE, xlab = "value", xlim = xlim, ylim = ylim)
		
		return(invisible(points))
	}
	library(animation)
	saveHTML({
				for (t in 0:38) {
					plotDensity(5, t, c(0, 1.2), c(0, 5.0))
				}
			})
}
plotBivariateDensity1 <- function(par1, par2, t, ...) {
	
	x <- particles[[t+1]][ , par1]
	y <- particles[[t+1]][ , par2]
	
	xlim <- range(particles[[1]][ , par1])
	ylim <- range(particles[[1]][ , par2])
	plot(x, y, main = paste("t =", t), xlab = graphNames[par1], 
			ylab = graphNames[par2], xlim = xlim, ylim = ylim, 
			col = rgb(1, 0, 0, 0.01), pch = 16, ...)
	grid()
	truePars <- c(0.4, 0.1, 0.7, 0.4, 0.6, 0.9, 0.5, 0.3, 1.0, 1.5)
	points(truePars[par1], truePars[par2], col = "green", pch = 4, cex = 2, lwd = 3)
	
	invisible(NULL)
}

plotBivariateDensity2 <- function(par1, par2, t) {
	
	require(ggplot2)
	require(RColorBrewer)
	require(gridExtra)
	
	xlim <- range(particles[[1]][ , par1])
	ylim <- range(particles[[1]][ , par2])
	
	truePars <- c(0.4, 0.1, 0.7, 0.4, 0.6, 0.9, 0.5, 0.3, 1.0, 1.5)
	dat <- particles[[t+1]]
	
	p <- ggplot(dat, aes(dat[, par1], dat[, par2]), environment = environment()) + 
			stat_binhex(bins = 25) + 
			theme_classic() + 
			scale_fill_gradientn(colours = brewer.pal(9, "Reds")[2:9]) + 
			theme(legend.position="none") + 
			xlab(graphNames[[par1]]) + ylab(graphNames[[par2]]) + 
			xlim(xlim) + ylim(ylim) + 
			ggtitle(paste("t =", t)) +
			geom_point(aes(truePars[par1], truePars[par2]), color= "green", size = 5, shape = 4) 
	
	return(p)
}

pairs <- list(a = c(1, 3), b = c(2, 4), c = c(5, 8), d = c(6, 10))
for (p in pairs) {
	times <- c(2, 5, 10, 15, 30, 35)
	plots <- list()
	for (i in seq_along(times)) {
		plots[[i]] <- plotBivariateDensity2(p[1], p[2], times[i])
	}
	multi <- do.call(arrangeGrob, plots)
	ggsave(paste0("sim-filter-bivariate-", p[1], "-", p[2], ".png"), multi, 
			width = 12, height = 18)
}
#	png(paste0("sim-filter-bivariate-", p[1], "-", p[2], ".png"), width = 12, height = 18)
#	dev.off()

if (FALSE) {
	saveHTML({
				for (t in 0:38) {
					plotBivariateDensity(6, 8, t)
				}
			})
}
