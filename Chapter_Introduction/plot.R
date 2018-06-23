library(ggplot2)
library(reshape2)

rm(list = ls())

##########
##########

events <- read.csv("/home/jeff/workspace/Latex-Main/Chapter_Adaptive/data/2011-2012-data.csv", stringsAsFactors = FALSE)

goalTimes <- ggplot(subset(events, event %in% 1:2), aes(x = end)) + geom_bar() + 
		xlab("time") + ylab("frequency") + scale_x_discrete(breaks = seq(0, 90, by = 5)) + theme_bw()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Introduction/ggFigures/goal-times.pdf", goalTimes)

##########
##########

odds <- read.csv("/home/jeff/workspace/FootballAnalysisAndData/data/E0-2011-2012.csv", stringsAsFactors = FALSE)
odds$match <- seq_len(nrow(odds))

#n <- 50
#odds <- odds[sample.int(nrow(odds), n), ]
odds <- odds[order(odds$B365H), ]

odds1 <- melt(odds[ , c("match", "B365H", "B365D", "B365A")], id.var= "match")
#odds1$match <- 1:n
odds1$match <- 1:380

oddsPlot <- ggplot(odds1, aes(match, 1 / value, fill = variable)) + 
		geom_bar(stat = "identity") + 
		geom_hline(yintercept = 1.0, linetype = "dashed") + 
		scale_fill_manual(values = c("blue", "red", "green"), guide = "none") +
		scale_y_continuous(breaks = seq(0, 1.1, by = 0.1), limits = c(0, 1.1)) + 
		scale_x_discrete(breaks = seq(0, 380, by = 20)) + theme_bw() + 
		xlab("m") + 
		ylab(substitute(list(1/d[H] + 1/d[D] + 1/d[A]), list(H = "m, H", D = "m, D", A = "m, A"))) + 
		coord_flip()

ggsave("/home/jeff/workspace/Latex-Main/Chapter_Introduction/ggFigures/odds.pdf", oddsPlot, height = 9)





