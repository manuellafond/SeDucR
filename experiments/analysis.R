library(ggplot2)
theme_set(theme_bw(base_size=20))

dat <- read.csv("stats_kowhai.csv", header=T)

#default values
nH <- 50
nP <- 20
rB <- 2
pJ <- 0.5
d <- 10

#may need to read from file
reps <- 10

###varying nH
nHdat <- dat[dat$np == nP & dat$rb == rB & dat$pj == pJ & dat$dup_cost == d,]
attach(nHdat)

nHmeans <- aggregate(nHdat, by=list(nh, dup_cost, method), FUN=mean)
nHsds <- aggregate(nHdat, by=list(nh, dup_cost, method), FUN=sd)

#pdf("figures/nH-cost.pdf")
ggplot(data=nHmeans) + geom_point(aes(nh,solution_cost,col=Group.3)) +
	geom_errorbar(aes(nh,ymin=solution_cost-2*nHsds$solution_cost/sqrt(reps),ymax=solution_cost+2*nHsds$solution_cost/sqrt(reps),col=Group.3), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	xlab("nH") + ylab("Cost")
#dev.off()

#pdf("figures/nH-time.pdf")
ggplot(data=nHmeans) + geom_point(aes(nh,time, col=Group.3)) +
	geom_errorbar(aes(nh,ymin=time-2*nHsds$time/sqrt(reps),ymax=time+2*nHsds$time/sqrt(reps),col=Group.3), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) + scale_y_continuous(trans="log") +
	xlab("nH") + ylab("Time")
#dev.off()
