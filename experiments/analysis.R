library(ggplot2)
library(scales)
theme_set(theme_bw(base_size=20))

#helper function for axis breaks
lb = function(maj, by=1, lower=0.1, radix=10) {
  function(x) {
    minx         = floor(max(-10,min(logb(x,radix), na.rm=T))) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(seq(by,radix-by,by=by),radix)
      breaks = rep(steps, times=n_major) +
               rep(major_breaks, each=floor(radix/by)-1)
    }
    radix^breaks
  }
}


dat <- read.csv("stats_kowhai.csv", header=T)

#default values
nH <- 50
nP <- 20
rB <- 2
pJ <- 0.5
d <- 10

#may need to read from file
reps <- 10

### 1. varying nH
nHdat <- dat[dat$np == nP & dat$rb == rB & dat$pj == pJ & dat$dup_cost == d,]
attach(nHdat)

nHmeans <- aggregate(nHdat, by=list(nh, method), FUN=mean)
nHsds <- aggregate(nHdat, by=list(nh, method), FUN=sd)

#pdf("figures/nH-cost.pdf")
ggplot(data=nHmeans) + geom_point(aes(nh,solution_cost,col=Group.3)) +
	geom_errorbar(aes(nh,ymin=solution_cost-2*nHsds$solution_cost/sqrt(reps),ymax=solution_cost+2*nHsds$solution_cost/sqrt(reps),col=Group.3), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	xlab("nH") + ylab("Cost")
#dev.off()

#compare lca cost to insider cost
nHcostlca <- nHdat[nHdat$method == "lca",]
nHcostinsider <- nHdat[nHdat$method=="insider",]
nHcostlca$ratio <- nHcostlca$solution_cost/nHcostinsider$solution_cost
nHcost <- aggregate(nHcostlca, by=list(nHcostlca$nh), FUN=mean)
nHcostsds <- aggregate(nHcostlca, by=list(nHcostlca$nh), FUN=sd)

pdf("figures/nH-cost-insidervlca.pdf")
ggplot(data=nHcost) + geom_point(aes(nh,ratio)) +
	geom_errorbar(aes(nh,ymin=ratio-2*nHcostsds$ratio/sqrt(reps),ymax=ratio+2*nHcostsds$ratio/sqrt(reps)), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	xlab(bquote(n[H])) + ylab("Cost ratio")
dev.off()

pdf("figures/nH-time.pdf")
ggplot(data=nHmeans) + geom_point(aes(nh,time, col=Group.2)) +
	geom_errorbar(aes(nh,ymin=time-2*nHsds$time/sqrt(reps),ymax=time+2*nHsds$time/sqrt(reps),col=Group.2), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual(values=c(3,1,2,4), labels=c("FastMultRec","inSiDeR","LCA","segdup")) +
	theme(legend.title=element_blank(), legend.position = c(.95, .75), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(n[H])) + ylab("Time (s)")
dev.off()


### 2. varying d
ddat <- dat[dat$np == nP & dat$rb == rB & dat$pj == pJ & dat$nh == nH,]
attach(ddat)

dmeans <- aggregate(ddat, by=list(dup_cost, method), FUN=mean)
dsds <- aggregate(ddat, by=list(dup_cost, method), FUN=sd)

#compare lca cost to insider cost
dcostlca <- ddat[ddat$method == "lca",]
dcostinsider <- ddat[ddat$method=="insider",]
dcostlca$ratio <- dcostlca$solution_cost/dcostinsider$solution_cost
dcost <- aggregate(dcostlca, by=list(dcostlca$dup_cost), FUN=mean)
dcostsds <- aggregate(dcostlca, by=list(dcostlca$dup_cost), FUN=sd)

pdf("figures/d-cost-insidervlca.pdf")
ggplot(data=dcost) + geom_point(aes(dup_cost,ratio)) +
	geom_errorbar(aes(dup_cost,ymin=ratio-2*dcostsds$ratio/sqrt(reps),ymax=ratio+2*dcostsds$ratio/sqrt(reps)), width=1) + 
	scale_x_continuous(breaks=seq(0,100,by=5)) +
	xlab(bquote(delta)) + ylab("Cost ratio")
dev.off()

pdf("figures/d-time.pdf")
ggplot(data=dmeans) + geom_point(aes(dup_cost,time, col=Group.2)) +
	geom_errorbar(aes(dup_cost,ymin=time-2*dsds$time/sqrt(reps),ymax=time+2*dsds$time/sqrt(reps),col=Group.2), width=1) + 
	scale_x_continuous(breaks=seq(0,100,by=5)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual(values=c(3,1,2,4), labels=c("FastMultRec","inSiDeR","LCA","segdup")) +
	theme(legend.title=element_blank(), legend.position = c(.95, .85), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(delta)) + ylab("Time (s)")
dev.off()


### 3. varying nP
nPdat <- dat[dat$dup_cost == d & dat$rb == rB & dat$pj == pJ & dat$nh == nH,]
attach(nPdat)

nPmeans <- aggregate(nPdat, by=list(np, method), FUN=mean)
nPsds <- aggregate(nPdat, by=list(np, method), FUN=sd)

#compare lca cost to insider cost
nPcostlca <- nPdat[nPdat$method == "lca",]
nPcostinsider <- nPdat[nPdat$method=="insider",]
nPcostlca$ratio <- nPcostlca$solution_cost/nPcostinsider$solution_cost
nPcost <- aggregate(nPcostlca, by=list(nPcostlca$np), FUN=mean)
nPcostsds <- aggregate(nPcostlca, by=list(nPcostlca$np), FUN=sd)

pdf("figures/nP-cost-insidervlca.pdf")
ggplot(data=nPcost) + geom_point(aes(np,ratio)) +
	geom_errorbar(aes(np,ymin=ratio-2*nPcostsds$ratio/sqrt(reps),ymax=ratio+2*nPcostsds$ratio/sqrt(reps)), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) +
	xlab(bquote(n[P])) + ylab("Cost ratio")
dev.off()

pdf("figures/nP-time.pdf")
ggplot(data=nPmeans) + geom_point(aes(np,time, col=Group.2)) +
	geom_errorbar(aes(np,ymin=time-2*nPsds$time/sqrt(reps),ymax=time+2*nPsds$time/sqrt(reps),col=Group.2), width=2) + 
	scale_x_continuous(breaks=seq(0,100,by=20)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual(values=c(3,1,2,4), labels=c("FastMultRec","inSiDeR","LCA","segdup")) +
	theme(legend.title=element_blank(), legend.position = c(.95, .75), legend.justification = c("right", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(n[P])) + ylab("Time (s)")
dev.off()


### 4. varying rB
rBdat <- dat[dat$dup_cost == d & dat$np == nP & dat$pj == pJ & dat$nh == nH,]
attach(rBdat)

rBmeans <- aggregate(rBdat, by=list(rb, method), FUN=mean)
rBsds <- aggregate(rBdat, by=list(rb, method), FUN=sd)

#compare lca cost to insider cost
rBcostlca <- rBdat[rBdat$method == "lca",]
rBcostinsider <- rBdat[rBdat$method=="insider",]
rBcostlca$ratio <- rBcostlca$solution_cost/rBcostinsider$solution_cost
rBcost <- aggregate(rBcostlca, by=list(rBcostlca$rb), FUN=mean)
rBcostsds <- aggregate(rBcostlca, by=list(rBcostlca$rb), FUN=sd)

pdf("figures/rB-cost-insidervlca.pdf")
ggplot(data=rBcost) + geom_point(aes(rb,ratio)) +
	geom_errorbar(aes(rb,ymin=ratio-2*rBcostsds$ratio/sqrt(reps),ymax=ratio+2*rBcostsds$ratio/sqrt(reps)), width=0.2) + 
	scale_x_continuous(breaks=seq(0,100,by=1)) +
	xlab(bquote(r[B])) + ylab("Cost ratio")
dev.off()

pdf("figures/rB-time.pdf")
ggplot(data=rBmeans) + geom_point(aes(rb,time, col=Group.2),position=position_dodge(width=0.2)) +
	geom_errorbar(aes(rb,ymin=pmax(0,time-2*rBsds$time/sqrt(reps)),ymax=time+2*rBsds$time/sqrt(reps),col=Group.2), width=0.2,position=position_dodge(width=0.2)) + 
	scale_x_continuous(breaks=seq(0,100,by=1)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual(values=c(3,1,2,4), labels=c("FastMultRec","inSiDeR","LCA","segdup")) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(r[B])) + ylab("Time (s)")
dev.off()


### 5. varying pJ
pJdat <- dat[dat$dup_cost == d & dat$np == nP & dat$rb == rB & dat$nh == nH,]
attach(pJdat)

pJmeans <- aggregate(pJdat, by=list(pj, method), FUN=mean)
pJsds <- aggregate(pJdat, by=list(pj, method), FUN=sd)

#compare lca cost to insider cost
pJcostlca <- pJdat[pJdat$method == "lca",]
pJcostinsider <- pJdat[pJdat$method=="insider",]
pJcostlca$ratio <- pJcostlca$solution_cost/pJcostinsider$solution_cost
pJcost <- aggregate(pJcostlca, by=list(pJcostlca$pj), FUN=mean)
pJcostsds <- aggregate(pJcostlca, by=list(pJcostlca$pj), FUN=sd)

pdf("figures/pJ-cost-insidervlca.pdf")
ggplot(data=pJcost) + geom_point(aes(pj,ratio)) +
	geom_errorbar(aes(pj,ymin=ratio-2*pJcostsds$ratio/sqrt(reps),ymax=ratio+2*pJcostsds$ratio/sqrt(reps)), width=0.02) + 
	scale_x_continuous(breaks=seq(0,100,by=0.2)) +
	xlab(bquote(p[J])) + ylab("Cost ratio")
dev.off()

pdf("figures/pJ-time.pdf")
ggplot(data=pJmeans) + geom_point(aes(pj,time, col=Group.2),position=position_dodge(width=0.05)) +
	geom_errorbar(aes(pj,ymin=pmax(0,time-2*pJsds$time/sqrt(reps)),ymax=time+2*pJsds$time/sqrt(reps),col=Group.2), width=0.1,position=position_dodge(width=0.05)) + 
	scale_x_continuous(breaks=seq(0,100,by=0.2)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual(values=c(3,1,2,4), labels=c("FastMultRec","inSiDeR","LCA","segdup")) +
	theme(legend.title=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(p[J])) + ylab("Time (s)")
dev.off()
