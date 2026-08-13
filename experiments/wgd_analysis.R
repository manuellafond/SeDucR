library(ggplot2)
library(scales)
library(reshape2)
theme_set(theme_bw(base_size=20))

#helper function for axis breaks
lb = function(maj, by=1, lower=0.1, radix=10) {
  function(x) {
    minx         = floor(max(-15,min(logb(x,radix), na.rm=T))) - 1
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

dat <- read.csv("stats_wgd_all.csv", header=T)

attach(dat)

### 1. insider only vs d, colour = duprate

subdat <- subset(dat, method=="insider")
attach(subdat)

means <- aggregate(time ~ duprate + dup_cost, data = subdat, FUN=mean)
means$duprate <- factor(means$duprate)
counts <- aggregate(time ~ duprate + dup_cost, data = subdat, FUN=length)
sds <- aggregate(time ~ duprate + dup_cost, data = subdat, FUN=sd)


pdf("figures/wgd-d.pdf")
ggplot(data=means) + geom_point(aes(dup_cost, time, col=duprate)) +
	geom_errorbar(aes(dup_cost,ymin=time-2*sds$time/sqrt(counts$time),ymax=time+2*sds$time/sqrt(counts$time), col=duprate), width=0.5) + 
	scale_x_continuous(breaks=seq(0,100,by=5)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual("Duplication/loss rate", values=1:5, labels=c(bquote(10^-8),bquote(10^-9),bquote(10^-10),bquote(10^-11),bquote(10^-15))) +
	theme(legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab(bquote(delta)) + ylab("Time (s)")
dev.off()


### 2. all methods vs duprate, d = 15, colour = method

subdat <- subset(dat, dup_cost == 15)
subdat$duprate <- 10^(-subdat$duprate)
attach(subdat)

means <- aggregate(time ~ duprate + method, data = subdat, FUN=mean)
counts <- aggregate(time ~ duprate + method, data = subdat, FUN=length)
sds <- aggregate(time ~ duprate + method, data = subdat, FUN=sd)

pdf("figures/wgd-duprate.pdf")
ggplot(data=means) + geom_point(aes(duprate, time, col=method)) +
	geom_errorbar(aes(duprate,ymin=time-2*sds$time/sqrt(counts$time),ymax=time+2*sds$time/sqrt(counts$time), col=method), width=0.2) + 
	scale_x_log10(minor_breaks=lb(TRUE)) + scale_y_log10(minor_breaks=lb(FALSE,by=2)) +
	scale_colour_manual(values=c(3,1,2,4), labels=c("FastMultRec","inSiDeR","LCA","segdup")) +
	theme(legend.position = c(.05, .85), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab("Duplication/loss rate") + ylab("Time (s)")
dev.off()


### 3. percentage of suboptimal instances

subopt <- read.csv("subopt.csv")
subopt$duprate <- 10^(-subopt$duprate)
subopt[,-1] <- subopt[,-1]/subopt$nbinstances
subopt <- melt(subopt, measure.vars = c("segdup","fastmultrec_greedy","lcamap"))

pdf("figures/wgd-subopt.pdf")
ggplot(data=subopt) + geom_point(aes(duprate,value, col=variable)) + 
	scale_x_log10(minor_breaks=lb(TRUE)) +
	ylim(c(0,1)) +
	scale_colour_manual(values=c(4,3,2), labels=c("segdup","FastMultRec","LCA")) +
	theme(legend.title=element_blank(),legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(6, 6, 6, 6)) +
	xlab("Duplication/loss rate") + ylab("Proportion of suboptimal instances")
dev.off()
