dat <- read.csv("stats.csv", header=T)

insider <- dat[dat$method=="insider",]
fmr <- dat[dat$method=="fastmultrec_greedy",]
lca <- dat[dat$method=="lcamap",]

#cost
par(mfrow=c(1,2))
plot(insider$solution_cost ~ fmr$solution_cost, main="inSiDeR vs FastMultRec: cost")
curve(1*x,add=T,col="red")
plot(insider$solution_cost ~ lca$solution_cost, main="inSiDeR vs LCA: cost")
curve(1*x,add=T,col="red")

#time
par(mfrow=c(1,2))
boxplot(time ~ factor(method) + factor(dup_cost), data=dat, ylim=c(0,3), main = "Time for dup cost")
boxplot(time ~ factor(method) + factor(duprate), data=dat, ylim=c(0,3), main = "Time for dup rate")
