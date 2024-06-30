#function to computer posterior probabilities at the nodes (L38-47) by liam revell from http://blog.phytools.org/2012/12/plotting-node-piecharts-on-top-of.html
#function to plot simmaps and node posterior probabilities (L50-59) by Dr. Michael May (UC Berkeley)

library(phytools)
library(grDevices)
setwd("~/Desktop/SapindalesWoodAnatomy-main/")
tree<-read.nexus("Tree.tre")

# upload data 
data  <- read.csv("DataTable.csv", row.names = 1) 
datum <- as.data.frame(cbind(rownames(data), data[,2]), row.names = FALSE) #select the column with the trait of interest
datum <- datum[complete.cases(datum), ]

# species
species <- datum$V1

# reformat data
datas <- datum$V2
names(datas) <- species

# drop tips
pruned_tree <- keep.tip(tree, species)

#DISCRETE TRAITS
er  <- fitMk(pruned_tree, datas, model="ER",  ordered=TRUE, pi="estimated")
ard <- fitMk(pruned_tree, datas, model="ARD", ordered=TRUE, pi="estimated")

# chi-square test
1 - pchisq(2*abs(ard$logLik- er$logLik), 1) # if p is not < .05 then select ER model

# stochastic mapping
simmap.trees <- make.simmap(pruned_tree, datas, model = "ER", pi = "estimated", nsim=1000)

# make colors
col_vec <- c("dodgerblue4", "orangered3")
states<- c("absent","present") # make colors

# function to compute the node states
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

XX<-sapply(simmap.trees,foo)
pies<-t(apply(XX,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=1000))

#generate summary of stochastic maps with pies of posterior at nodes..
source("plot_simmap.R")

plot_simmap(time_tree = simmap.trees[[1]], 
            tree = simmap.trees[[1]], 
            simmaps = simmap.trees, 
            states = states,
            show.tip.label = T,
            lwd = 4,
            label.cex = .6,
            colors = col_vec, edge.width=0, nt=10001)

title(main="Ray Composition", adj = .0004, line= -.01)
add.simmap.legend(colors=col_vec, prompt=FALSE,x=105,y=180, fsize =.9)
nodelabels(pie=pies,cex=0.36,piecol=col_vec, lwd=1)
legend(x = -7, y = 150, legend = states, col = col_vec, pch = 20, yjust = 0, bty = 'n', cex =.7, pt.cex = 4)
