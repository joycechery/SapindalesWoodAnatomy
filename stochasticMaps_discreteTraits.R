#function to computer posterior probabilities at the nodes (L44-52) by liam revell from http://blog.phytools.org/2012/12/plotting-node-piecharts-on-top-of.html

library(phytools)
library(grDevices)
setwd("~/Desktop/")
tree<-read.nexus("Tree.tre")

# upload data 
data  <- read.csv("DataTable.csv", row.names = 1) 
datum <- as.data.frame(cbind(rownames(data), data[,12]), row.names = FALSE) #select the column with the trait of interest
datum <- datum[complete.cases(datum), ]

# species
species <- datum$V1

# reformat data
datas <- datum$V2
names(datas) <- species

# drop tips
pruned_tree <- keep.tip(tree, species)

#POLYMORPHIC TRAITS
er  <- fitpolyMk(pruned_tree, datas, model="ER",  ordered=TRUE, pi="estimated")
ard <- fitpolyMk(pruned_tree, datas, model="ARD", ordered=TRUE, pi="estimated")

# chi-square test
er_num_params  <- 1
ard_num_params <- sum(ard$index.matrix != 0)
1 - pchisq(2*abs(ard$logLik- er$logLik), ard_num_params - er_num_params) # if p is not < .05 then selecr ER model

ard$index.matrix
ard$data

# stochastic mapping
simmap.trees <- make.simmap(pruned_tree, ard$data, model = ard$index.matrix, pi = "estimated", nsim=1000)

# make colors
col_vec <- c("dodgerblue4", "orangered3", "olivedrab4")
states<- c("hetero","hetero+homo", "homo")

# function to compute the states
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

XX<-sapply(simmap.trees,foo)
pies<-t(apply(XX,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=1000))

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
