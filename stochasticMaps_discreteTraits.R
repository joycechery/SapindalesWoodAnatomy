library(phytools)
library(grDevices)
setwd("~/Desktop/")
tree<-read.nexus("Tree.tre")

# upload data 
data  <- read.csv("DataTable.csv", row.names = 1) 
datum <- as.data.frame(cbind(rownames(data), data[,13]), row.names = FALSE) #select the column with the trait of interest
datum <- datum[complete.cases(datum), ]

# reformat data
species <- datum$V1
datas <- datum$V2
names(datas) <- species

# drop tips
pruned_tree <- keep.tip(tree, species)

# model fitting
er  <- fitMk(pruned_tree, datas, model="ER",  pi="estimated")
ard <- fitMk(pruned_tree, datas, model="ARD",  pi="estimated")

# chi-square test
1 - pchisq(2*abs(ard$logLik- er$logLik), 1) # if p is not < .05 then select ER model

# stochastic mapping
simmap.trees <- make.simmap(pruned_tree, datas, model = "ARD", pi = "estimated", nsim=100)

# make colors
colrs <- c("dodgerblue4", "orangered3",  "olivedrab4", "grey")
cols  <- colrs[1:length(ard$states)]
names(cols) <- ard$states

# plot the maps
densityTree(simmap.trees,
            method = "plotSimmap",
            lwd = 3,
            nodes = "intermediate",
            colors = cols,
            compute.consensus=FALSE,
            fsize=.50, show.axis = F, direction="lefttwards", alpha = 0.5)

title(main="Growth Rings", adj = .0004, line= -.01)
add.simmap.legend(colors=cols, prompt=FALSE,x=109,y=170, fsize =.9)