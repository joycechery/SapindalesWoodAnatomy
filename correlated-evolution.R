###correlation evolution

setwd("~/Desktop/")
tree<-read.nexus("Tree.tre")

# upload data
data  <- read.csv("IntervesselPitsxTyloses_discrete.csv", row.names = 1)

# reformat character 1
datum <- as.data.frame(cbind(rownames(data), data[,1]), row.names = FALSE)
datum <- datum[complete.cases(datum), ]
species <- datum$V1
char1 <- datum$V2
names(char1) <- species

# reformat character 2  
datum <- as.data.frame(cbind(rownames(data), data[,2]), row.names = FALSE)
datum <- datum[complete.cases(datum), ]
species <- datum$V1
char2 <- datum$V2
names(char2) <- species

#pagels 1994 test of correlated evolution
ardPAGEL<-fitPagel(pruned_tree, x=char1, y=char2, model = "ARD", dep.var = "y") # ard model
plot(ardPAGEL)

#erPAGEL<-fitPagel(pruned_tree, x=char1, y=char2, model = "ER", dep.var = "y")  #er model



