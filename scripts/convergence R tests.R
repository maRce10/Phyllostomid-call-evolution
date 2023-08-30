
# install.packages("convevol")
library(convevol)

phyl<-rtree(10)

phendata <- fastBM(phyl,nsim=2)


plot(phendata, col = "white")

text(phendata, rownames(phendata))


convtips <- c("t1","t2","t3")

convtips <- c("t2","t4","t5")


phendata[4,1] <- -3
phendata[1,1] <- -2

convrat(phyl,phendata,convtips)


convratsig(phyl,phendata,convtips, 100)$Pvals

par(mfrow = c(1, 2))

plot(phyl)

cols<-c(ifelse(rownames(phendata) %in% convtips, "red", "black"), rep("black",phyl$Nnode))

names(cols)<-1:(length(phyl$tip) + phyl$Nnode)

phylomorphospace(phyl, phendata, control=list(col.node=cols))


# install.packages("windex")

library(windex)



data(sample.data)

data(sample.tree)

windex(dat = sample.data, 
       tree = sample.tree, 
       traits = c("ou1", "ou2"), 
       focal = sample.data[ , 2], 
       SE = FALSE)


test.windex(dat = sample.data, 
       tree = sample.tree, 
       traits = c("ou1", "ou2"), 
       focal = sample.data[ , 2], 
       reps = 30)


test.windex(dat = sample.data, 
            tree = sample.tree, 
            traits = c("bm1", "bm2"), 
            focal = sample.data[ , 2], 
            reps = 30)


mat <- sample.data[, c("bm2", "bm1")]

rownames(mat) <- sample.data$species

cols <- c(ifelse(sample.data$focals == 1, "red", "black"), rep("black", sample.tree$Nnode))

names(cols) <- 1:(length(sample.tree$tip) + sample.tree$Nnode)


phylomorphospace(sample.tree, mat, control=list(col.node=cols))

