# Plot fvec file

Plot graphic for feature vectors


R colour palette ideas [here](https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf). 

```
ex <- read.table("fvec.example", header=T)
#we don't need the first four columns
ex2 <- as.data.frame(ex[5:ncol(ex)])

#rescale each column to be 0-1. We do this so that the different stats can all be plotted on the same colour scale
rescaled <- apply(ex2, MARGIN=2, FUN=function(X) (X-min(X))/diff(range(X)))

#transform to long data
ex.t <- as.data.frame(t(rescaled))

#Add a column with all window numbers
ex.t$windows <- rep(0:10, 12)

#Add a column with all stats named
ex.t$statsNames <- c(rep("pi",11), rep("thetaW", 11), rep("tajD", 11), rep("distVar", 11), rep("distSkew",11), rep
("distKurt", 11), rep("nDiplos", 11), rep("diplo_H1", 11), rep("diplo_H12", 11), rep("diplo_H2.H1", 11), rep("dipl
o_ZnS", 11), rep("diplo_Omega",11))

#Calculate the mean of all the individual simulations
ex.t$mean <- rowMeans(ex.t[1:9])


#plot

ggplot(ex.t, aes(windows, statsNames, fill=mean))+geom_tile()+scale_fill_continuous(type = "gradient")
```









