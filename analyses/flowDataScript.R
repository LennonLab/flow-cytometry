# Setup
setwd("~/Github/flow-cytometry/")

# Load source code and dependances
## for installation of bioconductor packages
source("http://bioconductor.org/biocLite.R")
source("./bin/flowPrep.R") 
source("./bin/support_functions.R")

# Install bioconductor packages. This process will take a bit of time, so be patient.
#biocLite(c("flowPeaks","flowCore","flowStats","flowViz",
#           "flowClust","flowQ", "flowUtils","flowMeans","flowDensity"))
#biocLite("GEOmap");biocLite("Logicle")
#biocLite("flowQ")
#biocLite("S4Vectors")
#biocLite("quantreg")

# Load bioconductor and stats packages
library(S4Vectors)
library(quantreg)
library(flowPeaks)  
library(flowCore)   
library(flowStats)  
library(flowViz)    
library(flowQ)      
library(flowClust)  
library(flowUtils)  
library(flowMeans) 
library(flowDensity)
library(GEOmap)     
#library(Logicle)
library(ggplot2)

##########################################
# File prep

# add the name of the batch to be processed
file.name <- "batch.name"

control.path <- paste("./data/",file.name,"-controls/", sep = "")
batch.path <- paste("./data/",file.name,"-batch/", sep = "")
  
# Read in control data
fs.controls <- read.flowSet(path = control.path, 
                            pattern = ".fcs", transformation = FALSE,
                            #alter.names = TRUE, column.pattern = ".A")
                            column.pattern = "-A")
fs.controls

# Control Compensation
fs.controls.comp <- fsApply(fs.controls,function(frame){
  #extract compensation matrix from keywords
  comp <- keyword(frame)$`SPILL`
  new_frame <- compensate(frame,comp)
  new_frame
})


# Control Transformation
chnls <- colnames(fs.controls.comp)[1:7]
lgcl <- estimateLogicle(samp, channels = chnls)
fs.trans <- lgcl %on% fs.controls.comp

# Control Gating
samp.beads <- fs.trans[[6]]
bead.chnl <- c("FSC-A","Alexa Fluor 488-A")
bead.cols <- c(1,4)
beads.plot <- plotDens(samp.beads, bead.chnl, devn = FALSE,
                       xlab = "FSC-A", ylab = "Alexa Fluor 488-A", las = 1)
bead.gate <- rectangleGate(filterId = "beads", 
                           "FSC-A" = c(3.0,3.2),
                           "Alexa Fluor 488-A" = c(0,1.25))
# Live/Dead gating from eFluor 660 - based on negative staining
## Must use unstained control and check with stained control
control.ld <- fs.trans[[5]]
test.ld <- fs.trans[[3]]
ld.chnl <- c("APC-A","SSC-A")
ld.cols <- c(2,7)

plot.new()

plotDens(control.ld, ld.chnl, devn = FALSE,
         xlab = "APC-A", ylab = "SSC-A", las = 1)

plotDens(test.ld, ld.chnl, devn = FALSE,
         xlab = "APC-A", ylab = "SSC-A", las = 1)

live.gate <- rectangleGate(filterId = "live",
                           "APC-A" = c(0,max(ld.cols[1])),
                           "SSC-A" = c(1.15,4.5))

dead.gate <- rectangleGate(filterId = "dead",
                           "APC-A" = c(max(ld.cols[[1]]),Inf))

# Subset flow data for control
actdorm <- Subset(fs.trans, live.gate)

results <- matrix(NA, nrow = length(sampleNames(fs.trans)), ncol = 9)
results <- as.data.frame(results)
colnames(results) <- c("sample","NA","ratio.min","ratio.max",
                       "live.dens","act.dens","dorm.dens",
                       "act.perc","dorm.perc")

for(i in 1:length(sampleNames(fs.trans))){
  DNA <- exprs(fs.trans[[i]][,"Pacific Blue-A"])
  RNA <- exprs(fs.trans[[i]][,"PI (B)-A"])
  RDratio <- RNA/DNA
  
  dat <- data.frame(DNA,RNA,RDratio)
  
  #make plot
  p <- ggplot(dat,aes(x = dat[,3])) + geom_density()
  print(p)
  
  range(RDratio)
  min <- 1-sd(RDratio)
  max <- 1+sd(RDratio)
  live.pop.dens <- length(RDratio[RDratio > min])
  dorm.pop.dens <- length(RDratio[RDratio > min & RDratio < max])
  act.pop.dens <- live.pop.dens-dorm.pop.dens
  per.act <- (live.pop.dens - dorm.pop.dens)/live.pop.dens 
  per.dorm <- (act.pop.dens)/live.pop.dens 
  
  results[i,1] <- sampleNames(fs.trans)[[i]]
  results[i,2] <- NA
  results[i,3] <- min
  results[i,4] <- max
  results[i,5] <- live.pop.dens
  results[i,6] <- act.pop.dens
  results[i,7] <- dorm.pop.dens
  results[i,8] <- round(per.act*100, digits = 3)
  results[i,9] <- round(per.dorm*100, digits = 3)
}

results


# Batch script
fs1 <- read.flowSet(path = batch.path, 
                    pattern = ".fcs", transformation = FALSE,
                    #alter.names = TRUE, column.pattern = ".A")
                    column.pattern = "-A")
fs1

fs1[[1]]@description$'SPILL'

fs1.comp <- fsApply(fs1,function(frame){
  #extract compensation matrix from keywords
  comp <- keyword(frame)$`SPILL`
  new_frame <- compensate(frame,comp)
  new_frame
})

fs1.trans <- lgcl %on% fs1.comp

# Rectangle gating results
bead.results <- filter(fs1.trans, bead.gate)
live.results <- filter(fs1.trans, live.gate)
dead.results <- filter(fs1.trans, dead.gate)

summary(bead.results)
summary(live.results)
summary(dead.results)

actdorm <- Subset(fs1.trans, live.gate)
sampleNames(actdorm)

# Gather results 
results.batch <- matrix(NA, nrow = length(sampleNames(actdorm)), ncol = 8)
results.batch <- as.data.frame(results.batch)
colnames(results.batch) <- c("sample","NA","ratio.min","ratio.max",
                        "live.dens","act.dens","dorm.dens",
                        "act.perc","dorm.perc")

# later iterations will need to account for bead count by sample
bead.count <- 10000

for(i in 1:length(sampleNames(actdorm))){
  DNA <- exprs(actdorm[[i]][,"Pacific Blue-A"])
  RNA <- exprs(actdorm[[i]][,"PI (B)-A"])
  RDratio <- RNA/DNA
  
  dat <- data.frame(DNA,RNA,RDratio)
  dat
  #make plot
  p <- ggplot(dat,aes(x = dat[,3]), main = i) + geom_density()
  print(p)
  
  range(RDratio)
  min <- 1-sd(RDratio)
  max <- 1+sd(RDratio)
  
  #densities must be divided by bead counts (10000)[# bacteria per 10^-6 mL of sample]
  live.pop.dens <- (length(RDratio[RDratio > min])/bead.count)*1000000
  dorm.pop.dens <- (length(RDratio[RDratio > min & RDratio < max])/bead.count)*1000000
  act.pop.dens <- (live.pop.dens-dorm.pop.dens)
  per.act <- (live.pop.dens - dorm.pop.dens)/live.pop.dens 
  per.dorm <- (dorm.pop.dens)/live.pop.dens 
  
  results.batch[i,1] <- sampleNames(actdorm)[[i]]
  results.batch[i,2] <- min
  results.batch[i,3] <- max
  results.batch[i,4] <- live.pop.dens
  results.batch[i,5] <- act.pop.dens
  results.batch[i,6] <- dorm.pop.dens
  results.batch[i,7] <- round(per.act*100, digits = 3)
  results.batch[i,8] <- round(per.dorm*100, digits = 3)
}

results.batch

# Export data
write.csv(file = paste("./data/",file.name,"PopDat.csv"), results.batch)