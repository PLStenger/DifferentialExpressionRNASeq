# DeSeq2 give you different table.
# Here an example with invented datas :
# We want to know the differential expression of genes between green, yellow and red inner shell in the Pinctada margaritifera species (oyster pearl)
# DeSeq2 give you 3 tables :
# 1 - Yellow vs Green
# 2 - Red vs Green
# 3 - Red vs Yellow
# Kept these only 3 columns : Name, Log2 (fold change) and the FDR (P-value adjusted)
# You must compare these tables against them like that :
# 1 - Red versus Green VS Red versus Yellow
# 2 - Yellow versus Green VS Red versus Green
# 3 - Yellow versus Green VS Red versus Yellow
# The fisrt colour is the normalisateur and it should be the same for the confrontation like 1 (here the normalisateur is the red and allows to compare the Green and Yellow colours.
# So, simply multiply the log2 column by -1 for those who are not in the right direction to have some thing like that :
# 1 - Red versus Green VS Red versus Yellow
# 2 - Green versus Yellow VS Green versus Red
# 3 - Yellow versus Green VS Yellow versus Red
# So you should have 6 tables. Convert in .txt :
# - 1 : Yellow vs green (yg.txt)
# - 2 : Green vs yellow (gy.txt)
# - 3 : Red vs yellow (ry.txt)
# - 4 : Yellow vs red (yr.txt)
# - 5 : Red vs green (rg.txt)
# - 6 : Red vs green (rg.txt)
# - 7 : Green vs red (gr.txt)

# Yellow vs green
yg <- read.csv("yg.txt", header = T, sep="\t")
yg$Log2 <- as.numeric(as.character(yg$Log2))
head(yg)

# Green vs yellow
gy <- read.csv("gy.txt", header = T, sep="\t")
gy$Log2 <- as.numeric(as.character(gy$Log2))
head(yg)

# Red vs yellow
ry <- read.csv("ry.txt", header = T, sep="\t")
ry$Log2 <- as.numeric(as.character(ry$Log2))
head(ry)

# Yellow vs red
yr <- read.csv("yr.txt", header = T, sep="\t")
yr$Log2 <- as.numeric(as.character(yr$Log2))
head(yr)

# Red vs green
rg <- read.csv("rg.txt", header = T, sep="\t")
rg$Log2 <- as.numeric(as.character(rg$Log2))
head(rg)

# Green vs red
gr <- read.csv("gr.txt", header = T, sep="\t")
gr$Log2 <- as.numeric(as.character(gr$Log2))
head(gr)

# Put a threshold for your FDR
THRESHOLD <- 0.05
syg <- subset(yg, FDR<THRESHOLD) ## result was generated from DESeq2 Galaxy 
sgy <- subset(gy, FDR<THRESHOLD) ## result was generated from DESeq2 Galaxy  
sry <- subset(ry, FDR<THRESHOLD) ## result was generated from DESeq2 Galaxy  
syr <- subset(yr, FDR<THRESHOLD) ## result was generated from DESeq2 Galaxy  
sgr <- subset(gr, FDR<THRESHOLD) ## result was generated from DESeq2 Galaxy  
srg <- subset(rg, FDR<THRESHOLD) ## result was generated from DESeq2 Galaxy  

################################################################
########################### syg vs syr #########################
################################################################

# how many genes are differentially expressed in Yellow vs Green ?
length(syg$Name)
# 84

# how many genes are differentially expressed in Yellow vs Red ?
length(syr$Name)
# 64

# How many genes are differentially expressed in common ?
data.frame1 <- data.frame(syg)
colnames(data.frame1) <- c("Name", "Log2a", "FDRa")
data.frame2 <- data.frame(syr)
colnames(data.frame2) <- c("Name", "Log2b", "FDRb")
merge(data.frame1,data.frame2, by.x = "Name")
# 7

# How many genes are inducted and repressed in common ? (Inducted = TRUE, Repressed = FALSE)
dat <- merge(data.frame1,data.frame2, by.x = "Name")
summary(dat$Log2a>0)
# Inducted = 3
# Repressed = 4

# How many genes are inducted and repressed in Yellow vs Green ? (Inducted = TRUE, Repressed = FALSE)
dat1 <- data.frame(data.frame1[!(data.frame1$Name %in% dat$Name),])
dat2 <- dat1[complete.cases(dat1), ]
summary(dat2$Log2a>0)
# Inducted = 50
# Repressed = 27

# How many genes are inducted and repressed in Yellow vs Red ? (Inducted = TRUE, Repressed = FALSE)
dat1 <- data.frame(data.frame2[!(data.frame2$Name %in% dat$Name),])
dat2 <- dat1[complete.cases(dat1), ]
summary(dat2$Log2b>0)
# Inducted = 35
# Repressed = 22

################################################################
########################### srg vs sry #########################
################################################################

# how many genes are differentially expressed in Red vs Green ?
length(srg$Name)
# 72

# how many genes are differentially expressed in Red vs Yellow ?
length(sry$Name)
# 64

# How many genes are differentially expressed in common ?
data.frame1 <- data.frame(srg)
colnames(data.frame1) <- c("Name", "Log2a", "FDRa")
data.frame2 <- data.frame(sry)
colnames(data.frame2) <- c("Name", "Log2b", "FDRb")
merge(data.frame1,data.frame2, by.x = "Name")
# 24

# How many genes are inducted and repressed in common ? (Inducted = TRUE, Repressed = FALSE)
dat <- merge(data.frame1,data.frame2, by.x = "Name")
summary(dat$Log2a>0)
# Inducted = 6
# Repressed = 18

# How many genes are inducted and repressed in Red vs Green ? (Inducted = TRUE, Repressed = FALSE)
dat1 <- data.frame(data.frame1[!(data.frame1$Name %in% dat$Name),])
dat2 <- dat1[complete.cases(dat1), ]
summary(dat2$Log2a>0)
# Inducted = 29
# Repressed = 19

# How many genes are inducted and repressed in Red vs Yellow ? (Inducted = TRUE, Repressed = FALSE)
dat1 <- data.frame(data.frame2[!(data.frame2$Name %in% dat$Name),])
dat2 <- dat1[complete.cases(dat1), ]
summary(dat2$Log2b>0)
# Inducted = 20
# Repressed = 20


################################################################
########################### sgy vs sgr #########################
################################################################

# how many genes are differentially expressed in Green vs Yellow ?
length(sgy$Name)
# 84

# how many genes are differentially expressed in Green  vs Red ?
length(sgr$Name)
# 72

# How many genes are differentially expressed in common ?
data.frame1 <- data.frame(sgy)
colnames(data.frame1) <- c("Name", "Log2a", "FDRa")
data.frame2 <- data.frame(sgr)
colnames(data.frame2) <- c("Name", "Log2b", "FDRb")
merge(data.frame1,data.frame2, by.x = "Name")
# 19

# How many genes are inducted and repressed in common ? (Inducted = TRUE, Repressed = FALSE)
dat <- merge(data.frame1,data.frame2, by.x = "Name")
summary(dat$Log2a>0)
# Inducted = 9
# Repressed = 10

# How many genes are inducted and repressed in Green vs Yellow ? (Inducted = TRUE, Repressed = FALSE)
dat1 <- data.frame(data.frame1[!(data.frame1$Name %in% dat$Name),])
dat2 <- dat1[complete.cases(dat1), ]
summary(dat2$Log2a>0)
# Inducted = 22
# Repressed = 43

# How many genes are inducted and repressed in Green vs red ? (Inducted = TRUE, Repressed = FALSE)
dat1 <- data.frame(data.frame2[!(data.frame2$Name %in% dat$Name),])
dat2 <- dat1[complete.cases(dat1), ]
summary(dat2$Log2b>0)
# Inducted = 28
# Repressed = 25



################################################################################################################################
############################## VennDiagram with the "VennDiagram" package (just to check) ######################################
################################################################################################################################


# install.packages("VennDiagram")


# yg vs yr
library(VennDiagram)
pdf("venn_diagram_yg_vs_yr.pdf")
venn.plot <- venn.diagram(list(syg$Name, syr$Name), NULL, fill=c("chartreuse", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Yellow vs green", "Yellow vs red"))
grid.draw(venn.plot)
dev.off()

# ry vs rg
library(VennDiagram)
pdf("venn_diagram_ry_vs_rg.pdf")
venn.plot <- venn.diagram(list(sry$Name, srg$Name), NULL, fill=c("yellow", "chartreuse"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Red vs yellow", "Red vs green"))
grid.draw(venn.plot)
dev.off()

# gy vs gr
library(VennDiagram)
pdf("venn_diagram_gy_vs_gr.pdf")
venn.plot <- venn.diagram(list(sgy$Name, sgr$Name), NULL, fill=c("yellow", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Green vs yellow", "Green vs red"))
grid.draw(venn.plot)
dev.off()


venn.plot1 <- venn.diagram(list(syg$Name, syr$Name), NULL, fill=c("chartreuse", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Yellow vs green", "Yellow vs red"))
venn.plot2 <- venn.diagram(list(sry$Name, srg$Name), NULL, fill=c("yellow", "chartreuse"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Red vs yellow", "Red vs green"))
venn.plot3 <- venn.diagram(list(sgy$Name, sgr$Name), NULL, fill=c("yellow", "red"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("Green vs yellow", "Green vs red"))

grid.draw(venn.plot1)
grid.draw(venn.plot2)
grid.draw(venn.plot3)



################################################################################################################################
###################### VennDiagram with the "venneuler" package (to product more aesthetic graph) ##############################
################################################################################################################################


# install.packages("venneuler")
library(venneuler)


# YGvsYR
YGvsYR <- venneuler(c(YG =84, YR =64, "YG&YR"=7))
YGvsYR$labels <- c(" "," ")

class(YGvsYR) 

YGvsYR$colors<- c(0.31, 1) 
plot(YGvsYR)

text(0.5,0.5,"7")
text(0.2,0.5,"57")
text(0.8,0.5,"77")
text(0.2,0.9,"Yellow vs Red")
text(0.8,0.9,"Yellow vs Green")

# RGvsRY
RGvsRY <- venneuler(c(RG =72, RY =64, "RG&RY"=24))
RGvsRY$labels <- c(" "," ")

class(RGvsRY) 

RGvsRY$colors<- c(0.31, 0.2) 
plot(RGvsRY)

text(0.5,0.5,"24")
text(0.2,0.5,"40")
text(0.8,0.5,"48")
text(0.2,0.9,"Red vs Yellow")
text(0.8,0.9,"Red vs Green")


# GYvsGR
GYvsGR <- venneuler(c(GY =84, GR =72, "GY&GR"=19))
GYvsGR$labels <- c(" "," ")

class(GYvsGR) 

GYvsGR$colors<- c(0.2, 1) 
plot(GYvsGR)

text(0.5,0.5,"19")
text(0.2,0.5,"53")
text(0.8,0.5,"65")
text(0.2,0.9,"Green vs Red")
text(0.8,0.9,"Green vs Yellow")

