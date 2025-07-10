###
#Summer school 2025
#
##
install.packages("promises")
install.packages("adegenet")

#If you have error:
#Error in install.packages : ERROR: failed to lock directory ‘C:\Users\TFerme\AppData\Local\R\win-library\4.4’ for modifying
#Try removing ‘C:\Users\TFerme\AppData\Local\R\win-library\4.4/00LOCK’
#Remove 00LOCK folder and install packages again



library("ggplot2")
library("adegenet")
library("RColorBrewer")
library("ggpubr")

#Set working directory-write your path to data
setwd("/")

#Upload data, write number of individuals, number of loci
samples_str <- read.structure("Mikrosateliti.str", n.ind=000, n.loc=00, col.lab=1, col.pop=2, col.others = NULL, row.marknames = 1, NA.char = "-9", ask = TRUE, quiet = FALSE)

pca <- dudi.pca(samples_str, nf = 10, scannf = FALSE)

cols <- brewer.pal(n = nPop(samples_str), name = "Dark2")

set.seed(9)
pca.scores <- as.data.frame(pca$scores)

p1 <- ggplot(pca$li, aes(x=Axis1, y=Axis2, colour=samples_str@pop)) 
p1 <- p1 + geom_point(size=2, alpha=0.5)
p1 <- p1 + stat_ellipse(level = 0.95, size = 1)
p1 <- p1 + scale_color_manual(values = cols) 
p1 <- p1 + geom_hline(yintercept = 0) 
p1 <- p1 + geom_vline(xintercept = 0) 
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("PC1") + ylab("PC2")
p1 <- p1 + theme(legend.position = "right")
p1


p2 <- ggplot(pca$li, aes(x=Axis3, y=Axis4, colour=samples_str@pop))
p2 <- p2 + geom_point(size=2, alpha=0.5)
p2 <- p2 + stat_ellipse(level = 0.95, size = 1)
p2 <- p2 + scale_color_manual(values = cols) 
p2 <- p2 + geom_hline(yintercept = 0) 
p2 <- p2 + geom_vline(xintercept = 0) 
p2 <- p2 + theme_bw()
p2 <- p2 + xlab("PC3") + ylab("PC4")
p2 <- p2 + theme(legend.position = "right")
p2

str_pca <- ggarrange(p1, p2, ncol = 2, nrow = 1)

str_pca
