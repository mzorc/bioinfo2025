

################################################################################
############################   Example 1    ####################################
################################################################################


############   WORKFLOW   ############ 

###  1. Import of .ped/.map files after QC

###  2. ROH determination 

###  3. Import of Function.RData

###  4. Creation of ROH Incidence Matrix

###  5. Froh estimation

###  6. Regional Froh estimation

###  7. Visualisation of Froh results

###  8. Creation of Genomic Relationship Matrix (GRM)

###  9. GRM Visualization

######################################


### 1. Import of .ped/.map files after QC


## First, set the working directory where the datasets are stored:
# a) Manually: Session -> Set Working Directory -> Choose Directory...
# b) Using code:
# setwd("C:/Users/mario/Desktop/MATE2024 Population Genomics Workshop/Day 2/3_Estimation of Inbreeding and Relatedness/2_Practical excercises/Example1/Directory")


## Then we import data (using data.table package) 
install.packages('data.table')
library(data.table)
map_file <- as.data.frame(fread(file = 'LabradorRetriever_FinalReport_afterQC.map',header = FALSE))
ped_file <- as.data.frame(fread(file = 'LabradorRetriever_FinalReport_afterQC.ped',header = FALSE))


## NOTE: You must not have any special signs in your Sample ID and Family ID (such as (, -, / or .)


## We also import Path to data (neccesary for ROH determination)
genotypeFilePath <- 'LabradorRetriever_FinalReport_afterQC.ped'
mapFilePath <- 'LabradorRetriever_FinalReport_afterQC.map'
# Only the file name is required if the files are located in used Directory

# Ensure the file paths include file name at the end
# Note: When copying paths from Windows Explorer, replace every '\' with '/' 
# For example:
# Windows-style path: "C:\Users\mario\Desktop\example.ped"
# Correct R-style path: "C:/Users/mario/Desktop/example.ped"


################################################################################


###  2. ROH determination 


## Package detectRUNS will be used for ROH determination
install.packages("detectRUNS")
library(detectRUNS)

# Package detectRUNS - 2 methods for ROH determination: 1. sliding-window (Plink method)
#                                                       2. consecutive runs (better, without window, SNP by SNP)

# Use consecutive runs method - more precise -> consecutiveRUNS.run function

# Arguments: 1. minSNP (minimum number of homozygous/heterozygous SNP in the run)
#            2. ROHet (default = FALSE)
#            3. maxGap (maximum gap between consecutive SNPs, in bps)
#            4. minLenghtBps (minimum length of the run, in bps)
#            5. maxOppRun (maximum number of opposite genotypes in the run)
#            6. maxMissRun (maximum number of missing genotypes in the run)

# Suggestion for the 50K SNP arrays: minSNP = 15; ROHet = FALSE; maxGap = 1000000 (in bps); minLenghtBps = 2000000 or 4000000 (in bps);
#                                    maxOppRun = depends on you (my suggestion is 1 because you do not work with classes - SVS software needed); maxMissRun = same as before (put it at 1)


# Suggestion for the HD SNP arrays:  minSNP = 15; ROHet = FALSE; maxGap = 1000000 (in bps); minLenghtBps = 1000000 (in bps);
#                                    maxOppRun = depends on you (my suggestion is 1 because you do not work with classes - SVS software needed); maxMissRun = same as before (put it at 1)


ROHs <- consecutiveRUNS.run(genotypeFile = genotypeFilePath, mapFile = mapFilePath, minSNP = 15, ROHet = FALSE, 
                            maxGap = 1000000, minLengthBps = 1000000, maxOppRun = 1, maxMissRun = 1)

str(ROHs)
ROHs


################################################################################


###  3. Import of Function.RData


## Function.RData file contains 4 functions essential for subsequent steps in Froh estimation procedure
## Start.End.Index()
## Incidence.Matrix()
## Froh.estimation()
## Regional.Froh.estimation()

load('Function.RData')


################################################################################


###  4. Creation of ROH Incidence Matrix


## ROH incidence matrix is a binary matrix where rows represent SNPs, 
## columns represent individuals, and each cell indicates whether a specific SNP is 
## inside ROH in an individual (1 for presence, 0 for absence)

## Start.End.Index() and Incidence.Matrix() functions are used:

# Arguments: chr, from, to
#            id, Start.Index, End.Index

# chr -> Column name of chromosomes
# from -> Column name of ROH start
# to ->  -> Column name of ROH end

# id -> Column name of individuals
# Start.Index -> Column name of index of SNP corresponding to 'from' position
# End.Index -> Column name of index of SNP corresponding to 'to' position

# First add Start.Index and End.Index columns -> Start.End.Index() function
Index_test <- Start.End.Index(map = map_file, ROH_file = ROHs, chr = "chrom", from = 'from',to = 'to')

I_matrix_test <- Incidence.Matrix(map = map_file, ROH_file = Index_test, id = "id", Start.Index = "Start.Index", End.Index = "End.Index")


################################################################################


###  5. Froh estimation


## Froh is estimated as proportion of genome information used covered by ROH


## Values = from 0 to 1


## Froh.estimation() function is used:

# Arguments: return, by.fid and by.chr

# Options:
# return = 1  - Returns Froh values for all individuals
# return = 2  - Returns mean Froh values

# by.fid = FALSE            - All individuals are treated as one group
# by.fid = TRUE             - You get results for all breeds separately

# by.chr = FALSE            - All Chromosomes are treated as one group  (Use it also if your dataset contains only one Chromosome)
# by.chr = TRUE             - You get results for all chromosomes separately


## 6 combinations:  a) return = 1, by.chr = FALSE,                 --- returns Froh of all individuals
#    (results)      b) return = 1, by.chr = TRUE,                  --- returns Froh of all individuals by Chr

#                   c) return = 2, by.chr = FALSE, by.fid = FALSE, --- returns mean Froh
#                   d) return = 2, by.chr = FALSE, by.fid = TRUE,  --- returns mean Froh by breed
#                   e) return = 2, by.chr = TRUE,  by.fid = FALSE, --- returns mean Froh by Chr
#                   f) return = 2, by.chr = TRUE,  by.fid = TRUE   --- returns mean Froh by breed by Chr


# Froh of all individuals
Froh_test_allId <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 1, by.chr = FALSE)

# Froh of all individuals by Chromosome
Froh_test__allId_byChr <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 1, by.chr = TRUE)

# Mean Froh
Froh_test_mean <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 2,  by.chr = FALSE)

# Mean Froh by Family ID
Froh_test_mean_byFid <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 2, by.fid = TRUE, by.chr = FALSE)

# Mean Froh by Chromosome
Froh_test_meanbyChr <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 2, by.chr = TRUE)

# Mean Froh by Family ID inside each Chromosome
Froh_test_meanbyChr_byFid <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 2, by.fid = TRUE, by.chr = TRUE)


###############


### Now we will classify inbreeding by its age -> split initial ROHs object
###                                               into length categories (1-2 Mb, 2-4 Mb, 4-8 Mb, 8-16 Mb, >16 Mb)
###                                               and repeat the previous procedure


## ROH from 1 to 2 Mb
ROHs_1to2Mb <- subset(ROHs, ROHs$lengthBps < 2000000)

## ROH from 2 to 4 Mb
ROHs_2to4Mb <- subset(ROHs, ROHs$lengthBps >= 2000000 & ROHs$lengthBps < 4000000)

## ROH from 4 to 8 Mb
ROHs_4to8Mb <- subset(ROHs, ROHs$lengthBps >= 4000000 & ROHs$lengthBps < 8000000)

## ROH from 8 to 16 Mb
ROHs_8to16Mb <- subset(ROHs, ROHs$lengthBps >= 8000000 & ROHs$lengthBps < 16000000)

## ROH > 16 Mb
ROHs_from16Mb <- subset(ROHs, ROHs$lengthBps >= 16000000)


###############


###  5.a) Froh estimation from 1 to 2 Mb (from approximately 50 to 25 generations ago)


## Add Start and End Index columns
Index_test_1to2Mb <- Start.End.Index(map = map_file, ROH_file = ROHs_1to2Mb, chr = "chrom", from = 'from',to = 'to')

## ROH Incidence matrix creation
I_matrix_test_1to2Mb <- Incidence.Matrix(map = map_file, ROH_file = Index_test_1to2Mb, id = "id", Start.Index = "Start.Index", End.Index = "End.Index")

## Froh of all individuals
Froh_test_allId_1to2Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_1to2Mb, return = 1, by.chr = FALSE)

## Froh of all individuals by Chromosome
# Froh_test__allId_byChr_1to2Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_1to2Mb, return = 1, by.chr = TRUE)

## Mean Froh
Froh_test_mean_1to2Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_1to2Mb, return = 2,  by.chr = FALSE)

## Mean Froh by Family ID
# Froh_test_mean_byFid_1to2Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_1to2Mb, return = 2, by.fid = TRUE, by.chr = FALSE)

## Mean Froh by Chromosome
# Froh_test_meanbyChr_1to2Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_1to2Mb, return = 2, by.chr = TRUE)

## Mean Froh by Family ID inside each Chromosome
# Froh_test_meanbyChr_byFid_1to2Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_1to2Mb, return = 2, by.fid = TRUE, by.chr = TRUE)


###############


###  5.b) Froh estimation from 2 to 4 Mb (from approximately 25 to 12.5 generations ago)


## Add Start and End Index columns
Index_test_2to4Mb <- Start.End.Index(map = map_file, ROH_file = ROHs_2to4Mb, chr = "chrom", from = 'from',to = 'to')

## ROH Incidence matrix creation
I_matrix_test_2to4Mb <- Incidence.Matrix(map = map_file, ROH_file = Index_test_2to4Mb, id = "id", Start.Index = "Start.Index", End.Index = "End.Index")

## Froh of all individuals
Froh_test_allId_2to4Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_2to4Mb, return = 1, by.chr = FALSE)

## Froh of all individuals by Chromosome
# Froh_test__allId_byChr_2to4Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_2to4Mb, return = 1, by.chr = TRUE)

## Mean Froh
Froh_test_mean_2to4Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_2to4Mb, return = 2,  by.chr = FALSE)

## Mean Froh by Family ID
# Froh_test_mean_byFid_2to4Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_2to4Mb, return = 2, by.fid = TRUE, by.chr = FALSE)

## Mean Froh by Chromosome
# Froh_test_meanbyChr_2to4Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_2to4Mb, return = 2, by.chr = TRUE)

## Mean Froh by Family ID inside each Chromosome
# Froh_test_meanbyChr_byFid_2to4Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_2to4Mb, return = 2, by.fid = TRUE, by.chr = TRUE)


###############


###  5.c) Froh estimation from 4 to 8 Mb (from approximately 12.5 to 6.25 generations ago)


## Add Start and End Index columns
Index_test_4to8Mb <- Start.End.Index(map = map_file, ROH_file = ROHs_4to8Mb, chr = "chrom", from = 'from',to = 'to')

## ROH Incidence matrix creation
I_matrix_test_4to8Mb <- Incidence.Matrix(map = map_file, ROH_file = Index_test_4to8Mb, id = "id", Start.Index = "Start.Index", End.Index = "End.Index")

## Froh of all individuals
Froh_test_allId_4to8Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_4to8Mb, return = 1, by.chr = FALSE)

## Froh of all individuals by Chromosome
# Froh_test__allId_byChr_4to8Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_4to8Mb, return = 1, by.chr = TRUE)

## Mean Froh
Froh_test_mean_4to8Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_4to8Mb, return = 2,  by.chr = FALSE)

## Mean Froh by Family ID
# Froh_test_mean_byFid_4to8Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_4to8Mb, return = 2, by.fid = TRUE, by.chr = FALSE)

## Mean Froh by Chromosome
# Froh_test_meanbyChr_4to8Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_4to8Mb, return = 2, by.chr = TRUE)

## Mean Froh by Family ID inside each Chromosome
# Froh_test_meanbyChr_byFid_4to8Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_4to8Mb, return = 2, by.fid = TRUE, by.chr = TRUE)


###############


###  5.d) Froh estimation from 8 to 16 Mb (from approximately 6.25 to 3.125 generations ago)


## Add Start and End Index columns
Index_test_8to16Mb <- Start.End.Index(map = map_file, ROH_file = ROHs_8to16Mb, chr = "chrom", from = 'from',to = 'to')

## ROH Incidence matrix creation
I_matrix_test_8to16Mb <- Incidence.Matrix(map = map_file, ROH_file = Index_test_8to16Mb, id = "id", Start.Index = "Start.Index", End.Index = "End.Index")

## Froh of all individuals
Froh_test_allId_8to16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_8to16Mb, return = 1, by.chr = FALSE)

## Froh of all individuals by Chromosome
# Froh_test__allId_byChr_8to16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_8to16Mb, return = 1, by.chr = TRUE)

## Mean Froh
Froh_test_mean_8to16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_8to16Mb, return = 2,  by.chr = FALSE)

## Mean Froh by Family ID
# Froh_test_mean_byFid_8to16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_8to16Mb, return = 2, by.fid = TRUE, by.chr = FALSE)

## Mean Froh by Chromosome
# Froh_test_meanbyChr_8to16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_8to16Mb, return = 2, by.chr = TRUE)

## Mean Froh by Family ID inside each Chromosome
# Froh_test_meanbyChr_byFid_8to16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_8to16Mb, return = 2, by.fid = TRUE, by.chr = TRUE)


###############


###  5.e) Froh estimation >16 Mb (from approximately 3.125 to 1 generations ago)


## Add Start and End Index columns
Index_test_from16Mb <- Start.End.Index(map = map_file, ROH_file = ROHs_from16Mb, chr = "chrom", from = 'from',to = 'to')

## ROH Incidence matrix creation
I_matrix_test_from16Mb <- Incidence.Matrix(map = map_file, ROH_file = Index_test_from16Mb, id = "id", Start.Index = "Start.Index", End.Index = "End.Index")

## Froh of all individuals
Froh_test_allId_from16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_from16Mb, return = 1, by.chr = FALSE)

## Froh of all individuals by Chromosome
# Froh_test__allId_byChr_from16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_from16Mb, return = 1, by.chr = TRUE)

## Mean Froh
Froh_test_mean_from16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_from16Mb, return = 2,  by.chr = FALSE)

## Mean Froh by Family ID
# Froh_test_mean_byFid_from16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_from16Mb, return = 2, by.fid = TRUE, by.chr = FALSE)

## Mean Froh by Chromosome
# Froh_test_meanbyChr_from16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_from16Mb, return = 2, by.chr = TRUE)

## Mean Froh by Family ID inside each Chromosome
# Froh_test_meanbyChr_byFid_from16Mb <- Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test_from16Mb, return = 2, by.fid = TRUE, by.chr = TRUE)


###############


### Merging results together


Froh_test <- cbind(Froh_test_allId,Froh_test_allId_1to2Mb[,3],Froh_test_allId_2to4Mb[,3],Froh_test_allId_4to8Mb[,3],Froh_test_allId_8to16Mb[,3],Froh_test_allId_from16Mb[,3])
colnames(Froh_test) <- c('fid','ID','Froh','Froh_1to2Mb','Froh_2to4Mb','Froh_4to8Mb','Froh_8to16Mb','Froh_from16Mb')
# If some object is not found -> delete it in both rows
#For example: object 'Froh_test_allId_from16Mb' not found -> Delete Froh_test_allId_from16Mb[,3] from first row and 'Froh_from16Mb' from second row


## Save Froh_test file
write.table(Froh_test,'Froh_byID_byCLASS.txt',quote = FALSE, row.names = FALSE, col.names = TRUE)


################################################################################


###  6. Regional Froh estimation


## Using same formula as before, estimation of inbreeding across genome


## Regional.Froh.estimation() function is used:

# Arguments: return, normalization, by.fid, start, slide, end, Min.Percent.Window.Size

# Options:
# return = 1  - Returns mean Regional Froh values (with LogPvalue for Froh Coldspot and Hotspot windows)
# return = 2  - Returns Regional Froh values for all individuals

# normalization = 'by.chr' - Normalisation for each Chromosome separately 
# normalization = 'all'    - Normalisation for all Chromosomes together

# by.fid = FALSE                               - All individuals are treated as one group
# by.fid = TRUE                                - You get results for all breeds separately

# start                    - start position in bp
# end                      - end position in bp
# slide                    - size of slide in bp

# Min.Percent.Window.Size = 0.75 (default) -> Your minimum window size for Regional Froh estimation is --> (end-start) * Min.Percent.Window.Size --> otherwise calculation is skipped [it will happen if gap between two SNPs in you dataset is higher than Min.Percent.Window.Size*(end-start)]


## RECOMENDATION - start keep at 0; unless you are working with specific region 

## 5 combinations:  a) return = 1, by.fid = FALSE, normalization = 'by.chr',         --- returns mean Regional Froh values (Normalization by Chr)
#    (results)      b) return = 1, by.fid = FALSE, normalization = 'all',            --- returns normalized mean Regional Froh values
#                   c) return = 1, by.fid = TRUE, normalization = 'by.chr',          --- returns mean Regional Froh values by breed (Normalization by Chr)
#                   d) return = 1, by.fid = TRUE, normalization = 'all',             --- returns normalized mean Regional Froh values by breed

#                   e) return = 2, by.fid = FALSE                                    --- returns Regional Froh values for all individuals


## Regional Froh estimation
## Recommended values: Window size = 10000000 (10 Mb), slide = 500000 (500 kb)
# Regional.Froh_test_meanFroh_NormAll <- Regional.Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 1, by.fid = FALSE, normalization = 'all', slide = 500000, end = 10000000)

# Regional.Froh_test_meanFroh_NormbyChr <- Regional.Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 1, by.fid = FALSE, slide = 500000, end = 10000000)
# Regional.Froh_test_meanFrohbyBreed_NormbyChr <- Regional.Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 1, by.fid = TRUE, slide = 500000, end = 10000000)
# Regional.Froh_test_meanFrohbyBreed_NormAll <- Regional.Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 1, by.fid = TRUE, normalization = 'all', slide = 500000, end = 10000000)
# Regional.Froh_test_byIdFroh <- Regional.Froh.estimation(map = map_file, ped = ped_file, I_matrix = I_matrix_test, return = 2, slide = 500000, end = 10000000)


# Example: Chromosome 30 only (using all chromosomes would take too long for the workshop)
Example_I_Index <- I_matrix_test
colnames(Example_I_Index) <- c(30,map_file$V1)
Regional.Froh_test_meanFroh_NormAll <- Regional.Froh.estimation(map = map_file[map_file$V1 == 30,], ped = ped_file, I_matrix = Example_I_Index[,colnames(Example_I_Index) == 30], return = 1, by.fid = FALSE, normalization = 'all', slide = 500000, end = 10000000)


################################################################################


### 7. Visualisation of Froh results


## There are many ways to visualize data -> see ggplot2 package


## In this workshop, you will first see boxplot and stacked plot visualizations 
## of Froh values, followed by plot of regional Froh values


## Following packages will be used for the visualizations: ggplot2 and patchwork
install.packages('ggplot2')
install.packaged('patchwork')
library(ggplot2)
library(patchwork)


## Preparing data for visualization

Froh_test_allId_1to2Mb$ROH <- '1-2 Mb'
Froh_test_allId_2to4Mb$ROH <- '2-4 Mb'
Froh_test_allId_4to8Mb$ROH <- '4-8 Mb'
Froh_test_allId_8to16Mb$ROH <- '8-16 Mb'
Froh_test_allId_from16Mb$ROH <- '>16 Mb'

# Combining all datasets into single object for plotting
Froh_forPlot <- rbind(Froh_test_allId_from16Mb, Froh_test_allId_8to16Mb, Froh_test_allId_4to8Mb,
                      Froh_test_allId_2to4Mb, Froh_test_allId_1to2Mb)
Froh_forPlot$ROH <- factor(Froh_forPlot$ROH, levels = c(">16 Mb", "8-16 Mb", "4-8 Mb", "2-4 Mb", "1-2 Mb"))


## Visualization 

# Creating boxplot
px <- ggplot(Froh_test_allId, aes(x = factor(fid, levels = c('LAB')), 
                                  y = Froh, fill = fid)) +
  geom_boxplot(alpha = 0.4) +
  scale_fill_manual(values = c("lightpink2")) +    # Single color for a single breed
  theme(legend.position = "none") +
  ylab(expression("F"["ROH"])) +
  xlab("") +
  ggtitle("Boxplot for Labrador Retriever")
px

# Creating stacked plot
py <- ggplot(Froh_forPlot, aes(x = factor(fid, levels = c('LAB')), 
                               fill = ROH, y = Froh)) +
  geom_bar(stat = "summary", fun = "mean") +
  scale_fill_manual(values = c("yellow", "gray80", "gray60", "gray40", "black")) +
  theme_classic() +
  xlab("") +
  ylab(expression("F"["ROH"])) +
  ggtitle("Stacked Plot for Labrador Retriever")
py


# Combining boxplot and stacked plot
combined_plot <- px + py + plot_layout(ncol = 1)
combined_plot


## Regional visualization
# Plotting regional Froh values
plot(Regional.Froh_test_meanFroh_NormAll$start_position / 1000000, 
     Regional.Froh_test_meanFroh_NormAll$mean_Froh, 
     col = "red", pch = 19, 
     ylim = c(0, 0.3), 
     xlab = 'Chromosome 30 (Mb)', 
     ylab = expression("F"['ROH']))


## Exporting plots in high quality

dpi=600;width.cm<-20;height.cm<-10;width.in<-width.cm/2.45;height.in<-height.cm/2.45
png(file="./Froh_Boxplot_And_StackedPlot.png",width = width.in*dpi,height=height.in*dpi,units="px",res=dpi)
par(mar = c(4,4, 4, 4))
combined_plot
dev.off()

png(file="./Regional_Froh.png",width = width.in*dpi,height=height.in*dpi,units="px",res=dpi)
par(mar = c(4,4, 4, 4))
plot(Regional.Froh_test_meanFroh_NormAll$start_position / 1000000, 
     Regional.Froh_test_meanFroh_NormAll$mean_Froh, 
     col = "red", pch = 19, 
     ylim = c(0, 0.3), 
     xlab = 'Chromosome 30 (Mb)', 
     ylab = expression("F"['ROH']))
dev.off()


################################################################################


###  8. Creation of Genomic Relationship Matrix (GRM)


## Used both for estimation of inbreeding and relatedness


## Values = from -1 to 1


## VanRaden Matrix will be used

#######################################################
# VanRaden (2008) J. Dairy Sci. 91:4414-4423
#######################################################


## First step is to load GRM.R -> Solomon Boison's script: calc_gnrm() function
source("GRM.R")


## Then recode Plink file -> present genotypes as numeric values (1 or 2, 0 = NA) -> for calc_gnrm()
system("plink.exe --file LabradorRetriever_FinalReport_afterQC --nonfounders --output-missing-genotype 0 --recode12 --out Recoded_for_VanRadenGRM --dog")


## VanRaden GRM calculation
# vanRaden <- calc_gnrm(genofile='Recoded_for_VanRadenGRM.ped',genoformat="ped",ana_type="vanRaden",ped_option=F,
#                       outputformat="matrix",outputname="GRM",nIID=60,missinggeno=T,plots=T)

# vanRaden

# Due to long runtime, for example we will use pruned data created with Plink:
system("plink.exe --file Recoded_for_VanRadenGRM --indep-pairwise 50 5 0.3 --out pruned_data --dog")
system("plink.exe --file Recoded_for_VanRadenGRM --extract pruned_data.prune.in --recode --out reduced_genotype --dog")
# This command uses a sliding window of 50 SNPs, shifts the window by 5 SNPs at a time and removes SNPs with an r^2 value above 0.3

vanRaden <- calc_gnrm(genofile='reduced_genotype.ped',genoformat="ped",ana_type="vanRaden",ped_option=F,
                      outputformat="matrix",outputname="GRM",nIID=60,missinggeno=T,plots=T)

vanRaden


## GRM is giving you information about relatedness between individual


## To get information about inbreeding (based on IBS status of alleles) substract value 1 from diagonal of GRM
diag(vanRaden) <- diag(vanRaden) - 1
diag(vanRaden)
summary(diag(vanRaden))


write.table(vanRaden,'GRM.grm', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


################################################################################


###  9. GRM Visualization


## Following package will be used for visualization: gplots
install.packages('gplots')
library(gplots)


## Heatmap -> for relatedness visualization
heatmap.2(vanRaden,
          key.xlab="Relatedness",                     # Add a custom label for the x-axis of the color key
          key.ylab="Count",                           # Add a custom label for the y-axis of the color key
          col=bluered(100),                           # Use a blue-red color palette with 100 levels
          trace="none",                               # Disable trace lines on the heatmap
          margins=c(6, 6),                            # Set margins around the heatmap for row and column labels
          cexRow=0.7,                                 # Adjust the font size for row labels
          cexCol=0.7,                                 # Adjust the font size for column labels
          key.title="Color Key and Histogram",        # Customize the title of the color key
          key.par=list(cex=0.8))                      # Reduce font size for the legend title

# Alternatively, inbreeding values can be set to zero
# diag(vanRaden) <- 0
# heatmap.2(vanRaden,
#          key.xlab="Relatedness",                     # Add a custom label for the x-axis of the color key
#          key.ylab="Count",                           # Add a custom label for the y-axis of the color key
#          col=bluered(100),                           # Use a blue-red color palette with 100 levels
#          trace="none",                               # Disable trace lines on the heatmap
#          margins=c(6, 6),                            # Set margins around the heatmap for row and column labels
#          cexRow=0.7,                                 # Adjust the font size for row labels
#          cexCol=0.7,                                 # Adjust the font size for column labels
#          key.title="Color Key and Histogram",        # Customize the title of the color key
#          key.par=list(cex=0.8))                      # Reduce font size for the legend title


## Histrogram -> for inbreeding visualization
hist(diag(vanRaden), main="GRM Inbreeding distribution", xlab=expression("F"["VanRaden"]))


## Exporting heatmap and histogram in high quality

png(file="./VanRaden_GRM_Heatmap_Relatedness.png", width = 10*300, height = 10*300, units="px", res=300)
par(mar = c(5,5, 5, 5))
heatmap.2(vanRaden,
          key.xlab="Relatedness",                     # Add a custom label for the x-axis of the color key
          key.ylab="Count",                           # Add a custom label for the y-axis of the color key
          col=bluered(100),                           # Use a blue-red color palette with 100 levels
          trace="none",                               # Disable trace lines on the heatmap
          margins=c(6, 6),                            # Set margins around the heatmap for row and column labels
          cexRow=0.7,                                 # Adjust the font size for row labels
          cexCol=0.7,                                 # Adjust the font size for column labels
          key.title="Color Key and Histogram",        # Customize the title of the color key
          key.par=list(cex=0.8))                      # Reduce font size for the legend title
dev.off()

png(file="./VanRaden_GRM_Histogram_Inbreeeding.png",width = width.in*dpi,height=height.in*dpi,units="px",res=dpi)
par(mar = c(4,4, 4, 4))
hist(diag(vanRaden), main="GRM Inbreeding distribution", xlab=expression("F"["VanRaden"]))
dev.off()


################################################################################




