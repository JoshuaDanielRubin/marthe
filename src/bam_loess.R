#!/usr/bin/env Rscript

suppressMessages(
library("stringr"))


args <- (commandArgs(TRUE))
directory <- args[3]
name_bam <- args[4]
file_directory = paste(directory,name_bam, sep="")

#Creating GC content dataframe
dataGC <- read.table(args[1], header = FALSE)




#Creating read counts dataframe
dataCounts <- read.table(args[2], header = FALSE)
dataCounts <- dataCounts[c("V2")]   #select second columns, non unique reads mapped to reference
#This has changed, in the previous version I was taking the second column




#Creating combined dataframe with GC and read counts
data <-  dataGC[c("V2")]        #Select only second column, it contains the GC percentages

data$count <- unlist(dataCounts)        #append dataCounts list

#data <- data [complete.cases(data), ]


#It removes rows with missing values in any column of data frame
data_new <- data[data[,1]!=-1,]

colnames(data_new)[1] = "gc"


#Application of loess
loessMod10 <- loess(count ~ gc, data=data_new, span = 0.70) # 70% smoothing span


#Loess regression is a non parametric technique that uses local weighted regression to fit
# a smooth curve through points in a scatter plot.

#Using gc content as predictor for the loess
smoothed10 <- predict(loessMod10, new_data=data_new$gc)

#Activate when needed to get the plots
pdf(file=paste(file_directory, "_Rplots1.pdf", sep= ""))
plot(data_new$gc, data_new$count, xlab="GC content", ylab="Read Count", pch=19, main='Loess Regression Models')

#plot to check if it is correct
plot(data_new$gc, data_new$count, xlab="GC content", ylab="Read Count", pch=19, main='Loess Regression Models')
lines(smoothed10, x=data_new$gc, col='red')


#actual read count - prediction
residues <- data_new$count - smoothed10
data_new$count_mod <- residues


#Plot the corrected read counts
#plot(data_new$gc, data_new$count_mod, xlab="gc content", ylab="corrected depth", pch=19, main='Corrected GC bias')

#Activate when needed to get the plots
notprint <- dev.off()

#It prints the values of gc and the predictions (read counts) in a table
cat(rbind(data_new$gc, smoothed10), sep = "\n") #--> this is the original version, recover if required

