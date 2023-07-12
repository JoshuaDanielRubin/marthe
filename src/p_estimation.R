#!/usr/bin/env Rscript


requiredPackages = c('EnvStats','comprehenr')
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
#install.packages('EnvStats')
#install.packages("comprehenr")
#library("EnvStats")
#library("comprehenr")

args <- (commandArgs(TRUE))


#Reading the file line by line
lines <- readLines(args[1])


#Estimate the probability of success for each line of the input file
#Each line of the input file correspond to the observed read counts for each GC bin
res_list = list()
tmp_list = list()
for (line in lines) {
    data <- strsplit(line, split = "\t")

    #Save the data as a vector of floats
    v <- unlist(data)
    v <- as.vector(v, "numeric")

    #Estimation of probability fo success given a fixed number of successes, r = 5
    res <- enbinom(v, size = 5)
    res_tosave <- res$parameters[2]

    #Save the results for the different GC bins in a list
    res_tosave <-  as.numeric(res_tosave)
    tmp_list <- list(p = res_tosave, r = 5, m = mean(v), sd = sd(v))

    res_list <- c(res_list, tmp_list)
    #res_list <- c(res_list, res_tosave)
    tmp_list = list()
}


vec <- paste0(unlist(res_list), collapse = "\n")

#Display the list of p values
cat(vec)