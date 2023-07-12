#!/usr/bin/env Rscript

requiredPackages = c('ggplot2')
suppressMessages(
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE, quietly = TRUE)
})

args <- (commandArgs(TRUE))
directory <- args[2]

name_bam <- args[1]


df <- read.csv(paste(directory, "all_info.csv", sep = ""))

file_directory = paste(directory,name_bam, sep="")
pdf(file=paste(file_directory, "_Rplots2.pdf", sep= ""))

#Plot of p_est and r_fixed against GC content
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(df$GC_bins, df$p_estimated, xlab="GC content bins", ylab="P estimated") # first plot
par(new = TRUE)
plot(df$GC_bins, df$r_fixed, type = "l", lty = 4, axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(df$r_fixed)))
mtext("R fixed", side=4, line=3)

# Plot lambda with error bars
plot(df$GC_bins, df$Lambda, xlab="GC content bins", ylab="Predicted read counts", pch=20, cex=2)
# Add error bars
arrows(x0=df$GC_bins, y0=df$Lambda-df$Standard_deviations, x1=df$GC_bins, y1=df$Lambda+df$Standard_deviations, code=3, angle=90, length=0.1)

notprint<-dev.off()

################################################################

data <- read.csv(paste(directory, "conditions_rc.csv", sep = ""))

plot_list = list()
for(i in 3:ncol(data)) {
  #variable <-names(data)[x+2]
  #print(predicted)
  condition <-names(data[, i])
  plot<-ggplot(data, aes(x=rau)) +
    geom_point(aes(y=Observed_rc), col="blue") +
    geom_line(aes(y=unlist(data[ ,i])), col="red") +
    scale_y_continuous(name="Observed", sec.axis=sec_axis(~./1, name="Predicted")) +
    ggtitle(condition) +
    theme(
      axis.title.y.left=element_text(color="blue"),
      axis.text.y.left=element_text(color="blue"),
      axis.title.y.right=element_text(color="red"),
      axis.text.y.right=element_text(color="red")
    )
  plot_list[[i]] <- plot

}


pdf(file=paste(file_directory, "_Rplots3.pdf", sep= ""))
for(i in 3:ncol(data)) {
    print(plot_list[[i]])
}
not_print<-dev.off()