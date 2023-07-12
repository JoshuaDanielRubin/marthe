#!/usr/bin/env Rscript

# Load required packages or install them if not present
requiredPackages <- c('ggplot2')
for (p in requiredPackages) {
  if (!require(p, character.only = TRUE)) install.packages(p)
  library(p, character.only = TRUE, quietly = TRUE)
}

# Get arguments from command line
args <- commandArgs(TRUE)
name_bam <- args[1]
directory <- args[2]

# Set file paths
all_info_path <- file.path(directory, "all_info.csv")
conditions_rc_path <- file.path(directory, "conditions_rc.csv")
pdf_file1 <- file.path(directory, paste0(name_bam, "_Rplots2.pdf"))
pdf_file2 <- file.path(directory, paste0(name_bam, "_Rplots3.pdf"))

# Load data
df <- read.csv(all_info_path)
data <- read.csv(conditions_rc_path)

# Plot p_est and r_fixed against GC content
pdf(pdf_file1)
par(mar = c(5, 4, 4, 4) + 0.3)  # Leave space for z axis
plot(df$GC_bins, df$p_estimated, xlab="GC content bins", ylab="P estimated") # first plot
par(new = TRUE)
plot(df$GC_bins, df$r_fixed, type = "l", lty = 4, axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(df$r_fixed)))
mtext("R fixed", side=4, line=3)

# Plot lambda with error bars
plot(df$GC_bins, df$Lambda, xlab="GC content bins", ylab="Predicted read counts", pch=20, cex=2)
arrows(x0=df$GC_bins, y0=df$Lambda-df$Standard_deviations, x1=df$GC_bins, y1=df$Lambda+df$Standard_deviations, code=3, angle=90, length=0.1)
dev.off()

# Loop over conditions to generate plots
plot_list <- list()
for (i in 3:ncol(data)) {
  condition <- names(data[, i])
  plot <- ggplot(data, aes(x=rau)) +
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

# Save plots to PDF
pdf(pdf_file2)
lapply(plot_list, print)
dev.off()

