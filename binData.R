### splits the trajectories in 10 bins according to mean particle intensity

arg = commandArgs(trailingOnly=TRUE)
prefix = arg[1]

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(data.table, Hmisc)

print("Importing tracking file...")
tracking_file = as.data.frame(read.csv(paste(prefix,"_tracking_noHeader.txt",sep=""),sep=" ",header = FALSE))
colnames(tracking_file) = c("track","frame","X","Y","mean")

print("Binning trajectories...")
tracking_file$group <- as.numeric(cut2(tracking_file$mean, g=10))

tracking_file = tracking_file[with(tracking_file, order(track, frame)),]

print("Writing to file...")
write.table(tracking_file, file = paste(prefix,"_tracking_noHeader.txt",sep=""),row.names=FALSE, col.names=FALSE, sep="\t")