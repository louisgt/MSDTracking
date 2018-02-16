### splits the trajectories in 10 bins according to mean particle intensity

arg = commandArgs(trailingOnly=TRUE)
prefix = arg[1]

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(data.table, Hmisc)

print("Importing tracking file...")
tracking_file = as.data.frame(read.csv(paste(prefix,"_tracking_noHeader.txt",sep=""),sep=" "))
colnames(tracking_file) = c("track","X","Y","intensity")
intensities = aggregate(tracking_file[,"intensity"], list(tracking_file$track), mean)
colnames(intensities) = c("track","mean")

print("Binning trajectories...")
D1 = data.table(tracking_file, key="track")
D2 = data.table(intensities, key = "track")
D_merge = D1[D2,nomatch=0]
tracking_file$group <- as.numeric(cut2(D_merge$mean, g=10))

print("Writing to file...")
write.table(tracking_file, file = paste(prefix,"_tracking_noHeader.txt",sep=""),row.names=FALSE, col.names=FALSE, sep="\t")