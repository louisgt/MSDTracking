### splits the trajectories in 10 bins according to mean particle intensity

arg = commandArgs(trailingOnly=TRUE)
prefix = arg[1]

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(data.table, Hmisc)

print("Importing tracking file...")
tracking_file = as.data.frame(read.csv(paste(prefix,"_tracking_noHeader.txt",sep=""),sep=" ",header = FALSE))
colnames(tracking_file) = c("track","frame","X","Y","intensity")
intensities = aggregate(tracking_file[,"intensity"], list(tracking_file$track), mean)
stdevs = aggregate(tracking_file[,"intensity"], list(tracking_file$track), sd)
colnames(intensities) = c("track","mean")
colnames(stdevs) = c("track","sd")

print("Binning trajectories...")
D1 = data.table(tracking_file, key="track")
D2 = data.table(intensities, key = "track")
D3 = data.table(stdevs, key = "track")
D_merge = D1[D2,nomatch=0]
D_merge = D_merge[D3,nomatch=0]
tracking_file$group <- as.numeric(cut2(D_merge$mean, g=10))
tracking_file$sd <- as.numeric(D_merge$sd)


print("Writing to file...")
write.table(tracking_file, file = paste(prefix,"_tracking_noHeader.txt",sep=""),row.names=FALSE, col.names=FALSE, sep="\t")