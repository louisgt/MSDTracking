### normalizes intensity for binning

arg = commandArgs(trailingOnly=TRUE)
prefix = arg[1]

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(data.table, Hmisc)

print("Importing tracking file...")
tracking_file = as.data.frame(read.csv(paste(prefix,"_tracking_noHeader.txt",sep=""),sep=" ",header = FALSE))
colnames(tracking_file) = c("track","frame","X","Y","intensity")
intensities = aggregate(tracking_file[,"intensity"], list(tracking_file$track), mean)
#stdevs = aggregate(tracking_file[,"intensity"], list(tracking_file$track), sd)
colnames(intensities) = c("track","mean")
#colnames(stdevs) = c("track","sd")

by_frame = aggregate(tracking_file[,"intensity"], list(tracking_file$frame), mean,na.rm=TRUE)
colnames(by_frame) = c("frame","normalized")
tracking_file = merge(tracking_file,by_frame,by.x="frame",by.y="frame")
tracking_file$normalized = tracking_file$intensity/tracking_file$normalized

mean = aggregate(tracking_file[,"normalized"], list(tracking_file$track), mean)
colnames(mean)=c("track","mean")
tracking_file = merge(tracking_file,mean,by.x="track",by.y="track")

tracking_file = tracking_file[,-c(5,6)]
tracking_file = tracking_file[with(tracking_file, order(track, frame)),]

print("Writing to file...")
write.table(tracking_file, file = paste(prefix,"_tracking_noHeader.txt",sep=""),row.names=FALSE, col.names=FALSE, sep="\t")