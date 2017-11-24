arg = commandArgs(trailingOnly=TRUE)

n_frames = as.numeric(arg[1])
prefix = arg[2]

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2, reshape2)

print("Importing MSD data...")

#### IMPORT MELTED DF
msd = as.data.frame(read.csv(paste(prefix,"_melt.txt",sep=""),sep="\t"))

## Order by frame
msd = msd[order(msd$FRAME),]

## Reconstruct dataframe
unmelt = dcast(msd,TRACK~FRAME)

frame_loss = 401 - rowSums(is.na(unmelt[,2:ncol(unmelt)]))
n_per_frame = rowSums(sapply(X=frame_loss,y=seq(1,400), FUN=function(x,y) x > y))
time_avg = colMeans(unmelt[,-1],na.rm = TRUE)
time_sd = apply(unmelt[,2:ncol(unmelt)], 2, function(x) sd(x, na.rm=TRUE))
tau=seq(from=0,to=(n_frames*2.048)-3,by=2.048)
time_avg_MSD = cbind(tau,time_avg,time_sd)
time_avg_MSD = as.data.frame(time_avg_MSD)

a = ggplot(msd, aes(x=factor(FRAME),y=VALUE,group=TRACK,color=factor(TRACK))) + geom_line() + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + scale_x_discrete(limits=0:400, breaks = seq(0,400,25)) +
labs(x = "Time",y="MSD") + ggtitle("MSD - All particles")
a$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_single_trajectories.png",sep=""), plot=a, width = 11, height = 8.5, dpi=300)

b = ggplot(time_avg_MSD, aes(x=tau,y=time_avg)) + geom_errorbar(aes(ymin=time_avg-time_sd, ymax=time_avg+time_sd, color="grey"), width=.05) + geom_line() + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD") + ggtitle("MSD - Ensemble average")
b$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average_SD.png",sep=""), plot=b, width = 11, height = 8.5, dpi=300)

c = ggplot(time_avg_MSD, aes(x=tau,y=time_avg)) + geom_line() + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD") + ggtitle("MSD - Ensemble average")
c$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average.png",sep=""), plot=c, width = 11, height = 8.5, dpi=300)

d = ggplot(time_avg_MSD, aes(x=tau)) + geom_line(aes(y=time_avg, colour="blue")) + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD") + ggtitle("MSD - Ensemble average vs. number of tracks per frame") + geom_line(aes(y = n_per_frame/25, colour = "red")) + scale_y_continuous(sec.axis = sec_axis(~., name = "Number of contributing trajectories (x25)"))
d$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average_nb_tracks.png",sep=""), plot=d, width = 11, height = 8.5, dpi=300)
