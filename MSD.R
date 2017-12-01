arg = commandArgs(trailingOnly=TRUE)

prefix = arg[1]

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2, reshape2)

print("Importing MSD data...")

#### IMPORT MELTED DF

msd = as.data.frame(read.csv(paste(prefix,"_MSD.txt",sep=""),sep="\t"))
msd_tau = as.data.frame(read.csv(paste(prefix,"_MSD_tau.txt",sep=""),sep="\t"))

n_frames = nrow(msd)

## Order by frame
msd_tau = msd_tau[order(msd_tau$FRAME),]
msd_single = msd_tau

## Reconstruct dataframe
unmelt = dcast(msd_tau,TRACK~FRAME)

# frame_loss = 401 - rowSums(is.na(unmelt[,2:ncol(unmelt)]))
# n_per_frame = rowSums(sapply(X=frame_loss,y=seq(1,401), FUN=function(x,y) x > y))

msd_tau = colMeans(unmelt[,-1],na.rm = TRUE)
time_sd = apply(unmelt[,2:ncol(unmelt)], 2, function(x) sd(x, na.rm=TRUE))
tau=seq(from=0,to=(n_frames*2.048)-1,by=2.048)
MELT_tau = cbind(tau,msd_tau,time_sd)
MELT_tau = as.data.frame(MELT_tau)

both_MSD = cbind(msd_tau,msd)
both_MSD = melt(both_MSD,id="TAU")

a = ggplot(msd_single, aes(x=factor(FRAME),y=VALUE,group=TRACK,color=factor(TRACK))) + geom_line() + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + scale_x_discrete(limits=0:400, breaks = seq(0,400,25)) +
labs(x = "Time",y="MSD tau") + ggtitle(paste("MSD tau (All particles) - ",prefix,sep=""))
a$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_single_trajectories.png",sep=""), plot=a, width = 11, height = 8.5, dpi=300)

b = ggplot(MELT_tau, aes(x=tau,y=msd_tau)) + geom_errorbar(aes(ymin=msd_tau-time_sd, ymax=msd_tau+time_sd, color="grey"), width=.05) + geom_line() + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau") + ggtitle(paste("MSD tau (Time-ensemble average) - ",prefix,sep=""))
b$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_time-ensemble_average_SD.png",sep=""), plot=b, width = 11, height = 8.5, dpi=300)

c = ggplot(MELT_tau, aes(x=tau,y=msd_tau)) + geom_line() + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau") + ggtitle(paste("MSD tau (Time-ensemble average) - ",prefix,sep=""))
c$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_time-ensemble_average.png",sep=""), plot=c, width = 11, height = 8.5, dpi=300)

# d = ggplot(MELT_tau, aes(x=tau)) + geom_line(aes(y=msd_tau, colour="blue")) + theme_bw(base_size = 16) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau") + ggtitle(paste("Time-ensemble average MSD vs Number of tracks per frame - ",prefix,sep="")) + geom_line(aes(y = n_per_frame/25, colour = "red")) + scale_y_continuous(sec.axis = sec_axis(~., name = "Number of contributing trajectories (x25)"))
# d$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
# ggsave(paste(prefix,"_MSD_time-ensemble_average_nb_tracks.png",sep=""), plot=d, width = 11, height = 8.5, dpi=300)

e = ggplot(msd) + geom_line(aes(x=TAU,y=MSD)) + theme_bw(base_size=16) + labs(x="Time (s)",y="MSD") + ggtitle(paste("Ensemble average MSD vs Number of tracks per frame - ",prefix,sep=""))
e$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average.png",sep=""), plot=e, width = 11, height = 8.5, dpi=300)

f = ggplot(both_MSD) + geom_line(aes(x=TAU,y=value,group=variable,color=variable),size=0.7) + theme_bw(base_size = 16) + labs(x="Time (s)",y="MSD") + scale_colour_manual(values=c("firebrick1","dodgerblue"))
f$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_both.png",sep=""), plot=f, width = 11, height = 8.5, dpi=300)