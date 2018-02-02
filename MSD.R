arg = commandArgs(trailingOnly=TRUE)

prefix = arg[1]
frame_int = as.numeric(arg[2])

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2, reshape2, Hmisc, viridis)

print("Importing MSD data...")

#### IMPORT MELTED DF

msd = as.data.frame(read.csv(paste(prefix,"_MSD_all.txt",sep=""),sep="\t"))
msd_bin = as.data.frame(read.csv(paste(prefix,"_MSD_binned.txt",sep=""),sep="\t"))
msd_tau = as.data.frame(read.csv(paste(prefix,"_MSD_tau.txt",sep=""),sep="\t"))

n_frames = nrow(msd)

## Order by frame
# msd_single = msd_tau

## Reconstruct dataframe
# unmelt = dcast(msd_tau,TRACK~FRAME)

# msd_tau = colMeans(unmelt[,-1],na.rm = TRUE)
# MELT_tau = cbind(unique(msd_single$FRAME),msd_tau)
# MELT_tau = as.data.frame(MELT_tau)
# colnames(MELT_tau)[1]="TAU"

# both_MSD = cbind(msd_tau,msd)
both_MSD = merge(msd_tau,msd,by="TAU",all=T)
both_MSD = melt(both_MSD,id="TAU")

# single = ggplot(msd_single, aes(x=as.numeric(factor(FRAME))*frame_int,y=VALUE)) + geom_line() + theme_bw(base_size = 18) + scale_color_discrete(guide=FALSE) + labs(x = "Time (s)",y="MSD tau (um²)") + ggtitle(paste("MSD tau (All particles) - ",prefix,sep=""))
# single$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
# ggsave(paste(prefix,"_MSD_single_trajectories.png",sep=""), plot=single, width = 11, height = 8.5, dpi=300)

bins = ggplot(msd_bin) + geom_line(aes(x=TAU*frame_int,y=MSD,group=factor(BIN),color=BIN),size=1.2) + theme_bw(base_size=20) + labs(x="Time (s)",y="MSD") + scale_color_viridis(guide=guide_legend(title="Bin number",keywidth = 2, keyheight = 2,override.aes = list(alpha = 1),byrow = TRUE),labels=c("2","4","6","8","10"),breaks=c(2,4,6,8,10)) + ggtitle(paste("Ensemble averaged MSD - Binned by intensity - ",prefix,sep=""))
bins$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_bins.png",sep=""), plot=bins, width = 11, height = 8.5, dpi=300)

ens = ggplot(msd) + geom_line(aes(x=TAU*frame_int,y=MSD)) + theme_bw(base_size=18) + labs(x="Time (s)",y="MSD") + ggtitle(paste("Ensemble average MSD - ",prefix,sep=""))
ens$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average.png",sep=""), plot=ens, width = 11, height = 8.5, dpi=300)

time_ens = ggplot(msd_tau, aes(x=TAU*frame_int,y=VALUE)) + geom_line() + theme_bw(base_size = 18) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau") + ggtitle(paste("MSD tau (Time-ensemble average) - ",prefix,sep=""))
ggsave(paste(prefix,"_MSD_time-ensemble_average.png",sep=""), plot=time_ens, width = 11, height = 8.5, dpi=300)

both = ggplot(both_MSD) + geom_line(aes(x=TAU*frame_int,y=value,group=variable,color=variable),size=0.7) + theme_bw(base_size = 18) + labs(x="Time (s)",y="MSD") + scale_colour_manual(values=c("firebrick1","dodgerblue"),labels=c("Time", "Ensemble"))
both$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_both.png",sep=""), plot=both, width = 11, height = 8.5, dpi=300)