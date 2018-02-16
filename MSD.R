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
msd_single = as.data.frame(read.csv(paste(prefix,"_MSD_single.txt",sep=""),sep="\t"))

n_frames = nrow(msd)

both_MSD = merge(msd_tau,msd,by="TAU",all=T)
both_MSD = melt(both_MSD,id="TAU")

single = ggplot(msd_single, aes(x=as.numeric(factor(TAU))*frame_int,y=MSD,color=factor(TRACK))) + geom_line() + theme_bw(base_size = 12) + scale_color_discrete(guide=FALSE) + labs(x = "Time (s)",y="MSD tau (um²)") + ggtitle("MSD tau (All particles)")
single$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_single_trajectories.png",sep=""), plot=single, width = 11, height = 8.5, dpi=300)

bins = ggplot(msd_bin) + geom_line(aes(x=TAU*frame_int,y=MSD,group=factor(BIN),color=BIN),size=1.1) + theme_bw(base_size=12) + labs(x="Time (s)",y="MSD") + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("2","4","6","8","10"),breaks=c(2,4,6,8,10)) + ggtitle("Ensemble averaged MSD - Binned by intensity")
bins$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_bins.png",sep=""), plot=bins, width = 11, height = 8.5, dpi=300)

ens = ggplot(msd) + geom_line(aes(x=TAU*frame_int,y=MSD)) + theme_bw(base_size=12) + labs(x="Time (s)",y="MSD") + ggtitle("Ensemble average MSD")
ens$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average.png",sep=""), plot=ens, width = 11, height = 8.5, dpi=300)

time_ens = ggplot(msd_tau, aes(x=TAU*frame_int,y=VALUE)) + geom_line() + theme_bw(base_size = 12) + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau") + ggtitle("MSD tau (Time-ensemble average)")
ggsave(paste(prefix,"_MSD_time-ensemble_average.png",sep=""), plot=time_ens, width = 11, height = 8.5, dpi=300)

both = ggplot(both_MSD) + geom_line(aes(x=TAU*frame_int,y=value,group=variable,color=variable),size=0.7) + theme_bw(base_size = 12) + labs(x="Time (s)",y="MSD") + scale_colour_manual(values=c("firebrick1","dodgerblue"),labels=c("Time", "Ensemble"))
both$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_both.png",sep=""), plot=both, width = 11, height = 8.5, dpi=300)

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

png(
  paste(prefix,"_4panel.png",sep=""),
  width     = 11.25,
  height    = 8.25,
  units     = "in",
  res       = 1000,
  pointsize = 0.01
)
par(
  mar      = c(5, 5, 2, 2),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 1,
  cex.lab  = 1
)
multi = multiplot(single, bins, ens, time_ens, cols=2)
dev.off()