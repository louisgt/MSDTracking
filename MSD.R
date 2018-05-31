arg = commandArgs(trailingOnly=TRUE)

prefix = arg[1]
frame_int = as.numeric(arg[2])

#### LOAD UTILITIES (INSTALL IF MISSING)
if (!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")
pacman::p_load(ggplot2, reshape2, Hmisc, viridis,lme4)

print("Importing MSD data...")

#"publication ready" plot theme
theme_Publication <- function(base_size=18, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size=1.5),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

#### IMPORT MELTED DF

msd = as.data.frame(read.csv(paste(prefix,"_MSD_all.txt",sep=""),sep="\t"))
msd_bin = as.data.frame(read.csv(paste(prefix,"_MSD_binned.txt",sep=""),sep="\t"))
msd_tau = as.data.frame(read.csv(paste(prefix,"_MSD_tau.txt",sep=""),sep="\t"))
msd_single = as.data.frame(read.csv(paste(prefix,"_MSD_single.txt",sep=""),sep="\t"))

n_frames = nrow(msd)

both_MSD = merge(msd_tau,msd,by="TAU",all=T)
both_MSD = melt(both_MSD,id="TAU")

msd$TAU = msd$TAU*frame_int
msd_bin$TAU = msd_bin$TAU*frame_int

###########

displacements <- as.data.frame(read.csv(paste(prefix,"_displacements.txt",sep=""),sep="\t",header = FALSE))
h=hist(displacements$V2,breaks = 300,plot=FALSE)
x=h$mids
y=h$density
power_law = cbind(x,y)
power_law = as.data.frame(power_law)

###########
## PLOT THE MSD OF EACH TRAJECTORY SEPARATELY

single = ggplot(msd_single, aes(x=as.numeric(factor(TAU))*frame_int,y=MSD,color=factor(TRACK))) + geom_line() + theme_Publication() + scale_color_discrete(guide=FALSE) + labs(x = "Time (s)",y="MSD tau (um²)") + ggtitle("MSD tau (All particles)")
single$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_single_trajectories.png",sep=""), plot=single, width = 11, height = 8.5, dpi=300)

###########
## PLOT ENSEMBLE averaged MSD binned by particle intensity

bins = ggplot(subset(msd_bin, TAU > 0 & TAU < 25),aes(x=TAU,y=MSD,group=factor(BIN),color=BIN)) + geom_point(size=2) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)")
bins = bins + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("1","3","5","7","9"),breaks=c(1,3,5,7,9),option='D',direction=-1)
bins = bins + ggtitle("Ensemble averaged MSD - Partitioned by intensity") + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10),labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_x_log10(breaks = c(0.1,1,10,100),labels = scales::trans_format("log10", scales::math_format(10^.x)))
bins = bins + annotation_logticks() + geom_smooth(method = "lm",se = FALSE)
bins$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_bins.png",sep=""), plot=bins, width = 11, height = 8.5, dpi=300)

###########
## PLOT ENSEMBLE averaged MSD (all particles)

ens = ggplot(msd) + geom_line(aes(x=TAU,y=MSD)) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)") + ggtitle("Ensemble average MSD")
ens$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average.png",sep=""), plot=ens, width = 11, height = 8.5, dpi=300)

###########
## PLOT TIME+ENSEMBLE averaged MSD (all particles)

time_ens = ggplot(msd_tau, aes(x=TAU*frame_int,y=VALUE)) + geom_line() + theme_Publication() + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau (um²)") + ggtitle("MSD tau (Time-ensemble average)")
ggsave(paste(prefix,"_MSD_time-ensemble_average.png",sep=""), plot=time_ens, width = 11, height = 8.5, dpi=300)

###########
## PLOT TWO PREVIOUS GRAPHS TOGETHER

both = ggplot(both_MSD) + geom_line(aes(x=TAU*frame_int,y=value,group=variable,color=variable),size=0.7) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)") + scale_colour_manual(values=c("firebrick1","dodgerblue"),labels=c("Time", "Ensemble"))
both$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_both.png",sep=""), plot=both, width = 11, height = 8.5, dpi=300)

###########
## PLOT DISTRIBUTION OF DISPLACEMENTS FOR ALL TRAJECTORIES

disp = ggplot(power_law,aes(x=x, y=y)) + geom_line(size=0.8) + scale_x_log10(breaks = c(0.01,0.1,1,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))
disp = disp + scale_y_log10(breaks = c(0.001,0.01,0.1,1,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))
disp = disp + theme_Publication() + annotation_logticks() + xlab("Displacement length (um)") + ylab("Frequency") + ggtitle("Frequency of displacements")
disp$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_displacements.png",sep=""), plot=disp, width = 11, height = 8.5, dpi=300)

###########
## PLOT DISTRIBUTION OF LOG(AVERAGE INTENSITY) for all particles

intens = ggplot(msd_single,aes(x=INTENSITY)) + geom_histogram(fill="blue",alpha=0.8, bins=100) + xlab("Particle intensity") + ylab("Count") + ggtitle("Distribution of particle intensities") + theme_Publication()
intens = intens + scale_x_log10(breaks = c(100,1000,10000,100000,1000000),labels = scales::trans_format("log10", scales::math_format(10^.x))) + annotation_logticks(sides = "b")
intens$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_log_intensity.png",sep=""), plot=intens, width = 11, height = 8.5, dpi=300)


#coef(logfit)[2] gives slope -> alpha
# D is obtained by taking the intercept to the power of 10 (since intercept is log(Da))

# D = 10^coef(logfit)[1]

######
## FIT LINE to first 8 seconds of binned MSDs
logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 8), formula = log10(MSD)~log10(TAU)|BIN)

# the coefficients (slopes) give alpha
bin_coef = as.data.frame(coefficients(logfit)[2])

# also get confidence interval for each
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))

# bind in dataframe
alpha_coef = cbind(bin_coef,bin_CI)
alpha_coef$bin = as.numeric(row.names(alpha_coef))
colnames(alpha_coef) = c("alpha","low","high","bin")

###########
## PLOT ALPHA FOR INCREASING PARTICLE INTENSITY

alpha = ggplot(alpha_coef,aes(x=bin,y=alpha,color=bin)) + geom_point(size=5) + geom_errorbar(aes(ymin=low,ymax=high),width=0.) + theme_Publication() + ylab("α") + xlab("Bin") + ggtitle(prefix)
alpha = alpha + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("1","3","5","7","9"),breaks=c(1,3,5,7,9),option='D',direction=-1)
ggsave(paste(prefix,"_alpha_by_partition.png",sep=""), plot=alpha, width = 11, height = 8.5, dpi=300)

######
## GET DIFFUSION COEFFICIENT from intercept of fit to first 8 seconds of binned MSDs

D_coef = 10^coef(logfit)[1]

# also get confidence interval for each
D_CI = as.data.frame(10^confint(logfit,"(Intercept)"))

# bind in dataframe
D_data = cbind(D_coef,D_CI)
D_data$bin = as.numeric(row.names(D_coef))
colnames(D_data) = c("D","low","high","bin")

###########
## PLOT D FOR INCREASING PARTICLE INTENSITY

diff = ggplot(D_data,aes(x=bin,y=D,color=bin)) + geom_point(size=5) + geom_errorbar(aes(ymin=low,ymax=high),width=0.) + theme_Publication() + ylab("D (um²/s)") + xlab("Bin") + ggtitle(prefix)
diff = diff + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("1","3","5","7","9"),breaks=c(1,3,5,7,9),option='D',direction=-1)
ggsave(paste(prefix,"_D_by_partition.png",sep=""), plot=diff, width = 11, height = 8.5, dpi=300)


lm_eqn_log <- function(df){
    m <- lm(data = subset(df, TAU > 0 & TAU < 8), formula = log10(MSD) ~ log10(TAU));
    eq <- as.expression(bquote(paste(alpha, " = ",.(format(coef(m)[2], digits = 4)), ", " ~ r^2 ~ " = " ,.(format(summary(m)$r.squared, digits = 3)))))
    as.character(as.expression(eq));
}

loglog = ggplot(msd,aes(x=TAU,y=MSD)) + geom_point(color="red") + theme_Publication() + labs(x="Time (s)",y="MSD (um²)") + ggtitle("Log-log ensemble average MSD") + geom_smooth(data=subset(msd, TAU > 0 & TAU < 8),method='lm',se=FALSE,color="black")
loglog = loglog + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_x_log10(limits = c(1,1e2),breaks = c(0.1,1,10,100),labels = scales::trans_format("log10", scales::math_format(10^.x)))
loglog = loglog + annotate("label", x = 13, y = 0.08, label = paste(lm_eqn_log(msd)),parse = TRUE, size = 6, color="red") + annotation_logticks()
loglog$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_log-log.png",sep=""), plot=loglog, width = 11, height = 8.5, dpi=300)

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
  paste(prefix,"_4panel1.png",sep=""),
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
multi1 = multiplot(single, ens, time_ens, both, cols=2)
dev.off()

png(
  paste(prefix,"_4panel2.png",sep=""),
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
multi2 = multiplot(bins, alpha, diff, disp, cols=2)
dev.off()