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

both_MSD = rbind(msd_tau,msd)

msd$TAU = msd$TAU*frame_int
msd_bin$TAU = msd_bin$TAU*frame_int

###########
## GET MEAN INTENSITY PER FRAME

mean_intensity = aggregate(msd_single$INTENSITY, by=list(frame=msd_single$TAU), mean)
sd_intensity = aggregate(msd_single$INTENSITY, by=list(frame=msd_single$TAU), sd)
intensity_profile = cbind(mean_intensity,sd_intensity)
intensity_profile = intensity_profile[,-3]
avg_intens = ggplot(intensity_profile) + geom_point(aes(x=frame,y=x),size=2) + theme_Publication() + ylab("average intensity of particles")
avg_vs_stdev = ggplot(intensity_profile) + geom_point(aes(x=x,y=x.1),size=2) + theme_Publication() + ylab("standard deviation") + xlab("mean")

###########

displacements <- as.data.frame(read.csv(paste(prefix,"_displacements.txt",sep=""),sep="\t",header = FALSE))
displacements$V1 = displacements$V1 + 1
displacements$V1 = as.factor(displacements$V1)
colnames(displacements)[1] = "Bin"
h=hist(displacements$V2,breaks = 300,plot=FALSE)
x=h$mids
y=h$density
power_law = cbind(x,y)
power_law = as.data.frame(power_law)

###########
## PLOT THE MSD OF EACH TRAJECTORY SEPARATELY

message("R: Plotting single trajectories ...")

single = ggplot(msd_single, aes(x=as.numeric(factor(TAU))*frame_int,y=MSD,color=factor(TRACK))) + geom_line() + theme_Publication() + scale_color_discrete(guide=FALSE) + labs(x = "Time (s)",y="MSD tau (um²)") + ggtitle("MSD tau (All particles)")
single$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_single_trajectories.png",sep=""), plot=single, width = 11, height = 8.5, dpi=300)

###########
## PLOT ENSEMBLE averaged MSD binned by particle intensity

message("R: Plotting binned MSD (log-log) ...")

logbins = ggplot(subset(msd_bin, TAU > 0 & TAU < 25),aes(x=TAU,y=MSD,group=factor(BIN),color=BIN)) + geom_point(size=2) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)")
logbins = logbins + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("1","3","5","7","9"),breaks=c(1,3,5,7,9),option='D',direction=-1)
logbins = logbins + ggtitle("Ensemble averaged MSD\nPartitioned by intensity") + scale_y_log10(breaks = c(0.0001,0.001,0.01,0.1,1,10),labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_x_log10(breaks = c(0.1,1,10,100),labels = scales::trans_format("log10", scales::math_format(10^.x)))
logbins = logbins + annotation_logticks() + geom_smooth(data=subset(msd_bin, TAU < 3.75),method = "lm",se = FALSE)
logbins$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_log_bins.png",sep=""), plot=logbins, width = 11, height = 8.5, dpi=300)


message("R: Plotting binned MSD (linear) ...")

bins = ggplot(msd_bin,aes(x=TAU,y=MSD,group=factor(BIN),color=BIN)) + geom_line(size=2) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)")
bins = bins + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("1","3","5","7","9"),breaks=c(1,3,5,7,9),option='D',direction=-1)
bins = bins + ggtitle("Ensemble averaged MSD\nPartitioned by intensity")
bins$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_bins.png",sep=""), plot=bins, width = 11, height = 8.5, dpi=300)

###########
## PLOT ENSEMBLE averaged MSD (all particles)

message("R: Plotting ensemble averaged MSD ...")

ens = ggplot(msd,aes(x=TAU,y=MSD)) + geom_line(size=2) + geom_errorbar(aes(ymin=MSD-SE, ymax=MSD+SE),width=0.1) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)") + ggtitle("Ensemble average MSD")
ens$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_ensemble_average.png",sep=""), plot=ens, width = 11, height = 8.5, dpi=300)

###########
## PLOT TIME+ENSEMBLE averaged MSD (all particles)

message("R: Plotting time-ensemble averaged MSD ...")

time_ens = ggplot(msd_tau, aes(x=TAU*frame_int,y=MSD)) + geom_line(size=2) + geom_errorbar(aes(ymin=MSD-SE, ymax=MSD+SE),width=0.1) + theme_Publication() + scale_color_discrete(guide=FALSE) + labs(x="Time (s)",y="MSD tau (um²)") + ggtitle("MSD tau (Time-ensemble average)")
ggsave(paste(prefix,"_MSD_time-ensemble_average.png",sep=""), plot=time_ens, width = 11, height = 8.5, dpi=300)

###########
## PLOT TWO PREVIOUS GRAPHS TOGETHER

message("R: Plotting overlay MSD ...")

overlay = ggplot(both) + geom_line(aes(x=TAU*frame_int,y=MSD,group=type,color=type),size=0.7) + theme_Publication() + labs(x="Time (s)",y="MSD (um²)") + scale_colour_manual(values=c("firebrick1","dodgerblue"))
overlay$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_MSD_both.png",sep=""), plot=overlay, width = 11, height = 8.5, dpi=300)

###########
## PLOT POWER LAW OF DISPLACEMENTS FOR ALL TRAJECTORIES

message("R: Plotting log-log distribution of single step displacements ...")

norm=hist(rnorm(n = 1000000,mean = 0.005,sd = 0.2),breaks = 1000,plot=FALSE)
normx=norm$mids
normy=norm$density
norm_law = cbind(normx,normy)
norm_law = as.data.frame(norm_law)

power = ggplot(power_law,aes(x=x, y=y)) + geom_line(size=0.8) + scale_x_log10(breaks = c(0.01,0.1,1,10),labels = scales::trans_format("log10", scales::math_format(10^.x)))
power = power + scale_y_log10(breaks = c(0.001,0.01,0.1,1,10),labels = scales::trans_format("log10", scales::math_format(10^.x))) + geom_line(data=norm_law,aes(x=normx,y=normy),color="red",size=0.8)
power = power + theme_Publication() + annotation_logticks() + xlab("Displacement length (um)") + ylab("Frequency") + ggtitle("Frequency of 1-frame displacements")
power$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_displacements_log.png",sep=""), plot=power, width = 11, height = 8.5, dpi=300)

###########
## PLOT DISTRIBUTION OF DISPLACEMENTS FOR ALL TRAJECTORIES

message("R: Plotting histogram of single step displacements ...")

disp = ggplot(displacements,aes(x=log10(V2))) + geom_histogram(bins = 300,fill="#E01A63") + theme_Publication()
disp = disp + xlab("Log(Displacement (um))") + ylab("Count") + ggtitle("Distribution of displacements for a single step")
disp$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_displacements_distribution.png",sep=""), plot=disp, width = 11, height = 8.5, dpi=300)

###########
## PLOT DISTRIBUTION OF DISPLACEMENTS BY BIN

message("R: Plotting displacements by bin...")

bindisp = ggplot(displacements,aes(x=log10(V2))) + geom_density(aes(y = ..count..,group=Bin,fill=Bin),position="stack") + theme_Publication() + scale_fill_viridis(discrete = TRUE)
bindisp = bindisp + xlab("Log(Displacement (um))") + ylab("Count") + ggtitle("Cumulative distribution of\n1-step displacements")
bindisp$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_displacements_binned.png",sep=""), plot=bindisp, width = 11, height = 8.5, dpi=300)

###########
## PLOT DISTRIBUTION OF LOG(AVERAGE INTENSITY) for all particles

message("R: Plotting distribution of particle intensities ...")

intens = ggplot(msd_single) + geom_histogram(aes(x=log(INTENSITY),fill=as.factor(BIN)), bins=100) + scale_fill_viridis(discrete = TRUE,direction=-1,option="D") + xlab("Log(normalized particle intensity)") + ylab("Count") + ggtitle("Distribution of particle intensities") + theme_Publication()
intens$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
ggsave(paste(prefix,"_log_intensity.png",sep=""), plot=intens, width = 11, height = 8.5, dpi=300)


#coef(logfit)[2] gives slope -> alpha
# D is obtained by taking the intercept to the power of 10 (since intercept is log(Da))

# D = 10^coef(logfit)[1]

message("R: Fitting log10(MSD)~log10(TAU) ...")

logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 3), formula = log10(MSD)~log10(TAU)|BIN)
bin_coef = as.data.frame(coefficients(logfit)[2])
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))
bin_R2 = as.data.frame(summary(logfit)$r.squared)
alpha_coef = cbind(bin_coef,bin_CI,bin_R2)
alpha_coef$bin = as.numeric(row.names(alpha_coef))
colnames(alpha_coef) = c("alpha","low","high","R2","bin")
alpha_coef$n_fit = 3

logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 3.75), formula = log10(MSD)~log10(TAU)|BIN)
bin_coef = as.data.frame(coefficients(logfit)[2])
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))
bin_R2 = as.data.frame(summary(logfit)$r.squared)
temp = cbind(bin_coef,bin_CI,bin_R2)
temp$bin = as.numeric(row.names(temp))
colnames(temp) = c("alpha","low","high","R2","bin")
temp$n_fit = 4
alpha_coef = rbind(alpha_coef,temp)

logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 4), formula = log10(MSD)~log10(TAU)|BIN)
bin_coef = as.data.frame(coefficients(logfit)[2])
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))
bin_R2 = as.data.frame(summary(logfit)$r.squared)
temp = cbind(bin_coef,bin_CI,bin_R2)
temp$bin = as.numeric(row.names(temp))
colnames(temp) = c("alpha","low","high","R2","bin")
temp$n_fit = 5
alpha_coef = rbind(alpha_coef,temp)

logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 5), formula = log10(MSD)~log10(TAU)|BIN)
bin_coef = as.data.frame(coefficients(logfit)[2])
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))
bin_R2 = as.data.frame(summary(logfit)$r.squared)
temp = cbind(bin_coef,bin_CI,bin_R2)
temp$bin = as.numeric(row.names(temp))
colnames(temp) = c("alpha","low","high","R2","bin")
temp$n_fit = 6
alpha_coef = rbind(alpha_coef,temp)

logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 5.5), formula = log10(MSD)~log10(TAU)|BIN)
bin_coef = as.data.frame(coefficients(logfit)[2])
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))
bin_R2 = as.data.frame(summary(logfit)$r.squared)
temp = cbind(bin_coef,bin_CI,bin_R2)
temp$bin = as.numeric(row.names(temp))
colnames(temp) = c("alpha","low","high","R2","bin")
temp$n_fit = 7
alpha_coef = rbind(alpha_coef,temp)

logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 6.5), formula = log10(MSD)~log10(TAU)|BIN)
bin_coef = as.data.frame(coefficients(logfit)[2])
bin_CI = as.data.frame(confint(logfit,"log10(TAU)"))
bin_R2 = as.data.frame(summary(logfit)$r.squared)
temp = cbind(bin_coef,bin_CI,bin_R2)
temp$bin = as.numeric(row.names(temp))
colnames(temp) = c("alpha","low","high","R2","bin")
temp$n_fit = 8
alpha_coef = rbind(alpha_coef,temp)

message("R: Plotting alpha vs number of points fitted ...")

alpha_fit = ggplot(alpha_coef,aes(x=bin,y=alpha,color=n_fit)) + geom_jitter(aes(size=R2),position = "dodge") + theme_Publication() + ylab("α") + xlab("Bin") + ggtitle("Alpha exponent by bin\nby number of fitted points")
alpha_fit = alpha_fit + scale_color_viridis(guide=guide_legend(title="# points fitted",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),option='C',direction=-1)
ggsave(paste(prefix,"_alpha_by_fit.png",sep=""), plot=alpha_fit, width = 11, height = 8.5, dpi=300)

message("R: Plotting R-squared with number of points fitted ...")

corr = ggplot(alpha_coef,aes(x=n_fit,y=R2,color=bin)) + geom_jitter(size=4,position = "dodge") + theme_Publication() + ggtitle("R-squared of fit\nby number of points used") + ylab("R-squared") + xlab("Number of fitted points") + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),option='D',direction=-1)
ggsave(paste(prefix,"_R2_by_fit.png",sep=""), plot=corr, width = 11, height = 8.5, dpi=300)

######
## FIT LINE to first 8 seconds of binned MSDs
logfit <- lmList(data = subset(msd_bin, TAU > 0 & TAU < 3.75), formula = log10(MSD)~log10(TAU)|BIN)

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

message("R: Plotting alpha by bin ...")

alpha = ggplot(alpha_coef,aes(x=bin,y=alpha,color=bin)) + geom_point(size=5) + geom_errorbar(aes(ymin=low,ymax=high),width=0.) + theme_Publication() + ylab("α") + xlab("Bin") + ggtitle("Alpha exponent by bin")
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

message("R: Plotting diffusion coefficient by bin ...")

diff = ggplot(D_data,aes(x=bin,y=D,color=bin)) + geom_point(size=5) + geom_errorbar(aes(ymin=low,ymax=high),width=0.) + theme_Publication() + ylab("D (um²/s)") + xlab("Bin") + ggtitle("Diffusion coefficient by bin")
diff = diff + scale_color_viridis(guide=guide_legend(title="Bin",title.theme = element_text(size = 8,face = "bold",colour = "black",angle = 0),label.theme = element_text(size = 8,colour = "black",angle = 0),keywidth = 0.1, keyheight = 0.1,default.unit="inch",override.aes = list(alpha = 1),byrow = TRUE),labels=c("1","3","5","7","9"),breaks=c(1,3,5,7,9),option='D',direction=-1)
ggsave(paste(prefix,"_D_by_partition.png",sep=""), plot=diff, width = 11, height = 8.5, dpi=300)


lm_eqn_log <- function(df){
    m <- lm(data = subset(df, TAU > 0 & TAU < 8), formula = log10(MSD) ~ log10(TAU));
    eq <- as.expression(bquote(paste(alpha, " = ",.(format(coef(m)[2], digits = 4)), ", " ~ r^2 ~ " = " ,.(format(summary(m)$r.squared, digits = 3)))))
    as.character(as.expression(eq));
}

# message("R: Plotting log-log ensemble averaged MSD ...")

# loglog = ggplot(msd,aes(x=TAU,y=MSD)) + geom_point(color="red") + theme_Publication() + labs(x="Time (s)",y="MSD (um²)") + ggtitle("Log-log ensemble average MSD") + geom_smooth(data=subset(msd, TAU > 0 & TAU < 8),method='lm',se=FALSE,color="black")
# loglog = loglog + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),labels = scales::trans_format("log10", scales::math_format(10^.x))) + scale_x_log10(limits = c(1,1e2),breaks = c(0.1,1,10,100),labels = scales::trans_format("log10", scales::math_format(10^.x)))
# loglog = loglog + annotate("label", x = 13, y = 0.08, label = paste(lm_eqn_log(msd)),parse = TRUE, size = 6, color="red") + annotation_logticks()
# loglog$theme$plot.margin = unit(c(0.5,1,0.5,0.5),"cm")
# ggsave(paste(prefix,"_MSD_log-log.png",sep=""), plot=loglog, width = 11, height = 8.5, dpi=300)

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

message("R: Plotting grids ...")

## rescale text for multiplot
single = single + theme_Publication(base_size=16)
ens = ens + theme_Publication(base_size=16)
time_ens = time_ens + theme_Publication(base_size=16)
overlay = overlay + theme_Publication(base_size=16)
logbins = logbins + theme_Publication(base_size=16)
bins = bins + theme_Publication(base_size=16)
bindisp = bindisp + theme_Publication(base_size=16)
alpha = alpha + theme_Publication(base_size=16)
diff = diff + theme_Publication(base_size=16)
power = power + theme_Publication(base_size=16)
alpha_fit = alpha_fit + theme_Publication(base_size=16)
corr = corr + theme_Publication(base_size=16)

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
multi1 = multiplot(single, ens, time_ens, overlay, cols=2)
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
multi2 = multiplot(logbins, bins, bindisp, power, cols=2)
dev.off()

png(
  paste(prefix,"_4panel3.png",sep=""),
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
multi3 = multiplot(alpha, diff, alpha_fit, corr, cols=2)
dev.off()