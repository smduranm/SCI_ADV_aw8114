#######################################################################################################
#######################################################################################################
#
# Code set #4
#
# Informing trait-based ecology by assessing remotely-sensed functional diversity across a 4 broad tropical temperature gradient
# SM Duran, RE Martin, S Diaz, BS Maitner, Y Malhi, N Salinas, A Shenkin, MR Silman, DJ Wieczynski, GP Asner, LP Bentley, VM Savage, BJ Enquist.
# 
# Code to reproduce the main figures and analyses in the manuscript
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#Load necessary packages 
library(ggplot2)
library(grid)
library(gridExtra)
library(ggridges)
library(ggpubr)
library(stringr)

##setting up directory 
#Peru_Plot_Master.data<- read.csv(file="Users/brianjenquist/Github/R/Peru_Analyses/Peru_Gradient_NPP_Merged7.csv")
setwd("C:/Users/smdur/Dropbox/PDF/Carnegie/FD_indices/RS_paper")


##########################################################################
##########################################################################
#PLOTTING FIGURE 2:Variation of functional diversity across elevation
##########################################################################

#Load necessary datasets 

#Read mean values of functional diversity indices per site
sum_stats<-read.csv("data/Summary_stats_FDindices_per_site.csv", header = T)
names(sum_stats)

#Read  functional diversity indices calculated using 4m2 pixel site per site
FD<-read.csv("data/FDindices_4m2_pixel_per_site.csv")
head(FD)

#Subsetting file for functional diversity indices (FRic,FDiv) values for plotting
FRic<-subset(FD,FD_index=="FRic") 
FDiv<-subset(FD,FD_index=="FDiv") 

#code to plot figure 2 and save it in figures folder

png("figures/Fig.2.png", units = "px", width=2400, height=1900, res=300)

panel_2A<- 
  ggplot(FRic, aes(y=as.factor(Elev), x=value, fill=site)) +
  geom_density_ridges(scale = 0.8, alpha=0.8) +
  scale_x_continuous(breaks=c(0.0, 0.5, 1.0)) +
  scale_fill_cyclical(values = c("dodgerblue4"))+
  theme_bw()+
  theme(strip.text = element_text(size=12),
        panel.grid = element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        plot.title =element_text(size=12,colour="black",face="bold")) +
  xlab(label='FRic') +
  ylab("Elevation (m)")+
  labs(title="A") 

panel_2B<- 
  ggplot(FDiv, aes(y=as.factor(Elev), x=value, fill=site)) +
  geom_density_ridges(scale = 0.7, alpha=0.8) +
  scale_y_discrete(expand = c(0.09, 0))+
  scale_x_continuous(breaks=c(0.5,0.7, 0.9))+
  xlab(label='FDiv') +
  scale_fill_cyclical(values = c("dodgerblue4"))+
  theme_bw()+
  theme(strip.text = element_text(size=12),
        panel.grid = element_blank(), 
        axis.text = element_text(size=10),
        axis.title.y = element_blank(),
        plot.title =element_text(size=12,colour="black",face="bold"),
        axis.title.x = element_text(size=12)) +
    labs(title="B") 


panel_2C <- 
  ggplot(sum_stats, aes(x=Elev, y=FRic)) +
  geom_point(alpha=0.8, size=3,col="black") + 
  geom_text(aes(x = 3000, y = 0.70, label = "R^2 == 0.53"), parse = TRUE)+
  scale_y_continuous(limits=c(0, 1), breaks = c(0.0, 0.5, 1.0))+
  geom_errorbar(aes(ymin=FRicMeanLower, ymax=FRicMeanUpper)) +
  geom_smooth(data=Mean_values, method="lm",formula = y~x, se=T)+
  labs(title="C") +
  xlab(label="Elevation (m)")+
  ylab(label='FRic')  +
  theme_bw()+
  theme(axis.text = element_text(size=10),
        plot.title =element_text(size=12,colour="black",face="bold"),
        axis.title = element_text(size=12))
  
panel_2D <- 
  ggplot(sum_stats, aes(x=Elev, y=FDiv)) +
  geom_point(alpha=0.8, size=3,col="black") + 
    geom_errorbar(aes(ymin=FDivSELow, ymax=FDivSEUp)) +
  scale_y_continuous(limits=c(0.5, 0.9), breaks = c(0.5, 0.7, 0.9))+
  ylab(label='FDiv')  +
  theme_bw()+
  theme(axis.text = element_text(size=10),
        plot.title =element_text(size=12,colour="black", face="bold"),
        axis.title = element_text(size=12))+
  xlab(label="Elevation (m)")+
  labs(title="D") 


grid.arrange(panel_2A,panel_2B,panel_2C,panel_2D, nrow=2)

dev.off()

###########################################################
###########################################################
#PLOTTING FIGURE 3: FUNCTIONAL DIVERSITY AREA CURVES
##########################################################

#The dataset below contain mean values of FRic and FDiv per area calculate after running two functions (spatial & non-spatial) to calculate
#functional diversity area curves in each site. Check FDAR_curves_from_raster.R and associated files for specifics in each function and model
#the file obtains values of convex hull volume calculations per site, which were standardized from 0 to 1 (column FRic)
#Area in hectares was calculated by multiplying the number of pixels by pixel size (4m2)

#Load necessary datasets, stats summary of FRic and FDiv per area
FDAR.values<-read.csv("data/fdar_values_per_site.csv")
head(FDAR.values)

#Subsetting data per site for plotting. Need to be done as each site has a total number of pixels (e.g., total area is variable across sites)

ACJ_01<-subset(FDAR.values,site=="ACJ_01") 
ESP_01<-subset(FDAR.values,site=="ESP_01") 
PAN_02<-subset(FDAR.values,site=="PAN_02") 
SPD_01<-subset(FDAR.values,site=="SPD_01") 
SPD_02<-subset(FDAR.values,site=="SPD_02") 
TAM_05<-subset(FDAR.values,site=="TAM_05") 
TAM_06<-subset(FDAR.values,site=="TAM_06") 
TRU_04<-subset(FDAR.values,site=="TRU_04") 
WAY_01<-subset(FDAR.values,site=="WAY_01") 

#Code to develop Figure A and save it to figures folder
#Plotting panel A (Figure 3A) first: functional richness (FRic) vs. area per site

png("figures/Figure_3A.png", units = "px", width=2400, height=2200, res=370)

panel_1<-
  ggplot(TAM_06, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 3.6), breaks = c(0,3.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ('green', 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") 

panel_2<-
  ggplot(TAM_05, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0,4.6), breaks = c(0,4.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank())+ 
        facet_wrap(~Elev, scales="free_x") 
 
panel_3<-
  ggplot(PAN_02, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 5.55), breaks = c(0,5.4))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") 

panel_4<-
  ggplot(SPD_02, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 2.6), breaks = c(0,2.5))+
  scale_color_manual(values=c ("green", 'blue'))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") 

panel_5<-
  ggplot(SPD_01, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  #expand_limits(x = 0, y = 0)+
  scale_x_continuous(expand = c(0, 0),limits = c(0, 1.7), breaks = c(0,1.6))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") 

panel_6<-
  ggplot(TRU_04, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=(y) ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 2.3), breaks = c(0,2.2))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position = c(0.6, 0.3),
        legend.title=element_blank(),
        legend.text=element_text(size=9.5,colour="black"),
        panel.grid = element_blank(), 
        axis.title = element_blank(),
        strip.text = element_text(size=12)) +
  facet_wrap(~Elev, scales="free_x") 

panel_7<-
  ggplot(ESP_01, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 1.27), breaks = c(0,1.2))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.x = element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") +
  xlab("Area (ha)") 

panel_8<-
  ggplot(WAY_01, aes(x = area_ha, y =  FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 3.6), breaks = c(0,3.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") +
  xlab("Area (ha)") 

panel_9<-
  ggplot(ACJ_01, aes(x = area_ha, y = FRic)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 0.52), breaks = c(0,0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1), breaks = c(0,1))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        legend.position="none",
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales="free_x") +
  xlab("Area (ha)")

grid.arrange(panel_1,panel_2,panel_3,panel_4,panel_5,panel_6,panel_7,panel_8,panel_9, 
             top= textGrob("A", gp = gpar(fontface = "bold", cex = 1.5), hjust=1, x=0.05), ncol=3)

dev.off()

#Plotting panel B (Figure 3B) second: functional divergence (FDiv) vs. area per site

png("figures/Figure_3B.png", units = "px", width=2400, height=2200, res=370)

panel_1<-
  ggplot(TAM_06, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 3.7), breaks = c(0,3.6))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.74), breaks = c(0.71, 0.74))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") 
 

panel_2<-
  ggplot(TAM_05, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0,4.7), breaks = c(0,4.6))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.74), breaks = c(0.71, 0.74))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") 

panel_3<-
  ggplot(PAN_02, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 5.5), breaks = c(0,5.4))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.75), breaks = c(0.71, 0.75))+
  scale_color_manual(values=c ('green', 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.grid.minor = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") 

panel_4<-
  ggplot(SPD_02, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 2.6), breaks = c(0,2.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.75), breaks = c(0.71, 0.75))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") 

panel_5<-
  ggplot(SPD_01, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 1.62), breaks = c(0,1.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.7, 0.74), breaks = c(0.70, 0.74))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") 
 
panel_6<-
  ggplot(TRU_04, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=(y) ~ log(x), se=T, size=0.7, aes(color=Model)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 2.23), breaks = c(0,2.1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.76), breaks = c(0.71, 0.76))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        strip.text = element_text(size=12),
        legend.title=element_blank(),
        legend.position = c(0.55,0.30),
        legend.text=element_text(size=9.5,colour="black")) +
  facet_wrap(~Elev, scales = "free_y")
  
panel_7<-
  ggplot(ESP_01, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 1.25), breaks = c(0,1.2))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.75), breaks = c(0.71, 0.75))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        strip.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        axis.title = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") +
  xlab("Area (ha)")

panel_8<-
  ggplot(WAY_01, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 3.5), breaks = c(0,3.4))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.75), breaks = c(0.71, 0.75))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(legend.position="none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        strip.text = element_text(size=12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  facet_wrap(~Elev, scales = "free_y") +
  xlab("Area (ha)") 

panel_9<-
  ggplot(ACJ_01, aes(x = area_ha, y = FDiv)) +
  geom_smooth(method ='lm', formula=y ~ log(x), se=T, size=0.7, aes(color=fdar_type)) +
  scale_x_continuous(expand = c(0, 0),limits = c(0, 0.565), breaks = c(0,0.5))+
  scale_y_continuous(expand = c(0, 0), limits = c(0.71, 0.75), breaks = c(0.71, 0.75))+
  scale_color_manual(values=c ("green", 'blue'))+
  theme_bw()  +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        legend.position="none",
        strip.text = element_text(size=12),
        panel.grid = element_blank()) +
  facet_wrap(~Elev,scales = "free_y") +
  xlab("Area (ha)") 

grid.arrange(panel_1,panel_2,panel_3,panel_4,panel_5,panel_6,panel_7,panel_8,panel_9,  
             top= textGrob("B", gp = gpar(fontface = "bold", cex = 1.5), hjust=1, x=0.05), ncol=3)

dev.off()

#####################################################################################################################
#####################################################################################################################
#PLOTTING FIGURE 4: Changes in functional diversity indices along elevation after controlling for scale dependency
####################################################################################################################

#Figures of the slope of functional diversity vs. area were developed by fitting logarithmic models to the FRic-Area and FDiv-area relationships
#Standardized effect sizes (SES) were calculated using the equation by Gotelli & McCabe (Ref #52). 
#Mean values, standard deviation and 95% confidence intervals of SES were calculated from 1000 bootstrapped replicate samples.

#Load necessary files
sum_stats<-read.csv("data/Summary_stats_FDindices_per_site.csv", header = T)
names(sum_stats)

png("figures/Figure_4.png", units = "px", width=2400, height=1800, res=350)

panel_A <- 
  ggplot(sum_stats, aes(x=(Elev), y=(SlopeFRicArea))) +
  geom_point(alpha=0.8, size=3,col="black") + 
  geom_errorbar(aes(ymin=FRicAreaLower, ymax=FRicAreaUpper)) +
  ggtitle("A") +
  xlim(200,3600)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        plot.title =element_text(size=12,colour="black",face="bold"),
        strip.text = element_text(size=11))+
  ylab(label= "Slope FRic vs Area") +
  xlab(label='Elevation (m)')

panel_B <- 
  ggplot(sum_stats, aes(x=(Elev), y=(SlopeFDivArea))) +
  geom_point(alpha=0.8, size=3,col="black") + 
  geom_errorbar(aes(ymin=FDivAreaLower, ymax=FDivAreaUpper)) +
  ggtitle("B") +
  xlim(200,3600)+
  ylim(0.48,0.51)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        plot.title =element_text(size=12,colour="black",face="bold"),
        strip.text = element_text(size=11))+
  ylab(label= "Slope FDiv vs Area") +
  xlab(label='Elevation (m)')

panel_C <- 
  ggplot(sum_stats, aes(x = Elev, y = SES_FRic)) +
  geom_point(alpha=0.8, shape=21,size=3, stat="identity", position="identity",aes(colour = factor(Wilcoxon), fill=factor(Wilcoxon)))+
  scale_fill_manual(values = c("white", "black"))+
  scale_colour_manual(values = c("black", "black"))+
  geom_errorbar(aes(ymin=SESFRicLower, ymax=SESFRicUpper)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1)+
  ggtitle("C")+
  xlim(200,3600)+
  ylim(-8,0.35)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        legend.position="none",
        legend.title=element_blank(),
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size=12,colour="black",face="bold"),
        strip.text = element_text(size=11))+
  ylab(label= "SES - FRic") +
  xlab(label='Elevation (m)')

panel_D <- 
  ggplot(sum_stats, aes(x = Elev, y = SES_FDiv)) +
  geom_point(alpha=0.8, shape=16,size=3, stat="identity", position="identity")+
  geom_errorbar(aes(ymin=SESFDivLower, ymax=SESFDivUpper)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=1)+
  ggtitle("D")+
  xlim(200,3600)+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        legend.position="none",
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        plot.title =element_text(size=12,colour="black",face="bold"),
        strip.text = element_text(size=11))+
  ylab(label= "SES - FDiv") +
  xlab(label='Elevation (m)')

ggarrange(panel_A, panel_B, panel_C, panel_D, ncol = 2, nrow = 2)

dev.off()

#############################################################################################
#############################################################################################
#PLOTTING FIGURE 5: Linkages between forest productivity and functional diversity indices
#############################################################################################



#Load necessary files
sum_stats<-read.csv("data/Summary_stats_FDindices_per_site.csv", header = T)
names(sum_stats)


png("figures/Figure_5.png", units = "px", width=2400, height=1800, res=300)

p1 <- 
  ggplot(sum_stats, aes(x=CMT_LMA, y=GPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  scale_y_continuous(breaks=c(20, 30, 40)) +
  scale_x_continuous(breaks=c(110, 130, 150)) +
  geom_errorbar(aes(ymin=GPPLower, ymax=GPPUpper), width=0.1)+
  xlab(label='LMA ('*g ~ m^-2*')') +
  ylab(label='GPP ('*Mg~C~y^-1~ha^-1*')') +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11)) +
  geom_text(aes(x = 147, y = 41, label = "italic(r)^2 == 0.63"),  size=3,parse = TRUE)

p2 <- 
  ggplot(sum_stats, aes(x=CMT_NSC, y=GPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  ylab(label='GPP ('*Mg~C~y^-1~ha^-1*')')+
  xlab(label='NSC (%)')+
  scale_y_continuous(breaks=c(20, 30, 40)) +
  scale_x_continuous(breaks=c(45, 60, 75)) +
  geom_errorbar(aes(ymin=GPPLower, ymax=GPPUpper), width=0.1)+
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11)) +
  geom_text(aes(x = 68, y = 41,label = "italic(r)^2 == 0.54"),  size=3, parse = TRUE)

p3 <- 
  ggplot(sum_stats, aes(x=FRic, y=GPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  geom_errorbar(aes(ymin=GPPLower, ymax=GPPUpper), width=0.01)+
  ylab(label='GPP ('*Mg~C~y^-1~ha^-1*')')+
  xlab(label='FRic')  +
  scale_y_continuous(breaks=c(20, 30, 40)) +
  scale_x_continuous(breaks=c(0.1,0.3,0.5)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11)) +
  geom_text(aes(x = 0.18, y = 43, label = "italic(r)^2 == 0.58"),  size=3, parse = TRUE)


p4 <- 
  ggplot(sum_stats, aes(x=CMT_Chloro, y=GPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  geom_errorbar(aes(ymin=GPPLower, ymax=GPPUpper), width=0.1)+
  ylab(label='GPP ('*Mg~C~y^-1~ha^-1*')')+
  xlab(label='Chl ('*mg ~ g^-1*')') +
  scale_y_continuous(breaks=c(20, 30, 40)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11)) +
  geom_text(aes(x = 4.2, y = 40,label = "italic(r)^2 == 0.35"), size=3, parse = TRUE)


p5 <- 
  ggplot(sum_stats, aes(x=CMT_Water, y=GPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue" ) + 
  geom_errorbar(aes(ymin=GPPLower, ymax=GPPUpper), width=0.1)+
  ylab(label='GPP ('*Mg~C~y^-1~ha^-1*')') +
  xlab(label='Water (%)')+
  scale_y_continuous(breaks=c(20, 30, 40)) +
  scale_x_continuous(breaks=c(50, 60)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11))

p6 <- 
  ggplot(sum_stats, aes(x=FDiv, y=GPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  geom_errorbar(aes(ymin=GPPLower, ymax=GPPUpper), width=.001)+
  xlab(label='FDiv')  +
  ylab(label='GPP ('*Mg~C~y^-1~ha^-1*')')+
  scale_y_continuous(breaks=c(20, 30, 40)) +
  scale_x_continuous(breaks=c(0.73, 0.75)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11))

p7 <- 
  ggplot(sum_stats, aes(x=CMT_LMA, y=NPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  geom_errorbar(aes(ymin=NPPLower, ymax=NPPUpper), width=0.45)+
  xlab(label='LMA ('*g ~ m^-2*')') +
  ylab(label='NPP ('*Mg~C~y^-1~ha^-1*')')+
  scale_y_continuous(breaks=c(5, 10, 15)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11)) +
  geom_text(aes(x = 147, y = 14.5, label = "italic(r)^2 == 0.61"),  size=3, parse = TRUE)

p8 <- 
  ggplot(sum_stats, aes(x=CMT_NSC, y=NPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  geom_errorbar(aes(ymin=NPPLower, ymax=NPPUpper), width=0.3)+
  scale_y_continuous(breaks=c(5, 10, 15)) +
  scale_x_continuous(breaks=c(45, 60, 75)) +
  xlab(label='NSC (%)') +
  ylab(label='NPP ('*Mg~C~y^-1~ha^-1*')')+
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11)) +
  geom_text(aes(x = 68, y = 14.5, label = "italic(r)^2 == 0.37"), size=3, parse = TRUE)

p9 <- 
  ggplot(sum_stats, aes(x=FRic, y=NPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  geom_errorbar(aes(ymin=NPPLower, ymax=NPPUpper), width=0.01)+
  scale_y_continuous(breaks=c(5, 10, 15)) +
  scale_x_continuous(breaks=c(0.1,0.3,0.5)) +
  ylab(label='NPP ('*Mg~C~y^-1~ha^-1*')')+
  xlab(label='FRic')  +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11)) +
  geom_text(aes(x = 0.19, y = 15.7, label = "italic(r)^2 == 0.47"), size=3, parse = TRUE)

p10 <- 
  ggplot(sum_stats, aes(x=CMT_Chloro, y=NPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  stat_smooth(method="lm", col = "black", se=T) +
  geom_errorbar(aes(ymin=NPPLower, ymax=NPPUpper), width=0.09)+  
  ylab(label='NPP ('*Mg~C~y^-1~ha^-1*')')+
  xlab(label='Chl ('*mg ~ g^-1*')') +
  scale_y_continuous(breaks=c(5, 10, 15)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11)) +
  geom_text(aes(x = 4.2, y = 14.5, label = "italic(r)^2 == 0.44"), size=3, parse = TRUE)

p11 <- 
  ggplot(sum_stats, aes(x=CMT_Water, y=NPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue" ) + 
  geom_errorbar(aes(ymin=NPPLower, ymax=NPPUpper), width=0.1)+
  ylab(label='NPP ('*Mg~C~y^-1~ha^-1*')') +
  xlab(label='Water (%)')+
  scale_y_continuous(breaks=c(5, 10, 15)) +
  scale_x_continuous(breaks=c(50, 60)) +
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 11))


p12 <- 
  ggplot(sum_stats, aes(x=FDiv, y=NPP)) +
  geom_point(alpha=0.8, size=2.5,col="blue") + 
  geom_errorbar(aes(ymin=NPPLower, ymax=NPPUpper), width=0.001)+
  scale_y_continuous(breaks=c(5, 10, 15))+
  scale_x_continuous(breaks=c(0.73, 0.75)) +
  xlab(label='FDiv')  +
  ylab(label='NPP ('*Mg~C~y^-1~ha^-1*')')+
  theme_bw()+
  theme(axis.text = element_text(size = 9),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 11))


grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6, top="         Gross primary productivity (GPP)"),
             arrangeGrob(p7,p8,p9,p10,p11,p12, top="      Net primary productivity (NPP)"), ncol=2)

dev.off()





