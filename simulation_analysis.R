library(cowplot)
library(grid)
library(gridExtra)
library(gtable)

#kind="nonprog"
kind="prog"
dir=paste0("/Users/Diego/Desktop/singlecrypt/simulation/final2/mcmc/simulation_",kind,"/results/")
#titletext="Nonprogressors (11 leaves)"
titletext="Progressors (54 leaves)"
data=read.csv(paste0(dir,"rf.csv"))
data=data[,-1]
data2=read.csv(paste0(dir,"parameters.out"),sep=" ")

#Reconstituted rates
data$rate=data$rate
data2$r=data2$r ###data2$r is the gain rate. We want it to be the total rate, and all of them are the same, so we just need to multiply by 3
data2$cnv.loss=data2$clock.rate*data2$cnv.loss
data2$cnv.conversion=data2$clock.rate*data2$cnv.conversion
data2$cnv.gain=data2$clock.rate
data2$clock.rate=data2$cnv.loss+data2$cnv.conversion+data2$cnv.gain

#Reconstituted relative errors
data2$cnv.loss.error=abs(data2$cnv.loss-(data2$r/3))/(data2$r/3)
data2$cnv.conversion.error=abs(data2$cnv.conversion-(data2$r/3))/(data2$r/3)
data2$cnv.gain.error=abs(data2$cnv.gain-(data2$r/3))/(data2$r/3)
data2$clock.rate.error=abs(data2$clock.rate-data2$r)/data2$r

data2_rel=data2
#data2_rel$cnv.loss=abs(data2$cnv.loss-1)
#data2_rel$cnv.conversion=abs(data2$cnv.conversion-1)
#data2_rel$clock.rate=abs(data2$clock.rate-data2$r)/data2$r

#Mean bars inspired by aosmith16 @ github https://github.com/tidyverse/ggplot2/issues/1259 

means=aggregate(formula= rf ~ rate+mod+gd,data=data,FUN=mean)

rf_gd0=ggplot(data=data[data$gd==0,],aes(y=rf,x=as.factor(rate)))+
  geom_violin(aes(fill=mod))+
  geom_boxplot(aes(group=interaction(mod,as.factor(rate))),alpha=0,position=position_dodge(0.9),width=0.2)+
  geom_errorbar(data=means[means$gd==0,],aes(y=rf,ymax=rf,ymin=rf,group=interaction(mod,as.factor(rate))),width=0.7,position=position_dodge(0.9),size=1,linetype=1)+
  scale_x_discrete(name="Acquisition rate")+
  #scale_y_continuous(name="RF distance")+
  scale_y_continuous(name="RF distance",limits=c(0,1))+
  labs(title="No genome doubling")+
  scale_fill_discrete(name="Strategy")

rf_gd=ggplot(data=data[data$gd=="0.001",],aes(y=rf,x=as.factor(rate)))+
  geom_violin(aes(fill=mod))+
  geom_boxplot(aes(group=interaction(mod,as.factor(rate))),alpha=0,position=position_dodge(0.9),width=0.2)+
  geom_errorbar(data=means[means$gd=="0.001",],aes(y=rf,ymax=rf,ymin=rf,group=interaction(mod,as.factor(rate))),width=0.7,position=position_dodge(0.9),size=1,linetype=1)+
  scale_x_discrete(name="Acquisition rate")+
  #scale_y_continuous(name="RF distance")+
  scale_y_continuous(name="RF distance",limits=c(0,1))+
  labs(title="Simulated genome doubling")+
  scale_fill_discrete(name="Strategy")

rf_r=ggplot(data=data[data$gd=="random",],aes(y=rf,x=as.factor(rate)))+
  geom_violin(aes(fill=mod))+geom_boxplot(aes(group=interaction(mod,as.factor(rate))),alpha=0,position=position_dodge(0.9),width=0.2)+
  geom_errorbar(data=means[means$gd=="random",],aes(y=rf,ymax=rf,ymin=rf,group=interaction(mod,as.factor(rate))),width=0.7,position=position_dodge(0.9),size=1,linetype=1)+
  scale_x_discrete(name="Acquisition rate")+
  #scale_y_continuous(name="RF distance")+
  scale_y_continuous(name="RF distance",limits=c(0,1))+
  labs(title="Random leaf genome doubling")+
  scale_fill_discrete(name="Strategy")

grid_rf=plot_grid(rf_gd0,rf_gd,rf_r,labels="AUTO",ncol=3)


means=aggregate(formula= clock.rate.error ~ r+mod+gd,data=data2_rel,FUN=mean)

rate_gd0=ggplot(data=data2_rel[data2_rel$gd==0,],aes(y=clock.rate.error,x=as.factor(r)))+
  geom_violin(aes(fill=mod))+geom_boxplot(aes(group=interaction(mod,as.factor(r))),alpha=0,position=position_dodge(0.9),width=0.2)+
  geom_errorbar(data=means[means$gd==0,],aes(y=clock.rate.error,ymax=clock.rate.error,ymin=clock.rate.error,group=interaction(mod,as.factor(r))),width=0.7,position=position_dodge(0.9),size=1,linetype=1)+
  scale_x_discrete(name="Acquisition rate")+
  #scale_y_log10(name="Acquisition rate relative error",limits=c(min(data2_rel$clock.rate.error),max(data2_rel$clock.rate.error)))+
  scale_y_log10(name="Acquisition rate relative error",limits=c(0.0001,400))+
  labs(title="No genome doubling")+
  scale_fill_discrete(name="Strategy")

rate_gd=ggplot(data=data2_rel[data2_rel$gd=="0.001",],aes(y=clock.rate.error,x=as.factor(r)))+
  geom_violin(aes(fill=mod))+geom_boxplot(aes(group=interaction(mod,as.factor(r))),alpha=0,position=position_dodge(0.9),width=0.2)+
  geom_errorbar(data=means[means$gd=="0.001",],aes(y=clock.rate.error,ymax=clock.rate.error,ymin=clock.rate.error,group=interaction(mod,as.factor(r))),width=0.7,position=position_dodge(0.9),size=1,linetype=1)+
  scale_x_discrete(name="Acquisition rate")+
  #scale_y_log10(name="Acquisition rate relative error",limits=c(min(data2_rel$clock.rate.error),max(data2_rel$clock.rate.error)))+
  scale_y_log10(name="Acquisition rate relative error",limits=c(0.0001,400))+
  labs(title="Simulated genome doubling")+
  scale_fill_discrete(name="Strategy")

rate_r=ggplot(data=data2_rel[data2_rel$gd=="random",],aes(y=clock.rate.error,x=as.factor(r)))+
  geom_violin(aes(fill=mod))+geom_boxplot(aes(group=interaction(mod,as.factor(r))),alpha=0,position=position_dodge(0.9),width=0.2)+
  geom_errorbar(data=means[means$gd=="random",],aes(y=clock.rate.error,ymax=clock.rate.error,ymin=clock.rate.error,group=interaction(mod,as.factor(r))),width=0.7,position=position_dodge(0.9),size=1,linetype=1)+
  scale_x_discrete(name="Acquisition rate")+
  scale_y_log10(name="Acquisition rate relative error",limits=c(0.0001,400))+
  #scale_y_log10(name="Acquisition rate relative error",limits=c(min(data2_rel$clock.rate.error),max(data2_rel$clock.rate.error)))+
  labs(title="Random leaf genome doubling")+
  scale_fill_discrete(name="Strategy")

grid_rate=plot_grid(rate_gd0,rate_gd,rate_r,labels=c("D","E","F"),ncol=3)

supergrid=plot_grid(grid_rf,grid_rate,labels="",ncol=1)
supergrid_gtable=cowplot:::ggplot_to_gtable(supergrid)

#By baptiste @ Stackoverflow: http://stackoverflow.com/questions/31640916/how-can-i-add-a-title-to-a-tablegrob-plot
titlegrob = textGrob(titletext,gp=gpar(fontsize=25,fontface="bold"))
padding <- unit(7,"mm")
supergrid_gtable <- gtable_add_rows(supergrid_gtable, 
                         heights = grobHeight(titlegrob) + padding,
                         pos = 0)
supergrid_gtable <- gtable_add_grob(supergrid_gtable, titlegrob, 1, 1, 1, ncol(supergrid_gtable))
##

finalplot=ggdraw(supergrid_gtable)
save_plot(paste0(dir,"panel_",kind,".pdf"),finalplot,base_height=12,base_aspect_ratio=1.5)
save_plot(paste0(dir,"panel_",kind,".png"),finalplot,base_height=12,base_aspect_ratio=1.5)

methods=unique(data$mod)
n_comb=length(methods)*(length(methods)-1)/2 #Bonferroni correction

for (dup in unique(data$gd)) {
  for (rate in unique(data$rate)) {
    for (a in 1:(length(methods)-1)){
      for (b in (a+1):length(methods)){
        datA=data[data$gd==dup & data$rate==rate & data$mod==methods[a],]$rf
        datB=data[data$gd==dup & data$rate==rate & data$mod==methods[b],]$rf
        if(length(datA)>length(datB)) {
          print("Removing samples from datA since there are missing data in datB") ##Samples without Acquisitions after baseline correction
          filt=data[data$gd==dup & data$rate==rate & data$mod==methods[a],]$rep %in% data[data$gd==dup & data$rate==rate & data$mod==methods[b],]$rep
          datA=(data[data$gd==dup & data$rate==rate & data$mod==methods[a],])[filt,]$rf
        }else if (length(datB)>length(datA)) {
          print("Removing samples from datA since there are missing data in datA")
          filt=data[data$gd==dup & data$rate==rate & data$mod==methods[b],]$rep %in% data[data$gd==dup & data$rate==rate & data$mod==methods[a],]$rep
          datB=(data[data$gd==dup & data$rate==rate & data$mod==methods[b],])[filt,]$rf
         }
        
        p_value=wilcox.test(datA,datB,paired = TRUE)$p.value
        add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *",
                                                     ifelse(p_value > 0.001/n_comb, " **", " ***")))
        print(paste(paste("tTop",dup,rate,methods[a],methods[b],sep=":"),paste0(signif(p_value,digits=2),add),mean(datA),sd(datA)/sqrt(length(datA)),mean(datB),sd(datB)/sqrt(length(datB)),sep=","))

        datA=data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[a],]$clock.rate.error
        datB=data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[b],]$clock.rate.error
        if(length(datA)>length(datB)) {
          print("Removing samples from data2_rel since there are missing data in datB") ##Samples without Acquisitions after baseline correction
          filt=data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[a],]$rep %in% data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[b],]$rep
          datA=(data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[a],])[filt,]$clock.rate.error
        }else if (length(datB)>length(datA)) {
          print("Removing samples from data2_rel since there are missing data in datA")
          filt=data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[b],]$rep %in% data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[a],]$rep
          datB=(data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[b],])[filt,]$clock.rate.error
        }
        p_value=wilcox.test(datA,datB,paired = TRUE)$p.value
        add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *",
                                                     ifelse(p_value > 0.001/n_comb, " **", " ***")))
        print(paste(paste("mRate",dup,rate,methods[a],methods[b],sep=":"),paste0(signif(p_value,digits=2),add),mean(datA),sd(datA)/sqrt(length(datA)),mean(datB),sd(datB)/sqrt(length(datB)),sep=","))
      }
    }
  }
} 


n_comb=1
#n_comb=12 #Bonferroni correction

for (dup in c("0.001","random")) {
  for (rate in unique(data$rate)) {
        a=3 ##None
        b=1 #Baseline
        datA=data[data$gd=="0" & data$rate==rate & data$mod==methods[a],]$rf
        datB=data[data$gd==dup & data$rate==rate & data$mod==methods[b],]$rf
        
        p_value=wilcox.test(datA,datB,paired = FALSE)$p.value ##This are not paired
        add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *",
                                                     ifelse(p_value > 0.001/n_comb, " **", " ***")))
        print(paste(paste("tTop",dup,rate,methods[a],methods[b],sep=":"),paste0(signif(p_value,digits=2),add),mean(datA),sd(datA)/sqrt(length(datA)),mean(datB),sd(datB)/sqrt(length(datB)),sep=","))
        
        datA=data2_rel[data2_rel$gd=="0" & data2_rel$r==rate & data2_rel$mod==methods[a],]$clock.rate.error
        datB=data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[b],]$clock.rate.error

        p_value=wilcox.test(datA,datB,paired = FALSE)$p.value
        add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *",
                                                     ifelse(p_value > 0.001/n_comb, " **", " ***")))
        print(paste(paste("mRate",dup,rate,methods[a],methods[b],sep=":"),paste0(signif(p_value,digits=2),add),mean(datA),sd(datA)/sqrt(length(datA)),mean(datB),sd(datB)/sqrt(length(datB)),sep=","))
  }
} 


for (dup in c("0.001","random")) {
  for (rate in unique(data$rate)) {
    a=3 ##None
    b=3 #none
    datA=data[data$gd=="0" & data$rate==rate & data$mod==methods[a],]$rf
    datB=data[data$gd==dup & data$rate==rate & data$mod==methods[b],]$rf
    
    p_value=wilcox.test(datA,datB,paired = FALSE)$p.value ##This are not paired
    add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *",
                                                 ifelse(p_value > 0.001/n_comb, " **", " ***")))
    print(paste(paste("tTop",dup,rate,methods[a],methods[b],sep=":"),paste0(signif(p_value,digits=2),add),mean(datA),sd(datA)/sqrt(length(datA)),mean(datB),sd(datB)/sqrt(length(datB)),sep=","))
    
    datA=data2_rel[data2_rel$gd=="0" & data2_rel$r==rate & data2_rel$mod==methods[a],]$clock.rate.error
    datB=data2_rel[data2_rel$gd==dup & data2_rel$r==rate & data2_rel$mod==methods[b],]$clock.rate.error
    
    p_value=wilcox.test(datA,datB,paired = FALSE)$p.value
    add=ifelse(p_value > 0.05/n_comb, "", ifelse(p_value > 0.01/n_comb, " *",
                                                 ifelse(p_value > 0.001/n_comb, " **", " ***")))
    print(paste(paste("mRate",dup,rate,methods[a],methods[b],sep=":"),paste0(signif(p_value,digits=2),add),mean(datA),sd(datA)/sqrt(length(datA)),mean(datB),sd(datB)/sqrt(length(datB)),sep=","))
  }
} 