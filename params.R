library(cowplot)
library(plyr)
library(birdring)
library(grid)
library(ggrepel)
setwd("/Users/Diego/Desktop/singlecrypt/results/newprior/data")
arrayPat=c("391","437","451","740","848","852")
allPat=c("256","391","437","451","740","848","852","911")
burnin=1000 #10%
#arrayPat=allPat
types=c("crypts","pool")

if (length(types) != 2)
{
  print("The number of types to compare must be two")
  exit()
}
param="clock.rate"
rm(alldata)
results=data.frame(patient=character(0),crypts_mean=numeric(0),biopsies_mean=numeric(0),p_diff=numeric(0),stringsAsFactors = FALSE)
for (npatient in 1:length(arrayPat))
{
  patient=arrayPat[npatient]
  rm(pdata)
  for (ntype in 1:length(types))
  {
    type=types[ntype]
    Sys.sleep(1) ##I get weird random errors if I do not wait for a little while here
    dirs=dir(path=".",pattern=paste0(patient,".*",type),recursive = FALSE)
    replicates=list.files(path=dirs,pattern=paste0(".*",type,"\\.r.*",".log$"),full.names = TRUE)
    prior_file=list.files(path=dirs,pattern=paste0(".*",type,"\\.p\\.log$"),full.names = TRUE)
    #print(paste(type,patient,dirs,replicates))
    rm(tdata)
    priordata=read.csv(prior_file,sep="\t",comment.char="#")
    priordata=priordata[seq(from=burnin+1,to=nrow(priordata),by=1),]
    priordata$prior=rep(1,times=nrow(priordata))
    for (nrep in 1:length(replicates))
    {
      rdata=read.csv(replicates[nrep],sep="\t",comment.char="#")
      rdata=rdata[seq(from=burnin+1,to=nrow(rdata),by=1),]
      if(!exists("tdata"))
      {
        tdata=rdata
      }
      else
      {
        tdata=rbind(tdata,rdata) 
      }
    }
    tdata$prior=rep(0,times=nrow(tdata))
    tdata=rbind(tdata,priordata)
    tdata$type=rep(type,times=nrow(tdata))
    tdata$prior=as.factor(tdata$prior)
    if(!exists("pdata"))
    {
      pdata=tdata
    }
    else
    {
      pdata=rbind(pdata,tdata)
    }
  }
  pdata$patient=rep(patient,times=nrow(pdata))
  p=1-overlap(pdata[pdata$type==types[1] & pdata$prior==0,param],pdata[pdata$type==types[2] & pdata$prior==0,param],from=min(pdata[pdata$prior==0,param]),to=max(pdata[pdata$prior==0,param]),edge.of.parameter.space=TRUE) ##We may want to calculate the overlap of the CIs, since the estimate close to the edges of the parameter space may be unreliable
  pp1=overlap(pdata[pdata$type==types[1] & pdata$prior==0,param],pdata[pdata$type==types[1] & pdata$prior==1,param],from=min(pdata[pdata$type==types[1],param]),to=max(pdata[pdata$type==types[1],param]),edge.of.parameter.space=TRUE) ##We may want to calculate the overlap of the CIs, since the estimate close to the edges of the parameter space may be unreliable
  pp2=overlap(pdata[pdata$type==types[2] & pdata$prior==0,param],pdata[pdata$type==types[2] & pdata$prior==1,param],from=min(pdata[pdata$type==types[2],param]),to=max(pdata[pdata$type==types[2],param]),edge.of.parameter.space=TRUE) ##We may want to calculate the overlap of the CIs, since the estimate close to the edges of the parameter space may be unreliable
  print(paste(pp1,pp2))
  d1=density(pdata[pdata$type==types[1] & pdata$prior==0,param])
  d2=density(pdata[pdata$type==types[2]& pdata$prior==0 ,param])
  f1=approxfun(d1$x,d1$y)
  f2=approxfun(d2$x,d2$y)
  mu <- ddply(pdata[pdata$prior==0,], "type", summarise, grp.mean=mean(get(param)))
  mu$density=c(f1(mu[1,]$grp.mean),f2(mu[2,]$grp.mean))
  text <- grobTree(textGrob(paste0("P different=",round(p,2)), x=0.7,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=13,fontface="bold")))
  plot=ggplot(data=pdata[pdata$prior==0,],aes_string(fill="type",x=param))+geom_density(alpha=0.3)+geom_vline(data=mu,aes(xintercept=grp.mean,color=type))+scale_x_continuous(name = "Mutation rate") + scale_y_continuous(name = "Density") + scale_fill_discrete(name="Type",labels=c("Crypts","Biopsies")) + scale_color_discrete(name="Type",labels=c("Crypts","Biopsies")) + labs(title=paste0("Patient ",patient))+annotation_custom(text)+geom_text(data=mu,aes(x=grp.mean,y=density,label=round(grp.mean,4)))
  plotwithp=ggplot(data=pdata,aes_string(fill="type",x=param,linetype="prior"))+geom_density(alpha=0.3)+geom_vline(data=mu,aes(xintercept=grp.mean,color=type))+scale_x_continuous(name = "Mutation rate") + scale_y_continuous(name = "Density") + scale_fill_discrete(name="Type",labels=c("Crypts","Biopsies")) + scale_color_discrete(name="Type",labels=c("Crypts","Biopsies")) + labs(title=paste0("Patient ",patient))+annotation_custom(text)+geom_text(data=mu,aes(x=grp.mean,y=density,label=round(grp.mean,4)))+coord_cartesian(xlim=c(0,max(pdata[pdata$prior==0,param])*10))
  save_plot(paste0(patient,".pdf"),plot=plot,base_height = 8)
  save_plot(paste0(patient,".withp.pdf"),plot=plotwithp,base_height = 8)
  results[npatient,]=c(patient,mean(pdata[pdata$type==types[1] & pdata$prior==0 ,param]),mean(pdata[pdata$type==types[2] & pdata$prior==0 ,param]),p)
  #print(paste(patient,p, mean(pdata[pdata$type==types[1],param]),mean(pdata[pdata$type==types[2],param])))
  if(!exists("alldata"))
  {
    alldata=pdata
  }
  else
  {
    alldata=rbind(alldata,pdata)
  }
}
write.csv("results.csv",x = results)