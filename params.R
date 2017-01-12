library(cowplot)
library(plyr)
library(birdring)
library(grid)
library(ggrepel)
setwd("/Users/Diego/Desktop/singlecrypt/results/500_prob/results")
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
    } else
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
  d1p=density(pdata[pdata$type==types[1] & pdata$prior==1,param])
  d2p=density(pdata[pdata$type==types[2]& pdata$prior==1 ,param])
  f1=approxfun(d1$x,d1$y)
  f2=approxfun(d2$x,d2$y)
  f1p=approxfun(d1p$x,d1p$y)
  f2p=approxfun(d2p$x,d2p$y)
  mu <- ddply(pdata[pdata$prior==0,], "type", summarise, grp.mean=mean(get(param)))
  mu$density=c(f1(mu[1,]$grp.mean),f2(mu[2,]$grp.mean))
  mu$prior=c(0,0)
  mup= ddply(pdata[pdata$prior==1,], "type", summarise, grp.mean=mean(get(param)))
  mup$density=c(f1p(mup[1,]$grp.mean),f2p(mup[2,]$grp.mean))
  mup$prior=c(1,1)
  mu=rbind(mu,mup)
  text <- grobTree(textGrob(paste0("P different=",round(p,2)), x=0.7,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=13,fontface="bold")))
  plot=ggplot(data=pdata[pdata$prior==0,],aes_string(fill="type",x=param))+geom_density(alpha=0.3)+geom_vline(data=mu,aes(xintercept=grp.mean,color=type))+scale_x_continuous(name = "Mutation rate") + scale_y_continuous(name = "Density") + scale_fill_discrete(name="Type",labels=c("Crypts","Biopsies")) + scale_color_discrete(name="Type",labels=c("Crypts","Biopsies")) + labs(title=paste0("Patient ",patient))+annotation_custom(text)+geom_text_repel(data=mu,aes(x=grp.mean,y=density,label=round(grp.mean,4)),force=2)
  #plotwithp=ggplot(data=pdata,aes_string(fill="type",x=param,linetype="prior"))+geom_density(alpha=0.3)+geom_vline(data=mu,aes(xintercept=grp.mean,color=type))+scale_x_continuous(name = "Mutation rate") + scale_y_continuous(name = "Density") + scale_fill_discrete(name="Type",labels=c("Crypts","Biopsies")) + scale_color_discrete(name="Type",labels=c("Crypts","Biopsies")) + labs(title=paste0("Patient ",patient))+annotation_custom(text)+geom_text(data=mu,aes(x=grp.mean,y=density,label=round(grp.mean,4)))+coord_cartesian(xlim=c(0,max(pdata[pdata$prior==0,param])*10))
  #save_plot(paste0(patient,".withp.pdf"),plot=plotwithp,base_height = 8)
  #The prior is mean="-4.0" stdev="2.5"
  save_plot(paste0(patient,".pdf"),plot=plot,base_height = 8)
  results[npatient,]=c(patient,mean(pdata[pdata$type==types[1] & pdata$prior==0 ,param]),mean(pdata[pdata$type==types[2] & pdata$prior==0 ,param]),p)
  #print(paste(patient,p, mean(pdata[pdata$type==types[1],param]),mean(pdata[pdata$type==types[2],param])))
  if(exists("alldata") == FALSE) {
    alldata=pdata
  } else {
    alldata=rbind(alldata,pdata)
  }
}
write.csv("results.csv",x = results)

# HPDindexes=function(sample,proportion=0.95) { ##Indexes of the highest posterior density interval (credible interval) ##Untested
#   sorted_indexes=sort(sample,index.return=TRUE)$ix
#   nsamples=round(length(sample)*proportion)
#   finalLindex=1
#   range=Inf
#   for (i in 1:(length(sample)-nsamples)) {
#     nrange=abs(sample[sorted_indexes[i]]-sample[sorted_indexes[i+nsamples]])
#     if (nrange< range) {
#       range=nrange
#       finalLindex=i
#     }
#   }
#   return (sorted_indexes[finalLindex:(finalLindex+nsamples-1)])
# }
# params=c("type","patient")
# listparamvalues=sapply(params,function(x){unique(alldata[,x])})
# 
# makecomb=function(datalist,id,inlist){
#   thisoutlist=list()
#   for (i in 1:length(datalist[[id]])) {
#     iinlist=c(inlist,datalist[[id]][[i]])
#     if(id<length(datalist)) {
#       thisoutlist=c(thisoutlist,makecomb(datalist,id+1,iinlist))
#     } else {
#       thisoutlist=c(thisoutlist,list(iinlist))
#     }
#   }
#   return(thisoutlist)
# }
# listcombs=makecomb(listparamvalues,1,NULL)

HPDmask=function(mydata,param,mask,proportion=0.95) { ##Indexes of the highest posterior density interval (credible interval)
  sorted_indexes=sort(mydata[,param],index.return=TRUE)$ix
  valid_indexes=sorted_indexes[mask[sorted_indexes]]
  nsamples=round(length(valid_indexes)*proportion)
  finalLindex=1
  range=Inf
  for (i in 1:(length(valid_indexes)-nsamples)) {
    nrange=abs(mydata[valid_indexes[i],param]-mydata[valid_indexes[i+nsamples],param])
    if (nrange < range) {
      range=nrange
      finalLindex=i
    }
  }
  outmask=rep(FALSE,nrow(mydata))
  outmask[valid_indexes[finalLindex:(finalLindex+nsamples-1)]]=TRUE
  return (outmask)
}

params=c(which(names(alldata)=="type"),which(names(alldata)=="patient"))
listparamvalues=sapply(params,function(x){unique(alldata[,x])})

makemask=function(mydata,paramlist,id,mask){
  thisoutlist=list()
  for (i in 1:length(paramlist[[id]])) {
    imask=mydata[,params[id]]==paramlist[[id]][[i]]
    omask=imask & mask
    if(id<length(paramlist)) {
      thisoutlist=c(thisoutlist,makemask(mydata,paramlist,id+1,omask))
    } else {
      thisoutlist=c(thisoutlist,list(omask))
    }
  }
  return(thisoutlist)
}

lmasks=makemask(alldata[alldata$prior==0,],listparamvalues,1,rep(TRUE,nrow(alldata[alldata$prior==0,])))

initp=0.95
finalmask=rep(FALSE,nrow(alldata[alldata$prior==0,]))
for (i in 1:length(lmasks))
{
  finalmask=finalmask|HPDmask(alldata[alldata$prior==0,],param,lmasks[[i]],initp) ##CI0.95 (BASE)
}

nlayers=4

myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type"))+geom_violin(alpha=1/nlayers,size=1.25)

for (j in 1:(nlayers-1)) {
  finalmask=rep(FALSE,nrow(alldata[alldata$prior==0,]))
  for (i in 1:length(lmasks))
  {
    finalmask=finalmask|HPDmask(alldata[alldata$prior==0,],param,lmasks[[i]],(initp)/nlayers*(nlayers-j))
  }
  myplot=myplot+geom_violin(data=(alldata[alldata$prior==0,])[finalmask,],alpha=1/nlayers*(j+1),size=0.5)
  
}