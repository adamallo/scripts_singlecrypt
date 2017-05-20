library(cowplot)
library(plyr)
library(birdring)
library(grid)
library(ggrepel)
setwd("/Users/Diego/Desktop/singlecrypt/results/newfinal/results/")
arrayPat=c("391","437","451","740","848","852")
allPat=c("256","391","437","451","740","848","852","911")
rainbow <- c("#CC0033", "#FF6633", "#FFCC33", "#99CC33", "#009933", "#009999", "#003399", "#330066")
rainbow10 <- c(rainbow, "gray", "black")
nprog <- c("256", "437", "451", "911")
nprognames=c("256-NP", "437-NP", "451-NP", "911-NP")
prog <- c("391", "740", "848", "852")
prognames = c("391-P","740-P","848-P","852-P")
pats <- c(nprog, prog)
patCol <- rainbow10[1:length(pats)]
typeCol=c(rainbow[3],rainbow[5])
names(typeCol)=c("crypts","pool")
rlc=FALSE
finalnames=character(0)
for (pat in pats) {
  if (pat %in% prog) {
    pname=paste0(pat,"-P")
  } else {
    pname=paste0(pat,"-NP")
  }
  finalnames=c(finalnames,pname)
}
names(patCol) <- finalnames
names(finalnames)=pats

burnin=1000 #10%
#arrayPat=allPat
types=c("crypts","pool")

ages=read.csv(file = "/Users/Diego/Desktop/singlecrypt/newfinal_data/data/ages.tsv",sep = " ")
ages$newpatient=finalnames[as.character(ages$patient)]


####My functions#####
#####################

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

###################

#MAIN
#####
if (length(types) != 2)
{
  print("The number of types to compare must be two")
  exit()
}

rm(alldata)
plotparams=c("gain.rate","conversion.rate","loss.rate","sga.rate","ne","ne2","luca_age") #HARDCODED
n=length(arrayPat)*length(plotparams)
results=data.frame(patient=vector("character",n),crypts_mean=vector("numeric",n),biopsies_mean=vector("numeric",n),p_diff=vector("numeric",n),ratekind=vector("character",n),stringsAsFactors = FALSE)
for (npatient in 1:(length(arrayPat)))
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
    priordata=priordata[,!(names(priordata) %in% c("treeLikelihood"))] ##Priordata had an extra reapeated column
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
  #convenient changes and calculation of absolute rates
  if (patient %in% prog) {
    pname=paste0(patient,"-P")
  } else {
    pname=paste0(patient,"-NP")
  }
  pdata$patient=rep(pname,times=nrow(pdata))
  if (rlc) {
    pdata$cnv.gain=1/(pdata$cnv.loss+pdata$cnv.conversion+1)
    pdata$loss.rate=pdata$cnv.loss*clock.rate
    pdata$conversion.rate=pdata$cnv.conversion*clock.rate
    pdata$gain.rate=pdata$cnv.gain*clock.rate
    pdata$sga.rate=clock.rate
  } else {
    pdata$loss.rate=pdata$cnv.loss*pdata$clock.rate
    pdata$conversion.rate=pdata$cnv.conversion*pdata$clock.rate
    pdata$gain.rate=pdata$clock.rate
    pdata$sga.rate=pdata$gain.rate+pdata$conversion.rate+pdata$loss.rate
  }

  pdata$ne=pdata$constant.popSize*7.3 #0.02gen/day => 50days/gen => 7.3gen/year
  pdata$ne2=pdata$constant.popSize*52.14 #7 days/gen
  pdata$luca_age=apply(pdata,1,function(x){ages[ages$newpatient==pname & as.character(ages$type)==as.character(x["type"]),]$age-as.numeric(x["luca_height"])})
  
  
  for (nparam in 1:length(plotparams)) {
    param=plotparams[nparam]
    name=sub(x=sub(x=param,pattern=".rate",replacement = "",fixed = TRUE),pattern="^(.)",replacement="\\U\\1",perl = TRUE)
    #Estimation of posterior probabilities of the 2 distributions being the same (overlap)
    ##p=Posterior probability of null hypothesis (same distribution)
    p=overlap(pdata[pdata$type==types[1] & pdata$prior==0,param],pdata[pdata$type==types[2] & pdata$prior==0,param],from=min(pdata[pdata$prior==0,param]),to=max(pdata[pdata$prior==0,param]),edge.of.parameter.space=TRUE) ##We may want to calculate the overlap of the CIs, since the estimate close to the edges of the parameter space may be unreliable
    #pp1=Posterior probability of prior=posterior type1
    #pp1=overlap(pdata[pdata$type==types[1] & pdata$prior==0,param],pdata[pdata$type==types[1] & pdata$prior==1,param],from=min(pdata[pdata$type==types[1],param]),to=max(pdata[pdata$type==types[1],param]),edge.of.parameter.space=TRUE) ##We may want to calculate the overlap of the CIs, since the estimate close to the edges of the parameter space may be unreliable
    #pp2=Posterior probability of prior=posterior type2
    #pp2=overlap(pdata[pdata$type==types[2] & pdata$prior==0,param],pdata[pdata$type==types[2] & pdata$prior==1,param],from=min(pdata[pdata$type==types[2],param]),to=max(pdata[pdata$type==types[2],param]),edge.of.parameter.space=TRUE) ##We may want to calculate the overlap of the CIs, since the estimate close to the edges of the parameter space may be unreliable
    #print(paste(pp1,pp2))
    
    #Estimation of the means for each category (Expected value)
    #We estimate also the density of the means in order to label
    #The plot (y axis)
    
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
    
    text <- grobTree(textGrob(paste0("P null hypothesis=",round(p,2)), x=0.7,  y=0.95, hjust=0,gp=gpar(col="black", fontsize=13,fontface="bold")))
    
    plot=ggplot(data=pdata[pdata$prior==0,],aes_string(fill="type",x=param))+geom_density(alpha=0.3)+geom_vline(data=mu[mu$prior==0,],aes(xintercept=grp.mean,color=type))+scale_x_continuous(name = paste0(name," rate")) + scale_y_continuous(name = "Density") + scale_fill_discrete(name="Type",labels=c("Crypts","Biopsies")) + scale_color_discrete(name="Type",labels=c("Crypts","Biopsies")) + labs(title=paste0("Patient ",pname))+annotation_custom(text)+geom_text_repel(data=mu[mu$prior==0,],aes(x=grp.mean,y=density,label=round(grp.mean,4)),force=2)
    #plotwithp=ggplot(data=pdata,aes_string(fill="type",x=param,linetype="prior"))+geom_density(alpha=0.3)+geom_vline(data=mu,aes(xintercept=grp.mean,color=type))+scale_x_continuous(name = "Mutation rate") + scale_y_continuous(name = "Density") + scale_fill_discrete(name="Type",labels=c("Crypts","Biopsies")) + scale_color_discrete(name="Type",labels=c("Crypts","Biopsies")) + labs(title=paste0("Patient ",patient))+annotation_custom(text)+geom_text(data=mu,aes(x=grp.mean,y=density,label=round(grp.mean,4)))+coord_cartesian(xlim=c(0,max(pdata[pdata$prior==0,param])*10))
    #save_plot(paste0(patient,".withp.pdf"),plot=plotwithp,base_height = 8)
    #The prior is mean="-4.0" stdev="2.5"
    save_plot(paste0(pname,param,".pdf"),plot=plot,base_height = 8)
    results[npatient*length(plotparams)-length(plotparams)+nparam,]=c(pname,mean(pdata[pdata$type==types[1] & pdata$prior==0 ,param]),mean(pdata[pdata$type==types[2] & pdata$prior==0 ,param]),p,param)
    #print(paste(patient,p, mean(pdata[pdata$type==types[1],param]),mean(pdata[pdata$type==types[2],param])))
  }
  
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

##make a list of masks, one for each parameter combination

alldata$patient=factor(alldata$patient,levels=c("256-NP","437-NP","451-NP","911-NP","391-P","740-P","848-P","852-P")) ##HARDCODED!!!!!!
alldata$progressor=alldata$patient %in% prognames

params=c(which(names(alldata)=="type"),which(names(alldata)=="patient"))
listparamvalues=sapply(params,function(x){unique(alldata[,x])})
lmasks=makemask(alldata[alldata$prior==0,],listparamvalues,1,rep(TRUE,nrow(alldata[alldata$prior==0,])))
resultsplot=results
resultsplot$p_diff=round(as.numeric(resultsplot$p_diff),3)
resultsplot$y=rep(0,nrow(resultsplot))
resultsplot$patient=factor(resultsplot$patient,levels=c("256-NP","437-NP","451-NP","911-NP","391-P","740-P","848-P","852-P")) ##HARDCODED!!!!!!
resultsplot$progressor=resultsplot$patient %in% prognames
resultsplot$diffrwhole=as.numeric(resultsplot$crypts_mean)/as.numeric(resultsplot$biopsies_mean)
#resultsplot$y=mapply(function(x,y){return(mean(c(x,y)))},as.numeric(resultsplot$crypts_mean),as.numeric(resultsplot$biopsies_mean))*2

for (nparam in 1:length(plotparams)) {
  param=plotparams[nparam]
  name=sub(x=sub(x=param,pattern=".rate",replacement = "",fixed = TRUE),pattern="^(.)",replacement="\\U\\1",perl = TRUE)
  initp=0.95
  finalmask=rep(FALSE,nrow(alldata[alldata$prior==0,]))
  for (i in 1:length(lmasks))
  {
    finalmask=finalmask|HPDmask(alldata[alldata$prior==0,],param,lmasks[[i]],initp) ##CI0.95 (BASE)
  }
  
  tempy=ddply((alldata[alldata$prior==0,])[finalmask,], .(patient), summarise, y=max(get(param)))
  resultsplot[results$ratekind==param,]$y=tempy[match(resultsplot[results$ratekind==param,]$patient,tempy$patient),2]
  nlayers=4

  if (name=="Ne"){
    #myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type",color="patient"))+geom_violin(alpha=1/nlayers,size=1)+scale_y_continuous(name="Effective population size")+scale_x_discrete(name="Patient")+scale_fill_discrete(name="Data type")+scale_color_manual(values=patCol,guide=FALSE)+scale_fill_manual(name="",values=typeCol,labels=c("Crypts","Whole epithelium"))+theme(legend.justification = c(0, 1), legend.position = c(0, 1))
    myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type"))+geom_violin(alpha=1/nlayers,size=1)+scale_y_continuous(name="Effective population size, generation time 50 days")+scale_x_discrete(name="Patient")+scale_fill_discrete(name="Data type")+scale_fill_manual(name="",values=typeCol,labels=c("Crypts","Whole epithelium"))+theme(legend.justification = c(0, 1), legend.position = c(0, 1))
    
  } else if(name=="Ne2"){
    myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type"))+geom_violin(alpha=1/nlayers,size=1)+scale_y_continuous(name="Effective population size, generation time 7 days")+scale_x_discrete(name="Patient")+scale_fill_discrete(name="Data type")+scale_fill_manual(name="",values=typeCol,labels=c("Crypts","Whole epithelium"))+theme(legend.justification = c(0, 1), legend.position = c(0, 1))
    
  } else if(name=="Luca_age"){
    myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type"))+geom_violin(alpha=1/nlayers,size=1)+scale_y_continuous(name="Last universal common ancestor age (years from birthday)")+scale_x_discrete(name="Patient")+scale_fill_discrete(name="Data type")+scale_fill_manual(name="",values=typeCol,labels=c("Crypts","Whole epithelium"))+theme(legend.justification = c(0, 1), legend.position = c(0, 1))
  } else {
    #myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type",color="patient"))+geom_violin(alpha=1/nlayers,size=1)+scale_y_continuous(name=paste0(name," rate (events/year/fragment/allele)"))+scale_x_discrete(name="Patient")+scale_fill_discrete(name="Data type")+scale_color_manual(values=patCol,guide=FALSE)+scale_fill_manual(name="",values=typeCol,labels=c("Crypts","Whole epithelium"))+theme(legend.justification = c(0, 1), legend.position = c(0, 1))
    myplot=ggplot(data=(alldata[alldata$prior==0,])[finalmask,],aes_string(y=param,x="patient",fill="type"))+geom_violin(alpha=1/nlayers,size=1)+scale_y_continuous(name=paste0(name," rate (events/year/fragment/allele)"))+scale_x_discrete(name="Patient")+scale_fill_discrete(name="Data type")+scale_fill_manual(name="",values=typeCol,labels=c("Crypts","Whole epithelium"))+theme(legend.justification = c(0, 1), legend.position = c(0, 1))
  }
  for (j in 1:(nlayers-1)) {
    finalmask=rep(FALSE,nrow(alldata[alldata$prior==0,]))
    for (i in 1:length(lmasks))
    {
      finalmask=finalmask|HPDmask(alldata[alldata$prior==0,],param,lmasks[[i]],(initp)/nlayers*(nlayers-j))
    }
    myplot=myplot+geom_violin(data=(alldata[alldata$prior==0,])[finalmask,],alpha=1/nlayers*(j+1),size=0.3,color="black",linetype=3)
  }
  myplot=myplot+stat_summary(fun.y=mean,geom="point",position = position_dodge(width = 0.9),aes(shape=progressor),color="black",size=3)+geom_text(data=resultsplot[resultsplot$ratekind==param,],aes(y=y,x=patient,label=p_diff),color="black")+scale_shape_discrete(name="",labels=c("Non progressor","Progressor"))
  save_plot(paste0(name,".pdf"),plot=myplot,base_height=12)
}

##Statistical tests
resultsplot[resultsplot$ratekind=="sga.rate",]
wilcox.test(as.numeric(resultsplot[resultsplot$ratekind=="sga.rate" & resultsplot$progressor==TRUE,]$crypts_mean),as.numeric(resultsplot[resultsplot$ratekind=="sga.rate" & resultsplot$progressor==FALSE,]$crypts_mean))
resultsplot[resultsplot$ratekind=="luca_age",]
wilcox.test(as.numeric(resultsplot[resultsplot$ratekind=="luca_age" & resultsplot$progressor==TRUE,]$crypts_mean),as.numeric(resultsplot[resultsplot$ratekind=="luca_age" & resultsplot$progressor==FALSE,]$crypts_mean))





