library(cowplot)
library(birdring)
setwd("/Users/Diego/Desktop/singlecrypt/simulation/final2")

burnin=1000 #10%

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


mcmcfiles=list.files(path="mcmc/compare",pattern=".*.log$",full.names = TRUE)
mc3files=list.files(path="mc3/compare",pattern=".*.log$",full.names = TRUE)

stopifnot(length(mcmcfiles)==length(mc3files))

results=data.frame(cond=vector(mode = "character",length = length(mcmcfiles)),overlap=vector(mode="numeric",length=length(mcmcfiles)),stringsAsFactors = FALSE)

for (filen in 1:length(mcmcfiles)) {
  name=gsub(pattern = "mcmc/compare/", replacement = "", x = mcmcfiles[filen])
  datamcmc=read.csv(paste0("mcmc/compare/",name),sep="\t",comment.char="#")
  datamc3=read.csv(paste0("mc3/compare/",name),sep="\t",comment.char="#")
  if (nrow(datamc3)!=10000) { ##Three replicates with problems for mc3
    results[filen,]=c(name,NA)
    next
  }
  print(paste0("Working on ",name))
  datamcmc=datamcmc[seq(from=burnin+1,to=nrow(datamcmc),by=1),]
  datamc3=datamc3[seq(from=burnin+1,to=nrow(datamc3),by=1),]
  
  datamcmc$loss.rate=datamcmc$cnv.loss*datamcmc$clock.rate
  datamcmc$conversion.rate=datamcmc$cnv.conversion*datamcmc$clock.rate
  datamcmc$gain.rate=datamcmc$clock.rate
  datamcmc$sga.rate=datamcmc$gain.rate+datamcmc$conversion.rate+datamcmc$loss.rate
  
  datamc3$loss.rate=datamc3$cnv.loss*datamc3$clock.rate
  datamc3$conversion.rate=datamc3$cnv.conversion*datamc3$clock.rate
  datamc3$gain.rate=datamc3$clock.rate
  datamc3$sga.rate=datamc3$gain.rate+datamc3$conversion.rate+datamc3$loss.rate
  p=overlap(datamcmc$sga.rate,datamc3$sga.rate,from=min(c(datamcmc$sga.rate,datamc3$sga.rate)),to=max(c(datamcmc$sga.rate,datamc3$sga.rate)),edge.of.parameter.space=TRUE)
  results[filen,]=c(name,p)
}
results[,2]=as.numeric(results[,2])
hist(results[,2])
dataresults=as.data.frame(results)
plot1=ggplot(data=dataresults,aes(x=overlap))+geom_histogram(fill="#F8766D")+scale_x_continuous(name="Posterior distribution overlap")+scale_y_continuous(name="Number of replicates")+labs(title="Overlap between MCMC and MC3")
mean(results[,2]>0.8,na.rm = TRUE)
mean(results[,2]>0.75,na.rm = TRUE)
save_plot(filename="overlap.pdf",plot=plot1,base_height = 5)
save_plot(filename="overlap.png",plot=plot1,base_height = 5)



