library(rwty)
library(RcmdrPlugin.KMggplot2)
library(cowplot)
library(methods)

rwty.processors=1

burnin_p=5

args <- commandArgs(trailingOnly = TRUE)

if (length(args) !=1) {
  print("Usage: script directory")
  q()
}

setwd(args[1])
my.trees=load.multi(path=".",format = "beast")
burnin=round(length(my.trees[[1]]$trees)*burnin_p/100)
results=analyze.rwty(my.trees,burnin=burnin,fill.color='likelihood')
results$pseudo.ess(my.trees,burnin=burnin) #Default number of replicates, 20
results$pairsb=makeplot.pairs(my.trees,burnin=burnin,params=c("posterior","likelihood","clock.rate","treeModel.rootHeight","luca_height","constant.popSize","cnv.conversion","cnv.loss"))
save(results,file="results.Rdata")

printplot = function(obj,file="out.pdf") {
  require(cowplot)
  require(RcmdrPlugin.KMggplot2)
  if(is(obj,"gg")) {
    save_plot(file,obj,base_height=12,base_aspect_ratio=1.4)
  } else if (is(obj,"recordedplot")) {
    ggsaveKmg2(filename = file, plot = obj,height=12,width = 12*1.4)
  } else if (is(obj,"citation")) {
    #We don't need to print this
  } else {
    stop(paste0("Unsuported class", class(obj)))
  }
}

for(name in names(results)) {
  if (is(results[[name]],"list")) {
    for (subname in names(results[[name]]))
    {
      printplot(results[[name]][[subname]],paste0(paste(name,subname,sep="_"),".pdf"))
    }
  }
  else {
    printplot(results[[name]],paste0(name,".pdf"))
  }
}