library(rwty)
library(cowplot)
library(methods)

rwty.processors=1
burnin_p=5

args <- commandArgs(trailingOnly = TRUE)

if (length(args) !=3) {
  print("Usage: script directory burnin_proportion n_cores")
  q()
} else {
  print(paste0("Running the script in the folder ",args[1]," with a burnin of ",args[2]," % of the trees and ",args[3],"processors"))
}

setwd(args[1])
burnin_p=as.numeric(args[2])
rwty.processors=as.numeric(args[3])

my.trees=load.multi(path=".",format = "beast")
burnin=round(length(my.trees[[1]]$trees)*burnin_p/100)
results=analyze.rwty(my.trees,burnin=burnin,fill.color='likelihood',params=c("posterior","likelihood","clock.rate","treeModel.rootHeight","luca_height","constant.popSize","cnv.conversion","cnv.loss"),)
results$pseudo.ess=makeplot.pseudo.ess(my.trees,burnin=burnin) #Default number of replicates, 20
save(results,file="results.Rdata")

printplot = function(obj,file="out.pdf") {
  require(cowplot)
  require(methods)
  if(is(obj,"gg")) {
    save_plot(file,obj,base_height=12,base_aspect_ratio=1.4)
    print(file)
  } else if (is(obj,"recordedplot")) {
    print(file)
    pdf(file=file,width=12*1.4,height=12)
    print(obj)
    dev.off()
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
  } else {
    printplot(results[[name]],paste0(name,".pdf"))
  }
}
