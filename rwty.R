library(rwty)
library(cowplot)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) !=1) {
        print("Usage: script directory")
        q()
}

setwd(args[1])
my.trees=load.multi(".",ext.tree="trees",ext.p="log")
results=analyze.rwty(my.trees,burnin=500,fill.color='likelihood')
save(results,file="results.Rdata")


for(name in names(results)) { 
	tryCatch({
					save_plot(paste0(name,".pdf"),results[[name]],base_height=12,base_aspect_ratio=1.4)
				},error=function(e){
					save_plot(paste0(name,".pdf"),results[[name]][[1]],base_height=12,base_aspect_ratio=1.4)
			})
	}
