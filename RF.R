library(phangorn)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

if (length(args) !=3) {
        print("Usage: script problem_dir true_trees outfile")
        q()
}
pdir=args[1]
tTreeFile=args[2]
outfile=args[3]
ttrees=read.tree(tTreeFile)
trees=list.files(path=pdir,pattern="*MCC.trees")
final_data=data.frame(rep=numeric(0),rate=numeric(0),gd=character(),mod=character(),rf=numeric(0),stringsAsFactors=FALSE)
for (i in 1:length(trees)){
    ptree=read.nexus(paste0(pdir,"/",trees[i]))
    ptree$tip.label=sub(x=ptree$tip.label,pattern="\\*",replacement="") ##Otherwise the tips won't be the same
    #tree8_r0.001_gdrandom_modbaseline_MCC.trees #Example
    replicate=as.numeric(sub(x=trees[i],pattern="tree(.*)\\_r.*\\_gd.*\\_mod.*\\_MCC.trees","\\1"))
    rate=as.numeric(sub(x=trees[i],pattern="tree.*\\_r(.*)\\_gd.*\\_mod.*\\_MCC.trees","\\1"))
    gd=sub(x=trees[i],pattern="tree.*\\_r.*\\_gd(.*)\\_mod.*\\_MCC.trees","\\1")
    mod=sub(x=trees[i],pattern="tree.*\\_r.*\\_gd.*\\_mod(.*)\\_MCC.trees","\\1")
    ttree=drop.tip(ttrees[[replicate+1]],"outgroup")
    dist=RF.dist(ttree,ptree,check.labels=TRUE,rooted=TRUE,normalize=TRUE)
    #dist=RF.dist(ttree, ptree, check.labels = TRUE,rooted=TRUE)/((ttree$Nnode+1)*2-6) ##s_tree$Nnode= Number of internal nodes. This tree is rooted, so internal nodes+1 = n_leaves. 2*(n-3) = number of internal branches/bipartitions in an unrooted tree * 2.
    final_data[nrow(final_data)+1,]=c(replicate,rate,gd,mod,dist)
}

write.csv(final_data,outfile)
