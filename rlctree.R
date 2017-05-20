library(cowplot)
library(ggtree)
library(lubridate)

##Pierre's code to rename tips
shortSampNames <- function(x, tp1, tp2, patient=FALSE, time=TRUE) {
  if (x == "Normal" || is.na(x)) ##DM: added isna
    return(x)
  p <- gsub("^([0-9]+)_.+$", "\\1", x)
  t1 <- ifelse(length(grep(tp1[p], x)) > 0, 1, 0)
  tp <- gsub("^_", "", ifelse(t1, tp1[p], tp2[p]))
  y <- gsub("^[0-9]+_", "", x)
  y <- gsub("_[0-9]+$", "", y)
  y <- gsub(paste("^", tp, sep=""), "", y)
  y <- paste(ifelse(t1, "TP1", "TP2"), "_", y, sep="")
  if (patient) {
    return(paste(p, gsub("__", "_", gsub("_$", "", y)), sep="_"))
  }
  label=gsub("__", "_", gsub("_$", "", y))
  if(!time) {
    label=gsub("TP._","",label)
  }
  return(label)
}

tp1 <- c("_106P_", "_67K|62K_", "_44P_", "_115T_", "_78M_", "_85O|Endoscopy_", "_113M_", "_5U_")
tp2 <- c("_42Y_", "_175R_", "_90W_", "_76Z_", "__", "_Surgical_", "__", "_97A_")
names(tp1) <- names(tp2) <- c("256", "391", "437", "451", "740", "848", "852", "911")
# 
# equiv848 <- c()
# equiv848["848_Endoscopy_SS9_1_22608"] <- "848_Endoscopy_SS9_1_1_22608"
# equiv848["848_Endoscopy_SS9_1_22609"] <- "848_Endoscopy_SS9_1_2_22609"
# equiv848["848_Endoscopy_SS9_4_22617"] <- "848_Endoscopy_SS9_4_1_22617"
# equiv848["848_Endoscopy_SS9_4_22618"] <- "848_Endoscopy_SS9_4_2_22618"

# lngoMat <- read.table(sprintf("%s/%s_phased_100_lngoMat.txt", phasedir, p), header=T, stringsAsFactors=F)
# lngoMat$Normal <- rep("N", nrow(lngoMat))
# colnames(lngoMat) <- gsub("^X", "", colnames(lngoMat))
# colnames(lngoMat) <- gsub("852_113M_SS2_whole_epi__22648", "852_113M_SS8_whole_epi__22648", colnames(lngoMat))
# colnames(lngoMat) <- gsub("_62K_", "_67K_", colnames(lngoMat))
# colnames(lngoMat)[which(colnames(lngoMat) %in% names(equiv848))] =
#   equiv848[colnames(lngoMat)[which(colnames(lngoMat) %in% names(equiv848))]]
# dnaID <- gsub("^.+_([0-9]+$)", "\\1", colnames(lngoMat)[-c(1:5)])
# oldnames <- colnames(lngoMat)[-c(1:5)]
# sampNames <- colnames(lngoMat)[-c(1:5)] <- unlist(lapply(oldnames, shortSampNames, tp1=tp1, tp2=tp2))

setwd("/Users/Diego/Desktop/singlecrypt/results/rlc2/")
##Tree input
tree <- read.beast("/Users/Diego/Desktop/singlecrypt/results/rlc2/MCC/852_phased_100crypts_mod_MCC.trees")

##HARDCODED
mrse=783673200
dobe=-1581033600
daysinayear=364.25
luca.branch=3.455
##################  
  
years=as.numeric(sub(pattern = "d.*", replacement = "",x = seconds_to_period(mrse-dobe)))/daysinayear
remainderdays=(years-floor(years))*daysinayear
months=remainderdays/(daysinayear/12)
days=(months-floor(months))*(daysinayear/12)
mrsd=sprintf("%04d-%02d-%02d",floor(years),floor(months),floor(days))

##HARDCODED, All B8 have genome doublings, and apparentyl only them
treeframe=fortify(tree)
nodegd=MRCA(tree,which(grepl("B8",treeframe$label)))
treegrouped=groupClade(tree,node=nodegd)


treegroupedframe=fortify(treegrouped)
toadd=years-max(treegroupedframe[,"x"])

treegroupedframe$label=unlist(lapply(treegroupedframe$label,shortSampNames,tp1=tp1,tp2=tp2,time=FALSE))

treegroupedframe$x=treegroupedframe$x+toadd
treegroupedframe$branch=treegroupedframe$branch+toadd

root=which(treegroupedframe$parent==treegroupedframe$node)
luca=nrow(treegroupedframe)+1
treegroupedframe[luca,]=treegroupedframe[root,]
treegroupedframe[root,]=treegroupedframe[luca,]
treegroupedframe[root,"parent"]=luca
treegroupedframe[root,"branch.length"]=luca.branch
treegroupedframe[luca,"x"]=treegroupedframe[root,"x"]-luca.branch
treegroupedframe[luca,"node"]=luca

treegroupedframe$posteriorasterisk=as.character(ifelse(treegroupedframe$posterior>=0.8,"*",""))

tree852=ggtree(treegroupedframe,aes(color=log10(rate)),as.Date = FALSE,size=1.5)+
  scale_color_continuous(name=expression("log"[10]*"(SCA rate)"),oob=scales::squish,low="blue", high="red")+
  theme_tree2(legend.position='right')+
  geom_text2(aes(label=posteriorasterisk),vjust=+.7,color="black",size=rel(8.5)) +
  geom_text2(aes(x=branch, subset=(indicator>0.3 & !is.na(indicator) | node == 63),label=sprintf("%.1E (%.1f)",rate,indicator)),vjust=-.5,hjust=+1,size=rel(5)) +
  geom_tiplab(align = TRUE,offset = 3,linesize = 0.5,size=rel(2.5)) +
  ##geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
##  scale_x_continuous(name="Years from DOB") + ##This doesn't work
  #ggplot2::xlim(25,years+15) +
  geom_cladelabel(node=nodegd,label = "GD",offset = 8,align = T,barsize=1.3) +
##  scale_linetype_manual(name="GD",values=c("0"=1,"1"=2),labels=c("No","Yes")) +
  labs(title="852-P") + theme(plot.title= element_text(size = rel(4),hjust = 0.5))
save_plot("852P.pdf",tree852,base_height = 10,base_aspect_ratio = 1.5)

tree <- read.beast("/Users/Diego/Desktop/singlecrypt/results/rlc2/MCC/391_phased_100crypts_MCC.trees")

##HARDCODED
mrse=836204400
dobe=-1258934400
daysinayear=364.25
luca.branch=10.574
##################  

years=as.numeric(sub(pattern = "d.*", replacement = "",x = seconds_to_period(mrse-dobe)))/daysinayear
remainderdays=(years-floor(years))*daysinayear
months=remainderdays/(daysinayear/12)
days=(months-floor(months))*(daysinayear/12)
mrsd=sprintf("%04d-%02d-%02d",floor(years),floor(months),floor(days))

##HARDCODED
dups=c("391_175R_A3_2_1_22315","391_175R_A3_2_2_22316","391_175R_A3_3_1_22318","391_175R_A3_3_2_22319","391_175R_A3_4_1_22321","391_175R_A3_4_2_22322","391_175R_A4_1_1_22276","391_175R_A4_3_1_22282","391_175R_A4_3_2_22283","391_175R_A5_1_2_22229","391_175R_A5_3_1_22234","391_175R_A5_3_2_22235","391_175R_B2_2_1_22303","391_175R_B2_2_2_22304","391_175R_B2_3_3_22308","391_175R_B4_2_1_22255","391_175R_B4_2_2_22256","391_175R_B4_3_1_22258","391_175R_B4_3_2_22259","391_175R_C2_1_2_22289","391_175R_C2_2_1_22291","391_175R_C2_2_2_22292","391_175R_C2_3_2_22295","391_175R_C3_1_2_22265","391_175R_C3_2_1_22267","391_175R_C3_2_2_22268","391_175R_C3_3_1_22270","391_175R_C3_4_2_22274","391_175R_C4_2_1_22243")
treeframe=fortify(tree)
nodegd=MRCA(tree,dups)
treegrouped=groupClade(tree,node=nodegd)

#treegroupedframe=treeframe

treegroupedframe=fortify(treegrouped)
toadd=years-max(treegroupedframe[,"x"])

treegroupedframe$label=unlist(lapply(treegroupedframe$label,shortSampNames,tp1=tp1,tp2=tp2,time=FALSE))

treegroupedframe$x=treegroupedframe$x+toadd
treegroupedframe$branch=treegroupedframe$branch+toadd

root=which(treegroupedframe$parent==treegroupedframe$node)
luca=nrow(treegroupedframe)+1
treegroupedframe[luca,]=treegroupedframe[root,]
treegroupedframe[root,]=treegroupedframe[luca,]
treegroupedframe[root,"parent"]=luca
treegroupedframe[root,"branch.length"]=luca.branch
treegroupedframe[luca,"x"]=treegroupedframe[root,"x"]-luca.branch
treegroupedframe[luca,"node"]=luca

treegroupedframe$printrate=treegroupedframe$node%in%ratenodes
treegroupedframe$posteriorasterisk=as.character(ifelse(treegroupedframe$posterior>=0.8,"*",""))

tree391=ggtree(treegroupedframe,aes(color=log10(rate)),as.Date = FALSE,size=1.5)+
  scale_color_continuous(name=expression("log"[10]*"(SCA rate)"),oob=scales::squish,low="blue", high="red")+
  theme_tree2(legend.position='right')+
  geom_text2(aes(label=posteriorasterisk),vjust=+.7,color="black",size=rel(8.5)) +
  geom_text2(aes(x=branch, subset=(indicator>0.3 & !is.na(indicator) | node == 54),label=sprintf("%.1e (%.1f)",rate,indicator)),vjust=-.5,hjust=+1,size=rel(5)) +
  geom_tiplab(align = TRUE,offset = 3,linesize = 0.5,size=rel(2.5)) +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  ##  scale_x_continuous(name="Years from DOB") + ##This doesn't work
  #ggplot2::xlim(30,NA) +
  geom_cladelabel(node=nodegd,label = "GD",offset = 4.5,align = T,barsize=1.3) +
  scale_linetype_manual(name="GD",values=c("0"=1,"1"=2),labels=c("No","Yes")) +
  labs(title="391-P") + theme(plot.title= element_text(size = rel(4),hjust = 0.5))
save_plot("391P.pdf",tree391,base_height = 10,base_aspect_ratio = 1.5)
library(cowplot)
grid=plot_grid(tree852,tree391,labels = "auto")
save_plot("figure.pdf",grid,base_height = 6,base_aspect_ratio = 3)

#########BIOPSY DATA#############

#Tree input
tree <- read.beast("/Users/Diego/Desktop/singlecrypt/results/rlc2/MCC/852_phased_100pool_MCC.trees")

##HARDCODED
mrse=783673200
dobe=-1581033600
luca.branch=7.516
##################  

years=as.numeric(sub(pattern = "d.*", replacement = "",x = seconds_to_period(mrse-dobe)))/daysinayear
remainderdays=(years-floor(years))*daysinayear
months=remainderdays/(daysinayear/12)
days=(months-floor(months))*(daysinayear/12)
mrsd=sprintf("%04d-%02d-%02d",floor(years),floor(months),floor(days))

treeframe=fortify(tree)
treegrouped=groupClade(tree,node=nodegd)
treegroupedframe=treeframe
toadd=years-max(treegroupedframe[,"x"])

treegroupedframe$label=unlist(lapply(treegroupedframe$label,shortSampNames,tp1=tp1,tp2=tp2,time=FALSE))

treegroupedframe$x=treegroupedframe$x+toadd
treegroupedframe$branch=treegroupedframe$branch+toadd

root=which(treegroupedframe$parent==treegroupedframe$node)
luca=nrow(treegroupedframe)+1
treegroupedframe[luca,]=treegroupedframe[root,]
treegroupedframe[root,]=treegroupedframe[luca,]
treegroupedframe[root,"parent"]=luca
treegroupedframe[root,"branch.length"]=luca.branch
treegroupedframe[luca,"x"]=treegroupedframe[root,"x"]-luca.branch
treegroupedframe[luca,"node"]=luca

treegroupedframe$posteriorasterisk=as.character(ifelse(treegroupedframe$posterior>=0.8,"*",""))

#which(treegroupedframe$parent==treegroupedframe[root,]$node) To get the additional rate
addnode=13

tree852pool=ggtree(treegroupedframe,aes(color=log10(rate)),as.Date = FALSE,size=1.5)+
  scale_color_continuous(name=expression("log"[10]*"(SCA rate)"),oob=scales::squish,low="blue", high="red")+
  theme_tree2(legend.position='right')+
  geom_text2(aes(label=posteriorasterisk),vjust=+.7,color="black",size=rel(8.5)) +
  geom_text2(aes(x=branch, subset=(indicator>0.3 & !is.na(indicator) | node == addnode),label=sprintf("%.1E (%.1f)",rate,indicator)),vjust=-.5,hjust=+1,size=rel(5)) +
  geom_tiplab(align = TRUE,offset = 3,linesize = 0.5,size=rel(2.5)) +
  ##geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  ##  scale_x_continuous(name="Years from DOB") + ##This doesn't work
  ggplot2::xlim(45,years+5) +
  #geom_cladelabel(node=nodegd,label = "GD",offset = 8,align = T,barsize=1.3) +
  ##  scale_linetype_manual(name="GD",values=c("0"=1,"1"=2),labels=c("No","Yes")) +
  labs(title="852-P") + theme(plot.title= element_text(size = rel(4),hjust = 0.5))
save_plot("852pool.pdf",tree852pool,base_height = 10,base_aspect_ratio = 1.5)


tree <- read.beast("/Users/Diego/Desktop/singlecrypt/results/rlc2/MCC/391_phased_100pool_MCC.trees")

##HARDCODED
mrse=836204400
dobe=-1258934400
daysinayear=364.25
luca.branch=11.919
##################  

years=as.numeric(sub(pattern = "d.*", replacement = "",x = seconds_to_period(mrse-dobe)))/daysinayear
remainderdays=(years-floor(years))*daysinayear
months=remainderdays/(daysinayear/12)
days=(months-floor(months))*(daysinayear/12)
mrsd=sprintf("%04d-%02d-%02d",floor(years),floor(months),floor(days))

##GD non-monophyletic, I will add it manually
treeframe=fortify(tree)
#nodegd=MRCA(tree,dups)
#treegrouped=groupClade(tree,node=nodegd)

treegroupedframe=treeframe

#treegroupedframe=fortify(treegrouped)
toadd=years-max(treegroupedframe[,"x"])

treegroupedframe$label=unlist(lapply(treegroupedframe$label,shortSampNames,tp1=tp1,tp2=tp2,time=FALSE))

treegroupedframe$x=treegroupedframe$x+toadd
treegroupedframe$branch=treegroupedframe$branch+toadd

root=which(treegroupedframe$parent==treegroupedframe$node)
luca=nrow(treegroupedframe)+1
treegroupedframe[luca,]=treegroupedframe[root,]
treegroupedframe[root,]=treegroupedframe[luca,]
treegroupedframe[root,"parent"]=luca
treegroupedframe[root,"branch.length"]=luca.branch
treegroupedframe[luca,"x"]=treegroupedframe[root,"x"]-luca.branch
treegroupedframe[luca,"node"]=luca

treegroupedframe$posteriorasterisk=as.character(ifelse(treegroupedframe$posterior>=0.8,"*",""))

#which(treegroupedframe$parent==treegroupedframe[root,]$node) To get the additional rate
addnode=12

#ggtree(treegroupedframe,aes(color=log2(rate),,linetype=group),as.Date = FALSE,mrsd=mrsd,size=1.5)+
tree391pool=ggtree(treegroupedframe,aes(color=log10(rate)),as.Date = FALSE,size=1.5)+
  scale_color_continuous(name=expression("log"[10]*"(SCA rate)"),oob=scales::squish,low="blue", high="red")+
  theme_tree2(legend.position='right')+
  #geom_text2(aes(label=round(posterior,2)),hjust=-.3,color="black")+
  geom_text2(aes(label=posteriorasterisk),vjust=+.7,color="black",size=rel(8.5)) +
  geom_text2(aes(x=branch, subset=(indicator>0.3 & !is.na(indicator) | node == addnode),label=sprintf("%.1e (%.1f)",rate,indicator)),vjust=-.5,hjust=+1,size=rel(5)) +
  geom_tiplab(align = TRUE,offset = 3,linesize = 0.5,size=rel(2.5)) +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  ##  scale_x_continuous(name="Years from DOB") + ##This doesn't work
  ggplot2::xlim(40,75) +
  #geom_cladelabel(node=nodegd,label = "GD",offset = 4.5,align = T,barsize=1.3) +
  #scale_linetype_manual(name="GD",values=c("0"=1,"1"=2),labels=c("No","Yes")) +
  labs(title="391-P") + theme(plot.title= element_text(size = rel(4),hjust = 0.5))
save_plot("391P.pdf",tree391pool,base_height = 10,base_aspect_ratio = 1.5)
library(cowplot)
grid_pool=plot_grid(tree852pool,tree391pool,labels = "auto")
save_plot("figure_pool.pdf",grid_pool,base_height = 6,base_aspect_ratio = 3)



