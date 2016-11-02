library(cowplot)
data=read.csv("/Users/Diego/Desktop/singlecrypt/simulation_nonprog/constant_infor.csv")
data=data[,-1]
data2=read.csv("/Users/Diego/Desktop/singlecrypt/simulation_nonprog/parameters.out",sep=" ")
data2_rel=data2
data2_rel$cnv.loss=abs(data2$cnv.loss-1)
data2_rel$cnv.conversion=abs(data2$cnv.conversion-1)
data2_rel$clock.rate=abs(data2$clock.rate-data2$r)/data2$r

rf_gd0=ggplot(data=data[data$gd==0,],aes(y=rf,x=as.factor(rate)))+geom_violin(aes(fill=mod))+scale_x_discrete(name="Mutation rate")+scale_y_continuous(name="RF distance")+labs(title="No whole-genome duplication")+scale_fill_discrete(name="Strategy")
rf_gd=ggplot(data=data[data$gd=="0.001",],aes(y=rf,x=as.factor(rate)))+geom_violin(aes(fill=mod))+scale_x_discrete(name="Mutation rate")+scale_y_continuous(name="RF distance")+labs(title="Simulated whole-genome duplication")+scale_fill_discrete(name="Strategy")
rf_r=ggplot(data=data[data$gd=="random",],aes(y=rf,x=as.factor(rate)))+geom_violin(aes(fill=mod))+scale_x_discrete(name="Mutation rate")+scale_y_continuous(name="RF distance")+labs(title="Random tip-only whole-genome duplication")+scale_fill_discrete(name="Strategy")
grid_rf=plot_grid(rf_gd0,rf_gd,rf_r,labels="AUTO",ncol=3)
rate_gd0=ggplot(data=data2_rel[data2_rel$gd==0,],aes(y=clock.rate,x=as.factor(r)))+geom_violin(aes(fill=mod))+scale_x_discrete(name="Mutation rate")+scale_y_log10(name="Mutation rate relative error")+labs(title="No whole-genome duplication")+scale_fill_discrete(name="Strategy")
rate_gd=ggplot(data=data2_rel[data2_rel$gd=="0.001",],aes(y=clock.rate,x=as.factor(r)))+geom_violin(aes(fill=mod))+scale_x_discrete(name="Mutation rate")+scale_y_log10(name="Mutation rate relative error")+labs(title="Simulated whole-genome duplication")+scale_fill_discrete(name="Strategy")
rate_r=ggplot(data=data2_rel[data2_rel$gd=="random",],aes(y=clock.rate,x=as.factor(r)))+geom_violin(aes(fill=mod))+scale_x_discrete(name="Mutation rate")+scale_y_log10(name="Mutation rate relative error")+labs(title="Random tip-only whole-genome duplication")+scale_fill_discrete(name="Strategy")
grid_rate=plot_grid(rate_gd0,rate_gd,rate_r,labels="AUTO",ncol=3)
supergrid=plot_grid(grid_rf,grid_rate,labels="",ncol=1)
save_plot("results_inforNe.pdf",supergrid,base_height=12,base_aspect_ratio=1.5)
