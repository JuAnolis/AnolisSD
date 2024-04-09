#SURFACE#

####Ejecutar esta parte 

data_females<-read.csv("female.csv",header = T)
data_females<-textshape::column_to_rownames(data_females, loc = 1)
f_SVL = as.numeric((data_females)[,1])
f_HlL =as.numeric((data_females)[,2])
f_HL = as.numeric((data_females)[,3])
f_LN = as.numeric((data_females)[,4])
data_males<-read.csv("male.csv",header = T)
data_males<-textshape::column_to_rownames(data_males, loc = 1)
m_SVL = as.numeric((data_males)[,1])
m_HlL =as.numeric((data_males)[,2])
m_HL = as.numeric((data_males)[,3])
m_LN = as.numeric((data_males)[,4])
log_f_SVL<-log(f_SVL)
log_f_HlL<-log(f_HlL)
log_f_HL<-log(f_HL)
log_f_LN<-log(f_LN)
log_m_SVL<-log(m_SVL)
log_m_HlL<-log(m_HlL)
log_m_HL<-log(m_HL)
log_m_LN<-log(m_LN)
f_names = row.names(data_females)
m_names = row.names(data_males)
f_regHlL = lm(log_f_SVL ~ log_f_HlL) 
f_regHL = lm(log_f_SVL ~ log_f_HL) 
f_regLN = lm(log_f_SVL ~ log_f_LN) 
f_HlL.res = f_regHlL$resid
f_HL.res = f_regHL$resid
f_LN.res = f_regLN$resid
f_data = data.frame(cbind(log_f_SVL, f_HlL.res, f_HL.res, f_LN.res),row.names = f_names)
m_regHlL = lm(log_m_SVL ~ log_m_HlL) 
m_regHL = lm(log_m_SVL ~ log_m_HL) 
m_regLN = lm(log_m_SVL ~ log_m_LN) 
m_HlL.res = m_regHlL$resid
m_HL.res = m_regHL$resid
m_LN.res = m_regLN$resid
m_names = row.names(data_males)
m_data = data.frame(cbind(log_m_SVL, m_HlL.res, m_HL.res, m_LN.res),row.names = m_names)
tree = read.nexus("dimorfismo.tree")
plot(tree)
dimorfismo.tree = nameNodes(tree)
write.csv(f_data,"f_data.csv")
write.csv(m_data,"m_data.csv")

##males
chk_f=name.check(dimorfismo.tree, f_data)
chk_f
olist = convertTreeData(dimorfismo.tree,m_data)
otree = olist[[1]]; odata = olist[[2]]
fwd_m = surfaceForward(otree_m, odata, aic_threshold = 0, exclude = 0,
                      verbose = FALSE, plotaic = FALSE)
k_m = length(fwd_m)
fsum_m = surfaceSummary(fwd_m)
names(fsum_m)
fsum_m$aics
fsum_m$n_regimes
fsum_m
bwd_m = surfaceBackward(otree, odata, starting_model = fwd_m[[k_m]], aic_threshold = 0,
                       only_best = TRUE, verbose = FALSE, plotaic = FALSE)
bsum_m = surfaceSummary(bwd_m)
kk_m = length(bwd_m)
print("Mis valores alfa para cada rasgo son:")
data.frame(bsum_m$alpha)
print("Mis valores sigma para cada rasgo son:")
data.frame(bsum_m$sigma_squared)
print("Mis valores de aic para cada iteración son:")
data.frame(bsum_m$aics)
print("Mis valores de theta son:")
bsum_m$theta
bsum_m$n_regimes
bsum_m$aics
tiff('surf_m_fullHD.tiff', pointsize=10, width=2800, height=3500, res=600)
par(lwd = 5)
surf_m<-surfaceTreePlot(dimorfismo.tree, bwd_m[[kk_m]], cols= NULL,  convcol=T, labelshifts = T)# grafica el árbol con los picos adaptativos generados por el SURFACE
#surf_m<-recordPlot()
#surf_m<-as.ggplot(surf_m)
dev.off()
tiff("surf_m.tiff",width= 800,heigh=400,res = 800)
plot(surf_m)

dpi=600
surfaceAICPlot(fwd_m, bwd_m)
abline(h=bm[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic)-2,c("BM","OU1"),cex=0.5)

fitBM_m<-fitContinuous(tree,m_data)
fitBM_m
fitOU_m<-fitContinuous(tree,m_data,model="OU")
fitOU_m
fitEB_m<-fitContinuous(tree,m_data,model="EB")
fitEB_m

#females
olist_f = convertTreeData(dimorfismo.tree,f_data)
otree_f = olist_f[[1]]; odata_f = olist_f[[2]]
fwd_f = surfaceForward(otree_f, odata_f, aic_threshold = 0, exclude = 0,
                       verbose = FALSE, plotaic = FALSE)
k_f = length(fwd_f)
fsum_f = surfaceSummary(fwd_f)
names(fsum_f)
fsum_f$aics
fsum_f$n_regimes
fsum_f
bwd_f = surfaceBackward(otree_f, odata_f, starting_model = fwd_f[[k_f]], aic_threshold = 0,
                        only_best = TRUE, verbose = FALSE, plotaic = FALSE)
bsum_f = surfaceSummary(bwd_f)
kk_f = length(bwd_f)
print("Mis valores alfa para cada rasgo son:")
data.frame(bsum_f$alpha)
print("Mis valores sigma para cada rasgo son:")
data.frame(bsum_f$sigma_squared)
print("Mis valores de aic para cada iteración son:")
data.frame(bsum_f$aics)
print("Mis valores de theta son:")
bsum_f$theta
bsum_f$n_regimes
bsum_f$aics
tiff('surf_f_fullHD.tiff', pointsize=10, width=2800, height=3500, res=600)
par(lwd = 5) 
surf_f<-surfaceTreePlot(dimorfismo.tree, bwd_f[[kk_f]], cols= NULL,  convcol=T, labelshifts = T)# grafica el árbol con los picos adaptativos generados por el SURFACE
surf_f<-recordPlot()

dev.off()
surfaceAICPlot(fwd_f, bwd_f)
abline(h=bm[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic)-2,c("BM","OU1"),cex=0.5)
bm
fitBM_f<-fitContinuous(tree,f_data)
fitBM_f
fitOU_f<-fitContinuous(tree,f_data,model="OU")
fitOU_f
fitEB_f<-fitContinuous(tree,f_data,model="EB")
fitEB_f

plot_surf<-ggarrange(plotlist = myPlots, 
                  labels = c("A.", "B."),
                  ncol = 2, nrow = 1,common.legend = TRUE, legend = "right")
myPlots = list(surf_f, surf_m)
myPlots

plot_surf

surf_f<-as.ggplot(surf_f)
plot(surf_f)

surf_f<-readTIFF("surf_f_fullHD.tiff")
surf_f<-rasterGrob(surf_f)
surf_m<-readTIFF("surf_m_fullHD.tiff")
surf_m<-rasterGrob(surf_m)
surf_f<-as.ggplot(surf_f)
surf_m<-as.ggplot(surf_m)
surf_plot<-ggarrange(surf_f,surf_m,
                     labels = c("A.", "B."),
                     ncol = 2,nrow = 1,common.legend = TRUE, legend = "right")
surf_plot
ggsave("Fig_6_newnew_.tiff", 
              dpi = 600, 
              width = 200, height = 100, unit = "mm")
