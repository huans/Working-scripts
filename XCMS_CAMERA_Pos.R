###############################################################################################################
#XCMS-CAMERA
##############################################################################################################
library(xcms)
library(Biobase)
library(multtest)
Rawfiles <- list.files("/public/home/hsi/metabolomics_data/liu_paper/Pos",recursive = TRUE,full.names = TRUE)
xset<-xcmsSet(Rawfiles,
	method="centWave",ppm=25,peakwidth=c(2,10),fitgauss=TRUE, snthresh=10,prefilter=c(3,1000),mzdiff=0.025,
	integrate=1,nSlaves=20,profparam = list(step=0.005))
xset1 <- group(xset,bw=3,minfrac=0.3,mzwid=0.015)#bw=5 or 2
xset2 <- retcor(xset1,method="loess",family ="symmetric",plottype = "none")
xset3 <- group(xset2,bw=1, mzwid=0.015, minfrac=0.3, minsamp=1)
xset4 <- fillPeaks(xset3,method="chrom",nSlaves=20)
#dat <- groupval(xset3, "medret", "into") #get peak intensity matrix
#dat <- rbind(group = as.character(phenoData(xset3)$class), dat) # add group labels (GroupA, GroupB)
#write.csv(dat,file="liu_0925+1011_data.csv")
#reporttab <- diffreport(xset3,"part_1","part_2")#metlin data can be used
#write.csv(reporttab,file="liu_0925+1011_reporttab.csv")
library(CAMERA)
#imports = parent.env(getNamespace("CAMERA"))
#unlockBinding("groups", imports)
#imports[["groups"]] = xcms::groups
#lockBinding("groups", imports)
#diffreport <- annotateDiffreport(xset3,sample=NA, sigma=6, perfwhm=0.6,
#                                 cor_eic_th=0.75, cor_exp_th = 0.75, graphMethod="hcs", pval=0.05, calcCiS=TRUE,
#                                 calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5,
#                                 ppm=5, mzabs=0.015, quick=FALSE, psg_list=NULL, rules=NULL,
#                                 polarity="positive", multiplier=3, max_peaks=100, intval="into",
#                                 pval_th = NULL, fc_th = NULL, sortpval=TRUE,metlin=TRUE)
an  <- xsAnnotate(xset4, nSlaves=4)
calcCasmatrix <- calcCaS(an)
write.csv(calcCasmatrix,"calcCasmatrix_before_groupFWHM.csv")
an1 <- groupFWHM(an, perfwhm=0.6, sigma=6)#perfwhm and sigma were default
an2 <- findIsotopes(an1, mzabs=0.01)
an3 <- groupCorr(an2,calcCiS=TRUE,calcCaS=TRUE)

an4 <- findAdducts(an3, polarity="positive")
##Get peaklist 
peaklist <- getPeaklist(an4, intval="into")
write.csv(peaklist, file='peaklist.csv')


#write.csv(diffreport,"liu_0925+1011_camera.csv")

