#------------------------------------------------------------------------------------------------
#  Purpose: check if there is correlation between Tcell and Bcell immune responses for vtn100 study
#  Creation Date: Oct 6, 2016
#  Lily Zhang
#-------------------------------------------------------------------------------------------------
rm(list=ls())
library(reshape2)
library(ggplot2)
library(plyr)
library(scales)
library(grid)
library(Hmisc)


dir<-"/home/yzhang2/workspace/Tcell_BCell_cor"
# dir<-"T:/vaccine/ImmuneDurability"

dataDir <- file.path(dir,"data")

# read in data
dat_bama_raw <- read.csv(file.path(dataDir,"v100_bama.csv"),stringsAsFactors=F)
dat_ics_raw <- read.csv(file.path(dataDir,"e100fcm_fh_39_2cyt_resp_p.csv"),stringsAsFactors=F)
dat_nab_raw <- read.csv(file.path(dataDir,"v100_nab.csv"),stringsAsFactors=F)
dat_adcc_raw <- read.csv(file.path(dataDir,"s100adf_gf_resp.csv"),stringsAsFactors=F)
dat_phago_raw <- read.csv(file.path(dataDir,"phago_100_rv144.csv"),stringsAsFactors=F)
dat_rx <- read.csv(file.path(dataDir,"vtn100_rv144_pt_stacked.csv"),stringsAsFactors=F)
dat_rx <- rename(dat_rx,c("PTID"="ptid"))
dat_rx <- subset(dat_rx,select=c(ptid,rx_code,pop_PP),protocol=="v100")
dat_rx$ptid <- gsub("\\-","",dat_rx$ptid)
############# immune data processing #######################
#---- BAMA ----------------
dat_bama <- subset(dat_bama_raw,select=c(ptid,isotype,antigen,visitno,dilution,delta,response))
#for IgG, there are three dilutions, 50, 100 and 200
dat_bama$group <- paste0(dat_bama$isotype," ",dat_bama$antigen," ","d",dat_bama$dilution)
dat_bama <- subset(dat_bama,select=c(ptid,group,visitno,delta,response))
dat_bama <- rename(dat_bama,c("delta"="readout"))
dat_bama$measurement <- "delta"
dat_bama$assay <- "BAMA"
#----- ICS -----------------
dat_ics <- subset(dat_ics_raw,select=c(ptid,tcellsub,antigen,cytokine,visitno,pctpos_adj,response))
dat_ics$group <- paste0(dat_ics$tcellsub," ",dat_ics$antigen)
dat_ics <- subset(dat_ics,select=c(ptid,group,visitno,pctpos_adj,response))
dat_ics <- rename(dat_ics,c("pctpos_adj"="readout"))
dat_ics$measurement <- "pctpos_adj"
dat_ics$assay <- "ICS"
#------ NAB ----------------
dat_nab <- subset(dat_nab_raw,select=c(PTID,isolate,visitno,log_titer,response))
dat_nab <- rename(dat_nab,c("PTID"="ptid","isolate"="group","log_titer"="readout"))
dat_nab$measurement <- "log_titer"
dat_nab$assay <- "NAB"
dat_nab$ptid <- gsub("\\-","",dat_nab$ptid)
#------ ADCC ---------------
dat_adcc <- subset(dat_adcc_raw,select=c(ptid,antigen,visitno,activity_peak,response,auc)
                   ,summary_row==1)
dat_adcc <- melt(dat_adcc,id.vars=c("ptid","antigen","visitno","response"))
dat_adcc <- rename(dat_adcc,c("variable"="measurement","antigen"="group","value"="readout"))
dat_adcc <- dat_adcc[,c("ptid","group","visitno","readout","response","measurement")]
dat_adcc$assay <- "ADCC"
dat_adcc$ptid <- gsub("\\-","",dat_adcc$ptid)
####
# note: 
# check the dilution range to make sure different studies used the same range
####
#------- phago --------------
dat_phago <- subset(dat_phago_raw,select=c(ptid,visitno,antigen,avg_phagocytosis_score,response)
                    ,protocol=="100")
dat_phago <- rename(dat_phago,c("antigen"="group","avg_phagocytosis_score"="readout"))
dat_phago$measurement <- "avg_phagocytosis_score"
dat_phago$assay <- "Phago"
dat_phago <- dat_phago[,c("ptid","group","visitno","readout","response","measurement","assay")]
####
# note: only 75 ppts were samples for phago assay, much fewer than other assays
#### 
#------- combine -----------
dat <- rbind(dat_bama,dat_ics,dat_nab,dat_adcc,dat_phago)
dat$group <- paste(dat$assay,dat$group,sep=" ")
dat <- merge(dat,dat_rx,by="ptid")
dat <- subset(dat,rx_code=="T1"&pop_PP==1)
dat$measurement_adcc <- ""
dat$measurement_adcc[dat$measurement=="activity_peak"] <- "peak"
dat$measurement_adcc[dat$measurement=="auc"] <- "auc"
dat$group <- paste(dat$group,dat$measurement_adcc,sep=" ")
#visit 10 
dat_v10 <- subset(dat,visitno==10)

############ start plot ######################
dat_heatmap <- dcast(dat_v10,ptid~group,value.var="readout")

  #calculate correlation
  dat_cor_up <- rcorr(as.matrix(dat_up),type="spearman")$r
  dat_cor_up <- melt(dat_cor_up)
  dat_cor_up <- rename(dat_cor_up,c("value"="cor"))
  #calculate p value
  dat_p_up <- rcorr(as.matrix(dat_up),type="spearman")$P
  dat_p_up <- melt(dat_p_up)
  dat_p_up <- rename(dat_p_up,c("value"="p"))
  #combine correlation with p
  dat_plot_up <- merge(dat_cor_up,dat_p_up,by=c("Var1","Var2"))
  #keep upper triangle
  dat_plot_up <- merge(com_up,dat_plot_up,by=c("Var1","Var2"))
  #-------- lower triangle ----------------
  #calculate correlation
  dat_cor_low <- rcorr(as.matrix(dat_low),type="spearman")$r
  dat_cor_low <- melt(dat_cor_low)
  dat_cor_low <- rename(dat_cor_low,c("value"="cor"))
  #calculate p value
  dat_p_low <- rcorr(as.matrix(dat_low),type="spearman")$P
  dat_p_low <- melt(dat_p_low)
  dat_p_low <- rename(dat_p_low,c("value"="p"))
  #combine correlation with p
  dat_plot_low <- merge(dat_cor_low,dat_p_low,by=c("Var1","Var2"))
  #keep lower triangle
  dat_plot_low <- merge(com_low,dat_plot_low,by=c("Var1","Var2"))
  #---------diagonal-------------
  dat_dia <- data.frame(Var1=group,Var2=group,cor=1,p=NA)
  #------combine placebo and vaccine------
  dat_plot <- rbind(dat_plot_up,dat_plot_low,dat_dia)
  dat_plot$Var1 <- factor(dat_plot$Var1,levels=group)
  dat_plot$Var2 <- factor(dat_plot$Var2,levels=group)
  dat_plot$p_cat <- cut(dat_plot$p, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))
  
  #-------start the plot-------------
  
  p1 <- ggplot(dat_plot,aes(x=Var1,y=Var2))+
    geom_tile(aes(fill=cor),color="black")+
    geom_text(aes(label=p_cat))+
    geom_abline(intercept=0,slope=1)+
    scale_fill_gradient2("",low="blue", mid="white", high="red")+
    xlab("")+
    ylab("")+
    theme(axis.text.x=element_text(angle=45,size=12,color="black",hjust=1)
          ,axis.text.y=element_text(size=12,color="black"))
  p2 <- textGrob("*** p <= 0.001,  ** 0.001 < p <=0.01, * 0.01 < p <= 0.05",x=0.5,y=0.8,gp=gpar(fontsize=11))
  grid.arrange(p1,p2,ncol=1,heights=c(0.95,0.05))
    
}

pdf(file.path(dir,"figure/heatmap_tcellAndAntibody_immuneResponseCorrelation_t1vst3.pdf"),width=10.5,height=9)
f_p(dat_up=dat_imm_t1,dat_low=dat_imm_t3)
dev.off()

pdf(file.path(dir,"figure/heatmap_tcellAndAntibody_immuneResponseCorrelation_t2vst4.pdf"),width=10.5,height=9)
f_p(dat_up=dat_imm_t2,dat_low=dat_imm_t4)
dev.off()
