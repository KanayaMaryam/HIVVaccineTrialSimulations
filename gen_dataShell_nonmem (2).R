##################################################################################
#Purpose:create 1000 data shells for simulations of 1000 datasets of concentrations 
#        by NONMEM accounting for: 
#        (1) the expected distributions of AMP participants' characteristics including
#            sex, weight, and age 
#        (2) missed and/or terminated infusions and visits, as well as dynamic visit schedules and 
#            visit windows.
#        (3) time to infection using Cox models with cyclic time-varying covariates

#Notes: (1) the whole cohort will be generated, concentrations will be simulated only for 
#       ppts selected for case/control 
#       (2) including two cohorts (male and female) and two dose levels (10 mg/kg and 30 mg/kg) 
#       (3) assuming no HIV diagnostic test at 5 days post 2nd infusion

#Author: Lily Zhang
#Date  :  Jan 18, 2017 
##################################################################################

rm(list=ls())
library(plyr)
library(reshape2)
library(parallel)
library(zoo)

###################
#simulation function
##################
#' Simulation
#' @param n_infu  : number of infusions
#' @param prob_misInfu   : probability of missing infusions at each infusion visit
#' @param prob_misVis : probability of missing measurement visit
#' @param r_infuDiscont : annual rate of permanent infusion discontinuation
#' @param n_ppt  : number of ppts
#' @param prob_outWin  : probability of out of target window but within 
#'                       allowable window for infusion visits 
#' @param prob_disContInfu_Win : probability of out of target window but within 
#'                       allowable window for visits after infusion permanent discontinuation
#' @param r_dropout  : drop out annual rate
#' @param beta.t  : time-dependent covariate Cox model coefficient
#' @param r_incidence : annual incidence rate 
#' @param B  : number of datasets simulated 


sim <- function(n_infu=10
                ,prob_outWin=0
#                 ,prob_outWin=0.2
                ,prob_disContInfu_Win=0
#                 ,prob_disContInfu_Win=0.2
                ,n_male
                ,n_female
                ,prob_misInfu
                ,prob_misVis
                ,r_infuDiscont
                ,r_dropout
                ,beta.t.30
                ,beta.t.10
                ,r_incidence = 0.04
                ,pdf_inf = "exponential"
                ,B=1000){
  #### redefine some parameters #####
  n_ppt <- n_male + n_female
  prob_misInfu <- rep(prob_misInfu,n_infu)
  prob_outWin <- rep(prob_outWin,n_infu)
  prob_disContInfu_Win <- c(prob_disContInfu_Win/2,1-prob_disContInfu_Win,prob_disContInfu_Win/2)
  #case to ctl ratios
  if(r_dropout==0){
    ctl_to_case <- 0
  }else{
    ctl_to_case <- 0
  }
  #directory to save the current scenario data
  folder <- paste0("m_",n_ppt
                   ,"_p1_",unique(prob_misInfu)
                   ,"_p2_",prob_misVis
                   ,"_r1_",r_infuDiscont
                   ,"_r2_",r_dropout
                   ,"_b30_",beta.t.30
                   ,"_b10_",beta.t.10
                   ,"_r3_",r_incidence
                   )                      
  ##### start the simulation #####
  ind.0cases <- lapply(1:B
           ,function(b){
    set.seed(b)
    #-------- infusion visits -------------------
    dat <- data.frame(ID=1:n_ppt)
    dat$infu1 <- 1
    for(i in 2:n_infu){
      var_name <- paste0("infu",i)
      dat[,var_name] <- rbinom(n=n_ppt,size=1,prob=1-prob_misInfu[i])
    }
    #long format
    dat <- melt(dat,id.vars="ID")
    dat <- dat[order(dat$ID),]
    dat <- rename(dat,c("value"="f_infu","variable"="infusion"))
    dat$infusion <- as.numeric(gsub("infu","",dat$infusion))
    #out of target window but within allowable window
    dat <- ddply(dat, .(infusion), function(df){
      if(unique(df$infusion)==1){
        df$out_win <- NA
      }else{
        infu <- unique(df$infusion)
        dat_infu <- subset(df,f_infu==1)
        dat_infu$out_win <- rbinom(n=nrow(dat_infu),size=1,prob=prob_outWin[infu])
        dat_misInfu <- subset(df,f_infu==0)
        if(nrow(dat_misInfu)==0){
          df <- dat_infu
        }else{
          dat_misInfu$out_win <- NA
          df <- rbind(dat_infu,dat_misInfu)
        }
      }
      return(df)
    })
    #random draw from windows
    dat_inWin <- subset(dat,out_win==0)
    dat_inWin$window <- round(runif(n=nrow(dat_inWin),min=-7,max=7),0)
#     dat_inWin$window <- 0 #test perfect aherence no window
    dat_outWin <- subset(dat,out_win==1)
    dat_outWin$window <- round(runif(n=nrow(dat_outWin),min=8,max=48),0)
    dat_naWin <- subset(dat,is.na(out_win))
    dat_naWin$window <- 0
    dat <- rbind(dat_inWin,dat_outWin,dat_naWin)
    dat <- dat[order(dat$ID,dat$infusion),]
    #infusion day
    dat <- ddply(dat,.(ID),function(df){
      #   browser()
      df$day_infu[df$infusion==1] <- 0
      for(i in 2:n_infu){
        df$day_infu[df$infusion==i] <- df$day_infu[df$infusion==(i-1)]+56+df$window[df$infusion==i]
      }
      return(df)
    })
    dat$day_infu <- ifelse(dat$f_infu==0,NA,dat$day_infu)
    dat$out_win <- NULL
    dat$window <- NULL
    #------------ post-infusion visits ------------------
#     browser()
    dat$window <- round(runif(n=nrow(dat),min=-7,max=7),0)
    dat$day_4wkPostInfu <- dat$day_infu+28+dat$window
    dat$day_8wkPostInfu10 <- dat$day_infu+56+dat$window
    dat$window_5dayPostInfu2 <- round(runif(n=nrow(dat),min=-2,max=2),0)
    dat$day_5dayPostInfu2 <- dat$day_infu+5+dat$window_5dayPostInfu2
    #missing visit 
    dat <- ddply(dat,.(infusion),function(df){
      nppt_infu <- length(df$f_infu[df$f_infu==1])
      df$f_4wkPostInfuVis <- ifelse(df$f_infu==1,rbinom(n=nppt_infu,size=1,prob=1-prob_misVis),0)
      
      if(unique(df$infusion)==2){
        df$f_5dayPostInfu2Vis <- ifelse(df$f_infu==1,rbinom(n=nppt_infu,size=1,prob=1-prob_misVis),0)
      }else{
        df$f_5dayPostInfu2Vis <- 0
      }
      
      if(unique(df$infusion)==10){
        df$f_8wkPostInfu10Vis <- rbinom(n=nrow(df),size=1,prob=1-prob_misVis)
      }else{
        df$f_8wkPostInfu10Vis <- 0
      }
      return(df)
      
    })
    dat$day_4wkPostInfu <- ifelse(dat$f_4wkPostInfuVis==0,NA,dat$day_4wkPostInfu)
    dat$day_8wkPostInfu10 <- ifelse(dat$f_8wkPostInfu10Vis==0,NA,dat$day_8wkPostInfu10)
    dat$day_5dayPostInfu2 <- ifelse(dat$f_5dayPostInfu2Vis==0,NA,dat$day_5dayPostInfu2)
    dat$window <- NULL
    dat$window_5dayPostInfu2 <- NULL
    #------------ change data to long format ----------------
    #infusion visit
    dat_infu <- subset(dat,select=c(ID,infusion,f_infu,day_infu))
    dat_infu <- rename(dat_infu,c("infusion"="visit_n","f_infu"="f_vis","day_infu"="day_vis"))
    dat_infu$visit_n <- ifelse(dat_infu$visit_n%in%c(1,2),dat_infu$visit_n*2,dat_infu$visit_n*2+1)
    dat_infu$f_infu <- 1
    #4 week post infusion
    dat_4wkPostInfu <- subset(dat,select=c(ID,infusion,f_4wkPostInfuVis,day_4wkPostInfu))
    dat_4wkPostInfu <- rename(dat_4wkPostInfu,c("infusion"="visit_n"
                                                ,"f_4wkPostInfuVis"="f_vis"
                                                ,"day_4wkPostInfu"="day_vis"))
    dat_4wkPostInfu$visit_n <- ifelse(dat_4wkPostInfu$visit_n==1,3,dat_4wkPostInfu$visit_n*2+2)
    dat_4wkPostInfu$f_infu <- 0
    #8 week post 10th infusion
    dat_8wkPostInfu10 <- subset(dat
                              ,infusion==10
                              ,select=c(ID,infusion,f_8wkPostInfu10Vis,day_8wkPostInfu10))
    dat_8wkPostInfu10 <- rename(dat_8wkPostInfu10,c("infusion"="visit_n"
                                                      ,"f_8wkPostInfu10Vis"="f_vis"
                                                      ,"day_8wkPostInfu10"="day_vis"))
    dat_8wkPostInfu10$visit_n <- 23
    dat_8wkPostInfu10$f_infu <- 0
    # 5 days post 2nd infusion 
    dat_5dayPostInfu2 <- subset(dat
                                ,infusion==2
                                ,select=c(ID,infusion,f_5dayPostInfu2Vis,day_5dayPostInfu2))
    dat_5dayPostInfu2 <- rename(dat_5dayPostInfu2,c("infusion"="visit_n"
                                                    ,"f_5dayPostInfu2Vis"="f_vis"
                                                    ,"day_5dayPostInfu2"="day_vis"))
    dat_5dayPostInfu2$visit_n <- 5
    dat_5dayPostInfu2$f_infu <- 0
    #combine
    dat_long <- rbind(dat_infu,dat_5dayPostInfu2,dat_4wkPostInfu,dat_8wkPostInfu10)
    dat_long <- subset(dat_long,f_vis==1)
    dat_long <- dat_long[order(dat_long$ID,dat_long$visit_n),]
    #----------- permanent infusion discontinuation -------------------------------------------------
    #time to discontinuation of infusion
    if(r_infuDiscont==0){
      r_infuDiscont <- 0.00000001
    }
    time_discontInfu <- data.frame(ID=1:n_ppt,day_discontInfu=round(rexp(n=n_ppt,rate=r_infuDiscont/365),0))
    #don't allow discontiune at time 0
    time_discontInfu$day_discontInfu[time_discontInfu$day_discontInfu==0] <- 1
    dat_long <- merge(dat_long,time_discontInfu,by="ID")
    dat_long$f_discontInfu <- ifelse(dat_long$day_vis<dat_long$day_discontInfu,0,1)
    dat_long$f_vis <- ifelse(dat_long$f_discontInfu==1,0,dat_long$f_vis)
    dat_long$day_vis <- ifelse(dat_long$f_discontInfu==1, NA,dat_long$day_vis)
    # visits after the ppts discontinue infusions
    id_discontInfu <- unique(dat_long$ID[dat_long$f_discontInfu==1])
    dat_discont <- unique(dat_long[dat_long$ID%in%id_discontInfu,c("ID","day_discontInfu")])
    dat_discont <- data.frame(ID=rep(dat_discont$ID,each=6)
                              ,day_discontInfu=rep(dat_discont$day_discontInfu,each=6)
                              ,wk_post_enroll = rep(c("20wk","32wk","44wk","56wk","68wk","80wk"),nrow(dat_discont))
    ) 
    dat_discont$lw[dat_discont$wk_post_enroll=="20wk"] <- 99
    dat_discont$lw[dat_discont$wk_post_enroll=="32wk"] <- 183
    dat_discont$lw[dat_discont$wk_post_enroll=="44wk"] <- 267
    dat_discont$lw[dat_discont$wk_post_enroll=="56wk"] <- 351
    dat_discont$lw[dat_discont$wk_post_enroll=="68wk"] <- 435
    dat_discont$lw[dat_discont$wk_post_enroll=="80wk"] <- 519
    
    dat_discont$day_tar[dat_discont$wk_post_enroll=="20wk"] <- 140
    dat_discont$day_tar[dat_discont$wk_post_enroll=="32wk"] <- 224
    dat_discont$day_tar[dat_discont$wk_post_enroll=="44wk"] <- 308
    dat_discont$day_tar[dat_discont$wk_post_enroll=="56wk"] <- 392
    dat_discont$day_tar[dat_discont$wk_post_enroll=="68wk"] <- 476
    dat_discont$day_tar[dat_discont$wk_post_enroll=="80wk"] <- 560
    
    dat_discont$f_discont <- ifelse(dat_discont$day_discontInfu<dat_discont$lw,1,0)
    dat_discont <- subset(dat_discont,f_discont==1)
    
    win_type <- data.frame(t(rmultinom(nrow(dat_discont), size = 1, prob=prob_disContInfu_Win)))
    names(win_type) <- c("win_type1","win_type2","win_type3")
    dat_discont <- cbind(dat_discont,win_type)
    dat_discont$win1 <- round(runif(n=nrow(dat_discont),min=-41,max=-14),0)
    dat_discont$win2 <- round(runif(n=nrow(dat_discont),min=-14,max=14),0)
    dat_discont$win3 <- round(runif(n=nrow(dat_discont),min=14,max=42),0)
    dat_discont$win[dat_discont$win_type1==1] <- dat_discont$win1[dat_discont$win_type1==1]
    dat_discont$win[dat_discont$win_type2==1] <- dat_discont$win2[dat_discont$win_type2==1]
    dat_discont$win[dat_discont$win_type3==1] <- dat_discont$win3[dat_discont$win_type3==1]
    dat_discont$day_vis <- dat_discont$day_tar+dat_discont$win
    #missing visit
    dat_discont$f_vis <- rbinom(n=nrow(dat_discont),size=1,prob=1-prob_misVis)
    dat_discont$day_vis <- ifelse(dat_discont$f_vis==0,NA,dat_discont$day_vis)
    #visits for ppt who discontinue infusions 
    dat_discont <- subset(dat_discont,select=c(ID,wk_post_enroll,f_vis,day_vis))
    dat_discont <- rename(dat_discont,c("wk_post_enroll"="visit_n"))
    dat_discont$visit_n  <- as.character(dat_discont$visit_n)
    dat_discont$visit_n[dat_discont$visit_n=="20wk"] <- 72
    dat_discont$visit_n[dat_discont$visit_n=="32wk"] <- 73
    dat_discont$visit_n[dat_discont$visit_n=="44wk"] <- 74
    dat_discont$visit_n[dat_discont$visit_n=="56wk"] <- 75
    dat_discont$visit_n[dat_discont$visit_n=="68wk"] <- 76
    dat_discont$visit_n[dat_discont$visit_n=="80wk"] <- 77
    dat_discont$f_infu[!is.na(dat_discont$ID)] <- 0 #for situation that no ppts who discontinue infusions
    #combine
    dat_long$day_discontInfu <- NULL
    dat_long$f_discontInfu <- NULL
    dat_long <- rbind(dat_long,dat_discont)
    dat_long <- subset(dat_long,f_vis==1)
    #------------------- drop out ----------------------
    if(r_dropout==0){
      r_dropout <- 0.000000001
    }
    time_dropout <- data.frame(ID=1:n_ppt, time_dropout=rexp(n=n_ppt,rate = r_dropout/365))
    dat_long <- merge(dat_long,time_dropout)
    dat_long$f_vis <- ifelse(dat_long$day_vis>=dat_long$time_dropout,0, dat_long$f_vis)
    dat_long$day_vis <- ifelse(dat_long$day_vis>=dat_long$time_dropout,NA, dat_long$day_vis)
    dat_long <- dat_long[order(as.numeric(dat_long$visit_n)),]
    #final data
    dat_long <- subset(dat_long,f_vis==1)
    dat_long$time_dropout <- NULL
  
    #------- generate covariates (age, gender, sex and dose levels)--------
    #get the distribution of age from 503 female and 502 male 
    dataDirAge <- "C:\\Users\\zzhang2\\Desktop\\simdat\\"
    #   dataDir <- "T:/vaccine/Ad5Meta/data/"
    data502MaleFile<-"v502_meta_males_16JUL2013.csv"
    data503File<-"503_meta_16Sep2013.csv"
    dat_502_male_raw<-read.csv(paste(dataDirAge,data502MaleFile, sep=""),header=T)
    dat_503_raw<-read.csv(paste(dataDirAge,data503File, sep=""),header=T)
    dat_503_female <- subset(dat_503_raw,DEMsex=="Female")
    #calculate mean and sd of age for male and female respectively
    m_age_male <- mean(dat_502_male_raw$AGE)
    sd_age_male <- sd(dat_502_male_raw$AGE)
    m_age_female <- mean(dat_503_female$age)
    sd_age_female <- sd(dat_503_female$age)
    #WT data
    dataDirWT <- "C:\\Users\\zzhang2\\Desktop\\simdat\\"
    dat_502_WT <- read.csv(file.path(dataDirWT,"v502_wt_list.csv"))
    dat_502_WT <- unique(dat_502_WT)
    dat_503_WT <- read.csv(file.path(dataDirWT,"503_blweight_women.csv"))
    #WT, age and dose
    dat_m <- data.frame(ID=1:n_male
                        ,WT=sample(dat_502_WT$WT_KLG,n_male,replace=T)
                        ,AGE = round(rnorm(n=n_male,mean=m_age_male,sd=sd_age_male),0)
                        ,sex=1
                        ,dose=c(rep(10,floor(n_male/2)),rep(30,ceiling(n_male/2)))
    )
    
    dat_m$AGE <- ifelse(dat_m$AGE>50,50,dat_m$AGE)
    dat_m$AGE <- ifelse(dat_m$AGE<18,18,dat_m$AGE)
    dat_f <- data.frame(ID=(n_male+1):(n_female+n_male)
                        ,WT=sample(dat_503_WT$weight,n_female,replace=T)
                        ,AGE = round(rnorm(n=n_female,mean=m_age_female,sd=sd_age_female),0)
                        ,sex=0
                        ,dose=c(rep(10,floor(n_female/2)),rep(30,ceiling(n_female/2)))
    )
    
    dat_f$AGE <- ifelse(dat_f$AGE>40,40,dat_f$AGE)
    dat_f$AGE <- ifelse(dat_f$AGE<18,18,dat_f$AGE)
    dat_cov <- rbind(dat_m,dat_f)
    dat_long <- merge(dat_cov,dat_long,by="ID")
    #---------------- infection status ------------------
    if(beta.t.30+beta.t.10==0){
      if(pdf_inf=="exponential"){
        t_event <- data.frame(ID=1:n_ppt,t_event=rexp(n=n_ppt,rate=r_incidence/365))
        dat_long <- merge(dat_long,t_event,by="ID")
        dat_long$RNA <- ifelse(dat_long$day_vis>=dat_long$t_event,1,0)
        ID_cases <- unique(dat_long$ID[dat_long$RNA==1])
        dat_long$t_event <- ifelse(dat_long$ID%in%ID_cases,dat_long$t_event,NA)
#         dat_long$t_event <- NULL
      }
      if(pdf_inf=="binomial"){
        ind_cases <- data.frame(ID=1:n_ppt,f_cases=rbinom(n=n_ppt,size=1,prob=0.055))
        dat_long <- merge(dat_long,ind_cases,by="ID")
        dat_cases <- subset(dat_long,f_cases==1)
        dat_cases <- ddply(dat_cases,.(ID),function(df){
          #       browser()
          b <- max(df$day_vis)
          df$t_event <- round(runif(1,min=1,max=b),0)
          df$RNA <- ifelse(df$t_event<=df$day_vis,1,0)
          return(df)
        })
#         dat_cases$t_event <- NULL
        dat_control <- subset(dat_long,f_cases==0)
        dat_control$RNA <- 0
        dat_long <- rbind(dat_cases,dat_control)
        dat_long$f_cases <- NULL
      }
    }else{
      f_infect_mulDose_imperf <- function(data, beta.t, ts){
        lambda <- r_incidence/365/exp(beta.t*56*2)
        data <- ddply(data,.(ID),function(df){
          dat <- df
          dat <- dat[order(dat$day_vis),]
          dat_infu <- subset(dat,f_infu==1)
          dat_infu$u <- runif(1)
          dat_infu$tDiff <- dat_infu$day_vis-dplyr::lag(dat_infu$day_vis)
          dat_infu$past_thresh <- ifelse(dat_infu$tDiff > ts, 1, 0)
          dat_infu$exp_beta_tDiff <- ifelse(dat_infu$past_thresh == 0, exp(beta.t*dat_infu$tDiff), exp(beta.t*ts) + beta.t*(dat_infu$tDiff-ts)*exp(beta.t*ts))
          dat_infu$exp_beta_tDiff[dat_infu$day_vis==0] <- 0
          dat_infu$sum_exp_beta_tDiff <- cumsum(dat_infu$exp_beta_tDiff)
          dat_infu$k <- 1:nrow(dat_infu)

          dat_infu$a <- lambda*(dat_infu$sum_exp_beta_tDiff-(dat_infu$k-1))/beta.t
          dat_infu$b <- lambda*(dat_infu$sum_exp_beta_tDiff + exp(beta.t*ts) - dat_infu$k)/beta.t
          dat_infu$c <- dplyr::lead(dat_infu$a)
          dat_infu$neg_log_u <- -log(dat_infu$u)
          dat_infu$event_k <- ifelse(dat_infu$neg_log_u<dat_infu$c&dat_infu$neg_log_u>=dat_infu$a,1,0)
          dat_infu$event_k <- ifelse(dat_infu$neg_log_u<dat_infu$c&dat_infu$neg_log_u>=dat_infu$b,2,dat_infu$event_k)
          dat_infu$event_k[dat_infu$neg_log_u>=dat_infu$a&is.na(dat_infu$c)] <- 1
          dat_infu$event_k[dat_infu$neg_log_u>=dat_infu$b&is.na(dat_infu$c)] <- 2
          
          
          dat_infu$t_event <- log(exp(beta.t*dat_infu$day_vis)*(beta.t*dat_infu$neg_log_u/lambda-dat_infu$sum_exp_beta_tDiff+dat_infu$k))/beta.t
          dat_infu$t_event <- ifelse(dat_infu$event_k == 2, (beta.t*dat_infu$neg_log_u/lambda-dat_infu$sum_exp_beta_tDiff - exp(beta.t*ts)+dat_infu$event_k)/beta.t/exp(beta.t*ts) + ts, dat_infu$t_event)

          t_event <- dat_infu$t_event[dat_infu$event_k==1 | dat_infu$event_k==2]

          dat$t_event <- t_event
          
          dat$RNA <- ifelse(dat$t_event<=dat$day_vis,1,0)
          #           dat$t_event <- NULL
          return(dat)
        })
        return(data)
      }
      #survival time is simulated by considering multiple dose infusions, limited to perfect adherence
      f_infect_mulDose <- function(data,beta.t){
        lambda <- r_incidence/365/exp(beta.t*56*2)
        
        data <- ddply(data,.(ID),function(df){
#                 browser()
          dat <- df
          dat <- dat[order(dat$day_vis),]
          
          dat_infu <- subset(dat,f_infu==1)
          dat_infu$u <- runif(1)
          dat_infu$tDiff <- dat_infu$day_vis-dplyr::lag(dat_infu$day_vis)
          dat_infu$exp_beta_tDiff <- exp(beta.t*dat_infu$tDiff)
          dat_infu$exp_beta_tDiff[dat_infu$day_vis==0] <- 0
          dat_infu$sum_exp_beta_tDiff <- cumsum(dat_infu$exp_beta_tDiff)
          dat_infu$k <- 1:nrow(dat_infu)
          dat_infu$a <- lambda*(dat_infu$sum_exp_beta_tDiff-(dat_infu$k-1))/beta.t
          dat_infu$b <- dplyr::lead(dat_infu$a)
          dat_infu$neg_log_u <- -log(dat_infu$u)
          dat_infu$event_k <- ifelse(dat_infu$neg_log_u<dat_infu$b&dat_infu$neg_log_u>=dat_infu$a,1,0)
          dat_infu$event_k[dat_infu$neg_log_u>=dat_infu$a&is.na(dat_infu$b)] <- 1
          dat_infu$t_event <- log(exp(beta.t*dat_infu$day_vis)*(beta.t*dat_infu$neg_log_u/lambda-dat_infu$sum_exp_beta_tDiff+dat_infu$k))/beta.t
          t_event <- dat_infu$t_event[dat_infu$event_k==1]
          
          dat$t_event <- t_event
          dat$RNA <- ifelse(dat$t_event<=dat$day_vis,1,0)
#           dat$t_event <- NULL
          return(dat)
        })
        return(data)
      }
      #survival time is simulated by considering single dose infusion
      f_infect_singleDose <- function(data,beta.t,ts){
        lambda <- r_incidence/365/exp(beta.t*ts)
        
        
        ddply(data,.(ID),function(df){
          #           browser()
          dat <- df
          dat <- dat[order(dat$day_vis),]
          #infection status
          dat_infu <- subset(dat,f_infu==1)
          if(dat$f_infu[nrow(dat)]==0){
            dat_lastVis <- dat[nrow(dat),]
            dat_infu <- rbind(dat_infu,dat_lastVis)
          }
          dat_infu$neg_log_u <- -log(runif(nrow(dat_infu)))
          dat_infu$cond <- lambda/beta.t*(exp(beta.t*ts)-1)
          dat_infu$t_event <- ifelse(dat_infu$neg_log_u<dat_infu$cond
                                     ,log(1+beta.t*dat_infu$neg_log_u/lambda)/beta.t
                                     ,dat_infu$neg_log_u/(lambda*exp(beta.t*ts))+(1-exp(beta.t*ts))/(beta.t*exp(beta.t*ts))+ts)
          
          dat_infu$t_interval <- dplyr::lead(dat_infu$day_vis)-dat_infu$day_vis
          dat_infu$f_event <- ifelse(dat_infu$t_interval>=dat_infu$t_event,1,0)
          #controls
          if(all(dat_infu$f_event!=1,na.rm=T)){
            dat$t_event <- NA
            dat$RNA <- 0
          }
          #cases
          if(any(dat_infu$f_event==1,na.rm=T)){
            #             browser()
            dat_event <- subset(dat_infu,f_event==1)
            dat_event <- dplyr::sample_n(dat_event,1)
            dat$t_event <- dat_event$t_event+dat_event$day_vis
            dat$RNA <- ifelse(dat$t_event<= dat$day_vis,1,0)
            #             dat$t_event <- NULL
          }
          return(dat)
        })
        
      }
      #dat_long_d10 <- f_infect_mulDose(subset(dat_long,dose==10),beta.t.10)
      #dat_long_d30 <- f_infect_mulDose(subset(dat_long,dose==30),beta.t.30)
      
      #dat_long_d10 <- f_infect_mulDose_imperf(subset(dat_long,dose==10),beta.t.10, ts=90)
      #dat_long_d30 <- f_infect_mulDose_imperf(subset(dat_long,dose==30),beta.t.30, ts=114)
      
      dat_long_d10 <- f_infect_singleDose(subset(dat_long,dose==10),beta.t.10, ts=90)
      dat_long_d30 <- f_infect_singleDose(subset(dat_long,dose==30),beta.t.30, ts=114)
      dat_long <- rbind(dat_long_d10,dat_long_d30)
      
    }
   
#       browser()
    ID_cases <- unique(dat_long$ID[dat_long$RNA==1])

    if (length(ID_cases)!=0) {
      # truncate visits after first RNA positive
      dat_long_cases <- subset(dat_long,ID%in%ID_cases)
      dat_long_cases <- ddply(dat_long_cases,.(ID),function(df){
#               browser()
        dat <- df
        dat <- dat[order(dat$day_vis),]
        dat$f_firstRNAP <- ifelse(dplyr::lag(dat$RNA)==0&dat$RNA==1,1,0)
        #assume no hiv diagnostic test at 5 days post 2nd infusion
        if(dat$visit_n[dat$f_firstRNAP==1]==5){
          dat$RNA[dat$f_firstRNAP==1] <- 0
          day_firstRNAP <- dat$day_vis[which(dat$f_firstRNAP==1)+1]
        }else{
          day_firstRNAP <- dat$day_vis[dat$f_firstRNAP==1]
        }
        dat <- subset(dat,day_vis<=day_firstRNAP)
        dat$f_firstRNAP <- NULL
        return(dat)
      })
      
      dat_long_ctl <- subset(dat_long,!ID%in%ID_cases)
      dat_long <- rbind(dat_long_cases,dat_long_ctl)
      dat_long$visit_n <- as.numeric(dat_long$visit_n)
      dat_long <- arrange(dat_long,ID,visit_n)
      #     browser()
      #------- caculate PE -----------------------------------------------
      n_cases <- length(ID_cases)
      #     browser()
      t_fu <- dplyr::summarise(dplyr::group_by(dat_long,ID),max_day=max(day_vis))
      dat_long$PE <- 1- (n_cases/(sum(t_fu$max_day)/360))/r_incidence
      dat_long$f_discontInfu <- NULL
              #---- save data for approach 1&3 of correlate paper without including concentration information----
      dat_long <- rename(dat_long,c("day_vis"="TIME"))
      dat_long <- subset(dat_long
                                           ,select=c(ID,sex,dose,visit_n,f_vis,TIME,f_infu,RNA,PE,t_event)
                                           ,visit_n!=5)
      
              day_max <- dplyr::summarise(dplyr::group_by(dat_long,ID),day_max=max(TIME))
              dat_com <-ddply(day_max,.(ID),function(df){
                #               browser()
                if(df$day_max!=0){
                  dat <- data.frame(ID=df$ID,TIME=1:df$day_max)
                  return(dat)
                }
              })
                      #merge visit information
                      dat_shell <- dat_long
                      dat_shell_infu <- subset(dat_shell,f_infu==1)
                      #trough measurement
                      dat_trgh <- subset(dat_shell_infu, visit_n!=2)
                      dat_shell_sample <- subset(dat_shell,f_infu==0)
                      dat_shell_sample <- rbind(dat_shell_sample,dat_trgh)
                      dat_shell_sample <- merge(dat_com,dat_shell_sample,by=c("ID","TIME"),all.x=T)
                      dat_shell_sample$f_vis[is.na(dat_shell_sample$f_vis)] <- 0
                      dat_shell_sample$f_infu[is.na(dat_shell_sample$f_infu)] <- 0
                      dat_shell_sample$RNA[is.na(dat_shell_sample$RNA)] <- 0

                      dat_shell <- rbind(dat_shell_infu,dat_shell_sample)

                      #     browser()
                      dat_shell <- dat_shell[order(dat_shell$ID,dat_shell$TIME),]


                      dat_shell$sex <- na.locf(dat_shell$sex)
                      dat_shell$dose <- na.locf(dat_shell$dose)
                      dat_shell$PE <- na.locf(dat_shell$PE)
                      dat_shell <- ddply(dat_shell, .(ID),function(df){
                        df$t_event <- df$t_event[df$TIME==0]
                        return(df)
                      })
                      dat_shell <- dat_shell[,c("ID","sex","dose","visit_n","f_vis","TIME","f_infu","RNA","PE","t_event")]
                                              
      dat_long = unique(dat_shell)
              
      folder <- paste0(folder,"_ratio_",0)
      dir_noConc <- file.path("C:\\Users\\zzhang2\\Desktop\\simdat\\",folder)
      csv.file <- paste0(b,".csv")
     if(file.exists(dir_noConc)==FALSE){
      dir.create(dir_noConc,recursive=T)   
      }
      write.csv(dat_long
                ,file=file.path(dir_noConc,csv.file)
                ,row.names = FALSE)
      
    }else{
        ind.0cases <- data.frame(b=b)
        return(ind.0cases)
    }
  })
}

###########################################################
# create different scenarios and run the simulation function
###########################################################
n_ppt <- list(c(1800,1000))
adherence <- list(
                 c(0.02,0.03,0.03,0.05)
               # ,c(0.1,0.15,0.15,0.15)
)
beta.t <- list(
  #c(0.02,0.025)
              # c(0.025,0.03)
               #,c(0.03,0.038)
               #,c(0.01,0.012)
               c(0,0) #currently running this one
)
scenario <- expand.grid(n_ppt,adherence,beta.t)
scenario <- apply(scenario,2,function(x)do.call("rbind",x))
scenario <- do.call("cbind",scenario)
scenario <- data.frame(scenario)
names(scenario) <- c("n_male","n_female","prob_misInfu","prob_misVis","r_infuDiscont","r_dropout","beta.t.30","beta.t.10")
print(scenario)

for(i in 1:nrow(scenario)){
  sim(n_male=scenario$n_male[i]
      ,n_female=scenario$n_female[i]
      ,prob_misInfu=scenario$prob_misInfu[i]
      ,prob_misVis=scenario$prob_misVis[i]
      ,r_infuDiscont=scenario$r_infuDiscont[i]
      ,r_dropout=scenario$r_dropout[i]
      ,beta.t.30=scenario$beta.t.30[i]
      ,beta.t.10=scenario$beta.t.10[i])
}

