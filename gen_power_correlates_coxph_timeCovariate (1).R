#-----------------------------------------------------------
# Erika Thommes
#
# Assessing PK Marker Correlates of HIV-Infection in AMP
# Program to Generate Simulation Results under 4 approaches
#
# 1. Approach 1 = Wald: mu_case - mu_control/2
# 2. Approach 2 = Wald: mu_s_controls - mu_s_cases
# 3. Approach 3 = PH Model with Time-Dependent Covariate
# 4. Approach 4 = Sign Test
#-----------------------------------------------------------

#-----------------------------------------------------------
library(zoo)
library(plyr)
library(reshape)
library(survival)

sysdate=tolower(gsub('([[:punct:]])|\\s+', '', format(Sys.Date(), "%d %b %y")))

in_dswd="C:\\Users\\zzhang2\\Desktop\\simdat\\"
out_dswd="C:\\Users\\zzhang2\\Desktop\\simdat\\"

ds1="m_2800_p1_0.02_p2_0.03_r1_0.03_r2_0.05_b30_0.02_b10_0.025_r3_0.04_ratio_0"
#ds2="m_2800_p1_0_p2_0_r1_0_r2_0_b30_0.03_b10_0.038_r3_0.04_ratio_0"
#ds3="m_2800_p1_0_p2_0_r1_0_r2_0_b30_0.025_b10_0.03_r3_0.04_ratio_0"
#dsn="m_2800_p1_0_p2_0_r1_0_r2_0_b30_0.02_b10_0.025_r3_0.04_ratio_0"
dsn=c(ds1)
#-----------------------------------------------------------

#-----------------------------------------------------------
# Define some useful functions
#-----------------------------------------------------------

week4_visits=c(3,6,8,10,12,14,16,18,20,22)
infu_visits=c(2,4,7,9,11,13,15,17,19,21,23)

calc_lag <- function(x){
  c(NA, x[1:length(x)-1])
}

calc_lead <- function(x){
  c(x[2:length(x)], NA)
}

calc_base_infu <- function(x){
  temp1=ifelse(x %in% infu_visits, x, 0)
  temp2=rep(0, length(x))
  for(i in 1:length(x)){
    temp2[i]=ifelse(temp1[i] == 0, temp2[i-1]+.01, temp1[i])
  }
  floor(temp2)
}

calc_grid <- function(t, tmin, tmax){
  temp=rep(0, length(t))
  for(i in 1:length(t)){
    temp[i]=ifelse(t[i] %in% seq(from=tmin[i], to=tmax[i], by=7), 1, 0)
  }
  temp
}

#-----------------------------------------------------------
# Define Approach 3 
#-----------------------------------------------------------

prep_approach3 <- function(dat, case_input){
  print("in prep approach")
  #format t_event: factor --> numeric (suppress warnings from missing values)
  dat$visit_n_num=suppressWarnings(as.numeric(as.character(dat$visit_n)))
  dat$t_event_num=suppressWarnings(as.numeric(as.character(dat$t_event)))
  # 
  # #replace estimated IPRE with observed concentration at visits
  # #replace values <= 0 with small non-negative value
  # dat$DV_obs=suppressWarnings(as.numeric(as.character(dat$DV_obs)))
  # dat$IPRE=ifelse(is.na(dat$DV_obs), dat$IPRE, dat$DV_obs)
  # dat$IPRE=ifelse(dat$IPRE <= 0, 0.00001, dat$IPRE)
  
  #exclude ppts who discontinued
  sched4ppts=unique(dat[dat$visit_n_num > 70,]$ID)
  dat_sched1ppts=dat[!dat$ID %in% sched4ppts,]
  dat_sched1ppts$temp <- dat_sched1ppts$TIME*dat_sched1ppts$f_infu
  dat_sched1ppts$temp <- ifelse(dat_sched1ppts$temp==0&dat_sched1ppts$TIME!=0,NA,dat_sched1ppts$temp)
  dat_sched1ppts$temp <- na.locf(dat_sched1ppts$temp)
  dat_sched1ppts$t_IPRE <- dat_sched1ppts$TIME-dat_sched1ppts$temp

  #subset data into cases and controls
  cases=unique(dat_sched1ppts[dat_sched1ppts$RNA==1,]$ID)
  dat1=dat_sched1ppts[dat_sched1ppts$ID %in% cases ,]
  dat0=dat_sched1ppts[!dat_sched1ppts$ID %in% cases,]
  
  # #----- cases ----- 
  # 
  # if(case_input=="est"){
  #   
  #   #exclude v5 data (5-day post infu #2) and calculate t_star
  #   dat1_subset=dat1[dat1$f_vis==1 & dat1$visit_n != 5,]
  #   dat1_subset=ddply(dat1_subset, "ID", transform, f_first_hiv_pos=ifelse(RNA==1 & calc_lag(RNA)==0, 1, 0))
  #   dat1_subset=ddply(dat1_subset, "ID", transform, f_last_hiv_neg=ifelse(RNA==0 & calc_lead(RNA)==1, 1, 0))
  #   dat1_subset_v3=subset(dat1_subset, 
  #                         subset=(f_first_hiv_pos==1 | f_last_hiv_neg==1),
  #                         select=c(ID, TIME,  f_first_hiv_pos, f_last_hiv_neg))
  #   dat1_melt=subset(melt(dat1_subset_v3, id=c("ID", "TIME")), value==1, select=c(ID, TIME, variable))
  #   dat1_wide=reshape(dat1_melt, idvar="ID", timevar="variable", direction="wide")
  #   dat1_wide$t_star=ceiling(dat1_wide$TIME.f_last_hiv_neg+(dat1_wide$TIME.f_first_hiv_pos-dat1_wide$TIME.f_last_hiv_neg)/2)
  #   
  #   #exclude those who missed the infusion visit preceding their first RNA+ visit
  # #   dat1_subset=ddply(dat1_subset, "ID", transform, lag_visit_n=calc_lag(visit_n))
  # #   dat1_subset=ddply(dat1_subset, "ID", transform, last_hiv_neg_infu_num=
  # #                       ifelse(f_first_hiv_pos==1 & lag_visit_n %in% infu_visits, lag_visit_n, 
  # #                              ifelse(f_first_hiv_pos==1 & calc_lag(lag_visit_n) %in% infu_visits, calc_lag(lag_visit_n), NA)))
  # #   dat1_subset=ddply(dat1_subset, "ID", transform, missed_prior_infu=
  # #                       ifelse(f_first_hiv_pos==1 & (visit_n-last_hiv_neg_infu_num) <= 3, 0,
  # #                              ifelse(f_first_hiv_pos != 1, 0, 1)))
  # #   cases_missed_prior_infu=unique(dat1_subset[dat1_subset$missed_prior_infu==1,]$ID)
  # #   dat1_subset_v2=dat1_subset[!dat1_subset$ID %in% cases_missed_prior_infu,]
  # #   dat1_subset_v3=subset(dat1_subset_v2, 
  # #                         subset=(f_first_hiv_pos==1 | f_last_hiv_neg==1),
  # #                         select=c(ID, TIME,  f_first_hiv_pos, f_last_hiv_neg))
  #       
  #   # for each case include only one record: time zero to t_star
  #   dat1_merge=merge(dat1_wide, dat1, by=c("ID"))
  #   dat1_final=subset(dat1_merge, TIME == t_star)
  #   dat1_final$status=ifelse(dat1_final$t_star==dat1_final$TIME, 1, 0)
  #   dat1_final$group=1
  #   dat1_final=dat1_final[,c("ID", "group", "TIME", "status", "IPRE","t_IPRE")]  
  # 
  # }
  
  if(case_input=="true"){
    
    #round t_true up to nearest integer
    
    dat1$t_event_round=ceiling(dat1$t_event_num)
    #dat1_subset=subset(dat1, TIME==t_event_round)
    # for each case include only one record: time zero to t_true
    dat1_final=subset(dat1, TIME == t_event_round)
    #dat1_final=subset(dat1, (is.na(calc_lead(TIME)) | (calc_lead(TIME) > t_event_round & TIME <= t_event_round)))
    dat1_final$status=ifelse(dat1_final$t_event_round==dat1_final$TIME, 1, 0)
    dat1_final$group=1
    dat1_final=dat1_final[,c("ID", "group", "TIME", "status","t_IPRE")]
    
  }
  
  #----- controls -----
  
  #calculate daily grid (exclude TIME=0)
  dat0_final=subset(dat0, TIME != 0)
  dat0_final$status=0
  dat0_final$group=0
  dat0_final=dat0_final[,c("ID", "group", "TIME", "status", "t_IPRE")]

  #calculate grid as all 4-wk post infusion visits for controls
  #this inherently subsets to 4-week visits for which the prior infusion visit was attended
#   dat0_final=subset(dat0, visit_n %in% week4_visits) 
#   dat0_final$status=0
#   dat0_final$group=0
#   dat0_final=dat0_final[,c("ID", "group", "TIME", "status", "IPRE")]
    
  #----- output analysis dataset -----

  data_out=rbind(dat0_final, dat1_final)
  data_out=data_out[order(data_out$ID, data_out$TIME),]
  #data_out$log_IPRE=log(data_out$IPRE)
  data_out=ddply(data_out, "ID", transform, tstart=ifelse(is.na(calc_lag(TIME)), 0, calc_lag(TIME)))
  data_out$tstop=data_out$TIME
  data_out=data_out[,c("ID", "group", "TIME", "status", "tstart", "tstop","t_IPRE")]
  return(data_out)
  
}

coxmodel <- function(dat2, alpha=0.025){
  print("in cox")
  # #define weights 
  # cases=unique(subset(dat2, group==1)$ID)
  # controls=unique(subset(dat2, group==0)$ID)
  # dat2$a=ifelse(dat2$group==1, 1, length(controls)/(2800-length(cases)))
  # dat2$w=1/dat2$a
  # dat2$w_star=ifelse(dat2$group==1, dat2$a*(1-dat2$a), 1-dat2$a)
  # dat2a=dat2[!duplicated(dat2$ID),]
  
  #modified self and prentice model 
  fit=coxph(Surv(tstart, tstop, status)~t_IPRE+cluster(ID), data=dat2)
 #print(head(print(fit), n=10))
     #print(head(dat2, n=10))
  # dfbeta=resid(fit, type='dfbeta', collapse=dat2$ID)
  # newvar=fit$naive.var + t(dfbeta) %*% diag(dat2a$w_star) %*% dfbeta
  # fit$var=newvar
  
  #retrieve results
  beta=coef(summary(fit))[,1]
  robust_se_beta=coef(summary(fit))[,4]
  beta_lo=confint(fit)[1]
  beta_hi=confint(fit)[2]
  test_statistic=coef(summary(fit))[,5]
  
  two_sided_test=ifelse(test_statistic > -qnorm(alpha) | test_statistic < qnorm(alpha), 1, 0) #to show beta != 0
  higher_tail_test=ifelse(test_statistic > -qnorm(alpha), 1, 0) #to show beta < 0, i.e. higher IPRE --> lower risk 
  
  cbind(beta, robust_se_beta, beta_lo, beta_hi, test_statistic, two_sided_test,higher_tail_test)
}

approach3 <- function(dat, alpha=0.025){
  print("in approach")
  
  dat=as.data.frame(dat)
  n_cases_total=length(unique(dat[dat$RNA==1,]$ID))
  
  #dat_est=prep_approach3(dat, case_input="est")
  dat_true=prep_approach3(dat, case_input="true")
  
  n_cases_appr=length(unique(dat_true[dat_true$group==1,]$ID))
  pe=unique(dat$PE)*100
  
  #dat_est_out=coxmodel(dat_est)
  #colnames(dat_est_out)= paste(colnames(dat_est_out), "est", sep="_")
  dat_true_out=coxmodel(dat_true)
  colnames(dat_true_out)= paste(colnames(dat_true_out), "true", sep="_")
  
  cbind(n_cases_total, n_cases_appr, alpha, pe, dat_true_out)
}  



#-----------------------------------------------------------
# Run all approaches
#-----------------------------------------------------------

result_list <- lapply(dsn,function(i){
  
  # Load 1000 datasets from each dsn directory as a list
  temp=list.files(path=paste0(in_dswd, i), pattern="*.csv", full.names=TRUE)
  n=length(temp)
  temp=temp[1:(n)]
  print(paste(n, "files uploaded"))
  # dat_i=lapply(temp
  #              # , mc.cores = 5
  #              , read.csv)
  
  dstemp <- vector("list", n)
  j <- 1
  while(j < n | j == n) {
    dat_i = read.csv(temp[j])
    temp3 <-approach3(dat_i)
    print(temp3)
    dstemp[[j]] = temp3
    j <- j+1
  }
  
  # 
  # # Approach 3 
  # dstemp=lapply(dat_i
  #                 # , mc.cores = 5
  #                 , approach3)
  
  dsresults=ldply(dstemp)
  print(dsresults)
  rm(n, temp)
  
  res <- data.frame(group=i
                    ,power_true=mean(dsresults$higher_tail_test_true))
                   # ,power_est=mean(dsresults$higher_tail_test_est))ins
  
}

)
result <- do.call(rbind,result_list)
folder <- "results"
dir <- file.path("C:\\Users\\zzhang2\\Desktop\\simdat\\",folder)
csv.file <- "results.csv"
if(file.exists(dir)==FALSE){
  dir.create(dir,recursive=T)   
}
write.csv(result
          ,file=file.path(dir,csv.file)
          ,row.names = FALSE)
print("ended")     


