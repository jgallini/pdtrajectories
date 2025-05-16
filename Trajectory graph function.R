
#This is a function that will create fitted graphs for trajectories
#including regression intervals for linear mixed effects models

#Plots will display different lines for the two groups of the chosen binary
#split variable (only designed for split variable to be binary 0/1)


traj_graph <- function (dat, baseline_dat, split,
                        title="Insert figure title",outcome,
                        seed=1046898,graph_groups=2) {
  

#libraries
library(lme4)
library(lmerTest)
library(table1)
library(dplyr)
library(tidyr)
library(merTools)
library(ggplot2)
library(rlang)

split_var<-sym(split)
outcome_var<-sym(outcome)
outcome_name<-ifelse(outcome=="PHC_MEM","Memory",
                     ifelse(outcome=="PHC_EXF","Executive Function",
                            ifelse(outcome=="PHC_LAN","Language","Unknown Outcome")))


#main mixed effects model using real data
#this model uses continuous versions of binary variables so that we can 
#fit the model using just proportions in the future when we're estimating
#average model trajectory (because R can't do math with factors)

formula_str <- paste0(outcome_var, " ~ time_fmbsl + NACCAGEB + RACE_cont + EDUC+
SEX_cont + Post_MCI + HYPERTEN_cont + DIABETES_cont + NACCTBI_cont + 
SMOKYRS + BMI_bin + PDb4CI + NACCNE4s_cont +
SEX_cont*time_fmbsl + NACCNE4s_cont*time_fmbsl + PDb4CI*time_fmbsl +
time_fmbsl*Post_MCI + SEX_cont*Post_MCI + NACCNE4s_cont*Post_MCI+ PDb4CI*Post_MCI +
SEX_cont*time_fmbsl*Post_MCI + NACCNE4s_cont*time_fmbsl*Post_MCI+
PDb4CI*time_fmbsl*Post_MCI + (time_fmbsl|NACCID)")

lmer_model <- lmer(as.formula(formula_str), data = dat)
fixef<-fixef(lmer_model)

#creating medians and proportions of variables to use in average model fit
#trajectory graphs
time_points<-6

#variables that we're adjusting for- checking to see if graph split by
#the variable first, if yes using all levels for that variable, if not
#using mean/median
#split option only for binary variables right now
NACCAGEB<-rep(median(baseline_dat$NACCAGEB),time_points*graph_groups)
time_fmbsl<-rep(seq(1:median(table(dat$NACCID))),graph_groups)
time_fmbsl<-time_fmbsl-1
RACE_cont<-rep(sum(baseline_dat$RACE_bin=="Non-White",na.rm=TRUE)/nrow(baseline_dat),
               time_points*graph_groups)
EDUC<-rep(median(baseline_dat$EDUC,na.rm=TRUE),12)
SMOKYRS<-rep(median(baseline_dat$SMOKYRS,na.rm = TRUE),time_points*graph_groups)

if (split == "SEX_cont") {
  SEX_cont <- c(rep(0, 6), rep(1, 6))
} else {
  SEX_cont <- rep(
    sum(baseline_dat$SEX == "Female", na.rm = TRUE) / nrow(baseline_dat),
    time_points * graph_groups
  )
}

HYPERTEN_cont<-rep(sum(baseline_dat$HYPERTEN=="Yes",na.rm=TRUE)/nrow(baseline_dat),
                   time_points*graph_groups)
DIABETES_cont<-rep(sum(baseline_dat$DIABETES=="Yes",na.rm=TRUE)/nrow(baseline_dat),
                   time_points*graph_groups)
NACCTBI_cont<-rep(sum(baseline_dat$NACCTBI=="Yes",na.rm=TRUE)/nrow(baseline_dat),
                  time_points*graph_groups)
BMI_bin<-rep(sum(baseline_dat$BMI_bin=="Obese",na.rm=TRUE)/nrow(baseline_dat),
             time_points*graph_groups)
if (split == "PDb4CI") {
  PDb4CI <- c(rep(0, 6), rep(1, 6))
} else {
  PDb4CI <- rep(
    sum(baseline_dat$PDb4CI == "Yes", na.rm = TRUE) / nrow(baseline_dat),
    time_points * graph_groups
  )
}
if (split == "NACCNE4s_cont") {
  NACCNE4s_cont <- c(rep(0, 6), rep(1, 6))
} else {
  NACCNE4s_cont <- rep(
    sum(baseline_dat$NACCNE4S_bin == "1 or 2 alleles", na.rm = TRUE) / nrow(baseline_dat),
    time_points * graph_groups
  )
}
#now creating interaction terms so we can solve for fitted mean
SEX_cont_time_fmbsl<-SEX_cont*time_fmbsl
NACCNE4s_cont_time_fmbsl<-NACCNE4s_cont*time_fmbsl
PDb4CI_time_fmbsl<-PDb4CI*time_fmbsl
time_fmbsl_Post_MCI<-time_fmbsl*Post_MCI
SEX_cont_Post_MCI<-SEX_cont*Post_MCI
NACCNE4s_cont_Post_MCI<-NACCNE4s_cont*Post_MCI
PDb4CI_Post_MCI<-PDb4CI*Post_MCI
SEX_cont_time_fmbsl_Post_MCI<-SEX_cont*time_fmbsl*Post_MCI
NACCNE4s_cont_time_fmbsl_Post_MCI<-NACCNE4s_cont*time_fmbsl*Post_MCI
PDb4CI_time_fmbsl_Post_MCI<-PDb4CI*time_fmbsl*Post_MCI

#Getting Post_MCI conversion median visit number to include spline
clean_dat2<-dat %>% group_by(NACCID) %>% 
  mutate(first_post=Post_MCI-lag(Post_MCI)==1)
med_switch<-median(clean_dat2$NACCVNUM[clean_dat2$first_post==TRUE],na.rm=TRUE)
Post_MCI<-rep(c(rep(0,med_switch-1),rep(1,median(table(dat$NACCID))-
                                          (med_switch-1))),graph_groups)

#creating data frame for fitted/predicted modeling
preddat<-as.data.frame(cbind(time_fmbsl,NACCAGEB,RACE_cont,EDUC,SEX_cont,Post_MCI,
               HYPERTEN_cont,DIABETES_cont,NACCTBI_cont,SMOKYRS,BMI_bin,PDb4CI,
               NACCNE4s_cont,SEX_cont_time_fmbsl, NACCNE4s_cont_time_fmbsl,
               PDb4CI_time_fmbsl,time_fmbsl_Post_MCI, SEX_cont_Post_MCI,
               NACCNE4s_cont_Post_MCI,PDb4CI_Post_MCI,
               SEX_cont_time_fmbsl_Post_MCI,NACCNE4s_cont_time_fmbsl_Post_MCI,
               PDb4CI_time_fmbsl_Post_MCI))

#solving for exact fitted estimates using fixed effects from this model
weighted_ests <- sweep(preddat, 2, fixef[-1], `*`)
preddat$exact_fit<-rowSums(weighted_ests)+fixef[1]

#adding NACCID as if this is one average person, this will ignore random effects
preddat$NACCID<-rep("NACC999",times=time_points*graph_groups)

set.seed(seed) #predictinterval uses parametric bootstrapping for estimates

#predicting new values of outcome
preds_model <- suppressWarnings(
  predictInterval(lmer_model, newdata = preddat, 
                  level = 0.95,  # 95% CI
                  n.sims = 1000,  # Number of simulations
                  stat = "mean",
                  type = "linear.prediction",
                  include.resid.var = TRUE,
                  returnSims=TRUE))
#getting this information out for bootstrapped variance

#These fitted values will be used to plot the actual lines
#graph_dat<-cbind(preddat,preds_model)

#now working on getting data to plot confidence intervals on graphs
#calculating bootstrapped standard deviation and merging back with data
boot_draws <- t(attr(preds_model, "sim.results"))
sd<-vector()
for (i in 1:(graph_groups*time_points)){
  sd[i]<-sd(boot_draws[,i])
}
graph_dat2<-cbind(preddat,sd)

#randomly drawing samples from appropriate distributions for
#each time point using actual empirical sample size at each time point
sample_sizes<-model.frame(lmer_model) %>% group_by(time_fmbsl,!!split_var) %>% 
  summarise(n=n(),.groups="drop")
sample_sizes$time_disc<-cut(sample_sizes$time_fmbsl, 
                            breaks=c(0,seq(1.5,time_points+0.5)))
sample_sizes<-sample_sizes %>% group_by(time_disc,!!split_var) %>% 
  summarise(n_disc=n(),.groups="drop")
sample_sizes <- sample_sizes[!is.na(sample_sizes$time_disc),] %>% 
  arrange(!!split_var)

draws<-matrix(NA,nrow=graph_groups*time_points, ncol=max(sample_sizes$n_disc))
for (i in 1:(graph_groups*time_points)){
  n <- sample_sizes$n_disc[i]
  vals <- rnorm(n, mean = graph_dat2$exact_fit[i], sd = graph_dat2$sd[i])
  draws[i, 1:n] <- vals 
}
graph_dat3<-as.data.frame(cbind(graph_dat2$time_fmbsl,graph_dat2$Post_MCI, 
                                graph_dat2[[split_var]],draws))

#transposing to graphing shape
start_col<-which(names(graph_dat3)=="V4")
graph_dat4 <- graph_dat3 %>%
  pivot_longer(
    cols = all_of("V4"):last_col(), 
    names_to = "Draw Number",         
    values_to = "fit"       
  )
colnames(graph_dat4)<-c("time_fmbsl","Post_MCI",split,"Draw Number","fit")

#subsetting data to only time points of interest for graphing 
#so we get straight lines for each spline
graph_dat5<-graph_dat4[graph_dat4$time_fmbsl %in% c(0,med_switch-2,time_points-1),]
spline1<-graph_dat5[graph_dat5$time_fmbsl<=(med_switch-2),]
spline2<-graph_dat5[graph_dat5$time_fmbsl>=(med_switch-2),]

#creating custom legend labels for each split by variable
legend_labs<-vector()
if (split == "PDb4CI") {
  legend_labs <- c("MCI before PD diagnosis",
                   "PD diagnosis before MCI")
} else if (split == "SEX_cont") {
  legend_labs <- c("Male", "Female")
} else {
  legend_labs <- c("0 e4 alleles", "1 or 2 e4 alleles")
}

  
#creating plot
#need to store actual mean from fitted models and use in ggplot portion,
#then do draws only for geom_smooth portion, still need to work on this
plot<-suppressWarnings(ggplot(graph_dat5,aes(x=time_fmbsl, y=exact_fit, color=as.factor(!!split_var),
              group=as.factor(!!split_var))) +
  #geom_point(data=graph_dat2,aes(x=time_fmbsl,y=exact_fit,color=as.factor(!!split_var)))+
  geom_smooth(data=spline1,aes(x=time_fmbsl,y=fit,color=as.factor(!!split_var)),
              formula = y ~ x,method=lm)+  
  geom_smooth(data=spline2,aes(x=time_fmbsl,y=fit,color=as.factor(!!split_var)),
              formula = y ~ x,method=lm)+
  ylab(paste0("Standardized ",outcome_name , " Factor Score"))+
  xlab("Time (years)")+ 
  ggtitle(title)+
  theme(legend.position = c(0.2, 0.3))+
  labs(color = NULL)+
  scale_color_manual(values=c("cornflowerblue","salmon"),
                     labels = legend_labs) +
  annotate("segment", x = 2, xend = 2, 
           y = quantile(graph_dat5$fit[graph_dat5$time_fmbsl==2],0.25,na.rm=TRUE)-0.2, 
           yend = quantile(graph_dat5$fit[graph_dat5$time_fmbsl==2],0.25,na.rm=TRUE), 
           arrow = arrow(type = "open", length = unit(0.2, "inches")),  
           color = "black", linewidth = 1) +
  annotate("text", x = 2, 
           y = quantile(graph_dat5$fit[graph_dat5$time_fmbsl==2],0.25,na.rm=TRUE)-0.24, 
           label = "MCI", color = "black",size = 5, hjust = 0.5))

return(plot)

}

