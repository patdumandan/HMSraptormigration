#===============================================================
#Project: Assessing temporal patterns of the migratory raptor assemblage structure 
#Analysis plan: Fit random change point models implemented in a Bayesian framework
#Objective: Fitting the random change point model to a time-series data
#Notes: removed BLVU (1934-1984) and TUVU (1969-1982), WWII(1942-1945) missing years#########
#===============================================================  
#  load required packages
library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)

rstan_options(auto_write = TRUE) #Avoid recompiling
options(mc.cores = parallel::detectCores()) #work on multiple cores, faster computation time
#=============================================================== 
hms=read.csv(file.choose(), h=T)
saveRDS(list(N=length(hms$year),count=hms$count, zeros4=zeros4,Nsp=16, years=hms$years,obs_days=hms$obs_days, spcode=hms$spcode), file="hms_dat.rds")
hms_dat=readRDS("hms_dat.rds") #dataset
zeros4=rep(0,4) #used for correlation matrix
#readRDS(list1_hms)
#str(list1_hms)
#in running Stan models, parameters need to be "real" (continuous) 
#so create a predictor variable centered on the mean year 
#average year=1982.671###
#stdev:21.20813#####

hms$years = (hms$year - mean(hms$year))/(2 *sd(hms$year)) #standardize/ center around mean year 
#range(hms$years) #lower=-0.86, upper=0.83, to be used as bounds for bkpoint years

hms$obs_days = log(hms$obs.hours/24) #transform obs.effort to have it around the same values as other predictors(year)

Nsp=16 #for ran_bkpoint (no.of species)


#===============================================================
#Model description:

#neg_binomial_2 distribution for count
#Stancode adapted from Brilleman et al. (2017)
#random effects correlated using unstructured variance-covariance matrix (Cholesky factor)
#priors: Cauchy, same as original breakpoint model

mod_c=stan (model_code= " 
            
            data {
            
            vector<lower=0,upper=0>[4] zeros4; // mean vector for random effects distribution
            int N; //no.of observations (length of hms$year)
            real years[N];//explanatory variable (centered around mean)
            int count[N];//no.of individuals per year
            real obs_days[N]; //explanatory variable (centered around mean)
            int Nsp; //no.of populations
            int spcode[N];// id for each population
            
            }
            
            parameters {
            
            vector[3] beta;            // fixed effects, intercept and slopes
            real<lower=-0.86,upper=0.83> bkpoint;// fixed effects, knotpoint (bounding specified to help convergence)
            vector<lower=0>[4] sigma_sp;   // level 2 error sd (sds of the random effects u[j])
            real<lower=0> phi;        // level 1 error sd
            vector[4] u[Nsp];         // random effects (level 2) errors
            cholesky_factor_corr[4] L_u_Corr; // cholesky factor for the random effects correlation matrix
            }
            
            transformed parameters {  
            vector[4] alpha[Nsp];     // random effects
            real y_mu[N];              // mean parameter based on regression equation
            real slope_difference; //difference between slope before and after
            
            
            //calculate slope_difference
            
            
            slope_difference= beta[3]-beta[2];
            
            
            //==========================
            // calculate random effects
            //==========================
            
            for (i in 1:Nsp) {
            for (k in 1:3) alpha[i,k] = beta[k] + u[i,k];
            alpha[i,4] = bkpoint + u[i,4];
            }
            
            //=====================
            // regression equation
            //=====================
            
            for (j in 1:N){
            
            if (years[j]< alpha[spcode[j],4]) {
            
            y_mu[j]= exp(alpha[spcode[j],1]+alpha[spcode[j],2]* (years[j]-alpha[spcode[j],4])+obs_days[j]);
            }
            
            else {
            y_mu[j]= exp(alpha[spcode[j],1]+alpha[spcode[j],3]* (years[j]-alpha[spcode[j],4])+obs_days[j]);
            }
            
            } 
            }
            
            model {
            
            //========
            // priors
            //========
            
            beta[1] ~ cauchy(0, 2.5); // prior: fixed effect, intercept
            beta[2] ~ cauchy(0, 2.5);   // prior: fixed effect, slope before knot
            beta[3] ~ cauchy(0, 2.5);   // prior: fixed effect, slope after knot
            bkpoint ~ uniform(0,2.5);
            // prior: fixed effect, knot point
            
            sigma_sp[1] ~ cauchy(0,5);    // prior: random effect sd, intercept
            sigma_sp[2] ~ cauchy(0,5);    // prior: random effect sd, slope before knot
            sigma_sp[3] ~ cauchy(0,5);    // prior: random effect sd, slope after knot
            sigma_sp[4] ~ cauchy(0,5);    // prior: random effect sd, knot point
            
            phi ~ cauchy(0,5);       // prior: level 1 error sd
            
            L_u_Corr ~ lkj_corr_cholesky(1);
            // prior: cholesky factor for random effects correlation matrix
            // NB. this prior is the lkj correlation distribution with shape parameter 1 
            // which is equivalent to a uniform distribution over the possible correlation 
            // matrices (where a shape parameter > 1 would have resulted in an upside down
            // U-shaped distribution with the mode being located at the identity matrix)
            
            //=============================
            // random effects distribution
            //=============================
            
            for (i in 1:Nsp) u[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(sigma_sp, L_u_Corr));
            // NB. the second parameter here is the cholesky factor L 
            // (for the correlation matrix). It only uses the sd rather 
            // than the variances since Sigma = L*L'
            
            //==================
            // model likelihood
            //==================
            
            count ~ neg_binomial_2(y_mu, phi); // likelihood for the observed data
            
            }
            
            generated quantities {
            
            real y_mu_pred[N];   // predicted mean
            corr_matrix[4] u_Corr;   // random effects correlation matrix
            matrix[4,4] u_Sigma;     // random effects covariance matrix
            vector[4] alpha_tosave[Nsp];
            // monitor random effects for a subset of patients only
            // (for plotting predictions) and do not monitor 'alpha' 
            // in the model above (since it consumes too much memory!)
            
            //==================================================
            // predicted mean outcome using regression equation
            //==================================================
            
            for (i in 1:Nsp) {  
            alpha_tosave[i] = alpha[i];
            }
            
            y_mu_pred = neg_binomial_2_rng(y_mu, phi);
            
            
            
            
            //=====================================================
            // recover the correlation and covariance matrices
            // using the cholesky factor of the correlation matrix
            //=====================================================
            
            u_Corr = multiply_lower_tri_self_transpose(L_u_Corr);    
            // correlation matrix: u_Corr = L_u_Corr * L_u_Corr'
            
            u_Sigma = quad_form_diag(u_Corr, sigma_sp);
            // covariance matrix: u_Sigma = diag(sigma_sp) * u_Corr * diag(sigma_sp)
            
            }", data=(list(N=length(hms$year),count=hms$count, zeros4=zeros4,Nsp=16, years=hms$years,obs_days=hms$obs_days, spcode=hms$spcode)),chain=4, warmup=4000, iter=5500, control=list(max_treedepth=15, adapt_delta=0.99))



#MODEL ASSESSMENT:
shinystan::launch_shinystan(mod_newcode_corr_prior_finali2)
print(mod_newcode_corr_prior_finali2,pars=c("beta[1]", "beta[2]", "beta[3]", "bkpoint", "slope_difference"))
print(mod_newcode_corr_prior_finali2,pars=c("alpha_tosave"))

posterior=as.data.frame(mod_newcode_corr_prior_finali2)

#plot coefficients:
require(bayesplot)
require(ggplot2)
gen1=mcmc_intervals(posterior, pars=c("beta[1]","beta[2]","beta[3]", "bkpoint"), prob=.8, prob_outer=0.95,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("intercept","pre-change slope","post-change slope", "breakpoint"))+vline_0(size = 0.25, color = "darkgray", linetype = 2)+
  labs(x="coefficient estimate")+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
gen1

#per random effects:
#intercept
int1=mcmc_intervals(posterior, pars=c("alpha[3,1]","alpha[16,1]","alpha[10,1]", "alpha[2,1]", "alpha[9,1]","alpha[15,1]","alpha[5,1]", "alpha[8,1]","alpha[13,1]","alpha[4,1]","alpha[14,1]", "alpha[12,1]","alpha[6,1]","alpha[1,1]","alpha[7,1]", "alpha[11,1]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  labs(x="intercept estimates")+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
int1
#slope1
slo1=mcmc_intervals(posterior, pars=c("alpha[16,2]","alpha[2,2]","alpha[3,2]", "alpha[11,2]", "alpha[6,2]","alpha[5,2]","alpha[7,2]", "alpha[13,2]","alpha[4,2]","alpha[9,2]","alpha[10,2]", "alpha[14,2]","alpha[15,2]","alpha[12,2]","alpha[1,2]", "alpha[8,2]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Turkey Vulture", "Bald Eagle", "Black Vulture", "Peregrine Falcon", "Golden Eagle", "Cooper's Hawk", "Merlin", "Red-shouldered Hawk", "Broad-winged Hawk","Northern Harrier", "Osprey", "Red-tailed Hawk", "Sharp-shinned Hawk", "Rough-legged Hawk", "American Kestrel", "Northern Goshawk"))+
  labs(x="pre-change slope estimates")+ggplot2::coord_flip()+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12),axis.text.x=element_text(angle=90))
slo1

#slope2
slo2=mcmc_intervals(posterior, pars=c("alpha[12,3]","alpha[8,3]","alpha[9,3]", "alpha[14,3]", "alpha[15,3]","alpha[1,3]","alpha[4,3]", "alpha[13,3]","alpha[10,3]","alpha[6,3]","alpha[5,3]", "alpha[11,3]","alpha[7,3]","alpha[3,3]","alpha[2,3]", "alpha[16,3]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Rough-legged Hawk","Northern Goshawk", "Northern Harrier","Red-tailed Hawk", "Sharp-shinned Hawk", "American Kestrel", "Broad-winged Hawk", "Red-shouldered Hawk", "Osprey", "Golden Eagle", "Cooper's Hawk","Peregrine Falcon", "Merlin", "Black Vulture", "Bald Eagle", "Turkey Vulture"))+
  labs(x="post-change slope estimates")+ggplot2::coord_flip()+scale_colour_manual(values = c("alpha[12,3]" = "#E08214", "alpha[4,3]" = "#E08214"))+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12), axis.text.x=element_text(angle=90))
slo2

#bkpoint
require(bayesplot)
bkp1=mcmc_intervals(posterior, pars=c("alpha[3,4]","alpha[16,4]","alpha[10,4]", "alpha[2,4]", "alpha[9,4]","alpha[15,4]","alpha[5,4]", "alpha[8,4]","alpha[13,4]","alpha[4,4]","alpha[14,4]", "alpha[12,4]","alpha[6,4]","alpha[1,4]","alpha[7,4]", "alpha[11,4]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  labs(x="breakpoint estimates")+vline_0(size = 0.25, color = "darkgray", linetype = 2)+theme_gray()+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
bkp1



#MODEL DIAGNOSTICS:
shinystan::launch_shinystan(mod_noeff)
print(mod_noeff)
hms=read.csv(file.choose(), h=T)
#check Rhat
#check neff
#check PPD

#PPD plots####

posterior=rstan::extract(mod_newcode_corr_prior_finali2)$y_mu
sp1=posterior[,which(hms$spcode==1)] #matrix of one species
sp2=posterior[,which(hms$spcode==2)] 
matplot(t(sp1),type="l",col="grey",lty=1,xaxt="n", ylab = "count (N)")
axis(1,c(1946:2018),at=c(0:72))
mean_sp2=apply(sp2,2,mean)
lines(mean_sp1~c(1:73),col="white", lwd=2)
amke=subset(hms, species=="AMKE")
points(amke$count, pch=21, cex=0.5, col="black", bg="black")
abline(v=25, col="black",lwd=3,lty=2)



#MODEL OUTPUT VISUALIZATION:

head(posterior)
gen1=mcmc_intervals(posterior, pars=c("beta[1]","beta[2]","beta[3]", "bkpoint"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("intercept","pre-change slope","post-change slope", "change-point"))+
  labs(x="coefficient estimate")+coord_flip()+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12), axis.text.x =element_text(angle=90))
gen1

#per random effects:
#intercept
require(ggplot2)
int1=mcmc_intervals(mod_newcode_corr_prior_finali2, pars=c("alpha[3,1]","alpha[16,1]","alpha[10,1]", "alpha[2,1]", "alpha[9,1]","alpha[15,1]","alpha[5,1]", "alpha[8,1]","alpha[13,1]","alpha[4,1]","alpha[14,1]", "alpha[12,1]","alpha[6,1]","alpha[1,1]","alpha[7,1]", "alpha[11,1]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  labs(x="intercept estimates")+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
int1
#slope1
slo1=mcmc_intervals(posterior, pars=c("alpha[3,2]","alpha[16,2]","alpha[10,2]", "alpha[2,2]", "alpha[9,2]","alpha[15,2]","alpha[5,2]", "alpha[8,2]","alpha[13,2]","alpha[4,2]","alpha[14,2]", "alpha[12,2]","alpha[6,2]","alpha[1,2]","alpha[7,2]", "alpha[11,2]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  labs(x="pre-change slope estimates")+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
slo1

#slope2
slo2=mcmc_intervals(posterior, pars=c("alpha[3,3]","alpha[16,3]","alpha[10,3]", "alpha[2,3]", "alpha[9,3]","alpha[15,3]","alpha[5,3]", "alpha[8,3]","alpha[13,3]","alpha[4,3]","alpha[14,3]", "alpha[12,3]","alpha[6,3]","alpha[1,3]","alpha[7,3]", "alpha[11,3]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  labs(x="post-change slope estimates")+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
slo2

#slo1+2
slo3=mcmc_intervals(posterior, pars=c("alpha[3,3]","alpha[3,2]","alpha[16,3]","alpha[10,3]", "alpha[2,3]", "alpha[9,3]","alpha[15,3]","alpha[5,3]", "alpha[8,3]","alpha[13,3]","alpha[4,3]","alpha[14,3]", "alpha[12,3]","alpha[6,3]","alpha[1,3]","alpha[7,3]", "alpha[11,3]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  labs(x="post-change slope estimates")+
  theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
        axis.line=element_blank(),axis.title=element_text(size=15),axis.text=element_text(size=12))
slo

#bkpoint
bkp1=mcmc_intervals(mod_newcode_corr_prior_finali2, pars=c("alpha[3,4]","alpha[16,4]","alpha[10,4]", "alpha[2,4]", "alpha[9,4]","alpha[15,4]","alpha[5,4]", "alpha[8,4]","alpha[13,4]","alpha[4,4]","alpha[14,4]", "alpha[12,4]","alpha[6,4]","alpha[1,4]","alpha[7,4]", "alpha[11,4]"), prob=.8,
                    point_est = "median")+
  ggplot2::scale_y_discrete(labels = c("Black Vulture", "Turkey Vulture", "Osprey", "Bald Eagle", "Northern Harrier", "Sharp-shinned Hawk", "Cooper's Hawk", "Northern Goshawk", "Red-shouldered Hawk", "Broad-winged Hawk","Red-tailed Hawk", "Rough-legged Hawk", "Golden Eagle", "American Kestrel", "Merlin", "Peregrine Falcon"))+
  scale_x_discrete()+scale_x_discrete(labels=yearlab)
#labs(x="change-point year")+
theme(plot.title=element_text(hjust=0,size=18),panel.border=element_rect(size=1,fill=NA,colour="black",linetype="solid"),
      axis.title=element_text(size=15),axis.text=element_text(size=12))

bkp1

#model interpretation
#1.probability of slope difference in global parameters to be significantly different
length(which(posterior$'beta[3]'>posterior$'beta[2]'))/length(posterior$'beta[3]')
quantile(posterior$'slope_difference',c(0.025,0.5,0.975))

#2.correlation between slope differences and traits
reference=read.csv(file.choose(), h=T)
cor.test(reference$value,reference$body.mass) #0.49
cor.test(reference$value,reference$Pop.Est) #0. 40

#3.contingency table for binary traits (functional)
ref=reference%>%group_by(variable)%>%mutate(n=row_number())
str(ref)
head(ref)

ref$n[ref$variable=="American Kestrel"]==ref$n[ref$variable=="Black Vulture"] #to check if the indexes are the same for all species

theta1=rep(NA, 6000)
theta2=rep(NA, 6000)

for (i in 1:6000) {
  subs=ref[ref$n==i,c("value","human.tolerance")]
  y1=sum(subs$value[subs$human.tolerance=="tolerant"]>0)
  y2=sum(subs$value[subs$human.tolerance=="intolerant"]>0)
  n1=10
  n2=6
  theta1[i]=rbeta(1,y1+1, (n1-y1)+1)
  theta2[i]=rbeta(1,y2+1, (n2-y2)+1)
  
}
str(ref)


diet_prob=length(which(theta1>theta2))/length(theta2) #0.2526
hab_prob=length(which(theta1>theta2))/length(theta2) #0.999
ddt_prob=length(which(theta1>theta2))/length(theta2) #0.893
beh_prob=length(which(theta1>theta2))/length(theta2) #0.9363
anth_prob=length(which(theta1>theta2))/length(theta2) #0.9475

#FIGURES USED IN MANSUCRIPT###############

#Behavior
behaviorplot=reference %>% 
  ggplot(aes(x= fct_reorder(variable,ordered.behavior), y=value, fill=Trend)) +
  geom_boxplot(outlier.size=0.5) +ggtitle("Migratory Behavior")+
  xlab("")+ylab("")+
  annotate("text", x=2.0, y=10.5, label= "complete", size=4)+
  annotate("text", x=10.0, y=10.5, label= "partial", size=4)+
  annotate("rect", xmin = 0, xmax = 3.5, ymin = -Inf, ymax = Inf,fill="grey20",alpha = 0.25, size=2)+
  scale_fill_manual(values=c("increasing"="sky blue", "decreasing"="orange", "no change"="pink"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

behaviorplot


#tolerance
str(reference)
tolerplot=reference %>% 
  ggplot(aes(x= fct_reorder(variable,ordered.tolerance), y=value, fill=Trend)) +
  geom_boxplot(outlier.size=0.5) +ggtitle("Anthropogenic Tolerance")+
  xlab("")+ylab("")+
  annotate("text", x=3.5, y=10.5, label= "intolerant", size=4)+
  annotate("text", x=10.0, y=10.5, label= "tolerant", size=4)+
  annotate("rect", xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf,fill="grey20",alpha = 0.25, size=2)+
  scale_fill_manual(values=c("increasing"="sky blue", "decreasing"="orange", "no change"="pink"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
tolerplot

#DIET
dietplot=reference %>% 
  ggplot(aes(x= fct_reorder(variable,ordered.diet_spec), y=value, fill=Trend)) + #use ordered.diet if sort by func.group
  geom_boxplot(outlier.size=0.5) +
  annotate("text", x=3.5, y=10.75, label= "specialist", size=4)+
  #annotate("text", x=4.5, y=10.75, label= "carnivore", size=4)+
  annotate("text", x=12, y=10.75, label= "generalist", size=4)+
  #annotate("text", x=13.5, y=10.75, label= "piscivore", size=4)+
  #annotate("text", x=16, y=10.75, label= "scavenger", size=4)+
  xlab("")+ylab("")+ggtitle("Diet Specialization")+
  scale_fill_manual(values=c("increasing"="sky blue", "decreasing"="orange", "no change"="pink"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  annotate("rect", xmin = 0, xmax = 8.5, ymin = -Inf, ymax = Inf,fill="grey20",alpha = 0.25, size=2)
dietplot

#DDT
ddtplot=reference %>% 
  ggplot(aes(x= fct_reorder(variable,DDT.order), y=value, fill=Trend)) +
  geom_boxplot(outlier.size=0.5) +ggtitle("DDT Susceptibility")+
  xlab("")+ylab("")+
  annotate("text", x=2.5, y=10.5, label= "non-susceptible", size=4)+
  annotate("text", x=10.0, y=10.5, label= "susceptible", size=4)+
  annotate("rect", xmin = 0, xmax = 3.5, ymin = -Inf, ymax = Inf,fill="grey20",alpha = 0.25, size=2)+
  scale_fill_manual(values=c("increasing"="sky blue", "decreasing"="orange", "no change"="pink"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ddtplot

#habitat########

habitataplot=reference %>% 
  ggplot(aes(x= fct_reorder(variable,ordered.habitat.spec), y=value, fill=Trend)) +
  geom_boxplot(outlier.size=0.5) +ggtitle("Habitat Specialization")+
  xlab("")+ylab("")+
  annotate("text", x=3.5, y=10.5, label= "specialist", size=4)+
  annotate("text", x=10.0, y=10.5, label= "generalist", size=4)+
  annotate("rect", xmin = 0, xmax = 6.5, ymin = -Inf, ymax = Inf,fill="grey20",alpha = 0.25, size=2)+
  scale_fill_manual(values=c("increasing"="sky blue", "decreasing"="orange", "no change"="pink"))+
  theme_bw()+theme(axis.text.x=element_text(angle=90, vjust=0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank())
