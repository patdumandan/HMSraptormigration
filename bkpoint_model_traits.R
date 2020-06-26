#===============================================================
#Project: Assessing temporal patterns of the migratory raptor assemblage structure 
#Analysis plan: Fit random change point model implemented in a Bayesian framework
#with slopes varying per functional trait
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

hms$diets = (hms$diet - mean(hms$diet))/(2 *sd(hms$diet)) 
hms$gentime = (hms$generation.time - mean(hms$generation.time))/(2 *sd(hms$generation.time)) 
hms$size = (hms$mass.kg - mean(hms$mass.kg))/(2 *sd(hms$mass.kg)) 
hms$habitat = (hms$habitat.types - mean(hms$habitat.types))/(2 *sd(hms$habitat.types)) 
hms$years = (hms$year - mean(hms$year))/(2 *sd(hms$year)) #standardize/ center around mean year 
hms$obs_days = log(hms$obs.hours/24) #transform obs.effort to have it around the same values as other predictors(year)
hms$diet=hms$diet.index
Nsp=16 #for ran_bkpoint (no.of species)
hms$migcode=as.integer(as.factor(hms$migration))
mod_d=stan (model_code= " 
            
            data {
            
            vector<lower=0,upper=0>[4] zeros4; // vector for random effects distribution
            int N; //no.of observations (length of hms$year)
            real years[N];//explanatory variable (centered around mean)
            int count[N];//no.of individuals per year
            real diet[N];
            real mass[N];
            real habitat[N];
            real gentime [N];
           // int migration [N];
            real obs_days[N]; //
            int Nsp; //no.of populations
            int spcode[N];// id for each population
            
            }
            
            parameters {
            
            vector[3] beta;            // fixed effects, intercept and slopes
            real<lower=-0.86,upper=0.83> bkpoint;// fixed effects, knotpoint (bounding specified to help convergence)
            vector<lower=0>[4] sigma_sp;   // level 2 error sd (sds of the random effects u[j])
            real<lower=0> phi;        // level 1 error sd
            vector[4] u[Nsp];         // 
            cholesky_factor_corr[4] L_u_Corr; // cholesky factor for the random effects correlation matrix
            real diet_effb;
            real mass_effb;
            real habitat_effb;
            real diet_effa;
            real mass_effa;
            real habitat_effa;
            real gentime_effa;
            real gentime_effb;
         //real migration_effa;
           // real migration_effb;
           
            }
            
            transformed parameters {  
            vector[4] alpha[Nsp];     // random effects
            real y_mu[N];              // mean parameter based on regression equation
            real slope_diff[N];
            
            for (i in 1:Nsp) {
            alpha[i,1] = beta[1] + u[i,1];
            alpha[i,2] = beta[2] +u[i,2];
            alpha[i,3] = beta[3] +u[i,3];
            alpha[i,4] = bkpoint + u[i,4];
            
            
            }
            for (k in 1:N) {
            slope_diff[k]= beta[3]-beta[2];};
            
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
            bkpoint ~ cauchy(0,2.5);
          diet_effb~ normal(0,1);
            mass_effb~ normal(0,1);
            habitat_effb~ normal(0,1);
           diet_effa~ normal(0,1);
            mass_effa~ normal(0,1);
            habitat_effa~ normal(0,0.1);
            gentime_effa~ normal(0,1);
            gentime_effb~ normal(0,1);
            
            for(i in 1:Nsp) 
            
            {
            slope_diff[i]~ normal((exp(alpha[i,1]+mass[i]*mass_effa)+(diet[i]*diet_effa)+(habitat[i]*habitat_effa)+(gentime[i]*gentime_effa)),1);
            };
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
            
            }", data=(list(N=length(hms$year),gentime=hms$gentime, habitat= hms$habitat ,count=hms$count,mass=hms$size, diet=hms$diets, zeros4=zeros4,Nsp=16, years=hms$years,obs_days=hms$obs_days, spcode=hms$spcode)),chain=2, iter=200)

print(mod_d,pars=c("beta[1]", "beta[2]", "beta[3]", "bkpoint", "diet_effa","diet_effb", "mass_effa","mass_effb", "habitat_effa","habitat_effb", "gentime_effa", "gentime_effb"))
print(mod_d, pars=c("slope_diff"))
print(mod_d,pars=c("alpha_tosave"))
posterior=as.data.frame(mod_c)
y_pred=rstan::extract(mod_d)
y_pred1=rstan::extract(mod_d)$y_mu
mean_y_mu1=apply(y_pred1,2, mean)
saveRDS(mod_d, "model_bkpttraits_v1.RDS")
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