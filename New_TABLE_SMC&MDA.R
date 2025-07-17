setwd("C:/Users/Administrator/Documents/Nouveau Projet R")
###Psi MDA
TheMDAStrategy<-function(t_begin_Camp,MDA_Number_of_cycle,
                         MDA_VectTime_between_cycles, MDA_Dur_cycle,
                         Gap,time,Nt){
  bar_psi=rep(0,Nt)
  if (MDA_Number_of_cycle >= 1) {
    for (t in 1:Nt) {
      for (k in 0:(MDA_Number_of_cycle - 1)) { 
        if (k == 0) { ta = t_begin_Camp}
        else {ta = t_begin_Camp + sum(MDA_Dur_cycle + MDA_VectTime_between_cycles[1:k]) }
        tb = ta + MDA_Dur_cycle - 1
        if (ta - Gap <= time[t] & time[t] < ta) { bar_psi[t] = (time[t]-ta+Gap)/Gap } 
        else if (ta <= time[t] & time[t] <= tb) { bar_psi[t] = 1 } 
        else if (tb < time[t] & time[t] <= tb + Gap){bar_psi[t] = (tb + Gap - time[t])/Gap}
      }
    }
  }
  
  FunList = list("bar_psiMDA" = bar_psi)
  return(FunList)
}
##phi SMC
TheSMCStrategy<-function(t_begin_Camp,SMC_Number_of_cycle,
                         SMC_VectTime_between_cycles,SMC_Dur_cycle,Gap,
                         time,Nt){
  bar_phi=rep(0,Nt)
  if (SMC_Number_of_cycle >= 1) {
    for (t in 1:Nt) {
      for (k in 0:(SMC_Number_of_cycle - 1)) { 
        if (k == 0) { ta = t_begin_Camp}
        else {ta = t_begin_Camp + sum(SMC_Dur_cycle + SMC_VectTime_between_cycles[1:k]) }
        tb = ta + SMC_Dur_cycle - 1
        if (ta - Gap <= time[t] & time[t] < ta) { bar_phi[t] = (time[t]-ta+Gap)/Gap } 
        else if (ta <= time[t] & time[t] <= tb) { bar_phi[t] = 1 } 
        else if (tb < time[t] & time[t] <= tb + Gap){bar_phi[t] = (tb + Gap - time[t])/Gap}
      }
    }
  }
  
  FunList = list("bar_phiSMC" = bar_phi)
  return(FunList)
}

# Efficacité du control SMC ou MDA
eta_fun <- function(sigma,T50,shape) {
  eta= 1 - 1 / (1 + (sigma / T50)^shape)
  return(eta)
}

# Taux de de perte Efficacité
kappa_fun <- function(sigma,T50,shape) {
  kappa= shape / T50 * ( (sigma / T50)^(shape - 1) ) / (1 + (sigma / T50)^shape)
  return(kappa)
}

Mainfunction = function(MDA_Dur_cycle, MDA_VectTime_between_cycles,
                        MDA_Number_of_cycle,t_begin_Camp,
                        SMC_Number_of_cycle,SMC_VectTime_between_cycles,
                        SMC_Dur_cycle,PropSMC, PropMDA, Gap,Nyr,LengYr, 
                        PropTreatment,Tmax,Lm, dt, p_f, a_lim,
                        initial_conditions = NULL){
  
  ##Variables 
  {
    Tmax = Nyr * LengYr;
    time = seq(0, Tmax, by = dt); Nt = length(time);
    da = 1; amax = 90; Va = seq(0, amax, by = dt); Na = length(Va);
    dsigmaSmc = 1; sigmaSmc_max=120#40;
    VsigmaSmc=seq(0,sigmaSmc_max, by = dsigmaSmc);
    NsigmaSmc=length(VsigmaSmc);
    dsigmaMda = 1; sigmaMda_max=120#60; 
    VsigmaMda=seq(0,sigmaMda_max, by = dsigmaMda);
    NsigmaMda=length(VsigmaMda);
    dsigmaR = 1; sigmaR_max=40; VsigmaR=seq(0,sigmaR_max, by = dsigmaR);
    NsigmaR=length(VsigmaR)
  } 
  
  ##Parametres fix
  {
    epsi_s=1/10; bar_deltam=0.55; bar_epss=1/10; bar_bm = 0.5; Lh = 30750;
    mum = 0.13; betam_s = 0.25 #0.25; #mum = 0.13; ref Hannah slater ou FranckDongmoEtAl
    deltam_s = 1/20; gammah_su_base = 0; gammah_ru_base = 0; deltah_st_base = 0;
    deltah_ru_base = 9 * 10^(-5); rhoh_base = 5.5 * 10^(-4)
    deltam_r = bar_deltam * deltam_s; betam_r = bar_bm * betam_s; PropTreatRate = 1; 
    betam = diag(c(betam_s, betam_r)); deltam = diag(c(deltam_s, deltam_r));
    epsi_r = bar_epss * epsi_s; InitPrevR = 0.001; 
  } 
  
  ## Initialization of parameters per age
  {
    deltah_st = rep(deltah_st_base, Na); deltah_ru = rep(deltah_ru_base, Na);
    gammah_su = rep(gammah_su_base, Na); rhoh = rep(rhoh_base, NsigmaR);
    gammah_ru = rep(gammah_ru_base, Na); 
  }    
  
  ## Définition des fonctions pour les paramètres dépendant de l'âge
  {
    ##recovery rate gammah(a)
    bar_gamma = 0.0;
    gammah_st = rep(3.704*10^(-2), Na); gammah_rt = bar_gamma * gammah_st
    
    ##betah_s(a) & betah_r(a) 
    alpha_1 = 0.04  #0.071[0.023;0.175]   
    alpha_2 = 0.18   #0.302[0.16; 0.475]
    betah_s = rep(NA, Na);betah_r = rep(NA, Na);
    
    # Define the function G(a)
    G = function(a) {
      return(22.7 * a * exp(-0.0934 * a))
    }
    # Compute beta_h for each age in Va
    hbr=0.5
    for (i in 1:Na) {
      betah_s[i] = hbr*alpha_1 * (G(Va[i]) ^ alpha_2)  
      betah_r[i] = betah_s[i]
    }
    
    ##nu_h(a) progression Asympto vers cas clinique
    {
      #p_f = 0.15 Fertility probability for women aged 15 to 40 in Bobo Dioulasso,
      nu_0=1/10; nu_1=1/150; nu_2=1/8;
      nuh_m=rep(0,Na) # disease mortality for men
      nuh_f = rep(0, Na) # disease mortality for female
      for (i in 1:Na) {
        if (Va[i] <= 5) { nuh_m[i] = nu_2; nuh_f[i] = nu_2 } 
        else if (Va[i] > 5 && Va[i] <= 15) {nuh_m[i] = nu_0; nuh_f[i] = nu_0 }
        else {nuh_m[i] = nu_1 }
        
        if (Va[i] > 15 && Va[i] <= 45) { nuh_f[i] = p_f * nu_0 + (1 - p_f) * nu_1 } 
        else if (Va[i] > 45) { nuh_f[i] = nu_1}
      }
      
      # plot(Va,nuh_m)
      # plot(Va,nuh_f)
    }
    
    ##mu_h(a) mortalité naturelle 
    muh_values = c(66.8, 7.6, 1.7, 0.9, 1.3, 1.9, 2.4, 2.8, 3.6, 4.7, 6.3, 8.9,
                   13.2, 19.8, 31.1, 47.7, 71.3, 110.5, 186.7) / 1000 
    age_labels0=c(0,1,seq(5,90,by=5))
    
    muh = rep(0, Na) 
    # Assign mu_h according to age_labels
    for (i in 1:Na) {
      for (IdAgeGroup in 1:(length(age_labels0)-1)) {
        if(age_labels0[IdAgeGroup]<=Va[i] & Va[i]<age_labels0[IdAgeGroup+1]) {
          muh[i]=muh_values[IdAgeGroup]
        }
      }
      if (i>=Na) {muh[i]=muh_values[length(muh_values)]}
    }
    
    ##delta_h(a), Mortalité induite maladie
    deltah_values = c(1.07 * 10^(-3), 7.02 * 10^(-4), 4.55 * 10^(-4),
                      5.73 * 10^(-5))
    bar_delta=0.5
    deltah_su = rep(0, Na); deltah_rt= rep(0, Na);
    
    for (i in 1:Na) {
      if(0<=Va[i] & Va[i]<=1) {
        deltah_su[i]=deltah_values[1]
        deltah_rt[i] = bar_delta * deltah_su[i]
      } else if(1<Va[i] & Va[i]<=5) {
        deltah_su[i]=deltah_values[2]
        deltah_rt[i] = bar_delta* deltah_su[i]
      } else if(5<Va[i] & Va[i]<=15) {
        deltah_su[i]=deltah_values[3]
        deltah_rt[i] = bar_delta * deltah_su[i]
      } else if(15<Va[i]) {
        deltah_su[i]=deltah_values[4]
        deltah_rt[i] = bar_delta * deltah_su[i];
      }
    }
    ##eta_M, k_M  
    eta_M = rep(0, NsigmaMda); k_M = rep(0, NsigmaMda)
    
    eta_M= eta_fun(VsigmaMda, T50=Eff_MDA_50, shape=shape_MDA)
    k_M = kappa_fun(VsigmaMda, T50=Eff_MDA_50, shape=shape_MDA)
    
    ##eta_S, k_S  
    eta_S =  rep(0,   NsigmaSmc); k_S = rep(0,   NsigmaSmc) 
    
    eta_S= eta_fun(VsigmaSmc, T50=Eff_SMC_50, shape=shape_SMC)
    k_S = kappa_fun(VsigmaSmc, T50=Eff_SMC_50, shape=shape_SMC)
    
    # Initialisation
    barbetahS_s = rep(0, Na)
    barbetahS_r = rep(0, Na)
    barbetahM_s = rep(0, Na)
    barbetahM_r = rep(0, Na)
    
    # Calcul du taux d’infection modifié par âge et par SMC et MDA
    for (i in 1:Na) {
      barbetahS_s[i] = 0.5*betah_s[i]
      barbetahS_r[i] = 0.5*betah_r[i]
      barbetahM_r[i]= 0.5*betah_r[i]
      barbetahM_s[i]= 0.5*betah_s[i]
    }
  }
  
  #######INITIALISATION DES MATRICES
  {
    betah = array(NA, dim = c(2, 2, Na))
    deltah = array(NA, dim = c(4, 4, Na))
    gammah = array(NA, dim = c(4, 4, Na))
    barbetahS = array(NA, dim = c(2, 2, Na))
    barbetahM = array(NA, dim = c(2, 2, Na))
    
    ######Matrix parameters with age
    for(a in 1:Na) {
      betah[,,a] = diag(c(betah_s[a], betah_r[a]))
      deltah[,,a] = diag(c(deltah_st[a], deltah_rt[a], deltah_su[a], deltah_ru[a]))
      gammah[,,a] = diag(c(gammah_st[a], gammah_rt[a], gammah_su[a], gammah_ru[a]))
      barbetahS[,,a] = diag(c(barbetahS_s[a], barbetahS_r[a])) 
      barbetahM[,,a] = diag(c(barbetahM_s[a], barbetahM_r[a]))
    } 
  }  
  
  epsi=matrix(c(1-epsi_s, epsi_r, 1, 0, epsi_s,1-epsi_r,0,1),ncol = 4,byrow = TRUE)
  f = matrix(1, ncol = 4, byrow = TRUE); m=matrix(1,ncol = 2,byrow = TRUE);
  PropTreatRate=1
  PropExpoI=rep(PropTreatment, Nt)
  
  #  Initialiser MatTheta 
  MatTheta = array(NA, dim = c(4, 2, Nt))
  for (t in 1:(Nt)) {
    MatTheta[,,t] = matrix(c(PropTreatRate*PropExpoI[t], 0,
                             0, PropTreatRate*PropExpoI[t],
                             1-PropTreatRate*PropExpoI[t], 0,
                             0, 1-PropTreatRate*PropExpoI[t]), 
                           ncol = 2, byrow = TRUE)
  }
  
  ## calcul des strategies SMC
  Modelphi=TheSMCStrategy(t_begin_Camp, SMC_Number_of_cycle,
                          SMC_VectTime_between_cycles, SMC_Dur_cycle,
                          Gap, time, Nt, nb_years = Nyr)
  
  bar_phi=Modelphi$bar_phiSMC
  # Calcul de phi pour SMC
  phi=matrix(0,nrow = Na,ncol = Nt) # SMC exposition rate for human 
  bar_phi=(PropSMC/SMC_Dur_cycle)*bar_phi
  for (t in 1:Nt) {
    for (a in 1:Na) {
      if(Va[a]<=a_lim){phi[a,t]=bar_phi[t]}
      else{phi[a,t]=0}
    }
  }
  
  ## calcul des strategies MDA
  Modelpsi=TheMDAStrategy(t_begin_Camp, MDA_Number_of_cycle,
                          MDA_VectTime_between_cycles, MDA_Dur_cycle,
                          Gap, time, Nt, nb_years = Nyr)
  bar_psi=Modelpsi$bar_psiMDA
  # calcul de psi pour MDA
  psi_f=matrix(0,nrow = Na,ncol = Nt)
  psi_m=matrix(0,nrow = Na,ncol = Nt)
  bar_psi=(PropMDA/MDA_Dur_cycle)*bar_psi
  
  for (t in 1:Nt) {
    for (a in 1:Na) {
      if(Va[a]>a_lim){psi_m[a,t]=bar_psi[t]}
      else{psi_m[a,t]=0}
      if(Va[a]>a_lim & (Va[a]<=15 | Va[a]>=45)){psi_f[a,t]=bar_psi[t]}
      else{psi_f[a,t]=0}
    }
  }
  
  ## Initialize Variables
  {
    # Initialize Variables unexposed to SMC and MDA 
    mSh = matrix(0, nrow = Na, ncol = Nt); mAh = array(0, dim = c(2, Na, Nt));
    mIh = array(0, dim = c(4, Na, Nt)); mRh = array(0, dim = c(NsigmaR, Na, Nt));
    fSh = matrix(0, nrow = Na, ncol = Nt); fAh = array(0, dim = c(2, Na, Nt));
    fIh = array(0, dim = c(4, Na, Nt)); fRh = array(0, dim = c(NsigmaR, Na, Nt));
    Sm=rep(0,Nt); Im = matrix(0,nrow =2,ncol = Nt); Dh=rep(0,Nt);
    
    #Initialize variables exposed to SMC and MDA 
    mSh_smc = array(0, dim = c(NsigmaSmc, Na, Nt));
    mAh_smc = array(0, dim = c(2, NsigmaSmc,Na, Nt));
    fSh_smc = array(0, dim = c(NsigmaSmc, Na, Nt));
    fAh_smc = array(0, dim = c(2, NsigmaSmc,Na, Nt));
    
    mSh_mda = array(0, dim = c(NsigmaMda, Na, Nt));
    mAh_mda = array(0, dim = c(2,NsigmaMda, Na, Nt));
    fSh_mda = array(0, dim = c(NsigmaMda, Na, Nt));
    fAh_mda = array(0, dim = c(2,NsigmaMda, Na, Nt));
    
    Nh = numeric(Nt); Nm = numeric(Nt); lambdam = matrix(0,nrow = 2,ncol = Nt);
    lambdah = matrix(0,nrow = 2,ncol = Nt);
  }
  
  if (is.null(initial_conditions)) {
    ## Initial Conditions of variables unexposed to SMC and MDA ############
    # Wedge_h=30750 #30750, recrutment rate Bobo dioulasso in 2012
    # Initial values (t=1), Ref Quentin Richard et al Bobo Dioulasso 2012
    
    {
      
      Int_values = 8136.10 * (c(12.9, 12.5, 11.5, 13.1, 11.9, 9.3, 7.3, 5.5, 4.4,
                                3.2, 2.6, 1.8, 1.4, 0.9, 0.7, 0.3, 0.2, 0.2) / 0.997)
      # apparently the data only count for 99.7 % of the population
      # Age labels corresponding to each range in mSh_values
      age_label_3 = c("0 to 5 years", "5 to 10 years", "10 to 15 years",
                      "15 to 20 years", "20 to 25 years", "25 to 30 years",
                      "30 to 35 years", "35 to 40 years", "40 to 45 years",
                      "45 to 50 years", "50 to 55 years", "55 to 60 years",
                      "60 to 65 years", "65 to 70 years", "70 to 75 years",
                      "75 to 80 years", "80 to 85 years", "85 to 90 years")
      #Prev, Prevalence des infected, PrevI Prevalence des symptomatics
      #Prev2 Propoertion initial des gueris
      Prev2=0.15
      Prev = 0.4; PrevI = 1/4
      for (i in seq_along(age_label_3)) {
        if (i == 1) {
          # For the first age range (0 to 5 years)
          n_ages_group = sum(Va <= 5)
          mSh[Va <= 5, 1] = (1-Prev -Prev2)*(1-p)*Int_values[i]/n_ages_group
          fSh[Va <= 5, 1] = (1-Prev -Prev2)*p*Int_values[i]/n_ages_group
          mAh[, Va <= 5, 1] = (1 - PrevI)* Prev*c((1-p)*Int_values[i]*(1 - InitPrevR),    
                                                  (1-p)*Int_values[i]*InitPrevR)/n_ages_group
          fAh[, Va <= 5, 1] = (1 - PrevI)* Prev*c(p*Int_values[i]*(1 - InitPrevR),
                                                  p*Int_values[i]*InitPrevR) /n_ages_group
          mIh[, Va <= 5, 1] = PrevI* Prev*c((1-p)*Int_values[i]*(1-InitPrevR),
                                            (1-p)*Int_values[i]*InitPrevR,
                                            (1-p)*Int_values[i]*(1-InitPrevR),
                                            (1-p)*Int_values[i]*InitPrevR) /n_ages_group
          fIh[, Va <= 5, 1] = PrevI* Prev*c(p*Int_values[i]*(1-InitPrevR),
                                            p*Int_values[i]*InitPrevR,
                                            p*Int_values[i]*(1-InitPrevR),
                                            p*Int_values[i]*InitPrevR) /n_ages_group
          mRh[1, Va <= 5, 1] = Prev2*(1-p)*Int_values[i] /n_ages_group
          fRh[1, Va <= 5, 1] = Prev2*p*Int_values[i] /n_ages_group
        } else if (i == length(age_label_3)) {
          # For the last age range (> 85 years)
          n_ages_group = sum(Va > 85)
          mSh[Va > 85, 1] = (1-Prev -Prev2)*(1-p)*Int_values[i] /n_ages_group
          fSh[Va > 85, 1] = (1-Prev -Prev2)*p*Int_values[i] /n_ages_group
          mAh[ ,Va > 85, 1] = (1 - PrevI)* Prev*c((1-p)*Int_values[i]*(1 - InitPrevR),
                                                  (1-p)*Int_values[i]*InitPrevR) /n_ages_group
          fAh[, Va > 85, 1] = (1 - PrevI)* Prev*c(p*Int_values[i]*(1 - InitPrevR),
                                                  p*Int_values[i]*InitPrevR) /n_ages_group
          mIh[, Va > 85, 1] = PrevI* Prev*c((1-p)*Int_values[i]*(1-InitPrevR),
                                            (1-p)*Int_values[i]*InitPrevR,
                                            (1-p)*Int_values[i]*(1-InitPrevR),
                                            (1-p)*Int_values[i]*InitPrevR) /n_ages_group
          fIh[, Va > 85, 1] = PrevI* Prev*c(p*Int_values[i]*(1-InitPrevR),
                                            p*Int_values[i]*InitPrevR,
                                            p*Int_values[i]*(1-InitPrevR),
                                            p*Int_values[i]*InitPrevR) /n_ages_group
          mRh[1, Va > 85, 1] = Prev2*(1-p)*Int_values[i] /n_ages_group
          fRh[1, Va > 85, 1] = Prev2*p*Int_values[i] /n_ages_group
          
        } else {
          # For intermediate age ranges
          age_range = which(Va > (i - 1) * 5 & Va <= i * 5)
          n_ages_group = length(age_range)
          mSh[age_range, 1] = (1-Prev -Prev2)*(1-p)*Int_values[i] /n_ages_group
          fSh[age_range, 1] = (1-Prev -Prev2)*p*Int_values[i] /n_ages_group
          mAh[, age_range, 1] = (1 - PrevI)* Prev*c((1-p)*Int_values[i]*(1 - InitPrevR),
                                                    (1-p)*Int_values[i]*InitPrevR) /n_ages_group
          fAh[, age_range, 1] = (1 - PrevI)* Prev*c(p*Int_values[i]*(1 - InitPrevR),
                                                    p*Int_values[i]*InitPrevR) /n_ages_group
          mIh[, age_range, 1] = PrevI* Prev*c((1-p)*Int_values[i]*(1-InitPrevR),
                                              (1-p)*Int_values[i]*InitPrevR,
                                              (1-p)*Int_values[i]*(1-InitPrevR),
                                              (1-p)*Int_values[i]*InitPrevR) /n_ages_group
          fIh[, age_range, 1] = PrevI* Prev*c(p*Int_values[i]*(1-InitPrevR),
                                              p*Int_values[i]*InitPrevR,
                                              p*Int_values[i]*(1-InitPrevR),
                                              p*Int_values[i]*InitPrevR) /n_ages_group
          mRh[1, age_range, 1] = Prev2*(1-p)*Int_values[i] /n_ages_group
          fRh[1, age_range, 1] = Prev2*p*Int_values[i] /n_ages_group 
        }
      }
      
      Dh[1]=1
      S0m=10^6
      Prev_m0=0.1
      Sm[1] = (1 -Prev_m0) * S0m
      Im[, 1] = c((1 - InitPrevR) * Prev_m0 * S0m, InitPrevR * Prev_m0 * S0m)
      
      mSh_smc[, , 1] = matrix(0,nrow=NsigmaSmc,ncol=Na);
      mAh_smc[, , , 1] = array(0, dim = c(2, NsigmaSmc, Na));
      mSh_mda[, ,1] = matrix(0,nrow=NsigmaMda,ncol=Na);
      mAh_mda[, , , 1] = array(0, dim = c(2, NsigmaMda, Na)); 
      fSh_smc[, , 1] = matrix(0,nrow=NsigmaSmc,ncol=Na); 
      fAh_smc[, , , 1] = array(0, dim = c(2, NsigmaSmc, Na));
      fSh_mda[, ,1] = matrix(0,nrow=NsigmaMda,ncol=Na); 
      fAh_mda[, , , 1] = array(0, dim = c(2, NsigmaMda, Na));
    }
    
  } else {
    # Utiliser les conditions initiales fournies 
    mSh[, 1] = initial_conditions$mSh_final
    fSh[, 1] = initial_conditions$fSh_final
    mAh[, , 1] = initial_conditions$mAh_final
    fAh[, , 1] = initial_conditions$fAh_final
    mIh[, , 1] = initial_conditions$mIh_final
    fIh[, , 1] = initial_conditions$fIh_final
    mRh[, , 1] = initial_conditions$mRh_final
    fRh[, , 1] = initial_conditions$fRh_final
    
    mSh_smc[, , 1] = initial_conditions$mSh_smc_final
    fSh_smc[, , 1] = initial_conditions$fSh_smc_final
    mAh_smc[, , , 1] = initial_conditions$mAh_smc_final
    fAh_smc[, , , 1] = initial_conditions$fAh_smc_final
    
    mSh_mda[, , 1] = initial_conditions$mSh_mda_final
    fSh_mda[, , 1] = initial_conditions$fSh_mda_final
    mAh_mda[, , , 1] = initial_conditions$mAh_mda_final
    fAh_mda[, , , 1] = initial_conditions$fAh_mda_final
    
    Sm[1] = initial_conditions$Sm_final
    Im[, 1] = initial_conditions$Im_final
    Dh[1] = initial_conditions$Dh_final
  }
  
  
  t=1
  Nh[t] =  sum(mSh[,t]+fSh[,t]) + sum(mAh[, , t]+fAh[ , , t]) + sum(mIh[, , t]+fIh[ , , t]) + 
    sum(mRh[, , t]+fRh[ , , t]) + sum(mSh_mda[,,t]+fSh_mda[,,t]) + sum(mSh_smc[,,t]+fSh_smc[,,t]) + 
    sum(mAh_mda[,,,t]+fAh_mda[,,,t])  + sum(mAh_smc[,,,t]+fAh_smc[,,,t])
  Nm[t] = Sm[t] + sum(Im[ , t])
  
  for (t in 1:(Nt - 1)) {
    
    lambdam[,t] = betam %*% Im[ , t] / Nh[t]
    
    ## Calcul de lambdah
    lambdah[,t] = matrix(0, nrow = 2, ncol = 1)
    lambdah_sum = matrix(0, nrow = 2, ncol = 1)
    for(a in 1:Na) {
      term1 = betah[,,a] %*% (mAh[, a, t] + fAh[, a, t] + 
                                epsi %*% (mIh[, a, t] + fIh[, a, t]))
      #Terme 2: barbetahS[,,a] %*% (somme des valeurs SMC)
      term2 = matrix(0, nrow = 2, ncol = 1)
      # Correction: convertir la somme en vecteur colonne
      for (sigmaSmc in 1:NsigmaSmc) {
        smc_sum = mAh_smc[,sigmaSmc,a,t] + fAh_smc[,sigmaSmc,a,t]
        term2 = term2 + eta_S[sigmaSmc]*barbetahS[,,a] %*% smc_sum 
      }
      # Terme 3: barbetahM[,,a] %*% (somme des valeurs MDA)
      term3 = matrix(0, nrow = 2, ncol = 1)
      for (sigmaMda in 1:NsigmaMda) {
        mda_sum = mAh_mda[,sigmaMda,a,t] + fAh_mda[,sigmaMda,a,t]
        term3 = term3+ eta_M[sigmaMda]*barbetahM[,,a] %*% mda_sum
      }
      # Addition des termes
      lambdah_sum = lambdah_sum + term1 + term2 + term3
    }   
    # Calcul final de lambdah
    lambdah[,t] = lambdah_sum / Nh[t]
    
    if(t >= 148 && t <= 220) {  
      cat("Temps t =", t, "\n")
      N_hTrait= (sum(mAh[,,t] + fAh[,,t])+sum(mIh[, , t]+fIh[ , , t])+sum(mSh[,t]+fSh[,t])+
                   sum(mRh[, , t]+fRh[ , , t]))/Nh[t] 
      cat("Population non traitée =", N_hTrait , "\n")
      SMC_Trait= (sum(mSh_smc[,,t]+fSh_smc[,,t]) +sum(mAh_smc[,,,t]+fAh_smc[,,,t]))/Nh[t]
      cat("Population SMC  =", SMC_Trait, "\n")
      MDA_Trait= (sum(mSh_mda[,,t]+fSh_mda[,,t]) +sum(mAh_mda[,,,t]+fAh_mda[,,,t]))/Nh[t]
      cat("Population MDA  =", MDA_Trait, "\n")
      cat("Taux d'infection total lambdah=", sum(lambdah[,t]), "\n")
      cat("Contribution normale =", sum(term1), "\n")
      cat("Contribution SMC =", sum(term2), "\n") 
      cat("Contribution MDA =", sum(term3), "\n") 
      cat("Taux d'infection total lambdam=", sum(lambdam[,t]), "\n")
      cat("prop SMC =", bar_phi[t], "\n")
      cat("prop MDA =", bar_psi[t], "\n")
      cat("-----\n")
    }
    
    if(t >= 148 && t <= 220) {
      # Flux entrants vers SMC/MDA
      flux_in_smc = sum(phi[,t] * (mSh[,t] + fSh[,t]) )+ sum(phi[,t] * (mAh[,,t] + fAh[,,t])) +
        sum(phi[,t] * (mRh[,,t] + fRh[,,t]))
      flux_in_mda = sum(psi_m[,t] * mSh[,t]) +sum(psi_m[,t] *(mAh[,,t])) + sum(psi_m[,t] *(mRh[,,t])) + 
        sum(psi_f[,t] * (fSh[,t])) +sum(psi_f[,t] * (fAh[,,t])) + sum(psi_f[,t] * (fRh[,,t]))
      
      # Flux sortants de SMC/MDA
      flux_out_smc = sum(k_S * (mSh_smc[,,t] + fSh_smc[,,t])) + sum(k_S * (mAh_smc[,,,t] + fAh_smc[,,,t]))
      flux_out_mda = sum(k_M* (mSh_mda[,,t] + fSh_mda[,,t])) + sum(k_M* (mAh_mda[,,,t] + fAh_mda[,,,t]))
      
      cat("Flux entrant SMC:", flux_in_smc, "- Flux sortant SMC:", flux_out_smc, "\n")
      cat("Flux entrant MDA:", flux_in_mda, "- Flux sortant MDA:", flux_out_mda, "\n")
    }
    
    ##Conditions aux bords 
    a=1;
    mSh[a,t+1]=((mSh[a,t]/dt)+((1-p)*Lh/da)+sum(rhoh[]*mRh[,a,t])+
                  sum(k_S[]*mSh_smc[,a,t])+sum(k_M[]*mSh_mda[,a,t]))/
      ((1/dt)+(1/da)+muh[a]+sum(lambdam[,t])+phi[a,t]+psi_m[a,t]);
    fSh[a,t+1]=((fSh[a,t]/dt)+(p*Lh/da)+sum(rhoh[]*fRh[,a,t])+
                  sum(k_S[]*fSh_smc[,a,t])+sum(k_M[]*fSh_mda[,a,t]))/
      ((1/dt)+(1/da)+muh[a]+sum(lambdam[,t])+phi[a,t]+psi_f[a,t]);
    mAh[,a,t+1]=c(0,0); mIh[,a,t+1]=c(0,0,0,0); mRh[,a,t+1]=rep(0,NsigmaR); 
    fAh[,a,t+1]=c(0,0); fIh[,a,t+1]=c(0,0,0,0); fRh[,a,t+1]=rep(0,NsigmaR);
    mSh_smc[,a,t+1]= rep(0,NsigmaSmc); fSh_smc[,a,t+1]= rep(0,NsigmaSmc);
    mAh_smc[,,a,t+1]= array(0, dim = c(2, NsigmaSmc));
    fAh_smc[,,a,t+1]= array(0, dim = c(2, NsigmaSmc));
    mSh_mda[,a,t+1]= rep(0,NsigmaMda); fSh_mda[,a,t+1]= rep(0,NsigmaMda);
    mAh_mda[,,a,t+1]= array(0, dim = c(2, NsigmaMda));
    fAh_mda[,,a,t+1]= array(0, dim = c(2, NsigmaMda));
    
    for (a in 2:Na) {
      mSh[a,t+1]=((mSh[a,t]/dt)+(mSh[a-1,t+1]/da)+sum(rhoh[]*mRh[,a,t])+
                    sum(k_S[]*mSh_smc[,a,t])+sum(k_M[]*mSh_mda[,a,t]))/((1/dt)+(1/da)+
                                                                          muh[a]+sum(lambdam[,t])+phi[a,t]+psi_m[a,t])
      fSh[a,t+1]=((fSh[a,t]/dt)+(fSh[a-1,t+1]/da)+sum(rhoh[]*fRh[,a,t])+
                    sum(k_S[]*fSh_smc[,a,t])+sum(k_M[]*fSh_mda[,a,t]))/((1/dt)+(1/da)+
                                                                          muh[a]+sum(lambdam[,t])+phi[a,t]+psi_f[a,t])
      mAh[,a,t+1]=(solve(diag(2)*((1/dt)+(1/da)+muh[a]+nuh_m[a]+ phi[a,t]+psi_m[a,t])))%*%
        ((mAh[,a,t]/dt)+(mAh[,a-1,t+1]/da)+lambdam[,t]*mSh[a,t]+ 
           rowSums(sweep(mAh_smc[,1:NsigmaSmc,a,t], 2, k_S[1:NsigmaSmc], "*")))
      fAh[,a,t+1]= (solve(diag(2)*((1/dt)+(1/da)+muh[a]+nuh_f[a]+ phi[a,t]+psi_f[a,t])))%*%
        ((fAh[,a,t]/dt)+(fAh[,a-1,t+1]/da)+lambdam[,t]*fSh[a,t]+ 
           rowSums(sweep(fAh_smc[,1:NsigmaSmc,a,t], 2, k_S[1:NsigmaSmc], "*")))
      mIh[,a,t+1]=(solve(diag(4)*((1/dt)+(1/da)+muh[a])+deltah[,,a]+gammah[,,a]))%*%
        ((mIh[,a,t]/dt)+(mIh[,a-1,t+1]/da)+nuh_m[a]*MatTheta[,,t]%*%mAh[,a,t])
      fIh[,a,t+1]=(solve(diag(4)*((1/dt)+(1/da)+muh[a])+deltah[,,a]+gammah[,,a]))%*%
        ((fIh[,a,t]/dt)+(fIh[,a-1,t+1]/da)+nuh_f[a]*MatTheta[,,t]%*%fAh[,a,t])
      
      ##Condition aux bords de sigma     
      sigmaR=1;
      mRh[sigmaR,a,t+1]=((mRh[sigmaR,a,t]/dt)+(mRh[sigmaR,a-1,t+1]/da)+(f%*%gammah[,,a]%*%mIh[, a,t]+
                                                                          m%*%rowSums(sweep(mAh_mda[,1:NsigmaMda,a,t], 2, k_M[1:NsigmaMda], "*")))/dsigmaR)/
        ((1/dt)+(1/da)+(1/dsigmaR)+muh[a]+rhoh[sigmaR]+phi[a,t]+psi_m[a,t]);
      fRh[sigmaR,a,t+1]=((fRh[sigmaR,a,t]/dt)+(fRh[sigmaR,a-1,t+1]/da)+(f%*%gammah[,,a]%*%fIh[, a,t]+
                                                                          m%*%rowSums(sweep(fAh_mda[,1:NsigmaMda,a,t], 2, k_M[1:NsigmaMda], "*")))/dsigmaR)/
        ((1/dt)+(1/da)+(1/dsigmaR)+ muh[a]+rhoh[sigmaR]+phi[a,t]+psi_f[a,t]);
      
      for (sigmaR in 2:NsigmaR) {
        mRh[sigmaR,a,t+1]=((mRh[sigmaR,a,t]/dt)+(mRh[sigmaR,a-1,t+1]/da)+
                             (mRh[sigmaR-1,a,t+1]/dsigmaR))/((1/dt)+(1/da)+(1/dsigmaR)+
                                                               muh[a]+rhoh[sigmaR]+phi[a,t]+psi_m[a,t])
        fRh[sigmaR,a,t+1]=((fRh[sigmaR,a,t]/dt)+(fRh[sigmaR,a-1,t+1]/da)+
                             (fRh[sigmaR-1,a,t+1]/dsigmaR))/((1/dt)+(1/da)+(1/dsigmaR)+
                                                               muh[a]+rhoh[sigmaR]+phi[a,t]+psi_f[a,t])
      } 
      
      sigmaSmc=1;
      mSh_smc[sigmaSmc,a,t+1]=((mSh_smc[sigmaSmc,a,t]/dt)+(mSh_smc[sigmaSmc,a-1,t+1]/da)+
                                 phi[a,t]*(mSh[a,t]+sum(mRh[,a,t]))/dsigmaSmc)/
        ((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+eta_S[sigmaSmc]*sum(lambdam[,t])+k_S[sigmaSmc]);
      fSh_smc[sigmaSmc,a,t+1]=((fSh_smc[sigmaSmc,a,t]/dt)+(fSh_smc[sigmaSmc,a-1,t+1]/da)+
                                 phi[a,t]*(fSh[a,t]+sum(fRh[,a,t]))/dsigmaSmc)/
        ((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+eta_S[sigmaSmc]*sum(lambdam[,t])+k_S[sigmaSmc]);
      mAh_smc[,sigmaSmc,a,t+1]=(solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+k_S[sigmaSmc])))%*%
        ((mAh_smc[,sigmaSmc,a,t]/dt)+(mAh_smc[,sigmaSmc,a-1,t+1]/da)+
           (phi[a,t]*mAh[,a,t]/dsigmaSmc)+ eta_S[sigmaSmc]*lambdam[,t]*mSh_smc[sigmaSmc,a,t])
      fAh_smc[,sigmaSmc,a,t+1]= (solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+k_S[sigmaSmc])))%*%
        ((fAh_smc[,sigmaSmc,a,t]/dt)+(fAh_smc[,sigmaSmc,a-1,t+1]/da)+
           (phi[a,t]*fAh[,a,t]/dsigmaSmc)+ eta_S[sigmaSmc]*lambdam[,t]*fSh_smc[sigmaSmc,a,t])
      
      for(sigmaSmc in 2:NsigmaSmc){
        mSh_smc[sigmaSmc,a,t+1]=((mSh_smc[sigmaSmc,a,t]/dt)+ (mSh_smc[sigmaSmc,a-1,t+1]/da)+
                                   (mSh_smc[sigmaSmc-1,a,t+1]/dsigmaSmc))/
          ((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+eta_S[sigmaSmc]*sum(lambdam[,t])+k_S[sigmaSmc])
        fSh_smc[sigmaSmc,a,t+1]= ((fSh_smc[sigmaSmc,a,t]/dt)+  (fSh_smc[sigmaSmc,a-1,t+1]/da)+
                                    (fSh_smc[sigmaSmc-1,a,t+1]/dsigmaSmc))/
          ((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+eta_S[sigmaSmc]*sum(lambdam[,t])+k_S[sigmaSmc])
        mAh_smc[,sigmaSmc,a,t+1]=(solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+k_S[sigmaSmc])))%*%
          ((mAh_smc[,sigmaSmc,a,t]/dt)+ (mAh_smc[,sigmaSmc,a-1,t+1]/da)+
             (mAh_smc[,sigmaSmc-1,a,t+1]/dsigmaSmc)+ eta_S[sigmaSmc]*lambdam[,t]*mSh_smc[sigmaSmc,a,t])
        fAh_smc[,sigmaSmc,a,t+1]=(solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaSmc)+muh[a]+k_S[sigmaSmc])))%*%
          ((fAh_smc[,sigmaSmc,a,t]/dt)+ (fAh_smc[,sigmaSmc,a-1,t+1]/da)+
             (fAh_smc[,sigmaSmc-1,a,t+1]/dsigmaSmc)+ eta_S[sigmaSmc]*lambdam[,t]*fSh_smc[sigmaSmc,a,t])
      }    
      sigmaMda=1;
      mSh_mda[sigmaMda,a,t+1]=((mSh_mda[sigmaMda,a,t]/dt)+ (mSh_mda[sigmaMda,a-1,t+1]/da)+
                                 (psi_m[a,t]*(mSh[a,t]+sum(mRh[,a,t]))/dsigmaMda))/
        ((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+eta_M[sigmaMda]*
           sum(lambdam[,t])+k_M[sigmaMda]);
      fSh_mda[sigmaMda,a,t+1]=((fSh_mda[sigmaMda,a,t]/dt)+(fSh_mda[sigmaMda,a-1,t+1]/da)+
                                 (psi_f[a,t]*(fSh[a,t]+sum(fRh[,a,t]))/dsigmaMda))/
        ((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+eta_M[sigmaMda]*
           sum(lambdam[,t])+k_M[sigmaMda]);
      mAh_mda[,sigmaMda,a,t+1]=(solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+k_M[sigmaMda])))%*%
        ((mAh_mda[,sigmaMda,a,t]/dt)+(mAh_mda[,sigmaMda,a-1,t+1]/da)+
           (psi_m[a,t]*mAh[,a,t]/dsigmaMda)+ eta_M[sigmaMda]*lambdam[,t]*mSh_mda[sigmaMda,a,t])
      fAh_mda[,sigmaMda,a,t+1]= (solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+k_M[sigmaMda])))%*%
        ((fAh_mda[,sigmaMda,a,t]/dt)+(fAh_mda[,sigmaMda,a-1,t+1]/da)+
           (psi_f[a,t]*fAh[,a,t]/dsigmaMda)+ eta_M[sigmaMda]*lambdam[,t]*fSh_mda[sigmaMda,a,t])
      
      ##With MDA 
      for(sigmaMda in 2:NsigmaMda){
        mSh_mda[sigmaMda,a,t+1]=((mSh_mda[sigmaMda,a,t]/dt)+(mSh_mda[sigmaMda,a-1,t+1]/da)+
                                   (mSh_mda[sigmaMda-1,a,t+1]/dsigmaMda))/
          ((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+eta_M[sigmaMda]* sum(lambdam[,t])+k_M[sigmaMda])
        fSh_mda[sigmaMda,a,t+1]=((fSh_mda[sigmaMda,a,t]/dt)+(fSh_mda[sigmaMda,a-1,t+1]/da)+
                                   (fSh_mda[sigmaMda-1,a,t+1]/dsigmaMda))/
          ((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+eta_M[sigmaMda]*sum(lambdam[,t])+k_M[sigmaMda])
        mAh_mda[,sigmaMda,a,t+1]=(solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+k_M[sigmaMda])))%*%
          ((mAh_mda[,sigmaMda,a,t]/dt)+ (mAh_mda[,sigmaMda,a-1,t+1]/da)+
             (mAh_mda[,sigmaMda-1,a,t+1]/dsigmaMda)+eta_M[sigmaMda]*lambdam[,t]* mSh_mda[sigmaMda,a,t]) 
        fAh_mda[,sigmaMda,a,t+1]= (solve(diag(2)*((1/dt)+(1/da)+(1/dsigmaMda)+muh[a]+k_M[sigmaMda])))%*%
          ((fAh_mda[,sigmaMda,a,t]/dt)+ (fAh_mda[,sigmaMda,a-1,t+1]/da)+
             (fAh_mda[,sigmaMda-1,a,t+1]/dsigmaMda)+eta_M[sigmaMda]*lambdam[,t]*fSh_mda[sigmaMda,a,t])
      }
    }
    
    Dh_increment = 0
    for(a in 1:Na) {
      Dh_increment = Dh_increment + sum(deltah[,,a] %*% (mIh[,a,t] + fIh[,a,t]))
    }
    Dh[t+1] = Dh[t] + dt * Dh_increment
    Sm[t+1]= (Sm[t]/dt+Lm)/(1/dt+mum+ sum(lambdah[,t]))
    Im[,t+1]= (solve(diag(2)/dt+deltam))%*%(Im[,t]/dt+ lambdah[,t]*Sm[t])
    
    Nh[t+1] = sum(mSh[,t+1]+fSh[,t+1]) + sum(mAh[, , t+1]+fAh[ , , t+1]) + sum(mIh[, , t+1]+fIh[ , , t+1]) + 
      sum(mRh[, , t+1]+fRh[ , , t+1]) + sum(mSh_mda[,,t+1]+fSh_mda[,,t+1]) + sum(mSh_smc[,,t+1]+fSh_smc[,,t+1]) + 
      sum(mAh_mda[,,,t+1]+fAh_mda[,,,t+1])  + sum(mAh_smc[,,,t+1]+fAh_smc[,,,t+1])
    
    Nm[t+1] = Sm[t+1] + sum(Im[ , t+1])
    
  }
  
  mShTot=numeric(Nt);mIhTot=numeric(Nt);mRhTot=numeric(Nt);mAhTot=numeric(Nt);
  fShTot=numeric(Nt);fIhTot=numeric(Nt);fRhTot=numeric(Nt);fAhTot=numeric(Nt);
  ShTot=numeric(Nt);AhTot=numeric(Nt);IhTot=numeric(Nt);RhTot=numeric(Nt);
  ResistantHumanCases=numeric(Nt);ImTot= numeric(Nt);
  TotPopSmc=numeric(Nt); TotPopMda=numeric(Nt);Totlambdah=numeric(Nt)
  
  for (t in 1: Nt) {
    
    ResistantHumanCases[t]= sum(mAh[2,,t]+fAh[2,,t]+mIh[2,,t]+mIh[4,,t]+
                                  fIh[2,,t]+fIh[4,,t])+ sum(mAh_smc[2,,,t])+
      sum(fAh_smc[2,,,t])+sum(mAh_mda[2,,,t])+sum(fAh_mda[2,,,t])
    ShTot[t] = sum(mSh[,t]) + sum(fSh[,t]) + sum(mSh_smc[,,t]) + sum(fSh_smc[,,t]) 
    +sum(mSh_mda[,,t]) + sum(fSh_mda[,,t])
    AhTot[t]=sum(mAh[1,,t]) + sum(mAh[2,,t]) + sum(fAh[1,,t]) + sum(fAh[2,,t]) +
      sum(mAh_smc[1,,,t]) + sum(mAh_smc[2,,,t]) +sum(fAh_smc[1,,,t]) +
      sum(fAh_smc[2,,,t]) + sum(mAh_mda[1,,,t]) + sum(mAh_mda[2,,,t]) + 
      sum(fAh_mda[1,,,t]) + sum(fAh_mda[2,,,t])
    IhTot[t]=sum(mIh[1,,t]) + sum(mIh[2,,t]) + sum(mIh[3,,t]) + sum(mIh[4,,t]) +
      sum(fIh[1,,t]) + sum(fIh[2,,t]) + sum(fIh[3,,t]) + sum(fIh[4,,t]);
    RhTot[t]= sum(mRh[,,t]) + sum(fRh[,,t])
    
    TotPopSmc[t]= sum(mAh_smc[1,,,t]) + sum(mAh_smc[2,,,t]) + sum(fAh_smc[1,,,t])
    + sum(fAh_smc[2,,,t])+sum(mSh_smc[,,t]) + sum(fSh_smc[,,t])
    TotPopMda[t]= sum(mAh_mda[1,,,t]) + sum(mAh_mda[2,,,t]) + sum(fAh_mda[1,,,t])
    + sum(fAh_mda[2,,,t])+sum(mSh_mda[,,t]) + sum(fSh_mda[,,t])
    Totlambdah[t]=sum(lambdah[1,t]+lambdah[2,t])
    
    ImTot[t]=Im[1,t]+Im[2,t]
  }
  
  #proportion
  PropSh=numeric(Nt); PropAh=numeric(Nt); PropIh=numeric(Nt); PropRh=numeric(Nt);
  PropSm=numeric(Nt); PropIm=numeric(Nt);
  #lambdam=numeric(Nt); lambdah=numeric(Nt);
  
  PropSh= ShTot/Nh; PropAh= AhTot/Nh; 
  PropIh= IhTot/Nh; PropRh= RhTot/Nh; 
  
  PropSm= Sm/Nm; PropIm= ImTot/Nm; 
  
  listOutput = list("mSh"=mSh, "mAh" =mAh, "mIh"=mIh, "mRh"=mRh, "fSh"=fSh,
                    "fAh" =fAh, "fIh"=fIh, "fRh"=fRh, "Sm"=Sm, "Im"=Im,
                    "mRh" =mRh,  "fRh" =fRh, "ShTot"= ShTot, "AhTot"=AhTot,
                    "IhTot"=IhTot, "RhTot"=RhTot, "ImTot"=ImTot, "Dh"=Dh,
                    "CumulDh"=Dh[Nt], "RhCases"=ResistantHumanCases,
                    "PrevRhCases"=ResistantHumanCases[Nt], "time"=time,
                    "Nh"=Nh, "Nm"=Nm, "TotPopSmc"=TotPopSmc, "TotPopMda"=TotPopMda,
                    "lambdam"=lambdam, "lambdah"=lambdah, "Totlambdah"=Totlambdah,
                    "PropSh"=PropSh,"PropAh"=PropAh, "PropIh"=PropIh,
                    "PropRh"=PropRh,"PropSm"=PropSm, "PropIm"=PropIm,
                    "mSh_final" = mSh[, Nt],"fSh_final" = fSh[, Nt], 
                    "mAh_final" = mAh[, , Nt],"fAh_final" = fAh[, , Nt],
                    "mIh_final" = mIh[, , Nt],"fIh_final" = fIh[, , Nt],
                    "mRh_final" = mRh[, , Nt],"fRh_final" = fRh[, , Nt],
                    "mSh_smc_final" = mSh_smc[, , Nt],"fSh_smc_final" = fSh_smc[, , Nt],
                    "mAh_smc_final" = mAh_smc[, , , Nt],"fAh_smc_final" = fAh_smc[, , , Nt],
                    "mSh_mda_final" = mSh_mda[, , Nt],"fSh_mda_final" = fSh_mda[, , Nt],
                    "mAh_mda_final" = mAh_mda[, , , Nt],"fAh_mda_final" = fAh_mda[, , , Nt],
                    "Sm_final" = Sm[Nt],"Im_final" = Im[, Nt],"Dh_final" = Dh[Nt])
  
  return(listOutput)  
}

###############TAB pour ANOVA#############################
{
  setwd("C:/Users/Administrator/Documents/Nouveau Projet R")
  
  VectPropTreatment=c(0.1,0.5,0.9)
  VectMDA_Number_of_cycle=c(3,4,5)
  VectPropMDA=c(0.1,0.5,0.9)
  VectPropSMC=c(0.1,0.5,0.9)
  Vectp_f=c(0.1,0.3)
  Vecta_lim=c(5,10)
  
  # CORRECTION: Remplacer 'a_lim' par 'Vecta_lim'
  SizeTab=length(VectPropTreatment)*length(VectPropMDA)*length(VectMDA_Number_of_cycle)*
    length(VectPropSMC)*length(Vectp_f)*length(Vecta_lim)
  
  SizeTab
  
  ColPropTreatment=rep(NA,SizeTab)
  ColPropMDA=rep(NA,SizeTab)
  ColMDA_Number_of_cycle=rep(NA,SizeTab)
  ColPropSMC=rep(NA,SizeTab)
  Colp_f=rep(NA,SizeTab)
  Cola_lim=rep(NA,SizeTab)
  
  ColDeath=rep(NA,SizeTab)
  ColPrevResisH=rep(NA,SizeTab)
  
  # Nyr = 5; LengYr = 365; Lm = 2.1*10^5; Tmax = Nyr * LengYr; t_begin_Camp=120;
   InitPrevR=1/1000; p=0.35
  
  RunNumber=0
  
  for (PropTreatment in VectPropTreatment) {
    for (PropMDA in VectPropMDA) {
      for (MDA_Number_of_cycle in VectMDA_Number_of_cycle) {
        for (PropSMC in VectPropSMC) {
          for (p_f in Vectp_f) {
            for (a_lim in Vecta_lim) {
              
              print(SizeTab-RunNumber)
              RunNumber=RunNumber+1  
              
              parameters=list(p_f = 0.15, a_lim=5, Nyr = 1, LengYr = 365,
                              Lm = 2.1*10^5, t_begin_Camp=120, PropTreatment=0.7, Gap=1, dt =1,
                              SMC_Number_of_cycle=1, SMC_VectTime_between_cycles=c(0,0,0), SMC_Dur_cycle=7,
                              PropSMC=0, MDA_Number_of_cycle=1, MDA_VectTime_between_cycles=c(0,0,0),
                              MDA_Dur_cycle=7, PropMDA=0,
                              initial_conditions = NULL)
              
              # Ajout des paramètres manquants probables
              parameters$Tmax = parameters$Nyr * parameters$LengYr
              
              Test_1= do.call(Mainfunction, parameters )
              
              initial_conditions_year2_6 = list(mSh_final = Test_1$mSh_final,fSh_final = Test_1$fSh_final,
                                                mAh_final =Test_1$mAh_final,fAh_final = Test_1$fAh_final,mIh_final = Test_1$mIh_final,
                                                fIh_final =Test_1$fIh_final, mRh_final = Test_1$mRh_final,fRh_final = Test_1$fRh_final,
                                                mSh_smc_final = Test_1$mSh_smc_final,fSh_smc_final = Test_1$fSh_smc_final,
                                                mAh_smc_final = Test_1$mAh_smc_final,fAh_smc_final = Test_1$fAh_smc_final,
                                                mSh_mda_final = Test_1$mSh_mda_final,fSh_mda_final = Test_1$fSh_mda_final,
                                                mAh_mda_final = Test_1$mAh_mda_final,fAh_mda_final = Test_1$fAh_mda_final,
                                                Sm_final = Test_1$Sm_final,Im_final = Test_1$Im_final, Dh_final = Test_1$Dh_final
              )
              
              parameters_2=list( p_f =  p_f, a_lim =a_lim, Nyr = 5, LengYr = 365,
                                 Lm = 2.1*10^5, t_begin_Camp=120, PropTreatment=PropTreatment, Gap=1, dt =1,
                                 SMC_Number_of_cycle=c(3,4,5), SMC_VectTime_between_cycles=c(30,30,30), SMC_Dur_cycle=7,
                                 PropSMC=PropSMC, MDA_Number_of_cycle=MDA_Number_of_cycle, MDA_VectTime_between_cycles=c(30,30,30),
                                 MDA_Dur_cycle=7, PropMDA=PropMDA,
                                 initial_conditions=initial_conditions_year2_6 )
              
              # Ajout des paramètres manquants
              parameters_2$Tmax = parameters_2$Nyr * parameters_2$LengYr
              
              Test_2= do.call(Mainfunction, parameters_2)
              
              ColPropTreatment[RunNumber]=PropTreatment
              ColPropMDA[RunNumber]=PropMDA
              ColMDA_Number_of_cycle[RunNumber]=MDA_Number_of_cycle
              ColPropSMC[RunNumber]=PropSMC
              Colp_f[RunNumber]=p_f
              Cola_lim[RunNumber]=a_lim
              
              ColDeath[RunNumber]=Test_2$CumulDh
              ColPrevResisH[RunNumber]=Test_2$PrevRhCases
              
              # Gestion des erreurs (optionnel)
              if (is.null(Test_2$CumulDh) || is.null(Test_2$PrevRhCases)) {
                warning(paste("Résultats manquants pour le run", RunNumber))
              }
            }
          }
        }
      }
    }
  }
  
  Tableau<-data.frame(PropTreatment=ColPropTreatment,PropMDA=ColPropMDA,
                      MDA_Number_of_cycle=ColMDA_Number_of_cycle, PropSMC=ColPropSMC,
                      p_f=Colp_f, a_lim=Cola_lim, Death=ColDeath,PrevR=ColPrevResisH)
  
  # Vérification des résultats
  print(paste("Nombre de runs complétés:", RunNumber))
  print(paste("Nombre de valeurs manquantes dans Death:", sum(is.na(ColDeath))))
  print(paste("Nombre de valeurs manquantes dans PrevR:", sum(is.na(ColPrevResisH))))
  
  write.table(Tableau, file = "TabMDA.csv", sep = ",", row.names = F)
  
  # Affichage des premières lignes pour vérification
  print("Premières lignes du tableau:")
  print(head(Tableau))
  
  # Retourner le tableau pour utilisation ultérieure
  return(Tableau)
  
}

getwd()