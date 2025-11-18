data {
  int<lower = 1> n_site; // number of random effect levels (sites) 
 
  // micro
  int<lower = 1> N; // Sample size for micro data 
  array[N] int<lower=1, upper=n_site> micro_site; // id of random effect (sites)
  
  vector[N] yeari;
  vector[N] yMicroi; // Observed micro
  
   // macro
  int<lower = 1> Nm; // Sample size for macro data
  array[Nm] int<lower=1, upper=n_site> macro_site;  
  
  vector[Nm] yMacroi; //
  
  vector[Nm] tempi; // predictor temp
  vector[Nm] LAIi; // predictor LAI
  // vector[Nm] macroAi; //  predictor macroA
  vector[Nm] year2i;
  
  vector[Nm] depthi;
//  vector[Nm] bedAi;


}

parameters{
  real mu_grand; // grand mean for micro value 
  array[n_site] real alphamuSt; 
// real muSt;
  real <lower = 0> sigmaSt;
//vector[n_site] b_muSt_ncp; // site offsets
  

  //array[n_site] real alphayear;
  real muyear;
  real<lower = 0> sigmayear;
  vector[n_site] b_year_ncp;
  
  real<lower = 0> sigma_microy; // sd general

  //macro
  real mu_grandM;
  
 // array[n_site] real alphayear2;
  real muyear2;
  real<lower = 0> sigmayear2;
  vector[n_site] b_year2_ncp;
  
  real muDepth;
  //real<lower = 0> sigmaDepth;
 // vector[n_site] b_depth_ncp;
  // 
  // real mubedA;
  // real<lower = 0> sigmabedA;
  // vector[n_site] b_bedA_ncp;
  
 
 //array[n_site] real alphatempSt;
 real mutempSt;
 real<lower = 0> sigmatempSt;
 vector[n_site] b_temp_ncp;
 
  // array[n_site] real alphaLAISt;
  real muLAISt;
  real<lower = 0> sigmaLAISt;
  vector[n_site] b_LAI_ncp;

 // array[n_site] real alphamacroASt;
 // real mumacroASt;
 // real<lower = 0> sigmamacroASt;
 // vector[n_site] b_macroA_ncp;

 
  array[n_site] real alphaMacroSt;
  //real muMacroSt;
  real<lower = 0> sigmaMacroSt;
  //vector[n_site] b_MacroSite_ncp;
  
  
  real betaMicroxtemp;
  real betaMicroxLAI;
  // real betaMicroxmacroA;
  
  real<lower = 0> sigmaMacro_y;
}

transformed parameters{
  // Micros
  vector[N] y_hat; 
  
  vector[n_site] b_year;
  b_year = muyear + sigmayear*b_year_ncp;
  // vector[n_site] b_site;
  //b_site = muSt + sigmaSt * b_muSt_ncp; 
  //macro
  // array[n_site] real betayear2;     //site level beta temp
  array[n_site] real betatempSt;     //site level beta temp
  // array[n_site] real betamacroASt;     //site level beta macroA
  array[n_site] real betaLAISt;     //site level beta LAI

  for (i in 1:N){
    y_hat[i] =  mu_grand + alphamuSt[micro_site[i]] + b_year[micro_site[i]]*yeari[i]//alphayear[micro_site[i]] * yeari[i] 
    ;
  }
  
  // macro
  vector[n_site] b_year2;
  vector[n_site] b_temp;
  vector[n_site] b_LAI;
  // vector[n_site] b_macroA;
  
 // vector[n_site] b_depth;
  // vector[n_site] b_bedA;

  // vector[n_site] b_MacroSite;

  
  b_year2 = muyear2 + sigmayear2 * b_year2_ncp;
  b_temp = mutempSt + sigmatempSt * b_temp_ncp;
  b_LAI = muLAISt + sigmaLAISt * b_LAI_ncp;
  // b_macroA = mumacroASt + sigmamacroASt * b_macroA_ncp;
    
//  b_depth = muDepth + sigmaDepth * b_depth_ncp;
  // b_bedA = mubedA + sigmabedA * b_bedA_ncp;

  // b_MacroSite = muMacroSt + sigmaMacroSt * b_MacroSite_ncp; 
   
  for (iSt in 1:n_site){
    betatempSt[iSt] = b_temp[iSt] + betaMicroxtemp * (mu_grand + alphamuSt[iSt]);
  }
  for (iSt in 1:n_site){
    betaLAISt[iSt] = b_LAI[iSt] + betaMicroxLAI * (mu_grand + alphamuSt[iSt]);
  }
  // for (iSt in 1:n_site){
  //   betamacroASt[iSt] = b_macroA[iSt] + betaMicroxmacroA * (mu_grand + alphamuSt[iSt]);
  // }

}

model{
  // Micros
  //// likelihood
  
  alphamuSt ~ normal(0, sigmaSt);
  // b_muSt ~ normal(muSt, sigma_St);
  //b_site ~ normal(0,1);
  mu_grand ~ normal(0, 5);
  // muSt ~ normal(0, 50);// beta(1,1); -3, 0.5
  // alphayear ~ normal(muyear, sigmayear);
  muyear ~ normal(0, 5);
  sigmayear ~ normal(0, 5);
  
  sigmaSt ~ normal(0, 5);
  sigma_microy ~ normal(0, 20);
  
  yMicroi ~ normal(y_hat, sigma_microy);
  
  b_year_ncp ~ normal(0,1);
  
 // macro
  // likelihood
  for (i in 1:Nm){
    yMacroi[i] ~ normal(mu_grandM + alphaMacroSt[macro_site[i]] +
    betatempSt[macro_site[i]] * tempi[i] + 
   // betamacroASt[macro_site[i]] * macroAi[i] +
    betaLAISt[macro_site[i]] * LAIi[i]  
    + b_year2[macro_site[i]] * year2i[i] 
     + muDepth * depthi[i] 
   //   + b_bedA[macro_site[i]] * bedAi[i] 
   , sigmaMacro_y);
  }
  
  alphaMacroSt ~ normal(0, sigmaMacroSt);
  // alphatempSt ~ normal(mutempSt, sigmatempSt);
  // alphaLAISt ~ normal(muLAISt, sigmaLAISt);
  // alphamacroASt ~ normal(mumacroASt, sigmamacroASt);
  
  //alphayear2 ~ normal(muyear2, sigmayear2);

  
  //// priors
  mu_grandM ~ normal(0, 5);
  
  // muMacroSt ~ normal(0,50);
  sigmaMacroSt ~ normal(0, 5);
  
  sigmaMacro_y ~ normal(0, 5);

  mutempSt ~ normal(0, 5);
  sigmatempSt ~ normal(0, 5);

  muLAISt ~ normal(0, 5);
  sigmaLAISt ~ normal(0, 5);
  
  // mumacroASt ~ normal(0,20);
  // sigmamacroASt ~ normal(0,20);

  muyear2 ~ normal(0, 5);
  sigmayear2 ~ normal(0, 5);
  
  muDepth ~ normal(0,5);
  //sigmaDepth ~ normal(0,20);
  
  // mubedA ~ normal(0,20);
  // sigmabedA ~ normal(0,20);

  betaMicroxtemp ~ normal(0,5);
  // betaMicroxmacroA ~ normal(0,20);
  betaMicroxLAI ~ normal(0,5);

  b_year2_ncp ~ normal(0,1);
  b_temp_ncp ~ normal(0,1);
  b_LAI_ncp ~ normal(0,1);
  // b_macroA_ncp ~ normal(0,1);
   // b_depth_ncp ~ normal(0,10);
  // b_bedA_ncp ~ normal(0,10);

  // b_muSt_ncp ~ normal(0,1);
  // b_MacroSite_ncp ~ normal(0,1);
  
}

generated quantities {
} 
