library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)
library(data.table)
library(readxl)
library(minpack.lm)
library(caret)


## Set TRUE if performing Rp nlmer regression with wi = 1, FALSE for wi = Slip_i/sum(Slip_i) ###

unit_flag = TRUE

### Set TRUE if performing regional nlmer regression (a0 correction for subduction zone name), FALSE for global nlmer regression ###

regional_flag = FALSE

###############################################
'Script to perform mixed effects regression for Rp-based GMM'
'Following the methodology of Stafford (2014)'
###############################################

path_functions = '/home/jtommy/Escritorio/Respaldo/functions/functions_regression.R'
source(path_functions)

# Load databases paths #
database_path = '/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/Drapela_database.csv'
site_response_path = '/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/PS21_Fs_DR_database.csv'



if (unit_flag == TRUE){
  Rp_path = "/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/Rp_unit_median_values.xlsx"
} else {
  Rp_path = "/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/Rp_median_values.xlsx"
}

# Output paths #
if (unit_flag == TRUE){
  if (regional_flag == TRUE){
    coef_Rp_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer_regional_unit.csv'
    events_performance_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/EVENTS_performance_nlmer_regional_Rp_unit.csv'
  } else {
    coef_Rp_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer_unit.csv'
    events_performance_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/EVENTS_performance_nlmer_Rp_unit.csv'
  }

} else {
  if (regional_flag == TRUE){
    coef_Rp_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer_regional.csv'
    events_performance_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/EVENTS_performance_nlmer_regional_Rp.csv'
    
  } else {
    coef_Rp_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer.csv'
    events_performance_path = '/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/EVENTS_performance_nlmer_Rp.csv'
    
  }
  
}

database = read.csv(database_path, check.names = FALSE)
Fs_DR = read.csv(site_response_path, check.names = FALSE)

columns_needed = c('NGAsubEQID','Earthquake_Magnitude','Station_Name','Vs30_Selected_for_Analysis_m_s',"DatabaseRegion")

p_values = seq(-50,50,0.5)
p_values = p_values[p_values!=0]

Periods = c(-1.0,0.0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.5,10.0)

coef_Rp = data.frame()
events_performance = data.frame()

for (i in seq_along(Periods)){
  for (p_value in p_values){

    ## Generate colnames array needed for obs and Fs csv
    obs_colname = get_column_name(Periods[i],Fs_Obs = "Obs")
    Fs_colname = get_column_name(Periods[i],Fs_Obs = "Fs")
    columns_needed_obs = copy(columns_needed)    
    columns_needed_obs[6] = obs_colname
    columns_needed_Fs = copy(columns_needed[-5])
    columns_needed_Fs[5] = Fs_colname

    print(paste0('p = ',p_value,' Intensity measure = ',obs_colname))

    ## Generate dataframe for Periods in actual iteration 

    Fs_database_T = copy(Fs_DR[,columns_needed_Fs])
    Sa_database_T = copy(database[,columns_needed_obs]) 
    Sa_database_T[,obs_colname] = log(Sa_database_T[,obs_colname])
    Sa_database_T = Sa_database_T[is.na(Sa_database_T[,obs_colname]) == FALSE,]    
    Sa_rock_database_T = merge(Sa_database_T,Fs_database_T,by = columns_needed[-5])
    Sa_rock_database_T$log_Sa_T_rock = Sa_rock_database_T[,obs_colname]-Sa_rock_database_T[,Fs_colname]    
    Sa_rock_database_T = subset(Sa_rock_database_T,Sa_rock_database_T$Vs30_Selected_for_Analysis_m_s >= 150)
    Sa_rock_database_T = df_Period_3reg(Sa_rock_database_T)

    Rp_median = read_excel(Rp_path, sheet = paste0("p = ",p_value))
    Rp_median_df = Rp_median[c("NGAsubEQID","Station_Name","Rp_median_km")]
    Sa_rock_database_T = merge(Sa_rock_database_T,Rp_median_df,by = c("NGAsubEQID","Station_Name")) 
    Rp = Sa_rock_database_T$Rp_median_km
    M = Sa_rock_database_T$Earthquake_Magnitude
    EQID = Sa_rock_database_T$NGAsubEQID
    Station_Name = Sa_rock_database_T$Station_Name
    Vs30 = Sa_rock_database_T$Vs30_Selected_for_Analysis_m_s

    regression_df_T = data.frame(log_obs_rock = Sa_rock_database_T$log_Sa_T_rock, R = Rp, M = M, EVENT = EQID,STATION = Station_Name,
                                REGION = Sa_rock_database_T$DatabaseRegion)

    start_values_ini = c(c0 = 0, c1 = -1,a0 = -0.001,c4 = 1)
    
    fitparam = nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048)
    fit_ini_1 = nls(log_obs_rock ~ GMM_rock_v3(c0,c1,R,M,a0,c4),data = regression_df_T, start = start_values_ini, control=fitparam,)
    

    start_values = c(c0 = coef(fit_ini_1)[['c0']], c1 = coef(fit_ini_1)[['c1']],a0 = coef(fit_ini_1)[['a0']],c4 = coef(fit_ini_1)[['c4']])

    # Mixed effects model fitting (Stafford 2014)

    if (regional_flag == TRUE){
      fit_1 <- nlmer(log_obs_rock ~ GMM_rock_v3(c0, c1, R, M, a0, c4) ~ c0 + c1 + a0 + c4 + (c0+0|EVENT) + (a0+0|REGION),data = regression_df_T,start = start_values,nAGQ = 1L)
      regional_coef = as.data.frame(t(ranef(fit_1)$REGION))
    } else {
    fit_1 <- nlmer(log_obs_rock ~ GMM_rock_v3(c0, c1, R, M, a0, c4) ~ c0 + c1 + a0 + c4 + (c0+0|EVENT),data = regression_df_T,start = start_values,nAGQ = 1L)
    }

    c0_T = fixef(fit_1)[['c0']]
    c1_T = fixef(fit_1)[['c1']]
    a0_T = fixef(fit_1)[['a0']]
    c4_T = fixef(fit_1)[['c4']]    
    tau_T = attr(VarCorr(fit_1)$EVENT, "stddev")[['c0']]
    phi_T = attr(VarCorr(fit_1),"sc")
    sigma_T = sqrt(tau_T^2 + phi_T^2)
    
    
    dBe_residual_T = ranef(fit_1)$EVENT
    dBe_residual_T$EVENT = rownames(dBe_residual_T)
    dBe_residual_T$Period = Periods[i]
    dBe_residual_T$p_value = p_value
    dBe_residual_T$dBe = dBe_residual_T$c0
    rownames(dBe_residual_T) = 1:nrow(dBe_residual_T)
    dBe_residual_T = dBe_residual_T[,-1]

    events_performance = rbind(events_performance,dBe_residual_T)
    

    coef_Rp_p_T = data.frame(Period = Periods[i],
                             c1 = c1_T,
                             Global_a0 = a0_T,
                             c4 = c4_T,
                             c0 = c0_T,
                             p_value = p_value,
                             logLik = logLik(fit_1)[1],
                             sigma = sigma_T,
                             tau = tau_T,
                             phi = phi_T
                             )
    if (regional_flag == TRUE){
      coef_Rp_p_T['a0_Alaska'] = coef_Rp_p_T['Global_a0'] + regional_coef['a0','Alaska']
      coef_Rp_p_T['a0_CentralAmerica&Mexico'] = coef_Rp_p_T['Global_a0'] + regional_coef['a0','CentralAmerica&Mexico']
      coef_Rp_p_T['a0_SouthAmerica'] = coef_Rp_p_T['Global_a0'] + regional_coef['a0','SouthAmerica']
      coef_Rp_p_T['a0_Japan'] = coef_Rp_p_T['Global_a0'] + regional_coef['a0','Japan']      
    } 
    coef_Rp = rbind(coef_Rp,coef_Rp_p_T)
  }
}

## Check path before save file

events_performance_sort = events_performance %>% arrange(EVENT)
write.csv(coef_Rp,coef_Rp_path)
write.csv(events_performance_sort,events_performance_path)

