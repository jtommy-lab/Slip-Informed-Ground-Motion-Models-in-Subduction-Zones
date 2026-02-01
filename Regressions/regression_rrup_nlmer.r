library(nlme)
library(lme4)
library(ggplot2)
library(dplyr)
library(data.table)
library(readxl)
library(minpack.lm)
library(caret)
library(here)

# =============================================================================
# 1. PATHS
# =============================================================================


database_path      <- here("Data", "Drapela_database.csv")
site_response_path <- here("Data", "PS21_Fs_DR_database.csv")


path_functions     <- here("Regressions", "functions_regression.R")


output_dir         <- here("Regressions", "Results")


if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message(paste("Carpeta creada:", output_dir))
}


if (file.exists(path_functions)) {
  source(path_functions)
} else {
  stop(paste("ERROR: No se encuentra 'functions_regression.R' en:", path_functions))
}


# =============================================================================
# 2. Model Parameters and Input/Output Files
# =============================================================================

### Set TRUE if performing regional nlmer regression ###
regional_flag = FALSE

# Define output file paths based on regional_flag
suffix <- if(regional_flag) "_regional" else ""

dBe_path       <- file.path(output_dir, paste0('dBe_residual_Rrup', suffix, '.csv'))
dWe_path       <- file.path(output_dir, paste0('dWe_residual_Rrup', suffix, '.csv'))
stats_path     <- file.path(output_dir, paste0('Stats_residual_Rrup', suffix, '.csv'))
coef_Rrup_path <- file.path(output_dir, paste0('coeficientes_Rrup_nlmer', suffix, '.csv'))

# =============================================================================
# 3. Processing and Regression
# =============================================================================

message("Loading database...")
database = read.csv(database_path, check.names = FALSE)
Fs_DR = read.csv(site_response_path, check.names = FALSE)

columns_needed = c('NGAsubEQID','Earthquake_Magnitude','Station_Name','Vs30_Selected_for_Analysis_m_s',"DatabaseRegion")

Periods = c(-1.0,0.0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.5,10.0)

dBe_residual_Rrup = data.frame()
dWe_residual_Rrup = data.frame()
Stats_residual_Rrup = data.frame()
coef_Rrup = data.frame()

for (i in seq_along(Periods)){

    ## Generate colnames array needed for obs and Fs csv
    obs_colname = get_column_name(Periods[i],Fs_Obs = "Obs")
    Fs_colname = get_column_name(Periods[i],Fs_Obs = "Fs")
    columns_needed_obs = copy(columns_needed)
    columns_needed_obs[6] = 'ClstD_km'
    columns_needed_obs[7] = obs_colname
    columns_needed_Fs = copy(columns_needed[-5])
    columns_needed_Fs[5] = Fs_colname

    print(paste0('Processing Period: ', Periods[i], 's | IM = ', obs_colname))

    ## Generate dataframe for Periods in actual iteration 
    Fs_database_T = copy(Fs_DR[,columns_needed_Fs])
    Sa_database_T = copy(database[,columns_needed_obs]) 
    Sa_database_T[,obs_colname] = log(Sa_database_T[,obs_colname])
    Sa_database_T = Sa_database_T[is.na(Sa_database_T[,obs_colname]) == FALSE,]    
    Sa_rock_database_T = merge(Sa_database_T,Fs_database_T,by = columns_needed[-5])
    Sa_rock_database_T[,obs_colname] = Sa_rock_database_T[,obs_colname]-Sa_rock_database_T[,Fs_colname]
    Sa_rock_database_T = subset(Sa_rock_database_T,Sa_rock_database_T$Vs30_Selected_for_Analysis_m_s >= 150)
    Sa_rock_database_T = df_Period_3reg(Sa_rock_database_T)

    ## Calculate predictions
    M = Sa_rock_database_T$Earthquake_Magnitude
    Rrup = Sa_rock_database_T$ClstD_km
    R = get_R(M,Rrup)$R
    EQID = Sa_rock_database_T[,'NGAsubEQID']
    Vs30 = Sa_rock_database_T[,'Vs30_Selected_for_Analysis_m_s']
    Station_Name = Sa_rock_database_T[,'Station_Name']

    regression_df_T = data.frame(log_obs_rock = Sa_rock_database_T[,obs_colname], R = R, M = M, EVENT = EQID,STATION = Station_Name,
                                 REGION = Sa_rock_database_T$DatabaseRegion)

    start_values_ini = c(c0 = 0, c1 = -1,a0 = -0.001,c4 = 1)
    fitparam = nls.control(maxiter = 1000, tol = 1e-05, minFactor = 1/2048)
    
    # Try/Catch para evitar que el loop se rompa si nls falla en un periodo
    tryCatch({
        fit_ini_1 = nls(log_obs_rock ~ GMM_rock_v3(c0,c1,R,M,a0,c4),data = regression_df_T, start = start_values_ini, control=fitparam)
        
        start_values = c(c0 = coef(fit_ini_1)[['c0']], c1 = coef(fit_ini_1)[['c1']],a0 = coef(fit_ini_1)[['a0']],c4 = coef(fit_ini_1)[['c4']])    

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
        dBe_residual_T$M =  regression_df_T %>%
                            distinct(EVENT, M, .keep_all = TRUE) %>%
                            filter(EVENT %in% dBe_residual_T$EVENT) %>%
                            arrange(match(EVENT, dBe_residual_T$EVENT)) %>%
                            pull(M)
        dBe_residual_T$Period = Periods[i]
        dBe_residual_T$dBe = dBe_residual_T$c0
        rownames(dBe_residual_T) = 1:nrow(dBe_residual_T)
        dBe_residual_T = dBe_residual_T[,-1]        
        
        dWe_residual_T = data.frame(NGAsubEQID = EQID, Station_Name = Station_Name,M = M,Rrup = Rrup,R = R,Period = Periods[i],dWe = residuals(fit_1))
       
        dBe_residual_Rrup = rbind(dBe_residual_Rrup,dBe_residual_T)
        dWe_residual_Rrup = rbind(dWe_residual_Rrup,dWe_residual_T)
        rownames(dWe_residual_Rrup) = NULL

        tau_T = attr(VarCorr(fit_1)$EVENT, "stddev")[['c0']]
        phi_T = attr(VarCorr(fit_1),"sc")
        sigma_T = sqrt(tau_T^2 + phi_T^2)
        logLik_T = logLik(fit_1)
        Stats_residual_T = data.frame(Period = Periods[i],sigma = sigma_T, tau = tau_T, phi = phi_T,logLik = logLik_T,AIC=AIC(fit_1),BIC=BIC(fit_1))
        Stats_residual_Rrup = rbind(Stats_residual_Rrup,Stats_residual_T)    

        coef_Rrup_T = data.frame(Period = Periods[i],
                                 c1 = c1_T,
                                 Global_a0 = a0_T,
                                 c4 = c4_T,
                                 Global_c0 = c0_T)
        if (regional_flag == TRUE){
            coef_Rrup_T['a0_Alaska'] = coef_Rrup_T['Global_a0'] + regional_coef['a0','Alaska']
            coef_Rrup_T['a0_CentralAmerica&Mexico'] = coef_Rrup_T['Global_a0'] + regional_coef['a0','CentralAmerica&Mexico']
            coef_Rrup_T['a0_SouthAmerica'] = coef_Rrup_T['Global_a0'] + regional_coef['a0','SouthAmerica']
            coef_Rrup_T['a0_Japan'] = coef_Rrup_T['Global_a0'] + regional_coef['a0','Japan']
        }    
        coef_Rrup = rbind(coef_Rrup,coef_Rrup_T)
        
    }, error = function(e){
        message(paste("Error en periodo:", Periods[i], "-", e$message))
    })
}

message(paste("Saving results in:", output_dir))
write.csv(dBe_residual_Rrup,dBe_path)
write.csv(dWe_residual_Rrup,dWe_path)
write.csv(Stats_residual_Rrup,stats_path)
write.csv(coef_Rrup,coef_Rrup_path)
message("Proceso finalizado.")
