library(readxl)
library(ggplot2)
library(dplyr)
library(abind)
library(openxlsx)
library(optimx)
library(data.table)
library(lme4)
library(patchwork)
library(pracma)
library(mgcv)
library(minpack.lm)
library(here) 

# =============================================================================
# 1. PATHS
# =============================================================================

## Set TRUE if performing Rp nlmer regression with wi = 1, FALSE for wi = Slip_i/sum(Slip_i) ###
unit_flag = FALSE

### Set TRUE if performing regional nlmer regression ###
regional_flag = FALSE

path_functions     <- here("Regressions", "functions_regression.R")
output_dir         <- here("Regressions", "Results")
data_dir           <- here("Data")
rp_data_dir        <- here("Data", "Rp")

# Crear carpeta de resultados si no existe
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source(path_functions)

database_path      <- file.path(data_dir, 'Drapela_database.csv')
site_response_path <- file.path(data_dir, 'PS21_Fs_DR_database.csv')

if (unit_flag == TRUE){
  Rp_path = file.path(rp_data_dir, "Rp_unit_median_values.xlsx")
} else {
  Rp_path = file.path(rp_data_dir, "Rp_median_values.xlsx")
}

suffix_unit <- if(unit_flag) "_unit" else ""
suffix_reg  <- if(regional_flag) "_regional" else ""

# Filenames
coef_filename     <- paste0('coeficientes_Rp_nlmer', suffix_reg, suffix_unit, '.csv')
p_values_filename <- paste0('p_selected_Rp', suffix_reg, suffix_unit, '.csv')
dBe_filename      <- paste0('dBe_residual_Rp', suffix_reg, suffix_unit, '.csv')
dWe_filename      <- paste0('dWe_residual_Rp', suffix_reg, suffix_unit, '.csv')
stats_filename    <- paste0('Stats_residual_Rp', suffix_reg, suffix_unit, '.csv')

#
coef_Rp_path  <- file.path(output_dir, coef_filename) 
p_values_path <- file.path(output_dir, p_values_filename)
dBe_path      <- file.path(output_dir, dBe_filename)
dWe_path      <- file.path(output_dir, dWe_filename)
stats_path    <- file.path(output_dir, stats_filename)

##### Checking paths ######
message("---- Checking paths ----")
message(paste("Input Coefficients:", coef_Rp_path))
message(paste("Input Rp Data:", Rp_path))
message(paste("Output Stats:", stats_path))


# =============================================================================
# 2. COMPUTE OPTIMAL EXPONENT P
# =============================================================================

database = read.csv(database_path, check.names = FALSE)
Fs_DR = read.csv(site_response_path, check.names = FALSE)
coeff_Rp = read.csv(coef_Rp_path, check.names = FALSE)
coeff_Rp = coeff_Rp[,-1] 

p_opt_values = data.frame()

for (T in unique(coeff_Rp[['Period']])){
    coeff_Rp_T = subset(coeff_Rp,coeff_Rp['Period'] == T)
    AIC_T = -2*coeff_Rp_T$logLik + 2*3
    AIC_min_T = min(AIC_T)
    AIC_diff_T = AIC_T - AIC_min_T
    p_min_T = coeff_Rp_T$p_value[which(AIC_diff_T <= 2)]
    AIC_opt_T = AIC_diff_T[which(AIC_diff_T <= 2)]
    weight = exp(-0.5*AIC_opt_T)/sum(exp(-0.5*AIC_opt_T))
    p_opt_T = sum(p_min_T*weight)
    p_weighted_sd_T = sqrt(sum(weight * (p_min_T - p_opt_T)^2))
    if (p_opt_T == 0){
        p_opt_T = 0.5
    }    
    p_opt_values = rbind(p_opt_values,data.frame(Period = T, p_opt = p_opt_T,p_sd = p_weighted_sd_T))
}

p_selected_Rp_curve = subset(p_opt_values, p_opt_values$Period > 0 )

#### If any p_sd == 0, to avoid problems with regression ###
min_sd_nonzero <- min(p_selected_Rp_curve$p_sd[p_selected_Rp_curve$p_sd > 0])
p_selected_Rp_curve$p_sd[p_selected_Rp_curve$p_sd == 0] <- min_sd_nonzero * 0.1

# p ~ A + (B*T)/(C+T)
tryCatch({
    mod_nls <- nlsLM(
      p_opt ~  A + (B * Period) / (C + Period),
      data  = p_selected_Rp_curve,
      start = list(A = mean(p_selected_Rp_curve$p_opt), B = -1,C = 3),
      weights = 1 / (p_sd^2),
      control = nls.lm.control(maxiter = 1000)
    )
    
    pred <- predict(mod_nls, newdata = data.frame(Period = p_opt_values$Period))
    p_opt_values$p_opt_pred = pred
    
    
    if(any(p_opt_values$Period == 3)) p_opt_values[p_opt_values$Period == 3, "p_opt_pred"] = subset(p_opt_values, Period == 3)$p_opt_pred
    if(any(p_opt_values$Period == 0.01)) p_opt_values[p_opt_values$Period == 0.01, "p_opt_pred"] = subset(p_opt_values, Period == 0.01)$p_opt_pred

}, error = function(e) {
    message("Advertencia: Falló el ajuste nlsLM para p_values. Usando valores crudos.")
    p_opt_values$p_opt_pred = p_opt_values$p_opt
})

write.csv(p_opt_values, p_values_path, row.names = FALSE)


p_opt_values = read.csv(p_values_path, check.names = FALSE)

# =============================================================================
# 3. COMPUTE RESIDUALS FOR OPTIMAL P EXPONENTS
# =============================================================================

columns_needed = c('NGAsubEQID','Earthquake_Magnitude','Station_Name','Vs30_Selected_for_Analysis_m_s',"DatabaseRegion")
Periods = c(-1.0,0.0,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.5,10.0)

dBe_residual_Rp = data.frame()
dWe_residual_Rp = data.frame()
Stats_residual_Rp = data.frame()

for (i in seq_along(Periods)){
    
    # Get p_value for maximum logLike model
    p_T_Rp = subset(p_opt_values,p_opt_values$Period == Periods[i])$p_opt_pred
    p_T_Rp = round(p_T_Rp/0.5)*0.5
    if (p_T_Rp == 0){
        p_T_Rp = -0.5
    }

    message(paste0("Procesando Periodo: ", Periods[i], " | p seleccionado: ", p_T_Rp))

    ## Generate colnames array needed for obs and Fs csv
    obs_colname = get_column_name(Periods[i],Fs_Obs = "Obs")
    Fs_colname = get_column_name(Periods[i],Fs_Obs = "Fs")
    columns_needed_obs = copy(columns_needed)    
    columns_needed_obs[6] = obs_colname
    columns_needed_Fs = copy(columns_needed[-5])
    columns_needed_Fs[5] = Fs_colname    

    ## Generate dataframe for Periods in actual iteration 
    Fs_database_T = copy(Fs_DR[,columns_needed_Fs])
    Sa_database_T = copy(database[,columns_needed_obs]) 
    Sa_database_T[,obs_colname] = log(Sa_database_T[,obs_colname])
    Sa_database_T = Sa_database_T[is.na(Sa_database_T[,obs_colname]) == FALSE,]    
    Sa_rock_database_T = merge(Sa_database_T,Fs_database_T,by = columns_needed[-5])
    Sa_rock_database_T$log_Sa_T_rock = Sa_rock_database_T[,obs_colname]-Sa_rock_database_T[,Fs_colname]    
    Sa_rock_database_T = subset(Sa_rock_database_T,Sa_rock_database_T$Vs30_Selected_for_Analysis_m_s >= 150)
    Sa_rock_database_T = df_Period_3reg(Sa_rock_database_T)

    sheet_name = paste0("p = ", p_T_Rp)
    tryCatch({
        Rp_median = suppressMessages(read_excel(Rp_path, sheet = sheet_name))
        
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
        } else {
          fit_1 <- nlmer(log_obs_rock ~ GMM_rock_v3(c0, c1, R, M, a0, c4) ~ c0 + c1 + a0 + c4 + (c0+0|EVENT),data = regression_df_T,start = start_values,nAGQ = 1L)
        }
    
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
    
        tau_T = attr(VarCorr(fit_1)$EVENT, "stddev")[['c0']]
        phi_T = attr(VarCorr(fit_1),"sc")
        sigma_T = sqrt(tau_T^2 + phi_T^2)
        Stats_residual_T = data.frame(Period = Periods[i],sigma = sigma_T, tau = tau_T, phi = phi_T)
        Stats_residual_Rp = rbind(Stats_residual_Rp,Stats_residual_T)     
        
        dWe_residual_T = data.frame(NGAsubEQID = EQID, Station_Name = Station_Name,M = M,Rp = Rp,Period = Periods[i],dWe = residuals(fit_1))    
    
        dBe_residual_Rp = rbind(dBe_residual_Rp,dBe_residual_T)
        dWe_residual_Rp = rbind(dWe_residual_Rp,dWe_residual_T)    

    }, error = function(e) {
        message(paste("Error en Periodo", Periods[i], ":", e$message))
    })
}


write.csv(dBe_residual_Rp,dBe_path)
write.csv(dWe_residual_Rp,dWe_path)
write.csv(Stats_residual_Rp,stats_path)
