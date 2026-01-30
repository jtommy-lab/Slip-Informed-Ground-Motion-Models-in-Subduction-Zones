import os
import sys
import numpy as np
import pandas as pd
import pygmt

# =============================================================================
# 1. PATH
# =============================================================================


script_dir = os.path.dirname(os.path.abspath(__file__))



repo_root = os.path.dirname(script_dir)


data_dir = os.path.join(repo_root, 'Data')
rp_dir = os.path.join(data_dir, 'Rp')
functions_path = os.path.join(repo_root, 'functions')
output_folder = os.path.join(repo_root, "Figuras_paper")


# Configurar salida
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


sys.path.append(functions_path)

import functions_py



# =============================================================================
# 2. LOAD DATA
# =============================================================================


db_path = os.path.join(data_dir, 'Drapela_database.csv')
fs_path = os.path.join(data_dir, 'PS21_Fs_DR_database.csv')
rp_excel_path = os.path.join(rp_dir, 'Rp_median_values.xlsx')


data = pd.read_csv(db_path, index_col=False)
Fs = pd.read_csv(fs_path, index_col=False)


Rp_min_df = pd.read_excel(rp_excel_path, sheet_name='p = -50', engine='openpyxl')
Rp_min_df = Rp_min_df.rename(columns={'Rp_median_km': 'Rp_median_km_minp'})

Rp_max_df = pd.read_excel(rp_excel_path, sheet_name='p = 50', engine='openpyxl')
Rp_max_df = Rp_max_df.rename(columns={'Rp_median_km': 'Rp_median_km_maxp'})


Rp_df = Rp_min_df.merge(Rp_max_df, how='inner')
data = data.merge(Rp_df, how='inner')

# =============================================================================
# 3. Graphics
# =============================================================================

p_values = [-9]
distance = np.logspace(0, 3.2, 3000)

projection = "X?/Y?"
projection_log = "X?l/Y?"

periods_plot = [0]  # [-1, 0, 0.5, 0.75, 3, 10]

for cont, T in enumerate(periods_plot):    

    if T == -1:
        column_obs = 'PGV_cm_sec'
        column_Fs = 'Fs - PGV_cm_sec'
    elif T == 0:
        column_obs = 'PGA_g'
        column_Fs = 'Fs - PGA_g'
    else:
        column_obs = 'T = ' + str(float(T))
        column_Fs = 'Fs - T = ' + str(float(T))
    
    p_T = p_values[cont]
    
    if float(p_T).is_integer():
        p_T = int(p_T)
        
    print(f"Procesando Periodo: {T} | p-value: {p_T}")

    
    sheet_name = 'p = ' + str(p_T)
    Rp_df_T = pd.read_excel(rp_excel_path, sheet_name=sheet_name, engine='openpyxl')
    
    
    data_T = data.copy()
    data_T = data_T.merge(Rp_df_T, how='inner')
    
    obs_Fs_df = pd.merge(data_T, Fs, how="inner", 
                         on=['NGAsubEQID', 'Earthquake_Magnitude', 'Station_Name', 'Vs30_Selected_for_Analysis_m_s'])
    
    Rrup = obs_Fs_df['ClstD_km']
    Rasp = obs_Fs_df['Rasp_median_km']
    Rehd = obs_Fs_df['Rehd_median_km']
    Rp_T = obs_Fs_df['Rp_median_km']
    M = round(obs_Fs_df['Earthquake_Magnitude'], 1)
    
    # --- AQUÍ ESTÁ EL CAMBIO PRINCIPAL ---
    R, _ = functions_py.get_R_Rref_interface(M, Rrup)
    # -------------------------------------
    
    obs_Fs_df_T = pd.DataFrame({
        'Rrup': Rrup,
        'Distance': R,
        'Rasp': Rasp,
        'Rehd': Rehd,
        'Rp': Rp_T,
        'Magnitude': M,
        'logObs_rock': np.log(obs_Fs_df[column_obs]) - obs_Fs_df[column_Fs]
    })
    
    
    fig1 = pygmt.Figure()
    pygmt.makecpt(cmap="jet", series=[6.6, 9.4, 0.4])
    pygmt.config(FONT_TITLE="14p,Helvetica,black")
    
   
    min_lim_obs = min(obs_Fs_df_T['logObs_rock']) - np.abs(min(obs_Fs_df_T['logObs_rock']) * 0.1)
    max_lim_obs = max(obs_Fs_df_T['logObs_rock']) + np.abs(max(obs_Fs_df_T['logObs_rock']) * 0.1)
    mid_y_point = max_lim_obs 
    mid_x_point = 10**(np.log10(10) + (np.log10(1400) - np.log10(10)) / 2)
    
    
    with fig1.subplot(nrows=2, ncols=1, figsize=("12c", "10c"), sharex=True):
        
        
        with fig1.set_panel(panel=[0, 0]):
            fig1.basemap(region=[10, 1400, min_lim_obs, max_lim_obs], projection=projection_log, 
                         frame=['Wsne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=obs_Fs_df_T['Distance'], y=obs_Fs_df_T['logObs_rock'], 
                      style='c0.1c', fill=obs_Fs_df_T['Magnitude'], cmap=True, pen='0.1p,black')
            fig1.text(text="a) R@-RUP@-", position='TL', font="10p,Helvetica,black", 
                      angle=0, pen="0.25p,black,solid", fill='white', offset='0.1/-0.1')
        
        
        with fig1.set_panel(panel=[1, 0]):
            fig1.basemap(region=[10, 1400, min_lim_obs, max_lim_obs], projection=projection_log, 
                         frame=['WSne', 'xa2fg3+lDistance [km]', 'ya2fg3'])
            fig1.plot(x=obs_Fs_df_T['Rp'], y=obs_Fs_df_T['logObs_rock'], 
                      style='c0.1c', fill=obs_Fs_df_T['Magnitude'], cmap=True, pen='0.1p,black')
            fig1.text(text="b) R@-P@-", position='TL', font="10p,Helvetica,black", 
                      angle=0, pen="0.25p,black,solid", fill='white', offset='0.1/-0.1')
            
            fig1.text(text=column_obs + " at V@-s30@- = 760 m/s [log units]", 
                      x=-110, y=mid_y_point, font="12p,Helvetica,black", 
                      projection=projection, no_clip=True, angle=90)
            
            fig1.colorbar(frame=["xa.4f.1+lMoment magnitude"], position="JMR+o0.5c/3c+w8c")

    
    filename = "Fig3.pdf"
    save_path = os.path.join(output_folder, filename)
    fig1.savefig(save_path, dpi=300)
    print(f"Figura guardada: {save_path}")

