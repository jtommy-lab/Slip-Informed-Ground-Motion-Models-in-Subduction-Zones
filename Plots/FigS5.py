import os
import sys
import numpy as np
import pandas as pd
import pygmt
import matplotlib.pyplot as plt
import xarray as xr
import random
import openquake.hazardlib.imt as IMT
from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SInter
from openquake.hazardlib.gsim.montalva_2017 import MontalvaEtAl2017SInter

# =============================================================================
# 1. PATH
# =============================================================================

try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.getcwd()

repo_root = os.path.dirname(script_dir)
data_dir = os.path.join(repo_root, 'Data')


rp_dir = os.path.join(data_dir, 'Rp')
models_dir = os.path.join(data_dir, 'Modelos_cosismicos')
locking_dir = os.path.join(data_dir, 'Modelos_locking')


results_dir = os.path.join(repo_root, 'Regressions', 'Results')

func_dir = os.path.join(repo_root, 'functions')
output_dir = script_dir 

if func_dir not in sys.path:
    sys.path.append(func_dir)
import functions_py

np.random.seed(42)
random.seed(42)

print(f"Directorio Data: {data_dir}")
print(f"Guardando figuras en: {output_dir}")

# =============================================================================
# 2. LOAD FILES
# =============================================================================

eqid = '4000068' # Tokachi


database = pd.read_csv(os.path.join(data_dir, 'Drapela_database.csv'), index_col=None)
Fs = pd.read_csv(os.path.join(data_dir, 'PS21_Fs_DR_database.csv'), index_col=None)


coeff_Rrup = pd.read_csv(os.path.join(results_dir, 'coeficientes_Rrup_nlmer.csv'), index_col=False)
coeff_Rp = pd.read_csv(os.path.join(results_dir, 'coeficientes_Rp_nlmer.csv'), index_col=False)


trench_tokachi = pd.read_csv(os.path.join(data_dir, 'trench-tohoku3.xyz'), sep='\s+', index_col=None, names=['Longitude', 'Latitude'])


modelos_cosismicos = pd.read_csv(os.path.join(models_dir, 'modelos_cosismicos.csv'), index_col=False)
modelos_cosismicos = modelos_cosismicos.loc[modelos_cosismicos['NGAsubEQID'] == eqid]


slab_tokachi = pygmt.load_dataarray(os.path.join(data_dir, 'kur_slab2_dep_02.24.18.grd'))
lock_tokachi_LL2016_grd = pygmt.load_dataarray(os.path.join(locking_dir, 'lock_tokachi_LL2016.grd'))


mask = np.isnan(slab_tokachi)
lock_tokachi_LL2016_grd = lock_tokachi_LL2016_grd.where(~mask)


Rp_sheets = pd.read_excel(os.path.join(rp_dir, 'Rp_median_values.xlsx'), sheet_name=None)
Rp_lock_sheets = pd.read_excel(os.path.join(rp_dir, 'Rp_lock_median_values.xlsx'), sheet_name=None)
Rp_lock_sheets_tokachi = pd.read_excel(os.path.join(rp_dir, '4000068_Rp_lock_values.xlsx'), sheet_name=None)


p_opt_Rp = pd.read_csv(os.path.join(results_dir, 'p_selected_Rp.csv'), index_col=None)


topo_grid = os.path.join(data_dir, 'tokachi.nc')
topo_gradient = os.path.join(data_dir, 'tokachi.int')
cpt_gray = os.path.join(data_dir, 'grayscale02.cpt')


# =============================================================================
# 3. OPENQUAKE CONFIG
# =============================================================================

gmpe_parker = ParkerEtAl2020SInter(region='JP', saturation_region='JP_Pac')
gmpe_montalva = MontalvaEtAl2017SInter()

tokachi_database = database.loc[database['NGAsubEQID'] == eqid]

R_gmm_Rrup = np.logspace(np.log10(40), np.log10(700), 100)
R_gmm_Rp = np.logspace(np.log10(60), np.log10(900), 100)
Rref_gmm = functions_py.get_R_Rref_interface(tokachi_database['Earthquake_Magnitude'].unique(), R_gmm_Rp)[1]

input_df = pd.DataFrame({
    "mag": len(R_gmm_Rrup) * [tokachi_database['Earthquake_Magnitude'].unique()[0]],
    "rrup": R_gmm_Rp,
    "vs30": len(R_gmm_Rp) * [760],
    "region": len(R_gmm_Rp) * ['JP'],
    "saturation_region": len(R_gmm_Rp) * ['JP_Pac'],
    "backarc": len(R_gmm_Rp) * [False]
})

n = input_df.shape[0]
ctx = np.recarray(n, dtype=np.dtype([("mag", float), ("rrup", float), ("vs30", float), ("region", str), ("saturation_region", str), ("backarc", bool)]))
ctx["rrup"] = input_df["rrup"]
ctx["mag"] = input_df["mag"]
ctx["vs30"] = input_df["vs30"]
ctx["region"] = input_df["region"]
ctx["saturation_region"] = input_df["saturation_region"]

# =============================================================================
# 4. PLOTS
# =============================================================================

periods_plot = [0]
R_gmm_v2 = functions_py.get_R_Rref_interface(tokachi_database['Earthquake_Magnitude'].unique(), R_gmm_Rrup)[0]

for T in periods_plot:
    print(f"Procesando Fig Tokachi para T={T}")
    
    if T == -1:
        column_obs, column_Fs = 'PGV_cm_sec', 'Fs - PGV_cm_sec'
        imt = IMT.PGV()
    elif T == 0:
        column_obs, column_Fs = 'PGA_g', 'Fs - PGA_g'
        imt = IMT.PGA()
    else:
        column_obs, column_Fs = f'T = {float(T)}', f'Fs - T = {float(T)}'
        imt = IMT.SA(T)

    
    raw_p = p_opt_Rp['p_opt_pred'].loc[p_opt_Rp['Period'] == T].values[0]
    p_opt_T = round(raw_p * 2) / 2
    print(f"  > P-value seleccionado: {p_opt_T}")

    
    coeff_Rp_T = coeff_Rp.loc[(coeff_Rp['Period'] == T) & (coeff_Rp['p_value'] == p_opt_T)]
    coeff_Rrup_T = coeff_Rrup.loc[coeff_Rrup['Period'] == T]

    
    tokachi_df = functions_py.get_event_period_df(
        database, Fs, eqid, T, coeff_Rp, coeff_Rrup, p_value=p_opt_T,
        Rp_sheets=Rp_sheets, Rp_sheets_lock=Rp_lock_sheets, region='Global', Rref=Rref_gmm
    )
    tokachi_df = tokachi_df.dropna(subset='log_obs_rock')
    obs_T = tokachi_df['log_obs_rock']
    mag = np.array(len(R_gmm_Rrup) * [tokachi_df['Earthquake_Magnitude'].unique()[0]])
    
    
    pred_Rrup = functions_py.DR_pred_GMM_Rp(R_gmm_v2, mag, coeff_Rrup_T, Rref=Rref_gmm)
    pred_Rp = functions_py.DR_pred_GMM_Rp(R_gmm_Rp, mag, coeff_Rp_T, Rref=Rref_gmm)

    Rp = tokachi_df['Rp']
    
    Rp_lock_tokachi_sheet = Rp_lock_sheets_tokachi['p = ' + str(p_opt_T)]
    Rp_lock_LL2016 = Rp_lock_tokachi_sheet.loc[:, Rp_lock_tokachi_sheet.columns.str.startswith('s2003')].median(axis=1).values
    
    R, Rrup = tokachi_df['R'], tokachi_df['Rrup']

    
    ln_mean_parker, sig_parker, tau_parker, phi_parker = [np.zeros([1, n]) for _ in range(4)]
    ln_mean_montalva, sig_montalva, tau_montalva, phi_montalva = [np.zeros([1, n]) for _ in range(4)]
    gmpe_parker.compute(ctx, [imt], ln_mean_parker, sig_parker, tau_parker, phi_parker)
    gmpe_montalva.compute(ctx, [imt], ln_mean_montalva, sig_montalva, tau_montalva, phi_montalva)
    ln_mean_parker = ln_mean_parker.reshape((len(mag),))
    ln_mean_montalva = ln_mean_montalva.reshape((len(mag),))

    fig1 = pygmt.Figure()
    projection_log = "X?l/Y?"
    projection_geo = "M5.5c"
    colores = [f"#{random.randint(0, 255):02x}{random.randint(0, 255):02x}{random.randint(0, 255):02x}" for _ in range(len(modelos_cosismicos))]
    cont = 0

   
    with fig1.subplot(nrows=1, ncols=1, figsize=("4c", "12c"), sharex=True):
        with fig1.set_panel(panel=[0, 0]):
            fig1.basemap(region=[139, 146, 37, 44], projection=projection_geo, frame=['WSne'])
            pygmt.makecpt(cmap=cpt_gray, series=[-8000, 3000, 10])
            
            
            fig1.grdimage(grid=topo_grid, shading=topo_gradient, projection=projection_geo, cmap=True)
            
            pygmt.makecpt(cmap='hot', series=[0, 1], reverse=True)
            fig1.grdimage(grid=lock_tokachi_LL2016_grd, projection=projection_geo, cmap=True, nan_transparent=True)
            
            fig1.plot(x=tokachi_df['Station_Longitude_deg'], y=tokachi_df['Station_Latitude_deg'], fill='green', style='t0.2c', projection=projection_geo, pen="0.1p,black")
            fig1.plot(x=trench_tokachi['Longitude'], y=trench_tokachi['Latitude'], style="f0.5i/0.1i+r+t", fill='black', pen="1.25p,black", projection=projection_geo)
            fig1.coast(shorelines=True, projection=projection_geo)
            fig1.text(text="a) LM16", position='TL', font="8p,Helvetica,black", projection=projection_geo, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            
            pygmt.config(FONT="10p,Helvetica,black")
            fig1.colorbar(frame="xa.25f.25+lLocking degree", position="g139/35.7+w5.5c/0.4c+h", projection=projection_geo, region=[139, 146, 37, 44])
            pygmt.config(FONT="8p,Helvetica,black")

            for i in modelos_cosismicos['event_tag']:
                model_path = os.path.join(models_dir, str(i) + '.xyz')
                model_i = pd.read_csv(model_path, sep='\s+', index_col=None, names=['Longitude', 'Latitude', 'Slip', 'Depth'])
                model_i = model_i[['Longitude', 'Latitude', 'Slip']]
                model_i['Slip'] = model_i['Slip'] / model_i['Slip'].max()
                
                resolution = np.sqrt((model_i['Longitude'].iloc[1] - model_i['Longitude'].iloc[0])**2 + (model_i['Latitude'].iloc[1] - model_i['Latitude'].iloc[0])**2)
                resolution_m = 3 * np.ceil(resolution * 60)
                region = [min(model_i['Longitude']) - 1, max(model_i['Longitude']) + 1, min(model_i['Latitude']) - 1, max(model_i['Latitude']) + 1]
                model_limits = [min(model_i['Longitude']), max(model_i['Longitude']), min(model_i['Latitude']), max(model_i['Latitude'])]
                model_grid = pygmt.nearneighbor(data=model_i, spacing="1m", region=model_limits, search_radius=str(resolution_m) + 'm', sectors="6")
                
                slab_model_grid = slab_tokachi.interp(y=model_grid.lat, x=model_grid.lon)
                model_grid = model_grid.where(slab_model_grid.notnull())
                fig1.grdcontour(grid=model_grid, limit=[0, 0.8], projection=projection_geo, annotation="n", levels=0.8, pen='1p,' + colores[cont])
                cont += 1

    
    fig1.shift_origin(xshift="w+3.8c", yshift="-2.8c")
    with fig1.subplot(nrows=3, ncols=1, figsize=("4c", "12c"), sharey=True):
        y_min = min(min(obs_T), min(pred_Rp), min(ln_mean_parker), min(ln_mean_montalva)) - 0.5
        y_max = max(max(obs_T), max(pred_Rp), max(ln_mean_parker), max(ln_mean_montalva)) + 1.8
        
        
        with fig1.set_panel(panel=[0, 0]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['Wsne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rp_lock_LL2016, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)
            
            fig1.plot(x=R_gmm_Rp, y=pred_Rp, pen='1p,red', label='This model')
            fig1.text(text="b) R@-p,lock@-", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.5c", box="+gwhite+p1p")

        
        with fig1.set_panel(panel=[1, 0]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['Wsne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rp, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)
            
            fig1.plot(x=R_gmm_Rp, y=pred_Rp, pen='1p,red', label='This model')
            fig1.text(text="c) R@-p@-", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.5c", box="+gwhite+p1p")

        
        with fig1.set_panel(panel=[2, 0]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['WSne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rrup, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)
            
            fig1.plot(x=R_gmm_Rrup, y=pred_Rrup, pen='1p,red', label='This model')
            fig1.plot(x=R_gmm_Rrup, y=ln_mean_parker, pen='1p,blue', label='Parker et al. (2021)')
            
            fig1.text(text="PGA at V@-s30@- = 760 m/s [log units]", position='ML', font="10p,Helvetica,black", projection=projection_log, offset='-1.1/0.8', no_clip=True, angle=90)
            fig1.text(text="d) R@-rup@-", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w3.5c", box="+gwhite+p1p")
            fig1.text(text="Distance [km]", position='BC', font="10p,Helvetica,black", projection=projection_log, offset='0/-0.9', no_clip=True)

    
    fig1.savefig(os.path.join(output_dir, 'FigS5.pdf'))

print("✅ Proceso completado. Figuras guardadas en 'Plots/'.")