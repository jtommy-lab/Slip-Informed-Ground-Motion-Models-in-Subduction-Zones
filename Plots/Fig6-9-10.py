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
# 1. PATHS
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
# -------------------------



results_dir = os.path.join(repo_root, 'Regressions', 'Results')
# ---------------------------

func_dir = os.path.join(repo_root, 'functions')
output_dir = script_dir 

if func_dir not in sys.path:
    sys.path.append(func_dir)
import functions_py

np.random.seed(42)
random.seed(42)

# =============================================================================
# 2. LOAD DATA
# =============================================================================

eqid = '4000001'


database = pd.read_csv(os.path.join(data_dir, 'Drapela_database.csv'), index_col=None)
Fs = pd.read_csv(os.path.join(data_dir, 'PS21_Fs_DR_database.csv'), index_col=None)


coeff_Rrup = pd.read_csv(os.path.join(results_dir, 'coeficientes_Rrup_nlmer.csv'), index_col=False)
coeff_Rp = pd.read_csv(os.path.join(results_dir, 'coeficientes_Rp_nlmer.csv'), index_col=False)

p_opt_Rp = pd.read_csv(os.path.join(results_dir, 'p_selected_Rp.csv'), index_col=None)


trench_tohoku = pd.read_csv(os.path.join(data_dir, 'trench-tohoku2.xyz'), delim_whitespace=True, names=['Longitude', 'Latitude'])


modelos_cosismicos = pd.read_csv(os.path.join(models_dir, 'modelos_cosismicos.csv'), index_col=False)

database = database[database['NGAsubEQID'].isin(modelos_cosismicos['NGAsubEQID'].unique())]
modelos_cosismicos = modelos_cosismicos[modelos_cosismicos['NGAsubEQID'].isin(database['NGAsubEQID'].unique())]
modelos_cosismicos = modelos_cosismicos.loc[modelos_cosismicos['NGAsubEQID'] == eqid]


lock_tohoku_LL2016_grd = pygmt.load_dataarray(os.path.join(locking_dir, 'lock_tohoku_LL2016.grd'))
lock_tohoku_LL2011_grd = pygmt.load_dataarray(os.path.join(locking_dir, 'lock_japan_LL2011.grd'))
slab_tohoku = pygmt.load_dataarray(os.path.join(data_dir, 'kur_slab2_dep_02.24.18.grd'))

mask = np.isnan(slab_tohoku)
lock_tohoku_LL2011_grd = lock_tohoku_LL2011_grd.where(~mask)
lock_tohoku_LL2016_grd = lock_tohoku_LL2016_grd.where(~mask)


Rp_sheets = pd.read_excel(os.path.join(rp_dir, 'Rp_median_values.xlsx'), sheet_name=None)
Rp_lock_sheets = pd.read_excel(os.path.join(rp_dir, 'Rp_lock_median_values.xlsx'), sheet_name=None)
Rp_sheets_tohoku = pd.read_excel(os.path.join(rp_dir, '4000001_Rp_values.xlsx'), sheet_name=None)
Rp_lock_sheets_tohoku = pd.read_excel(os.path.join(rp_dir, '4000001_Rp_lock_values.xlsx'), sheet_name=None)


topo_grid = os.path.join(data_dir, 'tohoku.nc')
topo_shade = os.path.join(data_dir, 'tohoku.int')
cpt_gray = os.path.join(data_dir, 'grayscale02.cpt')

# =============================================================================
# 3. OPENQUAKE CONFIG
# =============================================================================

bins = [6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
labels = ['6.5–7.0', '7.0–7.5', '7.5–8.0', '8.0–8.5', '8.5–9.0', '9.0–9.5']
database['Mw_bin'] = pd.cut(database['Earthquake_Magnitude'], bins=bins, labels=labels)
database['logRp/Rrup'] = np.log(database['Rasp_median_km'] / database['ClstD_km'])
database['Rp/Rrup'] = database['Rasp_median_km'] / database['ClstD_km']

gmpe_parker = ParkerEtAl2020SInter(region='JP', saturation_region='JP_Pac')
gmpe_montalva = MontalvaEtAl2017SInter()
tohoku_database = database.loc[database['NGAsubEQID'] == eqid]

R_gmm_Rrup = np.logspace(np.log10(40), np.log10(700), 100)
R_gmm_Rp = np.logspace(np.log10(60), np.log10(900), 100)
Rref_gmm = functions_py.get_R_Rref_interface(tohoku_database['Earthquake_Magnitude'].unique(), R_gmm_Rp)[1]

input_df = pd.DataFrame({
    "mag": len(R_gmm_Rrup) * [tohoku_database['Earthquake_Magnitude'].unique()],
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
# 4. FIG6
# =============================================================================

pygmt.config(FONT="8p,Helvetica,black")
periods_plot = [0]
R_gmm_v2 = functions_py.get_R_Rref_interface(tohoku_database['Earthquake_Magnitude'].unique(), R_gmm_Rrup)[0]

for T in periods_plot:
    print(f"Procesando Fig6 para T={T}")
    
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

    
    print(f"  > P-value crudo: {raw_p:.2f} | P-value redondeado: {p_opt_T}")

    coeff_Rp_T = coeff_Rp.loc[(coeff_Rp['Period'] == T) & (coeff_Rp['p_value'] == p_opt_T)]
    coeff_Rrup_T = coeff_Rrup.loc[coeff_Rrup['Period'] == T]

    tohoku_df = functions_py.get_event_period_df(
        database, Fs, eqid, T, coeff_Rp, coeff_Rrup, p_value=p_opt_T,
        Rp_sheets=Rp_sheets, Rp_sheets_lock=Rp_lock_sheets, region='Global', Rref=Rref_gmm
    )
    tohoku_df = tohoku_df.dropna(subset='log_obs_rock')
    obs_T = tohoku_df['log_obs_rock']
    mag = np.array(len(R_gmm_Rrup) * [tohoku_df['Earthquake_Magnitude'].unique()[0]])

    pred_Rrup = functions_py.DR_pred_GMM_Rp(R_gmm_v2, mag, coeff_Rrup_T, Rref=Rref_gmm)
    pred_Rp = functions_py.DR_pred_GMM_Rp(R_gmm_Rp, mag, coeff_Rp_T, Rref=Rref_gmm)

    Rp = tohoku_df['Rp']
    Rp_lock_tohoku = Rp_lock_sheets_tohoku['p = ' + str(p_opt_T)]
    Rp_lock_LL2011 = Rp_lock_tohoku.loc[:, Rp_lock_tohoku.columns.str.endswith('LL2011')].median(axis=1).values
    Rp_lock_LL2016 = Rp_lock_tohoku.loc[:, Rp_lock_tohoku.columns.str.endswith('epoch4')].median(axis=1).values
    R, Rrup = tohoku_df['R'], tohoku_df['Rrup']

    ln_mean_parker, sig_parker, tau_parker, phi_parker = [np.zeros([1, n]) for _ in range(4)]
    ln_mean_montalva, sig_montalva, tau_montalva, phi_montalva = [np.zeros([1, n]) for _ in range(4)]
    gmpe_parker.compute(ctx, [imt], ln_mean_parker, sig_parker, tau_parker, phi_parker)
    gmpe_montalva.compute(ctx, [imt], ln_mean_montalva, sig_montalva, tau_montalva, phi_montalva)
    ln_mean_parker = ln_mean_parker.reshape((len(mag),))
    ln_mean_montalva = ln_mean_montalva.reshape((len(mag),))

    fig1 = pygmt.Figure()
    projection_log = "X?l/Y?"
    projection_geo = "M6.5c"
    colores = [f"#{random.randint(0, 255):02x}{random.randint(0, 255):02x}{random.randint(0, 255):02x}" for _ in range(len(modelos_cosismicos))]
    cont = 0

    with fig1.subplot(nrows=2, ncols=1, figsize=("4c", "14c"), sharex=True):
        
        with fig1.set_panel(panel=[0, 0]):
            fig1.basemap(region=[139, 145, 35, 40], projection=projection_geo, frame=['Wsne'])
            pygmt.makecpt(cmap=cpt_gray, series=[-8000, 3000, 10])
            fig1.grdimage(grid=topo_grid, projection=projection_geo, shading=topo_shade, cmap=True)
            pygmt.makecpt(cmap='hot', series=[0, 1], reverse=True)
            fig1.grdimage(grid=lock_tohoku_LL2016_grd, projection=projection_geo, cmap=True, nan_transparent=True)
            fig1.plot(x=tohoku_df['Station_Longitude_deg'], y=tohoku_df['Station_Latitude_deg'], fill='green', style='t0.2c', projection=projection_geo, pen="0.1p,black")
            fig1.plot(x=trench_tohoku['Longitude'], y=trench_tohoku['Latitude'], style="f0.5i/0.1i+r+t", fill='black', pen="1.25p,black", projection=projection_geo)
            fig1.coast(shorelines=True, projection=projection_geo)
            fig1.text(text="a) LM16", position='TL', font="8p,Helvetica,black", projection=projection_geo, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')

            for i in modelos_cosismicos['event_tag']:
                model_path = os.path.join(models_dir, str(i) + '.xyz')
                model_i = pd.read_csv(model_path, delim_whitespace=True, index_col=None, names=['Longitude', 'Latitude', 'Slip', 'Depth'])
                model_i = model_i[['Longitude', 'Latitude', 'Slip']]
                model_i['Slip'] = model_i['Slip'] / model_i['Slip'].max()
                
                resolution = np.sqrt((model_i['Longitude'].loc[2]-model_i['Longitude'].loc[1])**2+(model_i['Latitude'].loc[2]-model_i['Latitude'].loc[1])**2)
                resolution_m = 3 * np.ceil(resolution * 60)
                region = [min(model_i['Longitude']) - 1, max(model_i['Longitude']) + 1, min(model_i['Latitude']) - 1, max(model_i['Latitude']) + 1]
                model_limits = [min(model_i['Longitude']), max(model_i['Longitude']), min(model_i['Latitude']), max(model_i['Latitude'])]
                model_grid = pygmt.nearneighbor(data=model_i, spacing="1m", region=model_limits, search_radius=str(resolution_m) + 'm', sectors="6")
                
                slab_model_grid = slab_tohoku.interp(y=model_grid.lat, x=model_grid.lon)
                model_grid = model_grid.where(slab_model_grid.notnull())
                fig1.grdcontour(grid=model_grid, limit=[0, 0.8], projection=projection_geo, annotation="n", levels=0.8, pen='1p,' + colores[cont])
                cont += 1
            cont = 0

        with fig1.set_panel(panel=[1, 0]):
            fig1.basemap(region=[139, 145, 35, 40], projection=projection_geo, frame=['WSne'])
            pygmt.makecpt(cmap=cpt_gray, series=[-8000, 3000, 10])
            fig1.grdimage(grid=topo_grid, projection=projection_geo, shading=topo_shade, cmap=True)
            pygmt.makecpt(cmap='hot', series=[0, 1], reverse=True)
            fig1.grdimage(grid=lock_tohoku_LL2011_grd, projection=projection_geo, cmap=True, nan_transparent=True)
            fig1.plot(x=tohoku_df['Station_Longitude_deg'], y=tohoku_df['Station_Latitude_deg'], fill='green', style='t0.2c', projection=projection_geo, pen="0.1p,black")
            fig1.plot(x=trench_tohoku['Longitude'], y=trench_tohoku['Latitude'], style="f0.5i/0.1i+r+t", fill='black', pen="1.25p,black", projection=projection_geo)
            fig1.coast(shorelines=True, projection=projection_geo)
            fig1.text(text="b) LM11", position='TL', font="8p,Helvetica,black", projection=projection_geo, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            
            pygmt.config(FONT="10p,Helvetica,black")
            fig1.colorbar(frame="xa.25f.25+lLocking degree", position="g139.5/34+w5.5c/0.4c+h", projection=projection_geo, region=[139, 145, 35, 40])
            pygmt.config(FONT="8p,Helvetica,black")
            
            for i in modelos_cosismicos['event_tag']:
                model_path = os.path.join(models_dir, str(i) + '.xyz')
                model_i = pd.read_csv(model_path, delim_whitespace=True, index_col=None, names=['Longitude', 'Latitude', 'Slip', 'Depth'])
                model_i = model_i[['Longitude', 'Latitude', 'Slip']]
                model_i['Slip'] = model_i['Slip'] / model_i['Slip'].max()
                
                resolution = np.sqrt((model_i['Longitude'].loc[2]-model_i['Longitude'].loc[1])**2+(model_i['Latitude'].loc[2]-model_i['Latitude'].loc[1])**2)
                resolution_m = 3 * np.ceil(resolution * 60)
                region = [min(model_i['Longitude']) - 1, max(model_i['Longitude']) + 1, min(model_i['Latitude']) - 1, max(model_i['Latitude']) + 1]
                model_limits = [min(model_i['Longitude']), max(model_i['Longitude']), min(model_i['Latitude']), max(model_i['Latitude'])]
                model_grid = pygmt.nearneighbor(data=model_i, spacing="1m", region=model_limits, search_radius=str(resolution_m) + 'm', sectors="6")
                
                slab_model_grid = slab_tohoku.interp(y=model_grid.lat, x=model_grid.lon)
                model_grid = model_grid.where(slab_model_grid.notnull())
                fig1.grdcontour(grid=model_grid, limit=[0, 0.8], projection=projection_geo, annotation="n", levels=0.8, pen='1p,' + colores[cont])
                cont += 1

    fig1.shift_origin(xshift="w+4.8c", yshift="-0.0c")
    with fig1.subplot(nrows=2, ncols=2, figsize=("14c", "14c"), sharey=True):
        y_min = min(min(obs_T), min(pred_Rp), min(ln_mean_parker), min(ln_mean_montalva)) - 0.5
        y_max = max(max(obs_T), max(pred_Rp), max(ln_mean_parker), max(ln_mean_montalva)) + 0.5
        
        with fig1.set_panel(panel=[0, 0]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['Wsne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rp_lock_LL2016, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)

            fig1.plot(x=R_gmm_Rp, y=pred_Rp, pen='1p,red', label='This model')
            fig1.text(text="c) R@-p,lock@- LM16", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.8c", box="+gwhite+p1p")

        with fig1.set_panel(panel=[1, 0]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['WSne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rp_lock_LL2011, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)

            fig1.plot(x=R_gmm_Rp, y=pred_Rp, pen='1p,red', label='This model')
            fig1.text(text="d) R@-p,lock@- LM11", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.8c", box="+gwhite+p1p")
            fig1.text(text="Distance [km]", position='BR', font="12p,Helvetica,black", projection=projection_log, offset='2/-1', no_clip=True)
            fig1.text(text="PGA at V@-s30@- = 760 m/s [log units]", position='ML', font="12p,Helvetica,black", projection=projection_log, offset='-1.2/0.5', no_clip=True, angle=90)

        with fig1.set_panel(panel=[0, 1]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['wsne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rp, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)

            fig1.plot(x=R_gmm_Rp, y=pred_Rp, pen='1p,red', label='This model')
            fig1.text(text="e) R@-p@-", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.8c", box="+gwhite+p1p")

        with fig1.set_panel(panel=[1, 1]):
            fig1.basemap(region=[40, 1000, y_min, y_max], projection=projection_log, frame=['wSne', 'xa2fg3', 'ya2fg3'])
            fig1.plot(x=Rrup, y=obs_T, style='c0.1c', fill='gray', pen='0.1p,black', intensity=0.3)

            fig1.plot(x=R_gmm_Rrup, y=pred_Rrup, pen='1p,red', label='This model')
            fig1.plot(x=R_gmm_Rrup, y=ln_mean_parker, pen='1p,blue', label='Parker et al. (2021)')
            fig1.text(text="f) R@-rup@-", position='TL', font="8p,Helvetica,black", projection=projection_log, fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w3.8c", box="+gwhite+p1p")

    fig1.savefig(os.path.join(output_dir, 'Fig6.pdf'))

# =============================================================================
# 5. FIG 9 AND FIG 10
# =============================================================================

n_stations = 12
periods_plot = [0]

for T in periods_plot:
    print(f"Procesando Fig9 y Fig10 para T={T}")
    
    if T == -1: column_obs = 'PGV_cm_sec'
    elif T == 0: column_obs = 'PGA_g'
    else: column_obs = f'T = {float(T)}'

    raw_p = p_opt_Rp['p_opt_pred'].loc[p_opt_Rp['Period'] == T].values[0]
    p_opt_T = round(raw_p * 2) / 2
    
    coeff_Rp_T = coeff_Rp.loc[(coeff_Rp['Period'] == T) & (coeff_Rp['p_value'] == p_opt_T)]
    
    Rp_sheets_tohoku_T = Rp_sheets_tohoku['p = ' + str(p_opt_T)]
    
    Rp_sheets_tohoku_T = Rp_sheets_tohoku_T.loc[Rp_sheets_tohoku_T.loc[:, Rp_sheets_tohoku_T.columns.str.contains('2011')].median(axis=1) < 300]
    Rp_sheets_tohoku_T = Rp_sheets_tohoku_T.sample(n=n_stations)
    Rp_T_models = Rp_sheets_tohoku_T.loc[:, Rp_sheets_tohoku_T.columns.str.contains('2011')]
    
    fig1 = pygmt.Figure()
    fig2 = pygmt.Figure()
    projection_log = "X?/Y?"
    row, col = 0, 0
    row1, col1 = 0, 0
    
    for i in range(n_stations):
        sta_name = Rp_sheets_tohoku_T['Station_Name'].iloc[i]
        obs_sta = np.log(tohoku_database[column_obs].loc[tohoku_database['Station_Name'] == sta_name].values[0])
        Rp_T_sta_i = Rp_T_models.iloc[i].values
        Rp_sta_median = np.median(Rp_T_sta_i)
        
        pred_Rp_T = functions_py.DR_pred_GMM_Rp(Rp_T_sta_i, 9.1, coeff_Rp_T, Rref=Rref_gmm)
        pred_Rp_median = functions_py.DR_pred_GMM_Rp(np.median(Rp_T_sta_i), 9.1, coeff_Rp_T, Rref=Rref_gmm)

        with fig1.subplot(nrows=3, ncols=4, figsize=("14c", "10c"), sharex=True):
            with fig1.set_panel(panel=[row, col]):
                frame_cfg = 'Wsne' if col == 0 and row == 0 else ('WSne' if col == 0 and row == 2 else ('wSne' if row == 2 else 'wsne'))
                fig1.basemap(region=[-8, 1, 0, 15], projection=projection_log, frame=[frame_cfg])
                
                if col == 0 and row == 1:
                    fig1.text(text='Counts', position='MC', font="10p,Helvetica,black", projection=projection_log, offset='-2.5/0', no_clip=True, angle=90)
                
                fig1.histogram(data=pred_Rp_T, series=[-10, 2, 0.5], fill="grey", pen="1p", histtype=0)
                
                label_obs = 'Observed+N2' if (col == 0 and row == 0) else None
                label_med = 'Median' if (col == 0 and row == 0) else None
                
                fig1.plot(x=[obs_sta, obs_sta], y=[-50, 50], projection=projection_log, pen='0.8p,red,-', label=label_obs)
                fig1.plot(x=[pred_Rp_median, pred_Rp_median], y=[-50, 50], projection=projection_log, pen='0.8p,blue,-', label=label_med)
                
                if col == 0 and row == 0:
                    fig1.legend(position="JBR+jBR+o-9.5c/-0.6c+w9c")
                if col == 1 and row == 2:
                    fig1.text(text='PGA at V@-s30@- = 760 m/s [log units]', position='BR', font="10p,Helvetica,black", projection=projection_log, offset='2.5/-1', no_clip=True)
            
            col += 1
            if col > 3: col, row = 0, row + 1

        with fig2.subplot(nrows=3, ncols=4, figsize=("14c", "10c"), sharex=True):
            with fig2.set_panel(panel=[row1, col1]):
                frame_cfg = 'Wsne' if col1 == 0 and row1 == 0 else ('WSne' if col1 == 0 and row1 == 2 else ('wSne' if row1 == 2 else 'wsne'))
                fig2.basemap(region=[40, 300, 0, 15], projection=projection_log, frame=[frame_cfg])
                
                if col1 == 0 and row1 == 1:
                    fig2.text(text='Counts', position='MC', font="10p,Helvetica,black", projection=projection_log, offset='-2.5/0', no_clip=True, angle=90)
                
                fig2.histogram(data=Rp_T_sta_i, series=[50, 300, 25], fill="grey", pen="1p", histtype=0)
                
                label_med = 'Median+N1' if (col1 == 0 and row1 == 0) else None
                fig2.plot(x=[Rp_sta_median, Rp_sta_median], y=[-50, 50], projection=projection_log, pen='0.8p,blue,-', label=label_med)
                
                if col1 == 0 and row1 == 0:
                    fig2.legend(position="JBR+jBR+o-11.8c/-0.6c+w9c")
                if col1 == 1 and row1 == 2:
                    fig2.text(text='R@-p@- [km]', position='BR', font="10p,Helvetica,black", projection=projection_log, offset='5.0/-1', no_clip=True)
            
            col1 += 1
            if col1 > 3: col1, row1 = 0, row1 + 1

    fig1.savefig(os.path.join(output_dir, 'Fig10.pdf'), dpi=300)
    fig2.savefig(os.path.join(output_dir, 'Fig9.pdf'), dpi=300)

