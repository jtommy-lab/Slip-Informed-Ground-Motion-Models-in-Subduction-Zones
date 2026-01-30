import numpy as np
import pandas as pd
import pygmt 
import math
import sys
sys.path.append('/home/jtommy/Escritorio/Respaldo/functions/')
import metricas_distancia
from skimage.feature import peak_local_max

def get_slab_filename(Model_Longitudes,Model_Latitudes):
    slab_info = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/slab2/grid/slab2_info.csv')
    Model_Longitudes[Model_Longitudes < 0] = Model_Longitudes[Model_Longitudes < 0] + 360
    slab_info = slab_info.loc[(slab_info['Minimum Longitude'] <= min(Model_Longitudes)) & (slab_info['Maximum Longitude'] >= max(Model_Longitudes)) &
                              (slab_info['Minimum Latitude'] <= min(Model_Latitudes)) & (slab_info['Maximum Latitude'] >= max(Model_Latitudes))]
    return(slab_info['Filename'].to_list()[0])

def get_slab_depth(Model_Longitudes,Model_Latitudes):
    grdfile = '/home/jtommy/Escritorio/Respaldo/base_de_datos/slab2/grid/'+get_slab_filename(Model_Longitudes,Model_Latitudes)
    Model_Longitudes[Model_Longitudes < 0] = Model_Longitudes[Model_Longitudes < 0] + 360
    Model_df = pd.DataFrame({'Longitudes': Model_Longitudes, 'Latitudes': Model_Latitudes})
    Model_depth = pygmt.grdtrack(grid = grdfile,points = Model_df, newcolname = 'Depth (km)',radius = True)
    Model_depth['Depth (km)'] = Model_depth['Depth (km)'] * -1
    return(Model_depth)

def get_R_Rref_interface(Magnitude,Rrup):
    h = 10**(-0.82+0.252*Magnitude)
    R = np.sqrt(Rrup**2+h**2)
    Rref = np.sqrt(1+h**2)
    return(R,Rref)

def PSHA20_Fp_interface(c1,b4,a0,Magnitude,Rmetric,h_flag = False):

    if h_flag == True:
        R,Rref = get_R_Rref_interface(Magnitude,Rmetric)
        R_Rref = R/Rref
        Fp = c1*np.log(R)+b4*Magnitude*np.log(R_Rref)+a0*R
    else:
        Fp = c1*np.log(Rmetric)+b4*Magnitude*np.log(Rmetric)+a0*Rmetric

    return(Fp)

def PSHA20_Fm_interface(c4,c5,c6,Magnitude,Mc = 7.9):
    
    if Magnitude<= Mc:
        Fm = c4*(Magnitude-Mc) + c5*(Magnitude-Mc)**2
    else:
        Fm = c6*(Magnitude-Mc)

    return(Fm)

def get_event_period_df(database,Fs,eqid,period,model_coeff_Rp,model_coeff_Rrup,p_value = 0,Rp_sheets = 0,Rp_sheets_lock = 0,region = 'Global',Rref = 1):
    
    model_coeff_Rp = model_coeff_Rp.loc[(model_coeff_Rp['Period'] == period) & (model_coeff_Rp['p_value'] == p_value) ]
    model_coeff_Rrup = model_coeff_Rrup.loc[model_coeff_Rrup['Period'] == period]
    data_eq = database.loc[database['NGAsubEQID'] == str(eqid)]
    Fs_eq = Fs.loc[Fs['NGAsubEQID'] == str(eqid)]
    data_eq_merge = pd.merge(data_eq,Fs_eq,how="inner", on=['NGAsubEQID','Station_Name','Earthquake_Magnitude','Vs30_Selected_for_Analysis_m_s'])
    if p_value != 0:
        Rp_df = Rp_sheets['p = '+str(p_value)]
        data_eq_merge = pd.merge(data_eq_merge,Rp_df,how="inner", on=['NGAsubEQID','Station_Name','Earthquake_Magnitude','Vs30_Selected_for_Analysis_m_s',
                                                                        'Station_Longitude_deg','Station_Latitude_deg'])
        if Rp_sheets_lock != 0:
            Rp_lock_df = Rp_sheets_lock['p = '+str(p_value)]
            data_eq_merge = pd.merge(data_eq_merge,Rp_lock_df,how="inner", on=['NGAsubEQID','Station_Name','Earthquake_Magnitude','Vs30_Selected_for_Analysis_m_s',
                                                                        'Station_Longitude_deg','Station_Latitude_deg'])

    if period == -1:
        column_obs = 'PGV_cm_sec'
        column_Fs = 'Fs - PGV_cm_sec'
    elif period == 0:
        column_obs = 'PGA_g'
        column_Fs = 'Fs - PGA_g'
    else:
        column_obs = 'T = '+str(float(period))
        column_Fs = 'Fs - T = '+str(float(period))
        
    R =  get_R_Rref_interface(data_eq_merge['Earthquake_Magnitude'],data_eq_merge['ClstD_km'])[0]
    data_eq_T = pd.DataFrame({'NGAsubEQID': data_eq_merge['NGAsubEQID'],'Station_Name': data_eq_merge['Station_Name'],'Station_Longitude_deg':data_eq_merge['Station_Longitude_deg'],
                            'Station_Latitude_deg':data_eq_merge['Station_Latitude_deg'],'Earthquake_Magnitude': data_eq_merge['Earthquake_Magnitude'],'Rp': data_eq_merge['Rp_median_km'],
                            'R':R,'Rrup': data_eq_merge['ClstD_km'],'log_obs': np.log(data_eq_merge[column_obs]),'log_obs_rock': np.log(data_eq_merge[column_obs])-data_eq_merge[column_Fs]})
    if Rp_sheets_lock != 0:
        data_eq_T['Rp_lock'] = data_eq_merge['Rp_lock_median_km']
        data_eq_T['pred_Rp_lock'] = model_coeff_Rp[str(region)+'_c0'].values[0] + model_coeff_Rp['c1'].values[0]*np.log(data_eq_T['Rp_lock']) + model_coeff_Rp[str(region)+'_a0'].values[0]*data_eq_T['Rp_lock'] + model_coeff_Rp['c4'].values[0]*data_eq_T['Earthquake_Magnitude']
    data_eq_T['pred_Rp'] = model_coeff_Rp[str(region)+'_c0'].values[0] + model_coeff_Rp['c1'].values[0]*np.log(data_eq_T['Rp']) + model_coeff_Rp[str(region)+'_a0'].values[0]*data_eq_T['Rp'] + model_coeff_Rp['c4'].values[0]*data_eq_T['Earthquake_Magnitude']
    if "b4" in model_coeff_Rp.columns:
        data_eq_T['pred_Rp'] = data_eq_T['pred_Rp'] + model_coeff_Rp["b4"].values[0] * data_eq_T['Earthquake_Magnitude'] * np.log(data_eq_T['Rp'] / Rref)
        data_eq_T['pred_Rp_lock'] = data_eq_T['pred_Rp'] + model_coeff_Rp["b4"].values[0] * data_eq_T['Earthquake_Magnitude'] * np.log(data_eq_T['Rp'] / Rref)
    data_eq_T['pred_Rrup'] = model_coeff_Rrup[str(region)+'_c0'].values[0] + model_coeff_Rrup['c1'].values[0]*np.log(data_eq_T['R']) + model_coeff_Rrup[str(region)+'_a0'].values[0]*data_eq_T['R'] + model_coeff_Rrup['c4'].values[0]*data_eq_T['Earthquake_Magnitude']  
    if "b4" in model_coeff_Rrup.columns:
        data_eq_T['pred_Rrup'] = data_eq_T['pred_Rrup'] + model_coeff_Rrup["b4"].values[0] * data_eq_T['Earthquake_Magnitude'] * np.log(data_eq_T['R'] / Rref)
    
    return(data_eq_T)

def DR_pred_GMM_Rp(Rp, Magnitude, coefficients, Rref = 1, region = 'Global'):    
    
    prediction = coefficients[str(region)+'_c0'].values[0] + coefficients['c1'].values[0]*np.log(Rp) + coefficients[str(region)+'_a0'].values[0]*Rp + coefficients['c4'].values[0]*Magnitude

    if "b4" in coefficients.columns:
        prediction = prediction + coefficients["b4"].values[0] * Magnitude * np.log(Rp / Rref)
   
    return(prediction)

def rotar_ff(lon_polo,lat_polo,lon_grilla,lat_grilla,rot_polo): 
    lat_buenas=[]
    lon_buenas=[]   
    for j in range(len(lat_grilla)):
        c=np.radians(lat_polo)
        d=np.radians(lon_polo)
        polo_60= np.array([np.cos(c)*np.cos(d), np.cos(c)*np.sin(d), np.sin(c)])
        a=np.radians(lat_grilla[j])
        b=np.radians(lon_grilla[j])
        ciudad=np.array([np.cos(a)*np.cos(b), np.cos(a)*np.sin(b),np.sin(a)])
            #generando matriz de rotacion
        co60=np.cos(rot_polo)
        se60=np.sin(rot_polo)
        r11=polo_60[0]*polo_60[0]*(1-co60)+co60
        r12=polo_60[0]*polo_60[1]*(1-co60)-polo_60[2]*se60
        r13=polo_60[0]*polo_60[2]*(1-co60)+polo_60[1]*se60
        r21=polo_60[1]*polo_60[0]*(1-co60)+polo_60[2]*se60
        r22=polo_60[1]*polo_60[1]*(1-co60)+co60
        r23=polo_60[1]*polo_60[2]*(1-co60)-polo_60[0]*se60
        r31=polo_60[2]*polo_60[0]*(1-co60)-polo_60[1]*se60
        r32=polo_60[2]*polo_60[1]*(1-co60)+polo_60[0]*se60
        r33=polo_60[2]*polo_60[2]*(1-co60)+co60    
            #matriz de rotacion
        r60= np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])       
        x60 = np.dot(r60,ciudad)
        lo_60= np.degrees(np.arctan2(x60[1],x60[0]))
        la_60=np.degrees(np.arctan(x60[2]/np.sqrt(x60[0]**2+x60[1]**2)))
        lon_buenas.append(lo_60); lat_buenas.append(la_60)        
    return np.array(lon_buenas), np.array(lat_buenas)

def calculate_s2_s1(region, coeff_T):
    """
    region: str
    coeff_T: pandas DataFrame con columnas como 'Global_s2', 'Taiwan_s1', etc.
             o pandas Series con esos nombres.
    """
    if region in ("global", "CAM"):
        s2 = coeff_T["Global_s2"]
        s1 = s2
    elif region in ("Taiwan", "Japan"):
        s2 = float(coeff_T[f"{region}_s2"])
        s1 = float(coeff_T[f"{region}_s1"])
    else:
        s2 = float(coeff_T[f"{region}_s2"])
        s1 = s2
    return {"s2": s2, "s1": s1}

import math

def linear_site_term_PS21(s1, s2, Vs30, V1, V2, Vref):
    """
    Calcula el término de sitio lineal del modelo PS21.

    Parámetros:
    s1, s2 : floats
        Coeficientes de pendiente para los tramos bajos y altos.
    Vs30 : float
        Velocidad de onda de corte en m/s.
    V1, V2 : float
        Límites de transición para Vs30.
    Vref : float
        Velocidad de referencia.

    Retorna:
    float : valor de Flin.
    """
    if Vs30 <= V1:
        Flin = s1 * math.log(Vs30 / V1) + s2 * math.log(V1 / Vref)
    elif Vs30 <= V2:
        Flin = s2 * math.log(Vs30 / Vref)
    else:
        Flin = s2 * math.log(V2 / Vref)
    return Flin

import math

def nonlinear_site_term_PS21(T, coeff_T, Vs30, Vref_Fnl, Vb, PGAr, f3):
    """
    Calcula el término no lineal de sitio (PS21).

    Parámetros:
    T : float
        Periodo.
    coeff_T : dict o pandas.Series
        Coeficientes, debe contener 'f4' y 'f5'.
    Vs30, Vref_Fnl, Vb, PGAr, f3 : float
        Parámetros del cálculo.

    Retorna:
    float : valor de Fnl.
    """
    if T >= 3:
        return 0.0
    f4 = coeff_T["f4"].values[0]
    f5 = coeff_T["f5"].values[0]
    vsmin = min(Vs30, Vref_Fnl)  # equivalente a pmin escalar
    f2 = f4 * (math.exp(f5 * (vsmin - Vb)) - math.exp(f5 * (Vref_Fnl - Vb)))
    Fnl = f2 * math.log((PGAr + f3) / f3)
    return Fnl





def PS21_Fs_SA(Rp,Magnitude,Vs30,Period,metric = 'Rp'):
    
    coefficients = pd.read_csv("/home/jtommy/Escritorio/Respaldo/base_de_datos/GMPE/Parkeretal2021/PSHAB20_Table_E1_Interface_Coefficients_OneRowHeader.csv")
    coeff_T = coefficients.loc[coefficients['Period'] == Period]
    coeff_Rp = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer_regional.csv')
    coeff_Rrup = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rrup/nlmer/coeficientes_Rrup_nlmer_regional.csv')
    coeff_Rrup_T = coeff_Rrup.loc[coeff_Rrup['Period'] == Period]
    p_value = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/p_selected_Rp_regional_smooth.csv')
    p_value_T = round(p_value.loc[p_value['Period'] == Period]['median_p_smooth']).values[0]
    coeff_Rp_T = coeff_Rp.loc[(coeff_Rp['Period'] == Period) & (coeff_Rp['p_value'] == p_value_T)]
    Vb = 200
    Vref_Fnl = 760
    f3 = 0.05
    R,Rref = get_R_Rref_interface(Magnitude,Rp)
    
    ### Linear site term
    
    V1 = coeff_T['V1_m_s'].values[0]
    V2 = coeff_T['V2_m_s'].values[0]
    Vref = coeff_T['Vref_m_s'].values[0]
    s1s2 = calculate_s2_s1('SA', coeff_T)
    s1 = s1s2['s1']
    s2 = s1s2['s2']
    Flin = linear_site_term_PS21(s1, s2, Vb, V1, V2, Vref_Fnl) 
    if metric == 'Rrup':
        PGAr = np.exp(DR_pred_GMM_Rp(R, Magnitude, coeff_Rrup_T, Rref = Rref, region = 'SA'))
    else:
        PGAr = np.exp(DR_pred_GMM_Rp(Rp, Magnitude, coeff_Rp_T, Rref = 1, region = 'SA'))
    Fnl = nonlinear_site_term_PS21(Period, coeff_T, Vs30, Vref_Fnl, Vb, PGAr, f3)
    Fs = Flin + Fnl
    
    return(Fs)

def local_peaks_km(da, min_sep_km=80, threshold_rel=0.2, max_peaks=None):
    '''
    Docstring for local_peaks_km
    
    :param da: model_grid
    :param min_sep_km: distance between maximums
    :param threshold_rel: minimum ratio between local maxima and global maxima
    :param max_peaks: Max number of peaks
    '''
    Z = da.values
    Zf = np.where(np.isfinite(Z), Z, -np.inf)
    # spacing = 1m -> 1 arc-minute por celda
    lat0 = float(np.nanmean(da["lat"].values))
    km_per_cell_lat = 111.32 / 60.0                  # ~1.855 km
    km_per_cell_lon = (111.32 * np.cos(np.deg2rad(lat0))) / 60.0
    # Para que SIEMPRE asegures al menos min_sep_km en cualquier dirección,
    # usa la escala más "fina" (la menor km/celda) -> distancia en celdas más grande.
    km_per_cell = min(km_per_cell_lat, km_per_cell_lon)
    min_dist_cells = max(1, int(np.ceil(min_sep_km / km_per_cell)))
    coords_rc = peak_local_max(
        Zf,
        min_distance=min_dist_cells,      # en celdas
        threshold_rel=threshold_rel,      # relativo al máximo global
        num_peaks=max_peaks,              # None = sin límite
        exclude_border=False
    )
    lats = da.lat.values[coords_rc[:, 0]]
    lons = da.lon.values[coords_rc[:, 1]]
    vals = Z[coords_rc[:, 0], coords_rc[:, 1]]
    out = pd.DataFrame({"lat": lats, "lon": lons, "z": vals})
    out = out.sort_values("z", ascending=False).reset_index(drop=True)
    return out, min_dist_cells, km_per_cell_lat, km_per_cell_lon
    
    


    
    