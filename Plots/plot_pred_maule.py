import numpy as np
import pandas as pd
import pygmt
import sys
sys.path.append('/home/jtommy/Escritorio/Respaldo/functions')  # Añade la ruta al sistema
sys.path.append('/home/jtommy/Escritorio/Respaldo/base_de_datos/GMPE/Montalvaetal2017')  # Añade la ruta al sistema

import functions_py
import random
from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SInter
from openquake.hazardlib.gsim.montalva_2017 import MontalvaEtAl2017SInter
import openquake.hazardlib.imt as IMT
import funcion_GMPE



eqid = '6000149'

database = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/Drapela_database.csv',index_col = None)
Fs = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/SS14_Fs_DR_database.csv',index_col = None)

coeff_Rrup = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rrup/nlmer/coeficientes_Rrup_nlmer.csv',index_col = False)
coeff_Rrup_regional = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rrup/nlmer/coeficientes_Rrup_nlmer_regional.csv',index_col = False)
coeff_Rp = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer.csv',index_col = False)
coeff_Rp_regional = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer_regional.csv',index_col = False)

trench_chile = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Trench/trench-chile',delim_whitespace=True,index_col = None,names=['Longitude','Latitude','Depth'])

modelos_cosismicos = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_cosismicos/modelos_cosismicos.csv',index_col = False)
modelos_cosismicos = modelos_cosismicos.loc[modelos_cosismicos['NGAsubEQID'] == eqid]

slab_soam = pygmt.load_dataarray('/home/jtommy/Escritorio/Respaldo/base_de_datos/slab2/grid/sam_slab2_dep_02.23.18.grd')
lock_chile = pygmt.load_dataarray('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/Chile_locking.grd')

Rp_sheets = pd.read_excel('/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/Rp_median_values.xlsx',sheet_name=None)
Rp_lock_sheets = pd.read_excel('/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/Rp_lock_median_values.xlsx',sheet_name=None)

p_opt_Rp = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/p_selected_Rp_smooth.csv',index_col = None)
p_opt_Rp_regional = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/p_selected_Rp_regional_smooth.csv',index_col = None)

gmpe_parker = ParkerEtAl2020SInter(region='SA',saturation_region='SA_S')
gmpe_montalva = MontalvaEtAl2017SInter()
maule_database = database.loc[database['NGAsubEQID'] == eqid]
maule_database = maule_database.sort_values(by='ClstD_km')

R_gmm_Rrup = np.logspace(np.log10(30),np.log10(300),100)
R_gmm_Rp = np.logspace(np.log10(40),np.log10(300),100)
Rref_gmm = functions_py.get_R_Rref_interface(maule_database['Earthquake_Magnitude'].unique(),R_gmm_Rrup)[1]
input_df = pd.DataFrame({"mag": len(R_gmm_Rrup) * [maule_database['Earthquake_Magnitude'].unique()[0]], "rrup" : R_gmm_Rrup, 
                            "vs30": len(R_gmm_Rrup)*[760], "region": len(R_gmm_Rrup)*['SA'], "saturation_region": len(R_gmm_Rrup)*['SA_S'],"backarc":len(R_gmm_Rrup)*[False]})
n = input_df.shape[0]
ctx = np.recarray(n,
                dtype=np.dtype([("mag", float),
                                ("rrup", float),
                                ("vs30", float),
                                ("region", str),
                                ("saturation_region", str),
                                ("backarc", bool)
                                ]),)
ctx["rrup"] = input_df["rrup"]
ctx["mag"] = input_df["mag"]
ctx["vs30"] = input_df["vs30"]
ctx["region"] = input_df["region"]
ctx["saturation_region"] = input_df["saturation_region"]

predictions = pd.DataFrame(columns=["mag", "rrup", "vs30", "period","ln_mean", "sig", "tau", "phi"])
#functions_py.get_R_Rref_interface(maule_database['Earthquake_Magnitude'],maule_database['ClstD_km'])[0].values
pred_HQ_maule = []
for Ri in R_gmm_Rrup:
    T,Sa,sig,phi,phiSS,phiS2S,tau, f_site = funcion_GMPE.GMPE_CHILE_MBR(R = Ri,M = 8.8,I = 1,Zh = 30.4055,Vs30 = 760,F_FABA = 0,epsilon = 0,HQ = 1)
    pred_HQ_maule.append(np.log(Sa[0]))

#gmt grdgradient '/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/maule.nc' -G'/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/maule.int' -A225 -Nt

periods_plot = [0]
##[-1,0,0.5,0.75,3,10]
R_gmm_v2 = functions_py.get_R_Rref_interface(maule_database['Earthquake_Magnitude'].unique(),R_gmm_Rrup)[0]
for T in periods_plot:
    if T == -1:
        column_obs = 'PGV_cm_sec'
        column_Fs = 'Fs - PGV_cm_sec'
        imt = IMT.PGV() 
    elif T == 0:
        column_obs = 'PGA_g'
        column_Fs = 'Fs - PGA_g'
        imt = IMT.PGA()
    else:
        column_obs = 'T = '+str(float(T))
        column_Fs = 'Fs - T = '+str(float(T))
        imt = IMT.SA(T) 
    p_opt_T = round(p_opt_Rp['median_p_smooth'].loc[p_opt_Rp['Period'] == T].values[0])
    p_opt_regional_T = round(p_opt_Rp_regional['median_p_smooth'].loc[p_opt_Rp_regional['Period'] == T].values[0])
    coeff_Rp_regional_T = coeff_Rp_regional.loc[(coeff_Rp_regional['Period'] == T) & (coeff_Rp_regional['p_value'] == p_opt_regional_T) ]
    coeff_Rp_T = coeff_Rp.loc[(coeff_Rp['Period'] == T) & (coeff_Rp['p_value'] == p_opt_T) ]
    coeff_Rrup_T = coeff_Rrup.loc[coeff_Rrup['Period'] == T]
    coeff_Rrup_regional_T = coeff_Rrup_regional.loc[coeff_Rrup_regional['Period'] == T]
    maule_df = functions_py.get_event_period_df(database,Fs,eqid,T,coeff_Rp,coeff_Rrup,p_value=p_opt_T,Rp_sheets=Rp_sheets,Rp_sheets_lock=Rp_lock_sheets,region = 'Global',Rref = Rref_gmm)   
    maule_df = maule_df.dropna(subset='log_obs_rock')
    maule_nearest = maule_df.nsmallest(4,'Rp')
    obs_T = maule_df['log_obs_rock']
    mag = np.array(len(R_gmm_Rp)*[maule_df['Earthquake_Magnitude'].unique()[0]])
    R = maule_df['R']
    Rrup = maule_df['Rrup']
    Rp = maule_df['Rp']
    Rp_lock = maule_df['Rp_lock']
    pred_Rrup = functions_py.DR_pred_GMM_Rp(R_gmm_v2,mag,coeff_Rrup_T,Rref = Rref_gmm)
    pred_Rrup_regional = functions_py.DR_pred_GMM_Rp(R_gmm_v2,mag,coeff_Rrup_regional_T,Rref = Rref_gmm,region = 'SA')
    pred_Rp = functions_py.DR_pred_GMM_Rp(R_gmm_Rp,mag,coeff_Rp_T,Rref = Rref_gmm)
    pred_Rp_regional = functions_py.DR_pred_GMM_Rp(R_gmm_Rp,mag,coeff_Rp_regional_T,Rref = Rref_gmm, region = 'SA')
    #### GMMs externas ###
    ln_mean_parker = np.zeros([1, n])
    sig_parker = np.zeros([1, n])
    tau_parker = np.zeros([1, n])
    phi_parker = np.zeros([1, n])
    ln_mean_montalva = np.zeros([1, n])
    sig_montalva = np.zeros([1, n])
    tau_montalva = np.zeros([1, n])
    phi_montalva = np.zeros([1, n])
    gmpe_parker.compute(ctx, [imt], ln_mean_parker, sig_parker, tau_parker, phi_parker)
    gmpe_montalva.compute(ctx,[imt], ln_mean_montalva, sig_montalva, tau_montalva, phi_montalva)     
    ln_mean_parker = ln_mean_parker.reshape((len(R_gmm_Rrup),))
    ln_mean_montalva = ln_mean_montalva.reshape((len(R_gmm_Rrup),))
    fig1 = pygmt.Figure()
    pygmt.config(FONT="8p,Helvetica,black",FONT_LABEL ="10p,Helvetica,black")
    projection_log = "X?l/Y?"
    projection_geo = "M6.5c"
    colores = [f"#{random.randint(0, 255):02x}{random.randint(0, 255):02x}{random.randint(0, 255):02x}" for _ in range(len(modelos_cosismicos))]
    cont = 0
    with fig1.subplot(nrows = 1, ncols = 1, figsize=("4c","12c")):
        with fig1.set_panel(panel = [0,0]):
            fig1.basemap(region=[-75, -71,-38,-34], projection=projection_geo, frame=['WSne','xa2f1','ya2f1'])        
            pygmt.makecpt(cmap='/home/jtommy/Escritorio/graficos_GMT/paletas/grayscale02.cpt',series = [-30000,6000,10])
            fig1.grdimage(grid ='/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/maule.nc',projection = projection_geo,
                        shading ='/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/maule.int', cmap=True )
            pygmt.makecpt(cmap='hot',series = [0,1],reverse=True)
            fig1.grdimage(grid = lock_chile,projection = projection_geo, cmap=True,nan_transparent = True ) 
            pygmt.config(FONT="10p,Helvetica,black")
            fig1.colorbar(frame="xaf+lLocking degree",position="g-74.7/-38.5+w5.5c/0.4c+h",projection=projection_geo,region=[-75, -71,-38,-34])
            pygmt.config(FONT="8p,Helvetica,black")   
            fig1.text(text="a) MM10", position='TL', font="8p,Helvetica,black", projection=projection_geo,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1')            
            fig1.coast(shorelines=True,projection = projection_geo)
            fig1.plot(x=maule_df['Station_Longitude_deg'], y=maule_df['Station_Latitude_deg'], fill='green',style='t0.2c',projection = projection_geo,pen="0.1p,black")
            fig1.plot(x=maule_nearest['Station_Longitude_deg'], y=maule_nearest['Station_Latitude_deg'], fill='blue',style='t0.2c',projection = projection_geo,pen="0.1p,black")
            for i in modelos_cosismicos['event_tag']:
                print(i)
                model_i = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_cosismicos/'+str(i)+'.xyz',delim_whitespace=True,index_col = None,names=['Longitude','Latitude','Slip','Depth'])
                model_i = model_i[['Longitude','Latitude','Slip']]
                model_i['Slip'] = model_i['Slip']/model_i['Slip'].max()
                model_i['Longitude'].loc[model_i['Longitude'] < 0] = model_i['Longitude'].loc[model_i['Longitude'] < 0] + 360 
                resolution = np.sqrt((model_i['Longitude'].loc[2]-model_i['Longitude'].loc[1])**2+(model_i['Latitude'].loc[2]-model_i['Latitude'].loc[1])**2)
                resolution_m = 3*np.ceil(resolution*60)
                slip_th30 = 0.8
                slip_th50 = 0.5
                max_slip = 0.05
                region = [min(model_i['Longitude'])-1, max(model_i['Longitude'])+1, min(model_i['Latitude'])-1, max(model_i['Latitude'])+1] 
                model_limits = [min(model_i['Longitude']), max(model_i['Longitude']), min(model_i['Latitude']), max(model_i['Latitude'])] 
                model_grid = pygmt.nearneighbor(data=model_i,spacing="1m",region=model_limits,search_radius=str(resolution_m)+'m',sectors="6")
                slab_model_grid = slab_soam.interp(y=model_grid.lat,x=model_grid.lon)
                model_grid = model_grid.where(slab_model_grid.notnull())
                #fig1.grdimage(grid =model_grid,projection = projection_geo, cmap=True ) 
                fig1.grdcontour(grid = model_grid,limit = [0,slip_th30],projection = projection_geo,annotation = "n",levels=slip_th30,pen='1p,'+colores[cont])
                #fig1.grdcontour(grid = model_grid,limit = [0,slip_th50],projection = projection_geo,annotation = "n",levels=slip_th50,pen='1p,'+colores[cont])
                #fig1.grdcontour(grid = grid,limit = [0,max_slip],projection = projection_geo,annotation = "n",levels=max_slip,pen='1.8p,darkgrey') 
                fig1.plot(x=trench_chile['Longitude'], y=trench_chile['Latitude'], style="f0.5i/0.1i+r+t",fill='black', pen="1.25p,black",projection = projection_geo) 
                cont = cont+1
    fig1.shift_origin(xshift="w+4.8c",yshift="-2.5c")
    with fig1.subplot(nrows = 3, ncols = 1, figsize=("4c","12c"),sharex=True):
        with fig1.set_panel(panel = [0,0]):
            fig1.basemap(region=[20, 400,min(min(obs_T),min(pred_Rp),min(ln_mean_parker),min(ln_mean_montalva)) -2.85,max(max(obs_T),max(pred_Rp),max(ln_mean_parker),max(ln_mean_montalva))+0.5], projection=projection_log, frame=['Wsne','xa2fg3','ya2fg3'])
            fig1.plot(x = Rp_lock,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)
            #fig1.plot(x = R_gmm_Rp,y = pred_Rp,pen='1p,red',label='This model (Global)') 
            fig1.plot(x = R_gmm_Rp,y = pred_Rp_regional,pen='1p,red',label='This model')                
            fig1.text(text="b) R@-p,lock@-", position='TL', font="8p,Helvetica,black", projection=projection_log,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.8c",box = "+gwhite+p1p")      
        with fig1.set_panel(panel = [1,0]):
            fig1.basemap(region=[20, 400,min(min(obs_T),min(pred_Rp),min(ln_mean_parker),min(ln_mean_montalva)) -2.85,max(max(obs_T),max(pred_Rp),max(ln_mean_parker),max(ln_mean_montalva))+0.5], projection=projection_log, frame=['Wsne','xa2fg3','ya2fg3'])
            fig1.plot(x = Rp,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)
            #fig1.plot(x = R_gmm_Rp,y = pred_Rp,pen='1p,red',label='This model (Global)') 
            fig1.plot(x = R_gmm_Rp,y = pred_Rp_regional,pen='1p,red',label='This model')
            fig1.text(text="c) R@-p@-", position='TL', font="8p,Helvetica,black", projection=projection_log,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.8c",box = "+gwhite+p1p") 
        with fig1.set_panel(panel = [2,0]):
            fig1.basemap(region=[20, 400,min(min(obs_T),min(pred_Rp),min(ln_mean_parker),min(ln_mean_montalva)) -2.85,max(max(obs_T),max(pred_Rp),max(ln_mean_parker),max(ln_mean_montalva))+0.5], projection=projection_log, frame=['WSne','xa2fg3','ya2fg3'])
            fig1.plot(x = Rrup,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)
            #fig1.plot(x = R_gmm_Rrup,y = pred_Rrup,pen='1p,red',label = 'This model')
            fig1.plot(x = R_gmm_Rrup,y = pred_Rrup_regional,pen='1p,red',label = 'This model')
            fig1.plot(x = R_gmm_Rrup,y = ln_mean_parker,pen='1p,blue', label='Parker et al. (2021)')
            fig1.plot(x = R_gmm_Rrup,y = pred_HQ_maule,pen='1p,black', label='Montalva et al. (2017)')
            fig1.text(text="d) R@-rup@-", position='TL', font="8p,Helvetica,black", projection=projection_log,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1') 
            fig1.text(text="PGA at V@-s30@- = 760 m/s [log units]", position = 'ML', font="10p,Helvetica,black", projection=projection_log,offset='-1.1/0.8',no_clip=True,angle=90) 
            fig1.text(text="Distance [km]", position = 'BC', font="10p,Helvetica,black", projection=projection_log,offset='0/-0.8',no_clip=True)
            fig1.legend(position="JBL+jBL+o0.1c+w3.8c",box = "+gwhite+p1p") 

fig1.savefig('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Predicciones/Maule/PGA_pred_maule_final.png',dpi=500)
fig1.savefig('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Predicciones/Maule/PGA_pred_maule_final.pdf')
        
fig1.show()

