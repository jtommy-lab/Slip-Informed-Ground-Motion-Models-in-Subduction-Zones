
import numpy as np
import pandas as pd
import pygmt

import xarray as xr
import shapely
from shapely.geometry import LineString,Point
import random
import sys
sys.path.append('/home/jtommy/Escritorio/Respaldo/functions')  # Añade la ruta al sistema
import functions_py
from openquake.hazardlib.gsim.parker_2020 import ParkerEtAl2020SInter
from openquake.hazardlib.gsim.montalva_2017 import MontalvaEtAl2017SInter
import openquake.hazardlib.imt as IMT


eqid = '4000068'

database = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/Drapela_database.csv',index_col = None)
Fs = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/SS14_Fs_DR_database.csv',index_col = None)

coeff_Rrup = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rrup/nlmer/coeficientes_Rrup_nlmer.csv',index_col = False)


coeff_Rrup_regional = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rrup/nlmer/coeficientes_Rrup_nlmer_regional.csv',index_col = False)
coeff_Rp = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer.csv',index_col = False)
coeff_Rp_regional = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Regresiones/Rp/nlmer/coeficientes_Rp_nlmer_regional.csv',index_col = False)
#coeff_Rp = coeff_Rp.rename(columns={"c0": "Global_c0","a0": "Global_a0"})



trench_tokachi = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Trench/trench-tohoku3.xyz',delim_whitespace=True,index_col = None,names=['Longitude','Latitude'])

modelos_cosismicos = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_cosismicos/modelos_cosismicos.csv',index_col = False)
modelos_cosismicos = modelos_cosismicos.loc[modelos_cosismicos['NGAsubEQID'] == eqid]

lock_tokachi_LL2016 = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/Japan_locking_epoch2.xyz',names = ['Longitude','Latitude','Locking'],delim_whitespace=True,index_col = None)
lock_tokachi_LL2016_grd = pygmt.load_dataarray('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/lock_tokachi_LL2016.grd')
slab_tokachi = pygmt.load_dataarray('/home/jtommy/Escritorio/Respaldo/base_de_datos/slab2/grid/kur_slab2_dep_02.24.18.grd')
mask = np.isnan(slab_tokachi)
lock_tokachi_LL2016_grd = lock_tokachi_LL2016_grd.where(~mask)
#gmt triangulate /home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/Japan_locking_epoch2.xyz -G/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/lock_tokachi_LL2016.grd -R139/147/35/45 -I1m -V

Rp_sheets = pd.read_excel('/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/Rp_median_values.xlsx',sheet_name=None)
Rp_lock_sheets = pd.read_excel('/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/Rp_lock_median_values.xlsx',sheet_name=None)
Rp_lock_sheets_tokachi = pd.read_excel('/home/jtommy/Escritorio/Respaldo/base_de_datos/Rp/4000068_Rp_lock_values.xlsx',sheet_name=None)

p_opt_Rp = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/p_selected_Rp_smooth.csv',index_col = None)
p_opt_Rp_regional = pd.read_csv('/home/jtommy/Escritorio/Respaldo/Paper1_v2/residuales/nlmer/p_selected_Rp_regional_smooth.csv',index_col = None)

gmpe_parker = ParkerEtAl2020SInter(region='JP',saturation_region='JP_Pac')
gmpe_montalva = MontalvaEtAl2017SInter()
tokachi_database = database.loc[database['NGAsubEQID'] == eqid]
R_gmm_Rrup = np.logspace(np.log10(40),np.log10(700),100)
R_gmm_Rp = np.logspace(np.log10(60),np.log10(900),100)
Rref_gmm = functions_py.get_R_Rref_interface(tokachi_database['Earthquake_Magnitude'].unique(),R_gmm_Rp)[1]
input_df = pd.DataFrame({"mag": len(R_gmm_Rrup)*[tokachi_database['Earthquake_Magnitude'].unique()], "rrup" : R_gmm_Rp, 
                            "vs30": len(R_gmm_Rp)*[760], "region": len(R_gmm_Rp)*['JP'], "saturation_region": len(R_gmm_Rp)*['JP_Pac'],"backarc":len(R_gmm_Rp)*[False]})
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

filas_prediction = []

periods_plot = [0]
#gmt img2grd /home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/topo_6.2.img -E -G/home/jtommy/Escritorio/Respaldo/Paper1_v2/Predicciones/Tokachi/tokachi.grd -R139/146/37/44 -I0.5m -T1 -N1 -V
#gmt grdgradient '/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/tokachi.nc' -G'/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/tokachi.int' -A225 -Nt


pygmt.config(FONT="8p,Helvetica,black")
R_gmm_v2 = functions_py.get_R_Rref_interface(tokachi_database['Earthquake_Magnitude'].unique(),R_gmm_Rrup)[0]
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
    # if p_opt_T.is_integer() == True:
    #     p_opt_T = int(p_opt_T)
    coeff_Rp_regional_T = coeff_Rp_regional.loc[(coeff_Rp_regional['Period'] == T) & (coeff_Rp_regional['p_value'] == p_opt_regional_T) ]
    coeff_Rp_T = coeff_Rp.loc[(coeff_Rp['Period'] == T) & (coeff_Rp['p_value'] == p_opt_T) ]
    coeff_Rrup_T = coeff_Rrup.loc[coeff_Rrup['Period'] == T]
    coeff_Rrup_regional_T = coeff_Rrup_regional.loc[coeff_Rrup_regional['Period'] == T]       
    tokachi_df = functions_py.get_event_period_df(database,Fs,eqid,T,coeff_Rp,coeff_Rrup,p_value=p_opt_T,Rp_sheets=Rp_sheets,Rp_sheets_lock=Rp_lock_sheets,region = 'Global',Rref = Rref_gmm)    
    tokachi_df = tokachi_df.dropna(subset='log_obs_rock')
    obs_T = tokachi_df['log_obs_rock']
    mag = np.array(len(R_gmm_Rrup)*[tokachi_df['Earthquake_Magnitude'].unique()[0]])
    pred_Rrup = functions_py.DR_pred_GMM_Rp(R_gmm_v2,mag,coeff_Rrup_T,Rref = Rref_gmm)
    pred_Rrup_regional = functions_py.DR_pred_GMM_Rp(R_gmm_v2,mag,coeff_Rrup_regional_T,Rref = Rref_gmm,region = 'Japan')
    pred_Rp = functions_py.DR_pred_GMM_Rp(R_gmm_Rp,mag,coeff_Rp_T,Rref = Rref_gmm)
    pred_Rp_regional = functions_py.DR_pred_GMM_Rp(R_gmm_Rp,mag,coeff_Rp_regional_T,Rref = Rref_gmm, region = 'Japan')
    Rp = tokachi_df['Rp']
    Rp_lock = tokachi_df['Rp_lock']
    Rp_lock_tokachi = Rp_lock_sheets_tokachi['p = '+str(p_opt_T)]
    Rp_lock_LL2016 = Rp_lock_tokachi.loc[:,Rp_lock_tokachi.columns.str.startswith('s2003')].median(axis=1).values
    R = tokachi_df['R']
    Rrup = tokachi_df['Rrup']
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
    ln_mean_parker = ln_mean_parker.reshape((len(mag),))
    ln_mean_montalva = ln_mean_montalva.reshape((len(mag),))
    # sig = sig.reshape((len(mag),))
    # tau = tau.reshape((len(mag),))
    # phi = phi.reshape((len(mag),))
    fig1 = pygmt.Figure()
    projection_log = "X?l/Y?"
    projection_geo = "M5.5c"
    colores = [f"#{random.randint(0, 255):02x}{random.randint(0, 255):02x}{random.randint(0, 255):02x}" for _ in range(len(modelos_cosismicos))]
    cont = 0    
    with fig1.subplot(nrows = 1, ncols = 1, figsize=("4c","12c"),sharex=True):
        with fig1.set_panel(panel = [0,0]):
            fig1.basemap(region=[139, 146,37,44], projection=projection_geo, frame=['WSne'])        
            pygmt.makecpt(cmap='/home/jtommy/Escritorio/graficos_GMT/paletas/grayscale02.cpt',series = [-8000,3000,10])
            ### GENERAR ARCHIVO DE ILUMINACION ###
            fig1.grdimage(grid ='/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/tokachi.nc',
                          shading = '/home/jtommy/Escritorio/Respaldo/base_de_datos/Topobatimetrias/tokachi.int',projection = projection_geo,cmap=True )
            pygmt.makecpt(cmap='hot',series = [0,1],reverse=True)
            # Crear una copia temporal de la grilla con valores menores al umbral como NaN
            
            fig1.grdimage(grid =lock_tokachi_LL2016_grd,projection = projection_geo, cmap=True,nan_transparent=True )
            #fig1.plot(x=lock_tokachi_LL2016['Longitude'], y=lock_tokachi_LL2016['Latitude'], fill=lock_tokachi_LL2016['Locking'],cmap=True,style='t0.1c',projection = projection_geo,pen="0.1p,black")
            fig1.plot(x=tokachi_df['Station_Longitude_deg'], y=tokachi_df['Station_Latitude_deg'], fill='green',style='t0.2c',projection = projection_geo,pen="0.1p,black")
            fig1.plot(x=trench_tokachi['Longitude'], y=trench_tokachi['Latitude'], style="f0.5i/0.1i+r+t",fill='black', pen="1.25p,black",projection = projection_geo)    
            fig1.coast(shorelines=True,projection = projection_geo)
            fig1.text(text="a) LM16", position='TL', font="8p,Helvetica,black", projection=projection_geo,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1') 
            pygmt.config(FONT="10p,Helvetica,black")
            fig1.colorbar(frame="xa.25f.25+lLocking degree",position="g139/35.7+w5.5c/0.4c+h",projection=projection_geo,region = [139, 146,37,44])  
            pygmt.config(FONT="8p,Helvetica,black")         
            for i in modelos_cosismicos['event_tag']:
                print(i)
                model_i = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_cosismicos/'+str(i)+'.xyz',delim_whitespace=True,index_col = None,names=['Longitude','Latitude','Slip','Depth'])
                model_i = model_i[['Longitude','Latitude','Slip']]
                model_i['Slip'] = model_i['Slip']/model_i['Slip'].max()
                resolution = np.sqrt((model_i['Longitude'].loc[2]-model_i['Longitude'].loc[1])**2+(model_i['Latitude'].loc[2]-model_i['Latitude'].loc[1])**2)
                resolution_m = 3*np.ceil(resolution*60)
                slip_th30 = 0.8
                max_slip = 0.05
                region = [min(model_i['Longitude'])-1, max(model_i['Longitude'])+1, min(model_i['Latitude'])-1, max(model_i['Latitude'])+1] 
                model_limits = [min(model_i['Longitude']), max(model_i['Longitude']), min(model_i['Latitude']), max(model_i['Latitude'])] 
                model_grid = pygmt.nearneighbor(data=model_i,spacing="1m",region=model_limits,search_radius=str(resolution_m)+'m',sectors="6")
                slab_model_grid = slab_tokachi.interp(y=model_grid.lat,x=model_grid.lon)
                model_grid = model_grid.where(slab_model_grid.notnull())
                #fig1.grdimage(grid =model_grid,projection = projection_geo, cmap=True ) 
                fig1.grdcontour(grid = model_grid,limit = [0,slip_th30],projection = projection_geo,annotation = "n",levels=slip_th30,pen='1p,'+colores[cont])
                #fig1.grdcontour(grid = grid,limit = [0,slip_th50],projection = projection_geo,annotation = "n",levels=slip_th50,pen='0.4p,blue')
                #fig1.grdcontour(grid = grid,limit = [0,max_slip],projection = projection_geo,annotation = "n",levels=max_slip,pen='1.8p,darkgrey') 
                fig1.plot(x=trench_tokachi['Longitude'], y=trench_tokachi['Latitude'], style="f0.5i/0.1i+r+t",fill='black', pen="1.25p,black",projection = projection_geo)    
                cont = cont+1            
    fig1.shift_origin(xshift="w+3.8c",yshift="-2.8c")    
    with fig1.subplot(nrows = 3, ncols = 1, figsize=("4c","12c"),sharey=True):
        with fig1.set_panel(panel = [0,0]):
            fig1.basemap(region=[40, 1000,min(min(obs_T),min(pred_Rp),min(ln_mean_parker),min(ln_mean_montalva)) -0.5,max(max(obs_T),max(pred_Rp),max(ln_mean_parker),max(ln_mean_montalva))+1.8], projection=projection_log, frame=['Wsne','xa2fg3','ya2fg3'])
            #fig1.plot(x = Rp_lock_LL2011,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)  
            #fig1.plot(x = np.sort(Rp_lock_LL2011),y = pred_LL2011,pen='1p,red',label='This model - LM11')
            fig1.plot(x = Rp_lock_LL2016,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)      
            #fig1.plot(x = R_gmm_Rp,y = pred_Rp,pen='1p,red',label='This model (Global)') 
            fig1.plot(x = R_gmm_Rp,y = pred_Rp_regional,pen='1p,red',label='This model')    
            fig1.text(text="b) R@-p,lock@-", position='TL', font="8p,Helvetica,black", projection=projection_log,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1') 
            fig1.legend(position="JBL+jBL+o0.1c+w2.5c",box = "+gwhite+p1p")        
        with fig1.set_panel(panel = [1,0]):
            fig1.basemap(region=[40, 1000,min(min(obs_T),min(pred_Rp),min(ln_mean_parker),min(ln_mean_montalva)) -0.5,max(max(obs_T),max(pred_Rp),max(ln_mean_parker),max(ln_mean_montalva))+1.8], projection=projection_log, frame=['Wsne','xa2fg3','ya2fg3'])
            fig1.plot(x = Rp,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)
            #fig1.plot(x = R_gmm_Rp,y = pred_Rp,pen='1p,red',label='This model (Global)')  
            fig1.plot(x = R_gmm_Rp,y = pred_Rp_regional,pen='1p,red',label='This model')    
            fig1.text(text="c) R@-p@-", position='TL', font="8p,Helvetica,black", projection=projection_log,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1')
            fig1.legend(position="JBL+jBL+o0.1c+w2.5c",box = "+gwhite+p1p")             
        with fig1.set_panel(panel = [2,0]):
            fig1.basemap(region=[40, 1000,min(min(obs_T),min(pred_Rp),min(ln_mean_parker),min(ln_mean_montalva)) -0.5,max(max(obs_T),max(pred_Rp),max(ln_mean_parker),max(ln_mean_montalva))+1.8], projection=projection_log, frame=['WSne','xa2fg3','ya2fg3'])
            fig1.plot(x = Rrup,y = obs_T,style='c0.1c',fill='gray',pen='0.1p,black',intensity = 0.3)
            #fig1.plot(x = R_gmm_Rrup,y = pred_Rrup,pen='1p,red',label = 'This model (Global)')
            fig1.plot(x = R_gmm_Rrup,y = pred_Rrup_regional,pen='1p,red',label = 'This model')
            fig1.plot(x = R_gmm_Rrup,y = ln_mean_parker,pen='1p,blue', label='Parker et al. (2021)')
            #fig1.plot(x = R_gmm,y = ln_mean_montalva,pen='1p,blue', label='Montalva et al. (2017)')
            fig1.text(text="PGA at V@-s30@- = 760 m/s [log units]", position = 'ML', font="10p,Helvetica,black", projection=projection_log,offset='-1.1/0.8',no_clip=True,angle=90)
            fig1.text(text="d) R@-rup@-", position='TL', font="8p,Helvetica,black", projection=projection_log,fill="white",pen = "0.25p,black,solid",offset='0.1/-0.1') 
            fig1.legend(position="JBL+jBL+o0.1c+w3.5c",box = "+gwhite+p1p") 
            fig1.text(text="Distance [km]", position = 'BC', font="10p,Helvetica,black", projection=projection_log,offset='0/-0.9',no_clip=True)
        
fig1.show()

fig1.savefig('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Predicciones/Tokachi/PGA_pred_tokachi_final.png',dpi=300)
fig1.savefig('/home/jtommy/Escritorio/Respaldo/Paper1_v2/Predicciones/Tokachi/PGA_pred_tokachi_final.pdf')





##[-1,0,0.5,0.75,3,10]
distance = np.logspace(0,3.2,3000)




fig1 = pygmt.Figure()
projection = "X?/Y?"
projection_log = "X?l/Y?"
projection_geo = "M?"







with fig1.subplot(nrows = 3, ncols = 1, figsize=("14c","12c"),sharex=True):
    with fig1.set_panel(panel = [0,0]):
        fig1.basemap(region=[30, 1200,min(obs_T)-2,max(obs_T)+2], projection=projection_log, frame=['Wsne+t'+str(column_obs),'xa2fg3','ya2fg3'])
        fig1.plot(x = Rehd_lock,y = obs_T,style='c0.1c',fill='red',pen='0.1p,black')
        fig1.plot(x = np.sort(Rehd_lock),y = pred_Rehd_lock_T,pen='1p,blue')        
    with fig1.set_panel(panel = [1,0]):
        fig1.basemap(region=[30, 1200,min(obs_T)-2,max(obs_T)+2], projection=projection_log, frame=['Wsne','xa2fg3','ya2fg3'])
        fig1.plot(x = Rehd,y = obs_T,style='c0.1c',fill='red',pen='0.1p,black')
        fig1.plot(x = np.sort(Rehd),y = pred_Rehd_T,pen='1p,black')
    with fig1.set_panel(panel = [2,0]):
        fig1.basemap(region=[30, 1200,min(obs_T)-2,max(obs_T)+2], projection=projection_log, frame=['WSne','xa2fg3','ya2fg3'])
        fig1.plot(x = R,y = obs_T,style='c0.1c',fill='red',pen='0.1p,black')
        fig1.plot(x = np.sort(R),y = pred_Rrup_T,pen='1p,green')

            #fig1.text(text="R@-rup@", x=35, y=max(obs_Fs_df_T['logObs_rock'])+0.8, font="10p,Helvetica,black", projection=projection_log,no_clip = True,angle = 0,pen = "0.25p,black,solid")
        # with fig1.set_panel(panel = [1,0]):
        #     fig1.basemap(region=[30, 800,min_lim_obs,max_lim_obs], projection=projection_log, frame=['WSne','xa2fg3','ya2fg3'])
        #     fig1.plot(x = obs_Fs_df_T['Rasp'],y = obs_Fs_df_T['logObs_rock'],style='c0.1c',fill=obs_Fs_df_T['Magnitude'], cmap = True,pen='0.1p,black')
        #     fig1.text(text="R@-asp@", x=35, y=max(obs_Fs_df_T['logObs_rock'])+0.8, font="10p,Helvetica,black", projection=projection_log,no_clip = True,angle = 0,pen = "0.25p,black,solid") 
        #     fig1.text(text="Sa at Vs30 = 760 m/s [log units]", x=-40, y=3.5, font="11p,Helvetica,black", projection=projection,no_clip = True,angle = 90)   
        #     fig1.text(text="Rupture distance [km]", x=150, y=-9.5, font="11p,Helvetica,black", projection=projection_log,no_clip = True,angle = 0)             
        #     fig1.colorbar(frame=["xa.4f.1+lMoment magnitude"], position="JMR+o0.5c/2.5c+w8c")