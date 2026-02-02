import numpy as np
import pandas as pd
from datetime import datetime
from scipy.spatial import cKDTree
import pygmt
import sys
sys.path.append('/home/jtommy/Escritorio/Respaldo/functions')
import functions_py
print(functions_py.__file__)


database = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/flatfiles/Drapelaetal/Drapela_database.csv', index_col=False)
database_region  = database[['NGAsubEQID','DatabaseRegion']].drop_duplicates()
modelos_cosismicos = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_cosismicos/modelos_cosismicos.csv',index_col = False)
modelos_cosismicos = pd.merge(modelos_cosismicos, database_region, on='NGAsubEQID', how='left')
for col in ['event_lat', 'event_lon']:
    modelos_cosismicos.loc[modelos_cosismicos[col].abs() > 360, col] /= 1000
modelos_locking = pd.read_csv('/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/modelos_locking.csv',index_col = False, encoding='ISO-8859-1')
modelos_locking = modelos_locking.replace({
    r'\+AF8-': '_',
    r'\+AC0-': '-',
    r'\+ACY-': '&',
    r'\+AF8': '_'
}, regex=True)
modelos_locking.columns = modelos_locking.columns.str.replace(r'\+AF8-', '_', regex=True)
locking_models_paths = '/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_locking/'
coseismic_models_paths = '/home/jtommy/Escritorio/Respaldo/base_de_datos/Modelos_cosismicos/'


eqid = ['4000001']

#eqid = modelos_cosismicos['NGAsubEQID'].loc[modelos_cosismicos['NGAsubEQID']!= '-999.0'].unique()


#### CORRER PARA GENERAR MODELOS COSISMICOS NUEVAMENTE, DEJE LA ZORRA CON LOS DECIMALESa
for EVENT in eqid:
    print(EVENT)
    modelos_lock_eq = []
    modelos_cosismicos_eq = modelos_cosismicos.loc[modelos_cosismicos['NGAsubEQID'] == EVENT]
    if len(modelos_cosismicos_eq) > 0:
        eq_date = datetime.strptime(modelos_cosismicos_eq['event_date'].unique()[0], '%m/%d/%Y')
        year_days = 366 if (eq_date.year% 4 == 0 and (eq_date.year % 100 != 0 or eq_date.year % 400 == 0)) else 365
        region = modelos_cosismicos_eq['DatabaseRegion'].unique()[0]
        decimal_year = eq_date.year + (eq_date.timetuple().tm_yday - 1) / year_days
        longitude_cosis = np.mean(modelos_cosismicos_eq['event_lon'])
        latitude_cosis = np.mean(modelos_cosismicos_eq['event_lat'])
        modelos_locking_selected = modelos_locking.loc[modelos_locking['DatabaseRegion'] == region]
        if len(modelos_locking_selected) > 0:
            for idx in modelos_locking_selected.index:
                lon_min_lock = float(modelos_locking_selected['lon_min'][idx])
                lon_max_lock = float(modelos_locking_selected['lon_max'][idx])
                lat_min_lock = float(modelos_locking_selected['lat_min'][idx])
                lat_max_lock = float(modelos_locking_selected['lat_max'][idx])
                time_min_lock = modelos_locking_selected['time_min'][idx]
                time_max_lock = modelos_locking_selected['time_max'][idx]
                if (longitude_cosis > lon_min_lock) and (longitude_cosis< lon_max_lock) and (latitude_cosis > lat_min_lock) and (latitude_cosis< lat_max_lock) and (decimal_year>time_min_lock) and (decimal_year<time_max_lock):
                    modelos_locking_elegido = modelos_locking_selected.loc[idx]
                    modelos_lock_eq.append(modelos_locking_elegido)
                    print('Modelo de locking elegido: ',modelos_locking_elegido['model_tag'])
                    #break
                else:
                    modelos_locking_elegido = pd.DataFrame()
                    continue
            if len(modelos_locking_elegido) == 0:
                print('El terremoto no tiene modelos de locking disponible',EVENT)
                continue
            else:
                del modelos_locking_elegido
                for modelos_locking_elegido in modelos_lock_eq:
                    if modelos_locking_elegido['extension'] == '.xyz':
                        locking_model_eq = np.loadtxt(locking_models_paths+str(modelos_locking_elegido['model_tag'])+'.xyz')
                        locking_coords = locking_model_eq[:,:2]
                        tree = cKDTree(locking_coords)
                        for i in modelos_cosismicos_eq.index:
                            coseismic_model_eq =  np.loadtxt(coseismic_models_paths+str(modelos_cosismicos_eq['event_tag'].loc[i])+'.xyz')
                            coseismic_coords = coseismic_model_eq[:, :2]                
                            distancias, indices_mas_cercanos = tree.query(coseismic_coords)
                            locking_values_match = locking_model_eq[indices_mas_cercanos,2]
                            locking_model_eq_i = coseismic_model_eq.copy()
                            locking_model_eq_i[:,2] = locking_values_match
                            nombre_archivo = str(EVENT)+'_'+str(modelos_cosismicos_eq['event_tag'].loc[i])+'_'+str(modelos_locking_elegido['model_tag'])+'_lock.xyz'
                            np.savetxt(locking_models_paths+nombre_archivo, locking_model_eq_i, fmt='%.6f')
                    elif modelos_locking_elegido['extension'] == '.grd':
                        for i in modelos_cosismicos_eq.index:
                            coseismic_model_eq =  np.loadtxt(coseismic_models_paths+str(modelos_cosismicos_eq['event_tag'].loc[i])+'.xyz')
                            coseismic_coords = coseismic_model_eq[:, :2]     
                            points_df = pd.DataFrame(coseismic_coords, columns=['Longitudes', 'Latitudes'])
                            #print(points_df['Longitudes'].max())
                            if points_df['Longitudes'].max() >= 180.0:    
                                #print('-180 - 180')
                                points_df['Longitudes'] = points_df['Longitudes']-360
                            #print(points_df['Longitudes'].max())
                            locking_i = pygmt.grdtrack(grid=locking_models_paths+modelos_locking_elegido['model_tag']+'.grd',points = points_df,
                                                        newcolname = 'Locking degree',radius = True)
                            locking_i['Depth (km)'] = functions_py.get_slab_depth(points_df['Longitudes'],points_df['Latitudes'])['Depth (km)']
                            nombre_archivo = str(EVENT)+'_'+str(modelos_cosismicos_eq['event_tag'].loc[i])+'_'+str(modelos_locking_elegido['model_tag'])+'_lock.xyz'
                            np.savetxt(locking_models_paths+nombre_archivo, locking_i.to_numpy(), fmt='%.6f')
                    else:
                        print('El terremoto no tiene modelos de locking disponibles en .xyz: ',EVENT)
                        break
        else:
            print('El terremoto no tiene modelos de locking disponible')
    else:
        print('el terremoto no tiene modelos cosismicos y/o locking disponibles')
