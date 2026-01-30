import numpy as np
import pandas as pd
import pygmt
import matplotlib.pyplot as plt
from geographiclib.geodesic import Geodesic
import sys
sys.path.append('/home/jtommy/Escritorio/Respaldo/functions/')  # Añade la ruta al sistema
import functions_py
import funciones_fuente as ff
import xarray as xr


# Load slab

sam_slab = pygmt.load_dataarray('/home/jtommy/Escritorio/Respaldo/base_de_datos/slab2/grid/sam_slab2_dep_02.23.18.grd')
sam_filt = sam_slab.sel(y=slice(-40,-35))



#################

Mw = 8.8
top_center = np.array([-74, -37.8])
angulo_rotacion = 0
Nasp = 2
N_fallas = 200
std_lon = 2
std_lat = 2

# Distancia entre nodos debe ser igual, corregir funcion
fuente_ff = ff.crea_fuente(Mw,top_center,angulo_rotacion,
                           Nasp,N_fallas,std_lon=std_lon,std_lat=std_lat)

fuente = np.zeros((len(fuente_ff),4))
fuente[:,0] = np.round(fuente_ff[:,0],5)
fuente[:,1] = np.round(fuente_ff[:,1],5)
fuente[:,2] = fuente_ff[:,2]/np.max(fuente_ff[:,2])

x_vals = np.unique(fuente[:,0])
y_vals = np.unique(fuente[:,1])

x_inc = round(float(x_vals[1] - x_vals[0]), 5)
y_inc = round(float(y_vals[1] - y_vals[0]), 5)


print(x_inc,y_inc)

depth_min = 5    # km
depth_max = 40   # km


lon_norm = (fuente[:,0] - fuente[:,0].min()) / (fuente[:,0].max() - fuente[:,0].min())

depth = depth_max - (depth_max - depth_min) * lon_norm

fuente[:,3] = -depth
grid_depth = pygmt.xyz2grd(x=fuente[:,0],y=fuente[:,1],z=-fuente[:,3],region=[fuente[:,0].min()-1,fuente[:,0].max()+1,fuente[:,1].min()-1,fuente[:,1].max()+1],
                           spacing=[x_inc,y_inc])

lon_sta = -74.5
lat_sta = -32.5
depth_sta = 30
# Región y rango de profundidades
region_xy = [fuente[:,0].min(),fuente[:,0].max(),fuente[:,1].min(),fuente[:,1].max()]  # xmin, xmax, ymin, ymax
zlim = [-50, 50]  # Profundidades desde -50 km hasta 0 km

#hyp_vec = [grid_depth.x.values[3],grid_depth.y.values[5],-grid_depth.values[5,3]]
hyp_vec = [fuente[3,0],fuente[3,1],-fuente[3,3]]

x0, y0, z0 = hyp_vec[0], hyp_vec[1], hyp_vec[2]   # lon, lat, depth (km, típicamente negativa)
x1, y1, z1 = lon_sta,   lat_sta,   depth_sta
info = Geodesic.WGS84.Inverse(y1, x1, y0, x0)  # (lat, lon, lat, lon)
az = info["azi2"]                 # azimut en grados (dirección hipocentro->estación)
dist_km = info["s12"] / 1000.0      # distancia horizontal en km

dist_deg = dist_km / 111.32

dist_min = np.sqrt((111.11*(grid_depth.x-lon_sta))**2 + (111.11*(grid_depth.y-lat_sta))**2)
idx_flat = dist_min.data.argmin()
ix, iy = np.unravel_index(idx_flat, dist_min.shape)
lon_min   = grid_depth.x.isel(x=ix).item()
lat_min   = grid_depth.y.isel(y=iy).item()
depth_min = -grid_depth.isel(y=iy, x=ix).item()

Rrup_km = []

for i in range(len(fuente)):
    info_rup = Geodesic.WGS84.Inverse(lat_sta, lon_sta, fuente[i,1], fuente[i,0])
    dist_rup = info_rup["s12"] / 1000.0
    dz_rup = depth_sta - (-fuente[i,3])
    Rrup = np.sqrt(dist_rup**2 + dz_rup**2)
    Rrup_km.append(Rrup)
    
lon_min, lat_min, depth_min = fuente[np.argmin(Rrup_km),0], fuente[np.argmin(Rrup_km),1], fuente[np.argmin(Rrup_km),3]

# Componente vertical y longitud total en "grados"
dz_km = z1 - z0
dz_deg = dz_km / 111.32  # para mantener proporción
elev = np.degrees(np.arctan2(dz_km, dist_km))
L_map = info["a12"]               # grados de arco en el mapa
scale_len = 3.0                   # prueba 2–6 según se vea
L = L_map * scale_len

vec_sta = [x1, y1, z1]
vec_hyp = [x0, y0, z0]
#vec_directions = [1000,1000]
vec = np.array([vec_sta+vec_hyp])
vec_Rrup = [lon_min, lat_min, depth_min]
vec_Rrup = np.array([vec_sta+vec_Rrup])

az_plot, elev_plot,zlevel = 130,80,0
fig = pygmt.Figure()

# fig.grdview(
#     grid=grid_depth,
#     # Sets the view azimuth as 130 degrees, and the view elevation as 30
#     # degrees
#     perspective=[az_plot, elev_plot],
#     # Sets the x- and y-axis labels, and annotates the west, south, and east axes
#     frame=["xa0", "ya0", "wsne"],
#     # Sets a Mercator projection on a 15-centimeter figure
#     projection="M6c",
#     region=region_xy,
#     # Sets the height of the three-dimensional relief at 1.5 centimeters
#     zsize="7c",
# )
pygmt.makecpt(cmap='hot',series = [0,1],reverse=True)
fig.plot3d(
    x=fuente[:,0],
    y=fuente[:,1],
    z=fuente[:,3],
    style="r0.69c/0.91c", 
    fill=fuente[:,2],
    pen="black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy+zlim,
    cmap=True,
    no_clip=True
)
lon_polo = [-77.7] #-77.1
lat_polo = [-35.5]

## Legend manual
fig.plot3d(
    x=lon_polo,
    y=lat_polo,
    z=[40],
    style="u0.6c",      # círculo 0.3 cm
    fill="green",
    pen="black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
)

lon1_rot,lat1_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.4],[-35.4],np.radians(50))

fig.text(x = lon1_rot,y = lat1_rot,text="Site",  font="16p,Helvetica,black", projection="M6c",angle=50, no_clip=True,perspective=[az_plot, elev_plot,zlevel])

lon1_rot,lat1_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.4],[-35.6],np.radians(50))
fig.text(x = lon1_rot,y = lat1_rot,text="R@-RUP@-",  font="16p,Helvetica,black", projection="M6c",angle=50, no_clip=True,perspective=[az_plot, elev_plot,zlevel])
lon1_rot,lat1_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.8],[-35.7],np.radians(50))

fig.plot3d(
    data=np.array([[lon1_rot[0],lat1_rot[0],40,40,1]]),
    style="V0.4c+e+a45",      # círculo 0.3 cm
    fill="black",
    pen="1.5p,black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
)

lon2_rot,lat2_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.3],[-35.8],np.radians(50))
fig.text(x = lon2_rot,y = lat2_rot,text="R@-i = 1:Nsub@-",  font="16p,Helvetica,black", projection="M6c",angle=50, no_clip=True,perspective=[az_plot, elev_plot,zlevel])
lon2_rot,lat2_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.8],[-35.9],np.radians(50))
fig.plot3d(
    data=np.array([[lon2_rot[0],lat2_rot[0],40,40,1]]),
    style="V0.4c+e+a45",      # círculo 0.3 cm
    fill="grey",
    pen="1.5p,grey",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
)


lon3_rot,lat3_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.4],[-36.0],np.radians(50))
fig.text(x = lon3_rot,y = lat3_rot,text="R@-HYP@-",  font="16p,Helvetica,black", projection="M6c",angle=50, no_clip=True,perspective=[az_plot, elev_plot,zlevel])
lon3_rot,lat3_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.8],[-36.1],np.radians(50))
fig.plot3d(
    data=np.array([[lon3_rot[0],lat3_rot[0],40,40,1]]),
    style="V0.4c+e+a45",      # círculo 0.3 cm
    fill="red",
    pen="1.5p,red",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
)

lon4_rot,lat4_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.2],[-36.2],np.radians(50))
fig.text(x = lon4_rot,y = lat4_rot,text="Hypocenter",  font="16p,Helvetica,black", projection="M6c",angle=50, no_clip=True,perspective=[az_plot, elev_plot,zlevel])
lon4_rot,lat4_rot = functions_py.rotar_ff(lon_polo[0],lat_polo[0],[-77.7],[-36.3],np.radians(50))

fig.plot3d(
    x = [lon4_rot[0]],
    y = [lat4_rot[0]],
    z = [40],
    style="a1c",      # círculo 0.3 cm
    fill="blue",
    pen="black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
)

############ Grafico vectores individuales ##############


for i in range(len(fuente)):
    vec_Rp_i = [fuente[i,0]+0.15, fuente[i,1], fuente[i,3]]
    vec_Rp = np.array([vec_sta+vec_Rp_i])
    fig.plot3d(
        data=vec_Rp,
        style="V0.2c+s+e+a45",     # V=vector 3D; cabeza 0.5 cm; +e cabeza al final; +a apertura 45°
        pen="1p,grey",
        fill="grey",
        perspective=[az_plot, elev_plot,zlevel],
        projection="M6c", 
        zsize="7c",
        region=region_xy+zlim,
        no_clip=True                                                                                                                                                                                                                                                                                                                                                       
    )




fig.plot3d(
    x=[lon_sta],
    y=[lat_sta],
    z=[depth_sta],
    style="u0.5c",      # círculo 0.3 cm
    fill="green",
    pen="black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    frame=["xa0", "ya0","za0","wsne"],
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
)



fig.plot3d(
    x = [hyp_vec[0]],
    y = [hyp_vec[1]],
    z = [hyp_vec[2]],
    style="a1c",      # círculo 0.3 cm
    fill="blue",
    pen="black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c",
    frame=["xa0", "ya0","za0","wsne"],
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True
   
)

fig.plot3d(
    x=[hyp_vec[0], lon_sta],
    y=[hyp_vec[1], lat_sta],
    z=[hyp_vec[2], depth_sta],
    style="v0.4c+e+a45",     # flecha (v), tamaño 0.4 cm, con punta al final (+e)
    pen="2p,red",
    fill="red"
)


fig.plot3d(
    data=vec,
    style="V0.6c+s+e+a45",     # V=vector 3D; cabeza 0.5 cm; +e cabeza al final; +a apertura 45°
    pen="3p,red",
    fill="red",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c", 
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True                                                                                                                                                                                                                                                                                                                                                      
)

fig.plot3d(
    data=vec_Rrup,
    style="V0.6c+s+e+a45",     # V=vector 3D; cabeza 0.5 cm; +e cabeza al final; +a apertura 45°
    pen="3p,black",
    fill="black",
    perspective=[az_plot, elev_plot,zlevel],
    projection="M6c", 
    zsize="7c",
    region=region_xy+zlim,
    no_clip=True                                                                                                                                                                                                                                                                                                                                                   
)
#fig.legend()

fig.savefig('/home/jtommy/Escritorio/Respaldo/Paper2_v2/Figuras_paper/modelos_sinteticos_metrics.pdf')
fig.savefig('/home/jtommy/Escritorio/Respaldo/Paper2_v2/Figuras_paper/modelos_sinteticos_metrics.png',dpi=300)
fig.show()