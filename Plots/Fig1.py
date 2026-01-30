import os
import sys
import numpy as np
import pygmt
from geographiclib.geodesic import Geodesic
import funciones_fuente as ff 


# =============================================================================
# 1. PATH SETUP
# =============================================================================


script_dir = os.path.dirname(os.path.abspath(__file__))

repo_root = os.path.dirname(script_dir)

functions_path = os.path.join(repo_root, 'functions')

sys.path.append(functions_path)

sys.path.append(script_dir)

import functions_py

# Definir ruta de salida
output_folder = os.path.join(repo_root, "Figuras_paper")

# =============================================================================
# 2. SOURCE PARAMETERS AND GENERATION
# =============================================================================

Mw = 8.8
top_center = np.array([-74, -37.8])
angulo_rotacion = 0
Nasp = 2
N_fallas = 200
std_lon = 2
std_lat = 2

# Generate source
fuente_ff = ff.crea_fuente(Mw, top_center, angulo_rotacion,
                           Nasp, N_fallas, std_lon=std_lon, std_lat=std_lat)

fuente = np.zeros((len(fuente_ff), 4))
fuente[:, 0] = np.round(fuente_ff[:, 0], 5)
fuente[:, 1] = np.round(fuente_ff[:, 1], 5)
fuente[:, 2] = fuente_ff[:, 2] / np.max(fuente_ff[:, 2])

x_vals = np.unique(fuente[:, 0])
y_vals = np.unique(fuente[:, 1])

x_inc = round(float(x_vals[1] - x_vals[0]), 5)
y_inc = round(float(y_vals[1] - y_vals[0]), 5)

depth_min = 5    # km
depth_max = 40   # km

# Normalización para profundidad
lon_norm = (fuente[:, 0] - fuente[:, 0].min()) / (fuente[:, 0].max() - fuente[:, 0].min())
depth = depth_max - (depth_max - depth_min) * lon_norm
fuente[:, 3] = -depth

# Create grid 
grid_depth = pygmt.xyz2grd(
    x=fuente[:, 0],
    y=fuente[:, 1],
    z=-fuente[:, 3],
    region=[fuente[:, 0].min()-1, fuente[:, 0].max()+1, fuente[:, 1].min()-1, fuente[:, 1].max()+1],
    spacing=[x_inc, y_inc]
)

lon_sta = -74.5
lat_sta = -32.5
depth_sta = 30

region_xy = [fuente[:, 0].min(), fuente[:, 0].max(), fuente[:, 1].min(), fuente[:, 1].max()]
zlim = [-50, 50]

# Hipocentro
hyp_vec = [fuente[3, 0], fuente[3, 1], -fuente[3, 3]]
x0, y0, z0 = hyp_vec[0], hyp_vec[1], hyp_vec[2]
x1, y1, z1 = lon_sta, lat_sta, depth_sta

# Calcular Rrup mínimo (distancia a la ruptura)
dist_min = np.sqrt((111.11 * (grid_depth.x - lon_sta))**2 + (111.11 * (grid_depth.y - lat_sta))**2)
idx_flat = dist_min.data.argmin()
ix, iy = np.unravel_index(idx_flat, dist_min.shape)

Rrup_km = []
for i in range(len(fuente)):
    info_rup = Geodesic.WGS84.Inverse(lat_sta, lon_sta, fuente[i, 1], fuente[i, 0])
    dist_rup_val = info_rup["s12"] / 1000.0
    dz_rup = depth_sta - (-fuente[i, 3])
    Rrup = np.sqrt(dist_rup_val**2 + dz_rup**2)
    Rrup_km.append(Rrup)

idx_min_rrup = np.argmin(Rrup_km)
lon_min, lat_min, depth_min = fuente[idx_min_rrup, 0], fuente[idx_min_rrup, 1], fuente[idx_min_rrup, 3]

# Vectores para plotear
vec = np.array([[x1, y1, z1, x0, y0, z0]]) # Vector Estación -> Hipocentro
vec_Rrup = np.array([[x1, y1, z1, lon_min, lat_min, depth_min]]) # Vector Estación -> Punto más cercano

# =============================================================================
# 3. GRAFICAR
# =============================================================================

az_plot, elev_plot, zlevel = 130, 80, 0
fig = pygmt.Figure()

pygmt.makecpt(cmap='hot', series=[0, 1], reverse=True)

# 3.1 Plotting source as 3D scatter
fig.plot3d(
    x=fuente[:, 0], y=fuente[:, 1], z=fuente[:, 3],
    style="r0.69c/0.91c",
    fill=fuente[:, 2],
    pen="black",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c",
    zsize="7c",
    region=region_xy + zlim,
    cmap=True,
    no_clip=True
)

# 3.2 Reference point
lon_polo = [-77.7]
lat_polo = [-35.5]

# Ficticious station (Legend)
fig.plot3d(
    x=lon_polo, y=lat_polo, z=[40],
    style="u0.6c", fill="green", pen="black",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c", zsize="7c", region=region_xy + zlim, no_clip=True
)

# Textos Rotados
labels = [
    ("Site", -77.4, -35.4),
    ("R@-RUP@-", -77.4, -35.6),
    ("R@-i = 1:Nsub@-", -77.3, -35.8),
    ("R@-HYP@-", -77.4, -36.0),
    ("Hypocenter", -77.2, -36.2)
]

for label, lx, ly in labels:
    l_rot, lat_rot = functions_py.rotar_ff(lon_polo[0], lat_polo[0], [lx], [ly], np.radians(50))
    fig.text(x=l_rot, y=lat_rot, text=label, font="16p,Helvetica,black",
             projection="M6c", angle=50, no_clip=True, perspective=[az_plot, elev_plot, zlevel])

# Símbolos de la leyenda
symbols = [
    ("V0.4c+e+a45", "black", -77.8, -35.7),
    ("V0.4c+e+a45", "grey", -77.8, -35.9),
    ("V0.4c+e+a45", "red", -77.8, -36.1),
    ("a1c", "blue", -77.7, -36.3) # Hipocentro
]

for style, color, lx, ly in symbols:
    l_rot, lat_rot = functions_py.rotar_ff(lon_polo[0], lat_polo[0], [lx], [ly], np.radians(50))
    fig.plot3d(
        x=[l_rot[0]], y=[lat_rot[0]], z=[40],
        style=style, fill=color, pen=f"1.5p,{color}" if "V" in style else "black",
        perspective=[az_plot, elev_plot, zlevel],
        projection="M6c", zsize="7c", region=region_xy + zlim, no_clip=True
    )

# 3.3 Individual vectors Rp (Rayos grises)
for i in range(len(fuente)):
    vec_data = np.array([[x1, y1, z1, fuente[i, 0] + 0.15, fuente[i, 1], fuente[i, 3]]])
    
    fig.plot3d(
        data=vec_data,
        style="V0.2c+s+e+a45",
        pen="1p,grey",
        fill="grey",
        perspective=[az_plot, elev_plot, zlevel],
        projection="M6c", zsize="7c", region=region_xy + zlim, no_clip=True
    )

# 3.4 Estación (Sitio)
fig.plot3d(
    x=[lon_sta], y=[lat_sta], z=[depth_sta],
    style="u0.5c", fill="green", pen="black",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c",
    frame=["xa0", "ya0", "za0", "wsne"],
    zsize="7c", region=region_xy + zlim, no_clip=True
)

# 3.5 Hipocentro en el plano
fig.plot3d(
    x=[hyp_vec[0]], y=[hyp_vec[1]], z=[hyp_vec[2]],
    style="a1c", fill="blue", pen="black",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c", zsize="7c", region=region_xy + zlim, no_clip=True
)

# Hypocentral vector line (legend)
fig.plot3d(
    x=[hyp_vec[0], lon_sta], 
    y=[hyp_vec[1], lat_sta],
    z=[hyp_vec[2], depth_sta],
    style="v0.4c+e+a45",
    pen="2p,red", fill="red",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c", zsize="7c", region=region_xy + zlim
)

# 3.6 Hypocentral vector
fig.plot3d(
    data=vec,
    style="V0.6c+s+e+a45",
    pen="3p,red", fill="red",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c", zsize="7c", region=region_xy + zlim, no_clip=True
)

# 3.7 Rrup vector
fig.plot3d(
    data=vec_Rrup,
    style="V0.6c+s+e+a45",
    pen="3p,black", fill="black",
    perspective=[az_plot, elev_plot, zlevel],
    projection="M6c", zsize="7c", region=region_xy + zlim, no_clip=True
)


output_file = os.path.join(output_folder, 'Fig1.pdf')
fig.savefig(output_file)
print(f"Figure saved at: {output_file}")
