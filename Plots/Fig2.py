import os
import sys
import pandas as pd
import numpy as np
import pygmt

# =============================================================================
# 1. PATHS
# =============================================================================

script_dir = os.path.dirname(os.path.abspath(__file__))

repo_root = os.path.dirname(script_dir)


data_dir = os.path.join(repo_root, 'Data')
rp_dir = os.path.join(data_dir, 'Rp')


db_path = os.path.join(data_dir, 'Drapela_database.csv')
events_path = os.path.join(data_dir, 'Rasp_events_info.csv')
rp_path = os.path.join(rp_dir, 'Rp_median_values.xlsx')


output_folder = os.path.join(repo_root, "Figuras_paper")



# =============================================================================
# 2. Load Data
# =============================================================================


data = pd.read_csv(db_path)

Rp_min = pd.read_excel(rp_path, sheet_name='p = -50')
Rp_max = pd.read_excel(rp_path, sheet_name='p = 50')

data_Rp_min = data.merge(Rp_min, how='inner', on=['NGAsubEQID', 'Station_Name'])
data_Rp_max = data.merge(Rp_max, how='inner', on=['NGAsubEQID', 'Station_Name'])

M = data['Earthquake_Magnitude']
Rrup = data['ClstD_km']

events_info = pd.read_csv(events_path)


mask_neg_ev = events_info['Hypocenter_Longitude_deg'] < 0
events_info.loc[mask_neg_ev, 'Hypocenter_Longitude_deg'] += 360

mask_neg_st = data['Station_Longitude_deg'] < 0
data.loc[mask_neg_st, 'Station_Longitude_deg'] += 360


events_info = events_info[events_info['NGAsubEQID'].isin(data['NGAsubEQID'])]

# Map parameters
lon_max = 320
lon_min = 120
lat_min = -55
lat_max = 65

Mw = events_info['Earthquake_Magnitude']

# =============================================================================
# 3. Graphics
# =============================================================================

fig = pygmt.Figure()


with fig.subplot(nrows=1, ncols=1, figsize=("12c", "12c")):
    with pygmt.config(MAP_FRAME_TYPE="plain"):        
        fig.basemap(region=[lon_min, lon_max, lat_min, lat_max], projection="M?", frame=["WSne", "xa40", "ya25"])
        fig.coast(land="black", water="skyblue")
    
    pygmt.makecpt(cmap="jet", series=[7, 9.2])
    
    fig.text(text="a)", position='TL', font="14p,Helvetica,black", projection="M?", 
             fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
    
    
    fig.plot(
        x=events_info['Hypocenter_Longitude_deg'],
        y=events_info['Hypocenter_Latitude_deg'],
        size=0.001 * 2**Mw,
        fill=Mw,
        cmap=True,
        style="cc",
        pen="black"
    )
    
    
    fig.plot(
        x=data["Station_Longitude_deg"],
        y=data["Station_Latitude_deg"],
        size=0.12 * np.ones(len(data)),
        style="t",
        fill="green",
        pen="black"
    )
    
    with pygmt.config(FONT_ANNOT="13p"):
        fig.colorbar(frame="a0.5f0.25+lM@-w@-", projection="M?")    


fig.shift_origin(xshift="w+2.5c", yshift="0.5c")


with fig.subplot(nrows=4, ncols=1, figsize=("10c", "10c")):
    
    
    fig.basemap(region=[10, 1400, 6, 9.3], projection="X?l/Y?", panel=[0, 0], frame=["xag3", "yag", "Wsne"])
    fig.plot(x=Rrup, y=M, style="t0.25c", pen='1p,black', fill='blue')
    fig.text(text="b) R@-rup@-", position='TL', font="10p,Helvetica,black", projection="X?l/Y?", 
             fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
    
    
    fig.basemap(region=[10, 1400, 6, 9.3], projection="X?l/Y?", panel=[1, 0], frame=["xag3", "yag", "Wsne"])
    fig.plot(x=data_Rp_min['Rp_median_km'], y=M, style="t0.25c", pen='1p,black', fill='grey')
    fig.text(text="c) R@-p = -50@-", position='TL', font="10p,Helvetica,black", projection="X?l/Y?", 
             fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
    
    
    fig.basemap(region=[10, 1400, 6, 9.3], projection="X?l/Y?", panel=[2, 0], frame=["xag3", "yag", "Wsne"])
    fig.plot(x=data['Rasp_median_km'], y=M, style="t0.25c", pen='1p,black', fill='red')
    fig.text(text="d) R@-p = 1@-", position='TL', font="10p,Helvetica,black", projection="X?l/Y?", 
             fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')
    
    
    fig.text(text="Moment magnitude (M@-w@-)", x=6, y=10, font="11p,Helvetica,black", 
             projection="X?l/Y?", no_clip=True, angle=90)
    
    
    fig.basemap(region=[10, 1400, 6, 9.3], projection="X?l/Y?", panel=[3, 0], 
                frame=["xag3+lRupture distance (km)", "yag", "WSne"])
    fig.plot(x=data_Rp_max['Rp_median_km'], y=M, style="t0.25c", pen='1p,black', fill='green')
    fig.text(text="e) R@-p = 50@-", position='TL', font="10p,Helvetica,black", projection="X?l/Y?", 
             fill="white", pen="0.25p,black,solid", offset='0.1/-0.1')

# =============================================================================
# 4. SAVE FIGURES
# =============================================================================


pdf_path = os.path.join(output_folder, "Fig2.pdf")


fig.savefig(pdf_path)

