import os
import sys
import numpy as np
import pandas as pd
import pygmt

# =============================================================================
# 1. PATHS
# =============================================================================

try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.getcwd()


repo_root = os.path.dirname(script_dir)

# Carpetas Clave
data_dir = os.path.join(repo_root, 'Data')
results_dir = os.path.join(repo_root, 'Regressions', 'Results')
output_folder = os.path.join(repo_root, "Figuras_paper") 

if not os.path.exists(output_folder):
    os.makedirs(output_folder) 
    print(f"Carpeta creada: {output_folder}")


# =============================================================================
# 2. LOAD FILES
# =============================================================================

files = {
    'database': os.path.join(data_dir, 'Drapela_database.csv'),
    'dBe_Rp': os.path.join(results_dir, 'dBe_residual_Rp.csv'),
    'dBe_Rp_unit': os.path.join(results_dir, 'dBe_residual_Rp_unit.csv'),
    'dBe_Rrup': os.path.join(results_dir, 'dBe_residual_Rrup.csv'),
    'dWe_Rp': os.path.join(results_dir, 'dWe_residual_Rp.csv'),
    'dWe_Rp_unit': os.path.join(results_dir, 'dWe_residual_Rp_unit.csv'),
    'dWe_Rrup': os.path.join(results_dir, 'dWe_residual_Rrup.csv'),
    'Stats_Rp': os.path.join(results_dir, 'Stats_residual_Rp.csv'),
    'Stats_Rp_unit': os.path.join(results_dir, 'Stats_residual_Rp_unit.csv'),
    'Stats_Rrup': os.path.join(results_dir, 'Stats_residual_Rrup.csv'),
    'p_unit': os.path.join(results_dir, 'p_selected_Rp_unit.csv'),
    'p_values': os.path.join(results_dir, 'p_selected_Rp.csv')
}

def safe_read(path):
    if os.path.exists(path):
        return pd.read_csv(path, index_col=False)
    return pd.DataFrame()

# Cargar DataFrames
database = safe_read(files['database'])
dBe_residual_Rp = safe_read(files['dBe_Rp'])
dBe_residual_Rrup = safe_read(files['dBe_Rrup'])
dBe_residual_Rp_unit = safe_read(files['dBe_Rp_unit'])
dWe_residual_Rrup = safe_read(files['dWe_Rrup'])
dWe_residual_Rp = safe_read(files['dWe_Rp'])
dWe_residual_Rp_unit = safe_read(files['dWe_Rp_unit'])
Stats_residual_Rrup = safe_read(files['Stats_Rrup'])
Stats_residual_Rp = safe_read(files['Stats_Rp'])
Stats_residual_Rp_unit = safe_read(files['Stats_Rp_unit'])
p_Rp = safe_read(files['p_values'])
p_Rp_unit = safe_read(files['p_unit'])


zero_line = np.zeros((1000, 2))
zero_line[:, 0] = np.linspace(0.1, 2000, len(zero_line))
periods_plot = [-1, 0, 0.1, 0.5, 1, 5]
labels = ["a) PGV", "b) PGA", "c) T = 0.1", "d) T = 0.5", "e) T = 1", "d) T = 5"]
projection = "X?/Y?"
projection_log = "X?l/Y?"

# Helper para bines
def get_binned_stats(df, x_col, y_col, bins, stat='mean'):
    df['bin'] = pd.cut(df[x_col], bins=bins)
    df['bin_mid'] = df['bin'].apply(lambda x: (x.left + x.right) / 2 if pd.notnull(x) else None)
    grouped = df.groupby('bin_mid')[y_col]
    
    if stat == 'mean':
        y_val = grouped.mean()
        y_err = grouped.std() / np.sqrt(grouped.count())
    elif stat == 'median':
        y_val = grouped.median()
        y_err = 0 # Placeholder for boxplots
        
    return pd.DataFrame({
        'x': y_val.index.values.to_numpy(),
        'y': y_val.values,
        'xerr': 0,
        'yerr': y_err
    })

# =============================================================================
# 4. FIG S2 (dBe Rp)
# =============================================================================

if not dBe_residual_Rp.empty:
    print("Generando FigS2 (dBe Rp)...")
    fig1 = pygmt.Figure()
    cont = 0
    M_bin = [6.5, 7.5, 8.1, 9.3]
    
    with fig1.subplot(nrows=2, ncols=3, figsize=("12c", "8c"), sharey='l', sharex='b'):
        row, col = 0, 0
        for T in periods_plot:
            df_T = dBe_residual_Rp[dBe_residual_Rp['Period'] == T].copy()
            stats = get_binned_stats(df_T, 'M', 'dBe', M_bin)
            
            with fig1.set_panel(panel=[row, col]):
                label_pts = "δB@-e@-+N2" if cont == 0 else None
                label_err = "Binned means" if cont == 0 else None
                
                fig1.basemap(region=[6.5, 9.5, -1.6, 1.6], projection=projection, frame=["ya.5f.25g.5", "xa1f.5g1"])
                fig1.text(text=labels[cont], position='TL', offset="0.1/-0.1", font="8p,Helvetica,black", fill="white", pen="0.25p,black,solid")
                fig1.plot(x=zero_line[:,0], y=zero_line[:,1], pen='1.2p,black,-')
                fig1.plot(x=df_T['M'], y=df_T['dBe'], style='t0.20c', fill='blue', pen="0.25p,black,solid", intensity=0.8, label=label_pts)
                fig1.plot(data=stats, error_bar='+p0.75p,blue', style='t0.25c', fill='blue', pen="0.25p,black,solid", label=label_err)
                
                if row == 1 and col == 1:
                    fig1.text(text="Moment magnitude", x=8, y=-2.2, font="10p,Helvetica,black", no_clip=True)
                if row == 0 and col == 0:
                    if cont == 0: fig1.legend(position="JBR+jBR+o-5c/-0.6c+w7.5c")
                    fig1.text(text="R@-P@--derived GMM δB@-e@-", x=5.5, y=-1.8, font="10p,Helvetica,black", no_clip=True, angle=90)
            
            col += 1
            if col > 2: col, row = 0, row + 1
            cont += 1
            
    fig1.savefig(os.path.join(output_folder, 'FigS2.pdf'))

# =============================================================================
# 5. FIG S1 (dBe Rrup)
# =============================================================================

if not dBe_residual_Rrup.empty:
    print("Generando FigS1 (dBe Rrup)...")
    fig1 = pygmt.Figure()
    cont = 0
    
    with fig1.subplot(nrows=2, ncols=3, figsize=("12c", "8c"), sharey='l', sharex='b'):
        row, col = 0, 0
        for T in periods_plot:
            df_T = dBe_residual_Rrup[dBe_residual_Rrup['Period'] == T].copy()
            stats = get_binned_stats(df_T, 'M', 'dBe', M_bin)
            
            with fig1.set_panel(panel=[row, col]):
                label_pts = "δB@-e@-+N2" if cont == 0 else None
                label_err = "Binned means" if cont == 0 else None
                
                fig1.basemap(region=[6.5, 9.5, -1.6, 1.6], projection=projection, frame=["ya.5f.25g.5", "xa1f.5g1"])
                fig1.text(text=labels[cont], position='TL', offset="0.1/-0.1", font="8p,Helvetica,black", fill="white", pen="0.25p,black,solid")
                fig1.plot(x=zero_line[:,0], y=zero_line[:,1], pen='1.2p,black,-')
                fig1.plot(x=df_T['M'], y=df_T['dBe'], style='c0.15c', fill='red', pen="0.2p,black,solid", label=label_pts, intensity=0.8)
                fig1.plot(data=stats, error_bar='+p0.75p,red', style='c0.15c', fill='red', label=label_err)

                if row == 1 and col == 1:
                    fig1.text(text="Moment magnitude", x=8, y=-2.2, font="10p,Helvetica,black", no_clip=True)
                if row == 0 and col == 0:
                    if cont == 0: fig1.legend(position="JBR+jBR+o-5c/-0.6c+w7.5c")
                    fig1.text(text="R@-RUP@--derived GMM δB@-e@-", x=5.5, y=-1.8, font="10p,Helvetica,black", no_clip=True, angle=90)
            
            col += 1
            if col > 2: col, row = 0, row + 1
            cont += 1
            
    fig1.savefig(os.path.join(output_folder, 'FigS1.pdf'))

# =============================================================================
# 6. FIG 4 (Boxplots Comparison)
# =============================================================================

if not dBe_residual_Rrup.empty and not dBe_residual_Rp.empty:
    print("Generando Fig4 (Boxplots)...")
    fig1 = pygmt.Figure()
    cont = 0
    
    def get_quantiles(df, col='dBe'):
        grouped = df.groupby('bin_mid')[col]
        q = grouped.quantile([0, 0.25, 0.75, 1]).unstack()
        return pd.DataFrame({
            'x': grouped.median().index.values.to_numpy(),
            'y': grouped.median().values,
            'y_q0': q[0.0].values, 'y_q0.25': q[0.25].values,
            'y_q0.75': q[0.75].values, 'y_q1': q[1.0].values
        })

    with fig1.subplot(nrows=2, ncols=3, figsize=("12c", "8c"), sharey='l', sharex='b'):
        row, col = 0, 0
        for T in periods_plot:
            df_Rrup = dBe_residual_Rrup[dBe_residual_Rrup['Period'] == T].copy()
            df_Rp = dBe_residual_Rp[dBe_residual_Rp['Period'] == T].copy()
            
            # Asignar bines antes de calcular
            df_Rrup['bin'] = pd.cut(df_Rrup['M'], bins=M_bin)
            df_Rrup['bin_mid'] = df_Rrup['bin'].apply(lambda x: (x.left + x.right)/2)
            df_Rp['bin'] = pd.cut(df_Rp['M'], bins=M_bin)
            df_Rp['bin_mid'] = df_Rp['bin'].apply(lambda x: (x.left + x.right)/2)

            stats_Rrup = get_quantiles(df_Rrup)
            stats_Rrup['x'] -= 0.1
            stats_Rp = get_quantiles(df_Rp)
            stats_Rp['x'] += 0.2

            with fig1.set_panel(panel=[row, col]):
                label_Rrup = "R@-RUP@-+N2" if cont == 0 else None
                label_Rp = "R@-P@-" if cont == 0 else None
                
                fig1.basemap(region=[6.5, 9.3, -1.6, 1.6], projection=projection, frame=["ya.5f.25g.5", "xa0f.5g1"])
                fig1.text(text=labels[cont], position='TL', offset="0.1/-0.1", font="8p,Helvetica,black", fill="white", pen="0.25p,black,solid")
                
                fig1.plot(data=stats_Rrup, error_bar='Y+p0.75p', fill='red')
                fig1.plot(data=stats_Rp, error_bar='Y+p0.75p', fill='blue')
                
                # Lineas dummy para leyenda
                fig1.plot(x=[10,10], y=[10,10], pen='6p,red', label=label_Rrup)
                fig1.plot(x=[10,10], y=[10,10], pen='6p,blue', label=label_Rp)

                if row == 1:
                    fig1.text(x=[7,8,9], y=[-1.75]*3, text=['(6.5-7.5]', '(7.5-8.1]', '(8.1-9.3]'], font="8p,Helvetica,black", no_clip=True, justify='CM')
                if row == 0 and col == 0 and cont == 0:
                    fig1.legend(position="JBR+jBR+o-6.2c/-0.6c+w6c")
                if row == 1 and col == 1:
                    fig1.text(text="Moment magnitude", x=7.8, y=-2.2, font="10p,Helvetica,black", no_clip=True)
                if row == 0 and col == 0:
                    fig1.text(text="Between-event residual", x=5.5, y=-1.8, font="10p,Helvetica,black", no_clip=True, angle=90)
            
            col += 1
            if col > 2: col, row = 0, row + 1
            cont += 1
            
    fig1.savefig(os.path.join(output_folder, 'Fig4.pdf'))

# =============================================================================
# 7. FIG S3 (dWe)
# =============================================================================

if not dWe_residual_Rrup.empty:
    print("Generando FigS3 (dWe)...")
    fig2 = pygmt.Figure()
    cont = 0
    nbins = 5
    
    # Pre-calcular bines globales
    Rrup_min, Rrup_max = dWe_residual_Rrup['R'].min(), dWe_residual_Rrup['R'].max()
    Rrup_bins_grid = np.logspace(np.log10(Rrup_min), np.log10(Rrup_max), nbins+1)
    
    Rp_min, Rp_max = dWe_residual_Rp['Rp'].min(), dWe_residual_Rp['Rp'].max()
    Rp_bins_grid = np.logspace(np.log10(Rp_min), np.log10(Rp_max), nbins+1)

    with fig2.subplot(nrows=2, ncols=3, figsize=("12c", "8c"), sharex='b', sharey='l'):
        row, col = 0, 0
        for T in periods_plot:
            df_Rrup = dWe_residual_Rrup[dWe_residual_Rrup['Period'] == T].copy()
            df_Rp = dWe_residual_Rp[dWe_residual_Rp['Period'] == T].copy()
            
            stats_Rrup = get_binned_stats(df_Rrup, 'Rrup', 'dWe', Rrup_bins_grid)
            stats_Rp = get_binned_stats(df_Rp, 'Rp', 'dWe', Rp_bins_grid)

            with fig2.set_panel(panel=[row, col]):
                label_Rrup = "R@-RUP@-+N2" if cont == 0 else None
                label_Rp = "R@-P@" if cont == 0 else None
                
                fig2.basemap(region=[10, 1400, -3, 3], projection=projection_log, frame=['xa1f1g3', 'ya1f.5g1'])
                fig2.text(text=labels[cont], position='TL', offset="0.1/-0.1", font="8p,Helvetica,black", fill="white", pen="0.25p,black,solid")
                
                fig2.plot(x=df_Rrup['R'], y=df_Rrup['dWe'], style='c0.08c', fill='red', intensity=0.8)
                fig2.plot(x=df_Rp['Rp'], y=df_Rp['dWe'], style='t0.08c', fill='black', intensity=0.8)
                fig2.plot(x=zero_line[:,0], y=zero_line[:,1], pen='1.2p,black')
                
                fig2.plot(data=stats_Rrup, error_bar='+p0.75p,red', style='c0.15c', fill='red', label=label_Rrup)
                fig2.plot(data=stats_Rp, error_bar='+p0.75p,blue', style='t0.25c', fill='blue', label=label_Rp, pen="0.25p,black,solid")

                if row == 0 and col == 0 and cont == 0:
                    fig2.legend(position="JBR+jBR+o-7.5c/-0.6c+w7.5c")
                if col == 0 and row == 1:
                    fig2.text(text="Within-event residual", x=2.8, y=3, font="10p,Helvetica,black", no_clip=True, angle=90)
                if col == 1 and row == 1:
                    fig2.text(text="Distance [km]", x=100, y=-4.3, font="10p,Helvetica,black", no_clip=True)
            
            col += 1
            if col > 2: col, row = 0, row + 1
            cont += 1
            
    fig2.savefig(os.path.join(output_folder, 'FigS3.pdf'))

# =============================================================================
# 8. FIG S4 (dWe Difference)
# =============================================================================

if not dWe_residual_Rrup.empty:
    print("Generando FigS4 (dWe Diff)...")
    fig2 = pygmt.Figure()
    cont = 0
    
    with fig2.subplot(nrows=2, ncols=3, figsize=("12c", "8c"), sharex='b', sharey='l'):
        row, col = 0, 0
        for T in periods_plot:
            df_Rrup = dWe_residual_Rrup[dWe_residual_Rrup['Period'] == T].copy()
            df_Rp = dWe_residual_Rp[dWe_residual_Rp['Period'] == T].copy()
            
            df_Rrup['bin'] = pd.cut(df_Rrup['Rrup'], bins=Rrup_bins_grid)
            df_Rrup['bin_mid'] = df_Rrup['bin'].apply(lambda x: (x.left + x.right)/2)
            df_Rp['bin'] = pd.cut(df_Rp['Rp'], bins=Rp_bins_grid)
            df_Rp['bin_mid'] = df_Rp['bin'].apply(lambda x: (x.left + x.right)/2)

            stats_Rrup = get_quantiles(df_Rrup, col='dWe')
            stats_Rp = get_quantiles(df_Rp, col='dWe')
            
            # Shifts
            stats_Rp['x'] = stats_Rp['x'] * 1.2 # Log shift hack visual

            with fig2.set_panel(panel=[row, col]):
                label_Rrup = "R@-RUP@-+N2" if cont == 0 else None
                label_Rp = "R@-P@-" if cont == 0 else None
                
                fig2.basemap(region=[10, 1400, -3.5, 3.5], projection=projection_log, frame=['xa1f1g2', 'ya1f.5g1'])
                fig2.text(text=labels[cont], position='TL', offset="0.1/-0.1", font="8p,Helvetica,black", fill="white", pen="0.25p,black,solid")
                
                fig2.plot(data=stats_Rrup, error_bar='Y+p0.75p', fill='red')
                fig2.plot(data=stats_Rp, error_bar='Y+p0.75p', fill='blue')
                
                # Lineas leyenda
                fig2.plot(x=[10,10], y=[10,10], pen='6p,red', label=label_Rrup)
                fig2.plot(x=[10,10], y=[10,10], pen='6p,blue', label=label_Rp)

                if row == 0 and col == 0 and cont == 0:
                    fig2.legend(position="JBR+jBR+o-6.2c/-0.6c+w6c")
                if col == 0 and row == 1:
                    fig2.text(text="Within-event residual", x=2.8, y=3.5, font="10p,Helvetica,black", no_clip=True, angle=90)
                if col == 1 and row == 1:
                    fig2.text(text="Distance [km]", x=100, y=-4.95, font="10p,Helvetica,black", no_clip=True)
            
            col += 1
            if col > 2: col, row = 0, row + 1
            cont += 1
            
    fig2.savefig(os.path.join(output_folder, 'FigS4.pdf'))

# =============================================================================
# 9. FIG 5 (Statistics Comparison)
# =============================================================================

if not Stats_residual_Rp.empty and not Stats_residual_Rrup.empty:
    print("Generando Fig5 (Stats)...")
    
    periods_all = Stats_residual_Rrup['Period'].unique()
    # Calcular diferencias porcentuales
    diff_df = pd.DataFrame({'Period': periods_all})
    
    # Vectorized calculation
    for col in ['sigma', 'tau', 'phi']:
        s_rup = Stats_residual_Rrup.set_index('Period')[col]
        s_rp = Stats_residual_Rp.set_index('Period')[col]
        diff = 100 * (np.exp(s_rup - s_rp) - 1)
        diff_df[col] = diff.values

    # Encontrar máximos
    max_sigma_idx = diff_df['sigma'].idxmax()
    max_sigma = diff_df.iloc[max_sigma_idx]
    
    max_tau_idx = diff_df['tau'].idxmax()
    max_tau = diff_df.iloc[max_tau_idx]
    
    max_phi_idx = diff_df['phi'].idxmax()
    max_phi = diff_df.iloc[max_phi_idx]

    # Plotting
    fig3 = pygmt.Figure()
    proj_lin = "X2.1c"
    proj_log_stat = "X?l/Y?"
    
    # Filtrar datos para plot
    def get_arrays(df, col):
        return df['Period'].values, df[col].values

    with fig3.subplot(nrows=3, ncols=2, figsize=("12c", "8c"), sharey=True, margins=["0.2c", "-6c", "0.2c", "0.2c"]):
        
        # --- SIGMA ---
        with fig3.set_panel(panel=[0, 0]):
            fig3.basemap(region=[-0.1, 0.1, 0.5, 1.1], projection=proj_lin, frame=["WSen", "xf0g0", "ya.2f.1+lσ"])
            # Plots dummy for first 2 points (PGA/PGV)
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rrup['sigma'].values[:2], style="d0.2c", pen='1.2p,red')
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rp['sigma'].values[:2], style="t0.2c", pen='1.2p,black')
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rp_unit['sigma'].values[:2], style="c0.2c", pen='1.2p,blue')
            fig3.text(text="a)", position='TL', offset='0.1/-0.1', font="10p,Helvetica,black", pen="0.25p,black,solid", fill='white')

        with fig3.set_panel(panel=[0, 1]):
            fig3.basemap(region=[0.009, 11, 0.5, 1.1], projection=proj_log_stat, frame=["wsen", "xa1g3", "ya.2f.1"])
            x, y = get_arrays(Stats_residual_Rrup, 'sigma')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,red', label='R@-RUP@-+N3')
            x, y = get_arrays(Stats_residual_Rp, 'sigma')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,black', label='R@-P@- w@-i@- = S@-i@-/ΣS@-i@-')
            x, y = get_arrays(Stats_residual_Rp_unit, 'sigma')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,blue', label='R@-P@- w@-i@- = 1')
            
            # Linea de diferencia
            # Need y-values at max_sigma period
            y_rup_max = Stats_residual_Rrup.loc[Stats_residual_Rrup['Period'] == max_sigma.Period, 'sigma'].values[0]
            y_rp_max = Stats_residual_Rp.loc[Stats_residual_Rp['Period'] == max_sigma.Period, 'sigma'].values[0]
            fig3.plot(x=[max_sigma.Period]*2, y=[y_rup_max, y_rp_max], pen='0.8p,black,-')
            fig3.text(x=max_sigma.Period, y=(y_rup_max + y_rp_max)/2, text=f"{round(max_sigma.sigma)}%", font="6p,Helvetica,black", pen="0.25p,black,solid", fill='white')
            fig3.legend(position="JBR+jBR+o1.5c/2.2c+w9c")

        # --- TAU ---
        with fig3.set_panel(panel=[1, 0]):
            fig3.basemap(region=[-0.1, 0.1, 0.3, 0.7], projection=proj_lin, frame=["WSen", "xf0g0", "ya.2f.1+lτ"])
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rrup['tau'].values[:2], style="d0.2c", pen='1.2p,red')
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rp['tau'].values[:2], style="t0.2c", pen='1.2p,black')
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rp_unit['tau'].values[:2], style="c0.2c", pen='1.2p,blue')
            fig3.text(text="b)", position='TL', offset='0.1/-0.1', font="10p,Helvetica,black", pen="0.25p,black,solid", fill='white')

        with fig3.set_panel(panel=[1, 1]):
            fig3.basemap(region=[0.009, 11, 0.3, 0.7], projection=proj_log_stat, frame=["wsen", "xa1g3", "ya.2f.1"])
            x, y = get_arrays(Stats_residual_Rrup, 'tau')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,red')
            x, y = get_arrays(Stats_residual_Rp, 'tau')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,black')
            x, y = get_arrays(Stats_residual_Rp_unit, 'tau')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,blue')
            
            y_rup_max = Stats_residual_Rrup.loc[Stats_residual_Rrup['Period'] == max_tau.Period, 'tau'].values[0]
            y_rp_max = Stats_residual_Rp.loc[Stats_residual_Rp['Period'] == max_tau.Period, 'tau'].values[0]
            fig3.plot(x=[max_tau.Period]*2, y=[y_rup_max, y_rp_max], pen='0.8p,black,-')
            fig3.text(x=max_tau.Period, y=(y_rup_max + y_rp_max)/2, text=f"{round(max_tau.tau)}%", font="6p,Helvetica,black", pen="0.25p,black,solid", fill='white')

        # --- PHI ---
        with fig3.set_panel(panel=[2, 0]):
            fig3.basemap(region=[-0.1, 0.1, 0.4, 0.8], projection=proj_lin, frame=["WSen", "xf0g0", "ya.2f.1+lφ"])
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rrup['phi'].values[:2], style="d0.2c", pen='1.2p,red')
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rp['phi'].values[:2], style="t0.2c", pen='1.2p,black')
            fig3.plot(x=[-0.05, 0.05], y=Stats_residual_Rp_unit['phi'].values[:2], style="t0.2c", pen='1.2p,blue')
            fig3.text(text="PGV", x=-0.05, y=0.34, font="8p,Helvetica,black", no_clip=True)
            fig3.text(text="PGA", x=0.05, y=0.34, font="8p,Helvetica,black", no_clip=True)
            fig3.text(text="c)", position='TL', offset='0.1/-0.1', font="10p,Helvetica,black", pen="0.25p,black,solid", fill='white')

        with fig3.set_panel(panel=[2, 1]):
            fig3.basemap(region=[0.009, 11, 0.4, 0.8], projection=proj_log_stat, frame=["wSen", "xa1g3", "ya.2f.1"])
            x, y = get_arrays(Stats_residual_Rrup, 'phi')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,red')
            x, y = get_arrays(Stats_residual_Rp, 'phi')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,black')
            x, y = get_arrays(Stats_residual_Rp_unit, 'phi')
            fig3.plot(x=x[2:], y=y[2:], pen='1.2p,blue')
            
            y_rup_max = Stats_residual_Rrup.loc[Stats_residual_Rrup['Period'] == max_phi.Period, 'phi'].values[0]
            y_rp_max = Stats_residual_Rp.loc[Stats_residual_Rp['Period'] == max_phi.Period, 'phi'].values[0]
            fig3.plot(x=[max_phi.Period]*2, y=[y_rup_max, y_rp_max], pen='0.8p,black,-')
            fig3.text(x=max_phi.Period, y=(y_rup_max + y_rp_max)/2, text=f"{round(max_phi.phi)}%", font="6p,Helvetica,black", pen="0.25p,black,solid", fill='white')
            fig3.text(text="Periods [s]", x=0.3, y=0.25, font="10p,Helvetica,black", no_clip=True)

    fig3.savefig(os.path.join(output_folder, 'Fig5.pdf'))

# =============================================================================
# 10. FIG 8 (p-values)
# =============================================================================

if not p_Rp.empty and not p_Rp_unit.empty:
    print("Generando Fig8 (p-values)...")
    
    # Preparar datos
    p_Rp_df = pd.DataFrame({
        'x': p_Rp['Period'][2:].values,
        'y': p_Rp['p_opt'][2:].values,
        'xerr': 0,
        'yerr': p_Rp['p_sd'][2:].values
    })
    
    p_Rp_unit_df = pd.DataFrame({
        'x': p_Rp_unit['Period'][2:].values,
        'y': p_Rp_unit['p_opt'][2:].values,
        'xerr': 0,
        'yerr': p_Rp_unit['p_sd'][2:].values
    })

    fig3 = pygmt.Figure()
    proj_lin = "X2.1c/8c"
    
    with fig3.subplot(nrows=1, ncols=2, figsize=("12c", "8c"), sharey=True, margins=["0.2c", "-6c", "0.2c", "0.2c"]):
        
        with fig3.set_panel(panel=[0, 0]):
            fig3.basemap(region=[-0.1, 0.1, -17, 2.5], projection=proj_lin, frame=["WSen", "xf0g0", "ya4f2+lPower p"])
            fig3.plot(x=[-0.05, 0.05], y=p_Rp['p_opt_pred'].values[:2], style="c0.1c", pen='0.1p,black', fill='red')
            fig3.plot(x=[-0.05, 0.05], y=p_Rp_unit['p_opt_pred'].values[:2], style="c0.1c", pen='0.1p,black', fill='blue')
            fig3.text(text="PGV", x=-0.05, y=-18.1, font="8p,Helvetica,black", no_clip=True)
            fig3.text(text="PGA", x=0.05, y=-18.1, font="8p,Helvetica,black", no_clip=True)
            
        with fig3.set_panel(panel=[0, 1]):
            fig3.basemap(region=[0.009, 11, -17, 2.5], projection=proj_log_stat, frame=["wSen", "xa1g3", "ya4f2"])
            fig3.plot(data=p_Rp_df, error_bar='+p0.75p,red', style='c0.1c', fill='red', pen="0.1p,black,solid", label='w@-i@- = S@-i@-/ΣS@-i@-+N4')
            fig3.plot(data=p_Rp_unit_df, error_bar='+p0.75p,blue', style='c0.1c', fill='blue', pen="0.1p,black,solid", label='w@-i@- = 1')
            fig3.plot(x=p_Rp['Period'][2:], y=p_Rp['p_opt_pred'][2:], pen='1.2p,red', label='Predicted')
            fig3.plot(x=p_Rp_unit['Period'][2:], y=p_Rp_unit['p_opt_pred'][2:], pen='1.2p,blue', label='Predicted')
            fig3.legend(position="JTR+jTR+o0.5c/-0.6c+w11c")
            fig3.text(text="Periods [s]", x=0.3, y=-19.4, font="12p,Helvetica,black", no_clip=True)

    fig3.savefig(os.path.join(output_folder, 'Fig8.pdf'))