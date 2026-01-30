import numpy as np
from scipy import signal

def crea_fuente(Mw, top_center, angulo_rotacion, Nasp, N_fallas,
                std_lon, std_lat, usar_celda_centrada=True):
    """
    Retorna un array (m, 3) con columnas [lon, lat, slip] para una fuente en Nasp asperezas.
    Se fuerza que el espaciamiento entre nodos sea el mismo en lon y lat (en grados)
    antes de la rotación.
    """

    # --- escalas sísmicas ---
    Mo = 10**(1.5*Mw + 9.1)
    mu = 30e9
    L = 10**(-2.9 + 0.63*Mw)    # km
    W = 10**(-0.86 + 0.35*Mw)   # km
    S = Mo / (mu * L * W * 1e6) # slip medio relativo (m)

    # conversión km -> grados (aprox igual en lon y lat para zona pequeña)
    km2deg = 1.0 / 111.11

    # --- tamaño de la malla por aspereza ---
    N_por_asp = max(1, int(round(N_fallas / Nasp)))
    Ls_km = L / Nasp                  # largo de cada strip (aspereza) en km
    lat_extent = Ls_km * km2deg       # extensión en grados de cada strip en lat
    lon_extent = W * km2deg           # extensión ideal en grados en lon

    # Tamaño de celda aproximado a partir del área / N_por_asp
    area_asp_km2 = Ls_km * W
    delta_km = np.sqrt(area_asp_km2 / N_por_asp)
    delta_deg_aprox = delta_km * km2deg  # paso angular "target"

    # Definimos ny a partir de la extensión en lat y delta_deg_aprox
    ny = max(1, int(round(lat_extent / delta_deg_aprox)))
    # Paso efectivo en lat
    delta_deg_lat = lat_extent / ny

    # Ahora forzamos que en lon el paso sea IGUAL que en lat
    delta_deg_lon = delta_deg_lat

    # nx a partir de la extensión ideal en lon y el delta forzado
    nx = max(1, int(round(lon_extent / delta_deg_lon)))
    # Ajustamos la extensión REAL en lon para que sea un múltiplo exacto
    lon_extent_real = nx * delta_deg_lon

    # --- preasignación ---
    m = Nasp * ny * nx
    fuente = np.zeros((m, 3), dtype=float)  # [lon, lat, slip]

    # --- utilidades para celdas centradas ---
    def linspace_centrado(a, b, n):
        if n == 1:
            return np.array([(a + b) / 2.0])
        if usar_celda_centrada:
            edges = np.linspace(a, b, n, endpoint=False)
            paso  = (b - a) / n
            return edges + 0.5 * paso
        else:
            return np.linspace(a, b, n)

    # --- bucle por aspereza ---
    idx0 = 0
    lon_pivot, lat_pivot = top_center

    for i in range(Nasp):
        # strip i a lo largo de L (en lat)
        y0 = lat_pivot + (i * Ls_km) * km2deg
        y1 = y0 + lat_extent  # y0 + Ls_km*km2deg

        # ancho W simétrico en lon alrededor del centro, pero usando lon_extent_real
        x0 = lon_pivot - lon_extent_real / 2.0
        x1 = lon_pivot + lon_extent_real / 2.0

        lat_vals = linspace_centrado(y0, y1, ny)
        lon_vals = linspace_centrado(x0, x1, nx)

        # ventana gaussiana separable (std_* en "número de muestras")
        win_lat = signal.windows.gaussian(ny, std=std_lat)[:, None]
        win_lon = signal.windows.gaussian(nx, std=std_lon)[None, :]
        grilla  = (win_lat * win_lon)
        grilla  = grilla / grilla.mean()
        grilla  = grilla.ravel()

        Lon, Lat = np.meshgrid(lon_vals, lat_vals, indexing="xy")
        Lon = Lon.ravel()
        Lat = Lat.ravel()

        # rotación (si tu coord_rotation acepta arrays, mejor vectorizar;
        # aquí lo dejo con bucle como en tu versión original)
        lon_rot = np.empty_like(Lon)
        lat_rot = np.empty_like(Lat)
        for j in range(Lon.size):
            lon_rot[j], lat_rot[j] = coord_rotation(
                lon_pivot, lat_pivot, Lon[j], Lat[j], angulo_rotacion
            )

        idx1 = idx0 + ny * nx
        fuente[idx0:idx1, 0] = lon_rot
        fuente[idx0:idx1, 1] = lat_rot
        fuente[idx0:idx1, 2] = grilla
        idx0 = idx1

    # escala de slip al promedio físico deseado
    fuente[:, 2] *= (S / fuente[:, 2].mean())

    return fuente





def coord_rotation(lon_polo,lat_polo,lon_punto,lat_punto,angulo_rotacion):
    
    '''
    Programa para rotar coordenadas geograficas con respecto a un eje definido por polo_lon
    y polo_lat, ademas del angulo de rotacion. Todo se debe expresar en grados
    '''

    lon_polo = np.radians(lon_polo)
    lat_polo = np.radians(lat_polo)
    c=lat_polo
    d=lon_polo
    rot_polo = np.radians(angulo_rotacion)

    polo_60 = np.array([np.cos(c)*np.cos(d), np.cos(c)*np.sin(d), np.sin(c)])

    a=np.radians(lat_punto)
    b=np.radians(lon_punto)
    punto=np.array([np.cos(a)*np.cos(b), np.cos(a)*np.sin(b),np.sin(a)])

    #Generamos la matriz de rotación

    co60=np.cos(np.radians(angulo_rotacion))
    se60=np.sin(np.radians(angulo_rotacion))

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

    x60 = np.dot(r60,punto)
    lo_60= np.degrees(np.arctan2(x60[1],x60[0]))
    la_60=np.degrees(np.arctan(x60[2]/np.sqrt(x60[0]**2+x60[1]**2)))

    return lo_60,la_60